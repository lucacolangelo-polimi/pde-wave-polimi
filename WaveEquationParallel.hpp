#ifndef WAVE_EQUATIONPARALLEL_H                     //last version
#define WAVE_EQUATIONPARALLEL_H

#include <deal.II/distributed/tria.h>               // deal.II distributed 
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/grid/grid_generator.h>            // grid
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>               // DoF handler
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_q.h>                        // FE Q1, Q2, ...
#include <deal.II/fe/fe_values.h>
//#include <deal.II/fe/fe_face_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>     // Trilinos sparse matrix
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/vector_tools.h>          // numerics
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/base/conditional_ostream.h>       // prints only on rank 0
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/timer.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <map>

using namespace dealii;


using TrilinosVector = TrilinosWrappers::MPI::Vector;               /// Shorthand for Trilinos MPI vector
using TrilinosMatrix = TrilinosWrappers::SparseMatrix;

// SIMULATION MODES
enum class SimulationMode
{
    PEBBLE_IN_POND,  // Gaussian pulse + Dirichlet BC + AMR
    MMS_CONVERGENCE, // Exact solution (standing wave), convergence study
    DAMPED_WAVE,     // Gaussian pulse + viscous damping
    ABSORBING_BC,    // Gaussian pulse + Sommerfeld absorbing BCs
    INTERFERENCE,    // Two Gaussian sources (superposition)
    REFRACTION,      // c = c(x): heterogeneous medium
    DIFFRACTION,     // Internal obstacle + slit
};

// TEMPORAL SCHEMES
enum class TimeScheme
{
    LEAPFROG,  // Explicit O(dt2), CFL. 
    NEWMARK,   // ImplicitoNewmark-β, unconditionally stable.
};

// ============================================================
//  Exact solution MMS:  u = cos(πt)·sin(πx)·sin(πy)[·sin(πz)]
//  Forcing:  f = π²·(dim·c²−1)·cos(πt)·sin(πx)·sin(πy)[·sin(πz)]
// ============================================================
template <int dim>
class ExactSolutionMMS : public Function<dim>
{
public:
    ExactSolutionMMS(double c_wave = 1.0)
        : Function<dim>(1), c(c_wave) {}

    double value(const Point<dim> &p, const unsigned int = 0) const override
    {
        double v = std::cos(M_PI * this->get_time())
                 * std::sin(M_PI * p[0])
                 * std::sin(M_PI * p[1]);
        if constexpr (dim == 3) v *= std::sin(M_PI * p[2]);
        return v;
    }

    Tensor<1, dim> gradient(const Point<dim> &p, const unsigned int = 0) const override
    {
        const double tf = std::cos(M_PI * this->get_time());
        Tensor<1, dim> g;
        g[0] = tf * M_PI * std::cos(M_PI * p[0]) * std::sin(M_PI * p[1]);
        g[1] = tf * M_PI * std::sin(M_PI * p[0]) * std::cos(M_PI * p[1]);
        if constexpr (dim == 3)
        {
            g[0] *= std::sin(M_PI * p[2]);
            g[1] *= std::sin(M_PI * p[2]);
            g[2]  = tf * std::sin(M_PI * p[0]) * std::sin(M_PI * p[1])
                       * M_PI * std::cos(M_PI * p[2]);
        }
        return g;
    }
private:
    double c;
};

template <int dim>
class InitialDisplacementMMS : public Function<dim>
{
public:
    double value(const Point<dim> &p, const unsigned int = 0) const override
    {
        double v = std::sin(M_PI * p[0]) * std::sin(M_PI * p[1]);
        if constexpr (dim == 3) v *= std::sin(M_PI * p[2]);
        return v;
    }
};

template <int dim>
class InitialVelocityMMS : public Function<dim>
{
public:
    double value(const Point<dim> &/*p*/, const unsigned int = 0) const override
    {
        return 0.0; // The derivative of cos(pi*t) evaluated at t=0 is ZERO!
    }
};

template <int dim>
class ForcingTermMMS : public Function<dim>
{
public:
    ForcingTermMMS(double c_wave = 1.0) : Function<dim>(1), c(c_wave) {}

    double value(const Point<dim> &p, const unsigned int = 0) const override
    {
        const double factor = M_PI * M_PI * (static_cast<double>(dim) * c * c - 1.0);
        double v = factor * std::cos(M_PI * this->get_time())
                          * std::sin(M_PI * p[0])
                          * std::sin(M_PI * p[1]);
        if constexpr (dim == 3) v *= std::sin(M_PI * p[2]);
        return v;
    }
private:
    double c;
};

template <int dim>
class ExactSolutionMMS_Decay : public Function<dim>                                             //new decay solution
{
public:
    ExactSolutionMMS_Decay(double c_wave = 1.0) : Function<dim>(1), c(c_wave) {}
    virtual double value(const Point<dim> &p, const unsigned int = 0) const override
    {
        return std::exp(-this->get_time()) * std::sin(M_PI * p[0]) * std::sin(M_PI * p[1]);
    }
    virtual Tensor<1,dim> gradient(const Point<dim> &p, const unsigned int = 0) const override
    {
        const double common = std::exp(-this->get_time());
        Tensor<1,dim> g;
        g[0] = common * M_PI * std::cos(M_PI*p[0]) * std::sin(M_PI*p[1]);
        g[1] = common * M_PI * std::sin(M_PI*p[0]) * std::cos(M_PI*p[1]);
        return g;
    }
private:
    double c;
};

template <int dim>
class InitialDisplacementMMS_Decay : public Function<dim>                                           //new decay solution
{
public:
    virtual double value(const Point<dim> &p, const unsigned int = 0) const override
    {
        return std::sin(M_PI * p[0]) * std::sin(M_PI * p[1]);
    }
};

template <int dim>
class InitialVelocityMMS_Decay : public Function<dim>                                               //new decay solution
{
public:
    virtual double value(const Point<dim> &p, const unsigned int = 0) const override
    {
        return -std::sin(M_PI * p[0]) * std::sin(M_PI * p[1]);
    }
};

template <int dim>
class ForcingTermMMS_Decay : public Function<dim>                                                   //new decay solution
{
public:
    ForcingTermMMS_Decay(double c_wave = 1.0) : Function<dim>(1), c(c_wave) {}
    virtual double value(const Point<dim> &p, const unsigned int = 0) const override
    {
        const double factor = 1.0 + 2.0 * M_PI * M_PI * c * c;
        return factor * std::exp(-this->get_time()) * std::sin(M_PI * p[0]) * std::sin(M_PI * p[1]);
    }
private:
    double c;
};


template <int dim>
class WaveEquation
{
public:

    // struct for dispersion analysis results
    struct DispersionResult
    {
        double k;              // Tested wave number
        double kh;             // Dimensionless wave number
        double c_numerical;    // Measured numerical phase velocity
        double c_exact;        // Exact phase velocity (= c)
        double relative_error; // (c_h - c) / c
    };

    struct ScalingResult
    {
        unsigned int n_procs;
        unsigned int n_dofs;
        unsigned int n_cells;
        double wall_time_total;    // Total elapsed time [s]
        double wall_time_assembly; // Assembly time only [s]
        double wall_time_solve;    // Solve time only [s]
        double wall_time_output;   // Output/I/O time only [s]
        double speedup;            // Speedup relative to 1 process
        double efficiency;         // Speedup / n_procs
    };
    
    double c             = 1.0;     // public parameters 
    double damping       = 0.0;
    double time_step     = 1e-3;
    double end_time      = 1.0;

    double newmark_beta  = 0.25;   // trapezoidal scheme (unconditionally stable, O(dt2))
    double newmark_gamma = 0.50;   // 0.5 for no numerical damping, >0.5 for some damping

    unsigned int initial_refinement   = 6;      // 8^6 = 262144 cells in 2D, 8^4 = 4096 cells in 3D
    unsigned int max_refinement_level = 8;      // maximum AMR level (relative to initial_refinement)
    unsigned int fe_degree            = 1;      // Q1 elements, can change to Q2 for more accuracy

    bool use_amr          = true;
    bool use_absorbing_bc = false;
    bool track_energy     = true;

    SimulationMode mode        = SimulationMode::PEBBLE_IN_POND;          
    TimeScheme     time_scheme = TimeScheme::NEWMARK;              

    // Rifraction parameters (only for mode == REFRACTION)
    double c_fast      = 2.0;
    double c_slow      = 0.8;
    double interface_y = 0.5;

    // AMR
    unsigned int amr_every_n_steps    = 20;
    double       amr_refine_fraction  = 0.30;
    double       amr_coarsen_fraction = 0.10;

    bool use_decay_mms = false;                               //new decay solution

    unsigned int output_every_n_steps = 10;

    // Constructor/ run method
    WaveEquation(MPI_Comm mpi_communicator);
    ~WaveEquation() = default;

    void run();
    void run_convergence_study();
    void prepare_for_analysis()         // useful for analysis of dispersion
    {
        make_grid();
        setup_system();
        assemble_matrices();
    }

    // Performs n_steps time steps and measures performance (disables disk output to focus on raw compute time)
    ScalingResult run_scaling_benchmark(unsigned int n_steps);

        // Prints and saves the formatted scaling table
    void print_scaling_table(const std::vector<ScalingResult> &results,
                             const std::string &label) const;           //const!!!

    // Launches the analysis for a list of wave numbers and Returns a table of results (one per wave number)
    std::vector<DispersionResult> run_dispersion_analysis(const std::vector<double> &wave_numbers);

private:
    // MPI
    MPI_Comm           mpi_comm;      // MPI communicator
    const unsigned int n_mpi_procs;   // # of processes
    const unsigned int this_mpi_proc; // rank of this process

    // printing only on rank 0
    ConditionalOStream pcout;

    // Timer fot profiling
    TimerOutput computing_timer;

    // PRIVATE METHODS
    void make_grid();
    void make_grid_with_obstacle();
    void setup_system();
    void assemble_matrices();
    void assemble_rhs(double t);
    void solve_time_step();
    void solve_time_step_newmark();
    void refine_mesh();
    void output_results(unsigned int step);
    void check_cfl_condition() const;
    void perform_single_mms_run(unsigned int ref, double dt);

    double compute_kinetic_energy(); //  const;
    double compute_potential_energy(); //const;
    std::pair<double, double> compute_errors(double t) const;

    double wave_speed_at(const Point<dim> &p) const;

    //NEWMR MATRIX AND PRECONDITIONER 
    bool newmark_matrix_is_current = false;
    //TrilinosWrappers::PreconditionAMG newmark_preconditioner;                 //amr preconditioner, good for elliptic problems but non-optimal per Newmark (non-elliptic)
    TrilinosWrappers::PreconditionILU newmark_preconditioner;                   // ILU preconditioner, good for strongly diagonally dominant matrices like those from Newmark

    // Constructs A = M + β·dt²·K and initializes the AMG preconditioner.
    // Called after assemble_matrices() and after refine_mesh() (mesh changed).
    void build_newmark_system_matrix();


    // Measures the numerical ω for a single wave number k ia phase correlation <u_h(T), u_exact(T)>
    double measure_numerical_phase_speed(double k, double n_periods);


    //FEM data structures


    // Distributed triangulation and DoF handler
    parallel::distributed::Triangulation<dim> triangulation;

    std::unique_ptr<FE_Q<dim>> fe_ptr;
    DoFHandler<dim>            dof_handler;

    // IndexSet: locally_owned_dofs (writing) e locally_relevant_dofs (reading)
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;  // owned + ghost

    AffineConstraints<double> constraints;
    // Vincoli solo per hanging nodes, usati per la costruzione della matrice
    // per garantire la simmetria.
    AffineConstraints<double> matrix_constraints;

    // Distributed matrices (sparse, MPI-aware)
    TrilinosMatrix laplace_matrix;
    TrilinosMatrix boundary_mass_matrix;
    TrilinosMatrix system_matrix_newmark; // A = M + β·dt²·K (Newmark)

    // MASS LUMPED We use a TrilinosVector because we need ghost communication
    TrilinosVector mass_matrix_diagonal;

    // vector solution (with ghost values for reading)
    TrilinosVector solution_u;          // u^n  (locally relevant)
    TrilinosVector solution_u_old;      // u^{n-1}
    TrilinosVector solution_u_new;      // u^{n+1}
    TrilinosVector velocity_u;          // v^n
    TrilinosVector acceleration_u;      // a^n (Newmark)
    TrilinosVector system_rhs;          // RHS (locally owned, no ghost)

    // Vettore per il profilo dello Sponge Layer
    TrilinosVector nodal_sponge_profile;

    // "owned only" vectors for the update (no ghost, for writing)
    TrilinosVector owned_solution_u;
    TrilinosVector owned_solution_u_old;
    TrilinosVector owned_velocity_u;
    TrilinosVector owned_acceleration_u;

    double       time;             // current simulation time
    unsigned int step_number;      // current time step number

    // Log energy (only rank 0)
    std::ofstream energy_log;

    ConvergenceTable convergence_table;

    //Energy tracking
    double energy_initial = 0.0;          // E0 at step t=0
    bool   energy_initialized = false;    // flag for the first step

    // NormL∞ solution 
    double compute_Linfty_norm(); //const;

    void write_energy_report() const;
};


template <int dim>
class PlaneWave : public Function<dim>
{
public:
    PlaneWave(double wave_number, double wave_speed, double phase = 0.0)
        : Function<dim>(1)
        , k(wave_number)
        , c(wave_speed)
        , phi0(phase)
    {}

    double value(const Point<dim> &p, const unsigned int = 0) const override
    {
        // Onda che si propaga nella direzione x
        return std::sin(k * p[0] - c * k * this->get_time() + phi0);
    }

    // Velocità iniziale: u_t(x,0) = -c·k·cos(k·x + phi0)
    double time_derivative_at_zero(const Point<dim> &p) const
    {
        return -c * k * std::cos(k * p[0] + phi0);
    }

private:
    double k, c, phi0;
};
#endif 
