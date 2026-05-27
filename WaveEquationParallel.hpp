#ifndef WAVE_EQUATION_H
#define WAVE_EQUATION_H

// ============================================================
// WaveEquation.hpp  —  FEM solver MPI-parallelo per l'equazione
// delle onde 2D/3D con deal.II + Trilinos
//
//   u_tt − c²(x)·Δu + d·u_t = f(x,t)    in Ω
//   u = g                                  su ∂Ω
//   u(0) = u0,   u_t(0) = u1              in Ω
//
// Parallelismo MPI:
//   - parallel::distributed::Triangulation  : mesh distribuita (p4est)
//   - DoFHandler su mesh distribuita
//   - TrilinosWrappers::SparseMatrix        : matrici distribuite
//   - TrilinosWrappers::MPI::Vector         : vettori distribuiti
//   - SolverCG con precondizionatore SSOR   : per Newmark
//   - Ogni processo assembla solo le celle locally_owned
//   - Comunicazione implicita via Trilinos/MPI
//
// Schema temporale: Leapfrog esplicito (default) o Newmark-β implicito
// Massa: Lumped via Gauss-Lobatto
// AMR:   KellyErrorEstimator + parallel::distributed::SolutionTransfer
// ABC:   Sommerfeld (FEFaceValues)
// Validazione: MMS con ConvergenceTable, energia discreta
// ============================================================

// ---- deal.II distributed ----
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/grid_refinement.h>

// ---- Grid ----
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>

// ---- DoF ----
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

// ---- FE ----
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_face_values.h>
#include <deal.II/fe/mapping_q1.h>

// ---- LAC Trilinos (vettori e matrici MPI) ----
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector.h>

// ---- Numerics ----
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

// ---- Base ----
#include <deal.II/base/conditional_ostream.h>   // stampa solo da rank 0
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

// ============================================================
//  Shorthand per i tipi Trilinos MPI usati ovunque
// ============================================================
using TrilinosVector = TrilinosWrappers::MPI::Vector;
using TrilinosMatrix = TrilinosWrappers::SparseMatrix;

// ============================================================
//  Modalità di simulazione
// ============================================================
enum class SimulationMode
{
    PEBBLE_IN_POND,  // Gaussiana + Dirichlet + AMR
    MMS_CONVERGENCE, // Soluzione esatta (standing wave), convergence study
    DAMPED_WAVE,     // Gaussiana + smorzamento viscoso
    ABSORBING_BC,    // Gaussiana + condizioni assorbenti di Sommerfeld
    INTERFERENCE,    // Due sorgenti gaussiane (sovrapposizione)
    REFRACTION,      // c = c(x): mezzo eterogeneo
    DIFFRACTION,     // Ostacolo interno + fenditura
};

// ============================================================
//  Schema temporale
// ============================================================
enum class TimeScheme
{
    LEAPFROG,  // Esplicito O(dt²), richiede CFL. Nessun solve lineare.
    NEWMARK,   // Implicito Newmark-β, incondizionatamente stabile.
};

// ============================================================
//  Soluzione esatta MMS:  u = cos(πt)·sin(πx)·sin(πy)[·sin(πz)]
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
    double value(const Point<dim> &p, const unsigned int = 0) const override
    {
        double v = -M_PI * std::sin(M_PI * p[0]) * std::sin(M_PI * p[1]);
        if constexpr (dim == 3) v *= std::sin(M_PI * p[2]);
        return v;
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

// ============================================================
//  Classe principale WaveEquation (MPI-parallel)
// ============================================================
template <int dim>
class WaveEquation
{
public:
    // ---- Parametri pubblici ----
    double c             = 1.0;
    double damping       = 0.0;
    double time_step     = 1e-3;
    double end_time      = 1.0;

    double newmark_beta  = 0.25;   // schema trapezoidale
    double newmark_gamma = 0.50;

    unsigned int initial_refinement   = 6;
    unsigned int max_refinement_level = 8;
    unsigned int fe_degree            = 1;

    bool use_amr          = true;
    bool use_absorbing_bc = false;
    bool track_energy     = true;

    SimulationMode mode        = SimulationMode::PEBBLE_IN_POND;
    TimeScheme     time_scheme = TimeScheme::LEAPFROG;

    // Rifrazione
    double c_fast      = 2.0;
    double c_slow      = 0.8;
    double interface_y = 0.5;

    // AMR
    unsigned int amr_every_n_steps    = 20;
    double       amr_refine_fraction  = 0.30;
    double       amr_coarsen_fraction = 0.10;

    unsigned int output_every_n_steps = 10;

    // ---- Costruttore / Run ----
    WaveEquation(MPI_Comm mpi_communicator);
    ~WaveEquation() = default;

    void run();
    void run_convergence_study();

private:
    // ---- MPI ----
    MPI_Comm           mpi_comm;      // comunicatore MPI
    const unsigned int n_mpi_procs;   // numero totale di processi
    const unsigned int this_mpi_proc; // rank del processo corrente

    // Stampa solo dal processo 0
    ConditionalOStream pcout;

    // Timer per il profiling
    TimerOutput computing_timer;

    // ---- Metodi privati ----
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

    double compute_kinetic_energy()   const;
    double compute_potential_energy() const;
    std::pair<double, double> compute_errors(double t) const;

    double wave_speed_at(const Point<dim> &p) const;

    // ---- Dati FEM (tutti distribuiti via MPI) ----

    // Triangolazione distribuita (usa p4est internamente)
    parallel::distributed::Triangulation<dim> triangulation;

    std::unique_ptr<FE_Q<dim>> fe_ptr;
    DoFHandler<dim>            dof_handler;

    // IndexSet: quali DoF appartengono a questo processo
    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;  // owned + ghost

    AffineConstraints<double> constraints;

    // Matrici Trilinos distribuite
    TrilinosMatrix laplace_matrix;
    TrilinosMatrix boundary_mass_matrix;
    TrilinosMatrix system_matrix_newmark; // A = M + β·dt²·K (Newmark)

    // Massa lumpata: ogni processo tiene solo la sua porzione
    // Usiamo un TrilinosVector perché serve la comunicazione ghost
    TrilinosVector mass_matrix_diagonal;

    // Vettori soluzione (con ghost per la lettura inter-processo)
    TrilinosVector solution_u;          // u^n  (locally relevant)
    TrilinosVector solution_u_old;      // u^{n-1}
    TrilinosVector solution_u_new;      // u^{n+1}
    TrilinosVector velocity_u;          // v^n
    TrilinosVector acceleration_u;      // a^n (Newmark)
    TrilinosVector system_rhs;          // RHS (locally owned, no ghost)

    // Vettori "owned only" per l'update (senza ghost, per la scrittura)
    TrilinosVector owned_solution_u;
    TrilinosVector owned_solution_u_old;
    TrilinosVector owned_velocity_u;
    TrilinosVector owned_acceleration_u;

    double       time;
    unsigned int step_number;

    // Log energia (solo rank 0 scrive)
    std::ofstream energy_log;

    ConvergenceTable convergence_table;

    //Energy tracking
    double energy_initial = 0.0;          // E0 al passo t=0
    bool   energy_initialized = false;    // flag primo step

    // Norma L∞ della soluzione (rileva blow-up numerici)
    double compute_Linfty_norm() const;

    // Scrive il report finale di energia su file
    void write_energy_report() const;
};

#endif // WAVE_EQUATION_H
