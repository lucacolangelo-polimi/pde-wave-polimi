#ifndef WAVE_EQUATION_H
#define WAVE_EQUATION_H

// ============================================================
// WaveEquation.hpp
// Finite Element solver for the 2D wave equation:
//   u_tt - c^2 * Delta(u) = f(x,t)   in Omega
//   u = g                              on d(Omega)
//   u(t=0) = u0                        in Omega
//   u_t(t=0) = u1                      in Omega
//
// Features:
//   - Explicit Leapfrog (Stormer-Verlet) time integration
//   - Lumped Mass Matrix via Gauss-Lobatto quadrature (no linear solve!)
//   - Adaptive Mesh Refinement (AMR) with SolutionTransfer
//   - Absorbing Boundary Conditions (Sommerfeld / first-order ABC)
//   - Damping term (Damped Wave Equation: u_tt + d*u_t - c^2*Delta u = f)
//   - Time-dependent source term f(x,t)
//   - Convergence Analysis via Method of Manufactured Solutions (MMS)
//   - Discrete Energy tracking (kinetic + potential) for stability check
//   - CFL condition check at runtime
// ============================================================

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_q.h>
//#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values.h> 
//#include <deal.II/fe/fe_face_values.h> // Verifica che non ci siano errori di battitura
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <map>

using namespace dealii;


// Simulation mode selector
enum class SimulationMode
{
    PEBBLE_IN_POND,    // Gaussian IC, zero forcing, Dirichlet BCs
    MMS_CONVERGENCE,   // Method of Manufactured Solutions (standing wave)
    DAMPED_WAVE,       // Gaussian IC with viscous damping
    ABSORBING_BC,      // Gaussian IC with Sommerfeld absorbing BCs
};



// ============================================================
// Exact solution for MMS (Method of Manufactured Solutions)
// u_exact(x,y,t) = cos(pi*t) * sin(pi*x) * sin(pi*y)                   //stationary wave
// => u_tt = -pi^2 * cos(pi*t) * sin(pi*x) * sin(pi*y)
// => -c^2*Delta(u) = 2*pi^2*c^2 * cos(pi*t)*sin(pi*x)*sin(pi*y)
// => f = u_tt - c^2*Delta(u) = pi^2*(2*c^2 - 1)*cos(pi*t)*sin(pi*x)*sin(pi*y)
// ============================================================
template <int dim>
class ExactSolutionMMS : public Function<dim>
{
public:
    ExactSolutionMMS(double c_wave = 1.0)
        : Function<dim>(1), c(c_wave)
    {}

    virtual double value(const Point<dim> &p,
                         const unsigned int /*component*/ = 0) const override
    {
        return std::cos(M_PI * this->get_time()) *
               std::sin(M_PI * p[0]) *
               std::sin(M_PI * p[1]);
    }

    virtual Tensor<1, dim> gradient(const Point<dim> &p,
                                    const unsigned int /*component*/ = 0) const override
    {
        Tensor<1, dim> grad;
        grad[0] = std::cos(M_PI * this->get_time()) *
                  M_PI * std::cos(M_PI * p[0]) *
                  std::sin(M_PI * p[1]);
        grad[1] = std::cos(M_PI * this->get_time()) *
                  std::sin(M_PI * p[0]) *
                  M_PI * std::cos(M_PI * p[1]);
        return grad;
    }

private:
    double c;
};

// ============================================================
// Initial displacement for MMS: u0 = sin(pi*x)*sin(pi*y)
// ============================================================
template <int dim>
class InitialDisplacementMMS : public Function<dim>
{
public:
    virtual double value(const Point<dim> &p,
                         const unsigned int /*component*/ = 0) const override
    {
        return std::sin(M_PI * p[0]) * std::sin(M_PI * p[1]);
    }
};

// ============================================================
// Initial velocity for MMS: u1 = -pi * sin(pi*x)*sin(pi*y)
// ============================================================
template <int dim>
class InitialVelocityMMS : public Function<dim>                                     ////////!!!!!!!!!
{
public:
    virtual double value(const Point<dim> &p,
                         const unsigned int /*component*/ = 0) const override
    {
        return -M_PI * std::sin(M_PI * p[0]) * std::sin(M_PI * p[1]);
    }
};

// ============================================================
// Forcing term for MMS:
//   f = pi^2*(2*c^2 - 1) * cos(pi*t) * sin(pi*x) * sin(pi*y)
// ============================================================
template <int dim>
class ForcingTermMMS : public Function<dim>
{
public:
    ForcingTermMMS(double c_wave = 1.0)
        : Function<dim>(1), c(c_wave)
    {}

    virtual double value(const Point<dim> &p,
                         const unsigned int /*component*/ = 0) const override
    {
        const double factor = M_PI * M_PI * (2.0 * c * c - 1.0);
        return factor *
               std::cos(M_PI * this->get_time()) *
               std::sin(M_PI * p[0]) *
               std::sin(M_PI * p[1]);
    }

private:
    double c;
};

// Main WaveEquation class
template <int dim>
class WaveEquation
{
public:
    // Physical parameters
    double c        = 1.0;    // Wave speed [m/s]
    double damping  = 0.0;    // Damping coefficient d (0 = undamped)           ovvero quanto è viscosa la membrana 

    // Time parameters
    double time_step = 1.0e-3;
    double end_time  = 1.0;

    // Mesh parameters
    unsigned int initial_refinement = 6;  // Global refinement level                    mesh dinamica 
    unsigned int max_refinement_level = 8;// Maximum level for AMR
    unsigned int fe_degree = 1;           // Polynomial degree (1=linear, 2=quadratic)

    // Feature flags
    bool use_amr              = true;   // Adaptive Mesh Refinement
    bool use_absorbing_bc     = false;  // Sommerfeld ABC (replaces Dirichlet)
    bool track_energy         = true;   // Print discrete energy at each output step
    SimulationMode mode       = SimulationMode::PEBBLE_IN_POND;

    // AMR parameters
    unsigned int amr_every_n_steps     = 20;
    double       amr_refine_fraction   = 0.30;
    double       amr_coarsen_fraction  = 0.10;

    // Output parameters
    unsigned int output_every_n_steps  = 10;

    // Constructor / Run 
    WaveEquation();
    void run();
    void run_convergence_study(); // Runs MMS for multiple refinements

private:
    void make_grid();
    void setup_system();
    void assemble_matrices();         // Build M_lumped and K
    void assemble_rhs(double t);      // Build RHS vector (forcing + ABC contributions)
    void solve_time_step();
    void refine_mesh();
    void output_results(unsigned int step);
    void check_cfl_condition() const;

    // Energy
    double compute_kinetic_energy() const;
    double compute_potential_energy() const;

    // L2 and H1 errors (for MMS)
    std::pair<double,double> compute_errors(double t) const;

    // FEM data 
    Triangulation<dim>        triangulation;
    FE_Q<dim>                 fe;
    DoFHandler<dim>           dof_handler;
    AffineConstraints<double> constraints;
    SparsityPattern           sparsity_pattern;

    // Matrices
    Vector<double>       mass_matrix_diagonal;  // Lumped mass (diagonal)
    SparseMatrix<double> laplace_matrix;        // Stiffness matrix K
    // For absorbing BC: boundary mass matrix B (FEFaceValues integral)
    SparseMatrix<double> boundary_mass_matrix;  // integral of phi_i * phi_j on Gamma_abs

    // Solution vectors (Leapfrog needs 3 time levels)
    Vector<double> solution_u;       // u^n
    Vector<double> solution_u_old;   // u^{n-1}
    Vector<double> solution_u_new;   // u^{n+1}
    Vector<double> velocity_u;       // velocity estimate (u^n - u^{n-1}) / dt  [for damping & ABC]
    Vector<double> system_rhs;       // global RHS

    // Time tracking 
    double       time;
    unsigned int step_number;

    // Energy log 
    std::ofstream energy_log;

    // Convergence table (MMS) 
    ConvergenceTable convergence_table;
};

#endif // WAVE_EQUATION_H