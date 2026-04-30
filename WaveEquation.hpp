#ifndef WAVE_EQUATION_H
#define WAVE_EQUATION_H

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <memory>
using namespace dealii;

template <int dim>
class WaveEquation
{
public:
    
    WaveEquation();         // Constructor

    void run();

private:

    // Setup mesh and assemly
    void make_grid();       // Generate square or read from file
    void setup_system();    // Initialize arrays, vectors, and DoFHandlers
    void assemble_system(); // Construct the Mass (M) and Stiffness (K) Matrix
    void refine_mesh();     // Adaptive Mesh Refinement (AMR)

    // Time evolution 
    void solve_time_step(); // Calculate u_new using Leapfrog (M*a = RHS)
    void output_results(unsigned int step_number); // Output results for visualization

    // Mesh e FEM
    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;            // CHANGED FE_Q instead of FE_SimplexP 
    DoFHandler<dim>      dof_handler;
    AffineConstraints<double> constraints; // (Dirichlet boundary conditions)

    SparsityPattern      sparsity_pattern;
    Vector<double>       mass_matrix_diagonal; // Lumped Mass Matrix (diagonal only)
    SparseMatrix<double> laplace_matrix;  // Matrix K (Stiffness)

    // Vector solution (Leapfrog requires 3 steps)
    Vector<double>       solution_u;      // u^n     (Current)
    Vector<double>       solution_u_old;  // u^{n-1} (old)
    Vector<double>       solution_u_new;  // u^{n+1} (New)
    Vector<double>       system_rhs;      // Force F (RHS)

    // Time-stepping parameters
    double time;
    double time_step;
    double end_time;
    const unsigned int degree = 1;    // polynomial degree (linear)
    const double       c      = 1.0;  // Velocity of wave propagation
};

#endif