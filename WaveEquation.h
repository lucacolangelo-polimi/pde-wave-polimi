#ifndef WAVE_EQUATION_H
#define WAVE_EQUATION_H

// ----------------------------------------------------------------------------
// BLOCK: Include
// let's include necessary deal.II headers
// ----------------------------------------------------------------------------
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>                //Added for FE_Q
// #include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/affine_constraints.h>

using namespace dealii;

// ----------------------------------------------------------------------------
// BLOCK: Classe Template
// <int dim> for 2D or 3D problems
// ----------------------------------------------------------------------------
template <int dim>
class WaveEquation
{
public:
    // Constructor
    WaveEquation();

    //"Driver" method: manages the entire simulation flow
    void run();

private:
    // ------------------------------------------------------------------------
    // SECTION1: Setup  and mesh 
    // ------------------------------------------------------------------------
    void make_grid();       // Generate square or read from file
    void setup_system();    // Initialize arrays, vectors, and DoFHandlers

    // ------------------------------------------------------------------------
    // SECTION2 : Assembly Matrices
    // ------------------------------------------------------------------------
    void assemble_system(); // Construct the Mass (M) and Stiffness (K) Matrix

    // ------------------------------------------------------------------------
    // SECTION 3: Time evolution 
    // ------------------------------------------------------------------------
    void solve_time_step(); // Calculate u_new using Leapfrog (M*a = RHS)
    void output_results(unsigned int step_number); // Output results for visualization

    // Mesh e FEM
    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;            // CHANGED FE_Q instead of FE_SimplexP 
    DoFHandler<dim>      dof_handler;
    AffineConstraints<double> constraints; // (Dirichlet boundary conditions)

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> mass_matrix;     // Matrix M
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
};

#endif