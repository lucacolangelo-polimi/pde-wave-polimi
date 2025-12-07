#include "../include/WaveEquation.h"

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>

#include <iostream>
#include <fstream>

// Constructor
template <int dim>
WaveEquation<dim>::WaveEquation()
  : fe(1) // TODO: P2 or P3?
{}

// ----------------------------------------------------------------------------
// BLOCK: Setup e Mesh
// ----------------------------------------------------------------------------
template <int dim>
void WaveEquation<dim>::make_grid()
{
    // TODO: Implement mesh generation.
    // Option A: GridGenerator::hyper_cube(...) for the square
    // Option B: GridIn to read an external .msh file
    // Option C: refine_global() to refine the mesh
}

template <int dim>
void WaveEquation<dim>::setup_system()
{
    // TODO: Distribute the DoFs (dof_handler.distribute_dofs)

    // TODO: Create the sparsity pattern (DynamicSparsityPattern)

    // TODO: Reinitialize the matrices (mass_matrix, laplace_matrix) using the pattern

    // TODO: Reinitialize the vectors (solution_u, u_old, u_new, rhs) to the appropriate size
}

// ----------------------------------------------------------------------------
// Matrices (M e K)
// ----------------------------------------------------------------------------
template <int dim>
void WaveEquation<dim>::assemble_system()
{
// TODO: Define quadrature (e.g., QGaussSimplex)
// TODO: Initialize FEValues ​​(with update_values, update_gradients, JxW)

// TODO: Loop through all cells (cell iterator)
// 1. Reinitialize fe_values ​​on the current cell
// 2. Loop through the quadrature points (q)
// 3. Loop through i (dofs) and j (dofs)
// - Calculate M_local: phi_i * phi_j * dx
// - Calculate K_local: grad_phi_i * grad_phi_j * dx
// 4. Add local contributions to the global matrices (mass_matrix, laplace_matrix)

// Note: The matrices are constant over time; they are assembled only once here.
}

// ----------------------------------------------------------------------------
// BLOCK: Time evolution (Solver)
// ----------------------------------------------------------------------------
template <int dim>
void WaveEquation<dim>::solve_time_step()
{
// TODO: Construct the known term (RHS)
// Formula: RHS = -c^2 * K * u_n (use vmult)

// TODO: Solve the linear system M * a = RHS
// Use SolverCG and PreconditionSSOR to find the acceleration 'a'

// TODO: Update u_new (Leapfrog)
// Formula: u_new = 2*u_n - u_old + dt^2 * a

// TODO: Apply Dirichlet constraints (if necessary, u=0 at the boundary)
}

template <int dim>
void WaveEquation<dim>::output_results(unsigned int step)
{
    // TODO: DataOut in order to create a file VTK
    // attach_dof_handler, add_data_vector(solution_u), build_patches, write_vtk
}

template <int dim>
void WaveEquation<dim>::run()
{
    // ------------------------------------------------------------------------
    // MAIN LOOP
    // ------------------------------------------------------------------------
    
    // 1. Preparation
    make_grid();
    setup_system();
    assemble_system(); // Assembly of M and K
    
    // 2. Initial Conditions
    // TODO: Set u_old and u_current (e.g., VectorTools::interpolate or manually)
    
    // 3. Temporal Loop
    time = 0.0;
    time_step = 0.001; /
    
    unsigned int step = 0;
    while (time < 1.0) // TODO: end_time
    {
        solve_time_step();
        
        // Shift 
        solution_u_old = solution_u;
        solution_u = solution_u_new;
        
        // Output every x steps (10 for now)
        if (step % 10 == 0) {
            output_results(step);
        }
        
        time += time_step;
        step++;
    }
}

// ----------------------------------------------------------------------------
// BLOCK: Template Instantiation
// Necessary because we declare in .h and implement in .cc
// ----------------------------------------------------------------------------
template class WaveEquation<2>;