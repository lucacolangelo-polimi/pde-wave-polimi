#include "WaveEquation.hpp"


#include <deal.II/grid/grid_generator.h>  
#include <deal.II/grid/grid_in.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/base/function_lib.h>

#include <iostream>
#include <fstream>

//Constructor
template <int dim>
WaveEquation<dim>::WaveEquation()
  : fe(1),                      // Initialization of FE_Q and FE_SimplexP of degree 1,  We are using first degree (linear) polynomials.
    dof_handler(triangulation), // DoFHandler to mesh
    time(0.0),                  // Initial time
    time_step(0.0)              // Temporal step 
{}


// ----------------------------------------------------------------------------
// BLOCK: Setup and Mesh
// ----------------------------------------------------------------------------
template <int dim>
void WaveEquation<dim>::make_grid()
{
    std::cout << "Generating grid..." << std::endl;

    // Option A: GridGenerator::hyper_cube(...) for the square
    // Generates a hypercube (square in 2D) from 0.0 to 1.0
    GridGenerator::hyper_cube(triangulation, 0.0, 1.0); 

    // Global refinement level
    // refine_global(6) in 2D means: 4^6 = 4096 cells.
    // Ideally suited for FE_Q elements.
    triangulation.refine_global(6);         

    std::cout << "   Number of active cells: " 
              << triangulation.n_active_cells() 
              << std::endl;
}

/*  ##WE NEED TO VERIFY THIS PART for void WaveEquation<dim>::make_grid()##
    // TODO: Implement mesh generation.
    // Option A: GridGenerator::hyper_cube(...) for the square
    // Option B: GridIn to read an external .msh file
    // Option C: refine_global() to refine the mesh
*/

template <int dim>
void WaveEquation<dim>::setup_system()
{
    std::cout << "Setting up system..." << std::endl;

    dof_handler.distribute_dofs(fe);

    std::cout << "   Number of degrees of freedom: " 
              << dof_handler.n_dofs() 
              << std::endl;
    constraints.clear();
    // lets apply Dirichlet BCs (u=0) on boundary_id = 0
    VectorTools::interpolate_boundary_values(dof_handler,
                                             0,
                                             Functions::ZeroFunction<dim>(),
                                             constraints);
    constraints.close(); 

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);

    sparsity_pattern.copy_from(dsp);

    mass_matrix_diagonal.reinit(dof_handler.n_dofs());
    laplace_matrix.reinit(sparsity_pattern);

    solution_u.reinit(dof_handler.n_dofs());
    solution_u_old.reinit(dof_handler.n_dofs());
    solution_u_new.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
}


// ----------------------------------------------------------------------------
// Matrices (M e K)
// ----------------------------------------------------------------------------
template <int dim>
void WaveEquation<dim>::assemble_system()
{
    std::cout << "Assembling matrices..." << std::endl;

    // Per FE_Q (quadrilaterals), we're using QGauss
    QGauss<dim> quadrature_formula(fe.degree + 1);
    FEValues<dim> fe_values(fe, quadrature_formula,
                            update_gradients | update_JxW_values);

    // For Mass Lumping, we use Gauss-Lobatto quadrature
    QGaussLobatto<dim> quadrature_formula_mass(fe.degree + 1);
    FEValues<dim> fe_values_mass(fe, quadrature_formula_mass,
                                 update_values | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points          = quadrature_formula.size();
    const unsigned int n_q_points_mass     = quadrature_formula_mass.size();

    FullMatrix<double> cell_laplace_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_mass_diagonal(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
        fe_values.reinit(cell);
        fe_values_mass.reinit(cell);
        
        cell_mass_diagonal = 0;
        cell_laplace_matrix = 0;

        // Assemble Stiffness Matrix using standard QGauss
        for (unsigned int q = 0; q < n_q_points; ++q)
        {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                    // Laplace Matrix: grad_phi_i * grad_phi_j
                    cell_laplace_matrix(i, j) += (fe_values.shape_grad(i, q) *
                                                  fe_values.shape_grad(j, q) *
                                                  fe_values.JxW(q));
                }
            }
        }

        // Assemble Lumped Mass Matrix using QGaussLobatto
        for (unsigned int q = 0; q < n_q_points_mass; ++q)
        {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                // With Gauss-Lobatto, shape functions are 1 at their node and 0 at others.
                // We only need to compute and store the diagonal entries.
                cell_mass_diagonal(i) += (fe_values_mass.shape_value(i, q) *
                                          fe_values_mass.shape_value(i, q) *
                                          fe_values_mass.JxW(q));
            }
        }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(cell_laplace_matrix, local_dof_indices, laplace_matrix);

        // Directly add local lumped mass to the global diagonal vector
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
            mass_matrix_diagonal(local_dof_indices[i]) += cell_mass_diagonal(i);
        }
    }
}

// ----------------------------------------------------------------------------
// BLOCK: Time evolution (Solver)
// ----------------------------------------------------------------------------
template <int dim>
void WaveEquation<dim>::solve_time_step()            //**let's verify other iterative methods but CG should be fine for symmetric positive definite matrices like M**
{
    // let's calculate the accelleration 'a': M * a = -c^2 * K * u_n
    // system_rhs = -c^2 * K * solution_u
    laplace_matrix.vmult(system_rhs, solution_u);
    system_rhs *= -(c * c);

    // Direct scalar division for lumped mass matrix: acc = system_rhs / M_diagonal
    Vector<double> acceleration(dof_handler.n_dofs());
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    {
        acceleration(i) = system_rhs(i) / mass_matrix_diagonal(i);
    }

    // Update u_new = 2*u - u_old + dt^2 * a
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    {
        solution_u_new(i) = 2.0 * solution_u(i) - solution_u_old(i) + 
                            (time_step * time_step) * acceleration(i);
    }

    //let's apply dirichlet boundary conditions
    constraints.distribute(solution_u_new);
}


template <int dim>
void WaveEquation<dim>::output_results(unsigned int step)
{
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution_u, "displacement");
    data_out.build_patches();

    std::ofstream output("solution-" + std::to_string(step) + ".vtk");
    data_out.write_vtk(output);
}

template <int dim>
void WaveEquation<dim>::run()
{
    std::cout << "Running simulation..." << std::endl;

    // 1. Setup
    make_grid();
    setup_system();
    assemble_system();

    // 2. Initial Conditions (The "Pebble in the Pond")
    std::cout << "Setting initial conditions..." << std::endl;
  
    const double amplitude = 1.0;
    const Point<dim> center(0.5, 0.5);
    const double width = 0.1;

    // Recuperiamo la posizione geometrica di ogni Grado di Libertà (DoF)
    std::vector<Point<dim>> support_points(dof_handler.n_dofs());
    DoFTools::map_dofs_to_support_points(MappingQ1<dim>(), dof_handler, support_points);

    // Calcoliamo il valore della Gaussiana per ogni nodo della mesh
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
    {
        double distance_sq = center.distance_square(support_points[i]);
        solution_u(i) = amplitude * std::exp(-distance_sq / (width * width));
    }

    // Velocità iniziale nulla: u_old = u_current
    solution_u_old = solution_u;

    // Applichiamo i vincoli al bordo (Dirichlet u=0)
    constraints.distribute(solution_u);
    constraints.distribute(solution_u_old);

    output_results(0);

    // 3. Temporal Loop
    time = 0.0;
    const double end_time = 1.0; 
    time_step = 0.001; // Assicurati che soddisfi CFL: dt < h/c
    
    unsigned int step = 0;
    while (time < end_time)
    {
        step++;
        time += time_step;

        solve_time_step();

        // Shift dei vettori per lo schema Leapfrog
        solution_u_old = solution_u;
        solution_u     = solution_u_new;
        
        if (step % 10 == 0) 
        {
            std::cout << "Step " << step << " at time " << time << std::endl;
            output_results(step);
        }
    }
    
    std::cout << "Simulation finished." << std::endl;
}
// ----------------------------------------------------------------------------
// BLOCK: Template Instantiation
// Necessary because we declare in .h and implement in .cc
// ----------------------------------------------------------------------------
template class WaveEquation<2>;