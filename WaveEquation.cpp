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
#include <limits>
#include <algorithm>

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
    if (dim == 2)
    {
        triangulation.refine_global(6);
    }
    else if (dim == 3)
    {
        triangulation.refine_global(5); // 8^5 = 32768 cells, to keep the problem size manageable in 3D
    }

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

    mass_matrix.reinit(sparsity_pattern);
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
                            update_values | update_gradients | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_laplace_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
        fe_values.reinit(cell);
        cell_mass_matrix = 0;
        cell_laplace_matrix = 0;

        for (unsigned int q = 0; q < n_q_points; ++q)
        {
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                {
                    // mass matrix : phi_i * phi_j
                    cell_mass_matrix(i, j) += (fe_values.shape_value(i, q) *
                                               fe_values.shape_value(j, q) *
                                               fe_values.JxW(q));

                    // Laplace Matrix: grad_phi_i * grad_phi_j
                    cell_laplace_matrix(i, j) += (fe_values.shape_grad(i, q) *
                                                  fe_values.shape_grad(j, q) *
                                                  fe_values.JxW(q));
                }
            }
        }
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(cell_mass_matrix, local_dof_indices, mass_matrix);
        constraints.distribute_local_to_global(cell_laplace_matrix, local_dof_indices, laplace_matrix);
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

    // Solving M * acc = system_rhs
    Vector<double> acceleration(dof_handler.n_dofs());
    SolverControl solver_control(1000, 1e-12 * system_rhs.l2_norm());           
    SolverCG<Vector<double>> solver(solver_control);
    
    // PreconditionIdentity because Mass matrix is well conditioned (**let's verify this assumption**)
    solver.solve(mass_matrix, acceleration, system_rhs, PreconditionIdentity());

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

   std::string filename = "solution-" + 
                           std::to_string(dim) + "d-" + 
                           std::to_string(step) + ".vtu";

    std::ofstream output(filename);
    data_out.write_vtu(output);
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

    Point<dim> center;                          //doing so, the generation of the starting point is independent of the dimension of the problem, and we can easily change the dimension without worrying about the initial condition generation.
    for (unsigned int d = 0; d < dim; ++d) {
    center(d) = 0.5; // Centro del dominio [0,1]^dim
}
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
    // ---------------------------------------------------------
    // CALCOLO DINAMICO DEL PASSO TEMPORALE (Condizione CFL)
    // ---------------------------------------------------------
    std::cout << "Calcolo del passo temporale per la stabilità CFL..." << std::endl;

    double h_min = std::numeric_limits<double>::max();
    
    // Iteriamo su tutte le celle attive per trovare la minima distanza tra i vertici
    for (const auto &cell : triangulation.active_cell_iterators())
    {
        h_min = std::min(h_min, cell->minimum_vertex_distance());
    }

    // Impostiamo il Numero di Courant. 
    // Valori tipici per Leapfrog + FE_Q(1) + Massa Consistente sono tra 0.1 e 0.5.
    const double courant_number = 0.2; 

    // Calcolo del time_step: dt = C * (h_min / c)
    time_step = courant_number * h_min / c;

    std::cout << "   Distanza minima h_min: " << h_min << std::endl;
    std::cout << "   Passo temporale dt: " << time_step << std::endl;
    // ---------------------------------------------------------
    
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
template class WaveEquation<3>; //If you want to change dimensions, just change this line and recompile
template class WaveEquation<2>;