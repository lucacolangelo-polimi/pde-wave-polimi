// ============================================================
// WaveEquation.cpp
//
// Finite Element solver for the 2D wave equation.
//
// Time scheme: Explicit Leapfrog (Stormer-Verlet)
//   u^{n+1} = 2*u^n - u^{n-1} + dt^2 * M^{-1} * RHS^n
// Where RHS^n = -c^2 * K * u^n + F^n - ABC terms - damping terms
//
// With lumped mass M is diagonal => no linear solve needed.
//
// Stability (CFL): dt < h / (c * sqrt(2))   [2D, Q1 elements]
// ============================================================

#include "WaveEquationOptimized.hpp"

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/base/utilities.h>

#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <numeric>


// Constructor
template <int dim>
WaveEquation<dim>::WaveEquation()
    : fe(fe_degree),
      dof_handler(triangulation),
      time(0.0),
      step_number(0)
{}

// make_grid
template <int dim>
void WaveEquation<dim>::make_grid()
{
    std::cout << "  Generating grid (hyper_cube [0,1]^" << dim
              << ", " << initial_refinement << " global refinements)..." << std::endl;

    GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
    triangulation.refine_global(initial_refinement);

    std::cout << "  Active cells: " << triangulation.n_active_cells() << std::endl;
}

// setup_system
template <int dim>
void WaveEquation<dim>::setup_system()
{
    std::cout << "  Setting up system..." << std::endl;

    dof_handler.distribute_dofs(fe);
    std::cout << "  Degrees of freedom: " << dof_handler.n_dofs() << std::endl;

    // Constraints (Dirichlet u=0 on all boundary, unless ABC is active)            vincoli sui bordi 
    constraints.clear();
    if (!use_absorbing_bc)
    {
        // Dirichlet BC on boundary_id = 0 (all faces of hyper_cube)
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 0,
                                                 Functions::ZeroFunction<dim>(),
                                                 constraints);
    }
    // For MMS mode: enforce exact solution on boundary
    if (mode == SimulationMode::MMS_CONVERGENCE)
    {
        // Will be applied at each time step during RHS assembly using the ExactSolution
        // Boundary condition for MMS is u=0 (sin(pi*0)*sin(pi*y)=0), so Dirichlet zero is correct.
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 0,
                                                 Functions::ZeroFunction<dim>(),
                                                 constraints);
    }
    constraints.close();

    // Sparsity pattern                                                         //creates a map with the non-zero entries of the matrix, based on the DoF connectivity and constraints. This is needed to efficiently allocate memory for the sparse matrix.
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    sparsity_pattern.copy_from(dsp);

    // Allocate matrices and vectors 
    mass_matrix_diagonal.reinit(dof_handler.n_dofs());
    laplace_matrix.reinit(sparsity_pattern);
    if (use_absorbing_bc)
        boundary_mass_matrix.reinit(sparsity_pattern);

    solution_u.reinit(dof_handler.n_dofs());
    solution_u_old.reinit(dof_handler.n_dofs());
    solution_u_new.reinit(dof_handler.n_dofs());
    velocity_u.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
}

// ============================================================
// assemble_matrices
//   - Lumped Mass Matrix (QGaussLobatto) -> mass_matrix_diagonal
//   - Stiffness Matrix K (QGauss)        -> laplace_matrix
//   - Boundary Mass Matrix (FEFaceValues) -> boundary_mass_matrix (only if ABC)
// ============================================================
template <int dim>
void WaveEquation<dim>::assemble_matrices()
{
    std::cout << "  Assembling matrices..." << std::endl;

    // Volume quadrature 
    QGauss<dim>        q_stiffness(fe.degree + 1);                  // For stiffness, we need to integrate gradients, so standard Gauss is fine.
    QGaussLobatto<dim> q_mass(fe.degree + 1);                       // For mass lumping, we use Gauss-Lobatto to get diagonal mass matrix.

    FEValues<dim> fev_stiff(fe, q_stiffness,                            // for stiffness K_ij = integral grad(phi_i) . grad(phi_j)
                            update_gradients | update_JxW_values);
    FEValues<dim> fev_mass(fe, q_mass,                                  // for mass lumping M_ii = sum_q phi_i(q)^2 * w_q
                           update_values | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;                // number of local DoFs per cell (e.g. 4 for Q1 in 2D)

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);       // local stiffness matrix
    Vector<double>     cell_mass_diag(dofs_per_cell);                   
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);  // For ABC: we need to compute the boundary mass matrix B_ij = integral_boundary phi_i * phi_j

    // Face quadrature for ABC 
    std::unique_ptr<FEFaceValues<dim>> fev_face;                            // Only needed if absorbing BC is active
    QGauss<dim - 1> q_face(fe.degree + 1);                                  // For ABC: we need to integrate phi_i * phi_j on the boundary faces, so we use FEFaceValues with a suitable quadrature.
    FullMatrix<double> cell_boundary_matrix(dofs_per_cell, dofs_per_cell);  // Local matrix for boundary integrals on each face

    if (use_absorbing_bc)                                                   // Initialize FEFaceValues for boundary mass matrix assembly (only if ABC is active)
    {
        fev_face = std::make_unique<FEFaceValues<dim>>(                     // for boundary integrals in ABC: B_ij = integral_boundary phi_i * phi_j
            fe, q_face,
            update_values | update_JxW_values);
    }

    for (const auto &cell : dof_handler.active_cell_iterators())            // Loop over all active cells to assemble local contributions to the global matrices
    {
        fev_stiff.reinit(cell);                                             // Compute local stiffness matrix K_ij = integral grad(phi_i) . grad(phi_j)
        fev_mass.reinit(cell);                                              // Compute local mass contributions for lumping: M_ii = sum_q phi_i(q)^2 * w_q

        cell_matrix    = 0.0;
        cell_mass_diag = 0.0;

        // Stiffness matrix (K_ij = integral grad(phi_i) . grad(phi_j)) 
        for (unsigned int q = 0; q < q_stiffness.size(); ++q)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    cell_matrix(i, j) +=
                        fev_stiff.shape_grad(i, q) *
                        fev_stiff.shape_grad(j, q) *
                        fev_stiff.JxW(q);

        //  Lumped mass (Gauss-Lobatto: phi_i(q_j) = delta_ij for Q1) 
        // We accumulate M_ii = sum_q phi_i(q)^2 * w_q
        // This is exact mass lumping for tensor-product elements.
        for (unsigned int q = 0; q < q_mass.size(); ++q)                // for integratio points
            for (unsigned int i = 0; i < dofs_per_cell; ++i)            // for nodes of the cell
                cell_mass_diag(i) +=
                    fev_mass.shape_value(i, q) *                        // value of shape function i at quadrature point q
                    fev_mass.shape_value(i, q) *                        // value of shape function i at quadrature point q (squared for lumping)
                    fev_mass.JxW(q);                                    // quadrature weight * det(Jacobian) at point q

        cell->get_dof_indices(local_dof_indices);                       // for mapping the local nodes

        // Distribute stiffness matrix (respects constraints)
        constraints.distribute_local_to_global(
            cell_matrix, local_dof_indices, laplace_matrix);

        // Accumulate lumped mass (global diagonal)
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
            mass_matrix_diagonal(local_dof_indices[i]) += cell_mass_diag(i);

        // Boundary mass matrix for ABC 
        // ABC (Sommerfeld): add  c * integral_boundary (phi_i * phi_j) * u_t
        // => introduces a flux term: -c * B * v  where v = velocity
        if (use_absorbing_bc)
        {
            for (const auto &face : cell->face_iterators())
            {
                if (face->at_boundary())            //calcola l'integrale solo nsulla faccia che tocca il bordo
                {
                    fev_face->reinit(cell, face);   // Compute local boundary mass matrix B_ij = integral_boundary phi_i * phi_j for this face
                    cell_boundary_matrix = 0.0;
                    for (unsigned int q = 0; q < q_face.size(); ++q)
                        for (unsigned int i = 0; i < dofs_per_cell; ++i)
                            for (unsigned int j = 0; j < dofs_per_cell; ++j)
                                cell_boundary_matrix(i, j) +=
                                    fev_face->shape_value(i, q) *
                                    fev_face->shape_value(j, q) *
                                    fev_face->JxW(q);

                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        for (unsigned int j = 0; j < dofs_per_cell; ++j)
                            boundary_mass_matrix.add(
                                local_dof_indices[i],
                                local_dof_indices[j],
                                cell_boundary_matrix(i, j));
                }
            }
        }
    }

    // Verify mass matrix is positive (sanity check)
    const double min_mass = *std::min_element(
        mass_matrix_diagonal.begin(), mass_matrix_diagonal.end());
    AssertThrow(min_mass > 0.0,
        ExcMessage("Lumped mass matrix has non-positive diagonal entry!"));
}


// assemble_rhs
//   Computes: RHS = -c^2 * K * u^n  +  F^n  -  ABC  -  damping
template <int dim>
void WaveEquation<dim>::assemble_rhs(double t)
{
    system_rhs = 0.0;

    // Term 1: -c^2 * K * u^n
    laplace_matrix.vmult(system_rhs, solution_u);
    system_rhs *= -(c * c);

    // Term 2: forcing term f(x,t) for MMS mode          // if we are in MMs mode we need tp add an artificial forcing term to the RHS to ensure that the manufactured solution is an exact solution of the PDE. 
    if (mode == SimulationMode::MMS_CONVERGENCE)
    {
        QGauss<dim> quadrature_formula(fe.degree + 1);
        FEValues<dim> fe_values(fe, quadrature_formula,
                                update_values | update_quadrature_points | update_JxW_values);

        ForcingTermMMS<dim> forcing_term(c);
        forcing_term.set_time(t);

        const unsigned int dofs_per_cell = fe.dofs_per_cell;
        Vector<double> cell_rhs(dofs_per_cell);
        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

        for (const auto &cell : dof_handler.active_cell_iterators())
        {
            fe_values.reinit(cell);
            cell_rhs = 0.0;

            for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
            {
                const double f_val = forcing_term.value(
                    fe_values.quadrature_point(q));
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    cell_rhs(i) += fe_values.shape_value(i, q) * f_val * fe_values.JxW(q);
            }

            cell->get_dof_indices(local_dof_indices);
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                system_rhs(local_dof_indices[i]) += cell_rhs(i);
        }
    }

    // Term 3: Absorbing BC  -c * B * v (Sommerfeld, first-order) 
    // velocity_u = (u^n - u^{n-1}) / dt  (centered difference)
    if (use_absorbing_bc)
    {
        Vector<double> abc_term(dof_handler.n_dofs());
        boundary_mass_matrix.vmult(abc_term, velocity_u);
        abc_term *= -c;
        system_rhs.add(1.0, abc_term);
    }

    // Term 4: Damping  -d * M * v  (viscous damping)
    // Adds: -damping * M * v  to RHS  (note: M is diagonal)
    if (damping > 0.0)
    {
        for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
            system_rhs(i) -= damping * mass_matrix_diagonal(i) * velocity_u(i);
    }
}

// ============================================================
// solve_time_step
//   Leapfrog update: u^{n+1} = 2*u^n - u^{n-1} + dt^2 * M^{-1} * RHS
//   With lumped mass: M^{-1} is just element-wise division.
// ============================================================
template <int dim>
void WaveEquation<dim>::solve_time_step()
{
    // Assemble RHS at current time
    assemble_rhs(time);

    // acceleration = M^{-1} * RHS  (element-wise, M is diagonal)
    Vector<double> acceleration(dof_handler.n_dofs());
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
        acceleration(i) = system_rhs(i) / mass_matrix_diagonal(i);

    // Leapfrog: u^{n+1} = 2*u^n - u^{n-1} + dt^2 * a
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
        solution_u_new(i) = 2.0 * solution_u(i)
                          - solution_u_old(i)
                          + (time_step * time_step) * acceleration(i);

    // Apply Dirichlet BCs (zero out constrained dofs)
    if (!use_absorbing_bc)
        constraints.distribute(solution_u_new);

    // Update velocity estimate: v^n = (u^{n+1} - u^{n-1}) / (2*dt)
    // Used for damping and ABC at next step.
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
        velocity_u(i) = (solution_u_new(i) - solution_u_old(i)) / (2.0 * time_step);

    // Shift solution vectors
    solution_u_old = solution_u;
    solution_u     = solution_u_new;
}


// Energy tracking
template <int dim>
double WaveEquation<dim>::compute_kinetic_energy() const
{
    // E_k = 0.5 * v^T * M * v  (with lumped M)
    double ek = 0.0;
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
        ek += 0.5 * mass_matrix_diagonal(i) * velocity_u(i) * velocity_u(i);
    return ek;
}

template <int dim>
double WaveEquation<dim>::compute_potential_energy() const
{
    // E_p = 0.5 * c^2 * u^T * K * u
    Vector<double> Ku(dof_handler.n_dofs());
    laplace_matrix.vmult(Ku, solution_u);
    double ep = 0.0;
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
        ep += 0.5 * c * c * solution_u(i) * Ku(i);
    return ep;
}

// Error computation for MMS
// Returns {L2_error, H1_error}
template <int dim>
std::pair<double, double> WaveEquation<dim>::compute_errors(double t) const
{
    ExactSolutionMMS<dim> exact_solution(c);
    exact_solution.set_time(t);

    Vector<double> difference_per_cell(triangulation.n_active_cells());

    // L2 error
    VectorTools::integrate_difference(dof_handler,
                                      solution_u,
                                      exact_solution,
                                      difference_per_cell,
                                      QGauss<dim>(fe.degree + 2),
                                      VectorTools::L2_norm);
    const double L2_error =
        VectorTools::compute_global_error(triangulation,
                                          difference_per_cell,
                                          VectorTools::L2_norm);

    // H1 (semi-norm gradient) error
    VectorTools::integrate_difference(dof_handler,
                                      solution_u,
                                      exact_solution,
                                      difference_per_cell,
                                      QGauss<dim>(fe.degree + 2),
                                      VectorTools::H1_seminorm);
    const double H1_error =
        VectorTools::compute_global_error(triangulation,
                                          difference_per_cell,
                                          VectorTools::H1_seminorm);

    return {L2_error, H1_error};
}

// CFL check                                                                // Courant-Friedrichs-Lewy, È un test matematico che serve a capire se il passo temporale scelto è troppo grande rispetto alla dimensione dei quadratini della griglia 
template <int dim>
void WaveEquation<dim>::check_cfl_condition() const
{
    double h_min = std::numeric_limits<double>::max();
    for (const auto &cell : dof_handler.active_cell_iterators())
        h_min = std::min(h_min, cell->minimum_vertex_distance());

    // CFL: dt < h / (c * sqrt(dim))
    const double cfl_limit = h_min / (c * std::sqrt(static_cast<double>(dim)));         // limite teorico per la stabilità del metodo esplicito, basato sulla velocità di propagazione c e sulla dimensione dei celle h. Se il passo temporale dt è maggiore di questo limite, la simulazione potrebbe diventare instabile.
    const double cfl_number = time_step / cfl_limit;

    std::cout << "  CFL check: dt=" << time_step
              << "  h_min=" << h_min
              << "  CFL_limit=" << cfl_limit
              << "  CFL_number=" << cfl_number;             //must be < 1 for stability

    if (cfl_number > 1.0)
        std::cout << "  *** WARNING: CFL VIOLATED! Simulation may be unstable! ***";
    else
        std::cout << "  [OK]";
    std::cout << std::endl;
}


// output_results
template <int dim>
void WaveEquation<dim>::output_results(unsigned int step)
{
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution_u, "displacement");
    data_out.add_data_vector(velocity_u, "velocity");
    data_out.build_patches();

    const std::string filename = "solution-" + std::to_string(step) + ".vtu";
    std::ofstream output(filename);
    data_out.write_vtu(output);
}

// refine_mesh  (AMR with SolutionTransfer)
template <int dim>
void WaveEquation<dim>::refine_mesh()
{
    std::cout << "  AMR: estimating errors and refining mesh..." << std::endl;

    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate(                 //to know if the cell must be refined or coarsened, we need an error estimator. Here we use the KellyErrorEstimator, which is a common choice for elliptic problems and can be adapted for wave equations by looking at the solution at the current time step. The estimator looks at the jump in the gradient of the solution across cell interfaces, which is a good indicator of where the solution is not well resolved.
        dof_handler,
        QGauss<dim - 1>(fe.degree + 1),
        std::map<types::boundary_id, const Function<dim> *>(),
        solution_u,
        estimated_error_per_cell);

    GridRefinement::refine_and_coarsen_fixed_fraction(
        triangulation, estimated_error_per_cell,
        amr_refine_fraction, amr_coarsen_fraction);

    // Enforce maximum refinement level (to keep CFL valid)
    for (const auto &cell : triangulation.active_cell_iterators())      //we are looking for a limit on the maximum refinement level to prevent the mesh from becoming too fine, which could violate the CFL condition and lead to instability. If a cell is marked for refinement but its current level is already at or above the specified maximum, we clear the refine flag for that cell, effectively preventing it from being refined further.
        if (cell->level() >= static_cast<int>(max_refinement_level))
            cell->clear_refine_flag();

    triangulation.prepare_coarsening_and_refinement();

    // Transfer solution_u and solution_u_old to the new mesh
    SolutionTransfer<dim, Vector<double>> solution_transfer(dof_handler);       //solution transfer si occupa di "proiettare" matematicamente l'onda dalla vecchia mesh a quella nuova
    std::vector<Vector<double>> x_vectors = {solution_u, solution_u_old, velocity_u};
    solution_transfer.prepare_for_coarsening_and_refinement(x_vectors);

    triangulation.execute_coarsening_and_refinement();

    // Re-distribute DoFs
    dof_handler.distribute_dofs(fe);
    std::cout << "  AMR: new DoFs = " << dof_handler.n_dofs()
              << "  cells = " << triangulation.n_active_cells() << std::endl;

    // Rebuild constraints
    constraints.clear();
    if (!use_absorbing_bc)
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 0,
                                                 Functions::ZeroFunction<dim>(),
                                                 constraints);
    constraints.close();

    // Rebuild sparsity and matrices
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    sparsity_pattern.copy_from(dsp);

    mass_matrix_diagonal.reinit(dof_handler.n_dofs());
    laplace_matrix.reinit(sparsity_pattern);
    if (use_absorbing_bc)
        boundary_mass_matrix.reinit(sparsity_pattern);

    // Interpolate old solutions onto new mesh
    std::vector<Vector<double>> tmp(3, Vector<double>(dof_handler.n_dofs()));
    solution_transfer.interpolate(x_vectors, tmp);

    solution_u     = tmp[0];
    solution_u_old = tmp[1];
    velocity_u     = tmp[2];
    solution_u_new.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.distribute(solution_u);
    constraints.distribute(solution_u_old);

    assemble_matrices();

    // Re-check CFL after refinement
    check_cfl_condition();
}

// run  (main time loop for standard simulation modes)
template <int dim>
void WaveEquation<dim>::run()
{
    std::cout << "=======================================" << std::endl;
    std::cout << "  WaveEquation<" << dim << "> simulation" << std::endl;
    std::cout << "  Mode: ";
    switch (mode)
    {
        case SimulationMode::PEBBLE_IN_POND:   std::cout << "Pebble in Pond";   break;
        case SimulationMode::MMS_CONVERGENCE:  std::cout << "MMS Convergence";  break;
        case SimulationMode::DAMPED_WAVE:      std::cout << "Damped Wave";      break;
        case SimulationMode::ABSORBING_BC:     std::cout << "Absorbing BC";     break;
    }
    std::cout << std::endl;
    std::cout << "  c=" << c << "  dt=" << time_step
              << "  T=" << end_time << std::endl;
    std::cout << "=======================================" << std::endl;

    // 1. Build mesh, system, matrices
    make_grid();
    setup_system();
    assemble_matrices();
    check_cfl_condition();

    // 2. Open energy log
    if (track_energy)
    {
        energy_log.open("energy_log.csv");
        energy_log << "step,time,kinetic_energy,potential_energy,total_energy\n";
    }

    // 3. Initial conditions
    std::cout << "  Setting initial conditions..." << std::endl;

    std::vector<Point<dim>> support_points(dof_handler.n_dofs());
    DoFTools::map_dofs_to_support_points(MappingQ1<dim>(), dof_handler, support_points);

    if (mode == SimulationMode::MMS_CONVERGENCE)
    {
        // MMS: u0 = sin(pi*x)*sin(pi*y), u1 = -pi*sin(pi*x)*sin(pi*y)
        InitialDisplacementMMS<dim> u0;
        InitialVelocityMMS<dim>     u1;
        VectorTools::interpolate(dof_handler, u0, solution_u);
        // u_old = u_current - dt * u1  (to encode initial velocity in Leapfrog)
        Vector<double> vel_vector(dof_handler.n_dofs());
        VectorTools::interpolate(dof_handler, u1, vel_vector);
        for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
            solution_u_old(i) = solution_u(i) - time_step * vel_vector(i);
        // Initial velocity estimate
        velocity_u = vel_vector;
    }
    else
    {
        // Default: Gaussian pulse (pebble in a pond)
        const double amplitude = 1.0;
        const Point<dim> center(0.5, 0.5);
        const double width = 0.05;

        for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
        {
            const double dist2 = center.distance_square(support_points[i]);
            solution_u(i) = amplitude * std::exp(-dist2 / (width * width));
        }
        // Zero initial velocity: u_old = u_current
        solution_u_old = solution_u;
        velocity_u     = 0.0;
    }

    constraints.distribute(solution_u);
    constraints.distribute(solution_u_old);

    // 4. Output step 0
    output_results(0);
    std::vector<std::pair<double, std::string>> times_and_names;
    times_and_names.push_back({0.0, "solution-0.vtu"});

    // 5. Time loop
    time        = 0.0;
    step_number = 0;

    while (time < end_time - 1e-12)
    {
        step_number++;
        time += time_step;

        // AMR every N steps (skip for MMS to keep uniform mesh)
        if (use_amr && mode != SimulationMode::MMS_CONVERGENCE &&
            step_number % amr_every_n_steps == 0)
        {
            refine_mesh();
        }

        solve_time_step();

        // Output and energy logging
        if (step_number % output_every_n_steps == 0)
        {
            std::cout << "Step " << step_number
                      << "  t=" << time << std::flush;

            if (track_energy)
            {
                const double Ek = compute_kinetic_energy();
                const double Ep = compute_potential_energy();
                const double Etot = Ek + Ep;
                std::cout << "  Ek=" << Ek
                          << "  Ep=" << Ep
                          << "  E_tot=" << Etot;
                energy_log << step_number << "," << time << ","
                           << Ek << "," << Ep << "," << Etot << "\n";
                energy_log.flush();
            }

            if (mode == SimulationMode::MMS_CONVERGENCE)
            {
                auto [L2, H1] = compute_errors(time);
                std::cout << "  L2=" << L2 << "  H1=" << H1;
            }

            std::cout << std::endl;

            output_results(step_number);
            times_and_names.push_back({time, "solution-" + std::to_string(step_number) + ".vtu"});

            // Update ParaView .pvd file
            std::ofstream pvd("solution.pvd");
            DataOutBase::write_pvd_record(pvd, times_and_names);
        }
    }

    if (track_energy) energy_log.close();
    std::cout << "Simulation finished. Total steps: " << step_number << std::endl;
}

// ============================================================
// run_convergence_study
//   Runs MMS on successively refined meshes to produce a
//   convergence table (L2 and H1 error vs. h).
//   Expected: O(h^2) for L2, O(h) for H1 with Q1 elements.
// ============================================================
template <int dim>
void WaveEquation<dim>::run_convergence_study()
{
    std::cout << "\n=== MMS Convergence Study ===" << std::endl;

    // Set MMS mode, disable AMR (we control the mesh manually)
    mode    = SimulationMode::MMS_CONVERGENCE;
    use_amr = false;

    // Run end_time = 1 cycle of the standing wave
    end_time = 1.0;

    const std::vector<unsigned int> refinement_levels = {3, 4, 5, 6}; //(livello 3 è una griglia 8x8 e livello 7  128x128)

    for (unsigned int ref : refinement_levels)
    {
        std::cout << "\n--- Refinement level " << ref << " ---" << std::endl;

        // Reset everything
        triangulation.clear();
        initial_refinement = ref;

        // Time step scaled to keep CFL: dt = 0.1 * h
        const double h = 1.0 / std::pow(2.0, ref);
        time_step = 0.4 * h / c;  // CFL number ~ 0.4  (safe for Q1 in 2D)
        // Round end_time to a multiple of dt to avoid drift
        const unsigned int n_steps = static_cast<unsigned int>(end_time / time_step);
        const double actual_end = n_steps * time_step;

        std::cout << "  h=" << h << "  dt=" << time_step
                  << "  n_steps=" << n_steps << std::endl;

        // Build
        make_grid();
        setup_system();
        assemble_matrices();
        check_cfl_condition();

        // Initial conditions (MMS)
        InitialDisplacementMMS<dim> u0;
        InitialVelocityMMS<dim>     u1;
        VectorTools::interpolate(dof_handler, u0, solution_u);
        Vector<double> vel_vector(dof_handler.n_dofs());
        VectorTools::interpolate(dof_handler, u1, vel_vector);
        for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
            solution_u_old(i) = solution_u(i) - time_step * vel_vector(i);
        velocity_u = vel_vector;
        constraints.distribute(solution_u);
        constraints.distribute(solution_u_old);

        // Time loop (no output to disk)
        time        = 0.0;
        step_number = 0;
        while (step_number < n_steps)
        {
            step_number++;
            time += time_step;
            solve_time_step();
        }

        // Compute errors at final time
        auto [L2, H1] = compute_errors(time);

        std::cout << "  L2 error = " << L2
                  << "  H1 error = " << H1 << std::endl;

        // Add to convergence table
        convergence_table.add_value("Level", ref);
        convergence_table.add_value("Cells", triangulation.n_active_cells());
        convergence_table.add_value("DoFs",  dof_handler.n_dofs());
        convergence_table.add_value("h",     h);
        convergence_table.add_value("L2",    L2);
        convergence_table.add_value("H1",    H1);
    }

    // Format and print convergence table
    convergence_table.set_precision("h",  4);
    convergence_table.set_precision("L2", 6);
    convergence_table.set_precision("H1", 6);
    convergence_table.set_scientific("h",  true);
    convergence_table.set_scientific("L2", true);
    convergence_table.set_scientific("H1", true);
    convergence_table.evaluate_convergence_rates(
        "L2", ConvergenceTable::reduction_rate_log2);
    convergence_table.evaluate_convergence_rates(
        "H1", ConvergenceTable::reduction_rate_log2);

    std::cout << "\n=== Convergence Table ===" << std::endl;
    convergence_table.write_text(std::cout);

    std::ofstream conv_file("convergence_table.txt");
    convergence_table.write_text(conv_file);
    std::cout << "\nConvergence table written to convergence_table.txt" << std::endl;
}

// Template instantiation
template class WaveEquation<2>;