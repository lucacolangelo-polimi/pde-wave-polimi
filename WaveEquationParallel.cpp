#include "WaveEquationParallel.hpp" //last version

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <chrono>
#include <iomanip> 

#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>

//Constructor
template <int dim>
WaveEquation<dim>::WaveEquation(MPI_Comm mpi_communicator)
    : mpi_comm(mpi_communicator)
    , n_mpi_procs(Utilities::MPI::n_mpi_processes(mpi_comm))
    , this_mpi_proc(Utilities::MPI::this_mpi_process(mpi_comm))
    , pcout(std::cout, this_mpi_proc == 0)   // print only on rank 0
    , computing_timer(mpi_comm,
                      pcout,
                      TimerOutput::never,
                      TimerOutput::wall_times)
    , triangulation(mpi_comm)                // distributed mesh p4est
    , fe_ptr(std::make_unique<FE_Q<dim>>(1))
    , dof_handler(triangulation)
    , time(0.0)
    , newmark_matrix_is_current(false)      
    , step_number(0)                         
{}


// wave_speed_at: c(x) — constant or eterogeneous (rifraction)
template <int dim>
double WaveEquation<dim>::wave_speed_at(const Point<dim> &p) const
{
    if (mode == SimulationMode::REFRACTION)
        return (p[1] > interface_y) ? c_fast : c_slow;
    return c;
}

// make_grid
template <int dim>
void WaveEquation<dim>::make_grid()
{
    TimerOutput::Scope t(computing_timer, "make_grid");

    pcout << "  Grid generation [0,1]^" << dim
          << "  (global refinements: " << initial_refinement << ")\n";

    GridGenerator::hyper_cube(triangulation, 0.0, 1.0);

    // refine_global on parallel::distributed::Triangulation
    //redistributes cells between processes
    triangulation.refine_global(initial_refinement);

    pcout << "  Active cells (global): "
          << triangulation.n_global_active_cells() << "\n";
}


// make_grid_with_obstacle  (DIFFRACTION)
template <int dim>
void WaveEquation<dim>::make_grid_with_obstacle()
{
    TimerOutput::Scope t(computing_timer, "make_grid_obstacle");
    pcout << "  grid generation with obstacle (diffraction)...\n";

    GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
    triangulation.refine_global(initial_refinement);

    //We mark material ONLY on local cells
    for (auto &cell : triangulation.active_cell_iterators())
    {
        if (!cell->is_locally_owned()) continue;

        const Point<dim> center = cell->center();
        // Wide wall in the center
        const bool in_wall_x = std::abs(center[0] - 0.5) < 0.05;
        // Narrow gap in the wall
        const bool in_gap_y  = (center[1] > 0.40 && center[1] < 0.60);

        if (in_wall_x && !in_gap_y)
            cell->set_material_id(1);
        else
            cell->set_material_id(0); 
    }

    pcout << "  Active cells (global): "
          << triangulation.n_global_active_cells() << "\n";
}

// setup_system
template <int dim>
void WaveEquation<dim>::setup_system()
{
    TimerOutput::Scope t(computing_timer, "setup_system");

    pcout << "  Setup system...\n";

    // Recreate FE with the chosen grade
    fe_ptr = std::make_unique<FE_Q<dim>>(fe_degree);
    dof_handler.distribute_dofs(*fe_ptr);

    // IndexSet: DoF owned e relevant 
    // locally_owned_dofs: those that this process "possesses"
    locally_owned_dofs    = dof_handler.locally_owned_dofs();
    // locally_relevant_dofs: owned + ghost (neighboring DoF needed for assembly)
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    pcout << "  DoF global: " << dof_handler.n_dofs()
          << "  (this process: " << locally_owned_dofs.n_elements() << ")\n";

    // CONSTRAINTS
    // for constraints:
    // 1. matrix_constraints: Contains only the constraints for the hanging nodes.
    // It is used to create the sparsity scheme and to assemble the
    // matrices, so as not to "condense" the matrix and keep it symmetric.
    // 2. constraints: Contains hanging nodes AND Dirichlet constraints. It is used
    // to apply the boundary values ​​to the vectors (solution, known terms).
    matrix_constraints.clear();
    matrix_constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, matrix_constraints);
    matrix_constraints.close();

    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    if (!use_absorbing_bc)
    {
        // Dirichlet u=0 on the boundary  (boundary_id=0)
        VectorTools::interpolate_boundary_values(
            dof_handler, 0,
            Functions::ZeroFunction<dim>(),
            constraints);
    }

    // Obstacle diffraction: binds the DoF of cells with material_id=1
    if (mode == SimulationMode::DIFFRACTION)
    {
        for (const auto &cell : dof_handler.active_cell_iterators())
        {
            if (cell->is_artificial()) continue;
            if (cell->material_id() != 1) continue;

            std::vector<types::global_dof_index> dof_ids(fe_ptr->dofs_per_cell);
            cell->get_dof_indices(dof_ids);
            
            for (auto idx : dof_ids)
            {
                if (locally_relevant_dofs.is_element(idx))
                {
                    if (!constraints.is_constrained(idx))
                    {
                        constraints.add_line(idx);
                        constraints.set_inhomogeneity(idx, 0.0);
                    }
                }
            }
        }
    }
    constraints.close();

    //Initialize the sponge vector
    nodal_sponge_profile.reinit(locally_owned_dofs, mpi_comm);
    
    if (mode == SimulationMode::DAMPED_WAVE)
    {
        // Map to get the physical coordinates of each DoF (node)
        std::map<types::global_dof_index, Point<dim>> support_points;
        MappingQ1<dim> mapping;
        DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

        for (const auto &pair : support_points)
        {
            const auto idx = pair.first;
            const auto &p  = pair.second;

            if (locally_owned_dofs.is_element(idx))
            {
                double sigma = 0.0;
                const double max_sigma = 50.0; // Maximum intensity of the damping
                const double sponge_thickness = 0.1; // Thickness of the sponge layer (from 0.9 to 1.0 and from 0 to 0.1)

                //Calculate the penetration of the node inside the frame
                double dist_x = std::max(sponge_thickness - p[0], p[0] - (1.0 - sponge_thickness));
                double dist_y = std::max(sponge_thickness - p[1], p[1] - (1.0 - sponge_thickness));
                
                double dist = std::max({0.0, dist_x, dist_y}); // Takes the greatest distance

                // If we are inside the sponge, sigma grows quadratically
                if (dist > 0.0)
                {
                    double normalized_dist = dist / sponge_thickness; //from 0 to 1
                    sigma = max_sigma * normalized_dist * normalized_dist;
                }
                
                nodal_sponge_profile(idx) = sigma;
            }
        }
    }

    // Sparsity pattern distributed 
    // DynamicSparsityPattern on locally_relevant_dofs, then We distribute it to remote processes with SparsityTools
    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, matrix_constraints, false);
    SparsityTools::distribute_sparsity_pattern(
        dsp,
        locally_owned_dofs,
        mpi_comm,
        locally_relevant_dofs);

    // Trilinos MATRICES
    // Each matrix is distributed: process i has the rows
    // matching locally_owned_dofs
    laplace_matrix.reinit(locally_owned_dofs,
                          locally_owned_dofs,
                          dsp, mpi_comm);
    if (use_absorbing_bc)
        boundary_mass_matrix.reinit(locally_owned_dofs,
                                    locally_owned_dofs,
                                    dsp, mpi_comm);
    if (time_scheme == TimeScheme::NEWMARK)
        system_matrix_newmark.reinit(locally_owned_dofs,
                                     locally_owned_dofs,
                                     dsp, mpi_comm);

    // ALLOCATION VECTORSD
    // "owned only" (no ghost): assembly (writings)
    mass_matrix_diagonal.reinit(locally_owned_dofs, mpi_comm);
    system_rhs.reinit(locally_owned_dofs, mpi_comm);
    owned_solution_u.reinit(locally_owned_dofs, mpi_comm);
    owned_solution_u_old.reinit(locally_owned_dofs, mpi_comm);
    owned_velocity_u.reinit(locally_owned_dofs, mpi_comm);
    owned_acceleration_u.reinit(locally_owned_dofs, mpi_comm);

    // "with ghost" (locally_relevant): for writing during assembly
    solution_u.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    solution_u_old.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    solution_u_new.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    velocity_u.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    acceleration_u.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
}

// ============================================================
// assemble_matrices
//
// Each process will only be on its "locally owned cells".
// Contributions are accumulated locally and then Trilinos
// performs communication with compresses(VectorOperation::add).
template <int dim>
void WaveEquation<dim>::assemble_matrices()
{
    TimerOutput::Scope t(computing_timer, "assemble_matrices");

    pcout << "  Assembly matrices (MPI, " << n_mpi_procs << " processes)...\n";

    // Reset
    laplace_matrix       = 0.0;
    mass_matrix_diagonal = 0.0;
    if (use_absorbing_bc)
        boundary_mass_matrix = 0.0;

    QGauss<dim>        q_stiff(fe_ptr->degree + 1);
    QGaussLobatto<dim> q_mass (fe_ptr->degree + 1);
    QGauss<dim - 1>    q_face (fe_ptr->degree + 1);

    FEValues<dim> fev_stiff(*fe_ptr, q_stiff,
        update_gradients | update_JxW_values | update_quadrature_points);
    FEValues<dim> fev_mass(*fe_ptr, q_mass,
        update_values | update_JxW_values);
    FEFaceValues<dim> fev_face(*fe_ptr, q_face,
        update_values | update_JxW_values);

    const unsigned int dpc = fe_ptr->dofs_per_cell;

    FullMatrix<double> cell_K(dpc, dpc);
    FullMatrix<double> cell_B(dpc, dpc);
    Vector<double>     cell_M(dpc);
    std::vector<types::global_dof_index> local_idx(dpc);

    // Iterates only on the locally owned cells 
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
        if (!cell->is_locally_owned()) continue;

        fev_stiff.reinit(cell);
        fev_mass.reinit(cell);

        cell_K = 0.0;
        cell_M = 0.0;

        // Obstacle cells: no physical contribution
        if (cell->material_id() == 1)
        {
            cell->get_dof_indices(local_idx);
            matrix_constraints.distribute_local_to_global(cell_K, local_idx, laplace_matrix); // FIX 1
            continue;
        }

        // matrix K with c²(x) variable
        for (unsigned int q = 0; q < q_stiff.size(); ++q)
        {
            const double cq  = wave_speed_at(fev_stiff.quadrature_point(q));
            const double JxW = fev_stiff.JxW(q);
            for (unsigned int i = 0; i < dpc; ++i)
                for (unsigned int j = 0; j < dpc; ++j)
                    cell_K(i, j) += cq * cq
                                  * fev_stiff.shape_grad(i, q)
                                  * fev_stiff.shape_grad(j, q)
                                  * JxW;
        }

        // Mass lumped (Gauss-Lobatto)
        for (unsigned int q = 0; q < q_mass.size(); ++q)
        {
            const double JxW = fev_mass.JxW(q);
            for (unsigned int i = 0; i < dpc; ++i)
            {
                const double phi = fev_mass.shape_value(i, q);
                cell_M(i) += phi * phi * JxW;
            }
        }

        cell->get_dof_indices(local_idx);

        // Distribuisce usando matrix_constraints (solo hanging nodes, NO Dirichlet)
        matrix_constraints.distribute_local_to_global(cell_K, local_idx, laplace_matrix); // FIX 2
        matrix_constraints.distribute_local_to_global(cell_M, local_idx, mass_matrix_diagonal); // FIX 3

       // Matrix for ABC (absorbing boundary conditions)
        if (use_absorbing_bc)
        {
            for (const auto &face : cell->face_iterators())
            {
                if (!face->at_boundary()) continue;
                fev_face.reinit(cell, face);
                cell_B = 0.0;
                for (unsigned int q = 0; q < q_face.size(); ++q)
                {
                    // Calcola c sul bordo per la condizione di Sommerfeld
                    const double cq = wave_speed_at(fev_face.quadrature_point(q));
                    const double JxW = fev_face.JxW(q);
                    for (unsigned int i = 0; i < dpc; ++i)
                        for (unsigned int j = 0; j < dpc; ++j)
                            // Aggiungiamo 'cq' all'integrale!
                            cell_B(i, j) += cq * fev_face.shape_value(i, q)
                                          * fev_face.shape_value(j, q)
                                          * JxW;
                }
                matrix_constraints.distribute_local_to_global(cell_B, local_idx, boundary_mass_matrix);
            }
        }
    }

    // MPI communication
    laplace_matrix.compress(VectorOperation::add);
    mass_matrix_diagonal.compress(VectorOperation::add);
    if (use_absorbing_bc)
        boundary_mass_matrix.compress(VectorOperation::add);

    // Force mass to 1 for constrained (Dirichlet) DoFs for algebraic safety
    for (const auto idx : locally_owned_dofs)
    {
        if (constraints.is_constrained(idx))
            mass_matrix_diagonal(idx) = 1.0;
    }
    mass_matrix_diagonal.compress(VectorOperation::insert);

    // mass positive check
    const double local_min = mass_matrix_diagonal.min();
    const double global_min = Utilities::MPI::min(local_min, mpi_comm);
    AssertThrow(global_min > 0.0,
        ExcMessage("Lumped mass matrix: entry diagonal non positive!"));

    pcout << "  Assembly completed. min(M_diag)=" << global_min << "\n";

    // FOR NEWMARK: build system matrix A = M + β·dt²·K
    if (time_scheme == TimeScheme::NEWMARK)
    {
        newmark_matrix_is_current = false; 
        build_newmark_system_matrix();
    }
}

// ============================================================
// assemble_rhs
//
// Calculate: RHS = −K·u + F(t) − c·B·v − d·M·v
// Note: laplace_matrix.vmult() and boundary_mass_matrix.vmult()
// use MPI_Allreduce internally to add up contributions
// of ghost processes. The result is already distributed.
template <int dim>
void WaveEquation<dim>::assemble_rhs(double t)
{
    system_rhs = 0.0;

    //We synchronize the owned vector by eliminating the ghosts before the Trilinos vmult
    owned_solution_u = solution_u; 
    
    // Now we use owned_solution_u (owned-only), NOT solution_u
    laplace_matrix.vmult(system_rhs, owned_solution_u);
    system_rhs *= -1.0;

    // Forcing term f(x,t) for MMS
    if (mode == SimulationMode::MMS_CONVERGENCE)
    {
        QGauss<dim> q(fe_ptr->degree + 1);
        FEValues<dim> fev(*fe_ptr, q, update_values | update_quadrature_points | update_JxW_values);
        
        std::unique_ptr<Function<dim>> forcing_ptr;
        if (use_decay_mms) {
            auto p = std::make_unique<ForcingTermMMS_Decay<dim>>(c);
            p->set_time(t);
            forcing_ptr = std::move(p);
        } else {
            auto p = std::make_unique<ForcingTermMMS<dim>>(c);
            p->set_time(t);
            forcing_ptr = std::move(p);
        }

        const unsigned int dpc = fe_ptr->dofs_per_cell;
        Vector<double> cell_rhs(dpc);
        std::vector<types::global_dof_index> local_idx(dpc);

        for (const auto &cell : dof_handler.active_cell_iterators())
        {
            if (!cell->is_locally_owned()) continue;
            fev.reinit(cell);
            cell_rhs = 0.0;
            for (unsigned int q_pt = 0; q_pt < q.size(); ++q_pt)
            {
                const double fval = forcing_ptr->value(fev.quadrature_point(q_pt));
                for (unsigned int i = 0; i < dpc; ++i)
                    cell_rhs(i) += fev.shape_value(i, q_pt) * fval * fev.JxW(q_pt);
            }
            cell->get_dof_indices(local_idx);
            matrix_constraints.distribute_local_to_global(cell_rhs, local_idx, system_rhs);
        }
        system_rhs.compress(VectorOperation::add);
    }

    // ABC: −c·B·v
    if (use_absorbing_bc)
    {
        // Let's convert velocity_u (ghosted) to owned before vmult
        owned_velocity_u = velocity_u; 
        
        TrilinosVector abc_contrib(locally_owned_dofs, mpi_comm);
        boundary_mass_matrix.vmult(abc_contrib, owned_velocity_u); // Uses owned_velocity_u
        system_rhs.add(-c, abc_contrib);
    }

    // Damping: −d·M·v
    if (damping > 0.0)
    {
        for (const auto idx : locally_owned_dofs)
            system_rhs(idx) -= damping * mass_matrix_diagonal(idx) * owned_velocity_u(idx);
        system_rhs.compress(VectorOperation::add);
    }
}

// ============================================================
// solve_time_step  —  Leapfrog esplicito
//
// u^{n+1} = 2·u^n − u^{n-1} + dt²·M^{-1}·RHS
// With lumpy mass M^{-1} is local scalar division:
// no MPI communication for the solve.
// Communication takes place only in assemble_rhs (vmult).
// ============================================================
template <int dim>
void WaveEquation<dim>::solve_time_step()
{
    TimerOutput::Scope t(computing_timer, "solve_leapfrog");

    assemble_rhs(time - time_step);

    // Update ghost values for the vectors we read
    solution_u.update_ghost_values();
    solution_u_old.update_ghost_values();

    // 1. Calculate the damping force at the boundary: C * (u^n - u^{n-1}) / dt
    TrilinosVector damping_forces(locally_owned_dofs, mpi_comm);
    if (use_absorbing_bc)
    {
        TrilinosVector vel_approx(locally_owned_dofs, mpi_comm);
        for (const auto idx : locally_owned_dofs)
        {
            vel_approx(idx) = (solution_u(idx) - solution_u_old(idx)) / time_step;
        }
        
        //Trilinos automatically communicates MPI edges during this multiplication
        boundary_mass_matrix.vmult(damping_forces, vel_approx);
    }

    // a_i = rhs_i / M_ii  — local operation (only owned dofs)
    for (const auto idx : locally_owned_dofs)
    {
        if (constraints.is_constrained(idx))
        {
            owned_acceleration_u(idx) = 0.0;
            continue;
        }

        const double m_ii = mass_matrix_diagonal(idx);

        AssertThrow(std::isfinite(m_ii) && m_ii > 1e-30,
                    ExcMessage("Invalid mass_matrix_diagonal entry"));

        // 2. Subtract the damping force from the known spring term
        double current_rhs = system_rhs(idx);
        if (use_absorbing_bc)
            current_rhs -= damping_forces(idx);

        owned_acceleration_u(idx) = current_rhs / m_ii;
    }

    // u^{n+1} = 2·u^n − u^{n-1} + dt²·a
    const double dt2 = time_step * time_step;
    for (const auto idx : locally_owned_dofs)
    {
        double sigma = 0.0;
        if (mode == SimulationMode::DAMPED_WAVE)
            sigma = nodal_sponge_profile(idx);

        // Damping factors derived from the centered discretization of the first derivative u_t
        const double damp_plus  = 1.0 + 0.5 * sigma * time_step;
        const double damp_minus = 1.0 - 0.5 * sigma * time_step;

        owned_solution_u(idx) = ( 2.0 * solution_u(idx)
                                - damp_minus * solution_u_old(idx)
                                + dt2 * owned_acceleration_u(idx) ) / damp_plus;
    }

    // Apply Dirichlet constraints: local operation
    if (!use_absorbing_bc)
        constraints.distribute(owned_solution_u);

    // Centered speed: v^n = (u^{n+1} − u^{n-1}) / (2·dt)
    for (const auto idx : locally_owned_dofs)
        owned_velocity_u(idx) = (owned_solution_u(idx) - solution_u_old(idx))
                              / (2.0 * time_step);

    // Shift: advances one step 
    // Copy owned → ghost (update_ghost_values propagates to neighbors)
    owned_solution_u_old = owned_solution_u;

    solution_u_old = solution_u;        // u^{n-1} <- u^n
    solution_u     = owned_solution_u;  // u^n <- u^{n+1}
    velocity_u     = owned_velocity_u;

    //Makes new ghost values available for the next step
    solution_u.update_ghost_values();
    solution_u_old.update_ghost_values();
    velocity_u.update_ghost_values();
}

// ------------------------------------------------------------
// solve_time_step_newmark OTTIMIZZATO
// uses AMG rather than SSOR (convergence in O(1) iterations)
template <int dim>
void WaveEquation<dim>::solve_time_step_newmark()
{
    TimerOutput::Scope timer(computing_timer, "solve_newmark");

    AssertThrow(newmark_matrix_is_current,
        ExcMessage("build_newmark_system_matrix() not called!"));

    const double dt  = time_step;
    const double b   = newmark_beta;
    const double gam = newmark_gamma;

    // 1 Update ghosts before starting
    solution_u.update_ghost_values();
    velocity_u.update_ghost_values();
    acceleration_u.update_ghost_values();

    // 2Calculate the predictors
    TrilinosVector u_pred(locally_owned_dofs, mpi_comm);
    TrilinosVector v_pred(locally_owned_dofs, mpi_comm);
    for (const auto idx : locally_owned_dofs)
    {
        u_pred(idx) = solution_u(idx)
                    + dt * velocity_u(idx)
                    + dt * dt * (0.5 - b) * acceleration_u(idx);
        v_pred(idx) = velocity_u(idx)
                    + dt * (1.0 - gam) * acceleration_u(idx);
    }
    u_pred.compress(VectorOperation::insert);

    TrilinosVector u_pred_ghosted(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    u_pred_ghosted = u_pred;
    u_pred_ghosted.update_ghost_values();

    // 3 RHS = F(t_{n+1}) - K * u_pred
    TrilinosVector rhs_newmark(locally_owned_dofs, mpi_comm);
    
    //--> Elastic term: -K * u_pred
    laplace_matrix.vmult(rhs_newmark, u_pred_ghosted);
    rhs_newmark *= -1.0;

    // --> Exact forcing term at t_{n+1}
    if (mode == SimulationMode::MMS_CONVERGENCE)
    {
        QGauss<dim> q(fe_ptr->degree + 1);
        FEValues<dim> fev(*fe_ptr, q, update_values | update_quadrature_points | update_JxW_values);
        
        ForcingTermMMS<dim> forcing(c); 
        forcing.set_time(time);      

        Vector<double> cell_rhs(fe_ptr->dofs_per_cell);
        std::vector<types::global_dof_index> local_idx(fe_ptr->dofs_per_cell);
        TrilinosVector F_ext(locally_owned_dofs, mpi_comm);

        for (const auto &cell : dof_handler.active_cell_iterators())
        {
            if (!cell->is_locally_owned()) continue;
            fev.reinit(cell);
            cell_rhs = 0.0;
            for (unsigned int q_pt = 0; q_pt < fev.get_quadrature().size(); ++q_pt)
            {
                const double fval = forcing.value(fev.quadrature_point(q_pt));
                for (unsigned int i = 0; i < fe_ptr->dofs_per_cell; ++i)
                    cell_rhs(i) += fev.shape_value(i, q_pt) * fval * fev.JxW(q_pt);
            }
            cell->get_dof_indices(local_idx);
            matrix_constraints.distribute_local_to_global(cell_rhs, local_idx, F_ext);
        }
        F_ext.compress(VectorOperation::add);
        rhs_newmark.add(1.0, F_ext); // Adds F(t_{n+1}) to the right side
    }

    // Reset the degrees of freedom on the Dirichlet edges
    constraints.set_zero(rhs_newmark);

    // 4. Solving the system A * a_new = rhs_newmark
    TrilinosVector a_new(locally_owned_dofs, mpi_comm);
    SolverControl solver_control(500, 1e-10 * rhs_newmark.l2_norm() + 1e-30);
    TrilinosWrappers::SolverCG cg_solver(solver_control);

    cg_solver.solve(system_matrix_newmark, a_new, rhs_newmark, newmark_preconditioner);
    constraints.distribute(a_new);

    if (step_number % output_every_n_steps == 0)
        pcout << "    Newmark CG: " << solver_control.last_step() << " iterazioni\n";

    // 5. Final corrections
    for (const auto idx : locally_owned_dofs)
    {
        owned_solution_u(idx)     = u_pred(idx) + b * dt * dt * a_new(idx);
        owned_velocity_u(idx)     = v_pred(idx) + gam * dt * a_new(idx);
        owned_acceleration_u(idx) = a_new(idx);
    }
    constraints.distribute(owned_solution_u);

    solution_u_old = solution_u;
    solution_u     = owned_solution_u;
    velocity_u     = owned_velocity_u;
    acceleration_u = owned_acceleration_u;

    solution_u.update_ghost_values();
    solution_u_old.update_ghost_values();
    velocity_u.update_ghost_values();
    acceleration_u.update_ghost_values();
}

// ============================================================
// compute_kinetic_energy
// E_k = 0.5 · v^T · M · v
// With diagonal M and v distributed, each process calculates
// the local sum, then MPI_Allreduce sum globally.
template <int dim>
double WaveEquation<dim>::compute_kinetic_energy() //const
{
    double local_ek = 0.0;
    for (const auto idx : locally_owned_dofs)
    {
        const double vi = velocity_u(idx);
        local_ek += 0.5 * mass_matrix_diagonal(idx) * vi * vi;
    }
    return Utilities::MPI::sum(local_ek, mpi_comm);
}

// ============================================================
// compute_potential_energy
//
// E_p = 0.5 · u^T · K · u
// laplace_matrix.vmult() already distrubuted
template <int dim>
double WaveEquation<dim>::compute_potential_energy()
{
    TrilinosVector Ku(locally_owned_dofs, mpi_comm);
    
    // Update ghosts before multiplying
    solution_u.update_ghost_values();
    laplace_matrix.vmult(Ku, solution_u); // Uses solution_u

    double local_ep = 0.0;
    for (const auto idx : locally_owned_dofs)
        local_ep += 0.5 * solution_u(idx) * Ku(idx); // Uses solution_u

    return Utilities::MPI::sum(local_ep, mpi_comm);
}

// ============================================================
// compute_Linfty_norm
// Calculate max|u_i| on locally owned DoF, then MPI_Allreduce.
// Useful for detecting numerical instability (violated CFL).
template <int dim>
double WaveEquation<dim>::compute_Linfty_norm()  //const
{
    double local_max = 0.0;
    for (const auto idx : locally_owned_dofs)
        local_max = std::max(local_max, std::abs(solution_u(idx)));

    return Utilities::MPI::max(local_max, mpi_comm);
}

// ============================================================
// compute_errors (MMS)
// VectorTools::integrate_difference already parallel:
// calculates the local error and then MPI_Allreduce
template <int dim>
std::pair<double,double> WaveEquation<dim>::compute_errors(double t) const
{
    std::unique_ptr<Function<dim>> exact_ptr;
    if (use_decay_mms) {
        auto p = std::make_unique<ExactSolutionMMS_Decay<dim>>(c);
        p->set_time(t);
        exact_ptr = std::move(p);
    } else {
        auto p = std::make_unique<ExactSolutionMMS<dim>>(c);
        p->set_time(t);
        exact_ptr = std::move(p);
    }

    Vector<float> diff(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler, solution_u, *exact_ptr, diff, QGauss<dim>(fe_ptr->degree + 2), VectorTools::L2_norm);
    const double L2 = VectorTools::compute_global_error(triangulation, diff, VectorTools::L2_norm);

    VectorTools::integrate_difference(dof_handler, solution_u, *exact_ptr, diff, QGauss<dim>(fe_ptr->degree + 2), VectorTools::H1_seminorm);
    const double H1 = VectorTools::compute_global_error(triangulation, diff, VectorTools::H1_seminorm);

    return {L2, H1};
}

// ============================================================
// check_cfl_condition
// h_min is calculated locally, then MPI_Allreduce takes the global min.
template <int dim>
void WaveEquation<dim>::check_cfl_condition() const
{
    if (time_scheme == TimeScheme::NEWMARK)
    {
        pcout << "  Newmark-β scheme:unconditionally stable.\n";
        return;
    }

    double local_h_min = std::numeric_limits<double>::max();
    for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
            local_h_min = std::min(local_h_min, cell->minimum_vertex_distance());

    // MPI_Allreduce: takes the global minimum
    const double h_min  = Utilities::MPI::min(local_h_min, mpi_comm);
    const double c_max  = (mode == SimulationMode::REFRACTION)
                           ? std::max(c_fast, c_slow) : c;
    const double cfl_lim = h_min / (c_max * std::sqrt(static_cast<double>(dim)));
    const double cfl_num = time_step / cfl_lim;

    pcout << "  CFL: dt=" << time_step
          << "  h_min=" << h_min
          << "  c_max=" << c_max
          << "  CFL=" << cfl_num;
    if (cfl_num > 1.0)
        pcout << "  *** WARNING: CFL VIOLATED! ***";
    else
        pcout << "  [OK]";
    pcout << "\n";
}

//============================================================
// output_results — parallel VTU output
//
// Strategy deal.II MPI:
// - Each process writes its own solution-NNNN-RRRR.vtu file  where RRRR = MPI rank
// - Process 0 writes the master file solution-NNNN.pvtu  which contains the list of all .vtu files
// - ParaView loads the .pvtu and rebuilds automatically the complete solution
template <int dim>
void WaveEquation<dim>::output_results(unsigned int step)
{
    TimerOutput::Scope t(computing_timer, "output");

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);

    // "with ghost" vectors already sinchronized
    data_out.add_data_vector(solution_u, "displacement");
    data_out.add_data_vector(velocity_u, "velocity");

    //ADDED WALL
    Vector<float> material_ids(triangulation.n_active_cells());
    int i = 0;
    for (auto &cell : triangulation.active_cell_iterators())
    {
        material_ids(i) = cell->material_id();
        ++i;
    }
    data_out.add_data_vector(material_ids, "Material_ID");

    // SPONGE LAYER (DAMPING)
    Vector<float> damping_profile(triangulation.n_active_cells());
    int j = 0;
    for (auto &cell : triangulation.active_cell_iterators())
    {
        const Point<dim> center = cell->center();
        if (center[0] > 0.9 || center[0] < 0.1 || center[1] > 0.9 || center[1] < 0.1)
            damping_profile(j) = 1.0; // we're in the sponge layer
        else
            damping_profile(j) = 0.0; // we're in the interior
            
        ++j;
    }
    data_out.add_data_vector(damping_profile, "Damping_Coefficient");

    // Adds the MPI rank as a field (useful for verifying partitioning)
    Vector<float> proc_id(triangulation.n_active_cells());
    proc_id = static_cast<float>(this_mpi_proc);
    data_out.add_data_vector(proc_id, "mpi_rank");

    data_out.build_patches();

    // file name base: solution-NNNN where NNNN is the time step number
    const std::string base = "solution-" + Utilities::int_to_string(step, 4);

    // for every process writes its own .vtu file with the local solution
    std::ofstream local_out(base + "-" +
        Utilities::int_to_string(this_mpi_proc, 4) + ".vtu");
    data_out.write_vtu(local_out);

    // Only rank 0 writes the .pvtu master file
    if (this_mpi_proc == 0)
    {
        std::vector<std::string> file_list;
        for (unsigned int r = 0; r < n_mpi_procs; ++r)
            file_list.push_back(base + "-" +
                Utilities::int_to_string(r, 4) + ".vtu");

        std::ofstream pvtu_out(base + ".pvtu");
        data_out.write_pvtu_record(pvtu_out, file_list);
    }
}

// ============================================================
// write_energy_report
// Call at the end of run() — print global statistics
// on stdout (only rank 0) and writes energy_report.txt
template <int dim>
void WaveEquation<dim>::write_energy_report() const
{
    if (this_mpi_proc != 0) return;

    std::ofstream report("energy_report.txt");
    report << "======================================\n";
    report << "  Energy Conservation Report\n";
    report << "  Time scheme: "
           << (time_scheme == TimeScheme::LEAPFROG ? "Leapfrog" : "Newmark-beta")
           << "\n";
    if (time_scheme == TimeScheme::NEWMARK)
        report << "  beta=" << newmark_beta << "  gamma=" << newmark_gamma << "\n";
    report << "  dt=" << time_step
           << "  T=" << end_time
           << "  fe_degree=Q" << fe_degree << "\n";
    report << "  MPI processes: " << n_mpi_procs << "\n";
    report << "======================================\n";
    report << "  E0 (initial total energy): " << energy_initial << "\n";
    report << "\n  See energy_log.csv for full time series.\n";
    report << "\n  Expected drift:\n";
    report << "    Leapfrog      : O(dt^2) oscillation, zero mean drift\n";
    report << "    Newmark b=0.25: near-zero drift (symplectic-like)\n";
    report << "    Newmark b>0.25: monotone decay (numerical dissipation)\n";
    report << "======================================\n";

    pcout << "\n  [Energy report salvato in energy_report.txt]\n";
}


//========================================================
// refine_mesh — parallel AMR with p4est
// parallel::distributed::GridRefinement distributes flags between processes and p4est manages load rebalancing.
// parallel::distributed::SolutionTransfer transfers the solution on the new distributed mesh.
template <int dim>
void WaveEquation<dim>::refine_mesh()
{
    TimerOutput::Scope t(computing_timer, "refine_mesh");

    pcout << "  AMR parallel...\n";

    // Estimate local error (only locally owned cells)
    Vector<float> err(triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate(
        dof_handler,
        QGauss<dim - 1>(fe_ptr->degree + 1),
        std::map<types::boundary_id, const Function<dim> *>(),
        solution_u, err);

    // parallel::distributed::GridRefinement: coordinates the flags between processes
    parallel::distributed::GridRefinement
        ::refine_and_coarsen_fixed_fraction(
            triangulation, err,
            amr_refine_fraction, amr_coarsen_fraction);

    // Limit the maximum level & PROTECT THE WALL
    for (auto &cell : triangulation.active_cell_iterators())
    {
        if (cell->is_locally_owned())
        {
            // Limite standard per il livello massimo
            if (cell->level() >= static_cast<int>(max_refinement_level))
            {
                cell->clear_refine_flag();
            }
            if (mode == SimulationMode::DIFFRACTION && cell->material_id() == 1)
            {
                cell->clear_coarsen_flag(); 
            }
        }
    }

    triangulation.prepare_coarsening_and_refinement();

    // SolutionTransfer distributed
    parallel::distributed::SolutionTransfer<dim, TrilinosVector>
        st_u(dof_handler), st_u_old(dof_handler), st_v(dof_handler);

    solution_u.update_ghost_values();
    solution_u_old.update_ghost_values();
    velocity_u.update_ghost_values();

    st_u.prepare_for_coarsening_and_refinement(solution_u);
    st_u_old.prepare_for_coarsening_and_refinement(solution_u_old);
    st_v.prepare_for_coarsening_and_refinement(velocity_u);


    // PERFORM GEOMETRIC REPARTITIONING
    triangulation.execute_coarsening_and_refinement();

    if (mode == SimulationMode::DIFFRACTION)
    {
        for (auto &cell : triangulation.active_cell_iterators())
        {
            if (!cell->is_locally_owned()) continue;

            const Point<dim> center = cell->center();
            const bool in_wall_x = std::abs(center[0] - 0.5) < 0.05;
            const bool in_gap_y  = (center[1] > 0.40 && center[1] < 0.60);

            if (in_wall_x && !in_gap_y)
                cell->set_material_id(1);
            else
                cell->set_material_id(0);
        }
    }

    // Redistributes DoF on new mesh
    dof_handler.distribute_dofs(*fe_ptr);

    pcout << "  AMR: DoF global=" << dof_handler.n_dofs()
          << "  cells=" << triangulation.n_global_active_cells() << "\n";

    setup_system(); 

    // Interpolates the solutions on the new mesh
    TrilinosVector interp_u(locally_owned_dofs, mpi_comm);
    TrilinosVector interp_u_old(locally_owned_dofs, mpi_comm);
    TrilinosVector interp_v(locally_owned_dofs, mpi_comm);

    st_u.interpolate(interp_u);
    st_u_old.interpolate(interp_u_old);
    st_v.interpolate(interp_v);

    constraints.distribute(interp_u);
    constraints.distribute(interp_u_old);

    solution_u     = interp_u;
    solution_u_old = interp_u_old;
    velocity_u     = interp_v;

    solution_u.update_ghost_values();
    solution_u_old.update_ghost_values();
    velocity_u.update_ghost_values();

    // Riassembly matrices on a new mesh
    assemble_matrices();
    check_cfl_condition();

    if (time_scheme == TimeScheme::NEWMARK)
    {
        assemble_rhs(time);

        for (const auto idx : locally_owned_dofs)
        {
            if (constraints.is_constrained(idx))
                owned_acceleration_u(idx) = 0.0;
            else
                owned_acceleration_u(idx) =
                    system_rhs(idx) / mass_matrix_diagonal(idx);
        }

        constraints.distribute(owned_acceleration_u);
        acceleration_u = owned_acceleration_u;
        acceleration_u.update_ghost_values();
    }

    if (time_scheme == TimeScheme::LEAPFROG)
    {
        double local_h = std::numeric_limits<double>::max();
        for (const auto &cell : dof_handler.active_cell_iterators())
        {
            if (cell->is_locally_owned())
                local_h = std::min(local_h, cell->minimum_vertex_distance());
        }
        const double h_min_new = Utilities::MPI::min(local_h, mpi_comm);
        const double c_max     = (mode == SimulationMode::REFRACTION)
                                ? std::max(c_fast, c_slow) : c;
        const double dt_safe = 0.45 * h_min_new / (c_max * std::sqrt((double)dim));
        
        if (time_step > dt_safe)
        {
            time_step = dt_safe;
            newmark_matrix_is_current = false;
            pcout << "  [AMR] dt riducted at " << time_step
                << " in order to satisfy CFL\n";
        }
    }
}

// ============================================================
// run  —  main time loop
template <int dim>
void WaveEquation<dim>::run()
{
    // I 
    energy_initialized = false;
    energy_initial     = 0.0;
    pcout << "\n=======================================\n"
          << "  WaveEquation<" << dim << ">  MPI  "
          << n_mpi_procs << " processi\n  Modes: ";
    switch (mode)
    {
        case SimulationMode::PEBBLE_IN_POND:  pcout << "Pebble in Pond\n";  break;
        case SimulationMode::MMS_CONVERGENCE: pcout << "MMS\n";             break;
        case SimulationMode::DAMPED_WAVE:     pcout << "Damped\n";        break;
        case SimulationMode::ABSORBING_BC:    pcout << "ABC Sommerfeld\n";  break;
        case SimulationMode::INTERFERENCE:    pcout << "Interference\n";    break;
        case SimulationMode::REFRACTION:      pcout << "Refraction\n";      break;
        case SimulationMode::DIFFRACTION:     pcout << "Diffraction\n";     break;
    }
    pcout << "  Scheme: "
          << (time_scheme == TimeScheme::LEAPFROG ? "Leapfrog" : "Newmark-β")
          << "  Q" << fe_degree
          << "  c=" << c << "  dt=" << time_step << "  T=" << end_time << "\n"
          << "=======================================\n";

    //grid and system setup
    if (mode == SimulationMode::DIFFRACTION)
        make_grid_with_obstacle();
    else
        make_grid();

    setup_system();
    assemble_matrices();
    check_cfl_condition();

    //Energy log (only rank 0)
    if (track_energy && this_mpi_proc == 0)
    {
        energy_log.open("energy_log.csv");
        energy_log << "step,time,kinetic,potential,total,drift_rel,Linfty\n";
    }

    // Initial conditions
    pcout << "  Initial conditions...\n";

    if (mode == SimulationMode::MMS_CONVERGENCE)
    {
        // IC from exact solution: u(x,0) = u0(x)  v(x,0) = u1(x)
        InitialDisplacementMMS<dim> u0;
        InitialVelocityMMS<dim>     u1;

        VectorTools::interpolate(dof_handler, u0, owned_solution_u);
        constraints.distribute(owned_solution_u);
        solution_u = owned_solution_u;

        VectorTools::interpolate(dof_handler, u1, owned_velocity_u);
        // u_old = u - dt*v 
        for (const auto idx : locally_owned_dofs)
            owned_solution_u_old(idx) = owned_solution_u(idx)
                                      - time_step * owned_velocity_u(idx);
        solution_u_old = owned_solution_u_old;
        velocity_u     = owned_velocity_u;
    }
    else if (mode == SimulationMode::INTERFERENCE)
    {
        // two gaussian wave packets centered at s1 and s2 (interference pattern)
        const double amp   = 1.0;
        const double width = 0.15;                                                              // 0.15 is betetr for 3d, 0.05 is suitable for 2d (more localized)
        const Point<dim> s1(0.3, 0.5), s2(0.7, 0.5);

        // We use interpolated with a lambda function via FunctionFromFunctionObjects (easier: we scroll through local support points)
        std::vector<Point<dim>> sp(dof_handler.n_dofs());
        MappingQ1<dim> mapping;
        DoFTools::map_dofs_to_support_points(mapping, dof_handler, sp);

        for (const auto idx : locally_owned_dofs)
        {
            const double d1 = s1.distance_square(sp[idx]);
            const double d2 = s2.distance_square(sp[idx]);
            owned_solution_u(idx) = amp * std::exp(-d1 / (width*width))
                                  + amp * std::exp(-d2 / (width*width));
        }
        constraints.distribute(owned_solution_u);
        solution_u     = owned_solution_u;
        solution_u_old = owned_solution_u;
        owned_velocity_u = 0.0;
        velocity_u       = owned_velocity_u;
    }
    else
    {
        // single gaussian wave packet centered at s1 (pebble in pond)
        std::vector<Point<dim>> sp(dof_handler.n_dofs());
        MappingQ1<dim> mapping;
        DoFTools::map_dofs_to_support_points(mapping, dof_handler, sp);

        const double amp   = 1.0;
        const double width = 0.05;
        const Point<dim> center = (dim == 2)
            ? Point<dim>(0.25, 0.5)
            : Point<dim>(0.25, 0.5, 0.5);

        for (const auto idx : locally_owned_dofs)
        {
            const double d2 = center.distance_square(sp[idx]);
            owned_solution_u(idx) = amp * std::exp(-d2 / (width*width));
        }
        constraints.distribute(owned_solution_u);
        solution_u     = owned_solution_u;
        solution_u_old = owned_solution_u;
        owned_velocity_u = 0.0;
        velocity_u       = owned_velocity_u;
    }

    solution_u.update_ghost_values();
    solution_u_old.update_ghost_values();
    velocity_u.update_ghost_values();

    if (time_scheme == TimeScheme::NEWMARK)
        {
            assemble_rhs(0.0);

            for (const auto idx : locally_owned_dofs)
            {
                if (constraints.is_constrained(idx))
                {
                    owned_acceleration_u(idx) = 0.0;
                    continue;
                }

                const double m_ii = mass_matrix_diagonal(idx);

                AssertThrow(std::isfinite(m_ii) && m_ii > 1e-30,
                            ExcMessage("Invalid mass_matrix_diagonal entry in initial acceleration"));

                owned_acceleration_u(idx) = system_rhs(idx) / m_ii;
            }

            constraints.distribute(owned_acceleration_u);
            acceleration_u = owned_acceleration_u;
            acceleration_u.update_ghost_values();
        }

    output_results(0);

    double linf_initial = compute_Linfty_norm();
    const double blowup_threshold = std::max(100.0, 100.0 * linf_initial);
    // .pvtu list for .pvd master (only rank 0)
    std::vector<std::pair<double, std::string>> pvd_list;
    if (this_mpi_proc == 0)
        pvd_list.push_back({0.0, "solution-0000.pvtu"});

    // time loop
    time        = 0.0;
    step_number = 0;

    while (time < end_time - 1e-12)
    {
        step_number++;
        time += time_step;

        // AMR parallel (not for MMS)
        if (use_amr && mode != SimulationMode::MMS_CONVERGENCE
                    && step_number % amr_every_n_steps == 0)
            refine_mesh();

        // Advance one time step
        if (time_scheme == TimeScheme::LEAPFROG)
            solve_time_step();
        else
            solve_time_step_newmark();

        // Output 
        if (step_number % output_every_n_steps == 0)
        {
            pcout << "  Step " << step_number << "  t=" << time;

            if (track_energy)
            {
                const double Ek    = compute_kinetic_energy();
                const double Ep    = compute_potential_energy();
                const double Etot  = Ek + Ep;
                const double Linf  = compute_Linfty_norm();

                // Initialized E0 at first output step (after ICs) to avoid startup transients
                if (!energy_initialized)
                {
                    energy_initial     = Etot;
                    energy_initialized = true;
                }

                // Relative drift: measures energy conservation
                const double drift = (energy_initial > 1e-30)
                                     ? std::abs(Etot - energy_initial) / energy_initial
                                     : 0.0;

                pcout << "  Ek=" << Ek
                      << "  Ep=" << Ep
                      << "  E=" << Etot
                      << "  drift=" << drift
                      << "  |u|_inf=" << Linf;

                if (this_mpi_proc == 0)
                    energy_log << step_number << ","
                               << time       << ","
                               << Ek         << ","
                               << Ep         << ","
                               << Etot       << ","
                               << drift      << ","
                               << Linf       << "\n";

                // Blow-up detection: whether |u|_inf > 100 * initial amplitude the simulation is probably unstable (CFL violated)
                if (Linf > blowup_threshold)
                    pcout << "\n  *** WARNING: numeric blow-up! (|u|_inf="
                        << Linf << " > " << blowup_threshold << ") ***\n";
            }

            if (mode == SimulationMode::MMS_CONVERGENCE)
            {
                auto [L2, H1] = compute_errors(time);
                pcout << "  L2=" << L2 << "  H1=" << H1;
            }
            pcout << "\n";

            output_results(step_number);

            if (this_mpi_proc == 0)
            {
                const std::string pvtu_name = "solution-"
                    + Utilities::int_to_string(step_number, 4) + ".pvtu";
                pvd_list.push_back({time, pvtu_name});
                std::ofstream pvd("solution.pvd");
                DataOutBase::write_pvd_record(pvd, pvd_list);
            }
        }
    }

    if (track_energy) 
    {
        if (this_mpi_proc == 0) energy_log.close();
        write_energy_report(); 
    }

    pcout << "  Simulation completed. Steps: " << step_number << "\n";
    computing_timer.print_summary();
}


template <int dim>
void WaveEquation<dim>::perform_single_mms_run(unsigned int ref, double dt)
{
    triangulation.clear();
    initial_refinement = ref;
    time_step = dt;

    //  Setup 
    make_grid();
    setup_system();
    assemble_matrices();

    // Setup IC
    if (use_decay_mms) 
    {
        InitialDisplacementMMS_Decay<dim> u0;
        InitialVelocityMMS_Decay<dim>     u1;
        VectorTools::interpolate(dof_handler, u0, owned_solution_u);
        VectorTools::interpolate(dof_handler, u1, owned_velocity_u);
    } 
    else 
    {
        InitialDisplacementMMS<dim> u0;
        InitialVelocityMMS<dim>     u1;
        VectorTools::interpolate(dof_handler, u0, owned_solution_u);
        VectorTools::interpolate(dof_handler, u1, owned_velocity_u);
    }

    constraints.distribute(owned_solution_u);
    solution_u = owned_solution_u;
    velocity_u = owned_velocity_u;

    solution_u.update_ghost_values();
    velocity_u.update_ghost_values();

    // Calculating the initial acceleration a0
    assemble_rhs(0.0);
    for (const auto idx : locally_owned_dofs)
    {
        if (constraints.is_constrained(idx))
        {
            owned_acceleration_u(idx) = 0.0;
            continue;
        }
        owned_acceleration_u(idx) = system_rhs(idx) / mass_matrix_diagonal(idx);
    }
    constraints.distribute(owned_acceleration_u);
    acceleration_u = owned_acceleration_u;
    acceleration_u.update_ghost_values();

    // DIAGNOSI (
    pcout << "    [diag a0] |a0|=" << owned_acceleration_u.l2_norm()
          << "  min(M)=" << mass_matrix_diagonal.min() << "\n";

    //  Initial Error Calculation (t=0)
    auto [L2_t0, H1_t0] = compute_errors(0.0);                  
    pcout << "  [t=0] L2=" << L2_t0 << "  H1=" << H1_t0 << "\n";

    // main Loop principale time-marching
    const unsigned int saved_output = output_every_n_steps;
    output_every_n_steps = std::numeric_limits<unsigned int>::max();

    time = 0.0; 
    step_number = 0;
    const unsigned int n_steps = static_cast<unsigned int>(end_time / time_step);
    
    while (step_number < n_steps) 
    {
        step_number++;
        time += time_step;
        if (time_scheme == TimeScheme::LEAPFROG) 
            solve_time_step();
        else 
            solve_time_step_newmark();
    }

    output_every_n_steps = saved_output;
    
    // final error
    auto [L2_end, H1_end] = compute_errors(time);

    convergence_table.add_value("Level", ref);
    convergence_table.add_value("h", 1.0 / std::pow(2.0, ref));
    convergence_table.add_value("L2_t0", L2_t0);
    convergence_table.add_value("H1_t0", H1_t0);
    convergence_table.add_value("L2_end", L2_end);
    convergence_table.add_value("H1_end", H1_end);
}

// ============================================================
// run_convergence_study  —  MMS on refining levels
template <int dim>
void WaveEquation<dim>::run_convergence_study()
{
    mode     = SimulationMode::MMS_CONVERGENCE;
    use_amr  = false;
    end_time = 0.05;                                        //1!

    // --- Studio 1: convergenza spaziale ---
    // dt = 0.1*h^2 => errore temporale O(dt) = O(0.1*h^2) << O(h^2)
    pcout << "\n=== SPATIAL ANALYSIS (dt = 0.1 * h^2) ===\n";
    for (unsigned int ref : {3u, 4u, 5u, 6u})
    {
        const double h  = 1.0 / std::pow(2.0, ref);
        const double dt = 0.1 * h * h;
        pcout << "  Level=" << ref << "  h=" << h
              << "  dt=" << dt
              << "  steps=" << (unsigned int)(end_time/dt) << "\n";
        perform_single_mms_run(ref, dt);
    }

    /*
    // --- Studio 2: convergenza temporale ---
    // h fissa (ref=5, h=0.03125), dt dimezza
    // Newmark trapezoidale: errore globale = O(dt^2) => rate atteso = 2
    pcout << "\n=== TEMPORAL ANALYSIS (ref=5 fisso) ===\n";
    for (double dt_val : {0.025, 0.0125, 0.00625, 0.003125})
    {
        pcout << "  dt=" << dt_val
              << "  steps=" << (unsigned int)(end_time/dt_val) << "\n";
        perform_single_mms_run(5, dt_val);
    }
        */
    
    
    //   convergenza temporale con T=0.5 e dt più piccoli
    end_time = 0.5; 
    pcout << "\n=== TEMPORAL ANALYSIS (ref=5 , T=0.5) ===\n";
    for (double dt_val : {0.05, 0.025, 0.0125, 0.00625})
    {
        pcout << "  dt=" << dt_val
              << "  steps=" << (unsigned int)(end_time/dt_val) << "\n";
        perform_single_mms_run(6, dt_val);
    }
    
    convergence_table.set_precision("h", 4);
    convergence_table.set_scientific("h", true);
    
    convergence_table.set_precision("L2_t0", 6);
    convergence_table.set_scientific("L2_t0", true);
    convergence_table.set_precision("H1_t0", 6);
    convergence_table.set_scientific("H1_t0", true);
    
    convergence_table.set_precision("L2_end", 6);
    convergence_table.set_scientific("L2_end", true);
    convergence_table.set_precision("H1_end", 6);
    convergence_table.set_scientific("H1_end", true);

    //  RATE FOR t=0 E t=end
    convergence_table.evaluate_convergence_rates("L2_t0", ConvergenceTable::reduction_rate_log2);
    convergence_table.evaluate_convergence_rates("H1_t0", ConvergenceTable::reduction_rate_log2);
    
    convergence_table.evaluate_convergence_rates("L2_end", ConvergenceTable::reduction_rate_log2);
    convergence_table.evaluate_convergence_rates("H1_end", ConvergenceTable::reduction_rate_log2);

    if (this_mpi_proc == 0)
    {
        pcout << "\n=== FINAL CONVERGENCE TABLE ===\n";
        convergence_table.write_text(std::cout);
    }

    computing_timer.print_summary();
}


//------------------------------------------------------------------------------------
// measure_numerical_phase_speed
//
// Launch a short simulation (n_wave cycle periods) with IC = flat wave sin(k·x).
// Measures the accumulated phase shift by comparing u_h(T) with  the shifted flat wave of various phase values, finding the maximum correlation (phase matching).
// Return c_h = ω_h / k (numeric phase speed).
template <int dim>
double WaveEquation<dim>::measure_numerical_phase_speed(
    double k,
    double n_periods)
{
    // Saves the current state (restored at the end)
    // Note: This method does NOT change the permanent status of the class, We use local copies of vectors for dispersion simulation.

    const double omega_exact = c * k;           // ω exact from dispersion relation of continuous wave equation
    const double T_period    = 2.0 * M_PI / omega_exact; // period
    const double T_sim       = n_periods * T_period;

    // dt for dispersion simulation: CFL with margin
    const double h_local = 1.0 / std::pow(2.0, initial_refinement);
    const double dt_disp = 0.4 * h_local / c;
    const unsigned int n_steps = static_cast<unsigned int>(T_sim / dt_disp) + 1;

    // Local vectors for dispersion simulation
    TrilinosVector d_u(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    TrilinosVector d_u_old(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    //TrilinosVector d_u_new(locally_owned_dofs, locally_relevant_dofs, mpi_comm);      --> non li usiamo alla fine 
    //TrilinosVector d_vel(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    TrilinosVector d_rhs(locally_owned_dofs, mpi_comm);

    //INITIAL CONDITIONS: palne wave sin(k·x)   
    PlaneWave<dim> pw(k, c, 0.0);
    pw.set_time(0.0);

    {
        TrilinosVector tmp(locally_owned_dofs, mpi_comm);
        VectorTools::interpolate(dof_handler, pw, tmp);
        constraints.distribute(tmp);
        d_u = tmp;
        d_u_old = d_u;
        //d_vel   = 0.0; !!!!!!
    }
    d_u.update_ghost_values();
    d_u_old.update_ghost_values();

    // Mini loop Leapfrog 
    for (unsigned int s = 0; s < n_steps; ++s)
    {
        d_rhs = 0.0;

        // RHS = -K·u
        laplace_matrix.vmult(d_rhs, d_u);
        d_rhs *= -1.0;

        // a = M^{-1} · RHS
        TrilinosVector d_acc(locally_owned_dofs, mpi_comm);
        for (const auto idx : locally_owned_dofs)
            d_acc(idx) = d_rhs(idx) / mass_matrix_diagonal(idx);

        // u_new = 2u - u_old + dt²·a
        TrilinosVector d_owned_new(locally_owned_dofs, mpi_comm);
        for (const auto idx : locally_owned_dofs)
            d_owned_new(idx) = 2.0 * d_u(idx)
                             - d_u_old(idx)
                             + dt_disp * dt_disp * d_acc(idx);

        constraints.distribute(d_owned_new);
        d_u_old = d_u;
        d_u     = d_owned_new;
        d_u.update_ghost_values();
        d_u_old.update_ghost_values();
    }

    //  Measurement of phase velocity by correlation 
    // We look for the value of φ in [0, 2π] that maximizes  C(φ) = <u_h(T), sin(k·x - ω_exact·T + φ)>
    // The φ optimal is the phase shift accumulated by the numerical solution.
    // The numerical phase speed is:
    // c_h = c · (1 - φ / (ω_exact · T_sim))
    //     = c · (1 - phase-out_per_cycle / (2π))

    const int n_phi = 360;  // Angular resolution: 1 degree
    double best_corr = -2.0;
    double best_phi  = 0.0;

    for (int ip = 0; ip < n_phi; ++ip)
    {
        const double phi = 2.0 * M_PI * ip / n_phi;

        PlaneWave<dim> ref(k, c, phi);
        ref.set_time(T_sim);

        // Calculate local correlation <u_h, ref>
        // by numerical integration to the nodes (lumbusted mass)
        TrilinosVector ref_vec(locally_owned_dofs, mpi_comm);
        VectorTools::interpolate(dof_handler, ref, ref_vec);

        double local_corr = 0.0;
        double local_norm_h = 0.0, local_norm_r = 0.0;
        for (const auto idx : locally_owned_dofs)
        {
            local_corr   += mass_matrix_diagonal(idx) * d_u(idx) * ref_vec(idx);
            local_norm_h += mass_matrix_diagonal(idx) * d_u(idx) * d_u(idx);
            local_norm_r += mass_matrix_diagonal(idx) * ref_vec(idx) * ref_vec(idx);
        }
        const double corr   = Utilities::MPI::sum(local_corr,   mpi_comm);
        const double norm_h = Utilities::MPI::sum(local_norm_h, mpi_comm);
        const double norm_r = Utilities::MPI::sum(local_norm_r, mpi_comm);

        const double normalized = corr / (std::sqrt(norm_h * norm_r) + 1e-30);
        if (normalized > best_corr)
        {
            best_corr = normalized;
            best_phi  = phi;
        }
    }

    // Total phase shift accumulated in T_sim
    // The numerical phase speed meets: ω_h · T_sim = ω_exact · T_sim - φ
    const double omega_numerical = omega_exact - best_phi / T_sim;
    const double c_numerical     = omega_numerical / k;

    return c_numerical;
}

//-------------------------------------------------------
// run_dispersion_analysis
//
// Head the dispersion for different wave numbers k.
// The useful range is kh ∈ [π/N, π] where N = 2^ref.
// kh = π is the Nyquist (the most dispersive) mode.
// kh → 0 is the continuous limit (error → 0).
template <int dim>
std::vector<typename WaveEquation<dim>::DispersionResult>
WaveEquation<dim>::run_dispersion_analysis(
    const std::vector<double> &wave_numbers)
{
    pcout << "\n=== Analisys Numerical Dispersion ===\n";
    pcout << "  Scheme: "
          << (time_scheme == TimeScheme::LEAPFROG ? "Leapfrog" : "Newmark")
          << "  Q" << fe_degree << "\n";
    pcout << "  Refinement: " << initial_refinement
          << "  (h = " << 1.0/std::pow(2.0,initial_refinement) << ")\n\n";

    const double h = 1.0 / std::pow(2.0, initial_refinement);

    std::vector<DispersionResult> results;
    results.reserve(wave_numbers.size());

    for (double k : wave_numbers)
    {
        const double kh = k * h;

        // Use 3 periods: enough to measure the phase shift, not too many to accumulate amplitude errors
        const double c_h = measure_numerical_phase_speed(k, 3.0);

        DispersionResult r;
        r.k             = k;
        r.kh            = kh;
        r.c_numerical   = c_h;
        r.c_exact       = c;
        r.relative_error = (c_h - c) / c;
        results.push_back(r);

        pcout << "  k=" << k
              << "  kh=" << kh
              << "  c_h=" << c_h
              << "  c_exact=" << c
              << "  relative_error=" << r.relative_error * 100.0 << "%\n";
    }

    // Write CSV (only rank 0)
    if (this_mpi_proc == 0)
    {
        std::ofstream f("dispersion.csv");
        f << "k,kh,c_numerical,c_exact,relative_error_pct\n";
        for (const auto &r : results)
            f << r.k << ","
              << r.kh << ","
              << r.c_numerical << ","
              << r.c_exact << ","
              << r.relative_error * 100.0 << "\n";

        pcout << "\n  Results saved in dispersion.csv\n";
        pcout << "  Use plot_dispersion.py to visualize.\n";
    }

    return results;
}

template <int dim>
void WaveEquation<dim>::build_newmark_system_matrix()
{
    TimerOutput::Scope t(computing_timer, "build_newmark_matrix");
    pcout << "  Building Newmark matrix A = M + β·dt²·K...\n";

    AssertThrow(time_scheme == TimeScheme::NEWMARK,
        ExcMessage("build_newmark_system_matrix without NEWMARK"));

    // A = β·dt²·K
    system_matrix_newmark.copy_from(laplace_matrix);
    system_matrix_newmark *= (newmark_beta * time_step * time_step);

    // A += M_lumped 
    for (const auto idx : locally_owned_dofs)
        system_matrix_newmark.add(idx, idx, mass_matrix_diagonal(idx));

    system_matrix_newmark.compress(VectorOperation::add);

    //  ILU PRECONDITIONER
    TrilinosWrappers::PreconditionILU::AdditionalData ilu_data;
    newmark_preconditioner.initialize(system_matrix_newmark, ilu_data);

    newmark_matrix_is_current = true;
    pcout << "  Matrix Newmark and ILU preconditioner ready.\n";
}



// ============================================================
// run_scaling_benchmark
// Runs n_steps time steps without output and measures wall-clock time
template <int dim>
typename WaveEquation<dim>::ScalingResult
WaveEquation<dim>::run_scaling_benchmark(unsigned int n_steps)

{
    ScalingResult result;
    result.n_procs = n_mpi_procs;
    result.n_dofs  = dof_handler.n_dofs();
    result.n_cells = triangulation.n_global_active_cells();

    // Timer reset 
    computing_timer.reset();

    // Warm-up: 1 inizialization cache and pipeline MPI
    {
        if (time_scheme == TimeScheme::LEAPFROG)
            solve_time_step();
        else
            solve_time_step_newmark();
    }
    // reset after the warmup
    computing_timer.reset();

    // Benchmark: n_steps 
    const auto t_start = std::chrono::high_resolution_clock::now();

    for (unsigned int s = 0; s < n_steps; ++s)
    {
        time += time_step;
        if (time_scheme == TimeScheme::LEAPFROG)
            solve_time_step();
        else
            solve_time_step_newmark();
    }

    const auto t_end = std::chrono::high_resolution_clock::now();
    result.wall_time_total =
        std::chrono::duration<double>(t_end - t_start).count();

    // Recover times per section from TimerOutput (labels must match those used in TimerOutput::Scope)
    const auto &summary = computing_timer.get_summary_data(
        TimerOutput::total_wall_time);

    result.wall_time_assembly = summary.count("assemble_rhs")
                              ? summary.at("assemble_rhs")
                              : 0.0;
    result.wall_time_solve =
            (time_scheme == TimeScheme::LEAPFROG)
            ? (summary.count("solve_leapfrog") ? summary.at("solve_leapfrog") : 0.0)
            : (summary.count("solve_newmark")  ? summary.at("solve_newmark")  : 0.0);
    result.wall_time_output   = 0.0;  // no output for benchmark

    // Speedup and efficiency are calculated off (require T(1))
    result.speedup    = 0.0;
    result.efficiency = 0.0;

    return result;
}

// ------------------------------------------------------------
// print_scaling_table
template <int dim>
void WaveEquation<dim>::print_scaling_table(
    const std::vector<ScalingResult> &results,
    const std::string &label) const
{
    if (this_mpi_proc != 0) return;

    pcout << "\n=== " << label << " ===\n";
    pcout << std::setw(8)  << "Procs"
          << std::setw(12) << "DoF"
          << std::setw(10) << "T_tot[s]"
          << std::setw(10) << "T_asm[s]"
          << std::setw(10) << "T_slv[s]"
          << std::setw(10) << "Speedup"
          << std::setw(12) << "Efficiency"
          << "\n";
    pcout << std::string(72, '-') << "\n";

    for (const auto &r : results)
    {
        pcout << std::setw(8)  << r.n_procs
              << std::setw(12) << r.n_dofs
              << std::setw(10) << std::fixed << std::setprecision(3) << r.wall_time_total
              << std::setw(10) << r.wall_time_assembly
              << std::setw(10) << r.wall_time_solve
              << std::setw(10) << std::setprecision(2) << r.speedup
              << std::setw(11) << r.efficiency * 100.0 << "%"
              << "\n";
    }

    // Saving
    std::ofstream f(label == "Strong Scaling" ? "strong_scaling.csv"
                                              : "weak_scaling.csv");
    f << "n_procs,n_dofs,n_cells,T_total,T_assembly,T_solve,speedup,efficiency\n";
    for (const auto &r : results)
        f << r.n_procs << ","
          << r.n_dofs  << ","
          << r.n_cells << ","
          << r.wall_time_total   << ","
          << r.wall_time_assembly << ","
          << r.wall_time_solve   << ","
          << r.speedup           << ","
          << r.efficiency        << "\n";

    pcout << "  Saved data in "
          << (label == "Strong Scaling" ? "strong_scaling.csv" : "weak_scaling.csv")
          << "\n";
}


// ============================================================
// Template instantiation for 2D and 3D
template class WaveEquation<2>;
template class WaveEquation<3>;
