#include "WaveEquationParallel.hpp"

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

    // Each process only marks its local cells
    for (auto &cell : triangulation.active_cell_iterators())
    {
        if (!cell->is_locally_owned()) continue;

        const Point<dim> center = cell->center();
        const bool in_wall_x = (center[0] > 0.48 && center[0] < 0.52);
        const bool in_gap_y  = (center[1] > 0.35 && center[1] < 0.65);
        if (in_wall_x && !in_gap_y)
            cell->set_material_id(1);
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
    // AffineConstraints uses locally_relevant_dofs such as IndexSet
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
        for (auto &cell : dof_handler.active_cell_iterators())
        {
            if (!cell->is_locally_owned()) continue;
            if (cell->material_id() != 1) continue;
            std::vector<types::global_dof_index> dof_ids(fe_ptr->dofs_per_cell);
            cell->get_dof_indices(dof_ids);
            for (auto idx : dof_ids)
                if (locally_relevant_dofs.is_element(idx)){
                    constraints.add_line(idx);
                    constraints.set_inhomogeneity(idx, 0.0);}
        }
    }
    constraints.close();

    // Sparsity pattern distributed 
    // DynamicSparsityPattern on locally_relevant_dofs, then We distribute it to remote processes with SparsityTools
    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
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

    pcout << "  Assembly matrici (MPI, " << n_mpi_procs << " processi)...\n";

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

        // Obstacle cells: no physical contribution, but we must still call distribute_local_to_global
        if (cell->material_id() == 1)
        {
            cell->get_dof_indices(local_idx);
            constraints.distribute_local_to_global(cell_K, local_idx, laplace_matrix);
            continue;
        }

        // matrix  K with c²(x) variable
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

        // distribute_local_to_global menage both communication and constraints application
        constraints.distribute_local_to_global(cell_K, local_idx, laplace_matrix);

        // Diagonal mass: we add directly to the global vector
        for (unsigned int i = 0; i < dpc; ++i)
            mass_matrix_diagonal(local_idx[i]) += cell_M(i);

        // Matric for ABC (absorbing boundary conditions):
        // matrix B = ∫_Γ φ_i φ_j ds 
        if (use_absorbing_bc)
        {
            for (const auto &face : cell->face_iterators())
            {
                if (!face->at_boundary()) continue;
                fev_face.reinit(cell, face);
                cell_B = 0.0;
                for (unsigned int q = 0; q < q_face.size(); ++q)
                {
                    const double JxW = fev_face.JxW(q);
                    for (unsigned int i = 0; i < dpc; ++i)
                        for (unsigned int j = 0; j < dpc; ++j)
                            cell_B(i, j) += fev_face.shape_value(i, q)
                                          * fev_face.shape_value(j, q)
                                          * JxW;
                }
                constraints.distribute_local_to_global(
                    cell_B, local_idx, boundary_mass_matrix);
            }
        }
    }

    // MPI communication : somma contributi di tutti i processi , it adds ocontributions from remote processes to the local rows and makes the result available on all processes
    // compress(add): every process send his "remote" rows to the others, and sums the contributions in the local rows
    laplace_matrix.compress(VectorOperation::add);
    mass_matrix_diagonal.compress(VectorOperation::add);
    if (use_absorbing_bc)
        boundary_mass_matrix.compress(VectorOperation::add);

    // mass positive check (global minimum via MPI_Allreduce)
    const double local_min = mass_matrix_diagonal.min();
    const double global_min = Utilities::MPI::min(local_min, mpi_comm);
    AssertThrow(global_min > 0.0,
        ExcMessage("Lumped mass matrix: entry diagonal non positive!"));

    pcout << "  Assembly completed. min(M_diag)=" << global_min << "\n";


    // FOR NEWMARK: build system matrix A = M + β·dt²·K
    if (time_scheme == TimeScheme::NEWMARK)
    {
        newmark_matrix_is_current = false;  // forza la ricostruzione
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

    // Update the ghost values of solution_u before vmult
    solution_u.update_ghost_values();

    // −K·u  (K already includes c²)
    laplace_matrix.vmult(system_rhs, owned_solution_u);
    system_rhs *= -1.0;

    /*
    (Se usi il damping, assicurati di usare owned_velocity_u(idx) anche nel ciclo for del damping poco sotto).
    */

    // Forcing term f(x,t) for MMS
    if (mode == SimulationMode::MMS_CONVERGENCE)
    {
        QGauss<dim> q(fe_ptr->degree + 1);
        FEValues<dim> fev(*fe_ptr, q,
            update_values | update_quadrature_points | update_JxW_values);

        ForcingTermMMS<dim> forcing(c);
        forcing.set_time(t);

        const unsigned int dpc = fe_ptr->dofs_per_cell;
        Vector<double> cell_rhs(dpc);
        std::vector<types::global_dof_index> local_idx(dpc);

        for (const auto &cell : dof_handler.active_cell_iterators())
        {
            if (!cell->is_locally_owned()) continue;
            fev.reinit(cell);
            cell_rhs = 0.0;
            for (unsigned int q = 0; q < fev.get_quadrature().size(); ++q)
            {
                const double fval = forcing.value(fev.quadrature_point(q));
                for (unsigned int i = 0; i < dpc; ++i)
                    cell_rhs(i) += fev.shape_value(i, q) * fval * fev.JxW(q);
            }
            cell->get_dof_indices(local_idx);
            // add_local_to_global on Trilinos vector 
            constraints.distribute_local_to_global(cell_rhs, local_idx, system_rhs);
        }
        system_rhs.compress(VectorOperation::add);
    }

    // ABC: −c·B·v
    if (use_absorbing_bc)
    {
        velocity_u.update_ghost_values();
        TrilinosVector abc_contrib(locally_owned_dofs, mpi_comm);
        boundary_mass_matrix.vmult(abc_contrib, velocity_u);
        system_rhs.add(-c, abc_contrib);
    }

    // Damping: −d·M·v  (M diagonal, local operation)
    if (damping > 0.0)
    {
        velocity_u.update_ghost_values();
        for (const auto idx : locally_owned_dofs)
            system_rhs(idx) -= damping * mass_matrix_diagonal(idx) * velocity_u(idx);
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

    assemble_rhs(time);

    // Update ghost values for the vectors we read
    solution_u.update_ghost_values();
    solution_u_old.update_ghost_values();

    // a_i = rhs_i / M_ii  — local operation (only owned dofs)
    for (const auto idx : locally_owned_dofs)
    {
        const double m_ii = mass_matrix_diagonal(idx);
        owned_acceleration_u(idx) = system_rhs(idx) / m_ii;
    }

    // u^{n+1} = 2·u^n − u^{n-1} + dt²·a
    const double dt2 = time_step * time_step;
    for (const auto idx : locally_owned_dofs)
    {
        owned_solution_u(idx) = 2.0 * solution_u(idx)
                              - solution_u_old(idx)
                              + dt2 * owned_acceleration_u(idx);
    }

    //Apply Dirichlet constraints: local operation
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

    // Security check: the matrix must be updated
    AssertThrow(newmark_matrix_is_current,
        ExcMessage("build_newmark_system_matrix() not called! "
                   "Call assemble_matrices() before the time loop."));

    const double dt  = time_step;
    const double b   = newmark_beta;
    const double gam = newmark_gamma;

    solution_u.update_ghost_values();
    velocity_u.update_ghost_values();
    acceleration_u.update_ghost_values();

    // predictors u_pred, v_pred (owned only, no ghost)
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

    // RHS_newmark = −K·u_pred 
    //TrilinosVector u_pred_g(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    //u_pred_g = u_pred;
    //u_pred_g.update_ghost_values();

    TrilinosVector rhs_newmark(locally_owned_dofs, mpi_comm);
    //laplace_matrix.vmult(rhs_newmark, u_pred_g);
    laplace_matrix.vmult(rhs_newmark, u_pred);
    rhs_newmark *= -1.0;

    //  Solve CG wuth AMG  preconditioner
    // system_matrix_newmark it's already built and the AMG preconditioner is initialized, so no overhead in the time loop
    TrilinosVector a_new(locally_owned_dofs, mpi_comm);

    SolverControl solver_control(500, 1e-10 * rhs_newmark.l2_norm() + 1e-30);
    TrilinosWrappers::SolverCG cg_solver(solver_control);

    // Uses the preconditioner that is already initialized — no overhead
    cg_solver.solve(system_matrix_newmark,
                    a_new,
                    rhs_newmark,
                    newmark_preconditioner);

    constraints.distribute(a_new);

    // Log iterations only every output_every_n_steps per not "fill" stdout
    if (step_number % output_every_n_steps == 0)
        pcout << "    Newmark CG: " << solver_control.last_step()
              << " iterazioni (AMG)\n";

    // CORRECTORS
    for (const auto idx : locally_owned_dofs)
    {
        owned_solution_u(idx)    = u_pred(idx) + b * dt * dt * a_new(idx);
        owned_velocity_u(idx)    = v_pred(idx) + gam * dt * a_new(idx);
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
double WaveEquation<dim>::compute_potential_energy() //const
{
    TrilinosVector Ku(locally_owned_dofs, mpi_comm);
    solution_u.update_ghost_values(); 
    laplace_matrix.vmult(Ku, owned_solution_u);

    double local_ep = 0.0;
    for (const auto idx : locally_owned_dofs)
        local_ep += 0.5 * owned_solution_u(idx) * Ku(idx);

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
std::pair<double, double> WaveEquation<dim>::compute_errors(double t) const
{
    ExactSolutionMMS<dim> exact(c);
    exact.set_time(t);

    Vector<float> diff(triangulation.n_active_cells());

    VectorTools::integrate_difference(
        dof_handler, solution_u, exact, diff,
        QGauss<dim>(fe_ptr->degree + 2),
        VectorTools::L2_norm);
    const double L2 = VectorTools::compute_global_error(
        triangulation, diff, VectorTools::L2_norm);

    VectorTools::integrate_difference(
        dof_handler, solution_u, exact, diff,
        QGauss<dim>(fe_ptr->degree + 2),
        VectorTools::H1_seminorm);
    const double H1 = VectorTools::compute_global_error(
        triangulation, diff, VectorTools::H1_seminorm);

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
// parallel::distributed::GridRefinement distributes flags  between processes and p4est manages load rebalancing.
// parallel::distributed::SolutionTransfer transfers the solution on the new distributed mesh.
template <int dim>
void WaveEquation<dim>::refine_mesh()
{
    TimerOutput::Scope t(computing_timer, "refine_mesh");

    pcout << "  AMR parallelo...\n";

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

    // Limit the maximum level
    for (auto &cell : triangulation.active_cell_iterators())
        if (cell->is_locally_owned() &&
            cell->level() >= static_cast<int>(max_refinement_level))
            cell->clear_refine_flag();

    triangulation.prepare_coarsening_and_refinement();

    // SolutionTransfer distrributed
    // Note: parallel::distributed::SolutionTransfer wants
    // only one vector per call. We use three transfers
    // separated or a chained vector (we use three calls).
    parallel::distributed::SolutionTransfer<dim, TrilinosVector>
        st_u(dof_handler), st_u_old(dof_handler), st_v(dof_handler);

    solution_u.update_ghost_values();
    solution_u_old.update_ghost_values();
    velocity_u.update_ghost_values();

    st_u.prepare_for_coarsening_and_refinement(solution_u);
    st_u_old.prepare_for_coarsening_and_refinement(solution_u_old);
    st_v.prepare_for_coarsening_and_refinement(velocity_u);

    triangulation.execute_coarsening_and_refinement();

    // Redistributes DoF on new mesh
    dof_handler.distribute_dofs(*fe_ptr);
    locally_owned_dofs    = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    pcout << "  AMR: DoF globali=" << dof_handler.n_dofs()
          << "  celle=" << triangulation.n_global_active_cells() << "\n";

    // Rebuilds constraints
    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    if (!use_absorbing_bc)
        VectorTools::interpolate_boundary_values(
            dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
    constraints.close();

    // Rebuilds sparsity and matrices
    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern(
        dsp, locally_owned_dofs, mpi_comm, locally_relevant_dofs);

    laplace_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_comm);
    if (use_absorbing_bc)
        boundary_mass_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_comm);
    if (time_scheme == TimeScheme::NEWMARK)
        system_matrix_newmark.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_comm);

    // Riallocation(?) vectors on the new partition
    mass_matrix_diagonal.reinit(locally_owned_dofs, mpi_comm);
    system_rhs.reinit(locally_owned_dofs, mpi_comm);
    owned_solution_u.reinit(locally_owned_dofs, mpi_comm);
    owned_solution_u_old.reinit(locally_owned_dofs, mpi_comm);
    owned_velocity_u.reinit(locally_owned_dofs, mpi_comm);
    owned_acceleration_u.reinit(locally_owned_dofs, mpi_comm);

    solution_u.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    solution_u_old.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    solution_u_new.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    velocity_u.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    acceleration_u.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);

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

    // Riassemly matrices on a new mesh
    assemble_matrices();
    check_cfl_condition();
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
          << n_mpi_procs << " processi\n  Modo: ";
    switch (mode)
    {
        case SimulationMode::PEBBLE_IN_POND:  pcout << "Pebble in Pond\n";  break;
        case SimulationMode::MMS_CONVERGENCE: pcout << "MMS\n";             break;
        case SimulationMode::DAMPED_WAVE:     pcout << "Smorzato\n";        break;
        case SimulationMode::ABSORBING_BC:    pcout << "ABC Sommerfeld\n";  break;
        case SimulationMode::INTERFERENCE:    pcout << "Interferenza\n";    break;
        case SimulationMode::REFRACTION:      pcout << "Rifrazione\n";      break;
        case SimulationMode::DIFFRACTION:     pcout << "Diffrazione\n";     break;
    }
    pcout << "  Schema: "
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
    pcout << "  Condizioni iniziali...\n";

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
        const double width = 0.05;
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

    output_results(0);

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
                if (Linf > 100.0)
                    pcout << "\n  *** WARNING: numeric blow-up! "
                             "Riduci dt o verifica CFL. ***\n";
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

// ============================================================
// run_convergence_study  —  MMS on refining levels
template <int dim>
void WaveEquation<dim>::run_convergence_study()
{
    pcout << "\n=== Convergence Study MMS ("
          << n_mpi_procs << " MPI processes, Q"
          << fe_degree << ", "
          << (time_scheme == TimeScheme::LEAPFROG ? "Leapfrog" : "Newmark-β")
          << ") ===\n";

    mode    = SimulationMode::MMS_CONVERGENCE;
    use_amr = false;
    end_time = 1.0;

    const std::vector<unsigned int> levels = {3, 4, 5, 6, 7};

    for (unsigned int ref : levels)
    {
        newmark_matrix_is_current = false;
        pcout << "\n--- Level " << ref << " ---\n";

        triangulation.clear();
        initial_refinement = ref;

        const double h = 1.0 / std::pow(2.0, ref);
        time_step = (time_scheme == TimeScheme::LEAPFROG)
                  ? 0.4 * h / c
                  : 0.5 * h;
        const unsigned int n_steps =
            static_cast<unsigned int>(end_time / time_step);

        pcout << "  h=" << h << "  dt=" << time_step
              << "  steps=" << n_steps << "\n";

        make_grid();
        setup_system();
        assemble_matrices();

        // IC MMS
        InitialDisplacementMMS<dim> u0;
        InitialVelocityMMS<dim>     u1;
        VectorTools::interpolate(dof_handler, u0, owned_solution_u);
        constraints.distribute(owned_solution_u);
        solution_u = owned_solution_u;

        VectorTools::interpolate(dof_handler, u1, owned_velocity_u);
        for (const auto idx : locally_owned_dofs)
            owned_solution_u_old(idx) = owned_solution_u(idx)
                                      - time_step * owned_velocity_u(idx);
        solution_u_old = owned_solution_u_old;
        velocity_u     = owned_velocity_u;
        owned_acceleration_u = 0.0;
        acceleration_u       = owned_acceleration_u;

        solution_u.update_ghost_values();
        solution_u_old.update_ghost_values();
        velocity_u.update_ghost_values();

        // Loop without output
        time = 0.0; step_number = 0;
        while (step_number < n_steps)
        {
            step_number++;
            time += time_step;
            if (time_scheme == TimeScheme::LEAPFROG)
                solve_time_step();
            else
                solve_time_step_newmark();
        }

        auto [L2, H1] = compute_errors(time);
        pcout << "  L2=" << L2 << "  H1=" << H1 << "\n";

        convergence_table.add_value("Level", ref);
        convergence_table.add_value("Cells",   (unsigned int)triangulation.n_global_active_cells());
        convergence_table.add_value("DoF",     dof_handler.n_dofs());
        convergence_table.add_value("h",       h);
        convergence_table.add_value("L2",      L2);
        convergence_table.add_value("H1",      H1);
    }

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

    // Only rank 0 prints the table and saves it in convergence_table.txt
    if (this_mpi_proc == 0)
    {
        pcout << "\n=== Convergence Table ===\n";
        convergence_table.write_text(std::cout);
        std::ofstream f("convergence_table.txt");
        convergence_table.write_text(f);
        pcout << "Table saved in convergence_table.txt\n";
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


// ------------------------------------------------------------
// build_newmark_system_matrix
// Builds A = M_lump + β·dt2·K and initializes the
// pre-conditioner AMG (Algebraic MultiGrid).
template <int dim>
void WaveEquation<dim>::build_newmark_system_matrix()
{
    TimerOutput::Scope t(computing_timer, "build_newmark_matrix");

    pcout << "  Building Newmark matrix A = M + β·dt²·K...\n";

    AssertThrow(time_scheme == TimeScheme::NEWMARK,
        ExcMessage("build_newmark_system_matrix called without NEWMARK mode"));

    // A = β·dt²·K  
    system_matrix_newmark.copy_from(laplace_matrix);
    system_matrix_newmark *= newmark_beta * time_step * time_step;

    // A += M  on the diagonal (M is a vector because we use mass lumping)
    for (const auto idx : locally_owned_dofs)
        system_matrix_newmark.add(idx, idx, mass_matrix_diagonal(idx));

    // MPI communication for the off-diagonal entries (if any)
    system_matrix_newmark.compress(VectorOperation::add);

    // AMG preconditioner for A
    // AMG is more expensive to build (O(n log n)) but much more
    // efficient for the solve (O(n) vs O(n^1.5) for SSOR).
    // Building it once amortizes the cost over all time steps.
    TrilinosWrappers::PreconditionAMG::AdditionalData amg_data;
    amg_data.elliptic              = true;   // A is elliptic (SPD)
    amg_data.higher_order_elements = (fe_degree > 1);
    amg_data.smoother_sweeps       = 2;
    amg_data.aggregation_threshold = 1e-4;

    newmark_preconditioner.initialize(system_matrix_newmark, amg_data);

    newmark_matrix_is_current = true;

    pcout << "  Matrix Newmark and AMG preconditioner ready.\n";
}





// ============================================================
// run_scaling_benchmark
// Runs n_steps Leapfrog steps without output and measures times
// with TimerOutput. Back to a ScalingResult.
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
        assemble_rhs(0.0);
        solve_time_step();
    }
    // reset after the warmup
    computing_timer.reset();

    // Benchmark: n_steps 
    const auto t_start = std::chrono::high_resolution_clock::now();

    for (unsigned int s = 0; s < n_steps; ++s)
    {
        time += time_step;
        solve_time_step();  // include assemble_rhs()
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
    result.wall_time_solve    = summary.count("solve_leapfrog")
                              ? summary.at("solve_leapfrog")
                              : 0.0;
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
