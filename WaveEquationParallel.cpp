// ============================================================
// WaveEquation.cpp  —  implementazione MPI-parallela
//
// ARCHITETTURA MPI IN DEAL.II:
// ─────────────────────────────────────────────────────────────
// 1. parallel::distributed::Triangulation
//    La mesh è partizionata automaticamente tra i processi MPI
//    tramite p4est. Ogni processo conosce:
//      - "locally owned cells": celle che elabora
//      - "ghost cells": celle dei vicini, necessarie per assembly
//
// 2. IndexSet
//    locally_owned_dofs  : DoF di proprietà di questo processo
//    locally_relevant_dofs: owned + ghost (per leggere i valori
//                           dei vicini durante assembly/solve)
//
// 3. TrilinosWrappers::MPI::Vector
//    Vettori distribuiti. Esistono in due varianti:
//      - "owned only" (no ghost): per scrivere (assembly RHS)
//      - "with ghost"            : per leggere (assembly K, output)
//    La sincronizzazione avviene con update_ghost_values() e
//    compress(VectorOperation::add o ::insert).
//
// 4. TrilinosWrappers::SparseMatrix
//    Matrice distribuita per righe. Ogni processo possiede le
//    righe corrispondenti ai suoi locally_owned_dofs.
//    L'assembly locale (celle owned+ghost) è thread-safe:
//    Trilinos gestisce internamente la comunicazione off-process.
//
// 5. Schema Leapfrog (esplicito, no solve lineare):
//    - assemble_rhs(): ogni processo calcola la sua porzione di RHS
//    - compress(add): somma i contributi tra processi
//    - divisione locale: a_i = rhs_i / M_ii  (solo owned dofs)
//    - update Leapfrog: u_new = 2u - u_old + dt²·a
//    - constraints.distribute(): sincronizza i DoF vincolati
//
// 6. Schema Newmark (implicito, solve CG con Trilinos):
//    - Risolve (M + β·dt²·K)·a = RHS con SolverCG Trilinos
//    - Il solver è già MPI-aware: comunicazione automatica
//
// 7. AMR parallelo:
//    - parallel::distributed::GridRefinement
//    - parallel::distributed::SolutionTransfer
//    - Dopo refine, repartiziona automaticamente con p4est
//
// 8. Output parallelo:
//    - Ogni processo scrive il suo file .vtu
//    - Il rank 0 scrive il file .pvtu master (indice dei file)
//    - ParaView carica il .pvtu e ricostruisce la soluzione
// ============================================================

#include "WaveEquation.hpp"

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

#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>

// ============================================================
// Costruttore
// ============================================================
template <int dim>
WaveEquation<dim>::WaveEquation(MPI_Comm mpi_communicator)
    : mpi_comm(mpi_communicator)
    , n_mpi_procs(Utilities::MPI::n_mpi_processes(mpi_comm))
    , this_mpi_proc(Utilities::MPI::this_mpi_process(mpi_comm))
    , pcout(std::cout, this_mpi_proc == 0)   // stampa solo rank 0
    , computing_timer(mpi_comm,
                      pcout,
                      TimerOutput::never,
                      TimerOutput::wall_times)
    , triangulation(mpi_comm)                // mesh distribuita p4est
    , fe_ptr(std::make_unique<FE_Q<dim>>(1))
    , dof_handler(triangulation)
    , time(0.0)
    , step_number(0)
{}

// ============================================================
// wave_speed_at: c(x) — costante o eterogeneo (rifrazione)
// ============================================================
template <int dim>
double WaveEquation<dim>::wave_speed_at(const Point<dim> &p) const
{
    if (mode == SimulationMode::REFRACTION)
        return (p[1] > interface_y) ? c_fast : c_slow;
    return c;
}

// ============================================================
// make_grid
// ============================================================
template <int dim>
void WaveEquation<dim>::make_grid()
{
    TimerOutput::Scope t(computing_timer, "make_grid");

    pcout << "  Generazione griglia [0,1]^" << dim
          << "  (raffinamenti globali: " << initial_refinement << ")\n";

    GridGenerator::hyper_cube(triangulation, 0.0, 1.0);

    // refine_global su parallel::distributed::Triangulation
    // redistribuisce automaticamente le celle tra i processi
    triangulation.refine_global(initial_refinement);

    pcout << "  Celle attive (globale): "
          << triangulation.n_global_active_cells() << "\n";
}

// ============================================================
// make_grid_with_obstacle  (diffrazione)
// ============================================================
template <int dim>
void WaveEquation<dim>::make_grid_with_obstacle()
{
    TimerOutput::Scope t(computing_timer, "make_grid_obstacle");

    pcout << "  Generazione griglia con ostacolo (diffrazione)...\n";

    GridGenerator::hyper_cube(triangulation, 0.0, 1.0);
    triangulation.refine_global(initial_refinement);

    // Ogni processo marca solo le sue celle locali
    for (auto &cell : triangulation.active_cell_iterators())
    {
        if (!cell->is_locally_owned()) continue;

        const Point<dim> center = cell->center();
        const bool in_wall_x = (center[0] > 0.48 && center[0] < 0.52);
        const bool in_gap_y  = (center[1] > 0.35 && center[1] < 0.65);
        if (in_wall_x && !in_gap_y)
            cell->set_material_id(1);
    }

    pcout << "  Celle attive (globale): "
          << triangulation.n_global_active_cells() << "\n";
}

// ============================================================
// setup_system
// ============================================================
template <int dim>
void WaveEquation<dim>::setup_system()
{
    TimerOutput::Scope t(computing_timer, "setup_system");

    pcout << "  Setup sistema...\n";

    // Ricrea FE con il grado scelto
    fe_ptr = std::make_unique<FE_Q<dim>>(fe_degree);
    dof_handler.distribute_dofs(*fe_ptr);

    // ---- IndexSet: DoF owned e relevant ----
    // locally_owned_dofs: quelli che questo processo "possiede"
    locally_owned_dofs    = dof_handler.locally_owned_dofs();
    // locally_relevant_dofs: owned + ghost (i vicini che ci servono)
    locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

    pcout << "  DoF globali: " << dof_handler.n_dofs()
          << "  (questo processo: " << locally_owned_dofs.n_elements() << ")\n";

    // ---- Constraints ----
    // AffineConstraints usa locally_relevant_dofs come IndexSet
    constraints.clear();
    constraints.reinit(locally_relevant_dofs);

    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    if (!use_absorbing_bc)
    {
        // Dirichlet u=0 su tutto il bordo (boundary_id=0)
        VectorTools::interpolate_boundary_values(
            dof_handler, 0,
            Functions::ZeroFunction<dim>(),
            constraints);
    }

    // Ostacolo diffrazione: vincola i DoF delle celle con material_id=1
    if (mode == SimulationMode::DIFFRACTION)
    {
        for (auto &cell : dof_handler.active_cell_iterators())
        {
            if (!cell->is_locally_owned()) continue;
            if (cell->material_id() != 1) continue;
            std::vector<types::global_dof_index> dof_ids(fe_ptr->dofs_per_cell);
            cell->get_dof_indices(dof_ids);
            for (auto idx : dof_ids)
                if (locally_relevant_dofs.is_element(idx))
                    constraints.add_constraint(idx, {}, 0.0);
        }
    }
    constraints.close();

    // ---- Sparsity pattern distribuito ----
    // DynamicSparsityPattern sui locally_relevant_dofs,
    // poi lo distribuiamo ai processi remoti con SparsityTools
    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern(
        dsp,
        locally_owned_dofs,
        mpi_comm,
        locally_relevant_dofs);

    // ---- Alloca matrici Trilinos ----
    // Ogni matrice è distribuita: processo i possiede le righe
    // corrispondenti a locally_owned_dofs
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

    // ---- Alloca vettori ----
    // "owned only" (no ghost): per assembly (scrittura)
    mass_matrix_diagonal.reinit(locally_owned_dofs, mpi_comm);
    system_rhs.reinit(locally_owned_dofs, mpi_comm);
    owned_solution_u.reinit(locally_owned_dofs, mpi_comm);
    owned_solution_u_old.reinit(locally_owned_dofs, mpi_comm);
    owned_velocity_u.reinit(locally_owned_dofs, mpi_comm);
    owned_acceleration_u.reinit(locally_owned_dofs, mpi_comm);

    // "with ghost" (locally_relevant): per lettura durante assembly
    solution_u.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    solution_u_old.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    solution_u_new.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    velocity_u.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    acceleration_u.reinit(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
}

// ============================================================
// assemble_matrices
//
// Ogni processo itera solo sulle sue "locally owned cells".
// I contributi vengono accumulati localmente e poi Trilinos
// esegue la comunicazione con compress(VectorOperation::add).
// ============================================================
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

    // ---- Itera solo sulle celle locally owned ----
    for (const auto &cell : dof_handler.active_cell_iterators())
    {
        if (!cell->is_locally_owned()) continue;

        fev_stiff.reinit(cell);
        fev_mass.reinit(cell);

        cell_K = 0.0;
        cell_M = 0.0;

        // Celle ostacolo: nessun contributo fisico
        if (cell->material_id() == 1)
        {
            cell->get_dof_indices(local_idx);
            constraints.distribute_local_to_global(cell_K, local_idx, laplace_matrix);
            continue;
        }

        // Matrice di rigidezza K con c²(x) variabile
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

        // Massa lumpata (Gauss-Lobatto)
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

        // distribute_local_to_global gestisce sia la comunicazione
        // inter-processo che l'applicazione dei constraints
        constraints.distribute_local_to_global(cell_K, local_idx, laplace_matrix);

        // Massa diagonale: sommiamo direttamente sul vettore globale
        for (unsigned int i = 0; i < dpc; ++i)
            mass_matrix_diagonal(local_idx[i]) += cell_M(i);

        // Matrice di bordo per ABC
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

    // ---- Comunicazione MPI: somma contributi di tutti i processi ----
    // compress(add): ogni processo invia le sue righe "remote" agli altri
    laplace_matrix.compress(VectorOperation::add);
    mass_matrix_diagonal.compress(VectorOperation::add);
    if (use_absorbing_bc)
        boundary_mass_matrix.compress(VectorOperation::add);

    // Controllo positività massa (su tutti i processi via MPI_Allreduce)
    const double local_min = mass_matrix_diagonal.min();
    const double global_min = Utilities::MPI::min(local_min, mpi_comm);
    AssertThrow(global_min > 0.0,
        ExcMessage("Massa lumpata: entry diagonale non positiva!"));

    pcout << "  Assembly completato. min(M_diag)=" << global_min << "\n";
}

// ============================================================
// assemble_rhs
//
// Calcola:  RHS = −K·u  +  F(t)  −  c·B·v  −  d·M·v
//
// Nota: laplace_matrix.vmult() e boundary_mass_matrix.vmult()
// usano internamente MPI_Allreduce per sommare i contributi
// dei processi ghost. Il risultato è già distribuito.
// ============================================================
template <int dim>
void WaveEquation<dim>::assemble_rhs(double t)
{
    system_rhs = 0.0;

    // Aggiorna i ghost values di solution_u prima di vmult
    solution_u.update_ghost_values();

    // −K·u  (K include già c²)
    laplace_matrix.vmult(system_rhs, solution_u);
    system_rhs *= -1.0;

    // Forcing term f(x,t) per MMS
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
            // add_local_to_global su vettore Trilinos
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

    // Smorzamento: −d·M·v  (M diagonale, operazione locale)
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
//
// Con massa lumpata M^{-1} è divisione scalare locale:
// nessuna comunicazione MPI per il solve.
// La comunicazione avviene solo in assemble_rhs (vmult).
// ============================================================
template <int dim>
void WaveEquation<dim>::solve_time_step()
{
    TimerOutput::Scope t(computing_timer, "solve_leapfrog");

    assemble_rhs(time);

    // Aggiorna ghost per i vettori che leggiamo
    solution_u.update_ghost_values();
    solution_u_old.update_ghost_values();

    // a_i = rhs_i / M_ii  — operazione locale (solo owned dofs)
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

    // Applica constraints Dirichlet: operazione locale
    if (!use_absorbing_bc)
        constraints.distribute(owned_solution_u);

    // Velocità centrata: v^n = (u^{n+1} − u^{n-1}) / (2·dt)
    for (const auto idx : locally_owned_dofs)
        owned_velocity_u(idx) = (owned_solution_u(idx) - solution_u_old(idx))
                              / (2.0 * time_step);

    // ---- Shift: avanza di un passo ----
    // Copia owned → ghost (update_ghost_values propaga ai vicini)
    owned_solution_u_old = owned_solution_u;

    solution_u_old = solution_u;    // u^{n-1} <- u^n
    solution_u     = owned_solution_u;  // u^n <- u^{n+1}
    velocity_u     = owned_velocity_u;

    // Rende i nuovi ghost values disponibili per il prossimo step
    solution_u.update_ghost_values();
    solution_u_old.update_ghost_values();
    velocity_u.update_ghost_values();
}

// ============================================================
// solve_time_step_newmark  —  schema Newmark-β implicito
//
// Risolve con SolverCG (Trilinos):
//   A·a^{n+1} = RHS_newmark
//   A = M_lumped + β·dt²·K
//
// Il solver CG di Trilinos è MPI-aware: le operazioni
// dot-product e axpy sono automaticamente ridotte su tutti
// i processi tramite MPI_Allreduce internamente.
// ============================================================
template <int dim>
void WaveEquation<dim>::solve_time_step_newmark()
{
    TimerOutput::Scope t(computing_timer, "solve_newmark");

    const double dt  = time_step;
    const double b   = newmark_beta;
    const double gam = newmark_gamma;

    solution_u.update_ghost_values();
    velocity_u.update_ghost_values();
    acceleration_u.update_ghost_values();

    // ---- Predittori ----
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

    // ---- RHS_newmark = −K·u_pred ----
    TrilinosVector u_pred_ghosted(locally_owned_dofs, locally_relevant_dofs, mpi_comm);
    u_pred_ghosted = u_pred;
    u_pred_ghosted.update_ghost_values();

    TrilinosVector rhs_newmark(locally_owned_dofs, mpi_comm);
    laplace_matrix.vmult(rhs_newmark, u_pred_ghosted);
    rhs_newmark *= -1.0;

    // ---- Matrice di sistema A = M_diag + β·dt²·K ----
    // Ricalcoliamo A solo se necessario (primo step o dopo AMR)
    // Per semplicità la costruiamo a ogni step (overhead accettabile
    // se il numero di step è grande rispetto ai DoF).
    system_matrix_newmark.copy_from(laplace_matrix);
    system_matrix_newmark *= b * dt * dt;
    // Aggiunge M sulla diagonale: A_ii += M_ii
    for (const auto idx : locally_owned_dofs)
        system_matrix_newmark.add(idx, idx, mass_matrix_diagonal(idx));
    system_matrix_newmark.compress(VectorOperation::add);

    // ---- Solve CG con precondizionatore Trilinos SSOR ----
    TrilinosVector a_new(locally_owned_dofs, mpi_comm);
    SolverControl solver_control(2000, 1e-10 * rhs_newmark.l2_norm() + 1e-30);
    TrilinosWrappers::SolverCG solver(solver_control);
    TrilinosWrappers::PreconditionSSOR preconditioner;
    preconditioner.initialize(system_matrix_newmark);
    solver.solve(system_matrix_newmark, a_new, rhs_newmark, preconditioner);
    constraints.distribute(a_new);

    pcout << "    Newmark CG: " << solver_control.last_step() << " iterazioni\n";

    // ---- Correttori ----
    for (const auto idx : locally_owned_dofs)
    {
        owned_solution_u(idx) = u_pred(idx) + b * dt * dt * a_new(idx);
        owned_velocity_u(idx) = v_pred(idx) + gam * dt * a_new(idx);
        owned_acceleration_u(idx) = a_new(idx);
    }
    constraints.distribute(owned_solution_u);

    // Shift
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
//
// E_k = 0.5 · v^T · M · v
// Con M diagonale e v distribuito, ogni processo calcola
// la somma locale, poi MPI_Allreduce somma globalmente.
// ============================================================
template <int dim>
double WaveEquation<dim>::compute_kinetic_energy() const
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
// laplace_matrix.vmult() già distribuito.
// ============================================================
template <int dim>
double WaveEquation<dim>::compute_potential_energy() const
{
    TrilinosVector Ku(locally_owned_dofs, mpi_comm);
    solution_u.update_ghost_values();  // const_cast non necessario: già chiamato
    laplace_matrix.vmult(Ku, solution_u);

    double local_ep = 0.0;
    for (const auto idx : locally_owned_dofs)
        local_ep += 0.5 * solution_u(idx) * Ku(idx);

    return Utilities::MPI::sum(local_ep, mpi_comm);
}

// ============================================================
// compute_Linfty_norm
// Calcola max|u_i| sui DoF locally owned, poi MPI_Allreduce.
// Utile per rilevare instabilità numerica (CFL violata).
// ============================================================
template <int dim>
double WaveEquation<dim>::compute_Linfty_norm() const
{
    double local_max = 0.0;
    for (const auto idx : locally_owned_dofs)
        local_max = std::max(local_max, std::abs(solution_u(idx)));

    return Utilities::MPI::max(local_max, mpi_comm);
}

// ============================================================
// compute_errors (MMS)
// VectorTools::integrate_difference è già parallelo:
// calcola l'errore locale e poi MPI_Allreduce.
// ============================================================
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
// h_min è calcolato localmente, poi MPI_Allreduce prende il min globale.
// ============================================================
template <int dim>
void WaveEquation<dim>::check_cfl_condition() const
{
    if (time_scheme == TimeScheme::NEWMARK)
    {
        pcout << "  Schema Newmark-β: incondizionatamente stabile.\n";
        return;
    }

    double local_h_min = std::numeric_limits<double>::max();
    for (const auto &cell : dof_handler.active_cell_iterators())
        if (cell->is_locally_owned())
            local_h_min = std::min(local_h_min, cell->minimum_vertex_distance());

    // MPI_Allreduce: prende il minimo globale
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
        pcout << "  *** ATTENZIONE: CFL VIOLATA! ***";
    else
        pcout << "  [OK]";
    pcout << "\n";
}

// ============================================================
// output_results  —  output VTU parallelo
//
// Strategia deal.II MPI:
//   - Ogni processo scrive il suo file  solution-NNNN-RRRR.vtu
//     dove RRRR = rank MPI
//   - Il processo 0 scrive il file master  solution-NNNN.pvtu
//     che contiene la lista di tutti i file .vtu
//   - ParaView carica il .pvtu e ricostruisce automaticamente
//     la soluzione completa
// ============================================================
template <int dim>
void WaveEquation<dim>::output_results(unsigned int step)
{
    TimerOutput::Scope t(computing_timer, "output");

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);

    // I vettori "with ghost" sono già sincronizzati
    data_out.add_data_vector(solution_u, "displacement");
    data_out.add_data_vector(velocity_u, "velocity");

    // Aggiunge il rank MPI come campo (utile per verificare il partizionamento)
    Vector<float> proc_id(triangulation.n_active_cells());
    proc_id = static_cast<float>(this_mpi_proc);
    data_out.add_data_vector(proc_id, "mpi_rank");

    data_out.build_patches();

    // Nome base del file
    const std::string base = "solution-" + Utilities::int_to_string(step, 4);

    // Ogni processo scrive il suo pezzo
    std::ofstream local_out(base + "-" +
        Utilities::int_to_string(this_mpi_proc, 4) + ".vtu");
    data_out.write_vtu(local_out);

    // Solo il rank 0 scrive il file .pvtu master
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
// Chiamata alla fine di run() — stampa statistiche globali
// su stdout (solo rank 0) e scrive energy_report.txt
// ============================================================
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

// ============================================================
// refine_mesh  —  AMR parallelo con p4est
//
// parallel::distributed::GridRefinement distribuisce i flag
// tra i processi e p4est gestisce il ribilanciamento del carico.
// parallel::distributed::SolutionTransfer trasferisce la
// soluzione sulla nuova mesh distribuita.
// ============================================================
template <int dim>
void WaveEquation<dim>::refine_mesh()
{
    TimerOutput::Scope t(computing_timer, "refine_mesh");

    pcout << "  AMR parallelo...\n";

    // Stima errore locale (solo celle locally owned)
    Vector<float> err(triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate(
        dof_handler,
        QGauss<dim - 1>(fe_ptr->degree + 1),
        std::map<types::boundary_id, const Function<dim> *>(),
        solution_u, err);

    // parallel::distributed::GridRefinement: coordina i flag tra processi
    parallel::distributed::GridRefinement
        ::refine_and_coarsen_fixed_fraction(
            triangulation, err,
            amr_refine_fraction, amr_coarsen_fraction);

    // Limita il livello massimo
    for (auto &cell : triangulation.active_cell_iterators())
        if (cell->is_locally_owned() &&
            cell->level() >= static_cast<int>(max_refinement_level))
            cell->clear_refine_flag();

    triangulation.prepare_coarsening_and_refinement();

    // ---- SolutionTransfer distribuito ----
    // Nota: parallel::distributed::SolutionTransfer vuole
    // un solo vettore per chiamata. Utilizziamo tre trasferimenti
    // separati oppure un vettore concatenato (usiamo tre chiamate).
    parallel::distributed::SolutionTransfer<dim, TrilinosVector>
        st_u(dof_handler), st_u_old(dof_handler), st_v(dof_handler);

    solution_u.update_ghost_values();
    solution_u_old.update_ghost_values();
    velocity_u.update_ghost_values();

    st_u.prepare_for_coarsening_and_refinement(solution_u);
    st_u_old.prepare_for_coarsening_and_refinement(solution_u_old);
    st_v.prepare_for_coarsening_and_refinement(velocity_u);

    triangulation.execute_coarsening_and_refinement();

    // Ridistribuisce DoF sulla nuova mesh
    dof_handler.distribute_dofs(*fe_ptr);
    locally_owned_dofs    = dof_handler.locally_owned_dofs();
    locally_relevant_dofs = DoFTools::extract_locally_relevant_dofs(dof_handler);

    pcout << "  AMR: DoF globali=" << dof_handler.n_dofs()
          << "  celle=" << triangulation.n_global_active_cells() << "\n";

    // Ricostruisce constraints
    constraints.clear();
    constraints.reinit(locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    if (!use_absorbing_bc)
        VectorTools::interpolate_boundary_values(
            dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
    constraints.close();

    // Ricostruisce sparsity e matrici
    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern(
        dsp, locally_owned_dofs, mpi_comm, locally_relevant_dofs);

    laplace_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_comm);
    if (use_absorbing_bc)
        boundary_mass_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_comm);
    if (time_scheme == TimeScheme::NEWMARK)
        system_matrix_newmark.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_comm);

    // Realloca vettori sulla nuova partizione
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

    // Interpola le soluzioni sulla nuova mesh
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

    // Riassembla matrici sulla nuova mesh
    assemble_matrices();
    check_cfl_condition();
}

// ============================================================
// run  —  loop temporale principale
// ============================================================
template <int dim>
void WaveEquation<dim>::run()
{
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

    // 1. Griglia e sistema
    if (mode == SimulationMode::DIFFRACTION)
        make_grid_with_obstacle();
    else
        make_grid();

    setup_system();
    assemble_matrices();
    check_cfl_condition();

    // 2. Log energia (solo rank 0)
    if (track_energy && this_mpi_proc == 0)
    {
        energy_log.open("energy_log.csv");
        energy_log << "step,time,kinetic,potential,total,drift_rel,Linfty\n";
    }

    // 3. Condizioni iniziali
    pcout << "  Condizioni iniziali...\n";

    if (mode == SimulationMode::MMS_CONVERGENCE)
    {
        // IC da soluzione esatta
        InitialDisplacementMMS<dim> u0;
        InitialVelocityMMS<dim>     u1;

        VectorTools::interpolate(dof_handler, u0, owned_solution_u);
        constraints.distribute(owned_solution_u);
        solution_u = owned_solution_u;

        VectorTools::interpolate(dof_handler, u1, owned_velocity_u);
        // u_old = u - dt*v  (encoda velocità iniziale nel Leapfrog)
        for (const auto idx : locally_owned_dofs)
            owned_solution_u_old(idx) = owned_solution_u(idx)
                                      - time_step * owned_velocity_u(idx);
        solution_u_old = owned_solution_u_old;
        velocity_u     = owned_velocity_u;
    }
    else if (mode == SimulationMode::INTERFERENCE)
    {
        // Due sorgenti gaussiane
        const double amp   = 1.0;
        const double width = 0.05;
        const Point<dim> s1(0.3, 0.5), s2(0.7, 0.5);

        // Usiamo interpolate con una funzione lambda tramite FunctionFromFunctionObjects
        // (più semplice: scorriamo i support points locali)
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
        // Gaussiana singola centrata a sx (funziona per tutti gli altri modi)
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

    // Lista .pvtu per il file .pvd master (solo rank 0)
    std::vector<std::pair<double, std::string>> pvd_list;
    if (this_mpi_proc == 0)
        pvd_list.push_back({0.0, "solution-0000.pvtu"});

    // 4. Loop temporale
    time        = 0.0;
    step_number = 0;

    while (time < end_time - 1e-12)
    {
        step_number++;
        time += time_step;

        // AMR parallelo (non per MMS)
        if (use_amr && mode != SimulationMode::MMS_CONVERGENCE
                    && step_number % amr_every_n_steps == 0)
            refine_mesh();

        // Avanzamento temporale
        if (time_scheme == TimeScheme::LEAPFROG)
            solve_time_step();
        else
            solve_time_step_newmark();

        // Output e statistiche
        if (step_number % output_every_n_steps == 0)
        {
            pcout << "  Step " << step_number << "  t=" << time;

            if (track_energy)
            {
                const double Ek    = compute_kinetic_energy();
                const double Ep    = compute_potential_energy();
                const double Etot  = Ek + Ep;
                const double Linf  = compute_Linfty_norm();

                // Inizializza E0 al primo step di output
                if (!energy_initialized)
                {
                    energy_initial     = Etot;
                    energy_initialized = true;
                }

                // Drift relativo: misura la conservazione dell'energia
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

                // Rilevamento blow-up: se |u|_inf > 100 * ampiezza iniziale
                // la simulazione è probabilmente instabile (CFL violata)
                if (Linf > 100.0)
                    pcout << "\n  *** ATTENZIONE: possibile blow-up numerico! "
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

    pcout << "  Simulazione completata. Steps: " << step_number << "\n";
    computing_timer.print_summary();
}

// ============================================================
// run_convergence_study  —  MMS su livelli di raffinamento
// ============================================================
template <int dim>
void WaveEquation<dim>::run_convergence_study()
{
    pcout << "\n=== Studio di Convergenza MMS ("
          << n_mpi_procs << " processi MPI, Q"
          << fe_degree << ", "
          << (time_scheme == TimeScheme::LEAPFROG ? "Leapfrog" : "Newmark-β")
          << ") ===\n";

    mode    = SimulationMode::MMS_CONVERGENCE;
    use_amr = false;
    end_time = 1.0;

    const std::vector<unsigned int> levels = {3, 4, 5, 6, 7};

    for (unsigned int ref : levels)
    {
        pcout << "\n--- Livello " << ref << " ---\n";

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

        // Loop senza output
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

        convergence_table.add_value("Livello", ref);
        convergence_table.add_value("Celle",   (unsigned int)triangulation.n_global_active_cells());
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

    // Solo rank 0 stampa e salva la tabella
    if (this_mpi_proc == 0)
    {
        pcout << "\n=== Tabella di Convergenza ===\n";
        convergence_table.write_text(std::cout);
        std::ofstream f("convergence_table.txt");
        convergence_table.write_text(f);
        pcout << "Tabella salvata in convergence_table.txt\n";
    }

    computing_timer.print_summary();
}

// ============================================================
// Template instantiation per 2D e 3D
// ============================================================
template class WaveEquation<2>;
template class WaveEquation<3>;
