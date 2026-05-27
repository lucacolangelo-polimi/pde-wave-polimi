# WaveEquation — Solver FEM MPI-parallelo

## Dipendenze

| Libreria | Versione | Ruolo |
|----------|----------|-------|
| deal.II | ≥ 9.4 | framework FEM |
| MPI | qualsiasi | comunicazione inter-processo |
| Trilinos | ≥ 13.0 | matrici/vettori distribuiti |
| p4est | ≥ 2.3 | mesh distribuita |

## Compilazione

```bash
mkdir build && cd build
cmake -DDEAL_II_DIR=/path/to/dealii ..
make -j$(nproc)
```

## Esecuzione

```bash
# 4 processi MPI
mpirun -np 4 ./wave_equation

# 8 processi MPI (consigliato per mesh fine, livello 7+)
mpirun -np 8 ./wave_equation

# Studio di convergenza MMS (decommenta in main.cpp)
mpirun -np 4 ./wave_equation
```

## Architettura MPI

```
parallel::distributed::Triangulation  →  mesh partizionata (p4est)
         ↓
DoFHandler::distribute_dofs()         →  DoF assegnati a ogni processo
         ↓
IndexSet: locally_owned_dofs          →  DoF "posseduti" da questo rank
IndexSet: locally_relevant_dofs       →  owned + ghost (vicini)
         ↓
TrilinosWrappers::SparseMatrix        →  righe distribuite tra processi
TrilinosWrappers::MPI::Vector         →  vettori distribuiti
         ↓
assemble_matrices()  [owned cells]    →  compress(add) → MPI_Allreduce
         ↓
solve_time_step() Leapfrog:
  a = M⁻¹·RHS   [locale, no MPI]     →  u_new = 2u - u_old + dt²·a
         ↓ oppure
solve_time_step_newmark():
  SolverCG Trilinos                   →  MPI_Allreduce interni al solver
         ↓
output_results():
  ogni rank → solution-NNNN-RRRR.vtu
  rank 0    → solution-NNNN.pvtu + solution.pvd
```

## Vettori: owned vs ghost

| Tipo | IndexSet | Uso |
|------|----------|-----|
| `owned_solution_u` | `locally_owned_dofs` | scrittura (assembly, update) |
| `solution_u` | `locally_relevant_dofs` | lettura (vmult, output) |

Dopo ogni aggiornamento owned, chiama `solution_u.update_ghost_values()`
per propagare i valori ai processi vicini.

## Visualizzazione con ParaView

```bash
# Apri il file master PVD (contiene tutti i timestep)
paraview solution.pvd

# Oppure apri un singolo PVTU
paraview solution-0010.pvtu
```

Il campo `mpi_rank` nell'output mostra la partizione della mesh tra i processi.
