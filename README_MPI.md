# WaveEquation — MPI-parallel FEM solver

## Dependencies

| Library | Version | Role |
|----------|----------|-------|
| deal.II | ≥ 9.4 | FEM framework |
| MPI | any | inter-process communication |
| Trilinos | ≥ 13.0 | distributed matrices/vectors |
| p4est | ≥ 2.3 | distributed mesh |

## Compilation

```bash
mkdir build && cd build
cmake -DDEAL_II_DIR=/path/to/dealii ..
make -j$(nproc)

# 4 MPI processes
mpirun -np 4 ./wave_equation

# 8 MPI processes (recommended for fine meshes, level 7+)
mpirun -np 8 ./wave_equation

# MMS convergence study (uncomment in main.cpp)
mpirun -np 4 ./wave_equation


parallel::distributed::Triangulation  →  partitioned mesh (p4est)
         ↓
DoFHandler::distribute_dofs()         →  DoFs assigned to each process
         ↓
IndexSet: locally_owned_dofs          →  DoFs "owned" by this rank
IndexSet: locally_relevant_dofs       →  owned + ghost (neighbors)
         ↓
TrilinosWrappers::SparseMatrix        →  rows distributed among processes
TrilinosWrappers::MPI::Vector         →  distributed vectors
         ↓
assemble_matrices()  [owned cells]    →  compress(add) → MPI_Allreduce
         ↓
solve_time_step() Leapfrog:
  a = M⁻¹·RHS   [local, no MPI]      →  u_new = 2u - u_old + dt²·a
         ↓ or
solve_time_step_newmark():
  Trilinos SolverCG                   →  internal solver MPI_Allreduce
         ↓
output_results():
  each rank → solution-NNNN-RRRR.vtu
  rank 0    → solution-NNNN.pvtu + solution.pvd

  # Open the master PVD file (contains all timesteps)
paraview solution.pvd

# Or open a single PVTU file
paraview solution-0010.pvtu