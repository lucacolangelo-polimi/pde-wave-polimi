# WaveEquation — MPI-parallel FEM Solver

A distributed-memory finite element solver for the 2D/3D wave equation,
built on deal.II with Trilinos and p4est.

## Dependencies

| Library | Version | Role |
|---------|---------|------|
| deal.II | ≥ 9.4   | FEM framework |
| MPI     | any     | inter-process communication |
| Trilinos| ≥ 13.0  | distributed matrices/vectors |
| p4est   | ≥ 2.3   | distributed mesh |

## Compilation

```bash
mkdir build && cd build
cmake -DDEAL_II_DIR=/path/to/dealii ..
make -j$(nproc)
```

## Running

```bash
# MMS convergence study (single process recommended)
mpirun -n 1 ./wave_equation --mode mms

# Physical scenarios
mpirun -n 4 ./wave_equation --mode pebble      --ref 6
mpirun -n 4 ./wave_equation --mode interference --ref 6
mpirun -n 4 ./wave_equation --mode refraction   --ref 6
mpirun -n 4 ./wave_equation --mode diffraction  --ref 6
mpirun -n 4 ./wave_equation --mode damping      --ref 6
mpirun -n 4 ./wave_equation --mode absorbing    --ref 6

# Scaling benchmarks
mpirun -n 1 ./wave_equation --mode strong --ref 8
mpirun -n 2 ./wave_equation --mode strong --ref 8
mpirun -n 4 ./wave_equation --mode strong --ref 8
```

## Available modes

| Mode | Description |
|------|-------------|
| `mms` | MMS convergence study (Newmark, single process) |
| `pebble` | Gaussian pulse, isotropic propagation (Leapfrog + AMR) |
| `interference` | Two Gaussian sources, interference pattern (Leapfrog + AMR) |
| `refraction` | Heterogeneous wave speed c(x) (Leapfrog + AMR) |
| `diffraction` | Rigid wall with narrow slit (Leapfrog + AMR) |
| `damping` | Viscous sponge layer (Leapfrog + AMR) |
| `absorbing` | Sommerfeld absorbing BC (Leapfrog + AMR) |
| `strong` | Strong scaling benchmark (Newmark, AMR off) |
| `weak` | Weak scaling benchmark (Newmark, AMR off) |

## MPI Architecture

```
parallel::distributed::Triangulation  →  partitioned mesh (p4est)
         ↓
DoFHandler::distribute_dofs()         →  DoFs assigned to each process
         ↓
IndexSet: locally_owned_dofs          →  DoFs written by this rank
IndexSet: locally_relevant_dofs       →  owned + ghost (neighbors)
         ↓
TrilinosWrappers::SparseMatrix        →  rows distributed among processes
TrilinosWrappers::MPI::Vector         →  distributed vectors
         ↓
assemble_matrices()  [owned cells]    →  compress(add) → MPI_Allreduce
         ↓
solve_time_step() — Leapfrog:
  a = (M^L)⁻¹·RHS  [local, no MPI]  →  u_new = 2u - u_old + dt²·a
         ↓  or
solve_time_step_newmark() — Newmark:
  Trilinos CG + ILU                  →  internal MPI_Allreduce
         ↓
output_results():
  each rank → solution-NNNN-RRRR.vtu
  rank 0    → solution-NNNN.pvtu + solution.pvd
```

## Visualization (ParaView)

```bash
# Open the master PVD file (all timesteps as animation)
paraview solution.pvd

# Or open a single timestep
paraview solution-0010.pvtu
```

## Output files

| File | Description |
|------|-------------|
| `solution-NNNN.pvtu` | Master file for timestep NNNN |
| `solution-NNNN-RRRR.vtu` | Local data for rank RRRR |
| `solution.pvd` | Time-series index for ParaView |
| `energy_log.csv` | Energy diagnostics (Ek, Ep, Etot, drift) |
| `energy_report.txt` | Summary report |