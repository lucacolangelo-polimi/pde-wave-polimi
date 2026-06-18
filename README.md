# Solving the Wave Equation — MPI Parallel FEM Solver

A distributed-memory finite element solver for the 2D/3D scalar wave equation,
developed as a final project for the Numerical Methods for PDEs course at
Politecnico di Milano.

## General formulation

$$
\frac{\partial^2 u}{\partial t^2} - \nabla \cdot (c^2(\mathbf{x}) \nabla u) = f \quad \text{in } \Omega \times (0,T]
$$

$$
u = 0 \quad \text{on } \partial\Omega \times (0,T]
$$

$$
u(\mathbf{x},0) = u_0(\mathbf{x}) \quad \text{in } \Omega
$$

$$
\frac{\partial u}{\partial t}(\mathbf{x},0) = u_1(\mathbf{x}) \quad \text{in } \Omega
$$

Where:
- $u(\mathbf{x},t)$ is the displacement field
- $c(\mathbf{x}) > 0$ is the local wave speed (possibly heterogeneous)
- $\Omega = [0,1]^d$ with $d = 2, 3$
- Homogeneous Dirichlet boundary conditions are imposed on $\partial\Omega$

# 1. Term-by-term analysis

### a. Partial Differential Equation (PDE)

$$
\frac{\partial^2 u}{\partial t^2} - \nabla \cdot (c^2(\mathbf{x}) \nabla u) = f
$$

**Hyperbolic term** — $\frac{\partial^2 u}{\partial t^2}$: acceleration of the field, makes the equation hyperbolic.

**Spatial term** — $\nabla \cdot (c^2(\mathbf{x}) \nabla u)$: restoring forces within the medium. The spatially varying wave speed $c(\mathbf{x})$ allows modeling heterogeneous media (refraction).

**Forcing term** — $f$: external energy source. If $f = 0$ the equation describes freely propagating waves.

### b. Boundary Condition

Homogeneous Dirichlet conditions $u = 0$ on $\partial\Omega$ model rigid boundaries.

### c. Initial Conditions

Two initial conditions are required (second-order in time):
- $u_0(\mathbf{x})$: initial displacement
- $u_1(\mathbf{x})$: initial velocity

# 2. Weak Form

Multiplying by a test function $v \in H_0^1(\Omega)$ and integrating over $\Omega$,
applying Green's first identity to eliminate the divergence term:

$$\int_\Omega \frac{\partial^2 u}{\partial t^2} v \, d\mathbf{x} + \int_\Omega c^2(\mathbf{x}) \nabla u \cdot \nabla v \, d\mathbf{x} = \int_\Omega f v \, d\mathbf{x} \quad \forall v \in H_0^1(\Omega)$$
The boundary integral vanishes because $v = 0$ on $\partial\Omega$.

# 3. Space Discretization — Galerkin FEM

The domain is discretized with $Q_1$ isoparametric elements (`FE_Q<dim>(1)`)
on a quadrilateral mesh. The discrete solution

$$u_h(\mathbf{x},t) = \sum_{j=1}^{N_h} U_j(t)\,\varphi_j(\mathbf{x})$$

leads to the semi-discrete system:

$$\mathbf{M}\,\ddot{\mathbf{U}}(t) + \mathbf{K}\,\mathbf{U}(t) = \mathbf{F}(t)$$

where:
- **$M$** — lumped mass matrix (assembled via Gauss-Lobatto quadrature, diagonal)
- **$K$** — stiffness matrix with $c^2(\mathbf{x})$ sampled at quadrature points
- **$F$** — load vector

Mass lumping via Gauss-Lobatto quadrature makes $M$ diagonal, enabling
trivial inversion $O(N_h)$ with no linear solver.

# 4. Time Discretization

Two schemes are implemented:

### Leapfrog (explicit, symplectic)

$$\mathbf{U}^{n+1} = 2\mathbf{U}^n - \mathbf{U}^{n-1} + \Delta t^2 (\mathbf{M}^L)^{-1}(\mathbf{F}^n - \mathbf{K}\mathbf{U}^n)$$

No linear solver required. Subject to the CFL condition:

$$\Delta t \leq \frac{h_{\min}}{c_{\max}\sqrt{d}}$$

Used for all physical simulations.

### Newmark-β (implicit, unconditionally stable)

$$(\mathbf{M}^L + \beta\Delta t^2 \mathbf{K})\,\ddot{\mathbf{U}}^{n+1} = \mathbf{F}^{n+1} - \mathbf{K}\tilde{\mathbf{U}}^{n+1}$$

with $\beta = 1/4$, $\gamma = 1/2$ (trapezoidal rule). Solved via CG + ILU preconditioner
(Trilinos). Used exclusively for the MMS convergence study, where very small
$\Delta t = 0.1h^2$ would violate the CFL condition.

# 5. MPI Parallelization

The solver is fully distributed via MPI using deal.II's
`parallel::distributed::Triangulation` and Trilinos for linear algebra.

Key design choices:
- **`locally_owned_dofs`**: DoFs written by each process
- **`locally_relevant_dofs`**: owned + ghost DoFs (needed for assembly)
- **Two constraint objects**: `matrix_constraints` (hanging nodes only,
  for matrix assembly) and `constraints` (hanging nodes + Dirichlet BC,
  for vectors)
- **AMR**: `KellyErrorEstimator` + `SolutionTransfer` on three vectors
  ($U^n$, $U^{n-1}$, $\dot{U}^n$)

# 6. Code Architecture

### `WaveEquation<dim>` class

Single templated class (`dim` ∈ {2,3}) encapsulating the full pipeline.

| Method | Description |
|--------|-------------|
| `make_grid()` | Generates $[0,1]^d$ mesh via `GridGenerator::hyper_cube` |
| `make_grid_with_obstacle()` | Mesh with internal wall for diffraction |
| `setup_system()` | Distributes DoFs, builds constraints and sparsity pattern |
| `assemble_matrices()` | Assembles $K$ and $M^L$ on locally owned cells |
| `solve_time_step()` | Leapfrog explicit update |
| `solve_time_step_newmark()` | Newmark predictor-corrector update |
| `refine_mesh()` | AMR with p4est + SolutionTransfer |
| `run()` | Main time loop with energy tracking and VTK output |
| `run_convergence_study()` | MMS verification |
| `run_scaling_benchmark()` | Strong/weak scaling benchmarks |

### Simulation modes

| Mode | Description |
|------|-------------|
| `mms` | MMS convergence study (Newmark, single process) |
| `pebble` | Gaussian pulse propagation (Leapfrog + AMR) |
| `interference` | Two-source interference pattern (Leapfrog + AMR) |
| `refraction` | Heterogeneous wave speed $c(\mathbf{x})$ (Leapfrog + AMR) |
| `diffraction` | Rigid wall with narrow slit (Leapfrog + AMR) |
| `damping` | Viscous sponge layer (Leapfrog + AMR) |
| `absorbing` | Sommerfeld absorbing BC (Leapfrog + AMR) |
| `strong` | Strong scaling benchmark (Newmark) |
| `weak` | Weak scaling benchmark (Newmark) |

# 7. Dependencies

| Library | Version | Role |
|---------|---------|------|
| deal.II | ≥ 9.4 | FEM framework |
| MPI | any | inter-process communication |
| Trilinos | ≥ 13.0 | distributed matrices/vectors |
| p4est | ≥ 2.3 | distributed mesh |

# 8. Compilation and Running

```bash
mkdir build && cd build
cmake -DDEAL_II_DIR=/path/to/dealii ..
make -j$(nproc)

# MMS convergence study
mpirun -n 1 ./wave_equation --mode mms

# Physical scenarios
mpirun -n 4 ./wave_equation --mode pebble --ref 6

# Scaling benchmarks
mpirun -n 1 ./wave_equation --mode strong --ref 8
mpirun -n 4 ./wave_equation --mode strong --ref 8
```

Output: `.pvtu` + `.pvd` files for ParaView, `energy_log.csv` for energy diagnostics.

# 9. Authors

Barrella Teresa, Chiari Stefano, Colangelo Luca, Sanjevic Milica  
Politecnico di Milano — June 2026
