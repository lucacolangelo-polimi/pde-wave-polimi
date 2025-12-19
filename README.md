# Solving the Wave Equation

General formula:

$$
\frac{\partial^2 u}{\partial t^2} - \Delta u = f \quad \text{in } \Omega \times (0,T]
$$

$$
u = g \quad \text{on } \partial\Omega \times [0,T]
$$

$$
u(x,0) = u_0(x) \quad \text{in } \Omega
$$

$$
\frac{\partial u}{\partial t}(x,0) = u_1(x) \quad \text{in } \Omega
$$

Where:

- $u(x,t)$ is the unknown function (the wave field).
- $\Omega$ is the spatial domain.
- $\partial\Omega$ is the boundary of the domain.
- The Laplace operator is defined as:

$$
\Delta u = \frac{\partial^2 u}{\partial x_1^2} + \frac{\partial^2 u}{\partial x_2^2}
$$

# 1. Term-by-term analysis

### a. Partial Differential Equation (PDE)

$$
\frac{\partial^2 u}{\partial t^2} - \Delta u = f \quad \text{in } \Omega \times (0,T]
$$

### Hyperbolic term

$\frac{\partial^2 u}{\partial t^2}$  
Represents the acceleration of the field over time. This is the key term that makes the equation hyperbolic and describes wave propagation.

### Spatial term

$\Delta u$  
Represents restoring forces or internal stresses within the medium. It is the source of spatial propagation of the motion.

### Forcing term

$f$  
Represents an external energy source or applied force acting on the system within the domain. If $f = 0$, the equation is homogeneous and describes freely propagating waves.

**Note:** The standard form of the wave equation includes a propagation constant $c^2$:

$$
\frac{\partial^2 u}{\partial t^2} - c^2 \Delta u = f
$$

In our problem, it is implicitly assumed that $c = 1$ (unit propagation speed).

---

### b. Boundary Condition

$$
u = g \quad \text{on } \partial\Omega \times [0,T]
$$

This is a **Dirichlet boundary condition**.  
It specifies the value of the wave field $u$ along the boundary $\partial\Omega$ of the domain for all times. $g$ is a known function.  
If $g = 0$, the boundary is “fixed,” or—depending on interpretation—an “ideal absorbing boundary.”

---

### c. Initial Conditions

Because the wave equation is second order in time, it requires **two initial conditions** to define the solution uniquely:

- $u_0$: the initial configuration (initial displacement) of the wave field at time $t = 0$.  
- $u_1$: the initial velocity of the wave field at time $t = 0$.

# 2. Weak Form

For the homogeneous equation (with $f = 0$ and $g = 0$), multiplying by a test function $v \in H_0^1(\Omega)$ and integrating over $\Omega$, we get:

$$
\int_\Omega \frac{\partial^2 u}{\partial t^2} v \, dx - \int_\Omega \Delta u \, v \, dx = 0
$$

Applying integration by parts (Green's theorem) to the Laplacian term:

$$
\int_\Omega \frac{\partial^2 u}{\partial t^2} v \, dx + \int_\Omega \nabla u \cdot \nabla v \, dx = 0
$$

The FEM implementation is based on the discretization of this weak form.

### Implications for the Numerical Solution

The hyperbolic nature of the equation is crucial for our implementation: numerical stability relies on the correct handling of the time-derivative term in conjunction with the spatial discretization.

---

# 3. Time Discretization

Since the problem is time-dependent, we use a time integration method. We chose the **Centered Finite Difference Method** (also known as the **Leapfrog scheme**).

The choice of the time-stepping method strongly affects **numerical stability** (CFL condition).

---

# 4. Space Discretization

This is where the Finite Element Method (FEM) comes into play. To use FEM, the problem is reformulated in its weak form within an appropriate functional space (typically $H^1(\Omega)$) and discretized using the **deal.II library**.

---

# 5. General Structure of the Project

The goal of this project is to solve the 2D wave equation using the **deal.II library** for the Finite Element Method (FEM) in space and an **explicit time-stepping scheme**. The solver follows these main steps:

## a. Mesh Generation and DoF Management
- The domain $\Omega$ is triangulated using deal.II's `Triangulation` class.
- We generate a structured grid (e.g., `GridGenerator::hyper_cube`) or import it from external formats (e.g., Gmsh).
- Degrees of Freedom (DoFs) are distributed using a **P1 finite element space** (`FE_SimplexP<dim>(1)`), representing linear polynomials on each triangle.

## b. System Assembly
Using deal.II's `FEValues` and quadrature formulas (`QGauss`), we compute the system matrices:
- **Mass Matrix ($M$)**: Represents the $L^2$ inner product of basis functions $\int \phi_i \phi_j$.
- **Stiffness Matrix ($K$)**: Represents the inner product of gradients $\int \nabla \phi_i \cdot \nabla \phi_j$.

## c. Boundary Conditions
- **Dirichlet boundaries** (where $u=g$) are handled by identifying boundary DoFs and applying constraints to the solution vectors or the system matrix.

## d. Time-Stepping Scheme (Leapfrog)
We implement the explicit central difference scheme:

$$
U^{n+1} = 2 U^n - U^{n-1} + \Delta t^2 M^{-1} (F^n - K U^n)
$$

In the deal.II implementation, solving $M^{-1}$ typically involves solving a linear system ($M a = RHS$) at each step using a Solver (like Conjugate Gradient) or using Mass Lumping techniques.

## e. Output and Analysis
- The solution is exported at regular time intervals using `DataOut` in **VTK format**.
- Results can be visualized in **Paraview** to observe wave propagation, reflection, and interference.

---

# 6. Code Architecture with deal.II

The project leverages the **deal.II** library structure. Instead of separating Mesh, Assembly, and Solver into different classes, the logic is encapsulated within a single template class that manages the entire simulation lifecycle.

### The `WaveEquation<dim>` Class
This is the core class of the project. It relies on C++ templates to be dimension-independent (working seamlessly in 2D or 3D).

#### 1. Setup and Grid Management (`make_grid`, `setup_system`)
- **`make_grid()`**: Generates the geometry using `GridGenerator` or reads a mesh file using `GridIn`. It handles global refinement.
- **`setup_system()`**: Initializes the `DoFHandler`, computes the sparsity pattern, and resizes the sparse matrices ($M$, $K$) and solution vectors ($U_{old}, U_{curr}, U_{new}$) to the correct size.

#### 2. Assembly (`assemble_system`)
- This method populates the **Mass Matrix** and **Stiffness Matrix**.
- It iterates over all active cells, initializes `FEValues`, computes local contributions via quadrature, and distributes them into the global sparse matrices.
- This is done **once** before the time loop (since the mesh is static).

#### 3. Time Evolution (`run`, `solve_time_step`)
- **`run()`**: The main driver. It orchestrates the setup, applies initial conditions, and executes the time loop.
- **`solve_time_step()`**: Implements the algebraic operations for the Leapfrog scheme. It computes the Right-Hand Side (RHS), solves the linear system involving the Mass Matrix (to find acceleration), and updates the displacement vector $U^{n+1}$.

#### 4. Output (`output_results`)
- Uses the `DataOut` class to attach the DoF handler and the solution vector.
- Writes `.vtk` files (e.g., `solution-001.vtk`) containing the wave field state at specific time steps for post-processing.

### Main Entry Point (`main.cc`)
A minimal file that initializes the `WaveEquation<2>` object and calls the `run()` method within a try-catch block to handle deal.II exceptions.

### Build System

Dependencies and build are handled through a Docker container made available from the course instructors at https://github.com/HPC-Courses/AMSC-Labs/blob/main/Labs/2025-26/00-environment_setup/README.md#2-linux-users.

Github will be used to keep track of the project files while the container will be handled locally be each partecipant of the project.
