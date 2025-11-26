# Solving the wave equation

General formula:

$$
\frac{\partial^2 u}{\partial t^2} - \Delta u = f \quad \text{in } \Omega \times (0,T]
$$

$$
u = g \quad \text{su } \partial\Omega \times [0,T]
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
\Delta u =
\frac{\partial^2 u}{\partial x_1^2}
+
\frac{\partial^2 u}{\partial x_2^2}.
$$

# 1. Term-by-term analysis

### a. Partial Differential Equation (PDE)

$$
\frac{\partial^2 u}{\partial t^2} - \Delta u = f 
\quad \text{in } \Omega \times (0,T]
$$

### Hyperbolic term

`∂²u/∂t²`  
Represents the acceleration of the field over time. This is the key term that makes the equation hyperbolic and describes wave propagation.

### Spatial term

`Δu`  
Represents restoring forces or internal stresses within the medium. It is the source of spatial propagation of the motion.

### Forcing term

`f`  
Represents an external energy source or applied force acting on the system within the domain. If `f = 0`, the equation is homogeneous and describes freely propagating waves.

**Note:**  
The standard form of the wave equation includes a propagation constant `c²`:

`∂²u/∂t² - c² Δu = f`

In our problem, it is implicitly assumed that `c = 1` (unit propagation speed).


---

### b. Boundary Condition

$$
u = g \quad \text{on } \partial\Omega \times [0,T]
$$

This is a **Dirichlet boundary condition**.  
It specifies the value of the wave field $u$ along the boundary $\partial\Omega$ of the domain for all times.   $g$ is a known function.  
If $g = 0$, the boundary is “fixed,” or—depending on interpretation—an “ideal absorbing boundary.”

---

### c. Initial Conditions

Because the wave equation is second order in time, it requires **two initial conditions** to define the solution uniquely:

- $u_0$: the initial configuration (initial displacement) of the wave field at time $t = 0$.  
- $u_1$: the initial velocity of the wave field at time $t = 0$.

# 2. Weak Form

For the homogeneous equation (with `f = 0` and `g = 0`), multiplying by a test function
`v \in H_0^1(\Omega)` and integrating over `Ω`, we get:

`∫_Ω (∂²u/∂t²) v dx - ∫_Ω Δu v dx = 0`

Applying integration by parts (Green's theorem) to the Laplacian term:

`∫_Ω (∂²u/∂t²) v dx + ∫_Ω ∇u ⋅ ∇v dx = 0`

The FEM implementation will be based on the discretization of this weak form.

### Implications for the Numerical Solution

The hyperbolic nature of the equation is crucial for our implementation:

---

# 3. Time Discretization

Since the problem is time-dependent (evolution in time), you will need to use a time integration method such as:

- the **Centered Finite Difference Method** (or **Leapfrog scheme**), which is classically used for the wave equation,  
- or more advanced methods such as **Runge–Kutta** schemes.

The choice of the time-stepping method will strongly affect **numerical stability**.

---

# 4. Space Discretization

This is where the Finite Element Method (FEM), which is required, comes into play.  
To use FEM, it is needed to reformulate the problem in its weak (or variational) form within an appropriate functional space (typically $H^1(\Omega)$).

---

# 5. General structure of the project
## 2D Wave Equation Solver

The goal of this project is to solve the 2D wave equation using the **Finite Element Method (FEM)** in space and an **explicit time-stepping scheme**. The solver follows these main steps:

## a. Mesh Generation and FEM Space

- The domain Ω is divided into triangles (2D mesh).  
- We define a **P1 function space**, i.e., linear polynomials on each triangle.  
  This allows us to represent the solution \(u(x, y, t)\) in a simple and linear way between mesh nodes.

## b. Assembly of Matrices

- **Mass matrix M**: represents the accumulation of “mass” at each node.  
- **Stiffness matrix K**: represents the “resistance” of the system to changes in the solution (spatial derivatives).  

For an explicit scheme, it is convenient to use **mass-lumping**: replace M with a diagonal approximation.  
This makes each node independent in the time update, avoiding the need to solve large linear systems.

## c. Application of Boundary Conditions

- The solution values on the domain boundaries (∂Ω) are fixed to \(g\).  
- These conditions are applied directly to the matrices or to the initial solution, depending on the implementation.

## d. Time-Stepping Scheme

We use a **centered explicit scheme** (central difference / Leap-Frog):

\[
U^{n+1} = 2 U^n - U^{n-1} + \Delta t^2 M^{-1} (F^n - K U^n)
\]

**Advantages:**
- Simple to implement.
- Stable if Δt is small enough (CFL condition).
- Very efficient if M is diagonal.

**Disadvantages:**
- Introduces slight numerical dispersion (waves may travel slightly faster or slower than the exact solution).
- Minimal numerical dissipation.

## e. Time Loop

- Starting from the initial conditions \(u_0\) and \(\partial_t u|_{t=0} = u_1\), the solution is updated step by step using the explicit scheme.  
- At each time step, the solution can be saved or visualized to analyze wave propagation.

## f. Output and Analysis

- Visualize the evolution of \(u(x, y, t)\) over time.  
- Discuss the effects of time step Δt, mesh resolution, mass-lumping, and numerical dispersion.

# 6. Code Architecture
## Project Architecture

The solver is organized in a modular way, with each module handling a specific part of the computation. The main modules are:

### A. Mesh Module
This module defines the 2D mesh and its properties:
- **Nodes, elements, and boundary nodes**  
- **Element areas and Jacobians**  

For small meshes, the mesh can be defined manually, or optionally read from a `.msh` file.

**Files:**
- `Mesh.hpp`: defines `Node` and `Element` structures and the `Mesh` class, which contains vectors of nodes, elements, and boundary nodes.
- `Mesh.cpp`: implements the constructor, computes element areas and Jacobians, and optionally parses a mesh file.

---

### B. P1 Shape Functions Module
This module handles the linear (P1) shape functions used in the FEM discretization:
- Evaluates `phi[i]` and `grad_phi[i]` for each element
- Computes **elementary integrals** for mass and stiffness matrices

**Files:**
- `ShapeFunctions.hpp`: declarations for P1 shape functions, gradients, and element-wise integration.
- `ShapeFunctions.cpp`: implements `computeLocalMassMatrix()` and `computeLocalStiffnessMatrix()` for each triangular element.

---

### C. Assembly Module
Responsible for building the global system:
- Iterates over all elements  
- Computes local contributions and assembles them into **sparse global matrices**  
- Applies Dirichlet boundary conditions  
- Performs **mass lumping** for explicit time-stepping  

**Files:**
- `Assembler.hpp`: declares functions or a class to assemble global sparse matrices, apply boundary conditions, and perform mass-lumping.
- `Assembler.cpp`: implements element loops and assembly using `Eigen::SparseMatrix`.

---

### D. Time-Stepping Module
Handles the temporal evolution of the solution:
- Initializes \(U^0\) and \(U^1\) using explicit Taylor expansion  
- Performs the explicit time-stepping loop using the **central difference scheme**  
- Updates the solution at each time step

**Files:**
- `TimeIntegrator.hpp`: contains a class (e.g., `WaveSolver`) with methods like `initialize()` and `stepForward()`.
- `TimeIntegrator.cpp`: implements initialization and explicit update rules.

---

### E. Output Module
Manages saving and visualizing the solution:
- Exports data in **VTK format** for visualization in Paraview  
- Saves snapshot files at user-defined time intervals

**Files:**
- `Output.hpp`: declarations for writing VTK files and saving snapshots.
- `Output.cpp`: implements file writing.

---

### F. Main Program
The main program ties all modules together:
- Creates the mesh and assigns simulation parameters  
- Initializes the solver  
- Runs the time loop  
- Saves outputs for visualization

**File:** `main.cpp`



