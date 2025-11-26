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

### Term-by-term analysis

---

## 1. Partial Differential Equation (PDE)

$$
\frac{\partial^2 u}{\partial t^2} - \Delta u = f 
\quad \text{in } \Omega \times (0,T]
$$

### Hyperbolic term
`$ \frac{\partial^2 u}{\partial t^2} $`  
Represents the acceleration of the field over time. This is the key term that makes the equation hyperbolic and describes wave propagation.

### Spatial term
`$ \Delta u $`  
Represents restoring forces or internal stresses within the medium. It is the source of spatial propagation of the motion.

### Forcing term
`$ f $`  
Represents an external energy source or applied force acting on the system within the domain. If `$ f = 0 $`, the equation is homogeneous and describes freely propagating waves.

**Note:**  
The standard form of the wave equation includes a propagation constant $c^2$:
$$
\frac{\partial^2 u}{\partial t^2} - c^2 \Delta u = f.
$$
In our problem, it is implicitly assumed that $c = 1$ (unit propagation speed).

---

## 2. Boundary Condition

$$
u = g \quad \text{on } \partial\Omega \times [0,T]
$$

This is a **Dirichlet boundary condition**.  
It specifies the value of the wave field $u$ along the boundary $\partial\Omega$ of the domain for all times.   $g$ is a known function.  
If $g = 0$, the boundary is “fixed,” or—depending on interpretation—an “ideal absorbing boundary.”

---

## 3. Initial Conditions

Because the wave equation is second order in time, it requires **two initial conditions** to define the solution uniquely:

- $u_0$: the initial configuration (initial displacement) of the wave field at time $t = 0$.  
- $u_1$: the initial velocity of the wave field at time $t = 0$.

### Weak Form

For the homogeneous equation (with `$ f = 0 $` and `$ g = 0 $`), multiplying by a test function  
`$ v \in H_0^1(\Omega) $` and integrating over `$\Omega$`, we get:

`$ \int_{\Omega} \frac{\partial^2 u}{\partial t^2} \, v \, dx - \int_{\Omega} \Delta u \, v \, dx = 0 $`

Applying integration by parts (Green's theorem) to the Laplacian term:

`$ \int_{\Omega} \frac{\partial^2 u}{\partial t^2} \, v \, dx + \int_{\Omega} \nabla u \cdot \nabla v \, dx = 0 $`

The FEM implementation will be based on the discretization of this weak form.

### Implications for the Numerical Solution

The hyperbolic nature of the equation is crucial for our implementation:

---

## 1. Time Discretization

Since the problem is time-dependent (evolution in time), you will need to use a time integration method such as:

- the **Centered Finite Difference Method** (or **Leapfrog scheme**), which is classically used for the wave equation,  
- or more advanced methods such as **Runge–Kutta** schemes.

The choice of the time-stepping method will strongly affect **numerical stability**.

---

## 2. Space Discretization

### Space Discretization

This is where the Finite Element Method (FEM), which is required, comes into play.  
To use FEM, it is needed to reformulate the problem in its weak (or variational) form within an appropriate functional space (typically $H^1(\Omega)$).



