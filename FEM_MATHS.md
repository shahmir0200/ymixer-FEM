# FEM maths explained

How the maths behind `ymixer.py` works, from the governing equations to the matrix the computer solves. Notation: **u** = velocity vector (u, v) in x and y, p = pressure, x_B = benzene mole fraction (written in mol % in the code, where it's called `B`).

---

## Background: linear vs non-linear

- **Linear** differential equation → the unknown (y) only appears as y, y′, y″… (not multiplied by itself)
- **Non-linear** → anything else, e.g. √y, y³, or y × y′

This matters because linear equations turn into a linear system that the computer can solve in one go. Non-linear ones need repeated solving (e.g. Newton's method).

---

## ① Navier–Stokes → Stokes (flow)

The **Navier–Stokes equations** describe how a fluid flows. There are two:

**1. Conservation of mass** → mass can't be created or destroyed. For an incompressible liquid:

$$\nabla\cdot\mathbf{u} = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0$$

**2. Momentum** → Newton's 2nd law for fluids (momentum = mass × velocity):

$$\rho\left(\frac{\partial \mathbf{u}}{\partial t} + \mathbf{u}\cdot\nabla\mathbf{u}\right) = -\nabla p + \mu\nabla^2\mathbf{u} + \rho\mathbf{g}$$

| Term | Name | Meaning |
|---|---|---|
| ρ(∂**u**/∂t + **u**·∇**u**) | inertial term | liquid's resistance to changing speed/direction |
| −∇p | pressure gradient | liquid flows from high → low pressure |
| μ∇²**u** | viscous term | friction & drag inside the liquid (layers dragging on each other and on the walls) |
| ρ**g** | body forces | external forces, e.g. gravity |

**Simplifying for my Y-junction (microfluidics):**

- **∂u/∂t = 0** → steady state assumed (nothing changes with time)
- **Inertial term → 0.** Why? → Reynolds number:

$$\mathrm{Re} = \frac{\rho U L}{\mu} = \frac{\text{inertial forces}}{\text{viscous forces}}$$

  With cyclohexane (ρ = 774 kg/m³, μ = 9 × 10⁻⁴ Pa·s), U = 1 mm/s and L = 100 µm → Re ≈ 0.09. If Re ≪ 1 → viscous forces >>> inertial forces → inertia neglected. Think: tiny channel + slow liquid → friction wins.
- **Body forces → 0.** Why? → the chip lies flat, so gravity points into the page, not along the flow. And for a liquid with constant density, gravity would only shift the pressure, not change the velocity.

Bonus: the inertial term **u**·∇**u** is non-linear (u × its own derivative), so dropping it also makes the equations **linear**.

**∴ Stokes equations (microfluidics):**

$$\nabla\cdot\mathbf{u} = 0, \qquad -\nabla p + \mu\nabla^2\mathbf{u} = 0$$

In words: the pressure push is balanced by viscous friction, and the liquid is incompressible.

---

## ② Mixing: flux of benzene + mass balance

**Flux of benzene** (N_B) = rate at which benzene passes through an area, per unit time:

**Flux = convection + diffusion**

- **convection** → benzene carried along by the flow: x_B **u**
- **diffusion** → benzene diffusing into cyclohexane, high → low concentration: **Fick's law**

**Fick's law** → 1D (x only): J_B = −D_B dC_B/dx. My project is 2D, so I need the x and y directions:

$$\mathbf{J}_B = -D_B\left(\frac{\partial x_B}{\partial x},\, \frac{\partial x_B}{\partial y}\right) = -D_B\nabla x_B$$

∇ ('del') = rate of change in the x and y directions. D_B = 2.1 × 10⁻⁹ m²/s (value from the Firedrake demo).

$$\therefore\; \mathbf{N}_B = x_B\,\mathbf{u} - D_B\nabla x_B$$

Fick's law works here because the benzene is **dilute** (1 mol %). For concentrated mixtures you'd need the full Onsager–Stefan–Maxwell equations, like Farrell's group.

**Mass balance** → in = out at every point. Why? → steady state (no build-up) and no reaction (no benzene made or destroyed). "Net flux out of a point" is the **divergence**, so:

$$\nabla\cdot\mathbf{N}_B = 0 \quad\Rightarrow\quad \nabla\cdot(x_B\,\mathbf{u}) - \nabla\cdot(D_B\nabla x_B) = 0$$

Simplify: ∇·(x_B **u**) = **u**·∇x_B + x_B(∇·**u**), and ∇·**u** = 0 from ① (incompressible). D_B is constant, so ∇·(D_B∇x_B) = D_B∇²x_B.

**∴ Convection–diffusion equation (aka mixing equation):**

$$\mathbf{u}\cdot\nabla x_B = D_B\nabla^2 x_B$$

benzene carried by the flow = benzene spreading by diffusion.

**We have now derived:**

| | Equations |
|---|---|
| Flow (Stokes, microfluidics) | ∇·**u** = 0 and μ∇²**u** − ∇p = 0 |
| Mixing (convection–diffusion) | **u**·∇x_B = D_B∇²x_B |

However, we're left with **PDEs** (because of the ∇'s), which can't be solved exactly for a Y-shaped channel. Next: turn them into linear equations the computer can solve, using FEM.

---

## ③ Boundary conditions

PDEs → solved by **setting boundary conditions**, then solving as **linear equations**. Without boundary conditions there's no unique answer.

**1) Dirichlet boundary condition** = fixed **value** set at the boundary, e.g.:
- inlets: parabolic laminar velocity profile, average 1 mm/s
- walls: **u** = 0 (no-slip: fluid in contact with a solid has no velocity relative to it)
- benzene: x_B = 1 mol % at inlet 1, 0 at inlet 2

**2) Neumann boundary condition** = fixed **rate of change across** the boundary (n = direction perpendicular to it), e.g.:
- walls and outlet: ∂x_B/∂n = 0 → no benzene diffuses through (Fick's law: no gradient = no diffusion)

**Outlet (flow):** liquid leaves freely (no condition imposed).

| Boundary | Flow | Benzene |
|---|---|---|
| Inlet 1 | parabolic profile, average 1 mm/s | x_B = 1 mol % (Dirichlet) |
| Inlet 2 | parabolic profile, average 1 mm/s | x_B = 0 (Dirichlet) |
| Walls | **u** = 0, no-slip (Dirichlet) | ∂x_B/∂n = 0 (Neumann) |
| Outlet | leaves freely | ∂x_B/∂n = 0 (Neumann) |

---

## ④ Weak form

**Problem:** how do we find values at every point on a complex shape?

- **Strong form** (our PDE) → must be exactly true at every single point
- **Weak form** → only requires the equation to **balance on average** (= 0) over a small area

**'True' equation** means it equals 0, e.g. x + y = 0 holds at every point. The opposite → an 'error' at a point, called the **residual**.

**Why we need it:** the mesh has ~33,000 nodes → 33,000 unknown benzene values (x_B) → ∴ we need 33,000 equations. Between nodes, x_B is approximated by simple polynomial pieces (joining the dots), which can't match the true smooth curve exactly → a **residual** is left over. So we can't make the equation true at every point, but we can make it true on average around each node.

**How:**

1. Write the solution as polynomial pieces with **unknown values at the nodes** (like a trial solution in ODEs)
2. **Multiply** the strong form by a weight w_i (**test function**) that's specific to each node: w_i = 1 at node i, 0 at every other node (a little 'hill' on node i)
3. **Integrate** over the whole shape (since w_i = 0 away from node i, this only adds up the residual around node i), then **integrate by parts**
4. **Set the result = 0**
5. Now we have **one linear equation per node** (w_1, w_2, w_3, … → equation 1, 2, 3, …)
6. Computer solves these all together in a matrix

**Why integrate by parts?** ∇² is a second derivative, which polynomial pieces with kinks between triangles can't handle. Integration by parts moves one derivative onto w (1D: ∫ f″w dx = [f′w] − ∫ f′w′ dx), so only first derivatives are left. Bonus: the boundary term is where the Neumann conditions (∂x_B/∂n = 0) come in automatically, which is why they need no code.

**Result:**

Weak form of the convection–diffusion equation ✓

$$\int (\mathbf{u}\cdot\nabla x_B)\, w \; dA + D_B \int \nabla x_B\cdot\nabla w \; dA = 0$$

Weak form of the Stokes equations ✓ (test functions **v** for velocity and q for pressure)

$$\int \mu\,\nabla\mathbf{u}:\nabla\mathbf{v}\; dA - \int p\,(\nabla\cdot\mathbf{v})\; dA + \int q\,(\nabla\cdot\mathbf{u})\; dA = 0$$

---

## ⑤ System of equations: [K]{x_B} = {F}

Each node's weak-form equation becomes **one row** of a matrix. Small 1D example (5 nodes, x_B = 1 at the left end and 0 at the right end, diffusion only). The weak form gives one equation for each middle node:

```
x_B0 − 2·x_B1 + x_B2 = 0
x_B1 − 2·x_B2 + x_B3 = 0
x_B2 − 2·x_B3 + x_B4 = 0
```

Each equation links a node to its neighbours, so they must be solved together. Put in the known boundary values (x_B0 = 1, x_B4 = 0) and move them to the right-hand side:

$$\begin{bmatrix} -2 & 1 & 0 \\ 1 & -2 & 1 \\ 0 & 1 & -2 \end{bmatrix} \begin{Bmatrix} x_{B1} \\ x_{B2} \\ x_{B3} \end{Bmatrix} = \begin{Bmatrix} -1 \\ 0 \\ 0 \end{Bmatrix}$$

- **[K]** = matrix of coefficients (the numbers from the weak form)
- **{x_B}** = vector of unknowns
- **{F}** = vector of known values (from the boundary conditions)

Computer solves [K]{x_B} = {F} for x_B → here x_B = 0.75, 0.5, 0.25 ✓

In my real 2D project it's the same, just bigger: K is about 33,000 × 33,000 for the mixing (and ~74,000 unknowns for the flow). Most entries are 0 (each node only connects to its neighbours), so it solves in seconds using LU decomposition (Gaussian elimination).
