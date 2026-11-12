# the following code creates a 2D visual representing how two components (dilute benzene in cyclohexane) mix and flow in a Y-Junction
# the visual shows a) flow speed and b) distribution of benzene in the Y-Junction
# maths behind it, the Finite Element Method is explained in my GitHub description
# essentially the Y-Junction shape is discretised into a mesh of triangles. the weak form then turns the PDE's into one linear equation per node,
# which are all solved together as a matrix system [K]{u}={F}. this is done twice: first for the flow, then for the mixing

import netgen.occ as ngocc #it uses netgen to create the Y-Junction mesh
import matplotlib.pyplot as plt #matplot used to make the visual
from firedrake import * #it uses firedrake to solve PDE's from my governing equations as a system of equations 
from firedrake.pyplot import tripcolor

#key parameters
U = 1e-3          # average speed of the liquid entering each inlet, m/s   (1 mm/s)
mu = 9e-4         # viscosity of cyclohexane, Pa s
D = 2.1e-9        # diffusivity of benzene in cyclohexane, m^2/s  (for Fick's law eqn.)
B_in = 1.0        # benzene mol % in inlet 1   (inlet 2 has 0 %). B = x_B in my notes


# making the mesh. x and y values describe the shape.
# corner points of the Y-junction in micrometres (channels are 100 um wide)
um = 1e-6 # um is micrometers
corners = [((1200, -50), "outlet"),
           ((1200, 50), "wall"),
           ((20.7, 50), "wall"),
           ((-176.8, 247.5), "inlet_1"),
           ((-247.5, 176.8), "wall"),
           ((-70.7, 0), "wall"),
           ((-247.5, -176.8), "inlet_2"),
           ((-176.8, -247.5), "wall"),
           ((20.7, -50), "wall")]

wp = ngocc.WorkPlane()
wp.MoveTo(corners[0][0][0] * um, corners[0][0][1] * um)
for i in range(len(corners)):
    next_point = corners[(i + 1) % len(corners)][0]
    wp.LineTo(next_point[0] * um, next_point[1] * um, name=corners[i][1])

h = 5 * um  # triangle size. mesh check: re-running with h = 10 um gives the same benzene value at the outlet, so 5 um is fine
ngmesh = ngocc.OCCGeometry(wp.Face(), dim=2).GenerateMesh(maxh=h)
mesh = Mesh(ngmesh)

# now work out which boundary number is at which edge, so boundary conditions can be applied
# firedrake doesn't read names like 'inlet_1', it reads the edges of each shape in order they were drawn (1,2,3..)
# python counts from 0 but firedrake counts edges from 1, hence the 'i + 1'
names = ngmesh.GetRegionNames(dim=1)
inlet_1 = [i + 1 for i, n in enumerate(names) if n == "inlet_1"] #assigning names to boundary numbers
inlet_2 = [i + 1 for i, n in enumerate(names) if n == "inlet_2"]
walls = [i + 1 for i, n in enumerate(names) if n == "wall"]
outlet = [i + 1 for i, n in enumerate(names) if n == "outlet"]
print(f"Mesh: {mesh.num_cells()} triangles") #prints how many triangles the mesh has

# Equations for flow: Stokes equations
#   -mu * lap(u) + grad(p) = 0   (pressure push balances viscous friction)
#    div(u) = 0                  (conservation of mass, liquid is incompressible)
# inertia neglected because Re ≈ 0.09 (viscous forces >>> inertial forces), body forces (gravity) neglected too
# "CG" = continuous Galerkin. quadratic velocity + linear pressure is called a Taylor-Hood pair,
# using the same order for both makes the system unstable and the solve fails
V = VectorFunctionSpace(mesh, "CG", 2)     # velocity: quadratic on each triangle
Q = FunctionSpace(mesh, "CG", 1)           # pressure: linear on each triangle
W = V * Q                                  # solve for velocity and pressure together

u, p = TrialFunctions(W)                   # unknowns
v, q = TestFunctions(W)                    # test (weight) functions

# weak form of the Stokes equations (Firedrake turns this into [K]{u} = {F})
a = (mu * inner(grad(u), grad(v)) - p * div(v) + q * div(u)) * dx
L = inner(Constant((0, 0)), v) * dx         # right hand side = 0, no body forces

# boundary conditions: the liquid enters each inlet with average speed U, doesn't move at the walls
def inlet_profile(centre_x, centre_y, dir_x, dir_y): 
    """laminar (parabolic) profile: average speed U, zero at the walls, fastest in the middle.

    (centre_x, centre_y) = middle of the inlet in um, (dir_x, dir_y) = direction of flow.
    """
    x, y = SpatialCoordinate(mesh)
    s = ((x - centre_x * um) * dir_y - (y - centre_y * um) * dir_x) / (100 * um)   # -0.5 to 0.5 across the inlet
    return 1.5 * U * (1 - 4 * s**2) * as_vector([dir_x, dir_y])

bcs = [DirichletBC(W.sub(0), inlet_profile(-212.15, 212.15, 0.7071, -0.7071), inlet_1),  # flowing down-right
       DirichletBC(W.sub(0), inlet_profile(-212.15, -212.15, 0.7071, 0.7071), inlet_2),  # flowing up-right
       DirichletBC(W.sub(0), as_vector([0, 0]), walls)]                                 # no slip condition
# these are all Dirichlet BCs (fixed values). the outlet has no BC, so the liquid leaves freely

flow = Function(W)
solve(a == L, flow, bcs=bcs)                # firedrake builds [K]{u}={F} and solves it with LU (gaussian elimination)
velocity, pressure = flow.subfunctions

# mixing equations: convection and diffusion (Fick's law)
#   u . grad(B) = D * lap(B)     (benzene carried by flow = benzene diffusing high -> low)
X = FunctionSpace(mesh, "CG", 2)           # benzene mol %: quadratic on each triangle
B = TrialFunction(X)
y = TestFunction(X)                        # test function (w in my notes), one 'hill' for each node

# weak form of mixing equations: multiply by test function, integrate, integrate by parts
a = (dot(velocity, grad(B)) * y + D * inner(grad(B), grad(y))) * dx
L = Constant(0) * y * dx

# boundary conditions
bcs = [DirichletBC(X, B_in, inlet_1),      # 1 mol% benzene enters through inlet 1
       DirichletBC(X, 0.0, inlet_2)]       # pure cyclohexane enters through inlet 2
# Neumann BCs on the walls and outlet (dB/dn = 0, no benzene diffuses through) need no code,
# they come out automatically from integrating by parts in the weak form

benzene = Function(X)
solve(a == L, benzene, bcs=bcs)

# mass balance check: everything that flows in must flow out, so the average
# benzene at the outlet should be halfway between the two inlets: 0.5 mol %
n = FacetNormal(mesh)
outflow = assemble(dot(velocity, n) * ds(tuple(outlet)))
benzene_out = assemble(benzene * dot(velocity, n) * ds(tuple(outlet)))
print(f"Average benzene at outlet: {benzene_out / outflow:.3f} mol %   (expected 0.500)")

# flow check: far from the junction the flow is laminar channel flow (a parabola), where peak speed = 1.5 x average speed.
# average = 2U (both inlets join into one channel of the same width), so peak should be 3U = 3 mm/s
peak = PointEvaluator(mesh, [[800 * um, 0]]).evaluate(velocity)[0][0]
print(f"Peak speed in main channel: {peak * 1000:.3f} mm/s   (expected {3 * U * 1000:.3f})")

# value to compare for the mesh check (benzene near the outlet, next to the top wall). should not change with h = 10 um
B_probe = PointEvaluator(mesh, [[1190 * um, 45 * um]]).evaluate(benzene)[0]
print(f"Benzene at outlet, near top wall: {B_probe:.4f} mol %")

# plotting Y-Junction for a) flow speed and b) mixing
fig, axes = plt.subplots(2, 1, figsize=(11, 7))
speed = Function(X).interpolate(sqrt(dot(velocity, velocity)) * 1000)        # in mm/s
c1 = tripcolor(speed, axes=axes[0], cmap="viridis")
fig.colorbar(c1, ax=axes[0], label="speed (mm/s)")
axes[0].set_title("Flow (Stokes equations)")
c2 = tripcolor(benzene, axes=axes[1], cmap="RdBu_r", vmin=0, vmax=B_in)
fig.colorbar(c2, ax=axes[1], label="benzene (mol %)")
axes[1].set_title("Mixing (Fick's law): benzene mol %, the rest is cyclohexane")
for ax in axes:
    ax.set_aspect("equal")
    ax.axis("off")
fig.savefig("ymixer.png", dpi=150, bbox_inches="tight")
print("Saved ymixer.png")
