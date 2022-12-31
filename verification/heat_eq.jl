using Pkg
Pkg.activate(".")

using Gridap
using WriteVTK

import Gridap:  ∇ 
# ∇(::typeof(u)) = ∇u


# domain size

LeftX = 0
RightX = 2

LeftY = 0
RightY = 2

x_c = (RightX - LeftX)/2
y_c = (RightY - LeftY)/2

domain = (LeftX, RightX, LeftY, RightY)
partition = (50, 50)

model = CartesianDiscreteModel(domain,partition)

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"bottom",[1,2,5])
add_tag_from_tags!(labels,"left",[7])
add_tag_from_tags!(labels,"right",[8])
add_tag_from_tags!(labels,"top",[3,4,6])
add_tag_from_tags!(labels,"inside",[9])
DC = ["left","right","top","bottom"]


degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

# Γ = BoundaryTriangulation(model)
# dΓ = Measure(Γ, degree)

order = 2
reffe = ReferenceFE(lagrangian, Float64, order)
V0 = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags=DC)

g(x) = 0.0

Ug = TrialFESpace(V0, g)

# force function
alpha = 40
f(x) =  exp(-alpha * (x[1] - x_c)^2 - alpha * (x[2] - y_c)^2)

# coefficient function
k1(x) = 1 + 1.5 * (x[1])
k2(x) = 1

# a1(u,v) = ∫( ∇(v)⋅∇(u) )*dΩ
# a(u,v) = ∫((∇ ⋅ ∇(u))*v)*dΩ
a1(u,v) = ∫( (k1*∇(u)) ⋅ (∇(v)))*dΩ
a2(u,v) = ∫( (k2*∇(u)) ⋅ (∇(v)))*dΩ
b(v) = ∫( v*f )*dΩ

op1 = AffineFEOperator(a1, b, Ug, V0)
op2 = AffineFEOperator(a2, b, Ug, V0)

ls = LUSolver()
solver = LinearFESolver(ls)

u_sol1 = solve(solver, op1)
u_sol2 = solve(solver, op2)

diff = abs(u_sol1 - u_sol2) / u_sol2

writevtk(Ω,"results",cellfields=["u1"=>u_sol1, "u2"=>u_sol2, "diff"=>diff, "coeff"=>k1])