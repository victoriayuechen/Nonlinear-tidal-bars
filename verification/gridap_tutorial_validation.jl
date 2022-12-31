using Gridap

u(x) = x[1] + x[2]
f(x) = 0

# We also need to define the gradient of $u$ since we will compute the $H^1$ error norm later. In that case, the gradient is simply defined as
#

# ∇u(x) = VectorValue(1,1)

# Note that we have used the constructor `VectorValue` to build the vector that represents the gradient. However, we still need a final trick. We need to tell the Gridap library that the gradient of the function `u` is available in the function `∇u` (at this moment `u` and `∇u` are two standard Julia functions without any connection between them). This is done by adding an extra method to the function `gradient` (aka `∇`) defined in Gridap:

import Gridap: ∇
∇(::typeof(u)) = ∇u

# Now, it is possible to recover function `∇u` from function `u` as `∇(u)`. You can check that the following expression evaluates to `true`.

∇(u) === ∇u

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)

domain3d = (0,1,0,1,0,1)
partition3d = (4,4,4)
model3d = CartesianDiscreteModel(domain3d,partition3d)


order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
V0 = TestFESpace(model,reffe,conformity=:H1,dirichlet_tags="boundary")
U = TrialFESpace(V0,u)

degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

a(u,v) = ∫( ∇(v)⊙∇(u) )*dΩ
b(v) = ∫( v*f )*dΩ

op = AffineFEOperator(a,b,U,V0)

uh = solve(op)

e = u - uh


writevtk(Ω,"error",cellfields=["e" => e, "man"=>u])