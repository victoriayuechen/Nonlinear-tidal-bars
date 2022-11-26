using Gridap

##NOT WORKING YET

domain = (0,1)
partition = (100)
model = CartesianDiscreteModel(domain,partition;isperiodic=(true,))

order = 1
refFE_u = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(model,refFE_u,dirichlet_tags = "boundary")

refFE_h = ReferenceFE(lagrangian,Float64,order)
Q = TestFESpace(model,refFE_h,dirichlet_tags ="boundary")


VT = TransientTrialFESpace(V)
QT = TransientTrialFESpace(Q)
Y = TransientMultiFieldFESpace([V,Q])
X = TransientMultiFieldFESpace([VT,QT])

degree = order
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)

res(u,v,h) = ∫(∂t(u)*v + v⋅∇(u) + v⋅∇(h))dΩ
op = TransientFEOperator(res,X,Y)

ls = LUSolver()
sf(x) = sin(π*x[1])
u_zero = interpolate_everywhere(0.0,VT(0.0))
h_zero = interpolate_everywhere(sf,QT(0.0))
uh_zero = interpolate_everywhere([u_zero,h_zero],X(0.0))
Δt = 0.05
θ = 0.5
t_zero = 0.0
t_end = 5.0
ode_solver = ThetaMethod(ls,Δt,θ)

uht,pht = solve(ode_solver,op,uh_zero,t_zero,t_end)
