using Gridap


T = CartesianDiscreteModel((0,1,0,1),(20,20))
Ω = Interior(T)
dΩ = Measure(Ω,2)


refFE = ReferenceFE(lagrangian,Float64,1)
V = TestFESpace(T,refFE,dirichlet_tags = "boundary")

g(x,t::Real) = 0.0
g(t::Real) = x -> g(x,t)

U = TransientTrialFESpace(V,g)


#Write weak form

k(t) = 1.0 + 0.95*sin(2π*t)
f(t) = sin(π*t)
#Now write the weak form with RHS being zero
res(t,u,v) = ∫(∂t(u)*v + k(t)*(∇(u)⋅∇(v)) - f(t)*v)dΩ

#This will perform automatic differentiation (make jacobian and stuff to do time-stepping)
op = TransientFEOperator(res,U,V)

linear_solver = LUSolver()

Δt = 0.05
θ = 0.5
ode_solver = ThetaMethod(linear_solver,Δt,θ)

u_zero = interpolate_everywhere(0.0,U(0.0))
t_zero = 0.0
t_end = 10.0
uht = solve(ode_solver,op,u_zero,t_zero,t_end)

createpvd("poisson_transient_solution") do pvd
    for (uₕ,t) in uht
      pvd[t] = createvtk(Ω,"poisson_transient_solution_$t"*".vtu",cellfields=["u"=>uₕ])
    end
  end

