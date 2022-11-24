using Gridap

dir = "examples/Gridap_Tutorials/output_17"

T = CartesianDiscreteModel((0,1,0,1),(20,20))
Ω = Interior(T)
dΩ = Measure(Ω,2)


refFE = ReferenceFE(lagrangian,Float64,1)
V = TestFESpace(T,refFE,dirichlet_tags = "boundary")

g(x,t::Real) = 0.0 #Boundary conditions
g(t::Real) = x -> g(x,t)

U = TransientTrialFESpace(V,g)


#Write weak form

k(t) = 1.0 + 0.95*sin(2π*t)
f(t) = sin(π*t)
#Now write the weak form with RHS being zero
res(t,u,v) = ∫(∂t(u)*v + k(t)*(∇(u)⋅∇(v)) - f(t)*v)dΩ

#This will perform automatic differentiation (make jacobian and stuff to do time-stepping) very nice that we do not even have to write it ourselves
op = TransientFEOperator(res,U,V)

linear_solver = LUSolver()#This problem is linear

Δt = 0.05
θ = 0.5
ode_solver = ThetaMethod(linear_solver,Δt,θ)#This determines what the ODE solver will be, can choose RK FE and some others. Will document this as this limits our scope.

u_zero = interpolate_everywhere(0.0,U(0.0))
t_zero = 0.0
t_end = 10.0
uht = solve(ode_solver,op,u_zero,t_zero,t_end)

if isdir(dir)
    createpvd(joinpath(dir,"poisson_transient_solution"*".pvd")) do pvd
        for (uₕ,t) in uht
          pvd[t] = createvtk(Ω,joinpath(dir,"poisson_transient_solution_$t"*".vtu"),cellfields=["u"=>uₕ])
        end
      end
else
    mkdir(dir)
    createpvd(joinpath(dir,"poisson_transient_solution"*".pvd")) do pvd
        for (uₕ,t) in uht
          pvd[t] = createvtk(Ω,joinpath(dir,"poisson_transient_solution_$t"*".vtu"),cellfields=["u"=>uₕ])
        end
      end
    end





