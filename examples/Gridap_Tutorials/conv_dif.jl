using Gridap


## This code solves the convection-diffusion equation on a 2D domain, with a initial condition and periodic boundary conditions on all sides. 
## Can be used to study the stability of some methods/play around with the transient FE solvers of Gridap.jl
dir = "examples/Gridap_Tutorials/output_conv_dif_eq"


#Making of the domain, this can very easily be changed to a 1D domain
domain = (0,1,0,1)
partition = (100,100)
model = CartesianDiscreteModel(domain,partition;isperiodic=(true,true))

order = 1
refFE = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(model,refFE,dirichlet_tags = "boundary")


# g(x,t::Real) = 0.0 #Boundary conditions
# g(t::Real) = x -> g(x,t)

U = TransientTrialFESpace(V)

degree = 2
Ω = Interior(model)
dΩ = Measure(Ω,degree)

velocity = VectorValue(1.0,0.0)

ff(t) = 0.0
res(t,u,v) = ∫(∂t(u)*v + 0.01*(∇(v)⋅∇(u)) + (velocity ⋅  ∇(u))*v - ff(t)*v)dΩ

op = TransientFEOperator(res,U,V)

ls = LUSolver()
sf(x) = sin(π*x[1])
Δt = 0.05
θ = 0.5
t_zero = 0.0
t_end = 5.0
ode_solver = ThetaMethod(ls,Δt,θ)

u_zero = interpolate_everywhere(sf,U(0.0))
uht = solve(ode_solver,op,u_zero,t_zero,t_end)

if isdir(dir)
    createpvd(joinpath(dir,"convdiv_transient_solution"*".pvd")) do pvd
        for (uₕ,t) in uht
          pvd[t] = createvtk(Ω,joinpath(dir,"convdiv_transient_solution_$t"*".vtu"),cellfields=["u"=>uₕ])
        end
      end
else
    mkdir(dir)
    createpvd(joinpath(dir,"convdiv_transient_solution"*".pvd")) do pvd
        for (uₕ,t) in uht
          pvd[t] = createvtk(Ω,joinpath(dir,"convdiv_transient_solution_$t"*".vtu"),cellfields=["u"=>uₕ])
        end
      end
    end

