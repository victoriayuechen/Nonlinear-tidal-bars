using Gridap
using Plots


LeftX = 0
RightX = 2

LeftY = 0
RightY = 2

x_c = (RightX - LeftX)/2
y_c = (RightY - LeftY)/2

domain = (LeftX, RightX, LeftY, RightY)
partition = (100, 100)

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

# ## FE space

# In this tutorial we will use linear Lagrangian Finite Elements.
refFE = ReferenceFE(lagrangian,Float64,2)

# The space of test functions is constant in time and is defined in steady problems:
V = TestFESpace(model,refFE,dirichlet_tags=DC)

# g(x,t::Real) = 0.0
# g(t::Real) = x -> g(x,t)

# g(x, t) = 0.0
um(x,t::Real) = sin(π*t/10) * (x[1]^2 + x[2]^2)
um(t::Real) = x -> um(x,t)

U = TransientTrialFESpace(V,um) # automatic BC (not sure if it works like that)

# manufactured time derivative
um_t(x,t::Real) = (π/10) * cos(π * t/10) * (x[1]^2 + x[2]^2)
um_t(t::Real) = x -> um_t(x,t)

#manufactured gradient
um_g(x,t::Real) = VectorValue(2 * sin(π*t/10) * x[1], 2 * sin(π*t/10) * x[2])
um_g(t::Real) = x -> um_g(x,t)


# corresponding force function
f(x,t::Real) = -((π/10) * cos(π * t/10) * (x[1]^2 + x[2]^2) + 4 * sin(π * t/10))
f(t::Real) = x -> f(x,t)
# ## Weak form

# κ(t) = 1.0 + 0.95*sin(2π*t)
κ(t) = 1.0

# f(t) = sin(π*t)
res(t,u,v) = ∫( ∂t(u)*v + κ(t)*(∇(u)⋅∇(v)) - f(t)*v )dΩ
jac(t,u,du,v) = ∫( κ(t)*(∇(du)⋅∇(v)) )dΩ
jac_t(t,u,duₜ,v) = ∫( duₜ*v )dΩ
op = TransientFEOperator(res,jac,jac_t,U,V)


linear_solver = LUSolver()

# Then, we define the ODE solver. That is, the scheme that will be used for the time integration. In this tutorial we use the `ThetaMethod` with $\theta = 0.5$, resulting in a 2nd order scheme. The `ThetaMethod` function receives the linear solver, the time step size $\Delta t$ (constant) and the value of $\theta $.
Δt = 0.05
θ = 0.5
ode_solver = ThetaMethod(linear_solver,Δt,θ)

# Finally, we define the solution using the `solve` function, giving the ODE solver, the FE operator, an initial solution, an initial time and a final time. To construct the initial condition we interpolate the initial value (in that case a constant value of 0.0) into the FE space $U(t)$ at $t=0.0$.
u₀ = interpolate_everywhere(0.0,U(0.0))
t₀ = 0.0
T = 10.0
uₕₜ = solve(ode_solver,op,u₀,t₀,T)

# ## Postprocessing

# We should highlight that `uₕₜ` is just an _iterable_ function and the results at each time steps are only computed when iterating over it, i.e., lazily. We can post-process the results and generate the corresponding `vtk` files using the `createpvd` and `createvtk` functions. The former will create a `.pvd` file with the collection of `.vtu` files saved at each time step by `createvtk`. The computation of the problem solutions will be triggered in the following loop:

el2 = Float64[] # L2 norm error array
eh1 = Float64[] # H2 norm eroor array
t_base = Float64[] # time base (I think Julia needs it to plot things...)

createpvd("poisson_transient_solution") do pvd
  for (uₕ,t) in uₕₜ
    e = um(t) - uₕ

    push!(el2, sqrt(sum( ∫( e*e )*dΩ )))
    push!(eh1, sqrt(sum( ∫( e*e + ∇(e)⋅∇(e) )*dΩ )))
    push!(t_base, t)
    pvd[t] = createvtk(Ω,"poisson_transient_solution_$t"*".vtu",cellfields=["u"=>uₕ, "force"=>f(t), "man"=>um(t), "error"=>(um(t) - uₕ), "e_grad"=>(um_g(t) - ∇(uₕ))])
  end
end

plot(t_base,[el2, eh1],
    label=["L2, H1"],
    shape=:auto,
    xlabel="time",ylabel="error norm")

# ![](../assets/poisson_transient/poisson_transient.gif)
