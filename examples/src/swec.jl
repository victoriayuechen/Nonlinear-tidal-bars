using DifferentialEquations, Setfield

includet("utils.jl")
includet("hydrodynamic_terms.jl")
includet("sediment_terms.jl")
includet("write_output.jl")


"""
Solve shallow water equations:
    ∂ζ/∂t + ∂(Du)/∂x + ∂(Dv)/∂y = 0
    ∂u/∂t + u∂u/∂x + v∂u/∂y + g∂ζ/∂x + cd/D |u⃗|u - fv = Fᵤ
    ∂v/∂t + u∂v/∂x + v∂v/∂y + g∂ζ/∂y + cd/D |u⃗|v + fu = Fᵥ
with,
    D = H+ζ-h
    Fᵤ = ∂u₀/∂t + cd √(u₀^2+v₀^2) u₀ - f v₀
    Fᵥ = ∂v₀/∂t + cd √(u₀^2+v₀^2) v₀ - f u₀
    u₀ = U cos(σ t) & v₀ = 0
    
Solve concentration equation:
    ∂c/∂t + ∇⋅q⃗ₛ  = αₛ H(uₑ²-u_crit^2)(uₑ²-u_crit^2) - γ c
with,
    q⃗ₛ = u⃗c - μ(∇C + c_top ∇ζ + c_bottom ∇h)
Calculate bed load sediment transport:
    q⃗b = αb H(uₑ²-u_crit^2)(uₑ²-u_crit^2)(u⃗ - Λuₑ∇h)

Grid:
    v   |v   v   v   v  | v
    c u |c u c u c u c u| c u
    v   |v   v   v   v  | v
    c u |c u c u c u c u| c u
    v   |v   v   v   v  | v
    c u |c u c u c u c u| c u
    v   |v   v   v   v  | v
where,
ζ: c grid
u: u grid
v: v grid
"""
function SWEC(A::MyArrays,Pars::Params;verbose=true)
    # unpack parameters
    @unpack Δt_num, t_span, save_moments_t = Pars
    @unpack output_residual = Pars

    # define ODEProblem
    p = (A,Pars,verbose)
    if output_residual
        prob = ODEProblem(rhs_swec!,A.ic,t_span,p,saveat=save_moments_t)
        savedvalues_swec = SavedValues(Float64, Tuple{Array{Float64,2},Array{Float64,2}})
        cb = SavingCallback(save_func_swec,saveat = save_moments_t,savedvalues_swec)
    else
        prob = ODEProblem(rhs_swec!,A.ic,t_span,p,save_everystep=false)
        cb = FunctionCallingCallback(save_func_swec,funcat = save_moments_t)
    end
    
    ∇h!(A,Pars); # calculate gradient of the bottom and save in A.dhdx, A.dhdy

    sol = solve(prob,RK4(),adaptive=false,dt=Δt_num,callback=cb)

    # update ζ,u,v,c with last value of time integration
    A.ic .= sol[:,:,:,end]
    
    #save residual ζ,u,v,c to residual_swec.h5
    if output_residual
        output_residual_swec(sol,savedvalues_swec.saveval,Pars)
        output_bottomheight(0,A,Pars)
    end
end

"""
rhs of ODEs
"""
function rhs_swec!(rhs,Ψ,p,t)
    # views for state variables and rhs
    ζ = @view Ψ[:,:,1]
    u = @view Ψ[:,:,2]
    v = @view Ψ[:,:,3]
    c = @view Ψ[:,:,4]
    dζdt = @view rhs[:,:,1]
    dudt = @view rhs[:,:,2]
    dvdt = @view rhs[:,:,3]
    dcdt = @view rhs[:,:,4]

    #unpack p
    A, Pars, _ = p

    # apply boundary conditions
    periodic_bc!(u,v,ζ,c)   # periodic boundary conditions
    noFlux_bc!(v)           # no water flux boundary condtion, note that no sediment flux boundary conditions of qsᵥ and qbᵥ are in q⃗s! and q⃗b!


    #precalculations
    precalculations!(A,ζ,u,v,c,Pars)    # calculate the depth, interpolate variables and calculate gradients

    # calculate terms on rhs for:
    # continuity equation
    continuity!(A,u,v,Pars)             # calculate ∇⋅(Du⃗)

    # momentum equation
    Fᵤ,Fᵥ = external_press_grad(t,Pars) # calculate (Fᵤ,Fᵥ)
    friction!(A,u,v,Pars)               # calculate cd|u⃗|u⃗/D
    advection!(A,u,v,Pars)              # calculate (u⃗⋅∇)u⃗
    pressgrad!(A,ζ,Pars)                # calculate g∇ζ
    coriolis!(A,Pars)                   # calculate -fv, fu

    # concentration of suspended sediment
    erosion!(A,Pars)                    # calculate erosion
    deposition!(A,Pars)                 # calculate deposition
    qs_advection!(A,u,v,Pars)           # calculate u⃗c
    qs_diff!(A,ζ,c,Pars)                # calculate -μ(∇C + c_top ∇ζ + c_bottom ∇h)
    q⃗s!(A,Pars)                         # calculate suspended sediment transport q⃗sᵤ and q⃗sᵥ
    div_q⃗s!(A,Pars)                     # calculate its divergence ∇⋅q⃗ₛ

    #bedload
    qb_advection!(A,u,v,Pars)           # calculate advective bed load transport
    qb_diff!(A,Pars)                    # calculate diffusive bed load transport (i.e., bed slope)
    q⃗b!(A,Pars)                         # calculate bed load sediment transport q⃗bᵤ and q⃗bᵥ
    div_q⃗b!(A,Pars)                     # calculate its divergence ∇⋅q⃗b

    # calculate rhs
    @unpack nx,ny = Pars
    @inbounds begin
        for j in 1:ny-1
            for i in 2:nx-1
                dζdt[i,j] = A.continuityᵤ[i,j] + A.continuityᵥ[i,j]
                dudt[i,j] = Fᵤ + A.advectionᵤ[i,j] + A.pressgradᵤ[i,j] + A.frictionᵤ[i,j] + A.coriolisᵤ[i,j]
                dcdt[i,j] = -A.div_q⃗s[i,j] + A.erosion[i,j] + A.deposition[i,j]
            end
        end
        for j in 2:ny-1
            for i in 2:nx-1
                dvdt[i,j] = Fᵥ + A.advectionᵥ[i,j] + A.pressgradᵥ[i,j] + A.frictionᵥ[i,j] + A.coriolisᵥ[i,j]
            end
        end
    end
    periodic_bc!(dζdt,dudt,dvdt,dcdt)
end

"""
output function to
    print progress in REPL 
    save ζ, u, v, c to hdf5 file
"""
function save_func_swec(Ψ,t,integrator)
    A, Pars, verbose = integrator.p
    @unpack t_end, output_hydro, output_residual = Pars

    if verbose
        print_progress(t,t_end) # print progress in REPL
    end

    if output_hydro
        output_swec(t,Ψ,A,Pars)   # save ζ, u, v, c to hdf5 file
    end
    if output_residual
        return (A.qsᵤ .+ A.qbᵤ, A.qsᵥ .+ A.qbᵥ)
    end
end

"""
output function to 
    return divergence of total sediment transport
"""
function save_div_q⃗(Ψ,t,integrator)
    A, Pars, verbose = integrator.p
    return A.div_q⃗b .+ A.div_q⃗s
end

"""
spinup function: integrate hydrodynamics for nr_of_tc_spinup tidal cycles.
Note that in SWEC the values of u,v,ζ and c at the last timestep are saved in A.
"""
function spinup(A,P,nr_of_tc_spinup;verbose=true)
    Pars_spinup = @set P.output_hydro = false
    Pars_spinup = @set Pars_spinup.output_residual = false
    if verbose
        println("spin-up for " * string(nr_of_tc_spinup) * " tidalcycle(s)..... ")
    end
    for TCn in 1:nr_of_tc_spinup
        if verbose
            println("spinup tidal cycle: " * string(TCn))
        end
        SWEC(A,Pars_spinup,verbose=false)
    end
    if verbose
        println("spin-up finished!")
    end
end
 