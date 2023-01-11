using DifferentialEquations

includet("swec.jl")

"""
solve ∂h/∂t = -1/(1-p̂) ⟨∇⋅q⃗⟩
where ⟨ ⟩ denotes the average over one tidalcycle
and q⃗ the total sediment transport
"""
function bed_evolution(A::MyArrays, Pars::Params; verbose=true)
    #prepare SWEC runs
    @unpack Δt_num, t_span, save_moments_t = Pars
    p_swec = (A,Pars,true)
    savedvalues_swec = SavedValues(Float64, Array{Float64,2})
    cb_swec          = SavingCallback(save_div_q⃗, saveat=save_moments_t, savedvalues_swec)
    prob_swec        = ODEProblem(rhs_swec!,A.ic,t_span,p_swec,save_everystep=false)
    integrator_swec  = init(prob_swec,RK4(),dt=Δt_num,callback=cb_swec)

    # unpack parameters
    @unpack Δτ_num,τ_span,save_moments_τ = Pars
    
    p = (A, Pars, verbose, integrator_swec, savedvalues_swec)

    #callback for output
    cb = FunctionCallingCallback(save_func_morpho,funcat = save_moments_τ)
    
    # define ODEProblem
    prob = ODEProblem(dhdτ!,A.h,τ_span,p,save_everystep=false)

    # Solve ODEProblem    
    sol = solve(
        prob,
        Euler(),
        adaptive=false,
        dt=Δτ_num,       #units: ε * seconds
        callback=cb,
    )
    return sol
end

"""
dhdτ
"""
function dhdτ!(dhdτ,h,p,τ)
    A, Pars, verbose, integrator_swec, savedvalues_swec = p
    
    #unpack parameters
    @unpack nx,ny,p̂,ε = Pars
   
    #save bottom height h to A
    A.h .= h
    ∇h!(A,Pars)
    
    #integrate swec for 1 tidal cycle
    reinit!(integrator_swec)
    solve!(integrator_swec)
    # @assert(maximum(abs.(A.ic .- integrator_swec.sol[:,:,:,end])) < 1e-4,"hydrodynamics not periodic (spin-up effect)")
    #update ic with last values of ζ, u, v, c
    A.ic .= integrator_swec.sol[:,:,:,end]

    #integrate swec for 2 tidal cycles (1 spinup)
    # A.ic .= zeros(nx,ny,4)
    # reinit!(integrator_swec)
    # solve!(integrator_swec)
    # # @assert(maximum(abs.(A.ic .- integrator_swec.sol[:,:,:,end])) < 1e-4,"hydrodynamics not periodic (spin-up effect)")
    # #update ic with last values of ζ, u, v, c
    # A.ic .= integrator_swec.sol[:,:,:,end]
    # reinit!(integrator_swec,A.ic)
    # solve!(integrator_swec)

    #calculate average divergence of total sediment transport A.avg_div_q⃗
    mean_sediment_transport!(A,savedvalues_swec.saveval,Pars)

    #calculate rhs
    @. dhdτ = -1.0/(1.0-p̂) * A.avg_div_q⃗ / ε
    periodic_bc!(dhdτ)
end

"""
calculate mean sediment transport on physcial domain
"""
function mean_sediment_transport!(A,saveval,Pars)
    @unpack nx,ny = Pars
    Nt = length(saveval)
    A.avg_div_q⃗ .= zeros(eltype(A.h),nx,ny)
    for j in 1:ny-1
        for i in 2:nx-1
            for tn in 1:Nt
                A.avg_div_q⃗[i,j] += saveval[tn][i,j]
            end
            A.avg_div_q⃗[i,j] /= Nt
        end
    end
    periodic_bc!(A.avg_div_q⃗)
end


"""
saving callback
"""
function save_func_morpho(Ψ,τ,integrator)
    A, Pars, verbose = integrator.p
    @unpack τ_end, output_morpho = Pars

    if verbose
        print_progress_morpho(τ,τ_end,Pars); # print progress in REPL
    end

    if output_morpho
        output_bottomheight(τ,A,Pars);   # save h to hdf5 file
    end
end