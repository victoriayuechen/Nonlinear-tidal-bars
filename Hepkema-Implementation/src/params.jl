using Parameters

# Parameters
@with_kw struct Params{T1 <: Real,T2 <: Integer} 

    # Physical parameters
    # --------------------------------------------------------------------
    # domain
    L::T1  = 10000.0           # length of channel
    B::T1  = 1000.0            # width of channel
    H::T1  = 3.0               # undisturbed depth of channel

    # hydrodynamics
    U::T1  = 0.5               # amplitude of background current
    cd::T1 = 2.5e-3            # drag coefficient
    r::T1  = 8/3π * cd * U     # linear friction coefficient (used if linear_friction_on = true)
    g::T1  = 9.81              # gravitational acceleration
    Ω::T1  = 2π/86164.0        # angular velocity of the Earth 
    φ::T1  = 50.0              # latitude in degrees N
    f::T1  = 2Ω*sin(φ*π/180)   # Coriolis parameter
    u_orbital::T1  = 0.25      # near bed wave orbital velocity
    u_critical::T1 = 0.3       # critical erosion velocity
    tidalcycle::T1 = 44700.0   # seconds in a tidalcyle
    σ::T1 = 2*π/tidalcycle     # tidal frequency
    
    # suspended sediment
    p̂::T1  = 0.4               # sediment porosity
    αs::T1 = 5e-6              # erosion constant
    μ::T1  = 10.0              # horizontal eddy diffusivity
    κᵥ::T1 = 1e-2              # vertical eddy diffusivity
    wₛ::T1 = 1.3e-2            # settling speed 
    
    # bedload
    αb::T1 = 3.0e-4            # advective bedload transport constant
    Λ::T1  = 2.0               # diffusive bedload transport constant
    # --------------------------------------------------------------------
    

    #Numerical parameters
    # --------------------------------------------------------------------
    #grid
    Δx::T1 = 100.0                  # along-channel spatial resolution
    Δy::T1 = 50.0                   # lateral spatial resolution
    nx_physical::T2 = Int(L/Δx)+1   # number of x gridpoints in physcial domain
    ny_physical::T2 = Int(B/Δy)+1   # number of y gridpoints in physcial domain
    nx::T2 = nx_physical+2          # number of x gridpoints in computational domain (+2 for periodic bc)
    ny::T2 = ny_physical            # number of y gridpoints in computational domain (one to many for U an C grid)
    
    # numerical timesteps
    Δt_num::T1 = 5.0                # in seconds; fixed Runge-Kutta timestep for swec
    Δτ_num::T1 = 6.048              # in ε sec;   fixed Runge-Kutta timestep for bed evolution; 6.048 is a week
    
    # rest
    δ::T1  = 1e-5                   # small parameter to make friction differentiable in h
    ε::T1  = 1e-5                   # scaling dhdt, so τ = εt
    sec_in_year::T1 = 3600*24*365.25 # number of seconds in a year

    # output, morphodynamics, units: ε * seconds
    τ_start::T1                 = 0.0         # start time of morpodynamical simulation (units ε * seconds)
    Δτ_out::T1                  = 5.0*Δτ_num  # timestep for output of morphodynamics; default = save output every 5 numerical timesteps (units ε * seconds)
    t_end_in_years::T1          = 1.0         # end time of morpodynamical simulation in years (it will simulate a bit longer to fit Δτ)
    nr_Δτ_steps::T1             = ceil(t_end_in_years  * sec_in_year * ε/ Δτ_num)  # number of numerical timesteps in morphodynamical simulation
    τ_end::T1                   = nr_Δτ_steps * Δτ_num          # actual end time of morpodynamical simulation (units ε * seconds)
    τ_span::Tuple{T1,T1}        = (τ_start,τ_end)               # morphodynamical simulation timespan
    save_moments_τ::Array{T1,1} = collect(τ_start:Δτ_out:τ_end) # moments to save/output the bottomheight
    
    # output, hydronamics
    t_start::T1                 = 0.0                           # start time of hydroodynamical simulation (units seconds)
    Δt_out::T1                  = 447.0                         # timestep for output of hydrodynamics (units seconds)
    t_end::T1                   = tidalcycle                    # end time of hydrodynamical simulation (units seconds)
    t_span::Tuple{T1,T1}        = (t_start,t_end)               # hydrodynamical simulation timespan
    save_moments_t::Array{T1,1} = collect(t_start:Δt_out:t_end) # moments to save/output the ζ, u, v, c, etc.
    # --------------------------------------------------------------------
    
    
    # Control parameters
    # ---------------------------------------------------------------------
    # output control
    output_path::String          = "output/test/"  # folder for output
    output_hydro::Bool           = false           # save ζ, u, v, c vs time and space to swec.h5
    output_residual::Bool        = false           # save ⟨ζ⟩, ⟨u⟩, ⟨v⟩, ⟨c⟩ vs space to residual_swec.h5
    output_morpho::Bool          = true            # save h vs time and space to bottomheight.h5

    # hydrodynamics
    continuity_on::Bool          = true            # calculate ∇⋅(Du⃗)
    advection_on::Bool           = true            # calculate (u⃗⋅∇)u⃗
    friction_on::Bool            = true            # calculate cd|u⃗|u⃗/D
    pressgrad_on::Bool           = true            # calculate g∇ζ
    coriolis_on::Bool            = true            # calculate -fv, fu
    external_forcing_on::Bool    = true            # calculate (Fᵤ,Fᵥ)
    linear_friction_on::Bool     = false           # if true: r/D u⃗, if false: quadratic friction: cd/D |u⃗| u⃗
    
    # suspended sediment
    concentration_on::Bool       = true            # calculate c and suspended sediment transport
    c_advection_on::Bool         = true            # calculate u⃗c
    c_diffusion_on::Bool         = true            # calculate -μ(∇C + c_top ∇ζ + c_bottom ∇h) or part of it
    c_diffusion_∇c_term_on::Bool = true            # calculate -μ∇C
    c_diffusion_∇h_term_on::Bool = true            # calculate -μ c_bottom ∇h 
    c_diffusion_∇ζ_term_on::Bool = true            # calculate -μ c_top ∇ζ
    c_top_bottom_dep_on_D::Bool  = true            # if false, replace D in the exp for H (saves calculations)
    
    # bedload
    bedload_on::Bool             = true            # calculate bedload sediment transport
    bl_advection_on::Bool        = true            # calculate  αb H(uₑ²-u_crit^2)(uₑ²-u_crit^2)u⃗
    bl_diffusion_on::Bool        = true            # calculate -αb H(uₑ²-u_crit^2)(uₑ²-u_crit^2)Λuₑ∇h

    # rest
    C∞_mode::Bool                = false           # smooth but slower mode, for jacobian calculations
    # ---------------------------------------------------------------------
end


