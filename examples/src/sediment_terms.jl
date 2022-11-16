using Parameters

includet("utils.jl")

"""
Erosion term
erosion
"""
function erosion!(A,Pars)
    @unpack αs, concentration_on = Pars
    if concentration_on
        @. A.erosion = αs * A.Hu
    end
end

"""
Deposition term
deposition
"""
function deposition!(A,Pars)
    @unpack wₛ, concentration_on = Pars
    if concentration_on
        @. A.deposition = -wₛ * A.c_bottom
    end
end

"""
suspended sediment advection term
qs_advectionᵤ and qs_advectionᵥ
"""
function qs_advection!(A,u,v,Pars)
    @unpack concentration_on, c_advection_on = Pars
    if concentration_on & c_advection_on
        @. A.qs_advectionᵤ = u * A.cᵤ
        @. A.qs_advectionᵥ = v * A.cᵥ
    end
end

"""
suspended sediment diffusion term
qs_diffᵤ and qs_diffᵥ
"""
function qs_diff!(A,ζ,c,Pars)
    @unpack Δx, Δy, μ, concentration_on, c_diffusion_on = Pars
    @unpack c_diffusion_∇c_term_on, c_diffusion_∇ζ_term_on, c_diffusion_∇h_term_on  = Pars
    if concentration_on & c_diffusion_on
        if c_diffusion_∇c_term_on
            ∂x_UC!(A.dcdx, c, Δx)   # calculate dcdx on U grid
            ∂y_VC!(A.dcdy, c, Δy)   # calculate dcdy on V grid
        end
        if c_diffusion_∇h_term_on
            ∂x_UC!(A.dhdx, A.h, Δx) # calculate dhdx on U grid
            ∂y_VC!(A.dhdy, A.h, Δy) # calculate dhdy on V grid
        end
        if c_diffusion_∇ζ_term_on
            ∂x_UC!(A.dζdx, ζ, Δx)   # calculate dζdx on U grid
            ∂y_VC!(A.dζdy, ζ, Δy)   # calculate dζdy on V grid
        end
        if c_diffusion_∇c_term_on && c_diffusion_∇h_term_on && c_diffusion_∇ζ_term_on
            @. A.qs_diffᵤ = -μ * (A.dcdx + A.c_bottomᵤ * A.dhdx + A.c_topᵤ * A.dζdx)
            @. A.qs_diffᵥ = -μ * (A.dcdy + A.c_bottomᵥ * A.dhdy + A.c_topᵥ * A.dζdy)
        elseif c_diffusion_∇c_term_on && !c_diffusion_∇h_term_on && !c_diffusion_∇ζ_term_on
            @. A.qs_diffᵤ = -μ * A.dcdx
            @. A.qs_diffᵥ = -μ * A.dcdy
        elseif !c_diffusion_∇c_term_on && c_diffusion_∇h_term_on && !c_diffusion_∇ζ_term_on
            @. A.qs_diffᵤ = -μ * A.c_bottomᵤ * A.dhdx
            @. A.qs_diffᵥ = -μ * A.c_bottomᵥ * A.dhdy
        elseif !c_diffusion_∇c_term_on && !c_diffusion_∇h_term_on && c_diffusion_∇ζ_term_on
            @. A.qs_diffᵤ = -μ * A.c_topᵤ * A.dζdx
            @. A.qs_diffᵥ = -μ * A.c_topᵥ * A.dζdy
        elseif c_diffusion_∇c_term_on && c_diffusion_∇h_term_on && !c_diffusion_∇ζ_term_on
            @. A.qs_diffᵤ = -μ * (A.dcdx + A.c_bottomᵤ * A.dhdx)
            @. A.qs_diffᵥ = -μ * (A.dcdy + A.c_bottomᵥ * A.dhdy)
        elseif c_diffusion_∇c_term_on && !c_diffusion_∇h_term_on && c_diffusion_∇ζ_term_on
            @. A.qs_diffᵤ = -μ * (A.dcdx + A.c_topᵤ * A.dζdx)
            @. A.qs_diffᵥ = -μ * (A.dcdy + A.c_topᵥ * A.dζdy)
        elseif !c_diffusion_∇c_term_on && c_diffusion_∇h_term_on && c_diffusion_∇ζ_term_on
            @. A.qs_diffᵤ = -μ * (A.c_bottomᵤ * A.dhdx + A.c_topᵤ * A.dζdx)
            @. A.qs_diffᵥ = -μ * (A.c_bottomᵥ * A.dhdy + A.c_topᵥ * A.dζdy)
        end
    end
end

"""
suspended sediment transport term
qsᵤ and qsᵥ
"""
function q⃗s!(A,Pars)
    @unpack concentration_on, c_advection_on, c_diffusion_on = Pars
    if concentration_on
        if c_advection_on && c_diffusion_on
            @. A.qsᵤ = A.qs_advectionᵤ + A.qs_diffᵤ
            @. A.qsᵥ = A.qs_advectionᵥ + A.qs_diffᵥ
        end
        if c_advection_on && !c_diffusion_on
            @. A.qsᵤ = A.qs_advectionᵤ
            @. A.qsᵥ = A.qs_advectionᵥ
        end
        if !c_advection_on && c_diffusion_on
            @. A.qsᵤ = A.qs_diffᵤ
            @. A.qsᵥ = A.qs_diffᵥ
        end
        noFlux_bc!(A.qsᵥ)
    end
end

"""
divergence of suspended sediment term
div_q⃗ₛ
"""
function div_q⃗s!(A,Pars)
    @unpack Δx,Δy, concentration_on = Pars
    if concentration_on
        ∂x_CU!(A.dqsᵤdx, A.qsᵤ, Δx) # calculate dqsᵤdx on C grid
        ∂y_CV!(A.dqsᵥdy, A.qsᵥ, Δy) # calculate dqsᵥdy on C grid
        @. A.div_q⃗s = A.dqsᵤdx + A.dqsᵥdy
    end
end

# # bedload transport -----------------------------------
"""
bedload sediment advection term
qb_advectionᵤ and qb_advectionᵥ
"""
function qb_advection!(A,u,v,Pars)
    @unpack αb, u_critical, bedload_on, bl_advection_on = Pars
    if bedload_on && bl_advection_on
        @. A.qb_advectionᵤ = αb * A.Huᵤ * u
        @. A.qb_advectionᵥ = αb * A.Huᵥ * v
    end
end

"""
bedload diffusion term
qb_diffᵤ and qb_diffᵥ
"""
function qb_diff!(A,Pars)
    @unpack Δx, Δy, Λ, αb, bedload_on, bl_diffusion_on = Pars
    if bedload_on && bl_diffusion_on
        @. A.qb_diffᵤ = -αb * A.Huᵤ * Λ * sqrt(A.u_effᵤ²) * A.dhdx
        @. A.qb_diffᵥ = -αb * A.Huᵥ * Λ * sqrt(A.u_effᵥ²) * A.dhdy
    end
end

"""
bedload sediment transport term
qbᵤ and qbᵥ
"""
function q⃗b!(A,Pars)
    @unpack bedload_on, bl_advection_on, bl_diffusion_on = Pars
    if bedload_on
        if bl_advection_on && bl_diffusion_on
            @. A.qbᵤ = A.qb_advectionᵤ + A.qb_diffᵤ
            @. A.qbᵥ = A.qb_advectionᵥ + A.qb_diffᵥ
        end
        if bl_advection_on && !bl_diffusion_on
            @. A.qbᵤ = A.qb_advectionᵤ
            @. A.qbᵥ = A.qb_advectionᵥ
        end
        if !bl_advection_on && bl_diffusion_on
            @. A.qbᵤ = A.qb_diffᵤ
            @. A.qbᵥ = A.qb_diffᵥ
        end
        noFlux_bc!(A.qbᵥ)
    end
end

"""
divergence of bedload sediment term
div_q⃗b
"""
function div_q⃗b!(A,Pars)
    @unpack Δx,Δy, bedload_on = Pars
    if bedload_on
        ∂x_CU!(A.dqbᵤdx, A.qbᵤ, Δx) # calculate dqbᵤdx on C grid
        ∂y_CV!(A.dqbᵥdy, A.qbᵥ, Δy) # calculate dqbᵥdy on C grid
        @. A.div_q⃗b = A.dqbᵤdx + A.dqbᵥdy
    end
end
#------------------------------------------------------------

"""
concentration at bottom and top of water column
c_bottom and c_top
"""
function c_top_bottom!(A,c,Pars)
    @unpack wₛ, κᵥ, nx, ny = Pars
    @unpack concentration_on, c_diffusion_on, c_diffusion_∇h_term_on, c_diffusion_∇ζ_term_on = Pars
    @unpack c_top_bottom_dep_on_D = Pars
    if concentration_on
        if c_top_bottom_dep_on_D
            # if the exp in c_b and c_t depend on D (and hence on x,y)
            @. A.c_bottom  = wₛ * c/κᵥ * (1 - exp(-wₛ*A.D/κᵥ))^(-1)
            if c_diffusion_on
                if c_diffusion_∇h_term_on
                    @. A.c_bottomᵤ = wₛ * A.cᵤ/κᵥ * (1 - exp(-wₛ*A.Dᵤ/κᵥ))^(-1)
                    @. A.c_bottomᵥ = wₛ * A.cᵥ/κᵥ * (1 - exp(-wₛ*A.Dᵥ/κᵥ))^(-1)
                end
                if c_diffusion_∇ζ_term_on
                    @. A.c_topᵤ    = wₛ * A.cᵤ/κᵥ * (exp( wₛ*A.Dᵤ/κᵥ) - 1)^(-1)
                    @. A.c_topᵥ    = wₛ * A.cᵥ/κᵥ * (exp( wₛ*A.Dᵥ/κᵥ) - 1)^(-1)
                end
            end
        else
            # if the exp in c_b and c_t depends only on constant H
            @unpack H = Pars
            βb_wₛdivκᵥ = wₛ/κᵥ * (1 - exp(-wₛ*H/κᵥ))^(-1)
            @. A.c_bottom  = c * βb_wₛdivκᵥ
            if c_diffusion_on
                if c_diffusion_∇h_term_on
                    @. A.c_bottomᵤ = A.cᵤ * βb_wₛdivκᵥ
                    @. A.c_bottomᵥ = A.cᵥ * βb_wₛdivκᵥ
                end
                if c_diffusion_∇ζ_term_on
                    βt_wₛdivκᵥ = wₛ/κᵥ * (exp( wₛ* H/κᵥ) - 1)^(-1)
                    @. A.c_topᵤ    = A.cᵤ * βt_wₛdivκᵥ
                    @. A.c_topᵥ    = A.cᵥ * βt_wₛdivκᵥ
                end
            end

        end

        # this is a little bit faster than using exp:
        # for j in 1:ny-1
        #     for i in 2:nx-1
        #         A.c_bottom[i,j]  = inv(4) * (A.c_bottomᵤ[i,j] + A.c_bottomᵤ[i+1,j] + A.c_bottomᵥ[i,j] + A.c_bottomᵥ[i,j+1])
        #     end
        # end
        # periodic_bc!(A.c_bottom)
    end
end

"""
"""
function Hu!(A,Pars)
    @unpack u_critical, C∞_mode, concentration_on,bedload_on = Pars
    if concentration_on || bedload_on
        if iszero(u_critical)
            @. A.Hu  = A.u_eff² 
            @. A.Huᵤ = A.u_effᵤ²
            @. A.Huᵥ = A.u_effᵥ²
        elseif C∞_mode
            @. A.Hu  = Heaviside_C∞(A.u_eff²  - u_critical^2) * (A.u_eff²  - u_critical^2)
            @. A.Huᵤ = Heaviside_C∞(A.u_effᵤ² - u_critical^2) * (A.u_effᵤ² - u_critical^2)
            @. A.Huᵥ = Heaviside_C∞(A.u_effᵥ² - u_critical^2) * (A.u_effᵥ² - u_critical^2)
        else
            @. A.Hu  = Heaviside(A.u_eff²  - u_critical^2) * (A.u_eff²  - u_critical^2)
            @. A.Huᵤ = Heaviside(A.u_effᵤ² - u_critical^2) * (A.u_effᵤ² - u_critical^2)
            @. A.Huᵥ = Heaviside(A.u_effᵥ² - u_critical^2) * (A.u_effᵥ² - u_critical^2)
        end
    end
end