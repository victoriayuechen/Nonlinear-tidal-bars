using Parameters

"""
spatially uniform external pressure gradient
Fᵤ and Fᵥ
"""
@inline function external_press_grad(t,Pars)
    @unpack U,H,σ,cd,f,r,δ = Pars
    @unpack external_forcing_on,friction_on,linear_friction_on,coriolis_on,C∞_mode = Pars
    
    u0   = U * sin(σ * t)
    ∂ₜu0 = σ * U * cos(σ * t)

    Fᵤ,Fᵥ = zero(U), zero(U)
    if external_forcing_on
        Fᵤ += ∂ₜu0
        if coriolis_on
            Fᵥ += f * u0
        end
        if friction_on
            if linear_friction_on
                Fᵤ += r/H * u0
            else
                if C∞_mode
                    Fᵤ += cd/H * sqrt(u0^2+δ) * u0  # add small parameter δ to keep it differentiable when u0 passes trough zero
                else 
                    Fᵤ += cd/H * abs(u0) * u0
                end
            end
        end
    end
    return Fᵤ,Fᵥ
end

"""
friction term
frictionᵤ and frictionᵥ
"""
@inline function friction!(A,u,v,Pars)
    @unpack friction_on,linear_friction_on, C∞_mode = Pars
    if friction_on
        if linear_friction_on
            @unpack cd,U,H,r = Pars
            @. A.frictionᵤ = -r/A.Dᵤ * u
            @. A.frictionᵥ = -r/A.Dᵥ * v
        else
            @unpack cd = Pars
            if C∞_mode    
                @unpack cd, δ = Pars
                @. A.frictionᵤ = -cd/A.Dᵤ * sqrt(A.Ûᵤ² + δ) * u
                @. A.frictionᵥ = -cd/A.Dᵥ * sqrt(A.Ûᵥ² + δ) * v
            else
                @unpack cd = Pars
                @. A.frictionᵤ = -cd/A.Dᵤ * sqrt(A.Ûᵤ²) * u
                @. A.frictionᵥ = -cd/A.Dᵥ * sqrt(A.Ûᵥ²) * v
            end

        end
    end
end

"""
pressure gradient term
pressgradᵤ and pressgradᵥ
"""
@inline function pressgrad!(A,ζ,Pars)
    @unpack g, Δx, Δy, pressgrad_on = Pars
    if pressgrad_on
        @. A.pressgradᵤ = -g * A.dζdx
        @. A.pressgradᵥ = -g * A.dζdy
    end
end

"""
advection term
advectionᵤ and advectionᵥ
"""
@inline function advection!(A,u,v,Pars)
    @unpack Δx, Δy, advection_on = Pars
    if advection_on
        @. A.advectionᵤ = - (  u  * A.dudx + A.vᵤ * A.dudy)
        @. A.advectionᵥ = - (A.uᵥ * A.dvdx +   v  * A.dvdy)
    end
end

"""
continuity term
continuityᵤ and continuityᵥ
"""
@inline function continuity!(A,u,v,Pars)
    @unpack Δx,Δy, continuity_on = Pars
    if continuity_on
        @. A.Dᵤu = -A.Dᵤ * u
        @. A.Dᵥv = -A.Dᵥ * v
        ∂x_CU!(A.continuityᵤ, A.Dᵤu , Δx) # calculate -d(Du)dx on C grid
        ∂y_CV!(A.continuityᵥ, A.Dᵥv , Δy) # calculate -d(Dv)dy on C grid
    end
end

"""
Coriolis term
coriolisᵤ and coriolisᵥ
"""
@inline function coriolis!(A,Pars)
    @unpack f,coriolis_on = Pars
    if coriolis_on
        @. A.coriolisᵤ =  f * A.vᵤ
        @. A.coriolisᵥ = -f * A.uᵥ
    end
end

"""
calculate depth
D
"""
@inline function D!(A,ζ,Pars)
    @unpack H = Pars
    @. A.D = H + ζ - A.h
    @assert(
        minimum(A.D) > zero(eltype(A.D)),
        "Negative depth D = " * string(minimum(A.D)) * " at " * string(argmin(A.D))
        )
end

"""
calculate speed² at U, V, C grid
Ûᵤ, Ûᵥ, Û
"""
@inline function Û2!(A,u,v,Pars)
    @unpack nx,ny = Pars
    @. A.Ûᵤ² = u^2    + A.vᵤ^2
    @. A.Ûᵥ² = A.uᵥ^2 + v^2

    
    @inbounds begin
        for j in 1:ny-1
            for i in 2:nx-1
                A.Û²[i,j]  = inv(4) * (A.Ûᵤ²[i,j] + A.Ûᵤ²[i+1,j] + A.Ûᵥ²[i,j] + A.Ûᵥ²[i,j+1])
            end
        end
    end
    periodic_bc!(A.Û²)
end

"""
calculate effective velocity at U, V and C grid
u_eff², u_effᵤ², u_effᵥ²
"""
@inline function u_effective!(A,Pars)
    @unpack u_orbital, H = Pars
    if iszero(u_orbital)
        @. A.u_eff²  = A.Û²  
        @. A.u_effᵤ² = A.Ûᵤ² 
        @. A.u_effᵥ² = A.Ûᵥ²
    else     
        @. A.u_eff²  = A.Û²  + 0.5 * (u_orbital*H/A.D )^2
        @. A.u_effᵤ² = A.Ûᵤ² + 0.5 * (u_orbital*H/A.Dᵤ)^2
        @. A.u_effᵥ² = A.Ûᵥ² + 0.5 * (u_orbital*H/A.Dᵥ)^2
    end
end

"""
calculate ∇h = (A.dhdx, A.dhdy), 
dhdx is calculated on the U grid, 
dhdy on the V grid.
"""
@inline function ∇h!(A,Pars)
    @unpack Δx, Δy = Pars
    @unpack concentration_on, c_diffusion_on, c_diffusion_∇h_term_on, bedload_on, bl_diffusion_on = Pars
    if (concentration_on && c_diffusion_on && c_diffusion_∇h_term_on) || (bedload_on && bl_diffusion_on)
        ∂x_UC!(A.dhdx, A.h, Δx) # calculate dhdx on U grid
        ∂y_VC!(A.dhdy, A.h, Δy) # calculate dhdy on V grid
    end
end


"""
Do precalculations for rhs swec
"""
function precalculations!(A,ζ,u,v,c,Pars)
    @unpack continuity_on, advection_on, friction_on, pressgrad_on, coriolis_on, linear_friction_on = Pars
    @unpack concentration_on, c_advection_on, c_diffusion_on = Pars
    @unpack c_diffusion_∇ζ_term_on, c_diffusion_∇c_term_on, c_diffusion_∇h_term_on = Pars
    @unpack bedload_on, bl_advection_on, bl_diffusion_on = Pars

    @unpack Δx, Δy = Pars

    D!(A,ζ,Pars)              # calculate depth
    interpolate_CU!(A.Dᵤ,A.D) # interpolate D to U grid
    interpolate_CV!(A.Dᵥ,A.D) # interpolate D to V grid
    interpolate_VU!(A.vᵤ,v)   # interpolate v to U grid
    interpolate_UV!(A.uᵥ,u)   # interpolate u to V grid
    interpolate_CU!(A.cᵤ,c)   # interpolate c to U grid
    interpolate_CV!(A.cᵥ,c)   # interpolate c to V grid
    Û2!(A,u,v,Pars)           # calculate current speed² on U and V and C grid
    u_effective!(A,Pars)      # effective velocity
    c_top_bottom!(A,c,Pars)   # calculate consentration at top and bottom of water column
    Hu!(A,Pars)               # calculate H(uₑ²-u_crit^2)(uₑ²-u_crit^2) on C, U and V grid

    if advection_on
        ∂x_UU!(A.dudx, u, Δx) # calculate dudx on U grid
        ∂y_UU!(A.dudy, u, Δy) # calculate dudy on U grid
        ∂x_VV!(A.dvdx, v, Δx) # calculate dvdx on V grid
        ∂y_VV!(A.dvdy, v, Δy) # calculate dvdy on V grid
    end

    if concentration_on && c_diffusion_on && c_diffusion_∇c_term_on
        ∂x_UC!(A.dcdx, c, Δx)   # calculate dcdx on U grid
        ∂y_VC!(A.dcdy, c, Δy)   # calculate dcdy on V grid
    end
    if pressgrad_on || (concentration_on && c_diffusion_on && c_diffusion_∇ζ_term_on)
        ∂x_UC!(A.dζdx, ζ, Δx)   # calculate dζdx on U grid
        ∂y_VC!(A.dζdy, ζ, Δy)   # calculate dζdy on V grid
    end
end


