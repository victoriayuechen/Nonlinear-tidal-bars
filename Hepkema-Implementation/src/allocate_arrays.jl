"""
Struct with all the allocated arrays
"""
struct MyArrays{T <: AbstractMatrix, T2<:AbstractArray}
    
    ic::T2              # initial conditions ζ,u,v,c
    h::T                # (initial) bottom height vs space

    # hydrodynamic terms (all on right hand side, hence the minuses)
    continuityᵤ::T      # -∂(Dᵤu)/∂x           on C grid
    continuityᵥ::T      # -∂(Dᵥv)/∂y           on C grid
    frictionᵤ::T        # cd|u⃗|u/Dᵤ            on U grid
    frictionᵥ::T        # cd|u⃗|v/Dᵥ            on V grid
    pressgradᵤ::T       # -g ∂ζ/∂x             on U grid 
    pressgradᵥ::T       # -g ∂ζ/∂y             on V grid 
    advectionᵤ::T       # -(u ∂u/∂x + vᵤ∂u/∂y) on U grid
    advectionᵥ::T       # -(uᵥ∂v/∂x + v ∂v/∂y) on V grid
    coriolisᵤ::T        #  fv                  on U grid
    coriolisᵥ::T        # -fu                  on V grid

    # pre calculations
    D::T       # Depth D = H + ζ - h
    Dᵤ::T      # Depth on U grid
    Dᵥ::T      # Depth on V grid
    Dᵤu::T     # Depth on U grid Dᵤ times u
    Dᵥv::T     # Depth on V grid Dᵥ times v
    uᵥ::T      # velocity in x direction on V grid
    vᵤ::T      # velocity in y direction on U grid
    Ûᵤ²::T     # current speed² u²+v² on U grid
    Ûᵥ²::T     # current speed² u²+v² on V grid
    Û²::T      # current speed² u²+v² on C grid
    
    u_eff²::T  # effective velocity on C grid
    u_effᵤ²::T # effective velocity on U grid
    u_effᵥ²::T # effective velocity on V grid
    Hu::T      # H(uₑ²-u_crit^2)(uₑ²-u_crit^2) on C grid
    Huᵤ::T     # H(uₑ²-u_crit^2)(uₑ²-u_crit^2) on U grid
    Huᵥ::T     # H(uₑ²-u_crit^2)(uₑ²-u_crit^2) on V grid

    dudx::T    #dudx on U grid    
    dudy::T    #dudy on U grid
    dvdx::T    #dvdx on V grid
    dvdy::T    #dvdy on V grid
    dcdx::T    #dcdx on U grid
    dcdy::T    #dcdy on V grid
    dζdx::T    #dζdx on U grid
    dζdy::T    #dζdy on V grid
    dhdx::T    #dhdx on U grid
    dhdy::T    #dhdy on V grid

    # suspended sediment concentration
    cᵤ::T             # concentration c on U grid
    cᵥ::T             # concentration c on V grid
    c_bottomᵤ::T      # concentration at bottom on U grid
    c_bottomᵥ::T      # concentration at bottom on V grid
    c_bottom::T       # concentration at bottom on C grid
    c_topᵤ::T         # concentration at top water column on U grid
    c_topᵥ::T         # concentration at top water column on V grid
    erosion::T        # erosion of suspended sediment
    deposition::T     # deposition of suspended sediment
    qs_advectionᵤ::T  # advection of suspended sediment in x direction
    qs_advectionᵥ::T  # advection of suspended sediment in y direction
    qs_diffᵤ::T       # diffusion of suspended sediment in x direction
    qs_diffᵥ::T       # diffusion of suspended sediment in y direction
    qsᵤ::T            # suspended sediment transport in x direction
    qsᵥ::T            # suspended sediment transport in y direction
    dqsᵤdx::T         # derivative wrt x of qsᵤ
    dqsᵥdy::T         # derivative wrt y of qsᵥ
    div_q⃗s::T         # divergence of suspended sediment
    
    # bedload
    qb_advectionᵤ::T  # advection of bedload sediment in x direction
    qb_advectionᵥ::T  # advection of bedload sediment in y direction
    qb_diffᵤ::T       # diffusion of bedload sediment in x direction
    qb_diffᵥ::T       # diffusion of bedload sediment in y direction
    qbᵤ::T            # bedload sediment transport in x direction
    qbᵥ::T            # bedload sediment transport in y direction
    dqbᵤdx::T         # derivative wrt x of qbᵤ
    dqbᵥdy::T         # derivative wrt y of qbᵥ
    div_q⃗b::T         # divergence of bedload sediment transport

    #total sediment transport
    avg_div_q⃗::T      # tidal average of total sediment transport divergence
end

"""
overwrite show for MyArrays
"""
Base.show(io::IO, myarrays::MyArrays) = dump(myarrays)

"""
function that returns the arrays with size nx,ny
"""
function allocateArrays(nx,ny)
    ns = length(fieldnames(MyArrays)) # length of struct
    MyArrays(
        zeros(nx,ny,4),
        [zeros(nx,ny) for _ in 2:ns]...
        )
end

"""
function that returns the arrays with size nx,ny of type T
"""
function allocateArrays(T,nx,ny)
    ns = length(fieldnames(MyArrays)) # length of struct
    MyArrays(
        zeros(T,nx,ny,4),
        [zeros(T,nx,ny) for _ in 2:ns]...
        )
end