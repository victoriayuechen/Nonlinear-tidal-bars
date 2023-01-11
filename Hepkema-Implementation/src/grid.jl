using Parameters

"""
struct for the grid
"""
struct Grid{T <: AbstractVector}
    x::T
    y::T
    xᵤ::T
    yᵤ::T
    xᵥ::T
    yᵥ::T
end

"""
return a Grid struct with x,y,xᵤ,yᵤ,xᵥ,yᵥ
"""
function Grid(Pars)
    @unpack Δx, Δy, L, B, nx, ny = Pars
    x  = collect(range(-Δx,  stop = L+Δx,   length = nx))
    y  = collect(range(Δy/2, stop = B-Δy/2, length = ny-1))
    xᵤ = collect(range(-3Δx/2,  stop = L+Δx/2,   length = nx))
    yᵤ = collect(range(Δy/2, stop = B-Δy/2, length = ny-1))
    xᵥ = collect(range(-Δx,  stop = L+Δx,   length = nx))
    yᵥ = collect(range(0, stop = B, length = ny))
    Grid(x,y,xᵤ,yᵤ,xᵥ,yᵥ)
end
