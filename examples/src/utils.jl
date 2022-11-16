using Statistics, Printf


# Boundary conditions --------------------------
"""
apply periodic boundary condition
"""
@inline function periodic_bc!(array)
    nx,ny = size(array)
    @inbounds for j in 1:ny
        array[1,j]   = array[end-1,j]
        array[end,j] = array[2,j]
    end
end
@inline function periodic_bc!(arrays...)
    for a in arrays
        periodic_bc!(a)
    end
end

"""
no flux boundary condition
at V grid
"""
@inline function noFlux_bc!(array)
    nx,ny = size(array)
    @inbounds for i in 1:nx
        array[i,1]   = zero(eltype(array))
        array[i,end] = zero(eltype(array))
    end
end
@inline function noFlux_bc!(arrays...)
    for a in arrays
        noFlux_bc!(a)
    end
end

"""
vanishing gradient boundary condition
setting the element closest to the boundaries (on C grid)
such that the derivative is zero (up to second order) at the boundaries (V grid)
"""
@inline function zeroGradient_bc!(array)
    nx,ny = size(array)
    @inbounds for i in 1:nx
        array[i,1]     = 1.5 * array[i,2]     - 0.5 * array[i,3]
        array[i,end-1] = 1.5 * array[i,end-2] - 0.5 * array[i,end-3]
    end
end
@inline function zeroGradient_bc!(array...)
    for a in array
        zeroGradient_bc!(a)
    end
end
#---------------------------------------------


# Dervatives ---------------------------------
"""
derivative wrt x, gradx on C grid, array on U grid
"""
@inline function ∂x_CU!(gradx, array, Δx)
    nx,ny = size(array)
    @inbounds begin
        for j in 1:ny-1
            for i in 2:nx-1
                gradx[i,j] = inv(Δx) * (array[i+1,j] - array[i,j])
            end
        end
    end
    periodic_bc!(gradx)
end

"""
derivative wrt y, gradx on C grid, array on V grid
"""
@inline function ∂y_CV!(grady, array, Δy)
    nx,ny = size(array)
    @inbounds begin
        for j in 1:ny-1
            for i in 2:nx-1
                grady[i,j] = inv(Δy) * (array[i,j+1] - array[i,j])
            end
        end
    end
    periodic_bc!(grady)
end

"""
derivative wrt x, gradx on U grid, array on C grid
"""
@inline function ∂x_UC!(gradx, array, Δx)
    nx,ny = size(array)
    @inbounds begin
        for j in 1:ny-1
            for i in 2:nx-1
                gradx[i,j] = inv(Δx) * (array[i,j] - array[i-1,j])
            end
        end
    end
    periodic_bc!(gradx)
end

"""
derivative wrt y, gradx on V grid, array on C grid
"""
@inline function ∂y_VC!(grady, array, Δy)
    nx,ny = size(array)
    @inbounds begin 
        for j in 2:ny-1
            for i in 2:nx-1
                grady[i,j] = inv(Δy) * (array[i,j] - array[i,j-1])
            end
        end
    end
    #not needed because v = 0 at the boundaries anyway
    @inbounds for i in 2:nx-1
        grady[i,1]   = 0
        grady[i,end] = 0
    end
    periodic_bc!(grady)
end

"""
derivative wrt x, gradx on U grid, array on U grid
"""
@inline function ∂x_UU!(gradx, array, Δx)
    nx,ny = size(array)
    @inbounds begin
        for j in 1:ny-1
            for i in 2:nx-1
                gradx[i,j] = inv(2*Δx) * (array[i+1,j] - array[i-1,j])
            end
        end
    end
    periodic_bc!(gradx)
end

"""
derivative wrt x, gradx on V grid, array on V grid
"""
@inline function ∂x_VV!(gradx, array, Δx)
    nx,ny = size(array)
    @inbounds begin
        for j in 1:ny
            for i in 2:nx-1
                gradx[i,j] = inv(2*Δx) * (array[i+1,j] - array[i-1,j])
            end
        end
    end
    periodic_bc!(gradx)
end

"""
derivative wrt y, grady on U grid, array on U grid
"""
@inline function ∂y_UU!(grady, array, Δy)
    nx,ny = size(array)
    @inbounds begin 
        for j in 2:ny-2
            for i in 2:nx-1
                grady[i,j] = inv(2*Δy) * (array[i,j+1] - array[i,j-1])
            end
        end
    end
    @inbounds for i in 2:nx-1
        grady[i,1]    = inv(2*Δy) * (-3*array[i,1]    + 4*array[i,2]    - array[i,3])
        grady[i,ny-1] = inv(2*Δy) * ( 3*array[i,ny-1] - 4*array[i,ny-2] + array[i,ny-3])
    end
    periodic_bc!(grady)
end


"""
derivative wrt y, gradx on V grid, array on V grid
"""
@inline function ∂y_VV!(grady, array, Δy)
    nx,ny = size(array)
    @inbounds begin
        for j in 2:ny-1
            for i in 2:nx-1
                grady[i,j] = inv(2*Δy) * (array[i,j+1] - array[i,j-1])
            end
        end
    end
    #not needed because v = 0 at the boundaries anyway
    @inbounds for i in 2:nx-1
        grady[i,1]  = 0
        grady[i,ny] = 0
    end
    periodic_bc!(grady)
end

"""
derivative wrt x, gradx on C grid, array on C grid
"""
@inline function ∂x_CC!(gradx, array, Δx)
    nx,ny = size(array)
    @inbounds begin
        for j in 1:ny-1
            for i in 2:nx-1
                gradx[i,j] = inv(2*Δx) * (array[i+1,j] - array[i-1,j])
            end
        end
    end
    periodic_bc!(gradx)
end

"""
derivative wrt y, gradx on C grid, array on C grid
"""
@inline function ∂y_CC!(grady, array, Δy)
    nx,ny = size(array)
    @inbounds begin
        for j in 2:ny-2
            for i in 2:nx-1
                grady[i,j] = inv(2*Δy) * (array[i,j+1] - array[i,j-1])
            end
        end
    end
    @inbounds for i in 2:nx-1
        grady[i,1]    = inv(2*Δy) * (-3*array[i,1]    + 4*array[i,2]    - array[i,3])
        grady[i,ny-1] = inv(2*Δy) * ( 3*array[i,ny-1] - 4*array[i,ny-2] + array[i,ny-3])
    end
    periodic_bc!(grady)
end
#---------------------------------------------



# interpolation ------------------------------
"""
interpolate array on C grid to V grid
with neuman bc
"""
@inline function interpolate_CV!(Aᵥ,A)
    nx,ny = size(A)
    @inbounds begin
        for j in 2:ny-1
            for i in 2:nx-1
                Aᵥ[i,j] = inv(2) * (A[i,j-1] + A[i,j])
            end
        end
    end
    # at closed boundaries Aᵥ is extrapolated
    @inbounds for i in 2:nx-1
        Aᵥ[i,1]   = 1.875 * A[i,1]     - 1.25 * A[i,2]     + 0.375 * A[i,3]
        Aᵥ[i,end] = 1.875 * A[i,end-1] - 1.25 * A[i,end-2] + 0.375 * A[i,end-3]
    end
    periodic_bc!(Aᵥ)
end

"""
interpolate array on C grid to U grid
"""
@inline function interpolate_CU!(Aᵤ,A)
    nx,ny = size(A)
    @inbounds begin
        for j in 1:ny
            for i in 2:nx-1
                Aᵤ[i,j] = inv(2) * (A[i,j] + A[i-1,j])
            end
        end
    end
    periodic_bc!(Aᵤ)
end

"""
interpolate array on V grid to U grid
"""
@inline function interpolate_VU!(Aᵤ,Aᵥ)
    nx,ny = size(Aᵥ)
    @inbounds begin
        for j in 1:ny
            for i in 2:nx-1
                Aᵤ[i,j] = inv(4) * (Aᵥ[i,j] + Aᵥ[i-1,j] + Aᵥ[i,j+1] + Aᵥ[i-1,j+1])
            end
        end
    end
    periodic_bc!(Aᵤ)
end

"""
interpolate array on U grid to V grid
"""
@inline function interpolate_UV!(Aᵥ,Aᵤ)
    nx,ny = size(Aᵤ)
    @inbounds begin
        for j in 2:ny-1
            for i in 2:nx-1
                Aᵥ[i,j] = inv(4) * (Aᵤ[i,j] + Aᵤ[i+1,j] + Aᵤ[i,j-1] + Aᵤ[i+1,j-1])
            end
        end
    end
    # at closed boundaries Aᵤ is extrapolated
    @inbounds for i in 2:nx-1
        Aᵥ[i,1]   = 0.5*(
            1.875 * Aᵤ[i,1]   - 1.25 * Aᵤ[i,2]   + 0.375 * Aᵤ[i,3]
          + 1.875 * Aᵤ[i+1,1] - 1.25 * Aᵤ[i+1,2] + 0.375 * Aᵤ[i+1,3]
          )

        Aᵥ[i,end] = 0.5*(
            1.875 * Aᵤ[i,end-1]   - 1.25 * Aᵤ[i,end-2]   + 0.375 * Aᵤ[i,end-3]
          + 1.875 * Aᵤ[i+1,end-1] - 1.25 * Aᵤ[i+1,end-2] + 0.375 * Aᵤ[i+1,end-3]
        )
    end
    periodic_bc!(Aᵥ)
end

"""
interpolate array on U grid to C grid
"""
@inline function interpolate_UC!(A,Aᵤ)
    nx,ny = size(A)
    @inbounds begin
        for j in 1:ny
            for i in 2:nx-1
                A[i,j] = inv(2) * (Aᵤ[i,j] + Aᵤ[i+1,j])
            end
        end
    end
    periodic_bc!(A)
end

"""
interpolate array on V grid to C grid
"""
@inline function interpolate_VC!(A,Aᵥ)
    nx,ny = size(A)
    @inbounds begin
        for j in 1:ny-1
            for i in 2:nx-1
                A[i,j] = inv(2) * (Aᵥ[i,j] + Aᵥ[i,j+1])
            end
        end
    end
    periodic_bc!(A)
end
#---------------------------------------------

# Rest ---------------------------------------
"""
Heaviside (step) function
"""
function Heaviside(x::T) where T<:AbstractFloat
    if x > zero(T)
        return one(T)
    elseif iszero(x)
        return one(T) / 2
    else
        return zero(T)
    end
end


"""
smooth version of Heaviside (step) function
"""
function Heaviside_C∞(x,k=10)
    # return @fastmath 0.5 * (1 + tanh(k*x))
    return 0.5 * (1 + tanh(k*x))
end

"""
subtract mean value over the physcial domain 
"""
function subtract_mean!(h,Pars)
    @unpack nx,ny = Pars
    avg = zero(eltype(h))
    for j in 1:ny-1
        for i in 2:nx-1
            avg += h[i,j]
        end
    end
    avg /= (ny-1) * (nx-2)
    h .-= avg
end



"""
print progres of integration of swec, for example:
t =  3/ 62h   5%
"""
function print_progress(t,tmax)
    per = floor(t/tmax * 100)
    if per % 2 == 0
        t_h    = @sprintf("%6.2f", t   /3600)
        tmax_h = @sprintf("%6.2f", tmax/3600)
        print("t = " * t_h * "/ " * tmax_h )
        println(" h \t" * @sprintf("%3.0f",per) * "%")
    end
end


"""
print progres of integration of bed evolution, for example:
τ =    2.0/ 200.00 yr    1%
"""
function print_progress_morpho(τ,τmax,Pars)
    @unpack sec_in_year, ε = Pars
    per = floor(τ/τmax * 100)
    if per % 2 == 0
        t    = @sprintf("%6.2f", τ    /(ε * sec_in_year)) 
        tmax = @sprintf("%6.2f", τmax /(ε * sec_in_year))
        print("t = " * t * "/ " * tmax )
        println(" yr \t" * @sprintf("%3.0f",per) * "%")
    end
end

#---------------------------------------------


