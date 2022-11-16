using Parameters, Random

includet("utils.jl")
includet("read_output.jl")

# initial bottomheights

"""
initialize bottom height with random perturbations between -ampl and +ampl
with seed
"""
function h₀_random(ampl,Pars,seed=1234)
    @unpack nx,ny = Pars
    Random.seed!(seed);
    h₀ = ampl .* (2 .* rand(Float64,(nx,ny)) .- 1);
    zeroGradient_bc!(h₀)
    periodic_bc!(h₀)
    subtract_mean!(h₀,Pars);
    return h₀
end

"""
initialize bottom height with:
    h₀(x,y) = ampl * cos(n π/B y) cos(m 2π/L x)
"""
function h₀_coscos(ampl,n,m,Pars,myGrid)
    @unpack nx,ny,B,L = Pars;
    h₀ = zeros(nx,ny);
    for j in 1:ny-1
        for i in 2:nx-1
            h₀[i,j] = ampl * cos(n * π/B * myGrid.y[j]) * cos(m * 2π/L * myGrid.x[i])
        end
    end
    periodic_bc!(h₀)
    return h₀
end

"""
initialize bottom height with bump in the middle:
    h₀(x,y) = ampl * exp(-(x -L/2)^2 /(2 st ^2) )
"""
function h₀_bump(ampl,st,Pars,myGrid)
    @unpack nx,ny, L = Pars;
    h₀ = zeros(nx,ny);
    for j in 1:ny-1
        for i in 2:nx-1
            h₀[i,j] = ampl * exp(-(myGrid.x[i]-L/2)^2 /(2*st^2) )
        end
    end
    periodic_bc!(h₀)
    return h₀
end

"""
initialize bottom height with output of an other run
"""
function h₀_output(path;tn=-1)
    t, h₀ = read_bottomheight(path,physical_domain=false);
    if tn == -1
        return h₀[:,:,end]
    else
        return h₀[:,:,tn]
    end
end
