using Gridap
using Parameters
includet("SWQ_Parameters.jl")

##''''''''''''''''Define initial solution''''''''''''''''##
function init_sol(X,Param)
    @unpack u_y,u_x,ζₙ = Param
    
    u₀(t) = VectorValue(u_y,u_x)
    ζ₀(t) = ζₙ

    uhn = interpolate_everywhere([u₀(0.0),ζ₀(0.0)],X(0.0))
    return uhn
end

function find_h(P,Param)
    @unpack B, L = Param
    h₀(x,t) = 0.1 * cos(π/B * x[1]) * cos(2π/L * x[2])      #Hepkema original function for h
    h₀(t::Real) = x->h₀(x,t)
    global h = interpolate_everywhere(h₀(0.0),P(0.0))
end
