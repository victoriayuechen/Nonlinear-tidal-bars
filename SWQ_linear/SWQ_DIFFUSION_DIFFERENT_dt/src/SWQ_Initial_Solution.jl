using Gridap
using Parameters
includet("SWQ_Parameters.jl")

##''''''''''''''''Define initial solution''''''''''''''''##
function init_sol(X,Param)
    @unpack u_y,u_x,ζₙ = Param
    
    u₀(t) = VectorValue(u_x,u_y)
    ζ₀(t) = ζₙ 
    # ζ₀(x,t) = 8.2e-3-x[1]*8.2e-3/500
    # ζ₀(t::Real) = x->ζ₀(x,t)
    # u₀(x,t) = VectorValue( 8.1e-3*(exp(-((x[2]-7750)^2*0.00000013+(x[1]-500)^2*0.000015)) - exp(-((x[2]-2750)^2*0.00000013+(x[1]-500)^2*0.000015))),0.86)
    # u₀(t::Real) = x->u₀(x,t)

    uhn = interpolate_everywhere([u₀(0.0),ζ₀(0.0)],X(0.0))
    return uhn
end

function find_h(P,Param)
    @unpack B, L = Param
    h₀(x,t) = 0.1 * cos(π/B * x[2]) * cos(2π/L * x[1])      #Hepkema original function for h
    # h₀(x,t) = 0.0
    h₀(t::Real) = x->h₀(x,t)
    global h = interpolate_everywhere(h₀(0.0),P(0.0))
end
