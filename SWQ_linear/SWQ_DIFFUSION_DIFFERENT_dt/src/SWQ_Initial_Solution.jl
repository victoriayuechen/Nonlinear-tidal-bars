using Gridap


##''''''''''''''''Define initial solution''''''''''''''''##
function init_sol(u₀,ζ₀,X)
    uhn = interpolate_everywhere([u₀(0.0),ζ₀(0.0)],X(0.0))
    return uhn
end

function find_h(h₀,P)
    h = interpolate_everywhere(h₀(0.0),P(0.0))
    return h
end
