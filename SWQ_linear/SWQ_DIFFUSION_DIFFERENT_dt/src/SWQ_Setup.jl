includet("SWQ_Make_Model.jl")
includet("SWQ_Initial_Solution.jl")


##"""""""""""""""Setup function""""""""""""""##
function setup()
    B = 1000
    L = 10000
    x_points = 20
    y_points = 100
    order = 1
    degree = 3
    
    @time Ω, dΩ, dΓ, Y, X, P = Make_model(B,L,x_points,y_points,order,degree)

    function u₀(t) 
        u = VectorValue(0.0,0.0)
        return u
    end

    function ζ₀(t)
        ζ = 0.0
        return ζ
    end

    h₀(x,t) = 0.1 * cos(π/B * x[1]) * cos(2π/L * x[2])      #Hepkema original function for h
    h₀(t::Real) = x->h₀(x,t)
    global h = find_h(h₀,P)
    # u₀(t) = VectorValue(0.0,0.0) evenly fast
    # ζ₀(t) = 0.0

    @time Initial_solutions = init_sol(u₀,ζ₀,X)

    return Initial_solutions, Ω, Y, X, dΩ, dΓ
end