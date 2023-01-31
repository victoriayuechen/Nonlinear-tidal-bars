using Parameters
includet("SWQ_Make_Model.jl")
includet("SWQ_Initial_Solution.jl")
includet("SWQ_Parameters.jl")


##"""""""""""""""Setup function""""""""""""""##
function setup(Param)    
    @time Ω, dΩ, dΓ, Y, X, P, nΓ = Make_model(Param)
    
    h = find_h(P,Param)

    @time Initial_solutions = init_sol(X,Param)

    return Initial_solutions, Ω, Y, X, dΩ, dΓ, P, h, nΓ
end