using Parameters
includet("SWQ_Make_Model_APVM.jl")
includet("SWQ_Initial_Solution_APVM.jl")
includet("SWQ_Parameters_APVM.jl")


##"""""""""""""""Setup function""""""""""""""##
function setup(Param)    
    @time Ω, dΩ, dΓ, Y, X, P, R, S, U, V = Make_model(Param)
    
    find_h(P,Param)

    @time Initial_solutions = init_sol(X,P,R,S,U,V,Y,dΩ,Param)

    return Initial_solutions, Ω, Y, X, dΩ, dΓ, P
end