
using Revise 
using Gridap
includet("../linear_SWE.jl")
using .Linear_SWE_solver
using GridapGmsh

function ζ₀((x,y))
    h =0.05*exp(-0.01*(x-25)^2 -0.01*(y-25)^2)
    h
end


function u₀((x,y))
    u = VectorValue(0.0,0.0)
    u
end

order = 0
degree = 4
model = GmshDiscreteModel("swe-solver/meshes/100x100periodic2.msh")

run_linear_SWE(order,degree,ζ₀,u₀,model)