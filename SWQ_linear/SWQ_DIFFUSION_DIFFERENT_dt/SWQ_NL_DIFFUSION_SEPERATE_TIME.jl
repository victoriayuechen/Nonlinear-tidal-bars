using Pkg
Pkg.activate(".")

#Add packages
using Gridap
using WriteVTK
using LineSearches: BackTracking
using Gridap.TensorValues: meas
using Parameters

includet("src/SWQ_Solver.jl")
includet("src/SWQ_Write_Output.jl")
includet("src/SWQ_Setup.jl")
includet("src/SWQ_Parameters.jl")

function SWES_Diff_dt(Initial_solutions, Ω, Y, X, dΩ, dΓ, dt1, dt2, Tend1, Tend2, dir,Param)
    x1 = solver(Initial_solutions, Y, X, dt1, 0.0, Tend1, dΩ, dΓ, true,Param)
    writing_output(dir, x1, Ω, Tend1)
    x2 = solver(xn,Y,X,dt2,Tend1,Tend2,dΩ,dΓ, true,Param)
    writing_output(dir, x2, Ω, Tend2)
end
Param = Parameter()
Initial_solutions, Ω, Y, X, dΩ, dΓ = setup(Param)

SWES_Diff_dt(Initial_solutions, Ω, Y, X, dΩ, dΓ, 1, 10, 800, 1600, "./test",Param)
