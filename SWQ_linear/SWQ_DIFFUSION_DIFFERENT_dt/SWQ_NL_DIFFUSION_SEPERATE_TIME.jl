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

function SWES_Diff_dt(Initial_solutions, Ω, Y, X, dΩ, dΓ, Param)
    @unpack dt1, dt2,dt3, Tstart, Tend1, Tend2,Tend3, theta1, theta2, theta3, dir = Param
    x1 = solver(Initial_solutions, Y, X, dt1, Tstart, Tend1, theta1, dΩ, dΓ, false,Param)
    @time writing_output(dir, x1, Ω, Tend1)
    

    # x2 = solver(xn,Y,X,dt2,Tend1,Tend2, theta2 ,dΩ,dΓ, false,Param)
    # @time writing_output(dir, x2, Ω, Tend2)
    # x3 = solver(xn,Y,X,dt3,Tend2,Tend2+5000, theta3 ,dΩ,dΓ, false,Param)
    # @time writing_output(dir, x3, Ω, Tend2+5000)
    # x4 = solver(xn,Y,X,40,Tend2+5000,Tend2+10000, theta3 ,dΩ,dΓ, false,Param)
    # @time writing_output(dir, x4, Ω, Tend2+10000)
end

# function SWES_Diff_dt(Initial_solutions, Ω, Y, X, dΩ, dΓ, Param)
#     @unpack dt1, dt2, Tstart, Tend1, Tend2, theta1, theta2, dir = Param
#     x1 = solver(Initial_solutions, Y, X, 5, Tstart, Tend1, theta1, dΩ, dΓ, false,Param)
#     @time writing_output(dir, x1, Ω, Tend1)
#     x2 = solver(Initial_solutions,Y,X,50,Tstart,Tend1, theta1 ,dΩ,dΓ, false,Param)
#     @time writing_output(dir, x2, Ω, Tend1)
#     x3 = solver(Initial_solutions,Y,X,10,Tstart,Tend1, theta1 ,dΩ,dΓ, false ,Param)
#     @time writing_output(dir, x3, Ω, Tend1)
#     x4 = solver(Initial_solutions,Y,X,20,Tstart,Tend1, theta1 ,dΩ,dΓ, false,Param)
#     @time writing_output(dir, x4, Ω, Tend1)
#     x5 = solver(Initial_solutions,Y,X,30,Tstart,Tend1, theta1 ,dΩ,dΓ, false,Param)
#     @time writing_output(dir, x5, Ω, Tend1)
#     x6 = solver(Initial_solutions,Y,X,40,Tstart,Tend1, theta1 ,dΩ,dΓ, false,Param)
#     @time writing_output(dir, x6, Ω, Tend1)
#     x7 = solver(Initial_solutions,Y,X,60,Tstart,Tend1, theta1 ,dΩ,dΓ, false,Param)
#     @time writing_output(dir, x7, Ω, Tend1)
#     x8 = solver(Initial_solutions,Y,X,100,Tstart,Tend1, theta1 ,dΩ,dΓ, false,Param)
#     @time writing_output(dir, x8, Ω, Tend1)
    
# end

# # Turns out 5 seconds is most optimal
# # After Tend1 = 10000 s 10 seconds is most optimal so dt2 = 10 s Speed increases by 17% Now Try to find it after 10000 seconds
## Na 35000 s is dt = 40 s most optimal

#=
Tend1 = 20000 --> Checken --> dt = 5
Tend2 = 35000/30000 --> Checken --> dt=10
Tend3 = undefined --> dt = 40
=#

Param = Parameter(dir = "test4",
    dt1 = 10,
    Tend1 = 44700
    )
Param2 = Parameter(dt1 = 5,
    dt2 = 30,
    x_points = 25,
    y_points = 5,
    Tend1 = 10
)

#Only first run
Initial_solutions, Ω, Y, X, dΩ, dΓ = setup(Param)
# SWES_Diff_dt(Initial_solutions, Ω, Y, X, dΩ, dΓ,Param2)


SWES_Diff_dt(Initial_solutions, Ω, Y, X, dΩ, dΓ,Param)

