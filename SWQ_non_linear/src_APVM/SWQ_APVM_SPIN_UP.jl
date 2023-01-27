using Pkg
Pkg.activate(".")

#Add packages
using Gridap
using Revise
using WriteVTK
using LineSearches: BackTracking
using Gridap.TensorValues: meas
using Parameters 
using CSV
using DataFrames
using DelimitedFiles
perp(u) = VectorValue(-u[2],u[1])

includet("SWQ_Solver_APVM.jl")
includet("SWQ_Write_Output_APVM.jl")
includet("SWQ_Setup_APVM.jl")
includet("SWQ_Parameters_APVM.jl")

function SWES_Diff_dt(Initial_solutions, Ω, Y, X, dΩ, dΓ, Param, P)
    @unpack dt1, dt2,dt3, Tstart, Tend1, Tend2,Tend3, theta1, theta2, theta3, dir = Param
    x1 = solver(Initial_solutions, Y, X, 100, Tstart, Tend1, theta1, dΩ, dΓ, false,Param)
    @time writing_output("APVM spinup 5010", "APVM dt_100",x1, Ω, Tend1, P, "zeta_NOspinup_APVM_order0_100_20_dt=5.csv")
    x2 = solver(xn, Y, X, dt1, Tstart, Tend1, theta1, dΩ, dΓ, false,Param)
    @time writing_output("APVM spin up dt_100 5010","APVM dt_100", x2, Ω, Tend1, P, "zeta_spinup_APVM_order0_100_20_dt=5.csv")
    # x3 = solver(xn, Y, X, 200, Tstart, Tend1, theta1, dΩ, dΓ, false,Param)
    # @time writing_output("APVM spin up dt_200","APVM dt_200", x3, Ω, Tend1+100)
    # x4 = solver(xn, Y, X, 400, Tstart, Tend1, theta1, dΩ, dΓ, false,Param)
    # @time writing_output("APVM spin up dt_400", "APVM dt_400",x4, Ω, Tend1+100)
    # x5 = solver(xn, Y, X, 800, Tstart, Tend1, theta1, dΩ, dΓ, false,Param)
    # @time writing_output("APVM spin up dt_800", "APVM dt_800",x5, Ω, Tend1+100)
    # x6 = solver(xn, Y, X, 1600, Tstart, Tend1, theta1, dΩ, dΓ, false,Param)
    # @time writing_output("APVM spin up dt_1600","APVM dt_1600", x6, Ω, Tend1+100)


    # x2 = solver(xn,Y,X,dt2,Tend1,Tend2, theta2 ,dΩ,dΓ, false,Param)
    # @time writing_output(dir, x2, Ω, Tend2)
    # x3 = solver(xn,Y,X,dt3,Tend2,Tend2+5000, theta3 ,dΩ,dΓ, false,Param)
    # @time writing_output(dir, x3, Ω, Tend2+5000)
    # x4 = solver(xn,Y,X,40,Tend2+5000,Tend2+10000, theta3 ,dΩ,dΓ, false,Param)
    # @time writing_output(dir, x4, Ω, Tend2+10000)
end


# # Turns out 5 seconds is most optimal
# # After Tend1 = 10000 s 10 seconds is most optimal so dt2 = 10 s Speed increases by 17% Now Try to find it after 10000 seconds
## Na 35000 s is dt = 40 s most optimal

#=
Tend1 = 20000 --> Checken --> dt = 5
Tend2 = 35000/30000 --> Checken --> dt=10
Tend3 = undefined --> dt = 40
=#

Param = Parameter(dir = "test4",
    dt1 = 5,
    Tend1 = 44700,
    dt2 = 447.0,
    Tend2 = 44700,
    order = 0,
    x_points = 100,
    y_points = 20
    )
Param2 = Parameter(dt1 = 20,
    dt2 = 30,
    x_points = 25,
    y_points = 5,
    Tend1 = 20,
    Tend2 = 10
)

#Only first run
Initial_solutions, Ω, Y, X, dΩ, dΓ, P = setup(Param)
# SWES_Diff_dt(Initial_solutions, Ω, Y, X, dΩ, dΓ,Param2)


SWES_Diff_dt(Initial_solutions, Ω, Y, X, dΩ, dΓ,Param, P)

