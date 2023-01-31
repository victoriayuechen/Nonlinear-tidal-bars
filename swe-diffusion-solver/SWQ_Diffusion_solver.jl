using Pkg
Pkg.activate(".")

#Add packages
using Gridap
using WriteVTK
using LineSearches: BackTracking
using Gridap.TensorValues: meas
using Parameters
using CSV, Tables
using DataFrames
using DelimitedFiles

includet("src/SWQ_Solver.jl")
includet("src/SWQ_Write_Output.jl")
includet("src/SWQ_Setup.jl")
includet("src/SWQ_Parameters.jl")

function SWES_Diff_dt(start_solutions, Ω, Y, X, dΩ, dΓ, Param, P, h, nΓ)
    @unpack dt, Tstart, Tend, theta, dir, name, show_iterations, save_CSV = Param
    x1 = solver(start_solutions, Y, X, dt, Tstart, Tend, dΩ, dΓ, show_iterations, Param, h, nΓ)
    @time writing_output(dir, name, x1, Ω, Tend, P, Param, false, h, save_CSV)
end


function Spinup(Initial_solutions, Ω, Y, X, dΩ, dΓ, Param, P, h, nΓ)
    @unpack dt_spinup, Tstart, Tend_spinup, theta = Param
    x1 = solver(Initial_solutions, Y, X, dt_spinup, Tstart, Tend_spinup, dΩ, dΓ, false, Param, h, nΓ)
    @time writing_output("SPINUP", "SPINUP", x1, Ω, Tend_spinup, P, Param, true, h, false)
end

Param = Parameter(dt = 5,
    tolerance = 1e-8,
    x_points = 50,
    y_points = 10,
    show_iterations = false,
    save_CSV = false
    )
Param_setup = Parameter(dt = 5,
    Tend = 10,
    Tend_spinup = 10
    )
function run(Param)
    Initial_solutions, Ω, Y, X, dΩ, dΓ,P, h, nΓ = setup(Param)
    Spinup(Initial_solutions, Ω, Y, X, dΩ, dΓ, Param, P, h, nΓ)
    SWES_Diff_dt(spinup_solution, Ω, Y, X, dΩ, dΓ, Param, P, h, nΓ)
    println("Done :)")
end
run(Param_setup)
run(Param)


