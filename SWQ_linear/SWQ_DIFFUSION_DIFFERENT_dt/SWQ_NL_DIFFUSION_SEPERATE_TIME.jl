using Pkg
Pkg.activate(".")

#Add packages
using Gridap
using WriteVTK
using LineSearches: BackTracking
using Gridap.TensorValues: meas

includet("src/SWQ_Solver.jl")
includet("src/SWQ_Write_Output.jl")
includet("src/SWQ_Setup.jl")

global η = 7.29e-5                             #angular speed of Earth rotation        (s^(-1))
latitude = 52
global f = 2*η*sin(latitude*(π/180))           #coriolis parameter                     (s^(-1))
global g = 9.81                                #Gravitational constant                 (ms^(-2))
global H = 3                                   #Constant layer depth at rest           (m)
global U_start = 0.5                           #Background current amplitude           (ms^(-1))
global σ = 2*pi/44700                          #Tidal frequency                        (s^(-1))
global cD = 0.0025                             #Drag coefficient                       ()

##''''''''''''''Stabilization Parameters''''''''''''''##
global α = 1e-6                                #Based on ν
global ν = 1e-6                                #From Anna Louka

function SWES_Diff_dt(Initial_solutions, Ω, Y, X, dΩ, dΓ, dt1, dt2, Tend1, Tend2, dir)
    x1 = solver(Initial_solutions, Y, X, dt1, 0.0, Tend1, dΩ, dΓ, false)
    writing_output(dir, x1, Ω, Tend1)
    x2 = solver(xn,Y,X,dt2,Tend1,Tend2,dΩ,dΓ, true)
    writing_output(dir, x2, Ω, Tend2)
end

Initial_solutions, Ω, Y, X, dΩ, dΓ = setup()

SWES_Diff_dt(Initial_solutions, Ω, Y, X, dΩ, dΓ, 5, 10, 400, 800, "./test")
