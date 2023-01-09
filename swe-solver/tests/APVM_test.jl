using Revise
includet("../nonlinear_APVM_transient.jl")
using Gridap



#Variable functions to be used to setup model, used for local tests
function h₀((x,y))
    hout = -topography((x,y)) +  0.5 + 0.05*exp(-0.01*(x-50)^2 -0.01*(y-25)^2)
    hout
end

function topography((x,y))
    p = 0.0
    p
end

function u₀((x,y))
    u = VectorValue(0.0,0.0)
    u
end

function forcfunc((x,y),t)
    f = VectorValue(0.0,0.0*0.5*cos(π*(1/50)*t))
    f
end

outputdir = "output_swe"
dir = "NL_SWE_APVM_test"
if !isdir(outputdir)
    mkdir(outputdir)
end

if !isdir(joinpath(outputdir,dir))
    mkdir(joinpath(outputdir,dir))
end

B = 100
L = 100
partition = (50,50)
domain = (0,B,0,L)

#=
Input:
order       = order of FE polynomials
degree      = lebesgue measure with quadruture rule of degree
h_0         = initial fluid depth h
u_0         = initial velocity vector field
topography  = bottom topography, passed as a function of x and Y
forcfunc    = The forcing function, passed as a function in x,y
outputdir   = the output directory of all output folders
dir         = the actual output folder NL_SWE_APVM_test
Periodic    = if true periodic in y-dir
Tend        = Total runtime
dt          = Time setup
=#
APVM_run(0,4,h₀,u₀,topography,forcfunc,joinpath(outputdir,dir),true,100.0,1.0,domain,partition)
