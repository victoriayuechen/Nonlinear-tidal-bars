
include("../nonlinear_APVM_transient.jl")
using .APVM_solver



# #Variable functions to be used to setup model, used for local tests
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

outputdir = "output_swe"
dir = "NL_SWE_APVM_test"
if !isdir(outputdir)
    mkdir(outputdir)
end

if !isdir(joinpath(outputdir,dir))
    mkdir(joinpath(outputdir,dir))
end

APVM_run(0,4,h₀,u₀,topography,joinpath(outputdir,dir),true,100.0,1.0)