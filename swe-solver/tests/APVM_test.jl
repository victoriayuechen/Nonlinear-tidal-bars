using Revise
includet("../nonlinear_APVM_transient.jl")
using Gridap
using GridapGmsh
using .MyNonlinearAPVM


#Variable functions to be used to setup model, used for local tests
function D₀((x,y))
    Dout = -topography((x,y)) +  3 #+ 0.01*exp(-0.1*(x-50)^2 -0.1*(y-30)^2)
    Dout
end

function topography((x,y))
    p =  0.1 * cos(π/B * y) * cos(2π/L * x)
    p
end

function u₀((x,y))
    u = VectorValue(0.0,0.0)
    u
end

function forcfunc((x,y),t)
    U_start =  0.5 
    η = 7.29e-5
    latitude = 50
    cD = 0.0025
    f = 2*η*sin(latitude*(π/180))
    σ = 2*pi/22350
    H = 3
    f = VectorValue(σ*U_start*cos(σ*t)+cD/H*abs(U_start*sin(σ*t))*U_start*sin(σ*t),f*U_start*sin(σ*t))
    f
end

outputdir = "output_swe"
dir = "NLAPVM_hepkema_400_order1_diffFESPACE"
if !isdir(outputdir)
    mkdir(outputdir)
end

if !isdir(joinpath(outputdir,dir))
    mkdir(joinpath(outputdir,dir))
end

B = 1000
L = 10000
partition = (100,20)
domain = (0,L,0,B)
# model = GmshDiscreteModel("swe-solver/meshes/1000x10000periodic.msh")
model = CartesianDiscreteModel(domain,partition;isperiodic=(true,false)) 
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"bottom",[1,2,5])
add_tag_from_tags!(labels,"left",[7])
add_tag_from_tags!(labels,"right",[8])
add_tag_from_tags!(labels,"top",[3,4,6])
add_tag_from_tags!(labels,"inside",[9])
DC = ["bottom","top"]
Tend = 44700
dt = 100
#=
Input:
order       = order of FE polynomials
degree      = lebesgue measure with quadruture rule of degree
D_0         = initial fluid depth h
u_0         = initial velocity vector field
topography  = bottom topography, passed as a function of x and Y
forcfunc    = The forcing function, passed as a function in x,y
outputdir   = the output directory of all output folders
dir         = the actual output folder NL_SWE_APVM_test
Periodic    = if true periodic in y-dir
Tend        = Total runtime
dt          = Time setup
=#

@time APVM_run(1,4,D₀,u₀,topography,forcfunc,joinpath(outputdir,dir),true,Tend,dt,model,DC)


