using Revise
includet("../swe-solver/nonlinear_APVM.jl")
using Gridap
using GridapGmsh
using .MyNonlinearAPVM

W = 1
L = 10

η = 7.29e-5
latitude = 50
cD = 0.0025
f = 2*η*sin(latitude*(π/180))
σ = 2*pi/Tend
H = 3.0

A = (π/W)
B = (π/L)
C = ((4 * π)/(1 * Tend))

D = (π/W)
E = (π/(7 * L))
F = A
G = B

top = 0.0
xref = W / 2
cd = cD
#Variable functions to be used to setup model, used for local tests
function D₀((x,y))
    Dout = -topography((x,y)) +  H + 0.01*exp(-0.01(x-0.5)^2 - 0.01*(y-5)^2)
    Dout
end

function topography((x,y))
    p =  0.0#0.1 * cos(π/B * y) * cos(2π/L * x)
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
    σ = 2*pi/Tend
    H = 3
    f = VectorValue(0.0,0.0)#VectorValue(σ*U_start*cos(σ*t)+cD/H*abs(U_start*sin(σ*t))*U_start*sin(σ*t),f*U_start*sin(σ*t))
    f
end

# function fu((x,y),t)
#     fu = VectorValue((cd*(0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)*sqrt((0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)^2+(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)^2))/(0.05*sin(F*x)*sin(G*y)*sin(C*t)-top+2*H)- 
#     0.1*E*sin(E*(x-xref))*sin(D*y)*sin(C*t)*(0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)+ 
#     0.1*D*cos(E*(x-xref))*cos(D*y)*sin(C*t)*(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)- 
#     f*(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)+ 
#     0.05*F*g*cos(F*x)*sin(G*y)*sin(C*t) + 
#     0.1*C*cos(E*(x-xref))*sin(D*y)*cos(C*t),
#     (cd*(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)*sqrt((0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)^2+(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)^2))/(0.05*sin(F*x)*sin(G*y)*sin(C*t)-top+2*H)+ 
#     0.1*A*cos(A*x)*sin(B*y)*sin(C*t)*(0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)+ 
#     f*(0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)+ 
#     0.1*B*sin(A*x)*cos(B*y)*sin(C*t)*(0.1*sin(A*x)*sin(B*y)*sin(C*z)+1.0)+ 
#     0.05*G*g*sin(F*x)*cos(G*y)*sin(C*t)+0.1*C*sin(A*x)*sin(B*y)*cos(C*t))
#     fu
# end

# function fh((x,y),t)
#     fh = -0.1*E*sin(E*(x-xref))*sin(D*y)*sin(C*t)*(0.05*sin(F*x)*sin(G*y)*sin(C*t)-top+2*H)+
#     0.1*B*sin(A*x)*cos(B*y)*sin(C*t)*(0.05*sin(F*x)*sin(G*y)*sin(C*t)-top+2*H)+
#     0.05*F*cos(F*x)*sin(G*y)*sin(C*t)*(0.1*cos(E*(x-xref))*sin(D*y)*sin(C*t)+1.0)+
#     0.05*G*sin(F*x)*cos(G*y)*sin(C*t)*(0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)+
#     0.05*C*sin(F*x)*sin(G*y)*cos(C*t)
#     fh
# end

# function mu((x,y),t)
#     u = VectorValue(0.1*sin(D*y)*cos(E*(x-xref))*sin(C*t)+1.0,
#     0.1*sin(A*x)*sin(B*y)*sin(C*t)+1.0)
#     u    
# end

# function mh((x,y),t)
#     h = 0.05*sin(F*x)*sin(G*y)*sin(C*t)+H
#     h
# end

outputdir = "output_swe"
dir = "NL_SWE_APVM_test_gc"
if !isdir(outputdir)
    mkdir(outputdir)
end

if !isdir(joinpath(outputdir,dir))
    mkdir(joinpath(outputdir,dir))
end

W = 1
L = 10
partition = (5,50)
domain = (0,L,0,W)
# model = GmshDiscreteModel("swe-solver/meshes/1000x10000periodic.msh")
model = CartesianDiscreteModel(domain,partition;isperiodic=(true,false)) 
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"bottom",[1,2,5])
add_tag_from_tags!(labels,"left",[7])
add_tag_from_tags!(labels,"right",[8])
add_tag_from_tags!(labels,"top",[3,4,6])
add_tag_from_tags!(labels,"inside",[9])
DC = ["bottom","top"]
Tend = 100
dt = 1
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
time_1 = time()
@time APVM_run(0,4,D₀,u₀,topography,forcfunc,joinpath(outputdir,dir),true,Tend,dt,model,DC)
println(time() - time_1)