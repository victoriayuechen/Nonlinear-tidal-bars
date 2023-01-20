
using Gridap
using GridapGmsh
using Revise

includet("../linear_SWE.jl")
using .MyLinearSWE
function ζ₀((x,y))
    h =0.01*exp(-0.1*(x-30)^2 -0.1*(y-50)^2)
    h
end

function u₀((x,y))
    u = VectorValue(0.0,0.0)
    u
end

function forcefunc((x,y),t)
    f = VectorValue(0.0,0.0*cos(π*(1/5)*t))
    f
end


order = 1
degree = 4
model = GmshDiscreteModel("swe-solver/meshes/ref.msh")
DC = ["bottom","top"]
filename = "test1"
dir = "output_swe/linear_SWE/test1"
H = 0.5
latitude = 50
Tend = 50
dt = 1
tcapture = 1.0

time_1 = time()
run_linear_SWE(order,degree,ζ₀,u₀,forcefunc,Tend,dt,model,H,DC,dir,latitude,filename,tcapture::Float64)
println(time() - time_1)