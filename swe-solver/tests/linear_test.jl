
using Gridap
using GridapGmsh
using Revise

includet("../linear_SWE.jl")
using .MyLinearSWE
function ζ₀((x,y))
    h =0.01*exp(-0.1*(x-50)^2 -0.1*(y-30)^2)
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
model = GmshDiscreteModel("swe-solver/meshes/100x100periodic_testing.msh")
DC = ["right","left"]
filename = "test_01"
dir = joinpath("output_swe/linear_SWE",filename)
H = 0.5
latitude = 52
Tend = 100
dt = 0.1
tcapture = 2.0

time_1 = time()
run_linear_SWE(order,degree,ζ₀,u₀,forcefunc,Tend,dt,model,H,DC,dir,latitude,filename,tcapture::Float64)
println(time() - time_1)