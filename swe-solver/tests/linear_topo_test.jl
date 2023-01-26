
using Gridap
includet("../linear_topo_SWE.jl")
using .MylinearSWETopo



function ζ₀((x,y))
    ```
    Function of initial ζ
    ``` 	
    h = 0.01*exp(-0.1*(x-50)^2 -0.1*(y-25)^2)
    h
end


function u₀((x,y))
    ```
    Function of initial u
    ```
    u = VectorValue(0.0,0.0)
    u
end

function forcefunc((x,y),t)
    ```
    Forcing function
    ```
    func = VectorValue(0.0,0.0)
    func
end

function topography((x,y))
    ```
    Function for topography
    ```
    b = 0.3*exp(-0.01*(y-50)^2)
    b
end


#Parameters
#################
order = 1
degree = 4
B = 100 #Height of rectangle
L = 100 #Width of rectangle
partition = (100,100)
model = CartesianDiscreteModel((0,B,0,L),partition;isperiodic=(true,false))
#Make labels
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"bottom",[1,2,5])
add_tag_from_tags!(labels,"left",[7])
add_tag_from_tags!(labels,"right",[8])
add_tag_from_tags!(labels,"top",[3,4,6])
add_tag_from_tags!(labels,"inside",[9])
DC = ["bottom","top"]
linear_dir = "output_swe/linear_SWE_topo"
if isdir(linear_dir)
    dir = joinpath(linear_dir,"test")
else
    mkdir(linear_dir)
    dir = joinpath(linear_dir,"test")
end
filename = "testt"
H = 0.5
latitude = 50
Tend = 20
dt = 0.1
#################

linear_SWE_topo(order,degree,ζ₀,u₀,topography,forcefunc,Tend,dt,model,H,DC,dir,latitude,filename)