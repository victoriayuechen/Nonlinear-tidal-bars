using Pkg
Pkg.activate(".")

using Gridap

#Parameters

B = 1000    #[m]
H = 3       #[m]
L = 10000   #[m]
dy = 50     #[m]
dx = 100    #[m]
dt = 5      #[s]
g = 9.81    #[m/s^2]
cd = 0.0025 #[-]


#Create model
domain = (B,L)
partition = (B/dy,L/dx)
model = CartesianDiscreteModel(domain,partition;isperiodic=(true,false))

labels = get_face_labeling(model)
add_tag_from_tags!(labels,"bottom",[1,2,5])
add_tag_from_tags!(labels,"left",[7])
add_tag_from_tags!(labels,"right",[8])
add_tag_from_tags!(labels,"top",[3,4,6])
add_tag_from_tags!(labels,"inside",[9])

order =  1
refFE = ReferenceFE(lagrangian,VectorValue(3,Float64),order)
V = TestFESpace(model,refFE,dirichlet_tags =["left","right"]) #We only want to set lateral velocity to zero at the sides, need to find 1) way to get only side and 2) Way to set only lateral to zero and the rest not


#Not correct boundary conditions
gl(x,t::Real) = VectorValue(0.0,0.0,0.0)
gr(x,t::Real) = VectorValue(0.0,0.0,0.0)
gl(t::Real) = x -> gl(x,t)
gr(t::Real) = x -> gr(x,t)


degree = 2
Ω = Interior(Model)
dΩ = Measure(Ω,degree)

