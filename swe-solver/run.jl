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
α = 7.29e-5 #[1/s]


#Create model
domain = (B,L)
partition = (B/dy,L/dx)
model = CartesianDiscreteModel(domain,partition;isperiodic=(true,false))#TODO: Check if this the correct side has been given periodic BC

order = 1
refFE_u = ReferenceFE(lagrangian,VectorValue(2,Float64),order)
V = TestFESpace(model,refFE,dirichlet_tags)
