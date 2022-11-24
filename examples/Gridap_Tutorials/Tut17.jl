using Gridap


#This will make the model and domain
T = CartesianDiscreteModel((0,1,0,1),(20,20))
Ω = Interior(T)
dΩ = Measure(Ω,2)#this is the lebesgue measure of the domain

refFE = ReferenceFE(lagrangian,Float64,1);
V = TestFESpace(T,refFE,dirichlet_tags= 'boundary')