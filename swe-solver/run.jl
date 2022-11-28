using Pkg
Pkg.activate(".")

using Gridap

function create_model(Nx::Int64, Ny::Int64, Xmin::Float64, Xmax::Float64, Ymin::Float64, Ymax::Float64)
    domain = (Xmin, Xmax, Ymin, Ymax)
    partition = (Nx,Ny)
    
    CartesianDiscreteModel(domain,partition)
end 
