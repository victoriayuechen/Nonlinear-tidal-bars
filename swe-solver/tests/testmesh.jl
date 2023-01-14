using Revise
includet("../mesh_generator.jl")
using .MyMeshGenerator


generate_rectangle_mesh(Float32(100), Float32(100), "1000x10000periodic.msh", "rectangle", Float32(100),true)

#generate_arc_mesh(Float32(100), Float32(100), "arc-larger.msh", "arc", Float32(1) ,true)

#generate_outlet_mesh(Float32(100), Float32(100), "canal.msh", "canal", Float32(1) ,true)

# Small test for whether the mesh generator works 