using Revise
includet("../swe-solver/mesh_generator.jl")
using .MyMeshGenerator


generate_rectangle_mesh(Float32(1), Float32(1), "1x1periodic01.msh", "rectangle", Float32(0.1),true)

#generate_arc_mesh(Float32(100), Float32(100), "arc-larger.msh", "arc", Float32(1) ,true)

#generate_outlet_mesh(Float32(100), Float32(100), "canal.msh", "canal", Float32(1) ,true)

# Small test for whether the mesh generator works 