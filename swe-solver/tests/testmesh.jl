using Revise
includet("../mesh_generator.jl")
using .MyMeshGenerator


generate_rectangle_mesh(Float32(100), Float32(100), "100x100periodic2.msh", "rectangle", Float32(4),true)

# Small test for whether the mesh generator works 