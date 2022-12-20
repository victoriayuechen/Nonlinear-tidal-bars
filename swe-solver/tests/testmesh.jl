using Revise
include("../mesh_generator.jl")
using .MyMeshGenerator


generate_rectangle_mesh(10.0, 10.0, "10x10periodic2.msh", "rectangle", 0.01,true)

# Small test for whether the mesh generator works 