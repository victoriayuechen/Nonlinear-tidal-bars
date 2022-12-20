using Revise
include("../mesh_generator.jl")
using .MyMeshGenerator


generate_rectangle_mesh(2.0, 2.0, "10x10periodic2.msh", "rectangle", 0.1, true)

# Small test for whether the mesh generator works 