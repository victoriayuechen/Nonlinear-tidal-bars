using Revise
include("../mesh_generator.jl")
using .MyMeshGenerator


generate_rectangle_mesh(100.0, 100.0, "100x100periodic2.msh", "rectangle", 2.0,true)

# Small test for whether the mesh generator works 