using Revise
include("../mesh_generator.jl")
using .MyMeshGenerator


generate_rectangle_mesh(5.0, 5.0, "testtest", "rectangle", 0.01)

# Small test for whether the mesh generator works 