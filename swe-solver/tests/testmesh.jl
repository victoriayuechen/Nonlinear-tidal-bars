using Revise
include("../mesh_generator.jl")
using .MyMeshGenerator


generate_rectangle_mesh(10.0, 10.0, "test.msh", "rectangle", 1.0, true)

# Small test for whether the mesh generator works 