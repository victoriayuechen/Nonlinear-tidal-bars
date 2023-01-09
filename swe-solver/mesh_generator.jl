module MyMeshGenerator
using Gmsh: gmsh

export generate_rectangle_mesh


"""
Generate rectangular mesh based on coordinates and triangle size. 
Parameters:
    Ly: Length of y-direction
    Lx: Length of in x-direction
    filename: Name of the mesh file 
    modelname: Name of the model 
    lc: The size of the triangles, default is 1e-2. 
Creates a <filename>.msh in the meshes folder.
"""
function generate_rectangle_mesh(Lx::Float32, Ly::Float32, filename::String, modelname::String, lc::Float32, periodic::Bool)
    # Initialise mesh generator 
    gmsh.initialize()
    gmsh.model.add(modelname)

    # Add the points 
    p1 = gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
    p2 = gmsh.model.geo.addPoint(Lx, 0, 0, lc, 2)
    p3 = gmsh.model.geo.addPoint(Lx, Ly, 0, lc, 3)
    p4 = gmsh.model.geo.addPoint(0, Ly, 0, lc, 4)

    # Add the lines 
    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(3, 2, 2)
    gmsh.model.geo.addLine(3, p4, 3)
    gmsh.model.geo.addLine(4, 1, p4)

    gmsh.model.geo.addCurveLoop([4, 1, -2, 3], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)

    gmsh.model.geo.synchronize()


    if periodic
        transformation_matrix = zeros(4, 4)
        transformation_matrix[1, 1] = 1
        transformation_matrix[2, 2] = 1 
        transformation_matrix[2,4] = -Ly
        transformation_matrix[3, 3] = 1
        transformation_matrix[4, 4] = 1
        transformation_matrix = vec(transformation_matrix')

        gmsh.model.mesh.set_periodic(1, [1], [3], transformation_matrix)
    end



    gmsh.model.addPhysicalGroup(0, [1, 2], 1, "bottom")
    gmsh.model.addPhysicalGroup(0, [3, 4], 2, "top")
    gmsh.model.addPhysicalGroup(1, [1], 1, "bottom")
    gmsh.model.addPhysicalGroup(1, [3], 2, "top")
    gmsh.model.addPhysicalGroup(1, [2], 3, "right")
    gmsh.model.addPhysicalGroup(1, [4], 4, "left")

    gmsh.model.addPhysicalGroup(2, [1], 1, "interior")
    gmsh.model.mesh.generate(2)
    gmsh.write(filename)

    gmsh.finalize()
end 
end



