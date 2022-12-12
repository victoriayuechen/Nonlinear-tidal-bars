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
function generate_rectangle_mesh(Lx::Float64, Ly::Float64, filename::String, modelname::String, lc::Float64)
    # Initialise mesh generator 
    gmsh.initialize()
    gmsh.model.add(modelname)

    # Add the points 
    gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
    gmsh.model.geo.addPoint(Lx, 0, 0, lc, 2)
    gmsh.model.geo.addPoint(Lx, Ly, 0, lc, 3)
    gmsh.model.geo.addPoint(0, Ly, 0, lc, 4)

    # Add the lines 
    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(2, 3, 2)
    gmsh.model.geo.addLine(3, 4, 3)
    gmsh.model.geo.addLine(4, 1, 4)

    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.setPeriodic(1, [3], [1], [1, 0, 0, 0, 0, 1, Ly, 0, 0, 0, 1, 0, 0, 0, 0, 1])

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
