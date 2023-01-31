module MyMeshGenerator
using Gmsh: gmsh

export generate_rectangle_mesh, generate_arc_mesh, generate_outlet_mesh

"""
Generate rectangular mesh based on coordinates and triangle size. 
Parameters:
    Ly: Length of y-direction
    Lx: Length of in x-direction
    filename: Name of the mesh file 
    modelname: Name of the model 
    lc: The size of the triangles, default is 1e-2. 
    periodic: Whether the top and bottom of rectangle is periodic 
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
        transformation_matrix[1 ,4] = Lx
        transformation_matrix[3, 3] = 1
        transformation_matrix[4, 4] = 1
        transformation_matrix = vec(transformation_matrix')

        gmsh.model.mesh.set_periodic(1, [2], [4], transformation_matrix)
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

"""
Generate arched mesh based on coordinates and triangle size. 
Parameters:
    inner radius: Inner Radius of the arc (similar to curvature)
    width: Width of in river
    filename: Name of the mesh file 
    modelname: Name of the model 
    lc: The size of the triangles, default is 1e-2. 
    periodic: Whether the top and bottom of arc is periodic 
Creates a <filename>.msh in the meshes folder.
"""
function generate_arc_mesh(innerradius::Float32, width::Float32, filename::String, modelname::String, lc::Float32, periodic::Bool)
    # Initialise mesh generator 
    gmsh.initialize()
    gmsh.model.add(modelname)
    outside = innerradius + width 

    center = gmsh.model.geo.add_point(0, 0, 0, lc)
    p1 = gmsh.model.geo.add_point(innerradius, 0, 0, lc)
    p2 = gmsh.model.geo.add_point(outside, 0.0, 0.0, lc)
    p3 = gmsh.model.geo.add_point(0.0, outside, 0.0, lc)
    p4 = gmsh.model.geo.add_point(0.0, innerradius, 0.0, lc)

    l1 = gmsh.model.geo.add_line(p1, p2)
    l2 = gmsh.model.geo.add_circle_arc(p2, center, p3)
    l3 = gmsh.model.geo.add_line(p3, p4)
    l4 = gmsh.model.geo.add_circle_arc(p4, center, p1)
 
    # Create the closed curve loop and the surface
    loop = gmsh.model.geo.add_curve_loop([l1, l2, l3, l4])
    surf = gmsh.model.geo.add_plane_surface([loop])

    gmsh.model.geo.synchronize()

    # Create the physical domains
    gmsh.model.add_physical_group(1, [l1], -1, "bottom")
    gmsh.model.add_physical_group(1, [l2], -1, "left")
    gmsh.model.add_physical_group(1, [l3], -1, "top")
    gmsh.model.add_physical_group(1, [l4], -1, "right")
    gmsh.model.add_physical_group(2, [surf])

    # Define Affine Operator 
    if periodic
        transformation_matrix = zeros(4, 4)
        transformation_matrix[1, 2] = 1  
        transformation_matrix[2, 1] = -1 
        transformation_matrix[3, 3] = 1
        transformation_matrix[4, 4] = 1
        transformation_matrix = vec(transformation_matrix')
        gmsh.model.mesh.set_periodic(1, [l1], [l3], transformation_matrix)
    end 

    # Generate a 2D mesh
    gmsh.model.mesh.generate(2)
    gmsh.write(filename)
    gmsh.finalize()
end


"""
Generate canal mesh based on length and width and triangle size. 
Parameters:
    Lx: Maximum width of canal
    Ly: Length of canal
    filename: Name of the mesh file 
    modelname: Name of the model 
    lc: The size of the triangles, default is 1e-2. 
    periodic: Whether the top and bottom of arc is periodic 
Creates a <filename>.msh in the meshes folder.
""" 
function generate_outlet_mesh(Lx::Float32, Ly::Float32, filename::String, modelname::String, lc::Float32, periodic::Bool)
    # Initialise mesh generator 
    gmsh.initialize()
    gmsh.model.add(modelname)

    # Add the points 
    p1 = gmsh.model.geo.add_point(0, 0, 0, lc, 1)
    p2 = gmsh.model.geo.add_point(Lx, 0, 0, lc, 2)
    p3 = gmsh.model.geo.add_point(Lx, Ly, 0, lc, 3)
    p4 = gmsh.model.geo.add_point(0, Ly, 0, lc, 4)

    # Between 1 and 4 
    center_l = gmsh.model.geo.add_point(-1*Lx, 0.5*Ly, 0, lc, 5)
    # Between 2 and 3
    center_r = gmsh.model.geo.add_point(2*Lx, 0.5*Ly, 0, lc, 6)

    # Add the lines 
    l1 = gmsh.model.geo.add_line(1, 2)
    l2 = gmsh.model.geo.add_circle_arc(2, center_r, 3)
    l3 = gmsh.model.geo.add_line(3, p4)
    l4 = gmsh.model.geo.add_circle_arc(4, center_l, 1)

    # Create the closed curve loop and the surface
    loop = gmsh.model.geo.add_curve_loop([l1, l2, l3, l4])
    surf = gmsh.model.geo.add_plane_surface([loop])

    gmsh.model.geo.synchronize()

    # Create the physical domains
    gmsh.model.add_physical_group(1, [l1], -1, "bottom")
    gmsh.model.add_physical_group(1, [l2], -1, "left")
    gmsh.model.add_physical_group(1, [l3], -1, "top")
    gmsh.model.add_physical_group(1, [l4], -1, "right")
    gmsh.model.add_physical_group(2, [surf])

    # Define Affine Operator 
    if periodic
        transformation_matrix = zeros(4, 4)
        transformation_matrix[1, 2] = 1  
        transformation_matrix[2, 1] = -1 
        transformation_matrix[3, 3] = 1
        transformation_matrix[4, 4] = 1
        transformation_matrix = vec(transformation_matrix')
        gmsh.model.mesh.set_periodic(1, [l1], [l3], transformation_matrix)
    end 

    # Generate a 2D mesh
    gmsh.model.mesh.generate(2)
    gmsh.write(filename)
    gmsh.finalize()

end

end 


