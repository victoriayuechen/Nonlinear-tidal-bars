using Gmsh: gmsh

#=
Tags:

top: y-max edge + two y-max corner points
bottom: y-min edge + two y-min corner points
left: x-min edge
right: x-max edge

=#

gmsh.initialize()

gmsh.model.add("p_rectangle")

size_x = 2
size_y = 10

lc = 1e-2

gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(size_x, 0, 0, lc, 2)
gmsh.model.geo.addPoint(size_x, size_y, 0, lc, 3)
gmsh.model.geo.addPoint(0, size_y, 0, lc, 4)

gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)

gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
gmsh.model.geo.addPlaneSurface([1], 1)

gmsh.model.geo.synchronize()

gmsh.model.mesh.setPeriodic(1, [3], [1], [1, 0, 0, 0, 0, 1, size_y, 0, 0, 0, 1, 0, 0, 0, 0, 1])

gmsh.model.addPhysicalGroup(0, [1, 2], 1, "bottom")
gmsh.model.addPhysicalGroup(0, [3, 4], 2, "top")
gmsh.model.addPhysicalGroup(1, [1], 1, "bottom")
gmsh.model.addPhysicalGroup(1, [3], 2, "top")
gmsh.model.addPhysicalGroup(1, [2], 3, "right")
gmsh.model.addPhysicalGroup(1, [4], 4, "left")

gmsh.model.addPhysicalGroup(2, [1], 1, "interior")

gmsh.model.mesh.generate(2)

gmsh.write("swe-solver/meshes/p_rectangle.msh")

if !("-nopopup" in ARGS)
    gmsh.fltk.run()
end

gmsh.finalize()