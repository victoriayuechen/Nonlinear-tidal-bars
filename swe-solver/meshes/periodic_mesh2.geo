// Gmsh project created on Thu Dec 08 15:19:19 2022
//+
SetFactory("OpenCASCADE");
//+
Rectangle(1) = {0, 0, 0, 1, 1, 0};
Periodic Curve {3} = {1} Translate {0, 1, 0};//+
Physical Curve("top") = {3};
//+
Physical Curve("bottom") = {1};
//+
Physical Curve("left") = {4};
//+
Physical Curve("right") = {2};
Physical Surface("surface") = {1};
Physical Point("top_points") = {3, 4};
Physical Point("bottom_points") = {1, 2};