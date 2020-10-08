/*
Author: T. Ulrich
Generates a surface mesh of a 3D box, to be used for generating a structural model with SimModeler
 */


mesh_size = 10e3;
Xmax = -160e3;
Xmin = 215e3;
Ymin = 1235e3;
Ymax = 1605e3;
Zmin=-200e3;
Zmax=5e3;

//Create the Volume
Point(1) = {Xmin, Ymin, Zmax, mesh_size};
Point(2) = {Xmin, Ymax, Zmax, mesh_size};
Point(3) = {Xmax, Ymax, Zmax, mesh_size};
Point(4) = {Xmax, Ymin, Zmax, mesh_size};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1,2,3,4};
Plane Surface(1) = {5};

Extrude {0,0, Zmin} { Surface{1}; }
