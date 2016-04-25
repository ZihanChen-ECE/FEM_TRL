// Gmsh project created on Thu Mar 31 21:19:22 2016
Point(1) = {0,0,0,0.01};
Point(2) = {0,1,0,0.01};
Point(3) = {1,1,0,0.01};
Point(4) = {1,0,0,0.01};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};

// Specify physical groups
Plane Surface(6) = {5};
Physical Line(666) = {4};
Physical Line(333) = {3,1};
Physical Line(444) = {2};
Physical Surface(10) = {6};

