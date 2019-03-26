// ---- Mesh Size options ----//
lc_gen = 50.000000;
// ---Point Coordinates--------//
Point(1) = {0.000000, 0.000000, 0.000000, lc_gen};
Point(2) = {0.000000, 5000.000000, 0.000000, lc_gen};
Point(3) = {5000.000000, 5000.000000, 0.000000, lc_gen};
Point(4) = {5000.000000, 0.000000, 0.000000, lc_gen};
// ---Lines indices--------//
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
// ---Polygons--------//
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
// ---Attractors and Thresholds--------//
// ---- Other Mesh options ----//
Mesh.ElementOrder = 1;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthMax = 50.000000;
Mesh.SecondOrderIncomplete = 0;
// ------------------------//
