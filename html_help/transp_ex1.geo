// ---- Mesh Size options ----//
lc_gen = 150.000000;
b1_DistMin = 5.000000;
b1_DistMax = 100.000000;
b1_LcMin = 5.000000;
b1_LcMax = 50.000000;
// ---Point Coordinates--------//
Point(1) = {0.000000, 0.000000, 0.000000, lc_gen};
Point(2) = {2000.000000, 0.000000, 0.000000, lc_gen};
Point(3) = {2000.000000, 1000.000000, 0.000000, lc_gen};
Point(4) = {0.000000, 1000.000000, 0.000000, lc_gen};
Point(5) = {339.000000, 777.000000, 0.000000, b1_LcMin};
Point(6) = {1896.000000, 660.000000, 0.000000, b1_LcMin};
Point(7) = {1615.000000, 459.000000, 0.000000, b1_LcMin};
Point(8) = {1247.000000, 570.000000, 0.000000, b1_LcMin};
Point(9) = {1136.000000, 308.000000, 0.000000, b1_LcMin};
Point(10) = {684.000000, 232.000000, 0.000000, b1_LcMin};
Point(11) = {735.000000, 494.000000, 0.000000, b1_LcMin};
Point(12) = {1044.000000, 806.000000, 0.000000, b1_LcMin};
Point(13) = {1325.000000, 488.000000, 0.000000, b1_LcMin};
Point(14) = {721.000000, 570.000000, 0.000000, b1_LcMin};
// ---Lines indices--------//
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
// ---Polygons--------//
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Point{1} In Surface{1};
Point{2} In Surface{1};
Point{3} In Surface{1};
Point{4} In Surface{1};
Point{5} In Surface{1};
Point{6} In Surface{1};
Point{7} In Surface{1};
Point{8} In Surface{1};
Point{9} In Surface{1};
Point{10} In Surface{1};
Point{11} In Surface{1};
Point{12} In Surface{1};
Point{13} In Surface{1};
Point{14} In Surface{1};
// ---Attractors and Thresholds--------//
Field[1] = Attractor;
Field[1].NodesList = {5,6,7,8,9,10,11,12,13,14};
Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = b1_LcMin;
Field[2].LcMax = b1_LcMax;
Field[2].DistMin = b1_DistMin;
Field[2].DistMax = b1_DistMax;
Field[3] = Min;
Field[3].FieldsList = {2};
Background Field = 3;
// ---- Other Mesh options ----//
Mesh.ElementOrder = 1;
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;
Mesh.CharacteristicLengthMax = 150.000000;
Mesh.SecondOrderIncomplete = 0;
// ------------------------//
