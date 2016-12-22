// Gmsh project created on Sun Dec 18 15:09:32 2016
lc = DefineNumber[ 0.0004, Name "Parameters/lc" ];
thk_adv = DefineNumber[ 0.00096, Name "Parameters/thk_adv" ];
thk_med = DefineNumber[ 0.00132, Name "Parameters/thk_med" ];
r = DefineNumber[ 0.01, Name "Parameters/r" ];
depth = DefineNumber[ 0.005, Name "Parameters/depth" ];
r_pla = DefineNumber[ 0.003, Name "Parameters/r_pla" ];

Point(1) = {0, 0, 0, lc};
Point(2) = {r, 0, 0, lc};
Point(3) = {0, r, 0, lc};
Point(4) = {-r, 0, 0, lc};
Point(5) = {0, -r, 0, lc};
Point(6) = {(r+thk_med), 0, 0, lc};
Point(7) = {0, (r+thk_med), 0, lc};
Point(8) = {-(r+thk_med), 0, 0, lc};
Point(9) = {0, -(r+thk_med), 0, lc};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};
Line Loop(9) = {6, 7, 8, 5};
Line Loop(10) = {2, 3, 4, 1};
Plane Surface(11) = {9, 10};
Extrude {0, 0, depth} {
  Surface{11};
}
Point(43) = {(r+thk_med+thk_adv), 0, 0, lc};
Point(44) = {0, (r+thk_med+thk_adv), 0, lc};
Point(45) = {-(r+thk_med+thk_adv), 0, 0, lc};
Point(46) = {0, -(r+thk_med+thk_adv), 0, lc};
Circle(54) = {43, 1, 44};
Circle(55) = {44, 1, 45};
Circle(56) = {45, 1, 46};
Circle(57) = {46, 1, 43};
Line Loop(58) = {57, 54, 55, 56};
Plane Surface(59) = {9, 58};
Extrude {0, 0, depth} {
  Surface{59};
}
Point(80) = {0, (-r+r_pla), 0, lc};
Ellipse(102) = {4, 1, 2, 80};
Ellipse(103) = {2, 1, 12, 80};
Line Loop(104) = {103, -102, 3, 4};
Plane Surface(105) = {104};
Extrude {0, 0, depth} {
  Surface{105};
}
