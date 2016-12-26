// Gmsh project created on Sun Dec 18 15:09:32 2016
lc = DefineNumber[ 0.0008, Name "Parameters/lc" ];
thk_adv = DefineNumber[ 0.00096, Name "Parameters/thk_adv" ];
thk_med = DefineNumber[ 0.00132, Name "Parameters/thk_med" ];
r = DefineNumber[ 0.01, Name "Parameters/r" ];
depth = DefineNumber[ 0.005, Name "Parameters/depth" ];
r_pla = DefineNumber[ 0.006, Name "Parameters/r_pla" ];
h_ca = (r_pla-0.002)/2;
w_ca = h_ca*2.5;

Point(1) = {0, 0, 0, lc};
Point(2) = {r, 0, 0, lc};
Point(3) = {0, r, 0, lc};
Point(4) = {-r, 0, 0, lc};
Point(5) = {0, -r, 0, lc};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Point(6) = {(r+thk_med), 0, 0, lc};
Point(7) = {0, (r+thk_med), 0, lc};
Point(8) = {-(r+thk_med), 0, 0, lc};
Point(9) = {0, -(r+thk_med), 0, lc};
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};

Point(10) = {(r+thk_med+thk_adv), 0, 0, lc};
Point(11) = {0, (r+thk_med+thk_adv), 0, lc};
Point(12) = {-(r+thk_med+thk_adv), 0, 0, lc};
Point(13) = {0, -(r+thk_med+thk_adv), 0, lc};
Circle(9) = {10, 1, 11};
Circle(10) = {11, 1, 12};
Circle(11) = {12, 1, 13};
Circle(12) = {13, 1, 10};

Point(14) = {0, (-r+r_pla), 0, lc};
Ellipse(13) = {4, 1, 2, 14};
Ellipse(14) = {2, 1, 4, 14};

Point(15) = {0, (-r+r_pla/2), 0, lc};
Point(16) = {w_ca,  (-r+r_pla/2), 0, lc};
Point(17) = {0, (-r+r_pla/2) +h_ca, 0, lc};
Point(18) = {-w_ca,  (-r+r_pla/2), 0, lc};
Point(19) = {0, (-r+r_pla/2) -h_ca, 0, lc};
Ellipse(15) = {16, 15, 18, 17};
Ellipse(16) = {17, 15, 19, 18};
Ellipse(17) = {18, 15, 16, 19};
Ellipse(18) = {19, 15, 17, 16};

Line Loop(19) = {15, 16, 17, 18};
Plane Surface(20) = {19};
Line Loop(21) = {3, 4, 14, -13};
Plane Surface(22) = {19, 21};

Line Loop(23) = {5, 6, 7, 8};
Line Loop(24) = {1, 2, 3, 4};
Plane Surface(25) = {23, 24};
Line Loop(26) = {12, 9, 10, 11};
Plane Surface(27) = {23, 26};



Extrude {0, 0, depth} {
  Surface{20};
}
//Physical Volume("cal", 4) = {1};
Extrude {0, 0, depth} {
  Surface{22};
}
//Physical Volume("plaque", 1) = {2};

Extrude {0, 0, depth} {
  Surface{25};
}
Extrude {0, 0, depth} {
  Surface{27};
}
//Physical Volume("med", 2) = {3};


//Physical Volume("adv", 3) = {4};



//Physical Surface("outer", 3) = {150, 146, 142, 154};
//Physical Surface("inner", 1) = {62, 58, 112, 108};
//Physical Surface("cross", 2) = {71, 49, 22, 20, 25, 27, 113, 155};

