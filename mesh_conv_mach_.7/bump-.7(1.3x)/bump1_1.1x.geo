// Gmsh project created on Mon Oct 25 09:30:34 2021

//+

Point(1) = {0, 0, -0, 1.0};
//+
Point(2) = {3, 0, 0, 1.0};
//+
Point(3) = {3, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {1, 0, 0, 1.0};
//+
Point(6) = {2, 0, 0, 1.0};
//+
Point(7) = {1.5, .1, 0, 1.0};
//+
Point(8) = {1.5, 1, 0, 1.0};
//+
Point(9) = {1.5, -1.2, 0, 1.0};


//+
Line(3) = {3, 2};
//+
Line(4) = {2, 6};
//+
Line(5) = {5, 1};
//+
Line(6) = {1, 4};
//+
Circle(7) = {7, 9, 5};
//+
Circle(8) = {7, 9, 6};
//+
Point(10) = {1, 1, 0, 1.0};
//+
Point(11) = {2, 1, 0, 1.0};
//+
Line(9) = {4, 10};
//+
Line(10) = {10, 8};
//+
Line(11) = {8, 11};
//+
Line(12) = {11, 3};
//+
Line(13) = {5, 10};
//+
Line(14) = {7, 8};
//+
Line(15) = {6, 11};
//+
Curve Loop(1) = {6, 9, -13, 5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {13, 10, -14, 7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, 15, -11, -14};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {4, 15, 12, 3};
//+
Plane Surface(4) = {4};
//+
Physical Curve("inlet", 16) = {6};
//+
Physical Curve("outlet", 17) = {3};
//+
Physical Curve("lower_wall", 18) = {8, 7, 5, 4};
//+
Physical Curve("upper_wall", 19) = {9, 10, 11, 12};
//+
Physical Surface("flow", 20) = {1, 2, 3, 4};
//+
Transfinite Surface {1} = {1, 4, 10, 5};
//+
Transfinite Curve {3} = 141 Using Progression .9804;
//+
Transfinite Curve {6, 13,14,15} = 141 Using Progression 1.02;
//+
Transfinite Curve {5} = 84 Using Progression 1.02;
//+
Transfinite Curve {9} = 84 Using Progression .9804;
//+
Transfinite Surface {2} = {5, 10, 8, 7};
//+
Transfinite Surface {3} = {7, 8, 11, 6};
//+
Transfinite Surface {4} = {6, 11, 3, 2};
//+
Transfinite Curve {10} = 60 Using Progression 1.02;
//+
Transfinite Curve {11} = 60 Using Progression .9804;
//+
Transfinite Curve {7} = 60 Using Progression .9804;
//+
Transfinite Curve {8} = 60 Using Progression .9804;
//+
Transfinite Curve {12} = 84 Using Progression 1.02;
//+
Transfinite Curve {4} = 84 Using Progression .9804;

//+
Recombine Surface {1, 2, 3, 4};
//+
Recombine Surface {1};
//+
Recombine Surface {2};
//+
Recombine Surface {4};
//+
Recombine Surface {3};
