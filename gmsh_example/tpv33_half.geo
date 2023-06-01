lc= 5000;
lc_fault = 100;

Point(1) = {-50000, 0, -50000, lc};
Point(2) = {-50000, 0, 0, lc};
Point(3) = {50000, 0, 0, lc};
Point(4) = {50000, 0, -50000, lc};
Point(100) = {-12000, 0, -10000, lc};
Point(101) = {4000, 0, -10000, lc};
Point(102) = {4000, 0, 0, lc};
Point(103) = {-12000, 0, 0, lc};
Point(200) = {-6000, 0, -6000, lc_fault};
Point(201) = {-6000, 0, -5450, lc_fault};
Point(202) = {-5450, 0, -6000, lc_fault};
Point(203) = {-6000, 0, -6550, lc_fault};
Point(204) = {-6550, 0, -6000, lc_fault};
Point(211) = {-6000, 0, -5200, lc_fault};
Point(212) = {-5200, 0, -6000, lc_fault};
Point(213) = {-6000, 0, -6800, lc_fault};
Point(214) = {-6800, 0, -6000, lc_fault};
Point(1001) = {-50000, 800, -50000, lc};
Point(1002) = {-50000, 800, 0, lc};
Point(1003) = {50000, 800, 0, lc};
Point(1004) = {50000, 800, -50000, lc};
Point(1011) = {-50000, 50000, -50000, lc};
Point(1012) = {-50000, 50000, 0, lc};
Point(1013) = {50000, 50000, 0, lc};
Point(1014) = {50000, 50000, -50000, lc};
Point(1100) = {-12000, 800, -10000, lc};
Point(1101) = {4000, 800, -10000, lc};
Point(1102) = {4000, 800, 0, lc};
Point(1103) = {-12000, 800, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 103};
Line(3) = {103, 100};
Line(4) = {100, 101};
Line(5) = {101, 102};
Line(6) = {102, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};
Line(9) = {103, 102};
Circle(200) = {201, 200, 202};
Circle(201) = {202, 200, 203};
Circle(202) = {203, 200, 204};
Circle(203) = {204, 200, 201};
Circle(210) = {211, 200, 212};
Circle(211) = {212, 200, 213};
Circle(212) = {213, 200, 214};
Circle(213) = {214, 200, 211};
Line(216) = {2, 1002};
Line(217) = {1002, 1012};
Line(218) = {3, 1003};
Line(219) = {1003, 1013};
Line(220) = {4, 1004};
Line(221) = {1004, 1014};
Line(222) = {1014, 1013};
Line(223) = {1004, 1003};
Line(224) = {1003, 1102};
Line(225) = {1103, 1102};
Line(226) = {1102, 1101};
Line(227) = {1101, 1100};
Line(228) = {1100, 1103};
Line(229) = {101, 1101};
Line(230) = {102, 1102};
Line(231) = {103, 1103};
Line(232) = {1103, 1002};
Line(233) = {1001, 1};
Line(234) = {1004, 1001};
Line(235) = {1011, 1012};
Line(236) = {1013, 1012};
Line(237) = {1011, 1014};
Line(238) = {1011, 1001};
Line(239) = {100, 1100};
Line(240) = {1001, 1002};
Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};
Line Loop(2) = {3, 4, 5, -9, 210, 211, 212, 213};
Plane Surface(2) = {2};
Line Loop(3) = {210, 211, 212, 213, 200, 201, 202, 203};
Plane Surface(3) = {3};
Line Loop(4) = {200, 201, 202, 203};
Plane Surface(4) = {4};
Line Loop(242) = {1, 216, -240, 233};
Plane Surface(242) = {242};
Line Loop(244) = {238, 240, 217, -235};
Plane Surface(244) = {244};
Line Loop(246) = {238, -234, 221, -237};
Plane Surface(246) = {246};
Line Loop(248) = {234, 233, -8, 220};
Plane Surface(248) = {248};
Line Loop(250) = {237, 222, 236, -235};
Plane Surface(250) = {250};
Line Loop(252) = {217, -236, -219, 224, -225, 232};
Plane Surface(252) = {252};
Line Loop(254) = {221, 222, -219, -223};
Plane Surface(254) = {254};
Line Loop(256) = {220, 223, -218, 7};
Plane Surface(256) = {256};
Line Loop(258) = {6, 218, 224, -230};
Plane Surface(258) = {258};
Line Loop(260) = {9, 230, -225, -231};
Plane Surface(260) = {260};
Line Loop(262) = {216, -232, -231, -2};
Plane Surface(262) = {262};
Line Loop(264) = {228, -231, 3, 239};
Plane Surface(264) = {264};
Line Loop(266) = {227, -239, 4, 229};
Plane Surface(266) = {266};
Line Loop(268) = {229, -226, -230, -5};
Plane Surface(268) = {268};
Line Loop(270) = {234, 240, -232, -228, -227, -226, -224, -223};
Plane Surface(270) = {270};
Line Loop(272) = {228, 225, 226, 227};
Plane Surface(272) = {272};
Line Loop(10000) = {3, 4, 5, -9};
Ruled Surface(10000) = {10000};
Surface Loop(274) = {244, 246, 254, 250, 252, 270, 272};
Volume(274) = {274};
Surface Loop(276) = {272, 268, 266, 264, 2, 260, 3, 4};
Volume(276) = {276};
Surface Loop(278) = {1, 242, 262, 248, 256, 258, 270, 268, 266, 264};
Volume(278) = {278};
Field[1] = Attractor;
Field[1].FacesList = {10000};
Field[1].FieldX = -1;
Field[1].FieldY = -1;
Field[1].FieldZ = -1;
Field[1].NNodesByEdge = 20;
Field[2] = MathEval;
//Field[2].F = Sprintf("0.05*F1 +(F1/2.5e3)^2 + %g", lc_fault);
Field[2].F = Sprintf("0.1*F1 +(F1/5.0e3)^2 + %g", lc_fault);
Background Field = 2;

Physical Surface(101) = {252, 258, 260, 262};
Physical Surface(105) = {242, 244, 246, 248, 250, 254, 256};
Physical Surface(103) = {2, 3, 4};
Physical Volume(2) = {276,278};
Physical Volume(3) = {274};

