lc  = 0.0005; // mesh size on the torch wall
lc2 = 0.005;   // mesh size on the far field
lc4 = 1.e-5; // mesh size at the coil (not used)

lengthTorch = 0.1;             // Length of the torch, along axial direction
heightTorch = 0.015;           // Radius of the torch
distance2Torch = 0.1;          // Distance between farfield and the torch
heightRefinementZone = 0.001;  // Radial size of the refinement zone near the torch wall
eps = 1.e-9;

// Coils // for display only
 // r1 = 0.019 ;
 // z1 = 0.04  ;
 // z2 = 0.053 ;
 // z3 = 0.066 ;
 // Point(5) = {z1 , r1 , 0, lc4}; 
 // Point(6) = {z2 , r1 , 0, lc4};
 // Point(7) = {z3 , r1 , 0, lc4};
 
// 4 points of the torch
Point(1) = {0, 0, 0, lc}; 
Point(2) = {0, heightTorch, 0, lc};
Point(3) = {lengthTorch, 0, 0, lc};
Point(4) = {lengthTorch, heightTorch, 0, lc};
 
// 4 points for the farfield
Point(28) = {-distance2Torch, 0, 0, lc2}; 
Point(29) = {lengthTorch+distance2Torch, 0, 0, lc2}; 
Point(30)= {lengthTorch+distance2Torch, heightTorch+distance2Torch, 0, lc2}; 
Point(31)= {-distance2Torch, heightTorch+distance2Torch, 0, lc2}; 

// Refinement zone at the torch wall
 Point(11) = {0,           heightTorch-heightRefinementZone,0,lc}; 
 Point(12) = {lengthTorch, heightTorch-heightRefinementZone,0,lc};
 Point(13) = {0,           heightTorch+heightRefinementZone,0,lc};
 Point(14) = {lengthTorch, heightTorch+heightRefinementZone,0,lc};

//***********
Line(1) = {2, 4};
Line(4) = {1, 3};
Line(10) = {29,30};
Line(11) = {30,31};
Line(12) = {31,28};
Line(6) = {28,1};
Line(7) = {3,29};
Line(15) = {1,11};
Line(16) = {11,2};
Line(17) = {2,13};
Line(18) = {13,14};
Line(19) = {14,4};
Line(20) = {4,12};
Line(21) = {12,3};
Line(22) = {11,12};

Line Loop(6) = {15,22,21,-4};
Plane Surface(6) = {6};
Line Loop(8) = {6,15,16,17,18,19,20,21,7,10,11,12};
Plane Surface(8) = {8}; // external domain
Line Loop(10) = {16,1,20,-22};
Plane Surface(10) = {10};
Line Loop(12) = {17,18,19,-1};
Plane Surface(12) = {12};

Transfinite Line {16,-17,19,-20} = 31 Using Progression 0.9; // number of nodes and progression in the refinement zone
Transfinite Line {1,4,18,22} = 401;                          // number of nodes on torch along axial direction
Transfinite Line {-15,21} = 81;                              // number of nodes on torch along radial direction
Transfinite Surface {10} = {11, 2, 4, 12} Alternate;
Transfinite Surface {12} = {2,13,14,4} Alternate;
Transfinite Surface {6} = {1,3,12,11} Alternate;

Physical Surface(1) = {8,12};   // external domain
Physical Surface(2) = {6,10};   // plasma
Physical Curve(1) = {10,11,12}; // far field
Physical Curve(2) = {6, 4, 7};  // axis
