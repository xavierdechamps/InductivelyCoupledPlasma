lc  = 0.0005; // ponderation pour le rectangle
lc2 = 0.005;   // cercle exterieur
lc3 = 1.e-4; // spire externe
lc4 = 1.e-5; // spire interne

largeurTorch = 0.1;
hauteurTorch = 0.015;
distance2torch = 0.1;
epaisseurRaffinement = 0.001;
eps = 1.e-9;

// Coils
 // r1 = 0.019 ;//+ 0.000000008;
 // z1 = 0.04  ;//+ 0.000000012 ;
 // z2 = 0.053 ;//+ 0.000000012 ;
 // z3 = 0.066 ;//+ 0.000000007 ;
 // Point(5) = {z1 , r1 , 0, lc4}; 
 // Point(6) = {z2 , r1 , 0, lc4};
 // Point(7) = {z3 , r1 , 0, lc4};
 
// 4 noeuds de la boîte rectangulaire
Point(1) = {0, 0, 0, lc}; 
Point(2) = {0, hauteurTorch, 0, lc};
Point(3) = {largeurTorch, 0, 0, lc};
Point(4) = {largeurTorch, hauteurTorch, 0, lc};
 
// points définissant le rectangle exterieur
Point(28) = {-distance2torch, 0, 0, lc2}; 
Point(29) = {largeurTorch+distance2torch, 0, 0, lc2}; 
Point(30)= {largeurTorch+distance2torch, hauteurTorch+distance2torch, 0, lc2}; 
Point(31)= {-distance2torch, hauteurTorch+distance2torch, 0, lc2}; 

// pour la zone de raffinement au niveau de l'interface air-métal
 Point(11) = {0,            hauteurTorch-epaisseurRaffinement,0,lc}; 
 Point(12) = {largeurTorch, hauteurTorch-epaisseurRaffinement,0,lc};
 Point(13) = {0,            hauteurTorch+epaisseurRaffinement,0,lc};
 Point(14) = {largeurTorch, hauteurTorch+epaisseurRaffinement,0,lc};

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
//Line(21) = {4,3};

Line Loop(6) = {15,22,21,-4};
Plane Surface(6) = {6};
Line Loop(8) = {6,15,16,17,18,19,20,21,7,10,11,12};
Plane Surface(8) = {8}; // external domain
Line Loop(10) = {16,1,20,-22};
Plane Surface(10) = {10};
Line Loop(12) = {17,18,19,-1};
Plane Surface(12) = {12};

//Point {10} In Surface {6};
Transfinite Line {16,-17,19,-20} = 31 Using Progression 0.9; //nombre de couches limites
Transfinite Line {1,4,18,22} = 401; // nombre d'éléments horizontalement
Transfinite Line {-15,21} = 81;
// Transfinite Line {65,66} = 21 Using Progression 0.8;
Transfinite Surface {10} = {11, 2, 4, 12} Alternate;
Transfinite Surface {12} = {2,13,14,4} Alternate;
Transfinite Surface {6} = {1,3,12,11} Alternate;

// Physical Surface(1) = {8,12,20,21,22,23,24,25}; // external domain
Physical Surface(1) = {8,12}; // external domain
Physical Surface(2) = {6,10}; // plasma
Physical Curve(1) = {10,11,12}; // far field
Physical Curve(2) = {6, 4, 7}; //axis
