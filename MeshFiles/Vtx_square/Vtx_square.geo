h=0.1;
//h=0.05;
//h=0.025;
//h=0.0125;
//h=0.00625;
//+
SetFactory("OpenCASCADE");
//+
Point(1)={0, 0, 0, h};
//+
Point(2)={1, 0, 0, h};
//+
Point(3)={1, 1, 0, h};
//+
Point(4)={0, 1, 0, h};
//+
Point(5)={0.5, 0.5, 0, h/8}; // Le point à ajouter au maillage
//+
Line(1) ={1,2};
//+
Line(2) ={2,3};
//+
Line(3) ={3,4};
//+
Line(4) ={4,1};
//+ 
Line Loop(7)={1,2,3,4};
//+ 
Plane Surface(8)= {7}; 
//+
Point{5} In Surface{8};
//+
Physical Line("bas")={1};
//+
Physical Line("Droit")={2};
//+
Physical Line("haut")={3};
//+
Physical Line("Gauche")={4};
//+
Physical Surface("dom") = {8};   // Setting a label to the Surface

