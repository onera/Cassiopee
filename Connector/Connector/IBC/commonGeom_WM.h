// PC: point corrige
// PW: point paroi
// PI: point interpole (ou exterieur)
// IN: (u,v,w) au point interpole
a0 = xPC[noind+ideb]-xPW[noind+ideb];
a1 = yPC[noind+ideb]-yPW[noind+ideb];
a2 = zPC[noind+ideb]-zPW[noind+ideb];

normb = sqrt(a0*a0+a1*a1+a2*a2);
normb = std::max(normb, 1.e-12);
n0 = a0/normb;
n1 = a1/normb;
n2 = a2/normb;
