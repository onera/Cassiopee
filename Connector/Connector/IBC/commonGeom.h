// PC: point corrige
// PW: point paroi
// PI: point interpole (ou exterieur)
// IN: (u,v,w) au point interpole
a0 = xPC[noind+ideb]-xPW[noind+ideb];
a1 = yPC[noind+ideb]-yPW[noind+ideb];
a2 = zPC[noind+ideb]-zPW[noind+ideb];

b0 = xPI[noind+ideb]-xPW[noind+ideb];
b1 = yPI[noind+ideb]-yPW[noind+ideb];
b2 = zPI[noind+ideb]-zPW[noind+ideb];

normb = sqrt(b0*b0+b1*b1+b2*b2);
normb = std::max(normb, 1.e-12);
n0 = b0/normb;
n1 = b1/normb;
n2 = b2/normb;

vn = u*n0+v*n1+w*n2;

alpha = a0*n0+a1*n1+a2*n2;
beta  = b0*n0+b1*n1+b2*n2;
if (K_FUNC::E_abs(beta)<1.e-12) beta = 1.e-12;
alphasbeta = alpha/beta;
