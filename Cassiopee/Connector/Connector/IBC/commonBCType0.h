# include "IBC/commonGeom.h"

vnc    = vn * alphasbeta;   // v.n du pt corrige
ucible = (u - vn*n0) + vnc*n0;// u du pt corrige
vcible = (v - vn*n1) + vnc*n1;// v du pt corrige
wcible = (w - vn*n2) + vnc*n2;// w du pt corrige

// ronutildeSA: lineaire
//if (nvars == 6) varSAOut[indR] = varSAIn[noind]*alphasbeta;
