# include "IBC/commonGeom.h"
ucible = u*alphasbeta;// u du pt corrige
vcible = v*alphasbeta;// v du pt corrige
wcible = w*alphasbeta;// w du pt corrige

# include "IBC/LBM/commonGeom11Dens.h"
rhocible = rho*(-gamma_local+alphasbeta+alphasbeta*gamma_local)/(1.-eta_local+alphasbeta*eta_local);// rho du pt corrige
