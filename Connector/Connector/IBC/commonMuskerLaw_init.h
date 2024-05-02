#     include "IBC/commonLaws1.h" 
// loi de Musker 


// Crocco-Busemann relationship
//Ta = Tb+0.5*pow(Pr,one_third)/(cv*gamma)*(ub^2-ua^2)
//In this case a=wall & b=ext(image)
twall  = text + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext*uext);               // Temperature a la paroi (assumes wall has 0 velocity)
rowall = pext/twall*cvgaminv;                                               // Densite a la paroi
// muwall = coefSuth * sqrt(K_FUNC::E_abs(twall)*Tsinv) / (1.+Cs/twall);    // Viscosite a la paroi en utilisant temperature absolue pour reference
muwall = muext*sqrt(twall/text)*(1+Cs/text)/(1+Cs/twall);                   // Viscosite a la paroi en utilisant temperature exterieure (PI) pour reference (Benjamin's formula)

utau0  = sqrt(muwall*uext/(yext*rowall));

press_vec[noind ]   = pext;

// yplus_vec[noind ]   = roext*yibc/muext;// yplus/utau
yplus_vec[noind ]   = rowall*yibc/muwall;// yplus/utau

// ro_vec[noind ]      = roext;
ro_vec[noind ]      = rowall;

utau_vec[noind ]    = utau0;

// aa_vec[noind]       = roext*yext/muext;
aa_vec[noind]       = rowall*yext/muwall;

uext_vec[noind]     = uext;
nutcible_vec[noind] = utau0;
sign_vec[noind]     = signibc/uext;
ucible_vec[noind]   = alphasbeta*un; // init : normal component of velocity is linearly reconstructed
vcible_vec[noind]   = alphasbeta*vn;
wcible_vec[noind]   = alphasbeta*wn;
ut_vec[noind]       = ut;
vt_vec[noind]       = vt;
wt_vec[noind]       = wt;

// mu_vec[noind]       = muext;
mu_vec[noind]       = muwall;

alpha_vec[noind]    = alpha;
tcible_vec[noind]   = text;
#  include "IBC/musker_vec.h"
// out= utau  et err 
