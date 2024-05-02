#     include "IBC/commonLaws1.h"

twall  = text + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext*uext);               // Temperature a la paroi
rowall = pext/twall*cvgaminv;                                               // Densite a la paroi
// muwall = coefSuth * sqrt(K_FUNC::E_abs(twall)*Tsinv) / (1.+Cs/twall);    // Viscosite a la paroi en utilisant temperature absolue pour reference
muwall = muext*sqrt(twall/text)*(1+Cs/text)/(1+Cs/twall);                   // Viscosite a la paroi en utilisant temperature exterieure (PI) pour reference (Benjamin's formula)

utau0  = sqrt(sqrt(pow(uext,3)*muwall/(xPW[noind+ideb]*rowall))*0.343);

press_vec[noind]      = pext;

ro_vec[noind]         = rowall;

utau_vec[noind]       = sqrt(sqrt(pow(uext,3)*muwall/(xPW[noind+ideb]*rowall))*0.343);
aa_vec[noind]         = rowall*yext/muwall;
yplus_vec[noind]      = rowall*yibc/muwall;

uext_vec[noind]       = uext;
nutcible_vec[noind]   = utau0;
sign_vec[noind]       = signibc;
ucible_vec[noind]     = alphasbeta*un;
vcible_vec[noind]     = alphasbeta*vn;
wcible_vec[noind]     = alphasbeta*wn;
ut_vec[noind]         = ut;
vt_vec[noind]         = vt;
wt_vec[noind]         = wt;

mu_vec[noind]         = muwall;

alpha_vec[noind]      = alpha;
tcible_vec[noind]     = text;
