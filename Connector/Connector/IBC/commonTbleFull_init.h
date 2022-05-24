#     include "IBC/commonLaws1.h" 

twall  = text + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext*uext);               // Temperature a la paroi
rowall = pext/twall*cvgaminv;                                               // Densite a la paroi
muwall = muext*sqrt(twall/text)*(1+Cs/text)/(1+Cs/twall);                   // Viscosite a la paroi en utilisant temperature exterieure (PI) pour reference (Benjamin's formula)

utau0  = sqrt(muwall*uext/(yext*rowall));

press_vec[noind ]   = pext;

yplus_vec[noind ]   = rowall*yibc/muwall;// yplus/utau

ro_vec[noind ]      = rowall;

utau_vec[noind ]    = utau0;

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

mu_vec[noind]       = muwall;

alpha_vec[noind]    = alpha;
tcible_vec[noind]   = text;

unext               = sqrt(un*un+vn*vn+wn*wn);
unext               = std::max(unext, 1.e-12);

// du/dt
tgradU =          (gradxUext*ut/uext  + gradyUext*vt/uext  + gradzUext*wt/uext) *ut/uext;
tgradU = tgradU + (gradxVext*ut/uext  + gradyVext*vt/uext  + gradzVext*wt/uext) *vt/uext;
tgradU = tgradU + (gradxWext*ut/uext  + gradyWext*vt/uext  + gradzWext*wt/uext) *wt/uext;

// du/dn
ngradU =          (gradxUext*un/unext + gradyUext*vn/unext + gradzUext*wn/unext)*ut/uext;
ngradU = ngradU + (gradxVext*un/unext + gradyVext*vn/unext + gradzVext*wn/unext)*vt/uext;
ngradU = ngradU + (gradxWext*un/unext + gradyWext*vn/unext + gradzWext*wn/unext)*wt/uext;

// Source term for TBLE
gradP_vec[noind]    = (gradxPext*ut/uext+gradyPext*vt/uext+gradzPext*wt/uext)/(rowall); // Pressure component : (1/rho)*dp/dt
conv_vec[noind]     = uext*tgradU + unext*ngradU; // convective terms : u*du/dt + v*du/dn

// // precaution trailing edge
// if (xPC[noind+ideb] > 0.99){
//   gradP_vec[noind] = 0.;
//   conv_vec[noind]  = 0.;
// }

#  include "IBC/musker_vec.h"
