#     include "IBC/commonLaws1.h"
// loi de Musker + correction afzal 

twall  = text + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext*uext);               // Temperature a la paroi
rowall = pext/twall*cvgaminv;                                               // Densite a la paroi
muwall = muext*sqrt(twall/text)*(1+Cs/text)/(1+Cs/twall);                   // Viscosite a la paroi en utilisant temperature exterieure (PI) pour reference (Benjamin's formula)

utau0  = sqrt(muwall*K_FUNC::E_abs(uext)/(yext*rowall));

press_vec[noind ]   = pext;

yplus_vec[noind ]   = rowall*yibc/muwall; // yplus/utau

ro_vec[noind ]      = rowall;

utauOri_vec[noind]  = utau0; // musker law
utau_vec[noind ]    = utau0; // mafzal law

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

unext               = sqrt(un*un+vn*vn+wn*wn);
unext               = std::max(unext, 1.e-12);

gradP_vec[noind]    = (gradxPext*ut/uext+gradyPext*vt/uext+gradzPext*wt/uext)*(muwall/(rowall*rowall));

if (MafzalMode < 3){
	if (gradP_vec[noind] < 0.){
		gradP_vec[noind] = 0.;
	}
}

mu_vec[noind]       = muwall;

alpha_vec[noind]    = alpha;
tcible_vec[noind]   = text;
// #  include "IBC/mafzal_vec.h"
#  include "IBC/musker_vec.h"
