#     include "IBC/commonLaws1.h" 
// loi de Musker 
press_vec[noind ]   = pext;
yplus_vec[noind ]   = roext*yibc/muext;// yplus/utau
ro_vec[noind ]      = roext;
utau_vec[noind ]    = utau0;
aa_vec[noind]       = roext*yext/muext;
uext_vec[noind]     = uext;
nutcible_vec[noind] = utau0;
sign_vec[noind]     = signibc/uext;
ucible_vec[noind]   = alphasbeta*un; // init : normal component of velocity is linearly reconstructed
vcible_vec[noind]   = alphasbeta*vn;
wcible_vec[noind]   = alphasbeta*wn;
ut_vec[noind]       = ut;
vt_vec[noind]       = vt;
wt_vec[noind]       = wt;
mu_vec[noind]       = muext;
alpha_vec[noind]    = alpha;
tcible_vec[noind]   = text;
#  include "IBC/musker_vec.h"
// out= utau  et err 
