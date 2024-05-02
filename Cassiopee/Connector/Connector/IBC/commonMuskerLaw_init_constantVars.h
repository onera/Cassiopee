b0 = xPI[noind]-xPW[noind];
b1 = yPI[noind]-yPW[noind];
b2 = zPI[noind]-zPW[noind];

normb = sqrt(b0*b0+b1*b1+b2*b2);
normb = std::max(normb, 1.e-12);
n0 = b0/normb;
n1 = b1/normb;
n2 = b2/normb;

vn = u*n0+v*n1+w*n2;
yext  = b0*n0+b1*n1+b2*n2; //distance yext du point interpole a la paroi

// uscaln: u au point interpole scalaire la normale au pt interpole
uscaln = u*n0 + v*n1 + w*n2;
  
//composante normale de la vitesse
un = uscaln*n0;
vn = uscaln*n1;
wn = uscaln*n2;

//composante tangentielle de la vitesse au pt interpole
ut = u-un;
vt = v-vn;
wt = w-wn;
// uext: norme de la composante tangentielle de la vitesse externe
uext = sqrt(ut*ut+vt*vt+wt*wt);
uext = std::max(uext, 1.e-12);

// loi de Musker 
rowall = roext;              
muwall = muext;
utau0  = sqrt(muwall*uext/(yext*rowall));

press_vec[noind ]   = pext;
yplus_vec[noind ]   = rowall*yext/muwall;
ro_vec[noind ]      = rowall;
utau_vec[noind ]    = utau0;
aa_vec[noind]       = rowall*yext/muwall;
uext_vec[noind]     = uext;
nutcible_vec[noind] = utau0;
ut_vec[noind]       = ut;
vt_vec[noind]       = vt;
wt_vec[noind]       = wt;
mu_vec[noind]       = muwall;

#  include "IBC/musker_vec.h"
// out= utau  et err 
