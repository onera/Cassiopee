// Pour Initialisation TBLE 

// PC: point linelets
// PW: point paroi
// PI: point interpole (ou exterieur)
// IN: (u,v,w) au point interpole
// a0 = xline[indl]-xPW[noind+ideb];
// a1 = yline[indl]-yPW[noind+ideb];
// a2 = zline[indl]-zPW[noind+ideb];
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
// if (K_FUNC::E_abs(beta)<1.e-12) beta = 1.e-12;
// alphasbeta = alpha/beta;


yext = beta; //distance yext du point interpole a la paroi
yibc = K_FUNC::E_abs(alpha); //si le pt IBC est interieur au corps cette distance est positive
if (yibc < 1.e-12)
{
  yibc = 1.e-12; signibc = 0;
}
else 
{
  signibc = alpha/yibc;//signe de la distance pt IBC : signibc=-1: pt IBC interieur au corps  
}

// Loi de Sutherland -> viscosite au point interpole
muext = coefSuth * sqrt(K_FUNC::E_abs(text)*Tsinv) / (1.+Cs/text);

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

// Calcul du frottement: loi "lineaire" (par defaut)
utau0 = sqrt(muext*uext/(yext*roext));
// if ( noind == 0 ) printf(" muext = %5.10f | uext = %5.10f \n", muext, uext); 


// loi de Musker 
press_vec[noind ]   = pext;
// yplus_vec[noind ]   = roext*yibc/muext;// yplus/utau
yplus_vec[noind ]   = roext*yline[indl]/muext;// yplus/utau
ro_vec[noind ]      = roext;
utau_vec[noind ]    = utau0;
aa_vec[noind]       = roext*yext/muext;
uext_vec[noind]     = uext;
nutcible_vec[noind] = utau0;
sign_vec[noind]     = signibc/uext;
ucible_vec[noind]   = yline[indl]/normb*un;//alphasbeta_line[noind + ideb]*un;//yline[indl]/normb*un; // init : normal component of velocity is linearly reconstructed
vcible_vec[noind]   = yline[indl]/normb*vn;//alphasbeta_line[noind + ideb]*vn;//yline[indl]/normb*vn;
wcible_vec[noind]   = yline[indl]/normb*wn;//alphasbeta_line[noind + ideb]*wn;//yline[indl]/normb*wn;
ut_vec[noind]       = ut;
vt_vec[noind]       = vt;
wt_vec[noind]       = wt;
mu_vec[noind]       = muext;
alpha_vec[noind]    = alpha;
tcible_vec[noind]   = text;
#  include "IBC/musker_vec.h"
// out= utau  et err 
