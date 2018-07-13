//out: alpha, beta, ni, vn
# include "IBC/commonGeom.h"
  
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
