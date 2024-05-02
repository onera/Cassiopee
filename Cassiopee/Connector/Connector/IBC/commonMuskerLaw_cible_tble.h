// Pour Initialisation TBLE

err = 1; skip =0;

//initialisation Newton SA  + vitesse cible

#ifdef _OPENM4
#pragma omp simd
#endif 

for (E_Int noind = 0; noind < ifin-ideb; noind++)
{
  // printf(" noind = %d | utau = %g \n", noind, utau_vec[noind]);
  yplus            = utau_vec[noind]*yplus_vec[noind];
  yplus_vec[noind] = yplus;
  
  denoml10 = yplus*yplus-8.15*yplus+86.;
  denoml10 = denoml10*denoml10;
  
  umod = utau_vec[noind ]*(5.424*atan((2.*yplus-8.15)/16.7) + log10(pow(yplus+10.6,9.6)/denoml10) - 3.52);
  umod = K_FUNC::E_abs(umod);

  ucible0 = sign_vec[noind] * umod;
  ucible_vec[noind] = ucible0 * ut_vec[noind]; // vitesse tangentielle pour le pt linelets
  vcible_vec[noind] = ucible0 * vt_vec[noind];
  wcible_vec[noind] = ucible0 * wt_vec[noind];


  //E_Int indR = rcvPts[noind+ideb];

  // uext: norme de la composante tangentielle de la vitesse externe
  uext = sqrt(ut_vec[noind]*ut_vec[noind]+vt_vec[noind]*vt_vec[noind]+wt_vec[noind]*wt_vec[noind]);
  uext = std::max(uext, 1.e-12);

  tcible_vec[noind] = tcible_vec[noind] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext*uext - umod*umod); // Crocco-Busemann

  // printf(" Doit valoir 1 : %g \n",  sqrt(mu_vec[noind]*umod/yibc/ro_vec[noind])/utau_vec[noind]);
  
  expy                = 1.-exp(-yplus/19.);// ranges 17 a 26
  nutcible_vec[noind] = (kappa * alpha_vec[noind])*utau_vec[noind ] * expy*expy;//negatif si pt ibc interieur aux corps

  nutilde             = K_FUNC::E_abs( nutcible_vec[noind] );
  aa_vec[noind]       = nutilde;

// # include "IBC/fnutilde_vec.h"

  // out= mut  et err 
  // ut_vec[noind]    = nutilde; // Sauvegarde du nutilde, au cas ou newton non convergent.
}
