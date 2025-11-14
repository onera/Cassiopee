err = 1; skip =0;

//initialisation Newton SA  + vitesse cible

#ifdef _OPENM4
#pragma omp simd
#endif

for (E_Int noind = 0; noind < ifin-ideb; noind++)
{
  utau_vec[noind] = std::max(utau_vec[noind], 1.e-12);

  if (gradP_vec[noind] == 0.) 
  {
    gradP_vec[noind] = 0.;
    utau_vec[noind] = utauOri_vec[noind];

    yplus            = utau_vec[noind]*yplus_vec[noind];
    yplus_vec[noind] = yplus;

    denoml10 = yplus*yplus-8.15*yplus+86.;
    denoml10 = denoml10*denoml10;

    px  = gradP_vec[noind]/pow(utau_vec[noind],3);
      
    l11 = pow(yplus+10.6,9.6);
    l12 = yplus*(yplus-8.15) + 86.;
    l13 = (2.*yplus-8.15)/16.7;

    l1 = 5.424*atan(l13) + log10(l11/(l12*l12)) - 3.52;
    l2 = -2.*kappainv*log((sqrt(1.+px*yplus)+1.)/2.);
    l3 = 2.*kappainv*(sqrt(1.+px*yplus)-1.);

    umod = utau_vec[noind]*(l1 + l2 + l3);
  }

  else 
  { 
    //Mafzal s'applique - FULL VERSION
    yplus = utau_vec[noind]*yplus_vec[noind];
    yplus_vec[noind] = yplus;

    denoml10 = yplus*yplus-8.15*yplus+86.;
    denoml10 = denoml10*denoml10;

    px  = gradP_vec[noind]/pow(utau_vec[noind],3);
    
    l11 = pow(yplus+10.6,9.6);
    l12 = yplus*(yplus-8.15) + 86.;
    l13 = (2.*yplus-8.15)/16.7;

    l1 = 5.424*atan(l13) + log10(l11/(l12*l12)) - 3.52;
    l2 = 0.;
    l3 = 0.;

    if (px > 0.){
      if (MafzalMode == 1){
        l2 = -2.*kappainv*log((sqrt(1.+px*yplus)+1.)/2.); //PRESSURE
        l3 = 2.*kappainv*(sqrt(1.+px*yplus)-1.); //PRESSURE
      }
      else{
        l2 = 2.*kappainv*log((sqrt(1.+px*yplus)+1.)/2.); //PRESSURE
        l3 = 0.; //PRESSURE
      }
    }
    else{
      if (MafzalMode == 3){
        px = -px;
        l2 = -2.*kappainv*log((sqrt(1.+px*yplus)+1.)/2.); //PRESSURE
        l3 = 0.; //PRESSURE
      }
      else{
        l2 = 0.; //PRESSURE
        l3 = 0.; //PRESSURE
      }
    }

    umod = utau_vec[noind]*(l1 + l2 + l3);
  }
 
  umod = K_FUNC::E_abs(umod);

  ucible0 = sign_vec[noind] * umod;
  ucible_vec[noind] += ucible0 * ut_vec[noind]; // vitesse tangentielle pour le pt IBC
  vcible_vec[noind] += ucible0 * vt_vec[noind];
  wcible_vec[noind] += ucible0 * wt_vec[noind];

  // uext: norme de la composante tangentielle de la vitesse externe
  uext = sqrt(ut_vec[noind]*ut_vec[noind]+vt_vec[noind]*vt_vec[noind]+wt_vec[noind]*wt_vec[noind]);
  uext = std::max(uext, 1.e-12);


  tcible_vec[noind] = tcible_vec[noind] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext*uext - umod*umod); // Crocco-Busemann

  expy                = 1.-exp(-yplus/19.);// ranges 17 a 26
  nutcible_vec[noind] = (kappa * alpha_vec[noind])*utau_vec[noind ] * expy*expy;//negatif si pt ibc interieur aux corps

  nutilde             = K_FUNC::E_abs( nutcible_vec[noind] );
  aa_vec[noind]       = nutilde;

# include "IBC/fnutilde_vec.h"

  // out= mut  et err
  ut_vec[noind]    = nutilde; // Sauvegarde du nutilde, au cas ou newton non convergent.
}
