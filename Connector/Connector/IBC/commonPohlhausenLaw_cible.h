err = 1; skip =0;

#ifdef _OPENM4
#pragma omp simd
#endif

for (E_Int noind = 0; noind < ifin-ideb; noind++)
{
  yibc             = mu_vec[noind]*yplus_vec[noind]/ro_vec[noind];
  yplus            = utau_vec[noind]*yplus_vec[noind];
  yplus_vec[noind] = yplus;

  denoml10 = yplus*yplus-8.15*yplus+86.;
  denoml10 = denoml10*denoml10;

  // uext: norme de la composante tangentielle de la vitesse externe
  uext = sqrt(ut_vec[noind]*ut_vec[noind]+vt_vec[noind]*vt_vec[noind]+wt_vec[noind]*wt_vec[noind]);
  uext = std::max(uext, 1.e-12);

  delta = 5.84/sqrt(uext*ro_vec[noind]/(mu_vec[noind]*xPW[noind+ideb]));
  eta   = yibc/delta;

  // if (eta < 0.){
  //   std::cout << "eta = "   << eta << std::endl;
  //   std::cout << "yibc = "  << yibc << std::endl;
  //   std::cout << "delta = " << delta << std::endl;
  // }
  
  if (eta < 1.){
    ut_vec[noind] *= (2*eta - 2*pow(eta, 3) + pow(eta, 4));
    vt_vec[noind] *= (2*eta - 2*pow(eta, 3) + pow(eta, 4));
    wt_vec[noind] *= (2*eta - 2*pow(eta, 3) + pow(eta, 4));
  }
  
  // mod(ut_cible) = mod(ut_image)*solutionBlasius_cible
  umod = sqrt(ut_vec[noind]*ut_vec[noind]+vt_vec[noind]*vt_vec[noind]+wt_vec[noind]*wt_vec[noind]);
  umod = K_FUNC::E_abs(umod);

  ucible0 = sign_vec[noind];
  ucible_vec[noind] += ucible0 * ut_vec[noind]; // vitesse tangentielle pour le pt IBC
  vcible_vec[noind] += ucible0 * vt_vec[noind];
  wcible_vec[noind] += ucible0 * wt_vec[noind];

  tcible_vec[noind] = tcible_vec[noind] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext*uext - umod*umod); // Crocco-Busemann
}
