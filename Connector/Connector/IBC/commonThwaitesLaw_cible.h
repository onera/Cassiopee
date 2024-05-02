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

  // angle d'attaque nul
  utinf = 1.;
  vtinf = 0.;
  wtinf = 0.;

  m_thwaites = acos(utinf*ut_vec[noind]/uext + vtinf*vt_vec[noind]/uext + wtinf*wt_vec[noind]/uext);

  if (yPW[noind+ideb] >= 0 &&  vt/uext < vtinf){
    m_thwaites = -m_thwaites;
  }
  else if (yPW[noind+ideb] <= 0 &&  vt/uext > vtinf){
    m_thwaites = -m_thwaites;
  }

  m_thwaites = m_thwaites/(M_PI - m_thwaites);

  lambda_thwaites = 0.45*m_thwaites/(5.*m_thwaites+1);

  if (lambda_thwaites >= 0.){
    t_thwaites = 0.22 + 1.57*lambda_thwaites - 1.8*lambda_thwaites*lambda_thwaites;
  }
  else{
    t_thwaites = 0.22 + 1.402*lambda_thwaites - 0.018*lambda_thwaites/(0.107+lambda_thwaites);
  }

  if (m_thwaites > 0.){
    a_thwaites = lambda_thwaites/36.;
    b_thwaites = 2.*lambda_thwaites/3. - pow(b_thwaites,2);
    c_thwaites = 4.*lambda_thwaites;
    c_thwaites = pow(b_thwaites,2) - 4*a_thwaites*c_thwaites;

    lambda_thwaites = (-b_thwaites-sqrt(c_thwaites))/(2.*a_thwaites);
    delta = sqrt((mu_vec[noind]*xPW[noind+ideb]*lambda_thwaites)/(uext*ro_vec[noind]*m_thwaites));
    utau_thwaites   = sqrt(sqrt(pow(uext,3)*mu_vec[noind]*(5*m_thwaites+1)/(xPW[noind+ideb]*ro_vec[noind]*0.45))*t_thwaites);
  }
  else{
    lambda_thwaites = 0.;
    delta = 5.84/sqrt(uext*ro_vec[noind]/(mu_vec[noind]*xPW[noind+ideb]));
    utau_thwaites   = sqrt(sqrt(pow(uext,3)*mu_vec[noind]/(xPW[noind+ideb]*ro_vec[noind]))*0.343);
  }

  eta   = yibc/delta;

  if (eta < 1.){
    ut_thwaites = ut_vec[noind]*(2*eta - 2*pow(eta, 3) + pow(eta, 4)) + 1./6.*lambda_thwaites*eta*pow((1-eta), 3);
    vt_thwaites = vt_vec[noind]*(2*eta - 2*pow(eta, 3) + pow(eta, 4)) + 1./6.*lambda_thwaites*eta*pow((1-eta), 3);
    wt_thwaites = wt_vec[noind]*(2*eta - 2*pow(eta, 3) + pow(eta, 4)) + 1./6.*lambda_thwaites*eta*pow((1-eta), 3);
  }
  else{
    ut_thwaites = ut_vec[noind];
    vt_thwaites = vt_vec[noind];
    wt_thwaites = wt_vec[noind];
  }

  // umod = sqrt(ut_vec[noind]*ut_vec[noind]+vt_vec[noind]*vt_vec[noind]+wt_vec[noind]*wt_vec[noind]);
  umod = sqrt(ut_thwaites*ut_thwaites+vt_thwaites*vt_thwaites+wt_thwaites*wt_thwaites);
  umod = K_FUNC::E_abs(umod);

  ucible0 = sign_vec[noind]*uext;
  ucible_vec[noind] += ucible0 * ut_thwaites; // vitesse tangentielle pour le pt IBC
  vcible_vec[noind] += ucible0 * vt_thwaites;
  wcible_vec[noind] += ucible0 * wt_thwaites;

  tcible_vec[noind] = tcible_vec[noind] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext*uext - umod*umod); // Crocco-Busemann
}
