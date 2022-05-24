#ifdef _OPENM4
#pragma omp simd
#endif

E_Float tol = 1.e-12;
E_Float Cv1 = 7.1;

for (E_Int noind = 0; noind < ifin-ideb; noind++)
{
  // printf(" noind = %d | utau = %g \n", noind, utau_vec[noind]);
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
  }
  else{
    lambda_thwaites = 0.;
    delta = 5.84/sqrt(uext*ro_vec[noind]/(mu_vec[noind]*xPW[noind+ideb]));
  }

  eta   = yibc/delta;

  if (eta < 1.){
    ut_vec[noind] = ut_vec[noind]*(2*eta - 2*pow(eta, 3) + pow(eta, 4)) + 1./6.*lambda_thwaites*eta*pow((1-eta), 3);
    vt_vec[noind] = vt_vec[noind]*(2*eta - 2*pow(eta, 3) + pow(eta, 4)) + 1./6.*lambda_thwaites*eta*pow((1-eta), 3);
    wt_vec[noind] = wt_vec[noind]*(2*eta - 2*pow(eta, 3) + pow(eta, 4)) + 1./6.*lambda_thwaites*eta*pow((1-eta), 3);
  }

  // mod(ut_cible) = mod(ut_image)*solutionBlasius_cible
  umod = sqrt(ut_vec[noind]*ut_vec[noind]+vt_vec[noind]*vt_vec[noind]+wt_vec[noind]*wt_vec[noind]);
  umod = K_FUNC::E_abs(umod);

  ucible0 = sign_vec[noind];
  ucible_vec[noind] += ucible0 * ut_vec[noind]; // vitesse tangentielle pour le pt IBC
  vcible_vec[noind] += ucible0 * vt_vec[noind];
  wcible_vec[noind] += ucible0 * wt_vec[noind];

  tcible_vec[noind] = tcible_vec[noind] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext*uext - umod*umod); // Crocco-Busemann

  //twall = tcible_vec[noind] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext*uext);
  //tcible_vec[noind] =  twall + (tcible_vec[noind] + 0.5*(uext*uext)/(cv*gamma) - twall)*(umod/uext) - 0.5*(umod*umod)/(cv*gamma); // Equations de Crocco, plus precises pour ecoulements compressibles

  // van dryst pour nut
  expy                = 1.-exp(-yplus/19.);// ranges 17 a 26
  nutcible_vec[noind] = (kappa * alpha_vec[noind])*utau_vec[noind] * expy*expy;//negatif si pt ibc interieur aux corps

  E_Float nutcible = K_FUNC::E_abs( nutcible_vec[noind] );

  // equation 4eme degre
  E_Float a = nutcible;
  E_Float b = nutcible*pow( (mu_vec[noind]/ro_vec[noind])*Cv1, 3.);
  //printf("x^4 + %g x^3 +%g = 0\n", -a,-b);

  // debug
  //a = 1.;
  //b = -12.;

  // changement de variable 4eme degre
  E_Float p = -3*a*a/8.;
  E_Float q = -a*a*a/8.;
  E_Float r = -3*pow(a/4.,4.)-b;

  // equation 3eme degre
  E_Float ap = 8.;
  E_Float bp = -4*p;
  E_Float cp = -8*r;
  //E_Float dp = (a*a)/2.*(pow(a/4.,4.)+3*b);
  E_Float dp = 4*p*r-q*q;
  //printf("cubique: %g x^3 + %g x^2 + %g x + %g\n",ap,bp,cp,dp);

  // racines du 3eme degre
  E_Float delta = 18*ap*bp*cp*dp-4*bp*bp*bp*dp+bp*bp*cp*cp-4*ap*cp*cp*cp-27*ap*ap*dp*dp;
  E_Float delta0 = bp*bp-3*ap*cp;
  E_Float delta1 = 2*bp*bp*bp-9*ap*bp*cp+27*ap*ap*dp;

  E_Float superdelta = -27.*ap*ap*delta;

  //printf("delta %g superdelta>0 = %g\n", delta, superdelta);


  E_Float y1 = -1.;
  E_Float y2 = -1.;
  if (fabs(delta) < tol && fabs(delta0) < tol)
  {
    y1 = -bp/(3*ap);  //printf("racine y1=%g\n", y1);
  }
  else if (fabs(delta) < tol)
  {
    y1 = (9*ap*dp - bp*cp)/(2*delta0);
    y2 = (4*ap*bp*cp-9*ap*ap*dp-bp*bp*bp)/(ap*delta0);
    //printf("racine y1=%g y2=%g\n", y1, y2);
  }
  else
  {
    /* version super delta */
    E_Float C1 = -1.; E_Float C2 = -1.;
    if (superdelta >= 0.)
    {
      E_Float root = sqrt(superdelta);
      if (delta1-root >= 0.)
        C1 = pow( (delta1 -root) /2., 1./3. );
      else
        C1 = -pow( (root -delta1) /2., 1./3. );
      if (delta1 +root >= 0)
        C2 = pow( (delta1 +root) /2., 1./3. );
      else
        C2 = -pow( -(delta1 +root) /2., 1./3. );
    }

    //printf("C1=%g, C2=%g\n", C1, C2);
    y1 = -1./(3*ap)*(bp + C1 + delta0/C1 );
    y2 = -1./(3*ap)*(bp + C2 + delta0/C2 );
    //printf("racine y1=%g y2=%g\n", y1, y2);
  }


  // racine de l'equation du 4eme degre
  E_Float c1 = 2*y1-p;
  E_Float c2 = 2*y2-p;
  //printf("c1 > 0 = %g, c2 > 0 = %g\n",c1,c2);

  E_Float z1 = -123456.;
  if (c1 >= tol)
  {
    E_Float p1 = -2*y1-p+2*q/(sqrt(c1));
    E_Float p2 = -2*y1-p-2*q/(sqrt(c1));
    //printf("1. p1=%g p2=%g\n", p1,p2);
    if (p1 >= tol) { z1 = 0.5*(sqrt(c1)+sqrt(p1));}
    else if (p2 >= tol) { z1 = 0.5*(sqrt(c1)+sqrt(p2));}
  }
  if (c2 >= tol && z1 == -123456)
  {
    E_Float p1 = -2*y2-p+2*q/(sqrt(c2));
    E_Float p2 = -2*y2-p-2*q/(sqrt(c2));
    //printf("2. p1=%g p2=%g\n", p1,p2);
    if (p1 >= tol) { z1 = 0.5*(sqrt(c2)+sqrt(fabs(p1)));}
    else if (p2 >= tol) { z1 = 0.5*(sqrt(c2)+sqrt(fabs(p2)));}
  }
  if (c1 <= tol && z1 == -123456)
  {
    E_Float b0 = y1*y1-r;
    if (b0 >= tol)
    {
      E_Float p1 = -2*y1-p+4.*sqrt(b0);
      E_Float p2 = -2*y1-p-4.*sqrt(b0);
      //printf("3. p1=%g p2=%g\n", p1,p2);
      if (p1 >= tol) { z1 = 0.5*(sqrt(fabs(c1))+sqrt(fabs(p1))); }
      else if (p2 >= tol) { z1 = 0.5*(sqrt(fabs(c1))+sqrt(fabs(p2))); }
    }
  }
  if (c2 <= tol && z1 == -123456)
  {
    E_Float b0 = y2*y2-r;
    if (b0 >= tol)
    {
      E_Float p1 = -2*y2-p+4.*sqrt(b0);
      E_Float p2 = -2*y2-p-4.*sqrt(b0);
      //printf("4. p1=%g p2=%g\n", p1,p2);
      if (p1 >= tol) { z1 = 0.5*(sqrt(fabs(c2))+sqrt(fabs(p1))); }
      else if (p2 >= tol) { z1 = 0.5*(sqrt(fabs(c2))+sqrt(fabs(p2))); }
    }
  }

  // nutile final
  E_Float nutilde1 = z1 + a/4.;
  aa_vec[noind] = nutilde1;
  //printf("nutilde final = %g\n", nutilde1);

}
