E_Float tol = 1.e-12;
E_Float Cv1 = 7.1;

// van driest pour nut
expy                = 1.-exp(-yplus/19.);// ranges 17 a 26
E_Float nutcible = (kappa * yline_ferrari)*utau_vec[noind] * expy*expy;//negatif si pt ibc interieur aux corps
nutcible = K_FUNC::E_abs(nutcible);
// equation 4eme degre
E_Float nu = (mu_vec[noind]/ro_vec[noind]);
E_Float a = nutcible/nu;
E_Float b = a*Cv1*Cv1*Cv1;

// changement de variable 4eme degre
E_Float p = -3*a*a/8.;
E_Float q = -a*a*a/8.;
E_Float r = -3*a*a*a*a/256.-b;

// equation 3eme degre
E_Float ap = 8.;
E_Float bp = -4*p;
E_Float cp = -8*r;
E_Float dp = 4*p*r-q*q;

// racines du 3eme degre
E_Float delta = 18*ap*bp*cp*dp-4*bp*bp*bp*dp+bp*bp*cp*cp-4*ap*cp*cp*cp-27*ap*ap*dp*dp;
E_Float delta0 = bp*bp-3*ap*cp;
E_Float delta1 = 2*bp*bp*bp-9*ap*bp*cp+27*ap*ap*dp;

E_Float superdelta = -27.*ap*ap*delta;

E_Float y1 = -1.;
E_Float y2 = -1.;
if (fabs(delta) < tol && fabs(delta0) < tol) 
{
  y1 = -bp/(3*ap);  
}
else if (fabs(delta) < tol)
{
  y1 = (9*ap*dp - bp*cp)/(2*delta0);
  y2 = (4*ap*bp*cp-9*ap*ap*dp-bp*bp*bp)/(ap*delta0);
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

  y1 = -1./(3*ap)*(bp + C1 + delta0/C1 );
  y2 = -1./(3*ap)*(bp + C2 + delta0/C2 );
}


E_Float c1 = 2*y1-p;
E_Float c2 = 2*y2-p;

E_Float z1 = -123456.;
if (c1 >= tol)
{
  E_Float p1 = -2*y1-p+2*q/(sqrt(c1));
  E_Float p2 = -2*y1-p-2*q/(sqrt(c1));
  if (p1 >= tol) { z1 = 0.5*(sqrt(c1)+sqrt(p1));}
  else if (p2 >= tol) { z1 = 0.5*(sqrt(c1)+sqrt(p2));}
}
if (c2 >= tol && z1 == -123456)
{
  E_Float p1 = -2*y2-p+2*q/(sqrt(c2));
  E_Float p2 = -2*y2-p-2*q/(sqrt(c2));
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

E_Float nutilde1 = (z1 + a/4.)*nu;
