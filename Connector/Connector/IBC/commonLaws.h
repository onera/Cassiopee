//==============================================================================
// loi log: estimation de utau
//==============================================================================
E_Float K_CONNECTOR::logf(E_Float x, E_Float a, E_Float b, E_Float c, E_Float d) 
{
  return a*log(b*x)-d/x+c;
}
// derivee
E_Float K_CONNECTOR::logfprime(E_Float x, E_Float a, E_Float d) 
{
  return a/x+d/(x*x);
}
//==============================================================================
// loi de musker : estimation de utau
//==============================================================================
E_Float K_CONNECTOR::musker(E_Float x, E_Float a, E_Float b)
{
  E_Float ax = a*x;
  E_Float ax2 = ax*ax;
  E_Float l1 = pow(ax+10.6,9.6);
  E_Float l2 = ax2-8.15*ax + 86.;
  E_Float l2s = l2*l2;
  return 5.424*atan((2.*ax-8.15)/16.7) + log10(l1/l2s) - b/x - 3.52;
}

//==============================================================================
E_Float K_CONNECTOR::muskerprime(E_Float x, E_Float a, E_Float b)
{
  E_Float t = pow((a*x + 10.6),9.6);
  E_Float tp = (a*9.6*pow((a*x + 10.6),8.6) * ((a*x)*(a*x) - 8.15*a*x + 86.)*((a*x)*(a*x) - 8.15*a*x + 86.)
                - t*(2.*(a*a)*x - 8.15*a)*2.*((a*x)*(a*x) - 8.15*a*x + 86.))/( ((a*x)*(a*x) - 8.15*a*x + 86.)*((a*x)*(a*x) - 8.15*a*x + 86.) );
  return 5.424*2./16.7*a/(1. + ((2.*a*x - 8.15)/16.7)*((2.*a*x - 8.15)/16.7) )+ tp/(t*log(10.)) + b/(x*x);
}
//==============================================================================
// fonctions pour estimer nutilde par Newton
//==============================================================================
E_Float K_CONNECTOR::fnutilde(E_Float nutilde, E_Float nut, E_Float rho, E_Float xmu)
{
  E_Float cv1cube = 7.1*7.1*7.1;   
  E_Float chi = (nutilde/xmu)*rho; 
  E_Float chicube = chi*chi*chi;
  E_Float fv1 = chicube/(chicube+cv1cube);
  return fv1*nutilde-nut;
}

E_Float K_CONNECTOR::fnutildeprime(E_Float nutilde, E_Float nut, E_Float rho, 
                                   E_Float xmu)
{
  E_Float cv1 = 7.1; 
  E_Float nutilde3 = nutilde*nutilde*nutilde;
  E_Float f = nutilde*nutilde3;//nutilde^4
  E_Float fnu = cv1*xmu/rho;
  E_Float g = nutilde3 + fnu*fnu*fnu;
  E_Float g2 = g*g;
  E_Float fp = 4*nutilde3;                              
  E_Float gp = 3*nutilde*nutilde;
  return (fp*g-f*gp)/g2;
}
