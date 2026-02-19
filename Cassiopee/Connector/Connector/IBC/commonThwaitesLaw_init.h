# include "IBC/commonLaws1.h"

twall  = text + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext*uext);               // Temperature a la paroi
rowall = pext/twall*cvgaminv;                                               // Densite a la paroi
// muwall = coefSuth * sqrt(K_FUNC::E_abs(twall)*Tsinv) / (1.+Cs/twall);    // Viscosite a la paroi en utilisant temperature absolue pour reference
muwall = muext*sqrt(twall/text)*(1+Cs/text)/(1+Cs/twall);                   // Viscosite a la paroi en utilisant temperature exterieure (PI) pour reference (Benjamin's formula)

// angle d'attaque nul
utinf = 1.;
vtinf = 0.;
wtinf = 0.;

m_thwaites = acos(utinf*ut/uext + vtinf*vt/uext + wtinf*wt/uext);

if (yPW[noind+ideb] >= 0 &&  vt/uext < vtinf)
{
  m_thwaites = -m_thwaites;
}
else if (yPW[noind+ideb] <= 0 &&  vt/uext > vtinf)
{
  m_thwaites = -m_thwaites;
}

m_thwaites = m_thwaites/(M_PI - m_thwaites);

lambda_thwaites = 0.45*m_thwaites/(5.*m_thwaites+1);

if (lambda_thwaites >= 0.)
{
  t_thwaites = 0.22 + 1.57*lambda_thwaites - 1.8*lambda_thwaites*lambda_thwaites;
}
else
{
  t_thwaites = 0.22 + 1.402*lambda_thwaites - 0.018*lambda_thwaites/(0.107+lambda_thwaites);
}

if (m_thwaites > 0.)
{
  a_thwaites = lambda_thwaites/36.;
  b_thwaites = 2.*lambda_thwaites/3. - pow(b_thwaites,2);
  c_thwaites = 4.*lambda_thwaites;
  c_thwaites = pow(b_thwaites,2) - 4*a_thwaites*c_thwaites;

  lambda_thwaites = (-b_thwaites-sqrt(c_thwaites))/(2.*a_thwaites);
  utau_vec[noind]   = sqrt(sqrt(pow(uext,3)*muwall*(5*m_thwaites+1)/(xPW[noind+ideb]*rowall*0.45))*t_thwaites);
  utau0             = sqrt(sqrt(pow(uext,3)*muwall*(5*m_thwaites+1)/(xPW[noind+ideb]*rowall*0.45))*t_thwaites);
}
else
{
  lambda_thwaites = 0.;
  utau_vec[noind]   = sqrt(sqrt(pow(uext,3)*muwall/(xPW[noind+ideb]*rowall))*0.343);
  utau0             = sqrt(sqrt(pow(uext,3)*muwall/(xPW[noind+ideb]*rowall))*0.343);
}

press_vec[noind]    = pext;

ro_vec[noind]       = rowall;

aa_vec[noind]       = rowall*yext/muwall;
yplus_vec[noind]    = rowall*yibc/muwall;

uext_vec[noind]     = uext;
nutcible_vec[noind] = utau0;
sign_vec[noind]     = signibc/uext;
ucible_vec[noind]   = alphasbeta*un;
vcible_vec[noind]   = alphasbeta*vn;
wcible_vec[noind]   = alphasbeta*wn;
ut_vec[noind]       = ut;
vt_vec[noind]       = vt;
wt_vec[noind]       = wt;

mu_vec[noind]       = muwall;

alpha_vec[noind]    = alpha;
tcible_vec[noind]   = text;
