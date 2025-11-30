E_Float Ro, Ux, Uy, Uz, ec;
E_Float coef1, coef2, coef3, coef4, coef5, coef6;

E_Float one_3 = 1./3.;
E_Float one_18= 1./18.;
E_Float one_36= 1./36.;

E_Int var_start=4;
if (nvars_loc == 36){ var_start=16;}

for (E_Int noind = pt_deb; noind < pt_fin; noind++)
{

  E_Int indR   = noind;

  Ro = vectOfRcvFields[0][indR]; //ro
  Ux = vectOfRcvFields[1][indR]; //vx
  Uy = vectOfRcvFields[2][indR]; //vy
  Uz = vectOfRcvFields[3][indR]; //vz

  ec   = Ux*Ux+ Uy*Uy +Uz*Uz;

  coef1 = Ux*Ux*Uy+Uy*Uz*Uz;
  coef2 = Ux*Uz*Uz+Ux*Uy*Uy;
  coef3 = Uy*Uy*Uz+Ux*Ux*Uz;
  coef4 = Ux*Ux*Uy-Uy*Uz*Uz;
  coef5 = Ux*Uz*Uz-Ux*Uy*Uy;
  coef6 = Uy*Uy*Uz-Ux*Ux*Uz;

  //p=1
  //psi = one_3*1.5*(corr_xx+corr_yy+corr_zz);
  E_Float psi = 0.;
  E_Float feq = Ro*one_3*(1.-1.5*ec);
  vectOfRcvFields[var_start+1][indR] = feq + 0.5*psi;

  //p=2
  //psi = -one_18*1.5*(2*corr_xx-corr_yy-corr_zz);
  E_Float vit      = Ux;
  feq      = Ro*one_18*(1.+vit*(3.+4.5*vit)-1.5*ec-9.*coef2);
  vectOfRcvFields[var_start+2][indR] = feq + 0.5*psi;

  //p=3
  feq      = feq - 6*Ro*one_18*(vit-3.*coef2);
  vectOfRcvFields[var_start+3][indR] = feq + 0.5*psi;

  //p=4
  // psi = -one_18*1.5*(2*corr_yy-corr_zz-corr_xx);
  vit      = Uy;
  feq      = Ro*one_18*(1.+vit*(3.+4.5*vit)-1.5*ec-9.*coef1);
  vectOfRcvFields[var_start+4][indR] = feq + 0.5*psi;

  //p=5
  feq      = feq - 6*Ro*one_18*(vit-3*coef1);
  vectOfRcvFields[var_start+5][indR] = feq + 0.5*psi;

  //p=6
  //psi = -one_18*1.5*(-corr_xx-corr_yy+2*corr_zz);
  vit      = Uz;
  feq      = Ro*one_18*(1.+vit*(3.+4.5*vit)-1.5*ec-9.*coef3);
  vectOfRcvFields[var_start+6][indR] = feq + 0.5*psi;

  //p=7
  feq      = feq - 6*Ro*one_18*(vit-3*coef1);
  vectOfRcvFields[var_start+7][indR] = feq + 0.5*psi;

  //p=8
  //psi_b = -one_36*1.5*(2*(corr_xx+corr_yy)-corr_zz);
  //psi   = psi_b - one_36*1.5*6*corr_z;
  vit           = Ux + Uy;
  E_Float herm  = -coef1-coef2-coef4+coef5;
  feq           = Ro*one_36*(1.+vit*(3.+4.5*vit)-1.5*ec -4.5*herm);
  vectOfRcvFields[var_start+8][indR] = feq + 0.5*psi;

  // p=9
  // psi_alt = psi_b + one_36*1.5*6*corr_z;
  E_Float  psi_alt = 0.;
  E_Float vit2     = -Ux + Uy;
  E_Float herm_alt = -coef1+coef2-coef4-coef5;
  E_Float feq_alt  = Ro*one_36*(1.+vit2*(3.+4.5*vit2)-1.5*ec-4.5*herm_alt);
  vectOfRcvFields[var_start+9][indR] = feq_alt + 0.5*psi_alt;

  //p=10
  feq_alt   = feq_alt-Ro*one_36*(6*vit2-9*herm_alt);
  vectOfRcvFields[var_start+10][indR] = feq_alt + 0.5*psi_alt;

  //p=11
  feq   = feq -Ro*one_36*(6*vit-9*herm);
  vectOfRcvFields[var_start+11][indR] = feq + 0.5*psi;


  //p=12
  //psi_b = -one_36*1.5*(2*(corr_xx+corr_zz)-corr_yy);
  //psi   = psi_b - one_36*1.5*6*corr_y;
  vit       = Ux + Uz;
  herm      = -coef2-coef3-coef5+coef6;
  feq       = Ro*one_36*(1.+ vit*(3.+4.5*vit) -1.5*ec-4.5*herm);
  vectOfRcvFields[var_start+12][indR] = feq + 0.5*psi;

  //p=13
  //psi_alt = psi_b + one_36*1.5*6*corr_y;
  vit2      = -Ux + Uz;
  herm_alt  = coef2-coef3+coef5+coef6;
  feq_alt   = Ro*one_36*(1.+vit2*(3.+4.5*vit2)-1.5*ec -4.5*herm_alt);
  vectOfRcvFields[var_start+13][indR] = feq_alt + 0.5*psi_alt;

  //p=14
  feq_alt   = feq_alt-Ro*one_36*(6*vit2-9*herm_alt);
  vectOfRcvFields[var_start+14][indR] = feq_alt + 0.5*psi_alt;

  //p=15
  feq   = feq -Ro*one_36*(6*vit-9*herm);
  vectOfRcvFields[var_start+15][indR] = feq + 0.5*psi;

  //p=16
  //psi_b = -one_36*1.5*(2*(corr_yy+corr_zz)-corr_xx);
  //psi   = psi_b - one_36*1.5*6*corr_x;
  vit       = Uy + Uz;
  herm      = -coef1-coef3+coef4-coef6;
  feq       = Ro*one_36*(1.+ vit*(3.+4.5*vit) -1.5*ec -4.5*herm);
  vectOfRcvFields[var_start+16][indR] = feq + 0.5*psi;

  // p=17
  //psi_alt = psi_b + one_36*1.5*6*corr_x;
  vit2      = -Uy + Uz;
  herm_alt  = coef1-coef3-coef4-coef6;
  feq_alt   = Ro*one_36*(1.+vit2*(3.+4.5*vit2)-1.5*ec -4.5*herm_alt);
  vectOfRcvFields[var_start+17][indR] = feq_alt + 0.5*psi_alt;

  //p=18
  feq_alt   = feq_alt-Ro*one_36*(6*vit2-9*herm_alt);
  vectOfRcvFields[var_start+18][indR] = feq_alt + 0.5*psi_alt;

  //p=19
  feq   = feq -Ro*one_36*(6*vit-9*herm);
  vectOfRcvFields[var_start+19][indR] = feq + 0.5*psi;
}

