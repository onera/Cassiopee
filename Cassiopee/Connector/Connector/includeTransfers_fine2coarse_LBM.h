//======================================================
// Transfert fine grid to coarse grid LBM
//======================================================

E_Float taug_f = ipt_param_realR[NoD][LBM_TAUG];
E_Float taug_c = ipt_param_realR[NoR][LBM_TAUG];

E_Int flag_psi = 0;
E_Int nb_dist    = ipt_param_intR[NoR][NEQ_LBM];
E_Int coll_model = ipt_param_intR[NoR][LBM_COLL_MODEL];
if ( coll_model==4 ) flag_psi = 1;

E_Int shift_macro = 0;
E_Int shift_sij   = 5;
E_Int shift_psi   = 11;
E_Int shift_q     = 17;

E_Float scale_f2c = 2*taug_c/taug_f;

// Adimensionnement masse volumique
E_Float rho_ref = ipt_param_realR[NoD][ROINF];
E_Float adim_rho = 1. / rho_ref;

// Adimensionnement vitesse
E_Float c0 = 1. / sqrt(3.);
E_Float gam  = ipt_param_realR[NoD][GAMMA];
E_Float rgp  = ipt_param_realR[NoD][CVINF] * (gam - 1);
E_Float Tref = ipt_param_realR[NoD][TINF];
E_Float cson = sqrt(gam * rgp * Tref);
E_Float adim_vit = c0 / cson;

// Adimensionnement corr
E_Float adim_corr = adim_rho * pow(adim_vit, 3.);

if (nb_dist == 19)
{

  E_Float w1 = 1. / 3.;
  E_Float w2 = 1. / 18.;
  E_Float w3 = 1. / 36.;

  E_Float vit, vit2, herm, herm_alt;
  E_Float Ro, Ux, Uy, Uz, ec;
  E_Float coef1, coef2, coef3, coef4, coef5, coef6;
  E_Float eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10;
  E_Float eq11, eq12, eq13, eq14, eq15, eq16, eq17, eq18, eq19;

  E_Float corr_xx, corr_yy, corr_zz, corr_x, corr_y, corr_z, psi, psi_b, psi_alt;

  E_Float offeq1, offeq2, offeq3, offeq4, offeq5, offeq6, offeq7, offeq8, offeq9, offeq10;
  E_Float offeq11, offeq12, offeq13, offeq14, offeq15, offeq16, offeq17, offeq18, offeq19;

  for (E_Int noind = pt_deb; noind < pt_fin; noind++)
  {

    indR = rcvPts[noind];

    // Calcul de l'equilibre
    Ro = vectOfRcvFields[0][indR]*adim_rho;
    Ux = vectOfRcvFields[1][indR]*adim_vit;
    Uy = vectOfRcvFields[2][indR]*adim_vit;
    Uz = vectOfRcvFields[3][indR]*adim_vit;
    ec = Ux * Ux + Uy * Uy + Uz * Uz;

    coef1 = Ux * Ux * Uy + Uy * Uz * Uz;
    coef2 = Ux * Uz * Uz + Ux * Uy * Uy;
    coef3 = Uy * Uy * Uz + Ux * Ux * Uz;
    coef4 = Ux * Ux * Uy - Uy * Uz * Uz;
    coef5 = Ux * Uz * Uz - Ux * Uy * Uy;
    coef6 = Uy * Uy * Uz - Ux * Ux * Uz;

    corr_xx = vectOfRcvFields[shift_psi    ][indR]*adim_corr;
    corr_yy = vectOfRcvFields[shift_psi + 1][indR]*adim_corr;
    corr_zz = vectOfRcvFields[shift_psi + 2][indR]*adim_corr;
    corr_x  = vectOfRcvFields[shift_psi + 3][indR]*adim_corr;
    corr_y  = vectOfRcvFields[shift_psi + 4][indR]*adim_corr;
    corr_z  = vectOfRcvFields[shift_psi + 5][indR]*adim_corr;

    // Q 1
    eq1 = Ro * w1 * (1. - 1.5 * ec);
    psi = 0.5 * w1 * 1.5 * (corr_xx + corr_yy + corr_zz);
    offeq1 = vectOfRcvFields[shift_q][indR] - eq1 - flag_psi * psi;
    vectOfRcvFields[shift_q][indR] = eq1 + offeq1 * scale_f2c + 2. * flag_psi * psi;

    // Q 2
    vit = Ux;
    eq2 = Ro * w2 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 9. * coef2);
    psi = -0.5 * w2 * 1.5 * (2 * corr_xx - corr_yy - corr_zz);
    offeq2 = vectOfRcvFields[shift_q + 1][indR] - eq2 - flag_psi * psi;
    vectOfRcvFields[shift_q + 1][indR] = eq2 + offeq2 * scale_f2c + 2. * flag_psi * psi;

    // Q 3
    eq3 = eq2 - 6 * Ro * w2 * (vit - 3. * coef2);
    offeq3 = vectOfRcvFields[shift_q + 2][indR] - eq3 - flag_psi * psi;
    vectOfRcvFields[shift_q + 2][indR] = eq3 + offeq3 * scale_f2c + 2. * flag_psi * psi;

    // Q 4
    vit = Uy;
    eq4 = Ro * w2 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 9. * coef1);
    psi = -0.5 * w2 * 1.5 * (2 * corr_yy - corr_zz - corr_xx);
    offeq4 = vectOfRcvFields[shift_q + 3][indR] - eq4 - flag_psi * psi;
    vectOfRcvFields[shift_q + 3][indR] = eq4 + offeq4 * scale_f2c + 2. * flag_psi * psi;

    // Q 5
    eq5 = eq4 - 6 * Ro * w2 * (vit - 3 * coef1);
    offeq5 = vectOfRcvFields[shift_q + 4][indR] - eq5 - flag_psi * psi;
    vectOfRcvFields[shift_q + 4][indR] = eq5 + offeq5 * scale_f2c + 2. * flag_psi * psi;

    // Q 6
    vit = Uz;
    eq6 = Ro * w2 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 9. * coef3);
    psi = -0.5 * w2 * 1.5 * (-corr_xx - corr_yy + 2 * corr_zz);
    offeq6 = vectOfRcvFields[shift_q + 5][indR] - eq6 - flag_psi * psi;
    vectOfRcvFields[shift_q + 5][indR] = eq6 + offeq6 * scale_f2c + 2. * flag_psi * psi;

    // Q 7
    eq7 = eq6 - 6 * Ro * w2 * (vit - 3 * coef3);
    offeq7 = vectOfRcvFields[shift_q + 6][indR] - eq7 - flag_psi * psi;
    vectOfRcvFields[shift_q + 6][indR] = eq7 + offeq7 * scale_f2c + 2. * flag_psi * psi;

    // Q 8
    vit = Ux + Uy;
    herm = -coef1 - coef2 - coef4 + coef5;
    eq8 = Ro * w3 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm);
    psi_b = -0.5 * w3 * 1.5 * (2 * (corr_xx + corr_yy) - corr_zz);
    psi = psi_b - 0.5 * w3 * 1.5 * 6 * corr_z;
    offeq8 = vectOfRcvFields[shift_q + 7][indR] - eq8 - flag_psi * psi;
    vectOfRcvFields[shift_q + 7][indR] = eq8 + offeq8 * scale_f2c + 2. * flag_psi * psi;

    // Q 9
    vit2 = -Ux + Uy;
    herm_alt = -coef1 + coef2 - coef4 - coef5;
    eq9 = Ro * w3 * (1. + vit2 * (3. + 4.5 * vit2) - 1.5 * ec - 4.5 * herm_alt);
    psi_alt = psi_b + 0.5 * w3 * 1.5 * 6 * corr_z;
    offeq9 = vectOfRcvFields[shift_q + 8][indR] - eq9 - flag_psi * psi_alt;
    vectOfRcvFields[shift_q + 8][indR] = eq9 + offeq9 * scale_f2c + 2. * flag_psi * psi_alt;

    // Q 10
    eq10 = eq9 - Ro * w3 * (6 * vit2 - 9 * herm_alt);
    offeq10 = vectOfRcvFields[shift_q + 9][indR] - eq10 - flag_psi * psi_alt;
    vectOfRcvFields[shift_q + 9][indR] = eq10 + offeq10 * scale_f2c + 2. * flag_psi * psi_alt;

    // Q 11
    eq11 = eq8 - Ro * w3 * (6 * vit - 9 * herm);
    offeq11 = vectOfRcvFields[shift_q + 10][indR] - eq11 - flag_psi * psi;
    vectOfRcvFields[shift_q + 10][indR] = eq11 + offeq11 * scale_f2c + 2. * flag_psi * psi;

    // Q 12
    vit = Ux + Uz;
    herm = -coef2 - coef3 - coef5 + coef6;
    eq12 = Ro * w3 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm);
    psi_b = -0.5 * w3 * 1.5 * (2 * (corr_xx + corr_zz) - corr_yy);
    psi = psi_b - 0.5 * w3 * 1.5 * 6 * corr_y;
    offeq12 = vectOfRcvFields[shift_q + 11][indR] - eq12 - flag_psi * psi;
    vectOfRcvFields[shift_q + 11][indR] = eq12 + offeq12 * scale_f2c + 2. * flag_psi * psi;

    // Q 13
    vit2 = -Ux + Uz;
    herm_alt = coef2 - coef3 + coef5 + coef6;
    eq13 = Ro * w3 * (1. + vit2 * (3. + 4.5 * vit2) - 1.5 * ec - 4.5 * herm_alt);
    psi_alt = psi_b + 0.5 * w3 * 1.5 * 6 * corr_y;
    offeq13 = vectOfRcvFields[shift_q + 12][indR] - eq13 - flag_psi * psi_alt;
    vectOfRcvFields[shift_q + 12][indR] = eq13 + offeq13 * scale_f2c + 2. * flag_psi * psi_alt;

    // Q 14
    eq14 = eq13 - Ro * w3 * (6 * vit2 - 9 * herm_alt);
    offeq14 = vectOfRcvFields[shift_q + 13][indR] - eq14 - flag_psi * psi_alt;
    vectOfRcvFields[shift_q + 13][indR] = eq14 + offeq14 * scale_f2c + 2. * flag_psi * psi_alt;

    // Q 15
    eq15 = eq12 - Ro * w3 * (6 * vit - 9 * herm);
    offeq15 = vectOfRcvFields[shift_q + 14][indR] - eq15 - flag_psi * psi;
    vectOfRcvFields[shift_q + 14][indR] = eq15 + offeq15 * scale_f2c + 2. * flag_psi * psi;

    // Q 16
    vit = Uy + Uz;
    herm = -coef1 - coef3 + coef4 - coef6;
    eq16 = Ro * w3 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm);
    psi_b = -0.5 * w3 * 1.5 * (2 * (corr_yy + corr_zz) - corr_xx);
    psi = psi_b - 0.5 * w3 * 1.5 * 6 * corr_x;
    offeq16 = vectOfRcvFields[shift_q + 15][indR] - eq16 - flag_psi * psi;
    vectOfRcvFields[shift_q + 15][indR] = eq16 + offeq16 * scale_f2c + 2. * flag_psi * psi;

    // Q 17
    vit2 = -Uy + Uz;
    herm_alt = coef1 - coef3 - coef4 - coef6;
    eq17 = Ro * w3 * (1. + vit2 * (3. + 4.5 * vit2) - 1.5 * ec - 4.5 * herm_alt);
    psi_alt = psi_b + 0.5 * w3 * 1.5 * 6 * corr_x;
    offeq17 = vectOfRcvFields[shift_q + 16][indR] - eq17 - flag_psi * psi_alt;
    vectOfRcvFields[shift_q + 16][indR] = eq17 + offeq17 * scale_f2c + 2. * flag_psi * psi_alt;

    // Q 18
    eq18 = eq17 - Ro * w3 * (6 * vit2 - 9 * herm_alt);
    offeq18 = vectOfRcvFields[shift_q + 17][indR] - eq18 - flag_psi * psi_alt;
    vectOfRcvFields[shift_q + 17][indR] = eq18 + offeq18 * scale_f2c + 2. * flag_psi * psi_alt;

    // Q 19
    eq19 = eq16 - Ro * w3 * (6 * vit - 9 * herm);
    offeq19 = vectOfRcvFields[shift_q + 18][indR] - eq19 - flag_psi * psi;
    vectOfRcvFields[shift_q + 18][indR] = eq19 + offeq19 * scale_f2c + 2. * flag_psi * psi;

  }
  // Fin cas D3Q19
}
else if (nb_dist == 27)
{

    E_Float w1 = 8. / 27.;
    E_Float w2 = 2. / 27.;
    E_Float w3 = 1. / 54.;
    E_Float w4 = 1. / 216.;

    E_Float vit, vit2, vit3, vit4;
    E_Float herm, herm2, herm3, herm4, herm_alt, herm2_alt, herm3_alt, herm4_alt;
    E_Float herm_alt2, herm2_alt2, herm3_alt2, herm4_alt2, herm_alt3, herm2_alt3, herm3_alt3, herm4_alt3;
    E_Float Ro, Ux, Uy, Uz, ec;
    E_Float axxy, axxz, axyy, axzz, ayzz, ayyz, axyz;
    E_Float axxyy, axxzz, ayyzz, axyzz, axyyz, axxyz;
    E_Float axxyzz, axxyyz, axyyzz, axxyyzz;
    E_Float eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10;
    E_Float eq11, eq12, eq13, eq14, eq15, eq16, eq17, eq18, eq19;
    E_Float eq20, eq21, eq22, eq23, eq24, eq25, eq26, eq27;
    E_Float corr_xx, corr_yy, corr_zz, corr_x, corr_y, corr_z, psi, psi_b, psi_alt;

    E_Float offeq1, offeq2, offeq3, offeq4, offeq5, offeq6, offeq7, offeq8, offeq9, offeq10;
    E_Float offeq11, offeq12, offeq13, offeq14, offeq15, offeq16, offeq17, offeq18, offeq19;
    E_Float offeq20, offeq21, offeq22, offeq23, offeq24, offeq25, offeq26, offeq27;

    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {

      indR = rcvPts[noind];

      // Calcul de l'equilibre
      Ro = vectOfRcvFields[0][indR]*adim_rho;
      Ux = vectOfRcvFields[1][indR]*adim_vit;
      Uy = vectOfRcvFields[2][indR]*adim_vit;
      Uz = vectOfRcvFields[3][indR]*adim_vit;
      ec = Ux * Ux + Uy * Uy + Uz * Uz;

      axxy = Ux * Ux * Uy;
      axxz = Ux * Ux * Uz;
      axyy = Ux * Uy * Uy;
      axzz = Ux * Uz * Uz;
      ayzz = Uy * Uz * Uz;
      ayyz = Uy * Uy * Uz;
      axyz = Ux * Uy * Uz;

      axxyy = Ux * Ux * Uy * Uy;
      axxzz = Ux * Ux * Uz * Uz;
      ayyzz = Uy * Uy * Uz * Uz;
      axyzz = Ux * Uy * Uz * Uz;
      axyyz = Ux * Uy * Uy * Uz;
      axxyz = Ux * Ux * Uy * Uz;

      axxyzz = Ux * Ux * Uy * Uz * Uz;
      axxyyz = Ux * Ux * Uy * Uy * Uz;
      axyyzz = Ux * Uy * Uy * Uz * Uz;

      axxyyzz = Ux * Ux * Uy * Uy * Uz * Uz;

      corr_xx = vectOfRcvFields[shift_psi    ][indR]*adim_corr;
      corr_yy = vectOfRcvFields[shift_psi + 1][indR]*adim_corr;
      corr_zz = vectOfRcvFields[shift_psi + 2][indR]*adim_corr;
      corr_x  = vectOfRcvFields[shift_psi + 3][indR]*adim_corr;
      corr_y  = vectOfRcvFields[shift_psi + 4][indR]*adim_corr;
      corr_z  = vectOfRcvFields[shift_psi + 5][indR]*adim_corr;

      // Q 1
      herm2 = axxyy + axxzz + ayyzz;
      herm4 = axxyyzz;
      eq1 = Ro * w1 * (1. - 1.5 * ec + 2.25 * herm2 + 3.375 * herm4);
      psi = w1 * 1.5 * (-corr_xx - corr_yy - corr_zz);
      offeq1 = vectOfRcvFields[shift_q][indR] - eq1 - flag_psi * psi;
      vectOfRcvFields[shift_q][indR] = eq1 + offeq1 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 2
      vit = Ux;
      herm = axyy + axzz;
      herm2 = -2 * axxyy - 2 * axxzz + ayyzz;
      herm3 = 3 * axyyzz;
      herm4 = 4 * axxyyzz;
      eq2 = Ro * w2 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm + 2.25 * (herm2 + herm3) + 3.375 * herm4);
      psi = w2 * 1.5 * (2 * corr_xx - corr_yy - corr_zz);
      offeq2 = vectOfRcvFields[shift_q + 1][indR] - eq2 - flag_psi * psi;
      vectOfRcvFields[shift_q + 1][indR] = eq2 + offeq2 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 3
      eq3 = eq2 - 6 * Ro * w2 * (vit - 1.5 * herm + 0.75 * herm3);
      offeq3 = vectOfRcvFields[shift_q + 2][indR] - eq3 - flag_psi * psi;
      vectOfRcvFields[shift_q + 2][indR] = eq3 + offeq3 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 4
      vit = Uy;
      herm = axxy + ayzz;
      herm2 = -2 * axxyy + axxzz - 2 * ayyzz;
      herm3 = 3 * axxyzz;
      eq4 = Ro * w2 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm + 2.25 * (herm2 + herm3) + 3.375 * herm4);
      psi = w2 * 1.5 * (-corr_xx + 2 * corr_yy - corr_zz);
      offeq4 = vectOfRcvFields[shift_q + 3][indR] - eq4 - flag_psi * psi;
      vectOfRcvFields[shift_q + 3][indR] = eq4 + offeq4 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 5
      eq5 = eq4 - 6 * Ro * w2 * (vit - 1.5 * herm + 0.75 * herm3);
      offeq5 = vectOfRcvFields[shift_q + 4][indR] - eq5 - flag_psi * psi;
      vectOfRcvFields[shift_q + 4][indR] = eq5 + offeq5 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 6
      vit = Uz;
      herm = axxz + ayyz;
      herm2 = axxyy - 2 * axxzz - 2 * ayyzz;
      herm3 = 3 * axxyyz;
      eq6 = Ro * w2 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm + 2.25 * (herm2 + herm3) + 3.375 * herm4);
      psi = w2 * 1.5 * (-corr_xx - corr_yy + 2 * corr_zz);
      offeq6 = vectOfRcvFields[shift_q + 5][indR] - eq6 - flag_psi * psi;
      vectOfRcvFields[shift_q + 5][indR] = eq6 + offeq6 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 7
      eq7 = eq6 - 6 * Ro * w2 * (vit - 1.5 * herm + 0.75 * herm3);
      offeq7 = vectOfRcvFields[shift_q + 6][indR] - eq7 - flag_psi * psi;
      vectOfRcvFields[shift_q + 6][indR] = eq7 + offeq7 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 8
      vit = Ux + Uy;
      herm = -2 * axxy - 2 * axyy + axzz + ayzz;
      herm2 = 4 * axxyy - 2 * axxzz - 2 * ayyzz - 3 * axyzz;
      herm3 = -6 * axxyzz - 6 * axyyzz;
      herm4 = -2 * axxyyzz;
      eq8 = Ro * w3 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm + 2.25 * (herm2 + herm3) + 3.375 * herm4);
      psi = w3 * 1.5 * (2 * corr_xx + 2 * corr_yy - corr_zz);
      offeq8 = vectOfRcvFields[shift_q + 7][indR] - eq8 - flag_psi * psi;
      vectOfRcvFields[shift_q + 7][indR] = eq8 + offeq8 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 9
      vit2 = -Ux + Uy;
      herm_alt = -2 * axxy + 2 * axyy - axzz + ayzz;
      herm2_alt = 4 * axxyy - 2 * axxzz - 2 * ayyzz + 3 * axyzz;
      herm3_alt = -6 * axxyzz + 6 * axyyzz;
      eq9 = Ro * w3 * (1. + vit2 * (3. + 4.5 * vit2) - 1.5 * ec - 4.5 * herm_alt + 2.25 * (herm2_alt + herm3_alt) + 3.375 * herm4);
      offeq9 = vectOfRcvFields[shift_q + 8][indR] - eq9 - flag_psi * psi;
      vectOfRcvFields[shift_q + 8][indR] = eq9 + offeq9 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 10
      eq10 = eq9 - 6. * Ro * w3 * (vit2 - 1.5 * herm_alt + 0.75 * herm3_alt);
      offeq10 = vectOfRcvFields[shift_q + 9][indR] - eq10 - flag_psi * psi;
      vectOfRcvFields[shift_q + 9][indR] = eq10 + offeq10 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 11
      eq11 = eq8 - 6 * Ro * w3 * (vit - 1.5 * herm + 0.75 * herm3);
      offeq11 = vectOfRcvFields[shift_q + 10][indR] - eq11 - flag_psi * psi;
      vectOfRcvFields[shift_q + 10][indR] = eq11 + offeq11 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 12
      vit = Ux + Uz;
      herm = -2 * axxz + axyy - 2 * axzz + ayyz;
      herm2 = -2 * axxyy + 4 * axxzz - 2 * ayyzz - 3 * axyyz;
      herm3 = -6 * axxyyz - 6 * axyyzz;
      eq12 = Ro * w3 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm + 2.25 * (herm2 + herm3) + 3.375 * herm4);
      psi = w3 * 1.5 * (2 * corr_xx - corr_yy + 2 * corr_yy);
      offeq12 = vectOfRcvFields[shift_q + 11][indR] - eq12 - flag_psi * psi;
      vectOfRcvFields[shift_q + 11][indR] = eq12 + offeq12 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 13
      vit2 = -Ux + Uz;
      herm_alt = -2 * axxz - axyy + 2 * axzz + ayyz;
      herm2_alt = -2 * axxyy + 4 * axxzz - 2 * ayyzz + 3 * axyyz;
      herm3_alt = -6 * axxyyz + 6 * axyyzz;
      eq13 = Ro * w3 * (1. + vit2 * (3. + 4.5 * vit2) - 1.5 * ec - 4.5 * herm_alt + 2.25 * (herm2_alt + herm3_alt) + 3.375 * herm4);
      offeq13 = vectOfRcvFields[shift_q + 12][indR] - eq13 - flag_psi * psi;
      vectOfRcvFields[shift_q + 12][indR] = eq13 + offeq13 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 14
      eq14 = eq13 - 6 * Ro * w3 * (vit2 - 1.5 * herm_alt + 0.75 * herm3_alt);
      offeq14 = vectOfRcvFields[shift_q + 13][indR] - eq14 - flag_psi * psi;
      vectOfRcvFields[shift_q + 13][indR] = eq14 + offeq14 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 15
      eq15 = eq12 - 6 * Ro * w3 * (vit - 1.5 * herm + 0.75 * herm3);
      offeq15 = vectOfRcvFields[shift_q + 14][indR] - eq15 - flag_psi * psi;
      vectOfRcvFields[shift_q + 14][indR] = eq15 + offeq15 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 16
      vit = Uy + Uz;
      herm = axxy + axxz - 2 * ayzz - 2 * ayyz;
      herm2 = -2 * axxyy - 2 * axxzz + 4 * ayyzz - 3 * axxyz;
      herm3 = -6 * axxyzz - 6 * axxyyz;
      eq16 = Ro * w3 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm + 2.25 * (herm2 + herm3) + 3.375 * herm4);
      psi = w3 * 1.5 * (-corr_xx + 2 * corr_yy + 2 * corr_zz);
      offeq16 = vectOfRcvFields[shift_q + 15][indR] - eq16 - flag_psi * psi;
      vectOfRcvFields[shift_q + 15][indR] = eq16 + offeq16 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 17
      vit2 = -Uy + Uz;
      herm_alt = -axxy + axxz + 2 * ayzz - 2 * ayyz;
      herm2_alt = -2 * axxyy - 2 * axxzz + 4 * ayyzz + 3 * axxyz;
      herm3_alt = 6 * axxyzz - 6 * axxyyz;
      eq17 = Ro * w3 * (1. + vit2 * (3. + 4.5 * vit2) - 1.5 * ec - 4.5 * herm_alt + 2.25 * (herm2_alt + herm3_alt) + 3.375 * herm4);
      offeq17 = vectOfRcvFields[shift_q + 16][indR] - eq17 - flag_psi * psi;
      vectOfRcvFields[shift_q + 16][indR] = eq17 + offeq17 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 18
      eq18 = eq17 - 6 * Ro * w3 * (vit2 - 1.5 * herm_alt + 0.75 * herm3_alt);
      offeq18 = vectOfRcvFields[shift_q + 17][indR] - eq18 - flag_psi * psi;
      vectOfRcvFields[shift_q + 17][indR] = eq18 + offeq18 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 19
      eq19 = eq16 - 6 * Ro * w3 * (vit - 1.5 * herm + 0.75 * herm3);
      offeq19 = vectOfRcvFields[shift_q + 18][indR] - eq19 - flag_psi * psi;
      vectOfRcvFields[shift_q + 18][indR] = eq19 + offeq19 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 20
      vit = Ux + Uy + Uz;
      herm = -2 * (axxy + axxz + axyy + axzz + ayzz + ayyz + 3 * axyz);
      herm2 = 4 * (axxyy + axxzz + ayyzz) + 6 * (axyzz + axyyz + axxyz);
      herm3 = 12 * (axxyzz + axxyyz + axyyzz);
      herm4 = 10 * axxyyzz;
      eq20 = Ro * w4 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm + 2.25 * (herm2 + herm3) + 3.375 * herm4);
      psi = w4 * 1.5 * (2 * corr_xx + 2 * corr_yy + 2 * corr_zz);
      offeq20 = vectOfRcvFields[shift_q + 19][indR] - eq20 - flag_psi * psi;
      vectOfRcvFields[shift_q + 19][indR] = eq20 + offeq20 * scale_f2c + 0.5 * flag_psi * psi;

      // Q21
      vit2 = -Ux + Uy + Uz;
      herm_alt = -2 * (axxy + axxz - axyy - axzz + ayzz + ayyz - 3 * axyz);
      herm2_alt = 4 * (axxyy + axxzz + ayyzz) + 6 * (axxyz - axyzz - axyyz);
      herm3_alt = 12 * (axxyzz + axxyyz - axyyzz);
      eq21 = Ro * w4 * (1. + vit2 * (3. + 4.5 * vit2) - 1.5 * ec - 4.5 * herm_alt + 2.25 * (herm2_alt + herm3_alt) + 3.375 * herm4);
      offeq21 = vectOfRcvFields[shift_q + 20][indR] - eq21 - flag_psi * psi;
      vectOfRcvFields[shift_q + 20][indR] = eq21 + offeq21 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 22
      vit3 = Ux - Uy + Uz;
      herm_alt2 = -2 * (-axxy + axxz + axyy + axzz - ayzz + ayyz - 3 * axyz);
      herm2_alt2 = 4 * (axxyy + axxzz + ayyzz) + 6 * (-axyzz + axyyz - axxyz);
      herm3_alt2 = 12 * (-axxyzz + axxyyz + axyyzz);
      eq22 = Ro * w4 * (1. + vit3 * (3. + 4.5 * vit3) - 1.5 * ec - 4.5 * herm_alt2 + 2.25 * (herm2_alt2 + herm3_alt2) + 3.375 * herm4);
      offeq22 = vectOfRcvFields[shift_q + 21][indR] - eq22 - flag_psi * psi;
      vectOfRcvFields[shift_q + 21][indR] = eq22 + offeq22 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 23
      vit4 = -Ux - Uy + Uz;
      herm_alt3 = -2 * (-axxy + axxz - axyy - axzz - ayzz + ayyz + 3 * axyz);
      herm2_alt3 = 4 * (axxyy + axxzz + ayyzz) + 6 * (axyzz - axyyz - axxyz);
      herm3_alt3 = 12 * (-axxyzz + axxyyz - axyyzz);
      eq23 = Ro * w4 * (1. + vit4 * (3. + 4.5 * vit4) - 1.5 * ec - 4.5 * herm_alt3 + 2.25 * (herm2_alt3 + herm3_alt3) + 3.375 * herm4);
      offeq23 = vectOfRcvFields[shift_q + 22][indR] - eq23 - flag_psi * psi;
      vectOfRcvFields[shift_q + 22][indR] = eq23 + offeq23 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 24
      eq24 = eq23 - 6 * Ro * w4 * (vit4 - 1.5 * herm_alt3 + 0.75 * herm3_alt3);
      offeq24 = vectOfRcvFields[shift_q + 23][indR] - eq24 - flag_psi * psi;
      vectOfRcvFields[shift_q + 23][indR] = eq24 + offeq24 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 25
      eq25 = eq22 - 6 * Ro * w4 * (vit3 - 1.5 * herm_alt2 + 0.75 * herm3_alt2);
      offeq25 = vectOfRcvFields[shift_q + 24][indR] - eq25 - flag_psi * psi;
      vectOfRcvFields[shift_q + 24][indR] = eq25 + offeq25 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 26
      eq26 = eq21 - 6 * Ro * w4 * (vit2 - 1.5 * herm_alt + 0.75 * herm3_alt);
      offeq26 = vectOfRcvFields[shift_q + 25][indR] - eq26 - flag_psi * psi;
      vectOfRcvFields[shift_q + 25][indR] = eq26 + offeq26 * scale_f2c + 0.5 * flag_psi * psi;

      // Q 27
      eq27 = eq20 - 6 * Ro * w4 * (vit - 1.5 * herm + 0.75 * herm3);
      offeq27 = vectOfRcvFields[shift_q + 26][indR] - eq27 - flag_psi * psi;
      vectOfRcvFields[shift_q + 26][indR] = eq27 + offeq27 * scale_f2c + 0.5 * flag_psi * psi;

    }
    // Fin cas D3Q27
}
