
// Coef d'adimensionnement masse volumique
E_Float rho_ref = ipt_param_realR[NoD][ROINF];
E_Float adim_rho = 1. / rho_ref;

// Coef d'adimensionnement vitesse
E_Float c0 = 1. / sqrt(3.);
E_Float gam = ipt_param_realR[NoD][GAMMA];
E_Float rgp = ipt_param_realR[NoD][CVINF] * (gam - 1);
E_Float Tref = ipt_param_realR[NoD][TINF];
E_Float cson = sqrt(gam * rgp * Tref);
E_Float adim_vit = c0 / cson;

// Recupere le nombre de vitesses dans le lattice
E_Int nb_dist = ipt_param_intR[NoR][NEQ_LBM];

// ==================================================
//  D3Q19
// ==================================================
if (nb_dist == 19)
{

  // Lattice weights
  E_Float w1 = 1. / 3.;
  E_Float w2 = 1. / 18.;
  E_Float w3 = 1. / 36.;

  // Shift pour acceder aux Qs dans vectOfRcvFields
  E_Int shift_q = 4;
  if (nvars_loc == 36) shift_q = 16;

  for (E_Int noind = pt_deb; noind < pt_fin; noind++)
  {

    // Indice du point receveur
    indR = rcvPts[noind];

    // On recupere Ro, U qu'on adimensionne
    E_Float Ro = vectOfRcvFields[0][indR] * adim_rho; // ro
    E_Float Ux = vectOfRcvFields[1][indR] * adim_vit; // vx
    E_Float Uy = vectOfRcvFields[2][indR] * adim_vit; // vy
    E_Float Uz = vectOfRcvFields[3][indR] * adim_vit; // vz
    E_Float ec = Ux * Ux + Uy * Uy + Uz * Uz;

    // Coefficients pour feq
    E_Float coef1 = Ux * Ux * Uy + Uy * Uz * Uz;
    E_Float coef2 = Ux * Uz * Uz + Ux * Uy * Uy;
    E_Float coef3 = Uy * Uy * Uz + Ux * Ux * Uz;
    E_Float coef4 = Ux * Ux * Uy - Uy * Uz * Uz;
    E_Float coef5 = Ux * Uz * Uz - Ux * Uy * Uy;
    E_Float coef6 = Uy * Uy * Uz - Ux * Ux * Uz;

    // Q 1
    E_Float feq = Ro * w1 * (1. - 1.5 * ec);
    vectOfRcvFields[shift_q + 1][indR] = feq;

    // Q 2 
    E_Float vit = Ux;
    feq = Ro * w2 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 9. * coef2);
    vectOfRcvFields[shift_q + 2][indR] = feq;

    // Q 3
    feq = feq - 6 * Ro * w2 * (vit - 3. * coef2);
    vectOfRcvFields[shift_q + 3][indR] = feq;

    // Q 4
    vit = Uy;
    feq = Ro * w2 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 9. * coef1);
    vectOfRcvFields[shift_q + 4][indR] = feq;

    // Q 5
    feq = feq - 6 * Ro * w2 * (vit - 3 * coef1);
    vectOfRcvFields[shift_q + 5][indR] = feq;

    // Q 6
    vit = Uz;
    feq = Ro * w2 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 9. * coef3);
    vectOfRcvFields[shift_q + 6][indR] = feq;

    // Q 7
    feq = feq - 6 * Ro * w2 * (vit - 3 * coef1);
    vectOfRcvFields[shift_q + 7][indR] = feq;

    // Q 8
    vit = Ux + Uy;
    E_Float herm = -coef1 - coef2 - coef4 + coef5;
    feq = Ro * w3 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm);
    vectOfRcvFields[shift_q + 8][indR] = feq;

    // Q 9
    E_Float psi_alt = 0.;
    E_Float vit2 = -Ux + Uy;
    E_Float herm_alt = -coef1 + coef2 - coef4 - coef5;
    E_Float feq_alt = Ro * w3 * (1. + vit2 * (3. + 4.5 * vit2) - 1.5 * ec - 4.5 * herm_alt);
    vectOfRcvFields[shift_q + 9][indR] = feq_alt;

    // Q 10
    feq_alt = feq_alt - Ro * w3 * (6 * vit2 - 9 * herm_alt);
    vectOfRcvFields[shift_q + 10][indR] = feq_alt;

    // Q 11
    feq = feq - Ro * w3 * (6 * vit - 9 * herm);
    vectOfRcvFields[shift_q + 11][indR] = feq;

    // Q 12
    vit = Ux + Uz;
    herm = -coef2 - coef3 - coef5 + coef6;
    feq = Ro * w3 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm);
    vectOfRcvFields[shift_q + 12][indR] = feq;

    // Q 13
    vit2 = -Ux + Uz;
    herm_alt = coef2 - coef3 + coef5 + coef6;
    feq_alt = Ro * w3 * (1. + vit2 * (3. + 4.5 * vit2) - 1.5 * ec - 4.5 * herm_alt);
    vectOfRcvFields[shift_q + 13][indR] = feq_alt;

    // Q 14
    feq_alt = feq_alt - Ro * w3 * (6 * vit2 - 9 * herm_alt);
    vectOfRcvFields[shift_q + 14][indR] = feq_alt;

    // Q 15
    feq = feq - Ro * w3 * (6 * vit - 9 * herm);
    vectOfRcvFields[shift_q + 15][indR] = feq;

    // Q 16
    vit = Uy + Uz;
    herm = -coef1 - coef3 + coef4 - coef6;
    feq = Ro * w3 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm);
    vectOfRcvFields[shift_q + 16][indR] = feq;

    // Q 17
    vit2 = -Uy + Uz;
    herm_alt = coef1 - coef3 - coef4 - coef6;
    feq_alt = Ro * w3 * (1. + vit2 * (3. + 4.5 * vit2) - 1.5 * ec - 4.5 * herm_alt);
    vectOfRcvFields[shift_q + 17][indR] = feq_alt;

    // Q 18
    feq_alt = feq_alt - Ro * w3 * (6 * vit2 - 9 * herm_alt);
    vectOfRcvFields[shift_q + 18][indR] = feq_alt;

    // Q 19
    feq = feq - Ro * w3 * (6 * vit - 9 * herm);
    vectOfRcvFields[shift_q + 19][indR] = feq;
  }
// END CASE D3Q19
}
// ==================================================
//  D3Q27
// ==================================================
else if (nb_dist == 27)
{

  // Lattice weights
  E_Float w1 = 8. / 27.;
  E_Float w2 = 2. / 27.;
  E_Float w3 = 1. / 54.;
  E_Float w4 = 1. / 216.;

  // Shift pour acceder aux Qs dans vectOfRcvFields
  E_Int shift_q = 4;
  if (nvars_loc == 44) shift_q = 16;

  for (E_Int noind = pt_deb; noind < pt_fin; noind++)
  {

    // Indice du point receveur
    indR = rcvPts[noind];

    // On recupere Ro, U qu'on adimensionne
    E_Float Ro = vectOfRcvFields[0][indR] * adim_rho;
    E_Float Ux = vectOfRcvFields[1][indR] * adim_vit;
    E_Float Uy = vectOfRcvFields[2][indR] * adim_vit;
    E_Float Uz = vectOfRcvFields[3][indR] * adim_vit;
    E_Float ec = Ux * Ux + Uy * Uy + Uz * Uz;

    // Coefficients pour feq
    E_Float axxy = Ux * Ux * Uy;
    E_Float axxz = Ux * Ux * Uz;
    E_Float axyy = Ux * Uy * Uy;
    E_Float axzz = Ux * Uz * Uz;
    E_Float ayzz = Uy * Uz * Uz;
    E_Float ayyz = Uy * Uy * Uz;
    E_Float axyz = Ux * Uy * Uz;

    E_Float axxyy = Ux * Ux * Uy * Uy;
    E_Float axxzz = Ux * Ux * Uz * Uz;
    E_Float ayyzz = Uy * Uy * Uz * Uz;
    E_Float axyzz = Ux * Uy * Uz * Uz;
    E_Float axyyz = Ux * Uy * Uy * Uz;
    E_Float axxyz = Ux * Ux * Uy * Uz;

    E_Float axxyzz = Ux * Ux * Uy * Uz * Uz;
    E_Float axxyyz = Ux * Ux * Uy * Uy * Uz;
    E_Float axyyzz = Ux * Uy * Uy * Uz * Uz;
    E_Float axxyyzz = Ux * Ux * Uy * Uy * Uz * Uz;

    // Q 1
    E_Float herm2 = axxyy + axxzz + ayyzz;
    E_Float herm4 = axxyyzz;
    E_Float eq1 = Ro * w1 * (1. - 1.5 * ec + 2.25 * herm2 + 3.375 * herm4);
    vectOfRcvFields[shift_q + 1][indR] = eq1;

    // Q 2
    E_Float vit = Ux;
    E_Float herm = axyy + axzz;
    herm2 = -2 * axxyy - 2 * axxzz + ayyzz;
    E_Float herm3 = 3 * axyyzz;
    herm4 = 4 * axxyyzz;
    E_Float eq2 = Ro * w2 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm + 2.25 * (herm2 + herm3) + 3.375 * herm4);
    vectOfRcvFields[shift_q + 2][indR] = eq2;

    // Q 3
    E_Float eq3 = eq2 - 6 * Ro * w2 * (vit - 1.5 * herm + 0.75 * herm3);
    vectOfRcvFields[shift_q + 3][indR] = eq3;

    // Q 4
    vit = Uy;
    herm = axxy + ayzz;
    herm2 = -2 * axxyy + axxzz - 2 * ayyzz;
    herm3 = 3 * axxyzz;
    E_Float eq4 = Ro * w2 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm + 2.25 * (herm2 + herm3) + 3.375 * herm4);
    vectOfRcvFields[shift_q + 4][indR] = eq4;

    // Q 5
    E_Float eq5 = eq4 - 6 * Ro * w2 * (vit - 1.5 * herm + 0.75 * herm3);
    vectOfRcvFields[shift_q + 5][indR] = eq5;

    // Q 6
    vit = Uz;
    herm = axxz + ayyz;
    herm2 = axxyy - 2 * axxzz - 2 * ayyzz;
    herm3 = 3 * axxyyz;
    E_Float eq6 = Ro * w2 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm + 2.25 * (herm2 + herm3) + 3.375 * herm4);
    vectOfRcvFields[shift_q + 6][indR] = eq6;

    // Q 7
    E_Float eq7 = eq6 - 6 * Ro * w2 * (vit - 1.5 * herm + 0.75 * herm3);
    vectOfRcvFields[shift_q + 7][indR] = eq7;

    // Q 8
    vit = Ux + Uy;
    herm = -2 * axxy - 2 * axyy + axzz + ayzz;
    herm2 = 4 * axxyy - 2 * axxzz - 2 * ayyzz - 3 * axyzz;
    herm3 = -6 * axxyzz - 6 * axyyzz;
    herm4 = -2 * axxyyzz;
    E_Float eq8 = Ro * w3 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm + 2.25 * (herm2 + herm3) + 3.375 * herm4);
    vectOfRcvFields[shift_q + 8][indR] = eq8;

    // Q 9
    E_Float vit2 = -Ux + Uy;
    E_Float herm_alt = -2 * axxy + 2 * axyy - axzz + ayzz;
    E_Float herm2_alt = 4 * axxyy - 2 * axxzz - 2 * ayyzz + 3 * axyzz;
    E_Float herm3_alt = -6 * axxyzz + 6 * axyyzz;
    E_Float eq9 = Ro * w3 * (1. + vit2 * (3. + 4.5 * vit2) - 1.5 * ec - 4.5 * herm_alt + 2.25 * (herm2_alt + herm3_alt) + 3.375 * herm4);
    vectOfRcvFields[shift_q + 9][indR] = eq9;

    // Q 10
    E_Float eq10 = eq9 - 6. * Ro * w3 * (vit2 - 1.5 * herm_alt + 0.75 * herm3_alt);
    vectOfRcvFields[shift_q + 10][indR] = eq10;

    // Q 11
    E_Float eq11 = eq8 - 6 * Ro * w3 * (vit - 1.5 * herm + 0.75 * herm3);
    vectOfRcvFields[shift_q + 11][indR] = eq11;

    // Q 12
    vit = Ux + Uz;
    herm = -2 * axxz + axyy - 2 * axzz + ayyz;
    herm2 = -2 * axxyy + 4 * axxzz - 2 * ayyzz - 3 * axyyz;
    herm3 = -6 * axxyyz - 6 * axyyzz;
    E_Float eq12 = Ro * w3 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm + 2.25 * (herm2 + herm3) + 3.375 * herm4);
    vectOfRcvFields[shift_q + 12][indR] = eq12;

    // Q 13
    vit2 = -Ux + Uz;
    herm_alt = -2 * axxz - axyy + 2 * axzz + ayyz;
    herm2_alt = -2 * axxyy + 4 * axxzz - 2 * ayyzz + 3 * axyyz;
    herm3_alt = -6 * axxyyz + 6 * axyyzz;
    E_Float eq13 = Ro * w3 * (1. + vit2 * (3. + 4.5 * vit2) - 1.5 * ec - 4.5 * herm_alt + 2.25 * (herm2_alt + herm3_alt) + 3.375 * herm4);
    vectOfRcvFields[shift_q + 13][indR] = eq13;

    // Q 14
    E_Float eq14 = eq13 - 6 * Ro * w3 * (vit2 - 1.5 * herm_alt + 0.75 * herm3_alt);
    vectOfRcvFields[shift_q + 14][indR] = eq14;

    // Q 15
    E_Float eq15 = eq12 - 6 * Ro * w3 * (vit - 1.5 * herm + 0.75 * herm3);
    vectOfRcvFields[shift_q + 15][indR] = eq15;

    // Q 16
    vit = Uy + Uz;
    herm = axxy + axxz - 2 * ayzz - 2 * ayyz;
    herm2 = -2 * axxyy - 2 * axxzz + 4 * ayyzz - 3 * axxyz;
    herm3 = -6 * axxyzz - 6 * axxyyz;
    E_Float eq16 = Ro * w3 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm + 2.25 * (herm2 + herm3) + 3.375 * herm4);
    vectOfRcvFields[shift_q + 16][indR] = eq16;

    // Q 17
    vit2 = -Uy + Uz;
    herm_alt = -axxy + axxz + 2 * ayzz - 2 * ayyz;
    herm2_alt = -2 * axxyy - 2 * axxzz + 4 * ayyzz + 3 * axxyz;
    herm3_alt = 6 * axxyzz - 6 * axxyyz;
    E_Float eq17 = Ro * w3 * (1. + vit2 * (3. + 4.5 * vit2) - 1.5 * ec - 4.5 * herm_alt + 2.25 * (herm2_alt + herm3_alt) + 3.375 * herm4);
    vectOfRcvFields[shift_q + 17][indR] = eq17;

    // Q 18
    E_Float eq18 = eq17 - 6 * Ro * w3 * (vit2 - 1.5 * herm_alt + 0.75 * herm3_alt);
    vectOfRcvFields[shift_q + 18][indR] = eq18;

    // Q 19
    E_Float eq19 = eq16 - 6 * Ro * w3 * (vit - 1.5 * herm + 0.75 * herm3);
    vectOfRcvFields[shift_q + 19][indR] = eq19;

    // Q 20
    vit = Ux + Uy + Uz;
    herm = -2 * (axxy + axxz + axyy + axzz + ayzz + ayyz + 3 * axyz);
    herm2 = 4 * (axxyy + axxzz + ayyzz) + 6 * (axyzz + axyyz + axxyz);
    herm3 = 12 * (axxyzz + axxyyz + axyyzz);
    herm4 = 10 * axxyyzz;
    E_Float eq20 = Ro * w4 * (1. + vit * (3. + 4.5 * vit) - 1.5 * ec - 4.5 * herm + 2.25 * (herm2 + herm3) + 3.375 * herm4);
    vectOfRcvFields[shift_q + 20][indR] = eq20;

    // Q21
    vit2 = -Ux + Uy + Uz;
    herm_alt = -2 * (axxy + axxz - axyy - axzz + ayzz + ayyz - 3 * axyz);
    herm2_alt = 4 * (axxyy + axxzz + ayyzz) + 6 * (axxyz - axyzz - axyyz);
    herm3_alt = 12 * (axxyzz + axxyyz - axyyzz);
    E_Float eq21 = Ro * w4 * (1. + vit2 * (3. + 4.5 * vit2) - 1.5 * ec - 4.5 * herm_alt + 2.25 * (herm2_alt + herm3_alt) + 3.375 * herm4);
    vectOfRcvFields[shift_q + 21][indR] = eq21;

    // Q 22
    E_Float vit3 = Ux - Uy + Uz;
    E_Float herm_alt2 = -2 * (-axxy + axxz + axyy + axzz - ayzz + ayyz - 3 * axyz);
    E_Float herm2_alt2 = 4 * (axxyy + axxzz + ayyzz) + 6 * (-axyzz + axyyz - axxyz);
    E_Float herm3_alt2 = 12 * (-axxyzz + axxyyz + axyyzz);
    E_Float eq22 = Ro * w4 * (1. + vit3 * (3. + 4.5 * vit3) - 1.5 * ec - 4.5 * herm_alt2 + 2.25 * (herm2_alt2 + herm3_alt2) + 3.375 * herm4);
    vectOfRcvFields[shift_q + 22][indR] = eq22;

    // Q 23
    E_Float vit4 = -Ux - Uy + Uz;
    E_Float herm_alt3 = -2 * (-axxy + axxz - axyy - axzz - ayzz + ayyz + 3 * axyz);
    E_Float herm2_alt3 = 4 * (axxyy + axxzz + ayyzz) + 6 * (axyzz - axyyz - axxyz);
    E_Float herm3_alt3 = 12 * (-axxyzz + axxyyz - axyyzz);
    E_Float eq23 = Ro * w4 * (1. + vit4 * (3. + 4.5 * vit4) - 1.5 * ec - 4.5 * herm_alt3 + 2.25 * (herm2_alt3 + herm3_alt3) + 3.375 * herm4);
    vectOfRcvFields[shift_q + 23][indR] = eq23;

    // Q 24
    E_Float eq24 = eq23 - 6 * Ro * w4 * (vit4 - 1.5 * herm_alt3 + 0.75 * herm3_alt3);
    vectOfRcvFields[shift_q + 24][indR] = eq24;

    // Q 25
    E_Float eq25 = eq22 - 6 * Ro * w4 * (vit3 - 1.5 * herm_alt2 + 0.75 * herm3_alt2);
    vectOfRcvFields[shift_q + 25][indR] = eq25;

    // Q 26
    E_Float eq26 = eq21 - 6 * Ro * w4 * (vit2 - 1.5 * herm_alt + 0.75 * herm3_alt);
    vectOfRcvFields[shift_q + 26][indR] = eq26;

    // Q 27
    E_Float eq27 = eq20 - 6 * Ro * w4 * (vit - 1.5 * herm + 0.75 * herm3);
    vectOfRcvFields[shift_q + 27][indR] = eq27;
  }
// END CASE D3Q19
}
