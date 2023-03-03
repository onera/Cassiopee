
E_Float dt_adim = ipt_param_realD[ NoD ][ DTC   ];
E_Float gam     = ipt_param_realD[ NoD ][ GAMMA ];
E_Float Rgp     = ipt_param_realD[ NoD ][ CVINF ]*(gam-1);
E_Float t_inf   = ipt_param_realD[ NoD ][ TINF  ];
E_Float ro_inf  = ipt_param_realD[ NoD ][ ROINF ];
E_Float c0_ref  = sqrt(gam*Rgp*t_inf);
E_Float u_scale = 1./(sqrt(3.)*c0_ref);

E_Float v1, v2, v3, v4;
E_Float g1, g2, g3, g4, g5, g6;

for (E_Int noind = pt_deb; noind < pt_fin; noind++)
{

  noind   = rcvPts[noind];

  v1 = vectOfRcvFields[0][noind]; //ro
  v2 = vectOfRcvFields[1][noind]; //vx
  v3 = vectOfRcvFields[2][noind]; //vy
  v4 = vectOfRcvFields[3][noind]; //vz
  // La temperature est modifiee par une BC
  g1 = vectOfRcvFields[5][noind]; //Sxx
  g2 = vectOfRcvFields[6][noind]; //Sxy
  g3 = vectOfRcvFields[7][noind]; //Sxz
  g4 = vectOfRcvFields[8][noind]; //Syy
  g5 = vectOfRcvFields[9][noind]; //Syz
  g6 = vectOfRcvFields[10][noind]; //Szz

  vectOfRcvFields[0][noind] = v1/ro_inf;
  vectOfRcvFields[1][noind] = v2*u_scale;
  vectOfRcvFields[2][noind] = v3*u_scale;
  vectOfRcvFields[3][noind] = 0.;//v4*u_scale;

  vectOfRcvFields[5][noind]  = g1*dt_adim;
  vectOfRcvFields[6][noind]  = g2*dt_adim;
  vectOfRcvFields[7][noind]  = 0.;//g3*dt_adim;
  vectOfRcvFields[8][noind]  = g4*dt_adim;
  vectOfRcvFields[9][noind]  = 0.;//g5*dt_adim;
  vectOfRcvFields[10][noind] = 0.;//g6*dt_adim;

}
