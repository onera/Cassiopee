
E_Float dt_adim = ipt_param_realR[ NoD ][ DTC   ];
E_Float gam     = ipt_param_realR[ NoD ][ GAMMA ];
E_Float Rgp     = ipt_param_realR[ NoD ][ CVINF ]*(gam-1);
E_Float t_inf   = ipt_param_realR[ NoD ][ TINF  ];
E_Float ro_inf  = ipt_param_realR[ NoD ][ ROINF ];
E_Float c0_ref  = sqrt(gam*Rgp*t_inf);
E_Float u_scale = 1./(sqrt(3.)*c0_ref);

E_Float v1, v2, v3, v4;
E_Float g1, g2, g3, g4, g5, g6;

for (E_Int noind = pt_deb; noind < pt_fin; noind++)
{

  indR   = rcvPts[noind];

  v1 = vectOfRcvFields[0][indR]; //ro
  v2 = vectOfRcvFields[1][indR]; //vx
  v3 = vectOfRcvFields[2][indR]; //vy
  v4 = vectOfRcvFields[3][indR]; //vz
  // La temperature est modifiee par une BC
  g1 = vectOfRcvFields[5][indR]; //Sxx
  g2 = vectOfRcvFields[6][indR]; //Sxy
  g3 = vectOfRcvFields[7][indR]; //Sxz
  g4 = vectOfRcvFields[8][indR]; //Syy
  g5 = vectOfRcvFields[9][indR]; //Syz
  g6 = vectOfRcvFields[10][indR]; //Szz

  //printf("dimNS %.14f %.14f %.14f %d \n", v1, ro_inf, u_scale, indR);
  vectOfRcvFields[0][indR] = v1/ro_inf;
  vectOfRcvFields[1][indR] = v2*u_scale;
  vectOfRcvFields[2][indR] = v3*u_scale;
  vectOfRcvFields[3][indR] = 0.;//v4*u_scale;

  vectOfRcvFields[5][indR]  = g1*dt_adim;
  vectOfRcvFields[6][indR]  = g2*dt_adim;
  vectOfRcvFields[7][indR]  = 0.;//g3*dt_adim;
  vectOfRcvFields[8][indR]  = g4*dt_adim;
  vectOfRcvFields[9][indR]  = 0.;//g5*dt_adim;
  vectOfRcvFields[10][indR] = 0.;//g6*dt_adim;

}
