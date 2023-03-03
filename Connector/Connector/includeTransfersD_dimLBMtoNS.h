
E_Float gam     = ipt_param_realD[ NoD ][ GAMMA ];
E_Float Rgp     = ipt_param_realD[ NoD ][ CVINF ]*(gam-1);
E_Float t_inf   = ipt_param_realD[ NoD ][ TINF  ];
E_Float ro_inf  = ipt_param_realD[ NoD ][ ROINF ];
E_Float c0_ref  = sqrt(gam*Rgp*t_inf);
E_Float u_scale = sqrt(3.)*c0_ref;

E_Float v1, v2, v3, v4;
E_Float pp, t_from_lbm;

for (E_Int noind = pt_deb; noind < pt_fin; noind++)
{

  noind   = rcvPts[noind];

  v1 = vectOfRcvFields[0][noind]; //ro
  v2 = vectOfRcvFields[1][noind]; //vx
  v3 = vectOfRcvFields[2][noind]; //vy
  v4 = vectOfRcvFields[3][noind]; //vz
  // La temperature est modifiee par une BC
  pp = ro_inf*gam*Rgp*t_inf*(-0.4/1.4 + v1 );
  t_from_lbm = pp/(Rgp*v1*ro_inf);
  //cout << t_from_lbm << endl;

  vectOfRcvFields[0][noind] = v1*ro_inf;
  vectOfRcvFields[1][noind] = v2*u_scale;
  vectOfRcvFields[2][noind] = v3*u_scale;
  vectOfRcvFields[3][noind] = 0.;//v4*u_scale;
  //vectOfRcvFields[4][noind] = t_from_lbm;

}
