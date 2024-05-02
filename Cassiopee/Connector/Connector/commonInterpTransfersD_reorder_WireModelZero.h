for (E_Int noind = pt_deb; noind < pt_fin; noind++)
  {
    vectOfRcvFields[0][noind] = -1e06;
    vectOfRcvFields[1][noind] = -1e06;
    vectOfRcvFields[2][noind] = -1e06;
    vectOfRcvFields[3][noind] = -1e06;
    vectOfRcvFields[4][noind] = -1e06;
    if (nvars_loc==6){
      vectOfRcvFields[5][noind] = -1e06;
    }
  }
