if(iscalcMacro_LBM==1 ){
  for (E_Int noind = pt_deb; noind < pt_fin; noind++){
    indR   = rcvPts[noind];
    vectOfmacroRcvFields[0][indR]=0.;
    vectOfmacroRcvFields[1][indR]=0.;
    vectOfmacroRcvFields[2][indR]=0.;
    vectOfmacroRcvFields[3][indR]=0.;
    for (E_Int ind_f = 0; ind_f<nvars_loc;ind_f++){
      vectOfmacroRcvFields[0][indR] += vectOfRcvFields[ind_f][noind];
      vectOfmacroRcvFields[1][indR] += vectOfRcvFields[ind_f][noind]*ipt_cvel[ind_f            ];
      vectOfmacroRcvFields[2][indR] += vectOfRcvFields[ind_f][noind]*ipt_cvel[ind_f+  nvars_loc];
      vectOfmacroRcvFields[3][indR] += vectOfRcvFields[ind_f][noind]*ipt_cvel[ind_f+2*nvars_loc];
    }
  }
 }

