E_Int bestFound = 1;
E_Int inddummy = -1;//a laisser a -1 absolument...pour k6compvolofstructcell
if (found > 0)
{      
  foundSav = 1;
  // determination du cellN de la cellule d'interpolation
  if (tmpType == 2 || tmpType == 3 || tmpType == 5)
  {
    firstCorner = tmpIndi[0]; 
    // calcul de cellvol
    k6compvolofstructcell_(ni, nj, nk, inddummy, firstCorner, oneField.begin(posx0),
                           oneField.begin(posy0), oneField.begin(posz0),vol);
    if (penalty == 1 && isBorder == 1) vol += 1.e3;
    if (vol < best+K_CONST::E_GEOM_CUTOFF) 
    {
      if ( isExtrapolation == 1 ) //si plusieurs donneurs par extrapolation possibles : extrap du plus proche
      {
        E_Float sumCf = 0.;
        for (E_Int nocf = 0; nocf < tmpCf.getSize(); nocf++)
          sumCf += K_FUNC::E_abs(tmpCf[nocf]);
        if ( sumCf < sumCoefMin ) sumCoefMin = sumCf; 
        else bestFound = 0;
      }

      if ( bestFound == 1)
      {
        best = vol;
        type = tmpType;
        donorIndices[0] = tmpIndi[0];
        for (E_Int nocf = 0; nocf < tmpCf.getSize(); nocf++)
          donorCoefs[nocf] = tmpCf[nocf];
        noblk = noz+1;          
      }
    }
  }//types 2, 3 et 5

  else if (tmpType == 4) 
  {
    if (meshtype == 1) 
    {
      //printf("getInterpolationCell: interpolation type 4 not implemented for structured donor zones.\n"); 
      return -1;
    }        
    noei = tmpIndi[0];      
    k6compvoloftetracell_(oneField.getSize(), 
                          (*cnloc)(noei,1)-1, (*cnloc)(noei,2)-1, (*cnloc)(noei,3)-1, (*cnloc)(noei,4)-1, 
                          oneField.begin(posx0), oneField.begin(posy0),oneField.begin(posz0), vol);
    if (penalty == 1 && isBorder == 1) vol += 1.e3;
    if (vol < best+K_CONST::E_GEOM_CUTOFF) 
    {
      if ( isExtrapolation == 1 )//si plusieurs donneurs par extrapolation possibles : extrap du plus proche
      {
        E_Float sumCf = 0.;
        for (E_Int nocf = 0; nocf < tmpCf.getSize(); nocf++)
          sumCf += K_FUNC::E_abs(tmpCf[nocf]);
        if ( sumCf < sumCoefMin ) sumCoefMin = sumCf; 
        else bestFound = 0;
      }
      if (bestFound == 1)
      {
        best = vol;
        type = tmpType;
        donorIndices[0] = tmpIndi[0];
        for (E_Int nocf = 0; nocf < tmpCf.getSize(); nocf++)
          donorCoefs[nocf] = tmpCf[nocf];
        noblk = noz+1;          
      }
    }
  }
  else if (tmpType == 22) // type 2 en 2D
  {
    firstCorner = tmpIndi[0]; 
    // calcul de cellvol
    vol = K_METRIC::compVolOfStructCell2D(ni,nj, oneField.begin(posx0),oneField.begin(posy0), oneField.begin(posz0),
                                          -1, firstCorner);
    if (penalty == 1 && isBorder == 1) vol += 1.e3;
    if (vol < best+K_CONST::E_GEOM_CUTOFF) 
    {
      if ( isExtrapolation == 1 )//si plusieurs donneurs par extrapolation possibles : extrap du plus proche
      {
        E_Float sumCf = 0.;
        for (E_Int nocf = 0; nocf < 4; nocf++)
          sumCf += K_FUNC::E_abs(tmpCf[nocf]);
        if ( sumCf < sumCoefMin ) sumCoefMin = sumCf; 
        else bestFound = 0;
      }

      if ( bestFound == 1)
      {
        best = vol;
        type = tmpType;
        donorIndices[0] = tmpIndi[0];
        //report des 8 coefs sur les 4 premiers sommets
        for (E_Int nocf = 0; nocf < 4; nocf++)
        {
          donorCoefs[nocf] = tmpCf[nocf]+tmpCf[nocf+4];
          donorCoefs[nocf+4] = 0.;
        }
        noblk = noz+1;          
      }
    }
  }
  else 
  {
    printf("getInterpolationCell: invalid interpolation type.\n");
    return -2;
  }
}// found > 0

