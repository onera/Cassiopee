type = types[noind];
SIZECF(type, meshtype, sizecoefs);
E_Float val;
switch (type)
{
  case 0:  //  nuage de pts quelconque
    ncfLoc = donorPts[noi];// nb de pts pour la formule
    for (E_Int eq = eq_deb; eq < eq_fin; eq++)
    {
      E_Float* fieldR = vectOfRcvFields[eq];
      E_Float* fieldD = vectOfDnrFields[eq];
      val = 0.;
      for (E_Int kk = 1; kk <= ncfLoc; kk++)
      {
        indD0 = donorPts[noi+kk];
        val += ptrCoefs[kk-1]*fieldD[indD0];
      }
      fieldR[indR] = val;
    }
    sizecoefs = ncfLoc;
    noi += ncfLoc+1;
    break;

  case 1:
    indD0 = donorPts[noi];
    for (E_Int eq = eq_deb; eq < eq_fin; eq++)
    {
      E_Float* fieldR = vectOfRcvFields[eq];
      E_Float* fieldD = vectOfDnrFields[eq];
      fieldR[indR] = ptrCoefs[0]*fieldD[indD0];      
    }
    sizecoefs = 1;
    noi+=1;
    break;
    
  case 2: // Structure Lineaire O2 par tetra
    indD0 = donorPts[noind];
    k = indD0/imdjmd;
    j = (indD0-k*imdjmd)/imd;
    i = (indD0-j*imd-k*imdjmd);
    for (E_Int eq = eq_deb; eq < eq_fin; eq++)
    {
      E_Float* fieldR = vectOfRcvFields[eq];
      E_Float* fieldD = vectOfDnrFields[eq];
      val = 0.; E_Int nocf = 0;
      for (E_Int kk=0; kk<2; kk++)
        for (E_Int jj=0; jj<2; jj++)
          for (E_Int ii=0; ii<2; ii++)
          {
            indD = (i+ii)+(j+jj)*imd+(k+kk)*imdjmd;
            val += ptrCoefs[nocf]*fieldD[indD];
            nocf++;
          }
      fieldR[indR] = val;
    }
    noi += 1;
    break;
    
  case 22:// O2CF 2D
    indD0 = donorPts[noind];
    j = indD0/imd;
    i = indD0-j*imd;
    for (E_Int eq = eq_deb; eq < eq_fin; eq++)
    {
      E_Float* fieldR = vectOfRcvFields[eq];
      E_Float* fieldD = vectOfDnrFields[eq];
      val = 0.; E_Int nocf = 0;
      for (E_Int jj=0; jj<2; jj++)
        for (E_Int ii=0; ii<2; ii++)
        {
          indD = (i+ii)+(j+jj)*imd;
          val += ptrCoefs[nocf]*fieldD[indD];
          nocf++;
        }
      fieldR[indR] = val;
    }
    noi += 1;
    break;

  case 3: // Lagrange O3
    indD0 = donorPts[noind];
    k = indD0/imdjmd;
    j = (indD0-k*imdjmd)/imd;
    i = (indD0-j*imd-k*imdjmd);
    for (E_Int eq = eq_deb; eq < eq_fin; eq++)
    {
      E_Float* fieldR = vectOfRcvFields[eq];
      E_Float* fieldD = vectOfDnrFields[eq];
      val = 0.; 
      for (E_Int kk=0; kk<3; kk++)
        for (E_Int jj=0; jj<3; jj++)
          for (E_Int ii=0; ii<3; ii++)
          {
            indD = (i+ii)+(j+jj)*imd+(k+kk)*imdjmd;
            val += ptrCoefs[ii]*ptrCoefs[jj+3]*ptrCoefs[kk+6]*fieldD[indD];               
          }
      fieldR[indR]=val;
    }
    noi += 1;
    break;
 
  case 44: // Lagrange O4
    indD0 = donorPts[noind];
    k = indD0/imdjmd;
    j = (indD0-k*imdjmd)/imd;
    i = (indD0-j*imd-k*imdjmd);
    for (E_Int eq = eq_deb; eq < eq_fin; eq++)
    {
      E_Float* fieldR = vectOfRcvFields[eq];
      E_Float* fieldD = vectOfDnrFields[eq];
      val = 0.; 
      for (E_Int kk=0; kk<4; kk++)
        for (E_Int jj=0; jj<4; jj++)
          for (E_Int ii=0; ii<4; ii++)
          {
            indD = (i+ii)+(j+jj)*imd+(k+kk)*imdjmd;
            val += ptrCoefs[ii]*ptrCoefs[jj+4]*ptrCoefs[kk+8]*fieldD[indD];               
          }
      fieldR[indR]=val;
    }
    noi += 1;
    break;
     
  case 4: // Tetra O2
    indD0 = donorPts[noind];
    // indD0 est le no de l elt, et les coefs sont aux noeuds
    for (E_Int eq = eq_deb; eq < eq_fin; eq++)
    {
      E_Float* fieldR = vectOfRcvFields[eq];
      E_Float* fieldD = vectOfDnrFields[eq];
      val = 0.;
      for (E_Int nov = 1; nov <= 4; nov++)
      {
        E_Int indv = ptrcnd[indD0*cnNfldD+nov-1]-1;
        val += fieldD[indv]*ptrCoefs[nov-1];
      }
      fieldR[indR]=val;
    }  
    noi += 1;      
    break;
      
  case 5: // Lagrange O5
    indD0 = donorPts[noind];
    k = indD0/imdjmd;
    j = (indD0-k*imdjmd)/imd;
    i = (indD0-j*imd-k*imdjmd);
    for (E_Int eq = eq_deb; eq < eq_fin; eq++)
    {
      E_Float* fieldR = vectOfRcvFields[eq];
      E_Float* fieldD = vectOfDnrFields[eq];
      val= 0.; 
      for (E_Int kk=0; kk<5; kk++)
        for (E_Int jj=0; jj<5; jj++)
          for (E_Int ii=0; ii<5; ii++)
          {
            indD = (i+ii)+(j+jj)*imd+(k+kk)*imdjmd;
            val += ptrCoefs[ii]*ptrCoefs[jj+5]*ptrCoefs[kk+10]*fieldD[indD];               
          }
      fieldR[indR] = val;
    }
    noi += 1;
    break;
      
  default: ;
}
