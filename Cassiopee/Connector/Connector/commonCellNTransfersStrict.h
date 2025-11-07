
type = types[noind];
SIZECF(type, meshtype, sizecoefs);

nocf = 0;
E_Float val = 1.;
E_Float threshold = 1.e-11; // thresold on coeff to ignore cellN

// - pas de cellule masquee, ni interpolee dans la molecule d'interpolation, 
// sauf si son cf associe est nul au sens de thresold.
// - somme des cf = 1.
// -----------------
// cellN*(2-cellN) renvoie 0 si cellN = 0 ou 2 (pt masque ou interpolee) et 1 si cellN =1 (pt calcule)
switch (type)
{
  case 0:  //  nuage de pts quelconque
    ncfLoc = donorPts[noi];// nb de pts pour la formule
    for (E_Int kk = 1; kk <= ncfLoc; kk++)
    {
      indD0 = donorPts[noi+kk];
      if (K_FUNC::E_abs(ptrCoefs[kk-1]) > threshold) { val *= cellND[indD0]*(2.-cellND[indD0]); }
    }
    if ( val > 0.5 ) cellNR[indR] = 1.;//egal a 1 
    else cellNR[indR] = 0.;
    sizecoefs = ncfLoc;
    noi += ncfLoc+1;
    break;
    
  case 1: // injection
    indD0 = donorPts[noi];
    val = cellND[indD0]*(2.-cellND[indD0]); 
    if ( val > 0.5 ) cellNR[indR] = 1.;//egal a 1 
    else cellNR[indR] = 0.;
    sizecoefs = 1;
    noi += 1;
    break;
    
  case 2: // Structure Lineaire O2 par tetra
    indD0 = donorPts[noind];
    for (E_Int kk=0; kk<2; kk++)
      for (E_Int jj=0; jj<2; jj++)
        for (E_Int ii=0; ii<2; ii++)
        {
          indD = indD0 + ii + jj*imd + kk*imdjmd;
          if (K_FUNC::E_abs(ptrCoefs[nocf]) > threshold) { val *= cellND[indD]*(2.-cellND[indD]); }
          nocf++;
        }    
    if ( val > 0.5 ) cellNR[indR] = 1.;//egal a 1 
    else cellNR[indR] = 0.;
    noi += 1;
    break;
    
  case 22:// O2CF 2D
    indD0 = donorPts[noind];
    for (E_Int jj=0; jj<2; jj++)
      for (E_Int ii=0; ii<2; ii++)
      {
        indD = indD0 + ii + jj*imd;
        if (K_FUNC::E_abs(ptrCoefs[nocf]) > threshold) { val *= cellND[indD]*(2.-cellND[indD]); }
        nocf++;
      }      
    if ( val > 0.5 ) cellNR[indR] = 1.;//egal a 1 
    else cellNR[indR] = 0.;
    noi += 1;
    break;
    
  case 3: // Lagrange O3
    indD0 = donorPts[noind];    
    for (E_Int kk=0; kk<3; kk++)
      for (E_Int jj=0; jj<3; jj++)
        for (E_Int ii=0; ii<3; ii++)
        {
          indD = indD0 + ii + jj*imd + kk*imdjmd;
          if (K_FUNC::E_abs(ptrCoefs[ii]*ptrCoefs[jj+3]*ptrCoefs[kk+6]) > threshold) { val *= cellND[indD]*(2.-cellND[indD]); }
        }
    if ( val > 0.5 ) cellNR[indR] = 1.;//egal a 1 
    else cellNR[indR] = 0.;
    noi += 1;
    break;
    
  case 4: // Tetra O2
    indD0 = donorPts[noind];
    // indD0 est le no de l elt, et les coefs sont aux noeuds
    for (E_Int nov = 1; nov <= 4; nov++)
    {
      indD = ptrcnd[indD0*cnNfldD+nov-1]-1;
      if (K_FUNC::E_abs(ptrCoefs[nov-1]) > threshold) { val *= cellND[indD]*(2.-cellND[indD]); }     
    }
    if ( val > 0.5 ) cellNR[indR] = 1.;//egal a 1 
    else cellNR[indR] = 0.;
    noi += 1;
    break;

  case 44: // Lagrange O4
    indD0 = donorPts[noind];  //car type 0 est toujour traité en dernier. Sinon noind pas valable
    for (E_Int kk=0; kk<4; kk++)
      for (E_Int jj=0; jj<4; jj++)
        for (E_Int ii=0; ii<4; ii++)
        {
          indD = indD0 + ii + jj*imd + kk*imdjmd;
          if (K_FUNC::E_abs(ptrCoefs[ii]*ptrCoefs[jj+4]*ptrCoefs[kk+8]) > threshold) { val *= cellND[indD]*(2.-cellND[indD]); }
        }
    if ( val > 0.5 ) cellNR[indR] = 1.;//egal a 1 
    else cellNR[indR] = 0.;
    noi += 1;
    break;

  case 5: // Lagrange O5
    indD0 = donorPts[noind];
    for (E_Int kk=0; kk<5; kk++)
      for (E_Int jj=0; jj<5; jj++)
        for (E_Int ii=0; ii<5; ii++)
        {
          indD = indD0 + ii + jj*imd + kk*imdjmd;
          if (K_FUNC::E_abs(ptrCoefs[ii]*ptrCoefs[jj+5]*ptrCoefs[kk+10]) > threshold) { val *= cellND[indD]*(2.-cellND[indD]); }
        }
    if ( val > 0.5 ) cellNR[indR] = 1.;//egal a 1 
    else cellNR[indR] = 0.;
    noi += 1;
    break;
    
  default:;
}
ptrCoefs += sizecoefs;
