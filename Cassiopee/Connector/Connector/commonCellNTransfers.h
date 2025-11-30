E_Float* ptrCoefs = donorCoefsF->begin();
E_Int indR, type, nocf;
E_Int indD0, indD, i, j, k, ncfLoc;
E_Int noi = 0; // compteur sur le tableau d indices donneur
E_Int sizecoefs = 0;
E_Float* cellNR = fieldsR[poscr];
E_Float* cellND = fieldsD[poscd];
for (E_Int noind = 0; noind < nbRcvPts; noind++)
{ 
  //
  // adressage indirect pour indR
  //
  indR = rcvPts[noind];
  type = types[noind];
  SIZECF(type, meshtype, sizecoefs);
  
  val = 0.;
  nocf = 0;

  // - pas de cellule masquee, ni interpolee dans la molecule d'interpolation, sauf si son cf associe est nul.
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
        val += ptrCoefs[kk-1]*cellND[indD0]*(2.-cellND[indD0]);        
      }
      cellNR[indR] = val;
      sizecoefs = ncfLoc;
      noi += ncfLoc+1;
      break;
      
    case 1: // injection
      indD0 = donorPts[noi];
      cellNR[indR] = ptrCoefs[0]*cellND[indD0]*(2.-cellND[indD0]);     
      sizecoefs = 1;
      noi += 1;
      break;
      
    case 2: // Structure Lineaire O2 par tetra
      indD0 = donorPts[noind];
      for (E_Int kk=0; kk<2; kk++)
        for (E_Int jj=0; jj<2; jj++)
          for (E_Int ii=0; ii<2; ii++)
          {
            indD = inD0 + ii + jj*imd + kk*imdjmd;
            val += ptrCoefs[nocf]*cellND[indD]*(2.-cellND[indD]);   
            nocf++;
          }
      cellNR[indR] = val;
      noi += 1;
      break;
    
    case 22:// O2CF 2D
      indD0 = donorPts[noind];
      for (E_Int jj=0; jj<2; jj++)
        for (E_Int ii=0; ii<2; ii++)
        {
          indD = indD0 + ii + jj*imd;
          val += ptrCoefs[nocf]*cellND[indD]*(2.-cellND[indD]);   
          nocf++;
        }      
      cellNR[indR] = val;
      noi += 1;
      break;
      
    case 3: // Lagrange O3
      indD0 = donorPts[noind];
      for (E_Int kk=0; kk<3; kk++)
        for (E_Int jj=0; jj<3; jj++)
          for (E_Int ii=0; ii<3; ii++)
          {
            indD = indD0 + ii + jj*imd + kk*imdjmd;
            val += ptrCoefs[ii]*ptrCoefs[jj+3]*ptrCoefs[kk+6]*cellND[indD]*(2.-cellND[indD]);             
          }
      cellNR[indR] = val;
      noi += 1;
      break;
      
    case 4: // Tetra O2
      indD0 = donorPts[noind];
      // indD0 est le no de l elt, et les coefs sont aux noeuds
      for (E_Int nov = 1; nov <= 4; nov++)
      {
        indD = ptrcnd[indD0*cnNfldD+nov-1]-1;
        val += cellND[indD]*(2.-cellND[indD])*ptrCoefs[nov-1];        
      }
      cellNR[indR] = val;      
      noi += 1;      
      break;
      
    case 44: // Lagrange O4
      indD0 = donorPts[noind];
      for (E_Int kk=0; kk<4; kk++)
        for (E_Int jj=0; jj<4; jj++)
          for (E_Int ii=0; ii<4; ii++)
          {
            indD = indD0 + ii + jj*imd + kk*imdjmd;
            val += ptrCoefs[ii]*ptrCoefs[jj+4]*ptrCoefs[kk+8]*cellND[indD]*(2.-cellND[indD]);         
          }
      cellNR[indR] = val;
      noi += 1;
      break;

    case 5: // Lagrange O5
      indD0 = donorPts[noind];
      for (E_Int kk=0; kk<5; kk++)
        for (E_Int jj=0; jj<5; jj++)
          for (E_Int ii=0; ii<5; ii++)
          {
            indD = indD0 + ii + jj*imd + kk*imdjmd;
            val += ptrCoefs[ii]*ptrCoefs[jj+5]*ptrCoefs[kk+10]*cellND[indD]*(2.-cellND[indD]);         
          }
      cellNR[indR] = val;
      noi += 1;
      break;
      
    default:;
  }
  ptrCoefs += sizecoefs;
} 
