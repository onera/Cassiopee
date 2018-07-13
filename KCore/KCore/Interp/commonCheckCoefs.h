// on verifie si tous les coefs sont nuls sauf 1. Si oui -> on passe en type 1
E_Int nozd = noDonorBlk-1;
E_Float cutOffCf = 1.e-10;//1.e-12;
if (a3[nozd] != NULL) //structure
{
  ni = *(E_Int*)a1[nozd]; nj = *(E_Int*)a2[nozd]; nk = *(E_Int*)a3[nozd]; 
  ninj = ni*nj;
}
else
{ 
  cnloc = (FldArrayI*)a1[nozd]; // non structure      
  ni = -1; nj = -1; nk = -1; ninj = -1;      
}

d = 0;
switch (type)
{
  case 2:
    indcell = donorIndices[0]; 
    for (E_Int kk = 0; kk < 2; kk++)
      for (E_Int jj = 0; jj < 2; jj++)
        for (E_Int ii = 0; ii < 2; ii++)
        {
          indLoc = indcell+ii+jj*ni+kk*ninj;
          coefLoc = donorCoefs[ii+jj*2+kk*4];
          d =  K_FUNC::E_abs(coefLoc) > cutOffCf;
          nbNonZeroCf += d;
          if (d == 1){ indNonZero = indLoc; coefNonZero = coefLoc;}
        }
    break;

  case 22:
    indcell = donorIndices[0]; 
      for (E_Int jj = 0; jj < 2; jj++)
        for (E_Int ii = 0; ii < 2; ii++)
        {
          indLoc = indcell+ii+jj*ni;
          coefLoc = donorCoefs[ii+jj*2];
          d =  K_FUNC::E_abs(coefLoc) > cutOffCf;
          nbNonZeroCf += d;
          if (d == 1){ indNonZero = indLoc; coefNonZero = coefLoc;}
        }
    break;

  case 3:
    indcell = donorIndices[0];  
    for (E_Int kk = 0; kk < 3; kk++)
      for (E_Int jj = 0; jj < 3; jj++)
        for (E_Int ii = 0; ii < 3; ii++)
        {
          indLoc = indcell+ii+jj*ni+kk*ninj; 
          coefLoc = donorCoefs[ii]*donorCoefs[jj+3]*donorCoefs[kk+6];
          d = K_FUNC::E_abs(coefLoc) > cutOffCf;
          nbNonZeroCf += d; 
          if (d == 1) {indNonZero = indLoc;coefNonZero = coefLoc;}
        }
    break;
   
  case 4:
    noei = donorIndices[0];
    for (E_Int nov = 1; nov <= 4; nov++)
    {
      indLoc = (*cnloc)(noei,nov)-1; 
      coefLoc = donorCoefs[nov-1];
      d = K_FUNC::E_abs(coefLoc) > cutOffCf;
      nbNonZeroCf += d; 
      if (d == 1) {indNonZero = indLoc; coefNonZero = coefLoc;}
    }
    break;
      
  case 5: 
    indcell = donorIndices[0]; 
    for (E_Int kk = 0; kk < 5; kk++)
      for (E_Int jj = 0; jj < 5; jj++)
        for (E_Int ii = 0; ii < 5; ii++)
        {
          indLoc =  indcell+ii+jj*ni+kk*ninj; 
          coefLoc = donorCoefs[ii]*donorCoefs[jj+5]*donorCoefs[kk+10];
          d = K_FUNC::E_abs(coefLoc) > cutOffCf;
          nbNonZeroCf += d; 
          if (d == 1) {indNonZero = indLoc; coefNonZero = coefLoc;}
        }
    break;
    
  default:
    break;      
}
if ( nbNonZeroCf == 1 )
{
  type = 1;
  donorIndices[0] = indNonZero;
  donorCoefs[0] = coefNonZero; 
}
