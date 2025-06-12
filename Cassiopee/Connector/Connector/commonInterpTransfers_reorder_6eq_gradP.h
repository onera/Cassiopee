E_Int ind000,ind100,ind010,ind110, ind001,ind101,ind011,ind111;
E_Int ind00 ,ind10 ,ind01 ,ind11, ind02,ind03;
E_Float val0, val1, val2, val3, val4, val5;

switch (type)
{
  case 0:  //  nuage de pts quelconque
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      indR   = rcvPts[noind];
      ncfLoc = donorPts[noi];// nb de pts pour la formule
      val0=0.; val1=0.; val2=0.; val3=0.; val4=0.; val5=0.;
      for (E_Int kk = 1; kk <= ncfLoc; kk++)
      {
        indD0         = donorPts[noi+kk];
        val0 += ptrCoefs[ indCoef + kk-1]*vectOfGradDnrFields[0][indD0];
        val1 += ptrCoefs[ indCoef + kk-1]*vectOfGradDnrFields[1][indD0];
        val2 += ptrCoefs[ indCoef + kk-1]*vectOfGradDnrFields[2][indD0];
        val3 += ptrCoefs[ indCoef + kk-1]*vectOfGradDnrFields[3][indD0];
        val4 += ptrCoefs[ indCoef + kk-1]*vectOfGradDnrFields[4][indD0];
        val5 += ptrCoefs[ indCoef + kk-1]*vectOfGradDnrFields[5][indD0];
      }
      vectOfGradRcvFields[0][indR] = val0;
      vectOfGradRcvFields[1][indR] = val1;
      vectOfGradRcvFields[2][indR] = val2;
      vectOfGradRcvFields[3][indR] = val3;
      vectOfGradRcvFields[4][indR] = val4;
      vectOfGradRcvFields[5][indR] = val5;
      sizecoefs = ncfLoc;
      noi      += ncfLoc+1;
      indCoef  += sizecoefs;
    }
    break;

  case 1:
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      indR   = rcvPts[noind];
      indD0  = donorPts[noind];  //car type 0 est toujour traite en dernier. Sinon noind pas valable
    
      vectOfGradRcvFields[0][indR] = vectOfGradDnrFields[0][indD0];
      vectOfGradRcvFields[1][indR] = vectOfGradDnrFields[1][indD0];
      vectOfGradRcvFields[2][indR] = vectOfGradDnrFields[2][indD0];
      vectOfGradRcvFields[3][indR] = vectOfGradDnrFields[3][indD0];
      vectOfGradRcvFields[4][indR] = vectOfGradDnrFields[4][indD0];
      vectOfGradRcvFields[5][indR] = vectOfGradDnrFields[5][indD0];
    }
    break;
    
  case 2: // Structure Lineaire O2 par tetra
  // #ifdef _OPENMP4
  //    #pragma omp simd
  // #endif
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      indR   = rcvPts[noind];
      ind000 = donorPts[noind];
      ind100 = ind000 +1; 
      ind010 = ind000 +imd; 
      ind110 = ind100 +imd; 
      ind001 = ind000+imdjmd;
      ind101 = ind100+imdjmd;
      ind011 = ind010+imdjmd;
      ind111 = ind110+imdjmd;

      val0  = ptrCoefs[ indCoef     ]*vectOfGradDnrFields[0][ind000];
      val0 += ptrCoefs[ indCoef + 1 ]*vectOfGradDnrFields[0][ind100];
      val0 += ptrCoefs[ indCoef + 2 ]*vectOfGradDnrFields[0][ind010];
      val0 += ptrCoefs[ indCoef + 3 ]*vectOfGradDnrFields[0][ind110];
      val0 += ptrCoefs[ indCoef + 4 ]*vectOfGradDnrFields[0][ind001];
      val0 += ptrCoefs[ indCoef + 5 ]*vectOfGradDnrFields[0][ind101];
      val0 += ptrCoefs[ indCoef + 6 ]*vectOfGradDnrFields[0][ind011];
      val0 += ptrCoefs[ indCoef + 7 ]*vectOfGradDnrFields[0][ind111];

      val1  = ptrCoefs[ indCoef     ]*vectOfGradDnrFields[1][ind000];
      val1 += ptrCoefs[ indCoef + 1 ]*vectOfGradDnrFields[1][ind100];
      val1 += ptrCoefs[ indCoef + 2 ]*vectOfGradDnrFields[1][ind010];
      val1 += ptrCoefs[ indCoef + 3 ]*vectOfGradDnrFields[1][ind110];
      val1 += ptrCoefs[ indCoef + 4 ]*vectOfGradDnrFields[1][ind001];
      val1 += ptrCoefs[ indCoef + 5 ]*vectOfGradDnrFields[1][ind101];
      val1 += ptrCoefs[ indCoef + 6 ]*vectOfGradDnrFields[1][ind011];
      val1 += ptrCoefs[ indCoef + 7 ]*vectOfGradDnrFields[1][ind111];

      val2  = ptrCoefs[ indCoef     ]*vectOfGradDnrFields[2][ind000];
      val2 += ptrCoefs[ indCoef + 1 ]*vectOfGradDnrFields[2][ind100];
      val2 += ptrCoefs[ indCoef + 2 ]*vectOfGradDnrFields[2][ind010];
      val2 += ptrCoefs[ indCoef + 3 ]*vectOfGradDnrFields[2][ind110];
      val2 += ptrCoefs[ indCoef + 4 ]*vectOfGradDnrFields[2][ind001];
      val2 += ptrCoefs[ indCoef + 5 ]*vectOfGradDnrFields[2][ind101];
      val2 += ptrCoefs[ indCoef + 6 ]*vectOfGradDnrFields[2][ind011];
      val2 += ptrCoefs[ indCoef + 7 ]*vectOfGradDnrFields[2][ind111];

      val3  = ptrCoefs[ indCoef     ]*vectOfGradDnrFields[3][ind000];
      val3 += ptrCoefs[ indCoef + 1 ]*vectOfGradDnrFields[3][ind100];
      val3 += ptrCoefs[ indCoef + 2 ]*vectOfGradDnrFields[3][ind010];
      val3 += ptrCoefs[ indCoef + 3 ]*vectOfGradDnrFields[3][ind110];
      val3 += ptrCoefs[ indCoef + 4 ]*vectOfGradDnrFields[3][ind001];
      val3 += ptrCoefs[ indCoef + 5 ]*vectOfGradDnrFields[3][ind101];
      val3 += ptrCoefs[ indCoef + 6 ]*vectOfGradDnrFields[3][ind011];
      val3 += ptrCoefs[ indCoef + 7 ]*vectOfGradDnrFields[3][ind111];

      val4  = ptrCoefs[ indCoef     ]*vectOfGradDnrFields[4][ind000];
      val4 += ptrCoefs[ indCoef + 1 ]*vectOfGradDnrFields[4][ind100];
      val4 += ptrCoefs[ indCoef + 2 ]*vectOfGradDnrFields[4][ind010];
      val4 += ptrCoefs[ indCoef + 3 ]*vectOfGradDnrFields[4][ind110];
      val4 += ptrCoefs[ indCoef + 4 ]*vectOfGradDnrFields[4][ind001];
      val4 += ptrCoefs[ indCoef + 5 ]*vectOfGradDnrFields[4][ind101];
      val4 += ptrCoefs[ indCoef + 6 ]*vectOfGradDnrFields[4][ind011];
      val4 += ptrCoefs[ indCoef + 7 ]*vectOfGradDnrFields[4][ind111];

      val5  = ptrCoefs[ indCoef     ]*vectOfGradDnrFields[5][ind000];
      val5 += ptrCoefs[ indCoef + 1 ]*vectOfGradDnrFields[5][ind100];
      val5 += ptrCoefs[ indCoef + 2 ]*vectOfGradDnrFields[5][ind010];
      val5 += ptrCoefs[ indCoef + 3 ]*vectOfGradDnrFields[5][ind110];
      val5 += ptrCoefs[ indCoef + 4 ]*vectOfGradDnrFields[5][ind001];
      val5 += ptrCoefs[ indCoef + 5 ]*vectOfGradDnrFields[5][ind101];
      val5 += ptrCoefs[ indCoef + 6 ]*vectOfGradDnrFields[5][ind011];
      val5 += ptrCoefs[ indCoef + 7 ]*vectOfGradDnrFields[5][ind111];

      vectOfGradRcvFields[0][indR] = val0;
      vectOfGradRcvFields[1][indR] = val1;
      vectOfGradRcvFields[2][indR] = val2;
      vectOfGradRcvFields[3][indR] = val3;
      vectOfGradRcvFields[4][indR] = val4;
      vectOfGradRcvFields[5][indR] = val5;
      indCoef  += 8;
    }
    break;
    
  case 22:// O2CF 2D
// #ifdef _OPENMP4
//     #pragma omp simd
// #endif
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      indR  = rcvPts[noind];
      ind00 = donorPts[noind];
      ind10 = ind00 +1; 
      ind01 = ind00 +imd; 
      ind11 = ind10 +imd; 

      val0  = ptrCoefs[ indCoef     ]*vectOfGradDnrFields[0][ind00];
      val0 += ptrCoefs[ indCoef + 1 ]*vectOfGradDnrFields[0][ind10];
      val0 += ptrCoefs[ indCoef + 2 ]*vectOfGradDnrFields[0][ind01];
      val0 += ptrCoefs[ indCoef + 3 ]*vectOfGradDnrFields[0][ind11];

      val1  = ptrCoefs[ indCoef     ]*vectOfGradDnrFields[1][ind00];
      val1 += ptrCoefs[ indCoef + 1 ]*vectOfGradDnrFields[1][ind10];
      val1 += ptrCoefs[ indCoef + 2 ]*vectOfGradDnrFields[1][ind01];
      val1 += ptrCoefs[ indCoef + 3 ]*vectOfGradDnrFields[1][ind11];

      val2  = ptrCoefs[ indCoef     ]*vectOfGradDnrFields[2][ind00];
      val2 += ptrCoefs[ indCoef + 1 ]*vectOfGradDnrFields[2][ind10];
      val2 += ptrCoefs[ indCoef + 2 ]*vectOfGradDnrFields[2][ind01];
      val2 += ptrCoefs[ indCoef + 3 ]*vectOfGradDnrFields[2][ind11];

      val3  = ptrCoefs[ indCoef     ]*vectOfGradDnrFields[3][ind00];
      val3 += ptrCoefs[ indCoef + 1 ]*vectOfGradDnrFields[3][ind10];
      val3 += ptrCoefs[ indCoef + 2 ]*vectOfGradDnrFields[3][ind01];
      val3 += ptrCoefs[ indCoef + 3 ]*vectOfGradDnrFields[3][ind11];

      val4  = ptrCoefs[ indCoef     ]*vectOfGradDnrFields[4][ind00];
      val4 += ptrCoefs[ indCoef + 1 ]*vectOfGradDnrFields[4][ind10];
      val4 += ptrCoefs[ indCoef + 2 ]*vectOfGradDnrFields[4][ind01];
      val4 += ptrCoefs[ indCoef + 3 ]*vectOfGradDnrFields[4][ind11];

      val5  = ptrCoefs[ indCoef     ]*vectOfGradDnrFields[5][ind00];
      val5 += ptrCoefs[ indCoef + 1 ]*vectOfGradDnrFields[5][ind10];
      val5 += ptrCoefs[ indCoef + 2 ]*vectOfGradDnrFields[5][ind01];
      val5 += ptrCoefs[ indCoef + 3 ]*vectOfGradDnrFields[5][ind11];

      vectOfGradRcvFields[0][indR] = val0;
      vectOfGradRcvFields[1][indR] = val1;
      vectOfGradRcvFields[2][indR] = val2;
      vectOfGradRcvFields[3][indR] = val3;
      vectOfGradRcvFields[4][indR] = val4;
      vectOfGradRcvFields[5][indR] = val5;
      indCoef  += 4;
    }
    break;

  case 3: // Lagrange O3
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      indR  = rcvPts[noind];
      indD0 = donorPts[noind];  //car type 0 est toujour traité en dernier. Sinon noind pas valable
      k     = indD0/imdjmd;
      j     = (indD0-k*imdjmd)/imd;
      i     = (indD0-j*imd-k*imdjmd);
      val0=0.; val1=0.; val2=0.; val3=0.; val4=0.; val5=0.;

      for (E_Int kk=0; kk<3; kk++)
        for (E_Int jj=0; jj<3; jj++)
          for (E_Int ii=0; ii<3; ii++)
          {
            indD = (i+ii)+(j+jj)*imd+(k+kk)*imdjmd;
            val0 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfGradDnrFields[0][indD];               
            val1 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfGradDnrFields[1][indD];               
            val2 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfGradDnrFields[2][indD];               
            val3 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfGradDnrFields[3][indD];               
            val4 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfGradDnrFields[4][indD];               
            val5 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfGradDnrFields[5][indD];               
          }
      vectOfGradRcvFields[0][indR] = val0;
      vectOfGradRcvFields[1][indR] = val1;
      vectOfGradRcvFields[2][indR] = val2;
      vectOfGradRcvFields[3][indR] = val3;
      vectOfGradRcvFields[4][indR] = val4;
      vectOfGradRcvFields[5][indR] = val5;
      noi      += 1;
      indCoef  += sizecoefs;
    }
    break;
    
  case 44: // Lagrange O4
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      indR  = rcvPts[noind];
      indD0 = donorPts[noind];  //car type 0 est toujour traité en dernier. Sinon noind pas valable
      k     = indD0/imdjmd;
      j     = (indD0-k*imdjmd)/imd;
      i     = (indD0-j*imd-k*imdjmd);
      val0=0.; val1=0.; val2=0.; val3=0.; val4=0.; val5=0.;

      for (E_Int kk=0; kk<4; kk++)
        for (E_Int jj=0; jj<4; jj++)
          for (E_Int ii=0; ii<4; ii++)
          {
            indD = (i+ii)+(j+jj)*imd+(k+kk)*imdjmd;
            val0 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+4]*ptrCoefs[ indCoef + kk+8]*vectOfGradDnrFields[0][indD];               
            val1 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+4]*ptrCoefs[ indCoef + kk+8]*vectOfGradDnrFields[1][indD];               
            val2 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+4]*ptrCoefs[ indCoef + kk+8]*vectOfGradDnrFields[2][indD];               
            val3 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+4]*ptrCoefs[ indCoef + kk+8]*vectOfGradDnrFields[3][indD];               
            val4 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+4]*ptrCoefs[ indCoef + kk+8]*vectOfGradDnrFields[4][indD];               
            val5 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+4]*ptrCoefs[ indCoef + kk+8]*vectOfGradDnrFields[5][indD];               
          }
      vectOfGradRcvFields[0][indR] = val0;
      vectOfGradRcvFields[1][indR] = val1;
      vectOfGradRcvFields[2][indR] = val2;
      vectOfGradRcvFields[3][indR] = val3;
      vectOfGradRcvFields[4][indR] = val4;
      vectOfGradRcvFields[5][indR] = val5;
      noi      += 1;
      indCoef  += sizecoefs;
    }
    break;
  
  case 4: // Tetra O2
// #ifdef _OPENM4
//     #pragma omp simd
// #endif
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      indR  = rcvPts[noind];
      indD0 = donorPts[noind];  //car type 0 est toujour traité en dernier. Sinon noind pas valable
      // indD0 est le no de l elt, et les coefs sont aux noeuds
    
      ind00 = ptrcnd[indD0*cnNfldD   ] -1;
      ind01 = ptrcnd[indD0*cnNfldD +1] -1;
      ind02 = ptrcnd[indD0*cnNfldD +2] -1;
      ind03 = ptrcnd[indD0*cnNfldD +3] -1;
      val0  = ptrCoefs[ indCoef   ]*vectOfGradDnrFields[0][ind00];
      val0 += ptrCoefs[ indCoef +1]*vectOfGradDnrFields[0][ind01];
      val0 += ptrCoefs[ indCoef +2]*vectOfGradDnrFields[0][ind02];
      val0 += ptrCoefs[ indCoef +3]*vectOfGradDnrFields[0][ind03];
      val1  = ptrCoefs[ indCoef   ]*vectOfGradDnrFields[1][ind00];
      val1 += ptrCoefs[ indCoef +1]*vectOfGradDnrFields[1][ind01];
      val1 += ptrCoefs[ indCoef +2]*vectOfGradDnrFields[1][ind02];
      val1 += ptrCoefs[ indCoef +3]*vectOfGradDnrFields[1][ind03];
      val2  = ptrCoefs[ indCoef   ]*vectOfGradDnrFields[2][ind00];
      val2 += ptrCoefs[ indCoef +1]*vectOfGradDnrFields[2][ind01];
      val2 += ptrCoefs[ indCoef +2]*vectOfGradDnrFields[2][ind02];
      val2 += ptrCoefs[ indCoef +3]*vectOfGradDnrFields[2][ind03];
      val3  = ptrCoefs[ indCoef   ]*vectOfGradDnrFields[3][ind00];
      val3 += ptrCoefs[ indCoef +1]*vectOfGradDnrFields[3][ind01];
      val3 += ptrCoefs[ indCoef +2]*vectOfGradDnrFields[3][ind02];
      val3 += ptrCoefs[ indCoef +3]*vectOfGradDnrFields[3][ind03];
      val4  = ptrCoefs[ indCoef   ]*vectOfGradDnrFields[4][ind00];
      val4 += ptrCoefs[ indCoef +1]*vectOfGradDnrFields[4][ind01];
      val4 += ptrCoefs[ indCoef +2]*vectOfGradDnrFields[4][ind02];
      val4 += ptrCoefs[ indCoef +3]*vectOfGradDnrFields[4][ind03];
      val5  = ptrCoefs[ indCoef   ]*vectOfGradDnrFields[5][ind00];
      val5 += ptrCoefs[ indCoef +1]*vectOfGradDnrFields[5][ind01];
      val5 += ptrCoefs[ indCoef +2]*vectOfGradDnrFields[5][ind02];
      val5 += ptrCoefs[ indCoef +3]*vectOfGradDnrFields[5][ind03];

      vectOfGradRcvFields[0][indR] = val0;
      vectOfGradRcvFields[1][indR] = val1;
      vectOfGradRcvFields[2][indR] = val2;
      vectOfGradRcvFields[3][indR] = val3;
      vectOfGradRcvFields[4][indR] = val4;
      vectOfGradRcvFields[5][indR] = val5;
      indCoef  += sizecoefs;
    }
    break;
      
  case 5: // Lagrange O5
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      indR  = rcvPts[noind];
      indD0 = donorPts[noind];  //car type 0 est toujour traité en dernier. Sinon noind pas valable
      k     = indD0/imdjmd;
      j     = (indD0-k*imdjmd)/imd;
      i     = (indD0-j*imd-k*imdjmd);
      val0=0.; val1=0.; val2=0.; val3=0.; val4=0.; val5=0.;
      for (E_Int kk=0; kk<5; kk++)
        for (E_Int jj=0; jj<5; jj++)
          for (E_Int ii=0; ii<5; ii++)
          {
            indD = (i+ii)+(j+jj)*imd+(k+kk)*imdjmd;
            val0 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfGradDnrFields[0][indD];               
            val1 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfGradDnrFields[1][indD];               
            val2 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfGradDnrFields[2][indD];               
            val3 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfGradDnrFields[3][indD];               
            val4 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfGradDnrFields[4][indD];               
            val5 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfGradDnrFields[5][indD];               
          }
      vectOfGradRcvFields[0][indR] = val0;
      vectOfGradRcvFields[1][indR] = val1;
      vectOfGradRcvFields[2][indR] = val2;
      vectOfGradRcvFields[3][indR] = val3;
      vectOfGradRcvFields[4][indR] = val4;
      vectOfGradRcvFields[5][indR] = val5;
      indCoef  += 15;
    }
    break;
      
  default: ;
}
