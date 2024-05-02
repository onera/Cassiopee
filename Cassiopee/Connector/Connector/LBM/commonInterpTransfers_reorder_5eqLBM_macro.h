//E_Int ind000,ind100,ind010,ind110, ind001,ind101,ind011,ind111;
//E_Int ind00 ,ind10 ,ind01 ,ind11, ind02,ind03;
//E_Float val0, val1, val2, val3, val4;
E_Int vel_temp;
E_Float valtmp;
switch (type)
{
  case 0:  //  nuage de pts quelconque
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      indR   = rcvPts[noind];
      ncfLoc = donorPts[noi];// nb de pts pour la formule
      val0=0.; val1=0.; val2=0.; val3=0.; val4=0.; 
      for (E_Int kk = 1; kk <= ncfLoc; kk++)
      {
        indD0         = donorPts[noi+kk];
        val0 += ptrCoefs[ indCoef + kk-1]*vectOfmacroDnrFields[0][indD0];
        val1 += ptrCoefs[ indCoef + kk-1]*vectOfmacroDnrFields[1][indD0];
        val2 += ptrCoefs[ indCoef + kk-1]*vectOfmacroDnrFields[2][indD0];
        val3 += ptrCoefs[ indCoef + kk-1]*vectOfmacroDnrFields[3][indD0];
        val4 += ptrCoefs[ indCoef + kk-1]*vectOfmacroDnrFields[4][indD0];
      }
      vectOfmacroRcvFields[0][indR] = val0;
      vectOfmacroRcvFields[1][indR] = val1;
      vectOfmacroRcvFields[2][indR] = val2;
      vectOfmacroRcvFields[3][indR] = val3;
      vectOfmacroRcvFields[4][indR] = val4;
      sizecoefs = ncfLoc;
      noi      += ncfLoc+1;
      indCoef  += sizecoefs;
    }
    break;

  case 1:
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      indR   = rcvPts[noind];
      indD0  = donorPts[noind];  //car type 0 est toujour traité en dernier. Sinon noind pas valable
    
      vectOfmacroRcvFields[0][indR] = vectOfmacroDnrFields[0][indD0];
      vectOfmacroRcvFields[1][indR] = vectOfmacroDnrFields[1][indD0];
      vectOfmacroRcvFields[2][indR] = vectOfmacroDnrFields[2][indD0];
      vectOfmacroRcvFields[3][indR] = vectOfmacroDnrFields[3][indD0];
      vectOfmacroRcvFields[4][indR] = vectOfmacroDnrFields[4][indD0];
    }
    break;
    
  case 2: // Structure Lineaire O2 par tetra
// #ifdef _OPENMP4
//     #pragma omp simd
// #endif
    //printf("nstep=%d\n",nstep);
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
      if(nstep ==1){
	vel_temp = 0;
#           include "nstep_1_macro.h"
	val0 = valtmp;

	vel_temp = 1;
#           include "nstep_1_macro.h"
	val1 = valtmp;

	vel_temp = 2;
#           include "nstep_1_macro.h"
	val2 = valtmp;

	vel_temp = 3;
#           include "nstep_1_macro.h"
	val3 = valtmp;

	vel_temp = 4;
#           include "nstep_1_macro.h"
	val4 = valtmp;
      }
      else if (nstep ==2){
	vel_temp = 0;
#           include "nstep_2_macro.h"
	val0 = valtmp;

	vel_temp = 1;
#           include "nstep_2_macro.h"
	val1 = valtmp;

	vel_temp = 2;
#           include "nstep_2_macro.h"
	val2 = valtmp;

	vel_temp = 3;
#           include "nstep_2_macro.h"
	val3 = valtmp;

	vel_temp = 4;
#           include "nstep_2_macro.h"
	val4 = valtmp;
      }
      
      vectOfmacroRcvFields[0][indR] = val0;
      vectOfmacroRcvFields[1][indR] = val1*val0;
      vectOfmacroRcvFields[2][indR] = val2*val0;
      vectOfmacroRcvFields[3][indR] = val3*val0;
      vectOfmacroRcvFields[4][indR] = val4*val0;

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
      val0  = ptrCoefs[ indCoef     ]*vectOfmacroDnrFields[0][ind00];
      val0 += ptrCoefs[ indCoef + 1 ]*vectOfmacroDnrFields[0][ind10];
      val0 += ptrCoefs[ indCoef + 2 ]*vectOfmacroDnrFields[0][ind01];
      val0 += ptrCoefs[ indCoef + 3 ]*vectOfmacroDnrFields[0][ind11];

      val1  = ptrCoefs[ indCoef     ]*vectOfmacroDnrFields[1][ind00];
      val1 += ptrCoefs[ indCoef + 1 ]*vectOfmacroDnrFields[1][ind10];
      val1 += ptrCoefs[ indCoef + 2 ]*vectOfmacroDnrFields[1][ind01];
      val1 += ptrCoefs[ indCoef + 3 ]*vectOfmacroDnrFields[1][ind11];

      val2  = ptrCoefs[ indCoef     ]*vectOfmacroDnrFields[2][ind00];
      val2 += ptrCoefs[ indCoef + 1 ]*vectOfmacroDnrFields[2][ind10];
      val2 += ptrCoefs[ indCoef + 2 ]*vectOfmacroDnrFields[2][ind01];
      val2 += ptrCoefs[ indCoef + 3 ]*vectOfmacroDnrFields[2][ind11];

      val3  = ptrCoefs[ indCoef     ]*vectOfmacroDnrFields[3][ind00];
      val3 += ptrCoefs[ indCoef + 1 ]*vectOfmacroDnrFields[3][ind10];
      val3 += ptrCoefs[ indCoef + 2 ]*vectOfmacroDnrFields[3][ind01];
      val3 += ptrCoefs[ indCoef + 3 ]*vectOfmacroDnrFields[3][ind11];

      val4  = ptrCoefs[ indCoef     ]*vectOfmacroDnrFields[4][ind00];
      val4 += ptrCoefs[ indCoef + 1 ]*vectOfmacroDnrFields[4][ind10];
      val4 += ptrCoefs[ indCoef + 2 ]*vectOfmacroDnrFields[4][ind01];
      val4 += ptrCoefs[ indCoef + 3 ]*vectOfmacroDnrFields[4][ind11];

      vectOfmacroRcvFields[0][indR] = val0;
      vectOfmacroRcvFields[1][indR] = val1;
      vectOfmacroRcvFields[2][indR] = val2;
      vectOfmacroRcvFields[3][indR] = val3;
      vectOfmacroRcvFields[4][indR] = val4;
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
      val0=0.; val1=0.; val2=0.; val3=0.; val4=0.; 

      for (E_Int kk=0; kk<3; kk++)
        for (E_Int jj=0; jj<3; jj++)
          for (E_Int ii=0; ii<3; ii++)
          {
            indD = (i+ii)+(j+jj)*imd+(k+kk)*imdjmd;
            val0 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfmacroDnrFields[0][indD];               
            val1 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfmacroDnrFields[1][indD];               
            val2 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfmacroDnrFields[2][indD];               
            val3 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfmacroDnrFields[3][indD];               
            val4 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfmacroDnrFields[4][indD];               
          }
      vectOfmacroRcvFields[0][indR] = val0;
      vectOfmacroRcvFields[1][indR] = val1;
      vectOfmacroRcvFields[2][indR] = val2;
      vectOfmacroRcvFields[3][indR] = val3;
      vectOfmacroRcvFields[4][indR] = val4;
      noi      += 1;
      indCoef  += sizecoefs;
    }
    break;
      
  case 4: // Tetra O2
// #ifdef _OPENMP4
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

      val0  = ptrCoefs[ indCoef   ]*vectOfmacroDnrFields[0][ind00];
      val0 += ptrCoefs[ indCoef +1]*vectOfmacroDnrFields[0][ind01];
      val0 += ptrCoefs[ indCoef +2]*vectOfmacroDnrFields[0][ind02];
      val0 += ptrCoefs[ indCoef +3]*vectOfmacroDnrFields[0][ind03];
      val1  = ptrCoefs[ indCoef   ]*vectOfmacroDnrFields[1][ind00];
      val1 += ptrCoefs[ indCoef +1]*vectOfmacroDnrFields[1][ind01];
      val1 += ptrCoefs[ indCoef +2]*vectOfmacroDnrFields[1][ind02];
      val1 += ptrCoefs[ indCoef +3]*vectOfmacroDnrFields[1][ind03];
      val2  = ptrCoefs[ indCoef   ]*vectOfmacroDnrFields[2][ind00];
      val2 += ptrCoefs[ indCoef +1]*vectOfmacroDnrFields[2][ind01];
      val2 += ptrCoefs[ indCoef +2]*vectOfmacroDnrFields[2][ind02];
      val2 += ptrCoefs[ indCoef +3]*vectOfmacroDnrFields[2][ind03];
      val3  = ptrCoefs[ indCoef   ]*vectOfmacroDnrFields[3][ind00];
      val3 += ptrCoefs[ indCoef +1]*vectOfmacroDnrFields[3][ind01];
      val3 += ptrCoefs[ indCoef +2]*vectOfmacroDnrFields[3][ind02];
      val3 += ptrCoefs[ indCoef +3]*vectOfmacroDnrFields[3][ind03];
      val4  = ptrCoefs[ indCoef   ]*vectOfmacroDnrFields[4][ind00];
      val4 += ptrCoefs[ indCoef +1]*vectOfmacroDnrFields[4][ind01];
      val4 += ptrCoefs[ indCoef +2]*vectOfmacroDnrFields[4][ind02];
      val4 += ptrCoefs[ indCoef +3]*vectOfmacroDnrFields[4][ind03];
      vectOfmacroRcvFields[0][indR] = val0;
      vectOfmacroRcvFields[1][indR] = val1;
      vectOfmacroRcvFields[2][indR] = val2;
      vectOfmacroRcvFields[3][indR] = val3;
      vectOfmacroRcvFields[4][indR] = val4;
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
      val0=0.; val1=0.; val2=0.; val3=0.; val4=0.; 
      for (E_Int kk=0; kk<5; kk++)
        for (E_Int jj=0; jj<5; jj++)
          for (E_Int ii=0; ii<5; ii++)
          {
            indD = (i+ii)+(j+jj)*imd+(k+kk)*imdjmd;
            val0 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfmacroDnrFields[0][indD];               
            val1 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfmacroDnrFields[1][indD];               
            val2 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfmacroDnrFields[2][indD];               
            val3 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfmacroDnrFields[3][indD];               
            val4 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfmacroDnrFields[4][indD];               
          }
      vectOfmacroRcvFields[0][indR] = val0;
      vectOfmacroRcvFields[1][indR] = val1;
      vectOfmacroRcvFields[2][indR] = val2;
      vectOfmacroRcvFields[3][indR] = val3;
      vectOfmacroRcvFields[4][indR] = val4;
      indCoef  += 15;
    }
    break;
      
  default: ;
}
