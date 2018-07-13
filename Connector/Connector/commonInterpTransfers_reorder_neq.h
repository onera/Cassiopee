E_Int ind000,ind100,ind010,ind110, ind001,ind101,ind011,ind111;
E_Int ind00 ,ind10 ,ind01 ,ind11, ind02,ind03;
E_Float val;
switch (type)
{
  case 0:  //  nuage de pts quelconque
   for (E_Int ne    = 0     ; ne    < nvars_loc ; ne++)
    {
     indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
     noi       = shiftDonor;
     for (E_Int noind = pt_deb; noind < pt_fin; noind++)
     {
      indR   = rcvPts[noind];
      ncfLoc = donorPts[noi];// nb de pts pour la formule

      val = 0.;
      for (E_Int kk = 1; kk <= ncfLoc; kk++)
      {
        indD0         = donorPts[noi+kk];
        val += ptrCoefs[ indCoef + kk-1]*vectOfDnrFields[ne][indD0];
      }
      vectOfRcvFields[ne][indR] = val;
      sizecoefs = ncfLoc;
      noi      += ncfLoc+1;
      indCoef  += sizecoefs;
     } //noind
    } //ne
    break;

  case 1:
    for (E_Int ne    = 0     ; ne    < nvars_loc ; ne++)
    {
     for (E_Int noind = pt_deb; noind < pt_fin; noind++)
     {
      indR   = rcvPts[noind];
      indD0  = donorPts[noind];  //car type 0 est toujour traite en dernier. Sinon noind pas valable
    
      vectOfRcvFields[ne][indR] = vectOfDnrFields[ne][indD0];
     }
    }
    break;
    
  case 2: // Structure Lineaire O2 par tetra
   for (E_Int ne    = 0     ; ne    < nvars_loc ; ne++)
    {
     indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
// #ifdef _OPENMP4
//      #pragma omp simd
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
    
      val  = ptrCoefs[ indCoef     ]*vectOfDnrFields[ne][ind000];
      val += ptrCoefs[ indCoef + 1 ]*vectOfDnrFields[ne][ind100];
      val += ptrCoefs[ indCoef + 2 ]*vectOfDnrFields[ne][ind010];
      val += ptrCoefs[ indCoef + 3 ]*vectOfDnrFields[ne][ind110];
      val += ptrCoefs[ indCoef + 4 ]*vectOfDnrFields[ne][ind001];
      val += ptrCoefs[ indCoef + 5 ]*vectOfDnrFields[ne][ind101];
      val += ptrCoefs[ indCoef + 6 ]*vectOfDnrFields[ne][ind011];
      val += ptrCoefs[ indCoef + 7 ]*vectOfDnrFields[ne][ind111];

      vectOfRcvFields[ne][indR] = val;
      indCoef  += 8;
     }
    }
    break;
    
  case 22:// O2CF 2D
    for (E_Int ne    = 0     ; ne    < nvars_loc ; ne++)
    {
    indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
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

      val  = ptrCoefs[ indCoef     ]*vectOfDnrFields[ne][ind00];
      val += ptrCoefs[ indCoef + 1 ]*vectOfDnrFields[ne][ind10];
      val += ptrCoefs[ indCoef + 2 ]*vectOfDnrFields[ne][ind01];
      val += ptrCoefs[ indCoef + 3 ]*vectOfDnrFields[ne][ind11];
      vectOfRcvFields[ne][indR] = val;

      indCoef  += 4;
    }
   }
    break;

  case 3: // Lagrange O3
   for (E_Int ne    = 0     ; ne    < nvars_loc ; ne++)
   { 
   indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
     for (E_Int noind = pt_deb; noind < pt_fin; noind++)
     {
      indR  = rcvPts[noind];
      indD0 = donorPts[noind];  //car type 0 est toujour traite en dernier. Sinon noind pas valable
      k     = indD0/imdjmd;
      j     = (indD0-k*imdjmd)/imd;
      i     = (indD0-j*imd-k*imdjmd);
    
      val=0.;
      for (E_Int kk=0; kk<3; kk++)
        for (E_Int jj=0; jj<3; jj++)
          for (E_Int ii=0; ii<3; ii++)
          {
            indD = (i+ii)+(j+jj)*imd+(k+kk)*imdjmd;
            val += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfDnrFields[ne][indD];
          }
      vectOfRcvFields[ne][indR] = val;

      noi      += 1;
      indCoef  += sizecoefs;
     }
   }
   break;
      
  case 4: // Tetra O2
   for (E_Int ne    = 0     ; ne    < nvars_loc ; ne++)
   {
     indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
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

      val  = ptrCoefs[ indCoef   ]*vectOfDnrFields[ne][ind00];
      val += ptrCoefs[ indCoef +1]*vectOfDnrFields[ne][ind01];
      val += ptrCoefs[ indCoef +2]*vectOfDnrFields[ne][ind02];
      val += ptrCoefs[ indCoef +3]*vectOfDnrFields[ne][ind03];
      vectOfRcvFields[ne][indR] = val;
      indCoef  += sizecoefs;
    }
   }
   break;
      
  case 5: // Lagrange O5
   for (E_Int ne    = 0     ; ne    < nvars_loc ; ne++)
   {
    indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      indR  = rcvPts[noind];
      indD0 = donorPts[noind];  //car type 0 est toujour traite en dernier. Sinon noind pas valable
      k     = indD0/imdjmd;
      j     = (indD0-k*imdjmd)/imd;
      i     = (indD0-j*imd-k*imdjmd);
      val=0.;
      for (E_Int kk=0; kk<5; kk++)
        for (E_Int jj=0; jj<5; jj++)
          for (E_Int ii=0; ii<5; ii++)
          {
            indD = (i+ii)+(j+jj)*imd+(k+kk)*imdjmd;
            val += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfDnrFields[ne][indD];
          }
      vectOfRcvFields[ne][indR] = val;
      indCoef  += 15;
    }
   }
   break;
      
  default: ;
}
