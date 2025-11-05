E_Int ind000,ind100,ind010,ind110, ind001,ind101,ind011,ind111;
E_Int ind00 ,ind10 ,ind01 ,ind11, ind02,ind03;
E_Float val0, val1, val2, val3, val4, val5;
switch (type)
{
  case 1:
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      indD0  = donorPts[noind];  
    
      vectOfRcvFields[0][noind] = vectOfDnrFields[0][indD0];
      vectOfRcvFields[1][noind] = vectOfDnrFields[1][indD0];
      vectOfRcvFields[2][noind] = vectOfDnrFields[2][indD0];
      vectOfRcvFields[3][noind] = vectOfDnrFields[3][indD0];
      vectOfRcvFields[4][noind] = vectOfDnrFields[4][indD0];
      if (nvars_loc==6){vectOfRcvFields[5][noind] = vectOfDnrFields[5][indD0];}
    }
    break;
  case 2: // Structure Lineaire O2 par tetra
// #ifdef _OPENMP4
//     #pragma omp simd
// #endif
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      ind000 = donorPts[noind];

      ind100 = ind000 +1; 
      ind010 = ind000 +imd; 
      ind110 = ind100 +imd; 
      ind001 = ind000+imdjmd;
      ind101 = ind100+imdjmd;
      ind011 = ind010+imdjmd;
      ind111 = ind110+imdjmd;
    
      val0  = ptrCoefs[ indCoef     ]*vectOfDnrFields[0][ind000];
      val0 += ptrCoefs[ indCoef + 1 ]*vectOfDnrFields[0][ind100];
      val0 += ptrCoefs[ indCoef + 2 ]*vectOfDnrFields[0][ind010];
      val0 += ptrCoefs[ indCoef + 3 ]*vectOfDnrFields[0][ind110];
      val0 += ptrCoefs[ indCoef + 4 ]*vectOfDnrFields[0][ind001];
      val0 += ptrCoefs[ indCoef + 5 ]*vectOfDnrFields[0][ind101];
      val0 += ptrCoefs[ indCoef + 6 ]*vectOfDnrFields[0][ind011];
      val0 += ptrCoefs[ indCoef + 7 ]*vectOfDnrFields[0][ind111];

      val1  = ptrCoefs[ indCoef     ]*vectOfDnrFields[1][ind000];
      val1 += ptrCoefs[ indCoef + 1 ]*vectOfDnrFields[1][ind100];
      val1 += ptrCoefs[ indCoef + 2 ]*vectOfDnrFields[1][ind010];
      val1 += ptrCoefs[ indCoef + 3 ]*vectOfDnrFields[1][ind110];
      val1 += ptrCoefs[ indCoef + 4 ]*vectOfDnrFields[1][ind001];
      val1 += ptrCoefs[ indCoef + 5 ]*vectOfDnrFields[1][ind101];
      val1 += ptrCoefs[ indCoef + 6 ]*vectOfDnrFields[1][ind011];
      val1 += ptrCoefs[ indCoef + 7 ]*vectOfDnrFields[1][ind111];

      val2  = ptrCoefs[ indCoef     ]*vectOfDnrFields[2][ind000];
      val2 += ptrCoefs[ indCoef + 1 ]*vectOfDnrFields[2][ind100];
      val2 += ptrCoefs[ indCoef + 2 ]*vectOfDnrFields[2][ind010];
      val2 += ptrCoefs[ indCoef + 3 ]*vectOfDnrFields[2][ind110];
      val2 += ptrCoefs[ indCoef + 4 ]*vectOfDnrFields[2][ind001];
      val2 += ptrCoefs[ indCoef + 5 ]*vectOfDnrFields[2][ind101];
      val2 += ptrCoefs[ indCoef + 6 ]*vectOfDnrFields[2][ind011];
      val2 += ptrCoefs[ indCoef + 7 ]*vectOfDnrFields[2][ind111];

      val3  = ptrCoefs[ indCoef     ]*vectOfDnrFields[3][ind000];
      val3 += ptrCoefs[ indCoef + 1 ]*vectOfDnrFields[3][ind100];
      val3 += ptrCoefs[ indCoef + 2 ]*vectOfDnrFields[3][ind010];
      val3 += ptrCoefs[ indCoef + 3 ]*vectOfDnrFields[3][ind110];
      val3 += ptrCoefs[ indCoef + 4 ]*vectOfDnrFields[3][ind001];
      val3 += ptrCoefs[ indCoef + 5 ]*vectOfDnrFields[3][ind101];
      val3 += ptrCoefs[ indCoef + 6 ]*vectOfDnrFields[3][ind011];
      val3 += ptrCoefs[ indCoef + 7 ]*vectOfDnrFields[3][ind111];

      val4  = ptrCoefs[ indCoef     ]*vectOfDnrFields[4][ind000];
      val4 += ptrCoefs[ indCoef + 1 ]*vectOfDnrFields[4][ind100];
      val4 += ptrCoefs[ indCoef + 2 ]*vectOfDnrFields[4][ind010];
      val4 += ptrCoefs[ indCoef + 3 ]*vectOfDnrFields[4][ind110];
      val4 += ptrCoefs[ indCoef + 4 ]*vectOfDnrFields[4][ind001];
      val4 += ptrCoefs[ indCoef + 5 ]*vectOfDnrFields[4][ind101];
      val4 += ptrCoefs[ indCoef + 6 ]*vectOfDnrFields[4][ind011];
      val4 += ptrCoefs[ indCoef + 7 ]*vectOfDnrFields[4][ind111];

      vectOfRcvFields[0][noind] = val0; //Density Pnt2
      vectOfRcvFields[1][noind] = val1; //VelocityX Pnt2
      vectOfRcvFields[2][noind] = val2; //VelocityY Pnt2
      vectOfRcvFields[3][noind] = val3; //VelocityZ Pnt2
      vectOfRcvFields[4][noind] = val4; //Temperature Pnt2

      if (nvars_loc==6)
      {
	      val5  = ptrCoefs[ indCoef     ]*vectOfDnrFields[5][ind000];
	      val5 += ptrCoefs[ indCoef + 1 ]*vectOfDnrFields[5][ind100];
	      val5 += ptrCoefs[ indCoef + 2 ]*vectOfDnrFields[5][ind010];
	      val5 += ptrCoefs[ indCoef + 3 ]*vectOfDnrFields[5][ind110];
	      val5 += ptrCoefs[ indCoef + 4 ]*vectOfDnrFields[5][ind001];
	      val5 += ptrCoefs[ indCoef + 5 ]*vectOfDnrFields[5][ind101];
	      val5 += ptrCoefs[ indCoef + 6 ]*vectOfDnrFields[5][ind011];
	      val5 += ptrCoefs[ indCoef + 7 ]*vectOfDnrFields[5][ind111];

	      vectOfRcvFields[5][noind] = val5; //TurbulentSANuTilde Pnt2
      }
      
      indCoef  += 8;
    }
    break;
    
  case 22:// O2CF 2D
// #ifdef _OPENMP4
//     #pragma omp simd
// #endif
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      ind00 = donorPts[noind];
      ind10 = ind00 +1; 
      ind01 = ind00 +imd; 
      ind11 = ind10 +imd; 
      val0  = ptrCoefs[ indCoef     ]*vectOfDnrFields[0][ind00];
      val0 += ptrCoefs[ indCoef + 1 ]*vectOfDnrFields[0][ind10];
      val0 += ptrCoefs[ indCoef + 2 ]*vectOfDnrFields[0][ind01];
      val0 += ptrCoefs[ indCoef + 3 ]*vectOfDnrFields[0][ind11];

      val1  = ptrCoefs[ indCoef     ]*vectOfDnrFields[1][ind00];
      val1 += ptrCoefs[ indCoef + 1 ]*vectOfDnrFields[1][ind10];
      val1 += ptrCoefs[ indCoef + 2 ]*vectOfDnrFields[1][ind01];
      val1 += ptrCoefs[ indCoef + 3 ]*vectOfDnrFields[1][ind11];

      val2  = ptrCoefs[ indCoef     ]*vectOfDnrFields[2][ind00];
      val2 += ptrCoefs[ indCoef + 1 ]*vectOfDnrFields[2][ind10];
      val2 += ptrCoefs[ indCoef + 2 ]*vectOfDnrFields[2][ind01];
      val2 += ptrCoefs[ indCoef + 3 ]*vectOfDnrFields[2][ind11];

      val3  = ptrCoefs[ indCoef     ]*vectOfDnrFields[3][ind00];
      val3 += ptrCoefs[ indCoef + 1 ]*vectOfDnrFields[3][ind10];
      val3 += ptrCoefs[ indCoef + 2 ]*vectOfDnrFields[3][ind01];
      val3 += ptrCoefs[ indCoef + 3 ]*vectOfDnrFields[3][ind11];

      val4  = ptrCoefs[ indCoef     ]*vectOfDnrFields[4][ind00];
      val4 += ptrCoefs[ indCoef + 1 ]*vectOfDnrFields[4][ind10];
      val4 += ptrCoefs[ indCoef + 2 ]*vectOfDnrFields[4][ind01];
      val4 += ptrCoefs[ indCoef + 3 ]*vectOfDnrFields[4][ind11];

      vectOfRcvFields[0][noind] = val0; //Density Pnt2d
      vectOfRcvFields[1][noind] = val1; //VelocityX Pnt2
      vectOfRcvFields[2][noind] = val2; //VelocityY Pnt2
      vectOfRcvFields[3][noind] = val3; //VelocityZ Pnt2
      vectOfRcvFields[4][noind] = val4; //Temperature Pnt2
      
      if (nvars_loc==6)
      {
	      val5  = ptrCoefs[ indCoef     ]*vectOfDnrFields[5][ind00];
	      val5 += ptrCoefs[ indCoef + 1 ]*vectOfDnrFields[5][ind10];
	      val5 += ptrCoefs[ indCoef + 2 ]*vectOfDnrFields[5][ind01];
	      val5 += ptrCoefs[ indCoef + 3 ]*vectOfDnrFields[5][ind11];
	      vectOfRcvFields[5][noind] = val5; //TurbulentSANuTilde Pnt2
      }
      indCoef += 4;
    }
    break;      
  default: ;
}
