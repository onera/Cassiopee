E_Int ind000,ind100,ind010,ind110, ind001,ind101,ind011,ind111;
E_Int ind00 ,ind10 ,ind01 ,ind11, ind02,ind03;
E_Float val0, val1, val2, val3, val4, val5, val6, val7, val8, val9;
E_Float val10, val11, val12, val13, val14, val15, val16, val17, val18;

switch (type)
{
  case 0:  //  nuage de pts quelconque
    for (E_Int noind = pt_deb; noind < pt_fin; noind++)
    {
      indR   = rcvPts[noind];
      ncfLoc = donorPts[noi];// nb de pts pour la formule
      val0=0.; val1=0.; val2=0.; val3=0.; val4=0.; val5=0.; val6=0.; val7=0.; val8=0.; val9=0.;
      val10=0.; val11=1.; val12=0.; val13=0.; val14=0.; val15=0.; val16=0.; val17=0.; val18=0.;
      for (E_Int kk = 1; kk <= ncfLoc; kk++)
      {
        indD0         = donorPts[noi+kk];
        val0  += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[0][indD0];
        val1  += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[1][indD0];
        val2  += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[2][indD0];
        val3  += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[3][indD0];
        val4  += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[4][indD0];
        val5  += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[5][indD0];
        val6  += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[6][indD0];
        val7  += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[7][indD0];
        val8  += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[8][indD0];
        val9  += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[9][indD0];
        val10 += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[10][indD0];
        val11 += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[11][indD0];
        val12 += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[12][indD0];
        val13 += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[13][indD0];
        val14 += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[14][indD0];
        val15 += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[15][indD0];
        val16 += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[16][indD0];
        val17 += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[17][indD0];
        val18 += ptrCoefs[ indCoef + kk-1]*vectOfQneqDnrFields[18][indD0];
      }
      vectOfQneqRcvFields[0][indR] = val0;
      vectOfQneqRcvFields[1][indR] = val1;
      vectOfQneqRcvFields[2][indR] = val2;
      vectOfQneqRcvFields[3][indR] = val3;
      vectOfQneqRcvFields[4][indR] = val4;
      vectOfQneqRcvFields[5][indR] = val5;
      vectOfQneqRcvFields[6][indR] = val6;
      vectOfQneqRcvFields[7][indR] = val7;
      vectOfQneqRcvFields[8][indR] = val8;
      vectOfQneqRcvFields[9][indR] = val9;
      vectOfQneqRcvFields[10][indR]= val10;
      vectOfQneqRcvFields[11][indR]= val11;
      vectOfQneqRcvFields[12][indR]= val12;
      vectOfQneqRcvFields[13][indR]= val13;
      vectOfQneqRcvFields[14][indR]= val14;
      vectOfQneqRcvFields[15][indR]= val15;
      vectOfQneqRcvFields[16][indR]= val16;
      vectOfQneqRcvFields[17][indR]= val17;
      vectOfQneqRcvFields[18][indR]= val18;


//      //macros
//#           include "commonInterpTransfers_reorder_macro.h"

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
    
      vectOfQneqRcvFields[0][indR] = vectOfQneqDnrFields[0][indD0];
      vectOfQneqRcvFields[1][indR] = vectOfQneqDnrFields[1][indD0];
      vectOfQneqRcvFields[2][indR] = vectOfQneqDnrFields[2][indD0];
      vectOfQneqRcvFields[3][indR] = vectOfQneqDnrFields[3][indD0];
      vectOfQneqRcvFields[4][indR] = vectOfQneqDnrFields[4][indD0];
      vectOfQneqRcvFields[5][indR] = vectOfQneqDnrFields[5][indD0];
      vectOfQneqRcvFields[6][indR] = vectOfQneqDnrFields[6][indD0];
      vectOfQneqRcvFields[7][indR] = vectOfQneqDnrFields[7][indD0];
      vectOfQneqRcvFields[8][indR] = vectOfQneqDnrFields[8][indD0];
      vectOfQneqRcvFields[9][indR] = vectOfQneqDnrFields[9][indD0];
      vectOfQneqRcvFields[10][indR]= vectOfQneqDnrFields[10][indD0];
      vectOfQneqRcvFields[11][indR]= vectOfQneqDnrFields[11][indD0];
      vectOfQneqRcvFields[12][indR]= vectOfQneqDnrFields[12][indD0];
      vectOfQneqRcvFields[13][indR]= vectOfQneqDnrFields[13][indD0];
      vectOfQneqRcvFields[14][indR]= vectOfQneqDnrFields[14][indD0];
      vectOfQneqRcvFields[15][indR]= vectOfQneqDnrFields[15][indD0];
      vectOfQneqRcvFields[16][indR]= vectOfQneqDnrFields[16][indD0];
      vectOfQneqRcvFields[17][indR]= vectOfQneqDnrFields[17][indD0];
      vectOfQneqRcvFields[18][indR]= vectOfQneqDnrFields[18][indD0];
//      //macros
//#           include "commonInterpTransfers_reorder_macro.h"

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

      val0  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[0][ind000];
      val0 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[0][ind100];
      val0 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[0][ind010];
      val0 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[0][ind110];
      val0 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[0][ind001];
      val0 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[0][ind101];
      val0 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[0][ind011];
      val0 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[0][ind111];

      val1  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[1][ind000];
      val1 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[1][ind100];
      val1 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[1][ind010];
      val1 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[1][ind110];
      val1 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[1][ind001];
      val1 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[1][ind101];
      val1 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[1][ind011];
      val1 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[1][ind111];

      val2  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[2][ind000];
      val2 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[2][ind100];
      val2 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[2][ind010];
      val2 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[2][ind110];
      val2 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[2][ind001];
      val2 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[2][ind101];
      val2 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[2][ind011];
      val2 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[2][ind111];

      val3  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[3][ind000];
      val3 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[3][ind100];
      val3 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[3][ind010];
      val3 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[3][ind110];
      val3 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[3][ind001];
      val3 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[3][ind101];
      val3 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[3][ind011];
      val3 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[3][ind111];

      val4  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[4][ind000];
      val4 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[4][ind100];
      val4 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[4][ind010];
      val4 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[4][ind110];
      val4 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[4][ind001];
      val4 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[4][ind101];
      val4 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[4][ind011];
      val4 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[4][ind111];

      val5  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[5][ind000];
      val5 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[5][ind100];
      val5 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[5][ind010];
      val5 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[5][ind110];
      val5 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[5][ind001];
      val5 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[5][ind101];
      val5 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[5][ind011];
      val5 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[5][ind111];

      val6  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[6][ind000];
      val6 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[6][ind100];
      val6 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[6][ind010];
      val6 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[6][ind110];
      val6 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[6][ind001];
      val6 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[6][ind101];
      val6 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[6][ind011];
      val6 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[6][ind111];

      val7  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[7][ind000];
      val7 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[7][ind100];
      val7 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[7][ind010];
      val7 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[7][ind110];
      val7 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[7][ind001];
      val7 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[7][ind101];
      val7 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[7][ind011];
      val7 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[7][ind111];

      val8  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[8][ind000];
      val8 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[8][ind100];
      val8 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[8][ind010];
      val8 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[8][ind110];
      val8 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[8][ind001];
      val8 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[8][ind101];
      val8 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[8][ind011];
      val8 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[8][ind111];

      val9  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[9][ind000];
      val9 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[9][ind100];
      val9 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[9][ind010];
      val9 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[9][ind110];
      val9 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[9][ind001];
      val9 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[9][ind101];
      val9 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[9][ind011];
      val9 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[9][ind111];

      val10  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[10][ind000];
      val10 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[10][ind100];
      val10 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[10][ind010];
      val10 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[10][ind110];
      val10 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[10][ind001];
      val10 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[10][ind101];
      val10 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[10][ind011];
      val10 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[10][ind111];

      val11  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[11][ind000];
      val11 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[11][ind100];
      val11 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[11][ind010];
      val11 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[11][ind110];
      val11 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[11][ind001];
      val11 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[11][ind101];
      val11 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[11][ind011];
      val11 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[11][ind111];

      val12  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[12][ind000];
      val12 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[12][ind100];
      val12 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[12][ind010];
      val12 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[12][ind110];
      val12 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[12][ind001];
      val12 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[12][ind101];
      val12 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[12][ind011];
      val12 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[12][ind111];

      val13  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[13][ind000];
      val13 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[13][ind100];
      val13 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[13][ind010];
      val13 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[13][ind110];
      val13 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[13][ind001];
      val13 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[13][ind101];
      val13 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[13][ind011];
      val13 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[13][ind111];

      val14  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[14][ind000];
      val14 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[14][ind100];
      val14 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[14][ind010];
      val14 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[14][ind110];
      val14 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[14][ind001];
      val14 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[14][ind101];
      val14 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[14][ind011];
      val14 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[14][ind111];

      val15  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[15][ind000];
      val15 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[15][ind100];
      val15 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[15][ind010];
      val15 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[15][ind110];
      val15 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[15][ind001];
      val15 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[15][ind101];
      val15 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[15][ind011];
      val15 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[15][ind111];

      val16  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[16][ind000];
      val16 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[16][ind100];
      val16 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[16][ind010];
      val16 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[16][ind110];
      val16 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[16][ind001];
      val16 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[16][ind101];
      val16 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[16][ind011];
      val16 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[16][ind111];

      val17  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[17][ind000];
      val17 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[17][ind100];
      val17 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[17][ind010];
      val17 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[17][ind110];
      val17 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[17][ind001];
      val17 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[17][ind101];
      val17 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[17][ind011];
      val17 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[17][ind111];

      val18  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[18][ind000];
      val18 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[18][ind100];
      val18 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[18][ind010];
      val18 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[18][ind110];
      val18 += ptrCoefs[ indCoef + 4 ]*vectOfQneqDnrFields[18][ind001];
      val18 += ptrCoefs[ indCoef + 5 ]*vectOfQneqDnrFields[18][ind101];
      val18 += ptrCoefs[ indCoef + 6 ]*vectOfQneqDnrFields[18][ind011];
      val18 += ptrCoefs[ indCoef + 7 ]*vectOfQneqDnrFields[18][ind111];


      vectOfQneqRcvFields[0][indR] = val0;
      vectOfQneqRcvFields[1][indR] = val1;
      vectOfQneqRcvFields[2][indR] = val2;
      vectOfQneqRcvFields[3][indR] = val3;
      vectOfQneqRcvFields[4][indR] = val4;
      vectOfQneqRcvFields[5][indR] = val5;
      vectOfQneqRcvFields[6][indR] = val6;
      vectOfQneqRcvFields[7][indR] = val7;
      vectOfQneqRcvFields[8][indR] = val8;
      vectOfQneqRcvFields[9][indR] = val9;
      vectOfQneqRcvFields[10][indR]= val10;
      vectOfQneqRcvFields[11][indR]= val11;
      vectOfQneqRcvFields[12][indR]= val12;
      vectOfQneqRcvFields[13][indR]= val13;
      vectOfQneqRcvFields[14][indR]= val14;
      vectOfQneqRcvFields[15][indR]= val15;
      vectOfQneqRcvFields[16][indR]= val16;
      vectOfQneqRcvFields[17][indR]= val17;
      vectOfQneqRcvFields[18][indR]= val18;
      
//      //macros
//#           include "commonInterpTransfers_reorder_macro.h"

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

      val0  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[0][ind00];
      val0 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[0][ind10];
      val0 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[0][ind01];
      val0 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[0][ind11];

      val1  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[1][ind00];
      val1 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[1][ind10];
      val1 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[1][ind01];
      val1 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[1][ind11];

      val2  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[2][ind00];
      val2 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[2][ind10];
      val2 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[2][ind01];
      val2 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[2][ind11];

      val3  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[3][ind00];
      val3 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[3][ind10];
      val3 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[3][ind01];
      val3 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[3][ind11];

      val4  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[4][ind00];
      val4 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[4][ind10];
      val4 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[4][ind01];
      val4 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[4][ind11];

      val5  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[5][ind00];
      val5 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[5][ind10];
      val5 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[5][ind01];
      val5 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[5][ind11];

      val6  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[6][ind00];
      val6 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[6][ind10];
      val6 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[6][ind01];
      val6 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[6][ind11];

      val7  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[7][ind00];
      val7 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[7][ind10];
      val7 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[7][ind01];
      val7 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[7][ind11];

      val8  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[8][ind00];
      val8 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[8][ind10];
      val8 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[8][ind01];
      val8 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[8][ind11];

      val9  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[9][ind00];
      val9 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[9][ind10];
      val9 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[9][ind01];
      val9 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[9][ind11];

      val10  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[10][ind00];
      val10 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[10][ind10];
      val10 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[10][ind01];
      val10 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[10][ind11];

      val11  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[11][ind00];
      val11 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[11][ind10];
      val11 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[11][ind01];
      val11 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[11][ind11];

      val12  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[12][ind00];
      val12 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[12][ind10];
      val12 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[12][ind01];
      val12 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[12][ind11];

      val13  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[13][ind00];
      val13 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[13][ind10];
      val13 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[13][ind01];
      val13 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[13][ind11];

      val14  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[14][ind00];
      val14 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[14][ind10];
      val14 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[14][ind01];
      val14 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[14][ind11];

      val15  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[15][ind00];
      val15 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[15][ind10];
      val15 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[15][ind01];
      val15 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[15][ind11];

      val16  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[16][ind00];
      val16 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[16][ind10];
      val16 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[16][ind01];
      val16 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[16][ind11];

      val17  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[17][ind00];
      val17 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[17][ind10];
      val17 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[17][ind01];
      val17 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[17][ind11];

      val18  = ptrCoefs[ indCoef     ]*vectOfQneqDnrFields[18][ind00];
      val18 += ptrCoefs[ indCoef + 1 ]*vectOfQneqDnrFields[18][ind10];
      val18 += ptrCoefs[ indCoef + 2 ]*vectOfQneqDnrFields[18][ind01];
      val18 += ptrCoefs[ indCoef + 3 ]*vectOfQneqDnrFields[18][ind11];

      vectOfQneqRcvFields[0][indR] = val0;
      vectOfQneqRcvFields[1][indR] = val1;
      vectOfQneqRcvFields[2][indR] = val2;
      vectOfQneqRcvFields[3][indR] = val3;
      vectOfQneqRcvFields[4][indR] = val4;
      vectOfQneqRcvFields[5][indR] = val5;
      vectOfQneqRcvFields[6][indR] = val6;
      vectOfQneqRcvFields[7][indR] = val7;
      vectOfQneqRcvFields[8][indR] = val8;
      vectOfQneqRcvFields[9][indR] = val9;
      vectOfQneqRcvFields[10][indR]= val10;
      vectOfQneqRcvFields[11][indR]= val11;
      vectOfQneqRcvFields[12][indR]= val12;
      vectOfQneqRcvFields[13][indR]= val13;
      vectOfQneqRcvFields[14][indR]= val14;
      vectOfQneqRcvFields[15][indR]= val15;
      vectOfQneqRcvFields[16][indR]= val16;
      vectOfQneqRcvFields[17][indR]= val17;
      vectOfQneqRcvFields[18][indR]= val18;

      
//      //macros
//#           include "commonInterpTransfers_reorder_macro.h"

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
      val0=0.; val1=0.; val2=0.; val3=0.; val4=0.; val5=0.; val6=0.; val7=0.; val8=0.; val9=0.;
      val10=0.; val11=1.; val12=0.; val13=0.; val14=0.; val15=0.; val16=0.; val17=0.; val18=0.;

      for (E_Int kk=0; kk<3; kk++)
        for (E_Int jj=0; jj<3; jj++)
          for (E_Int ii=0; ii<3; ii++)
          {
            indD = (i+ii)+(j+jj)*imd+(k+kk)*imdjmd;
            val0 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[0][indD];               
            val1 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[1][indD];               
            val2 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[2][indD];               
            val3 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[3][indD];               
            val4 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[4][indD];               
            val5 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[5][indD];               
            val6 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[6][indD];               
            val7 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[7][indD];               
            val8 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[8][indD];               
            val9 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[9][indD];               
            val10 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[10][indD];               
            val11 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[11][indD];               
            val12 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[12][indD];               
            val13 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[13][indD];               
            val14 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[14][indD];               
            val15 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[15][indD];               
            val16 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[16][indD];               
            val17 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[17][indD];               
            val18 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+3]*ptrCoefs[ indCoef + kk+6]*vectOfQneqDnrFields[18][indD];               
          }
      vectOfQneqRcvFields[0][indR] = val0;
      vectOfQneqRcvFields[1][indR] = val1;
      vectOfQneqRcvFields[2][indR] = val2;
      vectOfQneqRcvFields[3][indR] = val3;
      vectOfQneqRcvFields[4][indR] = val4;
      vectOfQneqRcvFields[5][indR] = val5;
      vectOfQneqRcvFields[6][indR] = val6;
      vectOfQneqRcvFields[7][indR] = val7;
      vectOfQneqRcvFields[8][indR] = val8;
      vectOfQneqRcvFields[9][indR] = val9;
      vectOfQneqRcvFields[10][indR]= val10;
      vectOfQneqRcvFields[11][indR]= val11;
      vectOfQneqRcvFields[12][indR]= val12;
      vectOfQneqRcvFields[13][indR]= val13;
      vectOfQneqRcvFields[14][indR]= val14;
      vectOfQneqRcvFields[15][indR]= val15;
      vectOfQneqRcvFields[16][indR]= val16;
      vectOfQneqRcvFields[17][indR]= val17;
      vectOfQneqRcvFields[18][indR]= val18;
//      //macros
//#           include "commonInterpTransfers_reorder_macro.h"

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
      val0  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[0][ind00];
      val0 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[0][ind01];
      val0 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[0][ind02];
      val0 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[0][ind03];
      val1  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[1][ind00];
      val1 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[1][ind01];
      val1 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[1][ind02];
      val1 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[1][ind03];
      val2  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[2][ind00];
      val2 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[2][ind01];
      val2 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[2][ind02];
      val2 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[2][ind03];
      val3  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[3][ind00];
      val3 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[3][ind01];
      val3 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[3][ind02];
      val3 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[3][ind03];
      val4  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[4][ind00];
      val4 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[4][ind01];
      val4 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[4][ind02];
      val4 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[4][ind03];
      val5  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[5][ind00];
      val5 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[5][ind01];
      val5 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[5][ind02];
      val5 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[5][ind03];
      val6  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[6][ind00];
      val6 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[6][ind01];
      val6 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[6][ind02];
      val6 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[6][ind03];

      val7  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[7][ind00];
      val7 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[7][ind01];
      val7 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[7][ind02];
      val7 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[7][ind03];

      val8  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[8][ind00];
      val8 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[8][ind01];
      val8 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[8][ind02];
      val8 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[8][ind03];

      val9  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[9][ind00];
      val9 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[9][ind01];
      val9 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[9][ind02];
      val9 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[9][ind03];

      val10  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[10][ind00];
      val10 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[10][ind01];
      val10 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[10][ind02];
      val10 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[10][ind03];

      val11  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[11][ind00];
      val11 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[11][ind01];
      val11 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[11][ind02];
      val11 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[11][ind03];

      val12  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[12][ind00];
      val12 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[12][ind01];
      val12 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[12][ind02];
      val12 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[12][ind03];

      val13  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[13][ind00];
      val13 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[13][ind01];
      val13 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[13][ind02];
      val13 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[13][ind03];

      val14  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[14][ind00];
      val14 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[14][ind01];
      val14 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[14][ind02];
      val14 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[14][ind03];

      val15  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[15][ind00];
      val15 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[15][ind01];
      val15 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[15][ind02];
      val15 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[15][ind03];

      val16  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[16][ind00];
      val16 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[16][ind01];
      val16 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[16][ind02];
      val16 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[16][ind03];

      val17  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[17][ind00];
      val17 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[17][ind01];
      val17 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[17][ind02];
      val17 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[17][ind03];

      val18  = ptrCoefs[ indCoef   ]*vectOfQneqDnrFields[18][ind00];
      val18 += ptrCoefs[ indCoef +1]*vectOfQneqDnrFields[18][ind01];
      val18 += ptrCoefs[ indCoef +2]*vectOfQneqDnrFields[18][ind02];
      val18 += ptrCoefs[ indCoef +3]*vectOfQneqDnrFields[18][ind03];

      vectOfQneqRcvFields[0][indR] = val0;
      vectOfQneqRcvFields[1][indR] = val1;
      vectOfQneqRcvFields[2][indR] = val2;
      vectOfQneqRcvFields[3][indR] = val3;
      vectOfQneqRcvFields[4][indR] = val4;
      vectOfQneqRcvFields[5][indR] = val5;
      vectOfQneqRcvFields[6][indR] = val6;
      vectOfQneqRcvFields[7][indR] = val7;
      vectOfQneqRcvFields[8][indR] = val8;
      vectOfQneqRcvFields[9][indR] = val9;
      vectOfQneqRcvFields[10][indR]= val10;
      vectOfQneqRcvFields[11][indR]= val11;
      vectOfQneqRcvFields[12][indR]= val12;
      vectOfQneqRcvFields[13][indR]= val13;
      vectOfQneqRcvFields[14][indR]= val14;
      vectOfQneqRcvFields[15][indR]= val15;
      vectOfQneqRcvFields[16][indR]= val16;
      vectOfQneqRcvFields[17][indR]= val17;
      vectOfQneqRcvFields[18][indR]= val18;
//      //macros
//#           include "commonInterpTransfers_reorder_macro.h"

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
      val0=0.; val1=0.; val2=0.; val3=0.; val4=0.; val5=0.; val6=0.; val7=0.; val8=0.; val9=0.;
      val10=0.; val11=1.; val12=0.; val13=0.; val14=0.; val15=0.; val16=0.; val17=0.; val18=0.;
      for (E_Int kk=0; kk<5; kk++)
        for (E_Int jj=0; jj<5; jj++)
          for (E_Int ii=0; ii<5; ii++)
          {
            indD = (i+ii)+(j+jj)*imd+(k+kk)*imdjmd;
            val0 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[0][indD];               
            val1 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[1][indD];               
            val2 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[2][indD];               
            val3 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[3][indD];               
            val4 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[4][indD];               
            val5 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[5][indD];
            val6 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[6][indD];               
            val7 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[7][indD];               
            val8 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[8][indD];               
            val9 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[9][indD];
            val10 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[10][indD];               
            val11 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[11][indD];               
            val12 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[12][indD];               
            val13 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[13][indD];               
            val14 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[14][indD];               
            val15 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[15][indD];               
            val16 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[16][indD];               
            val17 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[17][indD];               
            val18 += ptrCoefs[ indCoef + ii]*ptrCoefs[ indCoef + jj+5]*ptrCoefs[ indCoef + kk+10]*vectOfQneqDnrFields[18][indD];
          }
      vectOfQneqRcvFields[0][indR] = val0;
      vectOfQneqRcvFields[1][indR] = val1;
      vectOfQneqRcvFields[2][indR] = val2;
      vectOfQneqRcvFields[3][indR] = val3;
      vectOfQneqRcvFields[4][indR] = val4;
      vectOfQneqRcvFields[5][indR] = val5;
      vectOfQneqRcvFields[6][indR] = val6;
      vectOfQneqRcvFields[7][indR] = val7;
      vectOfQneqRcvFields[8][indR] = val8;
      vectOfQneqRcvFields[9][indR] = val9;
      vectOfQneqRcvFields[10][indR]= val10;
      vectOfQneqRcvFields[11][indR]= val11;
      vectOfQneqRcvFields[12][indR]= val12;
      vectOfQneqRcvFields[13][indR]= val13;
      vectOfQneqRcvFields[14][indR]= val14;
      vectOfQneqRcvFields[15][indR]= val15;
      vectOfQneqRcvFields[16][indR]= val16;
      vectOfQneqRcvFields[17][indR]= val17;
      vectOfQneqRcvFields[18][indR]= val18;
//      //macros
//#           include "commonInterpTransfers_reorder_macro.h"

      indCoef  += 15;
    }
    break;
      
  default: ;
}
