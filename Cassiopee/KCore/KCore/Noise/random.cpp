/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Noise/noise.h"
#include "Def/DefTypes.h"

//=============================================================================
// Minimum standard random generator
// Retourne un reel entre 0 et 1.
// Tirages correles. A initialiser avec un entier negatif.
//=============================================================================
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double K_NOISE::stdRand(E_LONG *idum)
{
  E_LONG k;
  double ans;

  *idum ^= MASK;
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;			//cycling
  
  ans=AM*(*idum);
  *idum ^= MASK;
  return ans;
}

//=============================================================================
// Minimum shuffle random generator
// Retourne un reel entre 0 et 1.
// Non correle pour un nbre de tirages < a 10e8.
// A initialiser avec un entier negatif.
//=============================================================================
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-10
#define RNMX (1.0-EPS)

double K_NOISE::shuffleRand(E_LONG *idum)
{
  int j;
  E_LONG k;
  static E_LONG iy=0;
  static E_LONG iv[NTAB];
  double temp;
  
  if (*idum <=0 || !iy)
  {
    if (-(*idum)<1) *idum=1;
    else *idum=-(*idum);
    for (j=NTAB+7;j>=0;j--)
    {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum<0) *idum += IM;
      if (j<NTAB) iv[j]=*idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum<0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j]=*idum;
  if ((temp=AM*iy)>RNMX) temp=RNMX;
  return temp;
}

//=============================================================================
// Minimum combined random generator
// Retourne un reel entre 0 et 1.
// Tres longue periode.
// A initialiser avec un entier negatif.
//=============================================================================
#define IM1 2147483563
#define IM2 2147483399
#define AM1 (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NDIV1 (1+IMM1/NTAB)

double K_NOISE::longRand(E_LONG *idum)
{
  int j;
  E_LONG k;
  static E_LONG idum2=123456789;
  static E_LONG iy=0;
  static E_LONG iv[NTAB];
  double temp;
  
  if (*idum <=0 || !iy)
  {
    if (-(*idum)<1) *idum=1;
    else *idum=-(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--)
    {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-IR1*k;
      if (*idum<0) *idum += IM1;
      if (j<NTAB) iv[j]=*idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-IR1*k;
  if (*idum<0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-IR2*k;
  if (idum2<0) idum2 += IM2;

  j=iy/NDIV1;
  iy=iv[j]-idum2;
  iv[j]=*idum;
  if (iy<1) iy += IMM1;
  if ((temp=AM1*iy)>RNMX) temp=RNMX;
  return temp;	
}
