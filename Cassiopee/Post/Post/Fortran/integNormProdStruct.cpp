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

# include "post.h"

//=============================================================================
// Integre les grandeurs de vect(F).vect(n)
// Retourne 1 si success, 0 si echec
//=============================================================================
E_Int K_POST::integNormProdStruct2D(E_Int ni, E_Int nj, E_Int nk, 
                                    E_Int center2node,
                                    E_Int posx, E_Int posy, E_Int posz,
                                    FldArrayF& coord, FldArrayF& F, 
                                    FldArrayF& ratio, E_Float& resultat)
{
  E_Int NI, NJ;
  E_Float result = 0.;
  
  if (nk == 1) { NI = ni; NJ = nj; }
  else if (nj == 1) { NI = ni; NJ = nk; }
  else if (ni == 1) { NI = nj; NJ = nk; }
  else return 0;
 
  // Compute surface of each "block" i cell, with coordinates coord
  FldArrayF nsurf((NI-1) * (NJ-1), 3);

  //E_Int npts = coord.getSize();
  K_METRIC::compNormStructSurf(
    NI, NJ, coord.begin(posx), coord.begin(posy), coord.begin(posz), 
    nsurf.begin(1), nsurf.begin(2), nsurf.begin(3));

  if (center2node == 1)
  {
    // Compute integral, coordinates defined in node 
    // and field F in center 
    integNormProdStructCellCenter2D(
      NI-1, NJ-1, ratio.begin(), 
      nsurf.begin(1), nsurf.begin(2), nsurf.begin(3),
      F.begin(1), F.begin(2), F.begin(3),
      result
    );
  }
  else
  {
    // Compute integral, coordinates and field have the same size
    integNormProdStructNodeCenter2D(
      NI, NJ, ratio.begin(),
      nsurf.begin(1), nsurf.begin(2),  nsurf.begin(3),  
      F.begin(1), F.begin(2), F.begin(3),
      result
    );
  }
  resultat += result;   
  return 1;
}

//==============================================================================
// Compute surface integral of the product vect(F).vect(n), coordinates
// are defined in nodes and F is defined in center
//==============================================================================
void K_POST::integNormProdStructNodeCenter2D(
  const E_Int ni, const E_Int nj,
  const E_Float* ratio,
  const E_Float* sx, const E_Float* sy, const E_Float* sz,
  const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float& result)
{
  E_Int ni1 = ni - 1;
  result = 0.0;
  E_Int nthreads = __NUMTHREADS__;
  E_Float* reti = new E_Float [nthreads];

  #pragma omp parallel
  {
    E_Int ind, ind1, ind2, ind3, ind4;
    E_Float r1, r2, r3, r4;
    E_Float fx, fy, fz;
    E_Int it = __CURRENT_THREAD__;
    reti[it] = 0.;

    #pragma omp for
    for (E_Int j = 0; j < nj - 1; j++)
    {
      for (E_Int i = 0; i < ni - 1; i++)
      {
        ind1 = i + j * ni;
        ind2 = ind1 + ni;
        ind3 = ind1 + 1;
        ind4 = ind3 + ni;

        ind = i + j * ni1;

        r1 = ratio[ind1];
        r2 = ratio[ind2];
        r3 = ratio[ind3];
        r4 = ratio[ind4];

        fx = r1 * vx[ind1] + r2 * vx[ind2] + r3 * vx[ind3] + r4 * vx[ind4];
        fy = r1 * vy[ind1] + r2 * vy[ind2] + r3 * vy[ind3] + r4 * vy[ind4];
        fz = r1 * vz[ind1] + r2 * vz[ind2] + r3 * vz[ind3] + r4 * vz[ind4];

        reti[it] += sx[ind] * fx + sy[ind] * fy + sz[ind] * fz;
      }
    }
  }

  for (E_Int it = 0; it < nthreads; it++)
  {
    result += reti[it];
  }

  result *= 0.25;

  delete [] reti;
}

//==============================================================================
// Compute surface integral of the product vect(F).vect(n), coordinates 
// are defined in nodes and F is defined in nodes (center-based formulation)
//==============================================================================
void K_POST::integNormProdStructCellCenter2D(
  const E_Int ni, const E_Int nj,
  const E_Float* ratio,
  const E_Float* sx, const E_Float* sy, const E_Float* sz,
  const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float& result)
{
  result = 0.0;
  E_Int nthreads = __NUMTHREADS__;
  E_Float* reti = new E_Float [nthreads];

  #pragma omp parallel
  {
    E_Int ind;
    E_Float sum;
    E_Int it = __CURRENT_THREAD__;
    reti[it] = 0.;

    #pragma omp for
    for (E_Int j = 0; j < nj; j++)
    {
      for (E_Int i = 0; i < ni; i++)
      {
        ind = i + j * ni;
        sum = sx[ind] * vx[ind] + sy[ind] * vy[ind] + sz[ind] * vz[ind];
        reti[it] += ratio[ind] * sum;
      }
    }
  }

  for (E_Int it = 0; it < nthreads; it++)
  {
    result += reti[it];
  }

  delete [] reti;
}