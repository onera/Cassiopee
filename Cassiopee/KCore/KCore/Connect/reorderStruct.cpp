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

#include "Connect/connect.h"
#include <stdio.h>
#include <stdlib.h>

using namespace K_FLD;
using namespace std;
 
// ============================================================================
/* Reorder la numerotation d'un array structure suivant oi, oj, ok */
// ============================================================================
void K_CONNECT::reorderStructField(
  E_Int& im, E_Int& jm, E_Int& km, 
  FldArrayF& f,
  E_Int oi, E_Int oj, E_Int ok)
{
  E_Int npts = im*jm*km;
  E_Int nfld = f.getNfld();
  FldArrayF fout(npts, nfld);

  reorderStructField(im, jm, km, 
                     f, fout, oi, oj, ok);

  #pragma omp parallel for collapse(2)
  for (E_Int ind = 0; ind < npts; ind++)
  for (E_Int n = 1; n <= nfld; n++) f(ind,n) = fout(ind,n);
}
//=============================================================================
void K_CONNECT::reorderStructField(
  E_Int& im, E_Int& jm, E_Int& km, 
  FldArrayF& f, FldArrayF& fout, 
  E_Int oi, E_Int oj, E_Int ok)
{
  E_Int nfld = f.getNfld();
  E_Int imjm = im*jm;

  // Validate oi, oj, ok
  if (K_FUNC::E_abs(oi) < 1 || K_FUNC::E_abs(oi) > 3 ||
      K_FUNC::E_abs(oj) < 1 || K_FUNC::E_abs(oj) > 3 ||
      K_FUNC::E_abs(ok) < 1 || K_FUNC::E_abs(ok) > 3)
  {
    printf("Error: reorder: bad value of oi, oj or ok.\n");
    exit(0);
  }

  // axisI/axisJ/axisK: which output axis (0=alpha,1=beta,2=gamma)
  // each input direction maps to
  E_Int axisI = K_FUNC::E_abs(oi)-1;
  E_Int axisJ = K_FUNC::E_abs(oj)-1;
  E_Int axisK = K_FUNC::E_abs(ok)-1;

  // Reject degenerate mappings (two inputs mapping to the same output axis)
  if (axisI == axisJ || axisJ == axisK || axisI == axisK)
  {
    printf("Error: reorder: bad value of oi, oj or ok.\n");
    exit(0);
  }

  // signI/signJ/signK: +1 if direction is kept as-is, -1 if reversed
  E_Int signI = (oi > 0) ? 1 : -1;
  E_Int signJ = (oj > 0) ? 1 : -1;
  E_Int signK = (ok > 0) ? 1 : -1;

  // New dimensions: dims[axis] = size of the input dimension mapped to that axis
  E_Int dims[3];
  dims[axisI] = im;
  dims[axisJ] = jm;
  dims[axisK] = km;
  E_Int in = dims[0], jn = dims[1], kn = dims[2];

  E_Int delta = in, epsilon = in*jn;

  #pragma omp parallel
  {
    E_Int i, j, k, ind2;
    E_Int coord[3];
    E_Int vi, vj, vk;

    #pragma omp for
    for (E_Int ind = 0; ind < imjm*km; ind++)
    { 
      k = ind / imjm;
      j = (ind-k*imjm)/im;
      i = ind-j*im-k*imjm;

      vi = (signI > 0) ? i : im-i-1;
      vj = (signJ > 0) ? j : jm-j-1;
      vk = (signK > 0) ? k : km-k-1;

      coord[axisI] = vi;
      coord[axisJ] = vj;
      coord[axisK] = vk;

      ind2 = coord[0] + coord[1]*delta + coord[2]*epsilon;
      for (E_Int n = 1; n <= nfld; n++) fout(ind2,n) = f(ind,n);
    }
  }

  im = in;
  jm = jn;
  km = kn;
}
