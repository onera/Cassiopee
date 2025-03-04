/*    
    Copyright 2013-2025 Onera.

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
  FldArrayF fout(im*jm*km, f.getNfld());
  reorderStructField(im, jm, km, 
                     f, fout, oi, oj, ok);
  f = fout;
}
//=============================================================================
void K_CONNECT::reorderStructField(
  E_Int& im, E_Int& jm, E_Int& km, 
  FldArrayF& f, FldArrayF& fout, 
  E_Int oi, E_Int oj, E_Int ok)
{
  // reordering
  E_Int delta, epsilon;
  E_Int nfld = f.getNfld();
  E_Int imjm = im*jm;
  E_Int in=0, jn=0, kn=0;

  switch (oi)
  {
    case 1:
      in = im; break;
    case -1:
      in = im; break;
    case 2:
      jn = im; break;
    case -2:
      jn = im; break;
    case 3:
      kn = im; break;
    case -3:
      kn = im; break;
    default:
      printf("Error: reorder: bad value of oi, oj or ok.\n"); 
      exit(0);
  }
  switch (oj)
  {
    case 1:
      in = jm; break;
    case -1:
      in = jm; break;
    case 2:
      jn = jm; break;
    case -2:
      jn = jm; break;
    case 3:
      kn = jm; break;
    case -3:
      kn = jm; break;
    default:
      printf("Error: reorder: bad value of oi, oj or ok.\n"); 
      exit(0);
  }
  switch (ok)
  {
    case 1:
      in = km; break;
    case -1:
      in = km; break;
    case 2:
      jn = km; break;
    case -2:
      jn = km; break;
    case 3:
      kn = km; break;
    case -3:
      kn = km; break;
    default:
      printf("Error: reorder: bad value of oi, oj or ok.\n"); 
      exit(0);
  }
  delta = in; epsilon = in*jn;
  
#pragma omp parallel default(shared)
  {
    E_Int i, j, k, ind2, alpha=0, beta=0, gamma=0;

#pragma omp for
    for (E_Int ind = 0; ind < imjm*km; ind++)
    { 
      k = ind / imjm;
      j = (ind-k*imjm)/im;
      i = ind-j*im-k*imjm;
      
      switch (oi)
      {
        case 1:
          alpha = i; break;
        case -1:
          alpha = im-i-1; break;
        case 2:
          beta = i; break;
        case -2:
          beta = im-i-1; break;
        case 3:
          gamma = i; break;
        case -3:
          gamma = im-i-1; break;
        default: //erreur deja testee
          exit(0);
      }
      switch (oj)
      {
        case 1:
          alpha = j; break;
        case -1:
          alpha = jm-j-1; break;
        case 2:
          beta = j; break;
        case -2:
          beta = jm-j-1; break;
        case 3:
          gamma = j; break;
        case -3:
          gamma = jm-j-1; break;
        default: //erreur deja testee
          exit(0);
      }
      switch (ok)
      {
        case 1:
          alpha = k; break;
        case -1:
          alpha = km-k-1; break;
        case 2:
          beta = k; break;
        case -2:
          beta = km-k-1; break;
        case 3:
          gamma = k; break;
        case -3:
          gamma = km-k-1; break;
        default: //erreur deja testee
          exit(0);
      }
      ind2 = alpha+beta*delta+gamma*epsilon;
      for (E_Int n = 1; n <= nfld; n++) fout(ind2,n) = f(ind,n);
    }  
  }

  im = in;
  jm = jn;
  km = kn;
}
