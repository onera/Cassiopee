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

# include "CompGeom/compGeom.h"

//===========================================================================
/* 
   Calcul le barycentre d'un ensemble de points.
   IN: n: nombre de pts dans le nuage
   IN: x, y, z: coord. du nuage de pts
   OUT: xb, yb, zb: coord. du barycentre
*/
//===========================================================================
void K_COMPGEOM::barycenter(
  E_Int n, 
  E_Float* x, E_Float* y, E_Float* z,
  E_Float& xb, E_Float& yb, E_Float& zb)
{
  E_Int nthreads = __NUMTHREADS__;
  E_Float* xp = new E_Float [3*nthreads];
  
  #pragma omp parallel
  {
    E_Int it = 3 * __CURRENT_THREAD__;
    xp[it] = 0.; xp[it+1] = 0.; xp[it+2] = 0.;
    #pragma omp for
    for (E_Int i = 0; i < n; i++) 
    { 
      xp[it] += x[i];
      xp[it+1] += y[i];
      xp[it+2] += z[i];
    }
  }
  E_Float onen = 1./n;
  xb = 0.; yb = 0.; zb = 0.;
  for (E_Int it = 0; it < nthreads; it++)
  {
    xb += xp[3*it]; yb += xp[3*it+1]; zb += xp[3*it+2];
  }
  xb = xb * onen; yb = yb * onen; zb = zb * onen;
  delete [] xp;
}

//=============================================================================
/* 
   Calcul le barycentre pondere d'un ensemble de points.
   IN: n: nombre de pts dans le nuage
   IN: x, y, z: coord. du nuage de pts
   IN: w: ponderation pour chaque point
   OUT: xb, yb, zb: coord. du barycentre
*/
//=============================================================================
void K_COMPGEOM::weightedBarycenter(E_Int n, 
E_Float* x, E_Float* y, E_Float* z, E_Float* w,
E_Float& xb, E_Float& yb, E_Float& zb)
{
  E_Int nthreads = __NUMTHREADS__;
  E_Float* xp = new E_Float [4*nthreads];
  
  #pragma omp parallel
  {
    E_Int it = 4 * __CURRENT_THREAD__;
    xp[it] = 0.; xp[it+1] = 0.; xp[it+2] = 0.; xp[it+3] = 0.;
    #pragma omp for
    for (E_Int i = 0; i < n; i++)
    { 
      xp[it] += x[i]*w[i]; 
      xp[it+1] += y[i]*w[i]; 
      xp[it+2] += z[i]*w[i];
      xp[it+3] += w[i];
    }
  }
  xb = 0.; yb = 0.; zb = 0.; E_Float wb = 0.;
  for (E_Int it = 0; it < nthreads; it++)
  {
    xb += xp[4*it]; yb += xp[4*it+1]; zb += xp[4*it+2];
    wb += xp[4*it+3];
  }
  wb = K_FUNC::E_max(wb, 1.e-10);
  E_Float onew = 1./wb;
  xb = xb * onew; yb = yb * onew; zb = zb * onew;
  delete [] xp;
}

//=============================================================================
/* Computes the barycenter of a NGON face 
   IN : posface : position of the face in NGon connectivity
   IN : ptrNG : pointer on the connectivity NGon
   OUT : coordinates of the barycenter of the face */
//=============================================================================
void K_COMPGEOM::getNGONFaceBarycenter(E_Int posface, E_Int* ptrNG, 
                                       E_Float* xt, E_Float* yt, E_Float* zt,
                                       E_Float& xbf, E_Float& ybf, E_Float& zbf)
{
  E_Int nvert = ptrNG[posface]; 
  E_Float invnvert = 1./nvert;
  xbf = 0.; ybf = 0.; zbf = 0.;
  E_Int ind;
  for (E_Int nv = 1; nv <= nvert; nv++)
  {
    ind = ptrNG[posface+nv]-1;
    xbf += xt[ind]; ybf += yt[ind]; zbf += zt[ind];
  }
  xbf = xbf*invnvert; ybf = ybf*invnvert; zbf = zbf*invnvert;
  return;
}
//=============================================================================
/* Computes the barycenter of a NGON  element
   IN : noet : number of the element
   IN : posEltsp : position of elts in connectivity NFace
   IN : posFacesp : position of faces in connectivity NGon
   IN : posface : position of the face in NGon connectivity
   IN : ptrNF : connectivity Elt/Faces (NFace connectivity)
   OUT : coordinates of the barycenter of the face */
//=============================================================================
void K_COMPGEOM::getNGONEltBarycenter(E_Int et, E_Int* posEltsp, 
                                      E_Int* posFacesp, 
                                      E_Int* ptrNF, E_Int* ptrNG,
                                      E_Float* xt, E_Float* yt, E_Float* zt,
                                      E_Float& xbg, E_Float& ybg, E_Float& zbg)
{
  E_Int pose = posEltsp[et];
  E_Int nf = ptrNF[pose];
  E_Int face, posf;
  E_Float xbf, ybf, zbf;
  xbg = 0.; ybg = 0.; zbg = 0.;
  for (E_Int nof = 1; nof <= nf; nof++)
  {
    face = ptrNF[pose+nof]-1;
    posf = posFacesp[face];  
    getNGONFaceBarycenter(posf, ptrNG, xt, yt, zt, xbf, ybf, zbf); 
    xbg += xbf; ybg += ybf; zbg += zbf;
  }
  xbg = xbg/nf; ybg = ybg/nf; zbg = zbg/nf;
}

