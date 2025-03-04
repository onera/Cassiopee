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
  E_Float x1, y1, z1;
  x1 = 0.; y1 = 0.; z1 = 0.;
#pragma omp parallel for shared(n,x,y,z) reduction(+:x1,y1,z1) 
  for (int i = 0; i < n; i++) 
  { x1 += x[i]; y1 += y[i]; z1 += z[i]; }

  E_Float onen = 1./n;
  xb = x1 * onen; yb = y1 * onen; zb = z1 * onen;
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
  E_Float wg;
  E_Float x1, y1, z1, wb;
  x1 = 0.; y1 = 0.; z1 = 0.; wb = 0.;
#pragma omp parallel for shared(n,x,y,z) private(wg) reduction(+:x1,y1,z1,wb) 
  for (int i = 0; i < n; i++)
  {  
    wg = w[i];
    x1 += wg*x[i]; y1 += wg*y[i]; z1 += wg*z[i];
    wb += wg;
  }

  wb = K_FUNC::E_max(wb, 1.e-10);
  E_Float onew = 1./wb;
  xb = x1 * onew; yb = y1 * onew; zb = z1 * onew;
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

