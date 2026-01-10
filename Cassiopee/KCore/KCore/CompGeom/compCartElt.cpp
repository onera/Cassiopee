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
# include "String/kstring.h"
# include <stdio.h>
# include <stdlib.h>

using namespace K_FLD;
using namespace std;

extern "C"
{
  void k6compcartelembox_(
    const E_Int& is1, const E_Int& is2,
    const E_Int& js1, const E_Int& js2,
    const E_Int& ks1, const E_Int& ks2,
    const E_Int& i2, const E_Int& j2, const E_Int& k2,
    const E_Int& im, const E_Int& jm, const E_Int& km,
    const E_Float* x, const E_Float* y,  const E_Float* z, 
    const E_Float* rotMat,
    const E_Float& du, const E_Float& dv,
    const E_Int& nbElts1, const E_Int& nbElts2,
    const E_Int& dir, const E_Int& area,
    const E_Float* nodeMin,
    const E_Float* nodeMax,
    const E_Int& szCartElt,
    E_Float* cartEltMin, E_Float* cartEltMax);
}
//=============================================================================
/* Calcul le tableau des elements cartesiens dans la direction dir 
   Comme il s'agit d elements cartesiens, sont stockes uniquemt  
   xmin, xmax de chq element : cartEltArray(elt,1) = xmin, ...
*/
//=============================================================================
short 
K_COMPGEOM::compCartEltsArray(E_Int dir, E_Int ni, E_Int nj, E_Int nk,
                              E_Float xmin, E_Float ymin, E_Float zmin,
                              E_Float xmax, E_Float ymax, E_Float zmax,
                              FldArrayF& coord, FldArrayF& cartEltArray)
{ 
  // Calcul des dimensions de la CEBB
  E_Int nbElts1, nbElts2;
  compDimOfCartElts(dir, ni, nj, nk, xmin, ymin, zmin,
                    xmax, ymax, zmax, coord, nbElts1, nbElts2);

  // Calcul de la hauteur des elemts cartesiens
  E_Int nCartElts = nbElts1*nbElts2;
  FldArrayF cartMin(nCartElts);
  FldArrayF cartMax(nCartElts);
  cartMin.setAllValuesAt(+K_CONST::E_MAX_FLOAT);
  cartMax.setAllValuesAt(-K_CONST::E_MAX_FLOAT);
  E_Float* cartMinp = cartMin.begin();
  E_Float* cartMaxp = cartMax.begin();
  short isOk = compCEBBHeightInDir(dir, nbElts1, nbElts2, ni, nj, nk, 
                                   xmin, ymin, zmin, xmax, ymax, zmax,
                                   coord, cartMin, cartMax);
  
  if (isOk == 0) return 0;
 
  cartEltArray.malloc(nCartElts, 6);
  E_Float dx = xmax-xmin;
  E_Float dy = ymax-ymin;
  E_Float dz = zmax-zmin;
  E_Float du = 1;
  E_Float dv = 1;
  E_Float xmin0 = xmin;
  E_Float ymin0 = ymin;
  E_Float zmin0 = zmin;
  E_Float xmin1 = 0.;
  E_Float ymin1 = 0.;
  E_Float zmin1 = 0.;
  E_Float xmax1 = 0.;
  E_Float ymax1 = 0.;
  E_Float zmax1 = 0.; 

  switch (dir)
  {
    case 1:
      du = dy/nbElts1;
      dv = dz/nbElts2;
      ymin0 = ymin0-du;
      ymin1 = ymin0;
      zmin0 = zmin0-dv;
      break;

    case 2:
      du = dz/nbElts1;
      dv = dx/nbElts2;
      xmin0 = xmin0-dv;
      zmin0 = zmin0-du;
      zmin1 = zmin0;
      break;

    case 3 :
      du = dx/nbElts1;
      dv = dy/nbElts2;
      xmin0 = xmin0-du;
      xmin1 = xmin0;
      ymin0 = ymin0-dv;
      break;
  }

  E_Int elt = 0;  
  for (E_Int i = 0 ; i < nbElts1; i++)
  {
    switch (dir)
    {
      case 1:
        ymin1 = ymin1+du;
        ymax1 = ymin1+du;
        zmin1 = zmin0;
        break;
      case 2:
        zmin1 = zmin1+du;
        zmax1 = zmin1+du;
        xmin1 = xmin0;
        break;
      case 3:
        xmin1 = xmin1+du;
        xmax1 = xmin1+du;
        ymin1 = ymin0;
        break;
    }
    for (E_Int j = 0; j < nbElts2; j++)
    {
      E_Int k = i+j*nbElts1;
      switch (dir)
      {
        case 1:
          xmin1 = cartMinp[k];
          xmax1 = cartMaxp[k];
          zmin1 = zmin1+dv;
          zmax1 = zmin1+dv;
          break;
        case 2:
          xmin1 = xmin1+dv;
          xmax1 = xmin1+dv;
          ymin1 = cartMinp[k];
          ymax1 = cartMaxp[k];
          break;
        case 3:
          ymin1 = ymin1+dv;
          ymax1 = ymin1+dv;
          zmin1 = cartMinp[k];
          zmax1 = cartMaxp[k];
          break;
      }

      // element has been set
      if (xmin1 < K_CONST::E_MAX_FLOAT-1. && xmax1> -K_CONST::E_MAX_FLOAT+1. &&
          ymin1 < K_CONST::E_MAX_FLOAT-1. && ymax1> -K_CONST::E_MAX_FLOAT+1. &&
          zmin1 < K_CONST::E_MAX_FLOAT-1. && zmax1> -K_CONST::E_MAX_FLOAT+1.) 
      {
        cartEltArray(elt,1) = xmin1;
        cartEltArray(elt,2) = ymin1;
        cartEltArray(elt,3) = zmin1;
        cartEltArray(elt,4) = xmax1;
        cartEltArray(elt,5) = ymax1;
        cartEltArray(elt,6) = zmax1;
        elt++;
      }
    }
  }
  cartEltArray.reAllocMat(elt,6);
  if (elt == 0) return 0;
  else return 1;
}
//=============================================================================
/* Calcul des dimensions dim1, dim2 des elements cartesiens */
//=============================================================================
void K_COMPGEOM::compDimOfCartElts(E_Int dir, E_Int ni, E_Int nj, E_Int nk,
                                   E_Float xmin, E_Float ymin, E_Float zmin,
                                   E_Float xmax, E_Float ymax, E_Float zmax,
                                   FldArrayF& coord,
                                   E_Int& nbElts1, E_Int& nbElts2)
{
  E_Float sizex = xmax - xmin;
  E_Float sizey = ymax - ymin;
  E_Float sizez = zmax - zmin;

  E_Float dxMax = 0.;
  E_Float dyMax = 0.;
  E_Float dzMax = 0.;

  E_Int ninj = ni*nj;
  E_Float dx, dy, dz;
  E_Int ind1, ind2;

  for (E_Int k = 0; k < nk; k++)
    for (E_Int j = 0; j < nj; j++)
      for (E_Int i = 0; i < ni-1; i++)
      {
        ind1 = i + j*ni + k *ninj;
        ind2 = ind1 + 1;
        dx = K_FUNC::E_abs(coord(ind2,1)-coord(ind1,1));
        dy = K_FUNC::E_abs(coord(ind2,2)-coord(ind1,2));
        dz = K_FUNC::E_abs(coord(ind2,3)-coord(ind1,3));

        dxMax = K_FUNC::E_max(dxMax,dx);
        dyMax = K_FUNC::E_max(dyMax,dy);
        dzMax = K_FUNC::E_max(dzMax,dz);   
      }

  for (E_Int k = 0; k < nk; k++)
    for (E_Int j = 0; j < nj-1; j++)
      for (E_Int i = 0; i < ni; i++)
      {
        ind1 = i + j*ni + k *ninj;
        ind2 = ind1 + ni;
        dx = K_FUNC::E_abs(coord(ind2,1)-coord(ind1,1));
        dy = K_FUNC::E_abs(coord(ind2,2)-coord(ind1,2));
        dz = K_FUNC::E_abs(coord(ind2,3)-coord(ind1,3));

        dxMax = K_FUNC::E_max(dxMax,dx);
        dyMax = K_FUNC::E_max(dyMax,dy);
        dzMax = K_FUNC::E_max(dzMax,dz);   
      }
  
  for (E_Int k = 0; k < nk-1; k++)
    for (E_Int j = 0; j < nj; j++)
      for (E_Int i = 0; i < ni; i++)
      {
        ind1 = i + j*ni + k *ninj;
        ind2 = ind1 + ninj;
        dx = K_FUNC::E_abs(coord(ind2,1)-coord(ind1,1));
        dy = K_FUNC::E_abs(coord(ind2,2)-coord(ind1,2));
        dz = K_FUNC::E_abs(coord(ind2,3)-coord(ind1,3));

        dxMax = K_FUNC::E_max(dxMax,dx);
        dyMax = K_FUNC::E_max(dyMax,dy);
        dzMax = K_FUNC::E_max(dzMax,dz);   
      }

  E_Int nx = K_FUNC::E_max(1, E_Int(sizex/dxMax));
  E_Int ny = K_FUNC::E_max(1, E_Int(sizey/dyMax));
  E_Int nz = K_FUNC::E_max(1, E_Int(sizez/dzMax));
  
  switch(dir)
  {
    case 1:
      nbElts1 = ny;
      nbElts2 = nz;
      break;
    case 2:
      nbElts1 = nz;
      nbElts2 = nx;
      break;
    case 3:
      nbElts1 = nx;
      nbElts2 = ny;
      break;
    default:
      printf("Error: compDimOfCartElts: bad dir value: " SF_D_ ".", dir);
      exit(0);
  }
}

//=============================================================================
/* Calcul de la hauteur des elements cartesiens ds la direction consideree */
//=============================================================================
short 
K_COMPGEOM::compCEBBHeightInDir(E_Int dir, E_Int nbElts1, E_Int nbElts2,
                                E_Int ni, E_Int nj, E_Int nk, 
                                E_Float xmin, E_Float ymin, E_Float zmin,
                                E_Float xmax, E_Float ymax, E_Float zmax,
                                FldArrayF& coord,
                                FldArrayF& cartMin, FldArrayF& cartMax)
{ 
  E_Float du, dv;
  E_Float sizex = xmax - xmin;
  E_Float sizey = ymax - ymin;
  E_Float sizez = zmax - zmin;

  switch (dir)
  {
    case 1:
      if (K_FUNC::fEqualZero(sizey) != true && 
          K_FUNC::fEqualZero(sizez) != true)
      {
        du = sizey/nbElts1;
        dv = sizez/nbElts2;
      }
      else return 0;     
      break;

    case 2:
      if (K_FUNC::fEqualZero(sizez) != true &&
          K_FUNC::fEqualZero(sizex) != true)
      {
        du = sizez/nbElts1;
        dv = sizex/nbElts2;
      }
      else return 0;
      break;

    case 3:
      if (K_FUNC::fEqualZero(sizex) != true &&
          K_FUNC::fEqualZero(sizey) != true)
      {
        du = sizex/nbElts1;
        dv = sizey/nbElts2;
      }
      else return 0;
      break;
  }
    
  FldArrayF nodeMin(3);
  FldArrayF nodeMax(4);
  nodeMin[0] = xmin; nodeMin[1] = ymin; nodeMin[2] = zmin;
  nodeMax[0] = xmax; nodeMax[1] = ymax; nodeMax[2] = zmax;
  FldArrayF rotMat(3,3);
  rotMat.setAllValuesAtNull();
  rotMat(0,1) = 1.; rotMat(1,2) = 1.; rotMat(2,3) = 1.;

  k6compcartelembox_(1, ni, 1, nj, 1, nk, ni, nj, nk, ni, nj, nk, 
                     coord.begin(1), coord.begin(2), coord.begin(3),
                     rotMat.begin(), 
                     du, dv, nbElts1, nbElts2, dir, 0,
                     nodeMin.begin(), nodeMax.begin(), nbElts1*nbElts2, 
                     cartMin.begin(), cartMax.begin());
                          
  return 1;
}
//============================ CompGeom/compCartElt.cpp =======================
