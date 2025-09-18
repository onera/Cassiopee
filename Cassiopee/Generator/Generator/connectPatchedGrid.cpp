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

#include "generator.h"

using namespace std;
using namespace K_FLD;
//=============================================================================
/* Determination des raccords match/nearmatch pour un maillage cartesien de 
   type patched grids 
   IN: now1: numero associe a la fenetre de la zone1 (now1 = 1 ; i = 1...)
   IN: array1: array de la zone1
   IN: array2: array de la zone2
   OUT: si raccords coincidents : retourne [1,[imin1,imax1,jmin1,jmax1,kmin1,kmax1],[imin2,imax2,jmin2,jmax2,kmin2,kmax2]]
        si raccords 1 point sur 2 : retourne [2,[imin1,imax1,jmin1,jmax1,kmin1,kmax1],[imin2,imax2,jmin2,jmax2,kmin2,kmax2]]
         sinon: retourne [] */
//=============================================================================
PyObject* K_GENERATOR::connectPatchGrid(PyObject* self, PyObject* args)
{
  E_Float eps = 1.e-10;
  PyObject *array1, *array2;
  E_Float xmin1, xmax1, ymin1, ymax1, zmin1, zmax1;
  E_Float xmin2, xmax2, ymin2, ymax2, zmin2, zmax2;
#ifdef E_DOUBLEREAL
  if (!PyArg_ParseTuple(args, "OO(dddddd)(dddddd)", &array1, &array2, 
                        &xmin1, &ymin1, &zmin1, &xmax1, &ymax1, &zmax1,
                        &xmin2, &ymin2, &zmin2, &xmax2, &ymax2, &zmax2))
#else
  if (!PyArg_ParseTuple(args, "OO(ffffff)(ffffff)", &array1, &array2, 
                        &xmin1, &ymin1, &zmin1, &xmax1, &ymax1, &zmax1,
                        &xmin2, &ymin2, &zmin2, &xmax2, &ymax2, &zmax2)) 
#endif
  return NULL;

  // identification du no de la win de la zone z1
  E_Int nowin1 = 0;

  if (K_FUNC::fEqualZero(xmin1-xmax2,eps) == true) nowin1 = 1;
  else if (K_FUNC::fEqualZero(xmax1-xmin2,eps) == true) nowin1 = 2;
  else if (K_FUNC::fEqualZero(ymin1-ymax2,eps) == true) nowin1 = 3;
  else if (K_FUNC::fEqualZero(ymax1-ymin2,eps) == true) nowin1 = 4;
  else if (K_FUNC::fEqualZero(zmin1-zmax2,eps) == true) nowin1 = 5;
  else if (K_FUNC::fEqualZero(zmax1-zmin2,eps) == true) nowin1 = 6;
  else return Py_None;

  // Check arrays
  E_Int ni1, nj1, nk1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray3(array1, varString1, f1, 
                                      ni1, nj1, nk1, cn1, eltType1); 
  if (res1 != 1)
  {
    RELEASESHAREDB(res1, array1, f1, cn1);
    PyErr_SetString(PyExc_TypeError,
                    "connectPatchGrid: 1st array must be structured.");
    return NULL;
  }
  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);

  if (posx1 == -1 || posy1 == -1 || posz1 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "connectPatchGrid: 1st array must contain (x,y,z).");
    RELEASESHAREDS(array1, f1); return NULL;
  }
  posx1++; posy1++; posz1++; 

  E_Int ni2, nj2, nk2;
  FldArrayF* f2; FldArrayI* cn2;
  char* varString2; char* eltType2;
  E_Int res2 = K_ARRAY::getFromArray3(array2, varString2, f2, 
                                      ni2, nj2, nk2, cn2, eltType2); 
  if (res2 != 1)
  {
    RELEASESHAREDS(array1, f1);
    RELEASESHAREDB(res2, array2, f2, cn2);
    PyErr_SetString(PyExc_TypeError,
                    "connectPatchGrid: 2nd array must be structured.");
    return NULL;
  }

  E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);

  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "connectPatchGrid: 2nd array must contain (x,y,z).");
    RELEASESHAREDS(array1, f1); RELEASESHAREDS(array2, f2);
    return NULL;
  }
  posx2++; posy2++; posz2++;
  
  // detection des raccords match/nearmatch
  E_Float* xt1 = f1->begin(posx1); E_Float* xt2 = f2->begin(posx2);
  E_Int im1, jm1, km1, im2, jm2, km2, ip1, jp1, kp1, ip2, jp2, kp2;
  E_Float dh1 = xt1[1]-xt1[0];//pas d espace sur les facettes
  E_Float dh2 = xt2[1]-xt2[0];
  E_Int dlevel = 0;
  E_Float res;
  if (K_FUNC::fEqualZero(dh1-dh2, eps) == true) dlevel = 1;
  else 
  {
    if (dh1 > dh2) {res = dh1/dh2; dlevel = E_Int(res+eps);}
    else { res = dh2/dh1; dlevel = E_Int(res+eps); }
  }
  E_Int found = 0;
  if (dlevel > 0) found = detectMatchCart(nowin1, ni1, nj1, nk1, 
                                          f1->begin(posx1),f1->begin(posy1),f1->begin(posz1),
                                          ni2, nj2, nk2, 
                                          f2->begin(posx2),f2->begin(posy2),f2->begin(posz2),
                                          im1, ip1, jm1, jp1, km1, kp1, 
                                          im2, ip2, jm2, jp2, km2, kp2);
  
  RELEASESHAREDS(array1, f1); RELEASESHAREDS(array2, f2);
  if (found != 0)
  {
    E_LONG imin1 = im1; E_LONG jmin1 = jm1; E_LONG kmin1 = km1;
    E_LONG imax1 = ip1; E_LONG jmax1 = jp1; E_LONG kmax1 = kp1;
    E_LONG imin2 = im2; E_LONG jmin2 = jm2; E_LONG kmin2 = km2;
    E_LONG imax2 = ip2; E_LONG jmax2 = jp2; E_LONG kmax2 = kp2;
    E_LONG found0 = dlevel*found;
    PyObject* info = Py_BuildValue("[l,[llllll],[llllll]]", found0,
                                   imin1, imax1, jmin1, jmax1, kmin1, kmax1,
                                   imin2, imax2, jmin2, jmax2, kmin2, kmax2);
    return info;
  }
  Py_INCREF(Py_None);
  return Py_None;
}
//=============================================================================
/* detection des raccords match pour les maillages cartesiens */
//=============================================================================
E_Int K_GENERATOR::detectMatchCart(
  E_Int nowin1, E_Int ni1, E_Int nj1, E_Int nk1,
  E_Float* xt1, E_Float* yt1, E_Float* zt1, 
  E_Int ni2, E_Int nj2, E_Int nk2,
  E_Float* xt2, E_Float* yt2, E_Float* zt2,
  E_Int& im1, E_Int& ip1, E_Int& jm1, E_Int& jp1, E_Int& km1, E_Int& kp1,
  E_Int& im2, E_Int& ip2, E_Int& jm2, E_Int& jp2, E_Int& km2, E_Int& kp2)
{
  E_Float eps = 1.e-10; E_Float eps2 = eps*eps;
  E_Int inc1, inc2, ind1, ind2;
  E_Int ni1nj1 = ni1*nj1; E_Int ni2nj2 = ni2*nj2;
  E_Float dx, dy, dz;
  im1 = ni1; ip1 =-1; jm1 = nj1; jp1 =-1; km1 = nk1; kp1 =-1;
  im2 = ni2; ip2 =-1; jm2 = nj2; jp2 =-1; km2 = nk2; kp2 =-1;
  E_Int dim = 2;
  if (nk1 != 1 && nj1 != 1) dim = 3;

  if (nowin1 == 1 || nowin1 == 2) //faces xmin1 avec xmax2
  {
    if (nowin1 == 1) 
    { inc1 = 0; inc2 = ni2-1; im1 = 1; ip1 = 1; im2 = ni2; ip2 = ni2; }
    else
    { inc1 = ni1-1; inc2 = 0; im1 = ni1; ip1 = ni1; im2 = 1; ip2 = 1; }
    for (E_Int k1 = 1; k1 <= nk1; k1++)
      for (E_Int j1 = 1; j1 <= nj1; j1++)
      {
        ind1 = inc1 + (j1-1)*ni1 + (k1-1)*ni1nj1;
        for (E_Int k2 = 1; k2 <= nk2; k2++)
          for (E_Int j2 = 1; j2 <= nj2; j2++)
          {       
            ind2 = inc2 + (j2-1)*ni2 + (k2-1)*ni2nj2;
            dx = xt1[ind1]-xt2[ind2]; dy = yt1[ind1]-yt2[ind2]; dz = zt1[ind1]-zt2[ind2];  
            if (dx*dx + dy*dy + dz*dz < eps2) 
            {
              jm1 = K_FUNC::E_min(j1, jm1); km1 = K_FUNC::E_min(k1,km1); 
              jp1 = K_FUNC::E_max(j1, jp1); kp1 = K_FUNC::E_max(k1,kp1); 
              jm2 = K_FUNC::E_min(j2, jm2); km2 = K_FUNC::E_min(k2,km2); 
              jp2 = K_FUNC::E_max(j2, jp2); kp2 = K_FUNC::E_max(k2,kp2);     
              goto next1;
            }
          }
        next1:;
      }
    if (jp1 == -1) return 0;
    if (jp1 != jm1 && kp1 != km1 && dim == 3) return 1;
    else if ((jp1 != jm1 || kp1 != km1) && dim == 2) return 1;
  }
  else if (nowin1 == 3 || nowin1 == 4) //faces ymin1 avec ymax2
  {
    if (nowin1 == 3)
    {inc1 = 0; inc2 = (nj2-1)*ni2; jm1 = 1; jp1 = 1; jm2 = nj2; jp2 = nj2;}
    else 
    {inc1 = (nj1-1)*ni1; inc2 = 0; jm1 = nj1; jp1 = nj1; jm2 = 1; jp2 = 1;}
    for (E_Int k1 = 1; k1 <= nk1; k1++)
      for (E_Int i1 = 1; i1 <= ni1; i1++)
      {
        ind1 = i1-1 + inc1 + (k1-1)*ni1nj1;
        for (E_Int k2 = 1; k2 <= nk2; k2++)
          for (E_Int i2 = 1; i2 <= ni2; i2++)
          {       
            ind2 = (i2-1) + inc2 + (k2-1)*ni2nj2;
            dx = xt1[ind1]-xt2[ind2]; dy = yt1[ind1]-yt2[ind2]; dz = zt1[ind1]-zt2[ind2];  
            if ( dx*dx + dy*dy + dz*dz < eps2 ) 
            {
              im1 = K_FUNC::E_min(i1,im1); km1 = K_FUNC::E_min(k1,km1); 
              ip1 = K_FUNC::E_max(i1,ip1); kp1 = K_FUNC::E_max(k1,kp1); 
              im2 = K_FUNC::E_min(i2,im2); km2 = K_FUNC::E_min(k2,km2); 
              ip2 = K_FUNC::E_max(i2,ip2); kp2 = K_FUNC::E_max(k2,kp2);     
              goto next2;
            }
          }
        next2:;
      }
    if (ip1 == -1) return 0;
    if (ip1 != im1 && kp1 != km1 && dim == 3) return 1;
    else if ( (ip1 != im1 || kp1 != km1) && dim == 2) return 1;
  }
  else //faces zmin1 avec zmax2
  {
    if (nowin1 == 5) 
    {inc1 = 0; inc2 = (nk2-1)*ni2nj2; km1 = 1; kp1 = 1; km2 = nk2; kp2 = nk2;}
    else
    {inc1 = (nk1-1)*ni1nj1; inc2 = 0; km1 = nk1; kp1 = nk1; km2 = 1; kp2 = 1;}
    for (E_Int j1 = 1; j1 <= nj1; j1++)
      for (E_Int i1 = 1; i1 <= ni1; i1++)
      {
        ind1 = i1-1 + (j1-1)*ni1 + inc1;
        for (E_Int j2 = 1; j2 <= nj2; j2++)
          for (E_Int i2 = 1; i2 <= ni2; i2++)
          {       
            ind2 = (i2-1) + (j2-1)*ni2 + inc2;
            dx = xt1[ind1]-xt2[ind2]; dy = yt1[ind1]-yt2[ind2]; dz = zt1[ind1]-zt2[ind2];  
            if (dx*dx + dy*dy + dz*dz < eps2) 
            {
              im1 = K_FUNC::E_min(i1,im1); jm1 = K_FUNC::E_min(j1,jm1); 
              ip1 = K_FUNC::E_max(i1,ip1); jp1 = K_FUNC::E_max(j1,jp1); 
              im2 = K_FUNC::E_min(i2,im2); jm2 = K_FUNC::E_min(j2,jm2); 
              ip2 = K_FUNC::E_max(i2,ip2); jp2 = K_FUNC::E_max(j2,jp2);  
              goto next3;
            }
          }
        next3:;
      }
    if (ip1 == -1) return 0;
    if (ip1 != im1 && jp1 != jm1 && dim == 3) return 1;
    else if ((ip1 != im1 || jp1 != jm1) && dim == 2) return 1;
  }
  return 0;
}
