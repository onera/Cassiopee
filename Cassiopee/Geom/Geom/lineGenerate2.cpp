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

// Generator of surface geometries

#include "geom.h"
#include <vector>
using namespace std;
using namespace K_FLD;

// ============================================================================
/* Generate a 2D surface mesh defined by an array from 
   a structured 1D line
   and a set of structured lines */
// ============================================================================
PyObject* K_GEOM::lineGenerate2(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* dCurves;
  if (!PyArg_ParseTuple(args, "OO", &array, &dCurves)) return NULL;
  
  // Driving curves
  vector<E_Int> res;
  vector<char*> structVarString;
  vector<char*> unstructVarString;
  vector<FldArrayF*> structF;
  vector<FldArrayF*> unstructF;
  vector<E_Int> ni; vector<E_Int> nj; vector<E_Int> nk;
  vector<FldArrayI*> c;
  vector<char*> eltType;
  vector<PyObject*> objs; vector<PyObject*> obju;
  K_ARRAY::getFromArrays(dCurves, res, structVarString,
                         unstructVarString, structF, unstructF,
                         ni, nj, nk, c, eltType, objs, obju,
                         true, true, false, true, true);

  E_Int structSize = structF.size();
  if (structSize == 0)
  {
    for (size_t i = 0; i < objs.size(); i++) RELEASESHAREDS(objs[i], structF[i]);
    for (size_t i = 0; i < obju.size(); i++) RELEASESHAREDU(obju[i], unstructF[i], c[i]);
    PyErr_SetString(PyExc_TypeError,
                    "lineGenerate2: no valid driving curves.");
    return NULL;
  }

  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int im1, jm1, km1;
  E_Int res1 = K_ARRAY::getFromArray3(array, varString1, f1, im1, jm1, km1, 
                                      cn1, eltType1);

  if (res1 == 2)
  {
    for (size_t i = 0; i < objs.size(); i++) RELEASESHAREDS(objs[i], structF[i]);
    for (size_t i = 0; i < obju.size(); i++) RELEASESHAREDU(obju[i], unstructF[i], c[i]);
    RELEASESHAREDU(array, f1, cn1);
    PyErr_SetString(PyExc_TypeError,
                    "lineGenerate2: array must be structured.");
    return NULL;
  }
  
  // Main curve
  E_Int nfld = f1->getNfld();
  E_Int n1 = f1->getSize();
  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
  if (posx1 == -1 || posy1 == -1 || posz1 == -1)
  {
    PyErr_SetString(PyExc_TypeError,"lineGenerate2: coordinates not found.");
    for (size_t i = 0; i < objs.size(); i++) RELEASESHAREDS(objs[i], structF[i]);
    for (size_t i = 0; i < obju.size(); i++) RELEASESHAREDU(obju[i], unstructF[i], c[i]);
    RELEASESHAREDS(array, f1);
    return NULL;
  }
  posx1++; posy1++; posz1++;
  E_Float* xt1 = f1->begin(posx1);
  E_Float* yt1 = f1->begin(posy1);
  E_Float* zt1 = f1->begin(posz1);
  E_Int im1jm1 = im1*jm1;

  // Driving curves
  E_Int nd = structF.size();
  
  E_Int im2 = ni[0]; // reference
  E_Int jm2 = nj[0];
  E_Int km2 = nk[0];
  E_Int n2 = structF[0]->getSize();

  E_Int posx2 = K_ARRAY::isCoordinateXPresent(structVarString[0]);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(structVarString[0]);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(structVarString[0]);
  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    PyErr_SetString(PyExc_TypeError,"lineGenerate2: coordinates not found.");
    for (size_t i = 0; i < objs.size(); i++) RELEASESHAREDS(objs[i], structF[i]);
    for (size_t i = 0; i < obju.size(); i++) RELEASESHAREDU(obju[i], unstructF[i], c[i]);
    RELEASESHAREDS(array, f1);
    return NULL;
  }
  posx2++; posy2++; posz2++;

  // L'idee est de remplacer la dimension egale a 1
  E_Int im3, jm3, km3, im3jm3;
  if (im1 == 1)
  {
    im3 = im2; jm3 = jm1; km3 = km1;
  }
  else if (jm1 == 1)
  {
    im3 = im1; jm3 = im2; km3 = km1;
  }
  else
  {
    im3 = im1; jm3 = jm1; km3 = im2;
  }
  
  im3jm3 = im3*jm3;
    
  FldArrayF* coord = new FldArrayF(im3jm3*km3, nfld);
  E_Float* xt = coord->begin(posx1);
  E_Float* yt = coord->begin(posy1);
  E_Float* zt = coord->begin(posz1);
    
  E_Int ind, ind2, indm;
  FldArrayI firstCurve(im3); FldArrayF firstDist(im3);
  FldArrayI secondCurve(im3); FldArrayF secondDist(im3);
  firstCurve.setAllValuesAt(-1); secondCurve.setAllValuesAt(-1);
  firstDist.setAllValuesAt(1.e6); secondDist.setAllValuesAt(1.e6);

  // Reordonne les driving curves si necessaire
  vector <E_Float*> dxt(nd);
  vector <E_Float*> dyt(nd);
  vector <E_Float*> dzt(nd);
  for (E_Int d = 0; d < nd; d++)
  {
    E_Int posx2 = K_ARRAY::isCoordinateXPresent(structVarString[d]); posx2++;
    E_Int posy2 = K_ARRAY::isCoordinateYPresent(structVarString[d]); posy2++;
    E_Int posz2 = K_ARRAY::isCoordinateZPresent(structVarString[d]); posz2++;
    E_Float* xt2 = structF[d]->begin(posx2);
    E_Float* yt2 = structF[d]->begin(posy2);
    E_Float* zt2 = structF[d]->begin(posz2);

    E_Float xa = xt2[0]; E_Float ya = yt2[0]; E_Float za = zt2[0];
    E_Float xb = xt2[n2-1]; E_Float yb = yt2[n2-1]; E_Float zb = zt2[n2-1];
    E_Float d1 = 1.e6; E_Float d2 = 1.e6;
    E_Float dx, dy, dz;
    for (E_Int i = 0; i < n1; i++)
    {
      dx = xt1[i]-xa; dy = yt1[i]-ya; dz = zt1[i]-za;
      d1 = K_FUNC::E_min(d1, dx*dx+dy*dy+dz*dz);
      dx = xt1[i]-xb; dy = yt1[i]-yb; dz = zt1[i]-zb;
      d2 = K_FUNC::E_min(d2, dx*dx+dy*dy+dz*dz);
    }
    if (d2 < d1)
      K_CONNECT::reorderStructField(im2, jm2, km2, *structF[d], -1, 2, 3);

    xt2 = structF[d]->begin(posx2);
    yt2 = structF[d]->begin(posy2);
    zt2 = structF[d]->begin(posz2);
    dxt[d] = xt2;  dyt[d] = yt2; dzt[d] = zt2;
    xa = xt2[0]; ya = yt2[0]; za = zt2[0];

    for (E_Int i = 0; i < n1; i++)
    {
      dx = xt1[i]-xa; dy = yt1[i]-ya; dz = zt1[i]-za;
      d1 = dx*dx+dy*dy+dz*dz;
      if (d1 < firstDist[i]) { secondDist[i] = firstDist[i]; secondCurve[i] = firstCurve[i]; firstDist[i] = d1; firstCurve[i] = d; }
      else if (d1 < secondDist[i]) { secondDist[i] = d1; secondCurve[i] = d; }
    }
  }

  // for (E_Int i = 0; i < n1; i++)
  // {
  //   printf("%d %d (%f) %d (%f)\n", i, firstCurve[i], firstDist[i], secondCurve[i], secondDist[i]);
  // }

  // Init (pour les champs)
  for (E_Int n = 1; n <= nfld; n++)
  {
    E_Float* ft = coord->begin(n);
    E_Float* ft1 = f1->begin(n);
    if (im1 == 1)
    {
      for (E_Int k = 0; k < km3; k++)
        for (E_Int j = 0; j < jm3; j++)
          for (E_Int i = 0; i < im3; i++)
          {
            ind = j*im1 + k*im1jm1;
            ind2 = i + j*im3 + k*im3jm3;
            ft[ind2] = ft1[ind];
          }
    }
    else if (jm1 == 1)
    {
      for (E_Int k = 0; k < km3; k++)
        for (E_Int j = 0; j < jm3; j++)
          for (E_Int i = 0; i < im3; i++)
          {
            ind = i + k*im1jm1;
            ind2 = i + j* im3 + k*im3jm3;
            ft[ind2] = ft1[ind];         
          }
    }
    else
    {
      for (E_Int k = 0; k < km3; k++)
        for (E_Int j = 0; j < jm3; j++)
          for (E_Int i = 0; i < im3; i++)
          {
            ind = i + j*im1;
            ind2 = i + j*im3 + k*im3jm3;
            ft[ind2] = ft1[ind];
          }
    }
  }
  
  E_Int c1, c2;
  E_Float alpha, beta, del1x, del1y, del1z, del2x, del2y, del2z, d1, d2;

  // Loop following line
  if (im1 == 1)
  {
    for (E_Int i = 1; i < im3; i++)
      for (E_Int k = 0; k < km3; k++)
        for (E_Int j = 0; j < jm3; j++)
        {
          ind = i + j*im3 + k*im3jm3;
          indm = j*im3 + k*im3jm3;

          c1 = firstCurve[j]; d1 = firstDist[j];
          c2 = secondCurve[j]; d2 = secondDist[j];
          alpha = d1/(d1+d2+1.e-10);
          beta = 1. - alpha;

          del1x = dxt[c1][i] - dxt[c1][0];
          del1y = dyt[c1][i] - dyt[c1][0];
          del1z = dzt[c1][i] - dzt[c1][0];
          
          del2x = dxt[c2][i] - dxt[c2][0];
          del2y = dyt[c2][i] - dyt[c2][0];
          del2z = dzt[c2][i] - dzt[c2][0];
 
          xt[ind]= xt1[indm] + alpha*del2x + beta*del1x;
          yt[ind]= yt1[indm] + alpha*del2y + beta*del1y;
          zt[ind]= zt1[indm] + alpha*del2z + beta*del1z;
        }
  }
  else if (jm1 == 1)
  {
    for (E_Int j = 1; j < jm3; j++)
      for (E_Int k = 0; k < km3; k++)
        for (E_Int i = 0; i < im3; i++)
        {
          ind = i + j*im3 + k*im3jm3;
          indm = i + k*im3jm3;

          c1 = firstCurve[i]; d1 = firstDist[i];
          c2 = secondCurve[i]; d2 = secondDist[i];
          alpha = d1/(d1+d2+1.e-10);
          beta = 1. - alpha;

          del1x = dxt[c1][j] - dxt[c1][0];
          del1y = dyt[c1][j] - dyt[c1][0];
          del1z = dzt[c1][j] - dzt[c1][0];
          
          del2x = dxt[c2][j] - dxt[c2][0];
          del2y = dyt[c2][j] - dyt[c2][0];
          del2z = dzt[c2][j] - dzt[c2][0];

          xt[ind]= xt1[indm] + alpha*del2x + beta*del1x;
          yt[ind]= yt1[indm] + alpha*del2y + beta*del1y;
          zt[ind]= zt1[indm] + alpha*del2z + beta*del1z;
        }
  }
  else
  {
    for (E_Int k = 1; k < km3; k++)
      for (E_Int j = 0; j < jm3; j++)
        for (E_Int i = 0; i < im3; i++)
        {
          ind = i +j*im3 + k*im3jm3;
          indm = i + j*im3;

          c1 = firstCurve[i]; d1 = firstDist[i];
          c2 = secondCurve[i]; d2 = secondDist[i];
          alpha = d1/(d1+d2+1.e-10);
          beta = 1. - alpha;

          del1x = dxt[c1][k] - dxt[c1][0];
          del1y = dyt[c1][k] - dyt[c1][0];
          del1z = dzt[c1][k] - dzt[c1][0];
          
          del2x = dxt[c2][k] - dxt[c2][0];
          del2y = dyt[c2][k] - dyt[c2][0];
          del2z = dzt[c2][k] - dzt[c2][0];

          xt[ind]= xt1[indm] + alpha*del2x + beta*del1x;
          yt[ind]= yt1[indm] + alpha*del2y + beta*del1y;
          zt[ind]= zt1[indm] + alpha*del2z + beta*del1z;
        }
  }
  
  for (size_t i = 0; i < objs.size(); i++) RELEASESHAREDS(objs[i], structF[i]);
  for (size_t i = 0; i < obju.size(); i++) RELEASESHAREDU(obju[i], unstructF[i], c[i]);
  RELEASESHAREDS(array, f1);
  PyObject* tpl = K_ARRAY::buildArray(*coord, varString1, 
                                      im3, jm3, km3);

  delete coord;
  return tpl; 
}
