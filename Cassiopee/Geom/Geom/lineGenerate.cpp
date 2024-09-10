/*    
    Copyright 2013-2024 Onera.

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

using namespace K_FLD;

// ============================================================================
/* Generate a 2D surface mesh defined by an array from 1D line mesh defined by
   array and a line generator arrayLine */
// ============================================================================
PyObject* K_GEOM::lineGenerateMesh(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* arrayLine;
  if (!PyArg_ParseTuple(args, "OO", &array, &arrayLine)) return NULL;
  
  // Check array
  E_Int im1, jm1, km1, im2, jm2, km2;
  E_Int im3, jm3, km3, im3jm3, im1jm1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  FldArrayF* f2; FldArrayI* cn2;
  char* varString2; char* eltType2;
  
  // Driving curve
  E_Int res2 = K_ARRAY::getFromArray(arrayLine, varString2, 
                                     f2, im2, jm2, km2, cn2, eltType2);

  if (res2 == 2)
  {
    delete f2; delete cn2;
    PyErr_SetString(PyExc_TypeError,
                    "lineGenerate: driving curve must be structured.");
    return NULL;
  }
  if (res2 != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "lineGenerate: driving curve array is invalid.");
    return NULL;
  }

  E_Int res1 = K_ARRAY::getFromArray(array, varString1, f1, im1, jm1, km1, 
                                     cn1, eltType1);

  E_Int nfld = f1->getNfld();

  if (res1 == 1) // structured
  {
    E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
    E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
    E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
    E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
    E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
    E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
    if (posx1 == -1 || posy1 == -1 || posz1 == -1 ||
        posx2 == -1 || posy2 == -1 || posz2 == -1)
    {
      PyErr_SetString(PyExc_TypeError,"lineGenerate: coordinates not found.");
      delete f1; delete f2; return NULL;
    }
    posx1++; posy1++; posz1++;
    posx2++; posy2++; posz2++;

    // L'idee est de remplacer la dimension egale a 1
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
    im1jm1 = im1*jm1;
    
    FldArrayF* coord = new FldArrayF(im3jm3*km3, nfld);
    E_Float* xt = coord->begin(posx1);
    E_Float* yt = coord->begin(posy1);
    E_Float* zt = coord->begin(posz1);
    E_Int ind, ind2, indm;

    E_Float* xt1 = f1->begin(posx1);
    E_Float* yt1 = f1->begin(posy1);
    E_Float* zt1 = f1->begin(posz1);
    E_Float* xt2 = f2->begin(posx2);
    E_Float* yt2 = f2->begin(posy2);
    E_Float* zt2 = f2->begin(posz2);

    // Reordonne la driving curve si necessaire
    E_Int n1 = f1->getSize();
    E_Int n2 = f2->getSize();
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
      K_CONNECT::reorderStructField(im2, jm2, km2, *f2, -1, 2, 3);

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

    // Loop following line
    if (im1 == 1)
    {
      for (E_Int i = 1; i < im3; i++)
        for (E_Int k = 0; k < km3; k++)
          for (E_Int j = 0; j < jm3; j++)
          {
            ind = i + j*im3 + k*im3jm3;
            indm = j*im3 + k*im3jm3; 
            xt[ind] = xt1[indm] + xt2[i] - xt2[0];
            yt[ind] = yt1[indm] + yt2[i] - yt2[0];
            zt[ind] = zt1[indm] + zt2[i] - zt2[0];
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
            xt[ind] = xt1[indm] + xt2[j] - xt2[0];
            yt[ind] = yt1[indm] + yt2[j] - yt2[0];
            zt[ind] = zt1[indm] + zt2[j] - zt2[0];
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
            xt[ind] = xt1[indm] + xt2[k] - xt2[0];
            yt[ind] = yt1[indm] + yt2[k] - yt2[0];
            zt[ind] = zt1[indm] + zt2[k] - zt2[0];
          }
    }
    delete f1; delete f2;
    PyObject* tpl = K_ARRAY::buildArray(*coord, varString1, 
                                        im3, jm3, km3);
    delete coord;
    return tpl; 
  }
  else if (res1 == 2) // unstructured
  {
    char eltType[8];
    if (strcmp(eltType1, "BAR") == 0) strcpy(eltType, "QUAD");
    else if (strcmp(eltType1, "QUAD") == 0) strcpy(eltType, "HEXA");
    else if (strcmp(eltType1, "TRI") == 0) strcpy(eltType, "PENTA");
    else
    {
      delete f2; delete f1; delete cn1;
      PyErr_SetString(PyExc_TypeError,
                      "lineGenerate: array must be structured, BAR, TRI or QUAD.");
      return NULL;
    }

    E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
    E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
    E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
    E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
    E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
    E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
    if ( posx1 == -1 || posy1 == -1 || posz1 == -1 ||
         posx2 == -1 || posy2 == -1 || posz2 == -1)
    {
      PyErr_SetString(PyExc_TypeError,"lineGenerate: coordinates not found.");
      delete f2; delete f1; delete cn1; return NULL;
    }
    posx1++; posy1++; posz1++;
    posx2++; posy2++; posz2++;
    
    FldArrayF* coord = new FldArrayF(f1->getSize()*im2, nfld);
    E_Float* xt = coord->begin(posx1);
    E_Float* yt = coord->begin(posy1);
    E_Float* zt = coord->begin(posz1);

    E_Float* xt1 = f1->begin(posx1);
    E_Float* yt1 = f1->begin(posy1);
    E_Float* zt1 = f1->begin(posz1);
    E_Float* xt2 = f2->begin(posx2);
    E_Float* yt2 = f2->begin(posy2);
    E_Float* zt2 = f2->begin(posz2);
    E_Int ind, ind2;

    // Reordonne la driving curve si necessaire
    E_Int n1 = f1->getSize();
    E_Int n2 = f2->getSize();
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
      K_CONNECT::reorderStructField(im2, jm2, km2, *f2, -1, 2, 3);

    // Init (pour les champs)
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* ft = coord->begin(n);
      E_Float* ft1 = f1->begin(n);
      
      for (E_Int i = 0; i < im2; i++)
      {
        for (E_Int j = 0; j < n1; j++)
        {
          ind = j;
          ind2 = j+n1*i;
          ft[ind2] = ft1[ind];
        }
      }
    }

    // Init connectivite
    E_Int nelt1 = cn1->getSize();
    E_Int net = cn1->getNfld();
    FldArrayI* connect = new FldArrayI(nelt1*(im2-1), net*2);
    FldArrayI& cn = *connect;
    FldArrayI& cnn1 = *cn1;
    for (E_Int i = 0; i < im2-1; i++)
    {
      for (E_Int e = 1; e <= net; e++)
      {
        for (E_Int j = 0; j < nelt1; j++)
        {
          cn(j+nelt1*i, e) = cnn1(j, e) + i*n1;
          cn(j+nelt1*i, net+e) = cnn1(j, e) + (i+1)*n1;
        }
      }
    }

    // Reorder pour les quads
    E_Int temp;
    if (strcmp(eltType, "QUAD") == 0)
    {
      for (E_Int j = 0; j < cn.getSize(); j++)
      {
        temp = cn(j, 3);
        cn(j, 3) = cn(j, 4);
        cn(j, 4) = temp;
      }
    }

    // Modification des coordonnees
    for (E_Int i = 0; i < im2; i++)
    {
      for (E_Int j = 0; j < n1; j++)
      {
        xt[j + i*n1] = xt1[j] + xt2[i] - xt2[0];
        yt[j + i*n1] = yt1[j] + yt2[i] - yt2[0];
        zt[j + i*n1] = zt1[j] + zt2[i] - zt2[0];
      }
    }

    delete f2; delete f1; delete cn1;
    PyObject* tpl = K_ARRAY::buildArray(*coord, varString1, 
                                        *connect, -1, eltType);
    delete coord; delete connect;
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "lineGenerate: invalid type of array.");
    return NULL;
  }
}
