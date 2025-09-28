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

# include "kcore.h"
# include "geom.h"
# include <stdio.h>

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

//=============================================================================
/* Add a separation line to a geometry */
//=============================================================================
PyObject* K_GEOM::addSeparationLineMesh(PyObject* self, PyObject* args)
{
  E_Float d, x, y, z, x2, y2, z2;
  E_Int i, ind;
  PyObject* array1; PyObject* array2;
  if (!PYPARSETUPLE_(args, OO_, &array1, &array2)) return NULL;

  // Check array
  E_Int im,im2,jm,jm2,km,km2;
  FldArrayF* f1; FldArrayF* f2;
  FldArrayI* cn1; FldArrayI* cn2;
  char* varString1; char* varString2;
  char* eltType1; char* eltType2;

  E_Int res1 = K_ARRAY::getFromArray3(array1, varString1, f1,
                                      im, jm, km, cn1, eltType1);
  if (res1 != 1 && res1 != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "addSeparationLine: unknown type of array.");
    return NULL;
  }
  else if (res1 == 2)
  {
    RELEASESHAREDB(res1, array1, f1, cn1);
    PyErr_SetString(PyExc_TypeError,
                    "addSeparationLine: array must be structured.");
    return NULL;
  }

  E_Int res2 = K_ARRAY::getFromArray3(array2, varString2, f2,
                                      im2, jm2, km2, cn2, eltType2);
  if (res2 != 1 && res2 != 2)
  {
    RELEASESHAREDB(res1, array1, f1, cn1);
    PyErr_SetString(PyExc_TypeError,
                    "addSeparationLine: unknown type of array.");
    return NULL;
  }
  else if (res2 == 2)
  {
    RELEASESHAREDB(res1, array1, f1, cn1);
    RELEASESHAREDB(res2, array2, f2, cn2);
    PyErr_SetString(PyExc_TypeError,
                    "addSeparationLine: array must be structured.");
    return NULL;
  }

  FldArrayF* coord0 = new FldArrayF();
  FldArrayF& coord = *coord0;
  FldArrayF* coord1 = new FldArrayF();
  FldArrayF& coord2 = *coord1;
  E_Int npt1, npt2;
  vector<E_Int> pos1;
  vector<E_Int> pos2;

  char* varString = new char [strlen(varString1)+strlen(varString2)+4];
  E_Int res0 = K_ARRAY::getPosition(varString1, varString2, pos1, pos2, varString);

  if (res0 == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "addSeparationLine: coordinates not found in array.");
    return NULL;
  }

  if (jm != 1 && km != 1 && jm2 != 1 && km2 != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "addSeparationLine: array must be i-array.");
    return NULL;
  }

  E_Int nfld = pos1.size();
  // Find the nearest point of two geometries
  E_Float dmax = 1.e6;
  E_Int isav = -1;
  E_Int isav2 = -1;

  // Extremite de la ligne
  x2 = (*f2)(0,pos2[0]);
  y2 = (*f2)(0,pos2[1]);
  z2 = (*f2)(0,pos2[2]);

  for (i = 0; i < im; i++)
  {
    x = (*f1)(i,pos1[0]);
    y = (*f1)(i,pos1[1]);
    z = (*f1)(i,pos1[2]);

    d = (x - x2)*(x - x2) + (y - y2)*(y - y2) + (z - z2)*(z - z2);
    if (d < dmax)
    {
      isav = i;
      dmax = d;
    }
  }

  x2 = (*f2)(im2-1,pos2[0]);
  y2 = (*f2)(im2-1,pos2[1]);
  z2 = (*f2)(im2-1,pos2[2]);

  for (i = 0; i < im; i++)
  {
    x = (*f1)(i,pos1[0]);
    y = (*f1)(i,pos1[1]);
    z = (*f1)(i,pos1[2]);

    d = (x - x2)*(x - x2) + (y - y2)*(y - y2) + (z - z2)*(z - z2);
    if (d < dmax)
    {
      isav2 = i;
      dmax = d;
    }
  }

  if (E_abs(dmax) > 1.e-4)
  {
    printf("Warning: addSeparationLine: msh2 is not stuck on mesh...");
    printf(" Sticking forced.\n");
  }

  // Formation du ou des nouveaux maillages
  coord.malloc(im+2*im2, nfld);
  coord2.malloc(im+2*im2, nfld);
  ind = 0;
  npt1 = 0; npt2 = 0;

  if (isav2 == -1 && (isav == 0 || isav == im-1)) // Extremite 1 est sur msh
  {
    for (i = im2-1; i >= 0; i--)
    {
      for (E_Int n = 1; n <= nfld; n++)
        coord(ind,n) = (*f2)(i,pos2[n-1]);
      ind++;
    }
    for (i = 0; i < im; i++)
    {
      for (E_Int n = 1; n <= nfld; n++)
        coord(ind,n) = (*f1)(i,pos1[n-1]);
      ind++;
    }
    for (i = 0; i < im2; i++)
    {
      for (E_Int n = 1; n <= nfld; n++)
        coord(ind,n) = (*f2)(i,pos2[n-1]);
      ind++;
    }
    npt1 = ind;
  }

  else if (isav2 == -1)  // Extremite 1 est sur msh
  {
    for (i = im2-1; i >= 0; i--)
    {
      for (E_Int n = 1; n <= nfld; n++)
        coord(ind,n) = (*f2)(i,pos2[n-1]);
      ind++;
    }
    for (i = isav; i < im; i++)
    {
      for (E_Int n = 1; n <= nfld; n++)
        coord(ind,n) = (*f1)(i,pos1[n-1]);
      ind++;
    }
    npt1 = ind;
    ind = 0;
    for (i = im2-1; i >= 0; i--)
    {
      for (E_Int n = 1; n <= nfld; n++)
        coord2(ind, n) = (*f2)(i,pos2[n-1]);
      ind++;
    }
    for (i = isav; i >= 0; i--)
    {
      for (E_Int n = 1; n <= nfld; n++)
        coord2(ind,n) = (*f1)(i,pos1[n-1]);
      ind++;
    }
    npt2 = ind;
  }

  else if (isav2 == 0 || isav2 == im-1) // Extremite 2 est sur msh
  {
    for (i = 0; i < im2; i++)
    {
      for (E_Int n = 1; n <= nfld; n++)
        coord(ind,n) = (*f2)(i,pos2[n-1]);
      ind++;
    }
    for (i = 0; i < im; i++)
    {
      for (E_Int n = 1; n <= nfld; n++)
        coord(ind, n) = (*f1)(i,pos1[n-1]);
      ind++;
    }
    for (i = im2-1; i >= 0; i--)
    {
      for (E_Int n = 1; n <= nfld; n++)
        coord(ind,n) = (*f2)(i,pos2[n-1]);
      ind++;
    }
    npt1 = ind;
  }
  else
  {
    for (i = 0; i < im2; i++)
    {
      for (E_Int n = 1; n <= nfld; n++)
        coord(ind,n) = (*f2)(i,pos2[n-1]);
      ind++;
    }
    for (i = isav; i < im; i++)
    {
      for (E_Int n = 1; n <= nfld; n++)
        coord(ind, n) = (*f1)(i,pos1[n-1]);
      ind++;
    }
    npt1 = ind;
    ind = 0;
    for (i = 0; i < im2; i++)
    {
      for (E_Int n = 1; n <= nfld; n++)
        coord2(ind,n) = (*f2)(i,pos2[n-1]);
      ind++;
    }
    for (i = isav; i < im; i++)
    {
      for (E_Int n = 1; n <= nfld; n++)
        coord2(ind,n) = (*f1)(i,pos1[n-1]);
      ind++;
    }
    npt2 = ind;
  }

  coord.reAllocMat(npt1, nfld);
  coord2.reAllocMat(npt2, nfld);

  K_CONNECT::supIdPoints(coord, pos1[0], pos1[1], pos1[2]);
  K_CONNECT::supIdPoints(coord2,pos1[0], pos1[1], pos1[2]);

  npt1 = coord.getSize();
  npt2 = coord2.getSize();

  PyObject* l = PyList_New(0);
  E_Int api = f1->getApi();

  if (npt1 > 0)
  {
    PyObject* tpl = K_ARRAY::buildArray3(*coord0, varString, npt1, 1, 1, api);
    delete coord0;
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }

  if (npt2 > 0)
  {
    PyObject* tpl = K_ARRAY::buildArray3(*coord1, varString, npt2, 1, 1, api);
    delete coord1;
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  delete [] varString;
  RELEASESHAREDS(array1, f1);
  RELEASESHAREDS(array2, f2);
  return l;
}
