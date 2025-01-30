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
#include "Array/Array.h"
#include "String/kstring.h"
#include <stdio.h>
#include <string.h>

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Build a structured array 
   IN: field: Fld structured field
   IN: varString: variables string
   IN: ni,nj,nk: number of points in field
   OUT: PyObject created. */
//=============================================================================
PyObject* K_ARRAY::buildArray(FldArrayF& field, const char* varString, 
                              E_Int ni, E_Int nj, E_Int nk)
{
  npy_intp dim[2];
  PyObject* tpl;
  IMPORTNUMPY;

  dim[1] = field.getSize();
  dim[0] = field.getNfld();
  PyArrayObject* a = (PyArrayObject*)PyArray_SimpleNew(2, dim, NPY_DOUBLE);
  //memcpy(PyArray_DATA(a), field.begin(), dim[0]*dim[1]*sizeof(double));

  E_Float* ptrf = field.begin();
  E_Float* ptra = (E_Float*)PyArray_DATA(a);
  E_Int s = dim[0]*dim[1];
#pragma omp parallel for default(shared)
  for (E_Int i = 0; i < s; i++) ptra[i] = ptrf[i];

  tpl = Py_BuildValue("[sOlll]", varString, a,  
                      (long)ni, (long)nj, (long)nk);
  Py_DECREF(a);
  return tpl;
}

//=============================================================================
/* Build an empty structured array 
   IN: nfld: nbre de champs
   IN: varString: variables string
   IN: ni,nj,nk: number of points in field
   OUT: PyObject created. */
//=============================================================================
PyObject* K_ARRAY::buildArray(E_Int nfld, const char* varString, 
                              E_Int ni, E_Int nj, E_Int nk)
{
  npy_intp dim[2];
  PyObject* tpl;
  IMPORTNUMPY;

  dim[1] = ni*nj*nk;
  dim[0] = nfld;
  PyArrayObject* a = (PyArrayObject*)PyArray_SimpleNew(2, dim, NPY_DOUBLE);
  tpl = Py_BuildValue("[sOlll]", varString, a,  
                      (long)ni, (long)nj, (long)nk);
  Py_DECREF(a);
  return tpl;
}

//=============================================================================
/* Build a unstructured array 
   IN: field: Fld unstructured field
   IN: varString: variable string
   IN: c: Fld connectivity
   IN: et: elements type: 0 (NODE), 1 (BAR), 2 (TRI), 3 (QUAD),
                          4 (TETRA), 5 (PYRA), 6 (PENTA), 7 (HEXA),
                          8 (NGON).
   IN: etString: if et=-1, etString is used instead.
   IN: center: set to true if field is localised in the centers of
   elements, otherwise let it to false.
   
   OUT: PyObject created. */
//=============================================================================
PyObject* K_ARRAY::buildArray(FldArrayF& field, const char* varString,
                              FldArrayI& c, E_Int et, const char* etString,
                              E_Boolean center)
{
  npy_intp dim[2];
  PyObject* tpl;
  char eltType[8];

  E_Int fSize = field.getSize();
  E_Int cSize = c.getSize(); // taille de la connectivite
  E_Int cEltSize = cSize; // nb d'elements dans la connectivite
  if (et == 8 || (et == -1 && K_STRING::cmp(etString, "NGON") == 0))
  {
    E_Int sizeFN = c[1]; cEltSize = c[sizeFN+2];
  }

  IMPORTNUMPY;

  if (et == -1)
  {
    if (etString != NULL)
    {
      strcpy(eltType, etString);
      E_Int pos = strlen(eltType)-1;
      if (cEltSize == fSize && eltType[pos] != '*' && center == true)
        strcat(eltType, "*");
      else if (cEltSize != fSize && eltType[pos] == '*')
        eltType[pos] = '\0';
    }
    else
    {
      printf("Warning: buildArray: element type is unknown.\n");
      strcpy(eltType, "UNKNOWN");
    }
  }
  else
  {
    switch (et)
    {
      case 0:
        strcpy(eltType, "NODE");
        if (cEltSize == fSize && center == true) strcat(eltType, "*");
        break;
      case 1:
        strcpy(eltType, "BAR");
        if (cEltSize == fSize && center == true) strcat(eltType, "*");
        break;
      case 2:
        strcpy(eltType, "TRI");
        if (cEltSize == fSize && center == true) strcat(eltType, "*");
        break;
      case 3:
        strcpy(eltType, "QUAD");
        if (cEltSize == fSize && center == true) strcat(eltType, "*");
        break;
      case 4:
        strcpy(eltType, "TETRA");
        if (cEltSize == fSize && center == true) strcat(eltType, "*");
        break;
      case 5:
        strcpy(eltType, "PYRA");
        if (cEltSize == fSize && center == true) strcat(eltType, "*");
        break;
      case 6:
        strcpy(eltType, "PENTA");
        if (cEltSize == fSize && center == true) strcat(eltType, "*");
        break;
      case 7:
        strcpy(eltType, "HEXA");
        if (cEltSize == fSize && center == true) strcat(eltType, "*");
        break;
      case 8:
        strcpy(eltType, "NGON");
        if (cEltSize == fSize && center == true) strcat(eltType, "*");
        break;
      default:
        printf("Warning: buildArray: element type is unknown.\n");
        strcpy(eltType, "UNKNOWN");
        break;
    }
  }
    
  // Build array of field 
  dim[1] = fSize;
  dim[0] = field.getNfld();
  E_Int s1 = dim[0]*dim[1];
  PyArrayObject* a = (PyArrayObject*)PyArray_SimpleNew(2, dim, NPY_DOUBLE);

  // Build array for connectivity   
  dim[1] = cSize;
  dim[0] = c.getNfld();
  E_Int s2 = dim[0]*dim[1];
  PyArrayObject* ac = (PyArrayObject*)PyArray_SimpleNew(2, dim, E_NPY_INT);

  E_Float* ptrf = field.begin();
  E_Float* ptra = (E_Float*)PyArray_DATA(a);
  E_Int* ptrc = c.begin();
  E_Int* ptrac = (E_Int*)PyArray_DATA(ac);
#pragma omp parallel default(shared)
  {
#pragma omp for nowait
    for (E_Int i = 0; i < s1; i++) ptra[i] = ptrf[i];

#pragma omp for
    for (E_Int i = 0; i < s2; i++) ptrac[i] = ptrc[i];
  }

  tpl = Py_BuildValue("[sOOs]", varString, a, ac, eltType);
  Py_DECREF(a); Py_DECREF(ac);

  return tpl;
}

//=============================================================================
/* Build an empty unstructured array 
   IN: nfld: number of fields
   IN: varString: variable string
   IN: nvertex: number of vertex
   IN: nelt: number of elements
   IN: et: elements type: 0 (NODE), 1 (BAR), 2 (TRI), 3 (QUAD),
                          4 (TETRA), 5 (PYRA), 6 (PENTA), 7 (HEXA),
                          8 (NGON).
   IN: etString: if et=-1, etString is used instead.
   IN: center: set to true if field is localised in the centers of
   elements, otherwise let it to false.
   IN: sizeConnect: connectivity size. Used only for NGONS. 

   OUT: PyObject created. */
//=============================================================================
PyObject* K_ARRAY::buildArray(E_Int nfld, const char* varString,
                              E_Int nvertex, E_Int nelt,
                              E_Int et, const char* etString,
                              E_Boolean center, E_Int sizeConnect)
{
  npy_intp dim[2];
  PyObject* a;
  PyObject* ac;
  PyObject* tpl;
  char eltType[8];

  E_Int fSize;
  if (center == true) fSize = nelt;
  else fSize = nvertex;
  E_Int cSize = nelt;
  E_Int nvpe = 0;
 
  IMPORTNUMPY;

  if (et == -1)
  {
    if (etString != NULL)
    {
      strcpy(eltType, etString);
      E_Int pos = strlen(eltType)-1;
      pos = 0;
      if (cSize == fSize && eltType[pos] != '*' && center == true)
        strcat(eltType, "*");
      else if (cSize != fSize && eltType[pos] == '*')
        eltType[pos] = '\0';
    }
    else
    {
      printf("Warning: buildArray: element type is unknown.\n");
      strcpy(eltType, "UNKNOWN");
    }
    if (K_STRING::cmp(eltType, "NODE") == 0) nvpe = 1;
    else if (K_STRING::cmp(eltType, "TRI") == 0) nvpe = 3;
    else if (K_STRING::cmp(eltType, "QUAD") == 0) nvpe = 4;
    else if (K_STRING::cmp(eltType, "TETRA") == 0) nvpe = 4;
    else if (K_STRING::cmp(eltType, "HEXA") == 0) nvpe = 8;
    else if (K_STRING::cmp(eltType, "BAR") == 0) nvpe = 2;
    else if (K_STRING::cmp(eltType, "PYRA") == 0) nvpe = 5;
    else if (K_STRING::cmp(eltType, "PENTA") == 0) nvpe = 6;
    else if (K_STRING::cmp(eltType, "NGON") == 0) { nvpe = 1; cSize = sizeConnect; }
    else if (K_STRING::cmp(eltType, "NODE*") == 0) nvpe = 1;
    else if (K_STRING::cmp(eltType, "TRI*") == 0) nvpe = 3;
    else if (K_STRING::cmp(eltType, "QUAD*") == 0) nvpe = 4;
    else if (K_STRING::cmp(eltType, "TETRA*") == 0) nvpe = 4;
    else if (K_STRING::cmp(eltType, "HEXA*") == 0) nvpe = 8;
    else if (K_STRING::cmp(eltType, "BAR*") == 0) nvpe = 2;
    else if (K_STRING::cmp(eltType, "PYRA*") == 0) nvpe = 5;
    else if (K_STRING::cmp(eltType, "PENTA*") == 0) nvpe = 6;
    else if (K_STRING::cmp(eltType, "NGON*") == 0) { nvpe = 1; cSize = sizeConnect; }
  }
  else
  {
    switch (et)
    {
      case 0:
        strcpy(eltType, "NODE");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 1;
        break;
      case 1:
        strcpy(eltType, "BAR");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 2;
        break;
      case 2:
        strcpy(eltType, "TRI");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 3;
        break;
      case 3:
        strcpy(eltType, "QUAD");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 4;
        break;
      case 4:
        strcpy(eltType, "TETRA");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 4;
        break;
      case 5:
        strcpy(eltType, "PYRA");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 5;
        break;
      case 6:
        strcpy(eltType, "PENTA");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 6;
        break;
      case 7:
        strcpy(eltType, "HEXA");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 8;
        break;
      case 8:
        strcpy(eltType, "NGON");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 1; cSize = sizeConnect;
        break;
      default:
        printf("Warning: buildArray: element type is unknown.\n");
        strcpy(eltType, "UNKNOWN");
        break;
    }
  }
  
  // Build array of field 
  dim[1] = fSize; dim[0] = nfld;
  a = PyArray_SimpleNew(2, dim, NPY_DOUBLE);

  // Build array for connectivity
  dim[1] = cSize; dim[0] = nvpe;
  ac = PyArray_SimpleNew(2, dim, E_NPY_INT);
  tpl = Py_BuildValue("[sOOs]", varString, a, ac, eltType);
  Py_DECREF(a); Py_DECREF(ac);

  return tpl;
}
