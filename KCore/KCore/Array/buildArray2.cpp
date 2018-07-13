/*    
    Copyright 2013-2018 Onera.

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
/* Build an empty structured array 
   IN: nfld: nbre de champs
   IN: varString: variables string
   IN: ni,nj,nk: number of points in field
   IN: api (1: array, 2: array2)
   OUT: PyObject created. */
//=============================================================================
PyObject* K_ARRAY::buildArray2(E_Int nfld, const char* varString, 
                               E_Int ni, E_Int nj, E_Int nk, E_Int api)
{
  PyObject* tpl;
  IMPORTNUMPY;

  if (api == 1) // Array
  {
    npy_intp dim[2];  
    dim[1] = ni*nj*nk;
    dim[0] = nfld;
    PyArrayObject* a = (PyArrayObject*)PyArray_SimpleNew(2, dim, NPY_DOUBLE);
    tpl = Py_BuildValue("[sOlll]", varString, a,  
                        (long)ni, (long)nj, (long)nk);
    Py_DECREF(a);
  }
  else // Array2
  {
    npy_intp dim[3]; int ndim=3;
    dim[0] = ni; dim[1] = nj; dim[2] = nk;
    if (nk == 1) ndim--;
    if (nj == 1) ndim--;
    PyObject* rake = PyList_New(0);
    for (E_Int n=0; n < nfld; n++)
    {
      PyArrayObject* a = (PyArrayObject*)PyArray_EMPTY(3, dim, NPY_DOUBLE, 1);
      PyList_Append(rake, (PyObject*)a); Py_DECREF(a);
    }
    tpl = Py_BuildValue("[sOlll]", varString, rake,  (long)ni, (long)nj, (long)nk);
    Py_DECREF(rake);
  }
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
   IN: sizeNGon, sizeNFace, nface: connectivity size. Used only for NGONS. 

   OUT: PyObject created. */
//=============================================================================
PyObject* K_ARRAY::buildArray2(E_Int nfld, const char* varString,
                               E_Int nvertex, E_Int nelt,
                               E_Int et, const char* etString,
                               E_Boolean center, E_Int sizeNGon, 
                               E_Int sizeNFace, E_Int nface, E_Int api)
{
  npy_intp dim[2];
  PyObject* a; PyObject* ac; PyObject* tpl;
  char eltType[8];

  E_Int fSize;
  if (center == true) fSize = nelt;
  else fSize = nvertex;
  E_Int cSize = nelt;
  E_Int nvpe = 0;
  E_Int ngon = 0;

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
    else if (K_STRING::cmp(eltType, "NGON") == 0) { nvpe = 1; cSize = 4+sizeNGon+sizeNFace; ngon = 1; }
    else if (K_STRING::cmp(eltType, "NODE*") == 0) nvpe = 1;
    else if (K_STRING::cmp(eltType, "TRI*") == 0) nvpe = 3;
    else if (K_STRING::cmp(eltType, "QUAD*") == 0) nvpe = 4;
    else if (K_STRING::cmp(eltType, "TETRA*") == 0) nvpe = 4;
    else if (K_STRING::cmp(eltType, "HEXA*") == 0) nvpe = 8;
    else if (K_STRING::cmp(eltType, "BAR*") == 0) nvpe = 2;
    else if (K_STRING::cmp(eltType, "PYRA*") == 0) nvpe = 5;
    else if (K_STRING::cmp(eltType, "PENTA*") == 0) nvpe = 6;
    else if (K_STRING::cmp(eltType, "NGON*") == 0) { nvpe = 1; cSize = 4+sizeNGon+sizeNFace; ngon = 1; }
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
        nvpe = 1; cSize = 4+sizeNGon+sizeNFace; ngon = 1;
        break; 
      default:
        printf("Warning: buildArray: element type is unknown.\n");
        strcpy(eltType, "UNKNOWN");
        break;
    }
  }
  
  // Build array of field
  if (api == 1) // Array1
  { 
    dim[1] = fSize; dim[0] = nfld;
    a = PyArray_SimpleNew(2, dim, NPY_DOUBLE);
  }
  else // Array2
  {
    dim[0] = fSize;
    a = PyList_New(0);
    for (E_Int n=0; n < nfld; n++)
    {
      PyArrayObject* ar = (PyArrayObject*)PyArray_EMPTY(1, dim, NPY_DOUBLE, 1);
      PyList_Append(a, (PyObject*)ar); Py_DECREF(ar);
    }
  } 

  // Build array for connectivity
  if (api == 1) // Array 1
  {
    dim[1] = cSize; dim[0] = nvpe;
    ac = PyArray_SimpleNew(2, dim, NPY_INT);
    E_Int* data = (E_Int*)PyArray_DATA((PyArrayObject*)ac);
    if (ngon == 1)
    {
      data[0] = nface;
      data[1] = sizeNGon;
      data[sizeNGon+2] = nelt;
      data[sizeNGon+3] = sizeNFace;
    }
  }
  else // Array2
  {
    ac = PyList_New(0);
    if (ngon == 0)
    {
      // autre que NGON - connectivite elt-noeuds
      dim[0] = cSize; dim[1] = nvpe;
      PyObject* ar = PyArray_EMPTY(2, dim, NPY_INT, 0);
      PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
    }
    else
    {
      // ngons - NGON - sizeNGon
      dim[0] = sizeNGon;
      PyObject* ar = PyArray_EMPTY(1, dim, NPY_INT, 0);
      PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
      // ngons - NFACE - sizeNFace
      dim[0] = sizeNFace;
      ar = PyArray_EMPTY(1, dim, NPY_INT, 0);
      PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
      // ngons - indPG - nfaces
      dim[0] = nface;
      ar = PyArray_EMPTY(1, dim, NPY_INT, 0);
      PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
      // ngons - indPH - nelts
      dim[0] = nelt;
      ar = PyArray_EMPTY(1, dim, NPY_INT, 0);
      PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
      // Eventuellement PE - 2*nface
      //dim[0] = nface; dim[1] = 2;
      //PyObject* ar = PyArray_EMPTY(2, dim, NPY_INT, 0);
      //PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
    }
  }
  tpl = Py_BuildValue("[sOOs]", varString, a, ac, eltType);
  Py_DECREF(a); Py_DECREF(ac);

  return tpl;
}
