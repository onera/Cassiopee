/*    
    Copyright 2013-2019 Onera.

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
/* Build a structured array from a FldArrayF
   IN: varString: variables string
   IN: ni,nj,nk: number of points in field
   IN: api (1: array, 2: array2)
   OUT: PyObject created. */
//=============================================================================
/*
PyObject* K_ARRAY::buildArray2(FldArrayF& field, const char* varString, 
                               E_Int ni, E_Int nj, E_Int nk, E_Int api)
{
  PyObject* o = buildArray2(field.getNfld(), varString, 
                            ni, nj, nk, api);
  FldArrayF* fp; char* varString2;
  getFromArray2(o, varString2, fp);
  delete fp;
  return o;
}
*/
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
  char eltType[12];

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
    E_Int l = strlen(eltType);
    if (eltType[l-1] == '*') l = l-1;

    if (K_STRING::cmp(eltType, l, "NODE") == 0) nvpe = 1;
    else if (K_STRING::cmp(eltType, l, "TRI") == 0) nvpe = 3;
    else if (K_STRING::cmp(eltType, l, "QUAD") == 0) nvpe = 4;
    else if (K_STRING::cmp(eltType, l, "TETRA") == 0) nvpe = 4;
    else if (K_STRING::cmp(eltType, l, "HEXA") == 0) nvpe = 8;
    else if (K_STRING::cmp(eltType, l, "BAR") == 0) nvpe = 2;
    else if (K_STRING::cmp(eltType, l, "PYRA") == 0) nvpe = 5;
    else if (K_STRING::cmp(eltType, l, "PENTA") == 0) nvpe = 6;
    else if (K_STRING::cmp(eltType, l, "NGON") == 0) { nvpe = 1; cSize = 4+sizeNGon+sizeNFace; ngon = 1; }
    // quadratic
    else if (K_STRING::cmp(eltType, l, "TRI_6") == 0) nvpe = 6;
    else if (K_STRING::cmp(eltType, l, "QUAD_8") == 0) nvpe = 8;
    else if (K_STRING::cmp(eltType, l, "TETRA_10") == 0) nvpe = 10;
    else if (K_STRING::cmp(eltType, l, "BAR_3") == 0) nvpe = 3;
    else if (K_STRING::cmp(eltType, l, "HEXA_20") == 0) nvpe = 20;
    else if (K_STRING::cmp(eltType, l, "QUAD_9") == 0) nvpe = 9;
    else if (K_STRING::cmp(eltType, l, "PYRA_14") == 0) nvpe = 14;
    else if (K_STRING::cmp(eltType, l, "PENTA_15") == 0) nvpe = 15;
    else if (K_STRING::cmp(eltType, l, "PENTA_18") == 0) nvpe = 18;
    else if (K_STRING::cmp(eltType, l, "HEXA_27") == 0) nvpe = 27;
    else if (K_STRING::cmp(eltType, l, "PYRA_13") == 0) nvpe = 13;
    // cubic
    else if (K_STRING::cmp(eltType, l, "TRI_9") == 0) nvpe = 9;
    else if (K_STRING::cmp(eltType, l, "TRI_10") == 0) nvpe = 10;
    else if (K_STRING::cmp(eltType, l, "QUAD_12") == 0) nvpe = 12;
    else if (K_STRING::cmp(eltType, l, "QUAD_16") == 0) nvpe = 16;
    else if (K_STRING::cmp(eltType, l, "BAR_4") == 0) nvpe = 4;
    else if (K_STRING::cmp(eltType, l, "TETRA_16") == 0) nvpe = 16;
    else if (K_STRING::cmp(eltType, l, "TETRA_20") == 0) nvpe = 20;
    else if (K_STRING::cmp(eltType, l, "PYRA_21") == 0) nvpe = 21;
    else if (K_STRING::cmp(eltType, l, "PYRA_29") == 0) nvpe = 29;
    else if (K_STRING::cmp(eltType, l, "HEXA_32") == 0) nvpe = 32;
    else if (K_STRING::cmp(eltType, l, "HEXA_56") == 0) nvpe = 56;
    else if (K_STRING::cmp(eltType, l, "HEXA_64") == 0) nvpe = 64;
    }
  else
  {
    switch (et)
    {
      // Elements lineaires
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

      // Elements quadratiques
      case 10:
        strcpy(eltType, "BAR_3");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 3;
        break;
      case 11:
        strcpy(eltType, "TRI_6");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 6;
        break;
      case 12:
        strcpy(eltType, "QUAD_8");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 8;
        break;
      case 13:
        strcpy(eltType, "QUAD_9");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 9;
        break;
      case 14:
        strcpy(eltType, "TETRA_10");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 10;
        break;
      case 15:
        strcpy(eltType, "PYRA_14");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 14;
        break;
      case 16:
        strcpy(eltType, "PENTA_15");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 15;
        break;
      case 17:
        strcpy(eltType, "PENTA_18");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 18;
        break;
      case 18:
        strcpy(eltType, "HEXA_20");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 20;
        break;
      case 19:
        strcpy(eltType, "HEXA_27");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 27;
        break;
      case 20:
        strcpy(eltType, "PYRA_13");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 13;
        break;

      // Elements cubiques
      case 30:
        strcpy(eltType, "BAR_4");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 4;
        break;
      case 31:
        strcpy(eltType, "TRI_9");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 9;
        break;
      case 32:
        strcpy(eltType, "TRI_10");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 10;
        break;
      case 33:
        strcpy(eltType, "QUAD_12");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 12;
        break;
      case 34:
        strcpy(eltType, "QUAD_16");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 16;
        break;
      case 35:
        strcpy(eltType, "TETRA_16");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 16;
        break;
      case 36:
        strcpy(eltType, "TETRA_20");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 20;
        break;
      case 37:
        strcpy(eltType, "PYRA_21");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 21;
        break;
      case 38:
        strcpy(eltType, "PYRA_29");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 29;
        break;
      case 39:
        strcpy(eltType, "PYRA_30");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 30;
        break;
      case 40:
        strcpy(eltType, "HEXA_32");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 32;
        break;
      case 41:
        strcpy(eltType, "HEXA_56");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 56;
        break;
      case 42:
        strcpy(eltType, "HEXA_64");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 64;
        break;
                    
      // Elements quartiques
      case 51:
        strcpy(eltType, "BAR_5");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 5;
        break;
      case 52:
        strcpy(eltType, "TRI_12");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 12;
        break;
      case 53:
        strcpy(eltType, "TRI_15");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 15;
        break;
      case 54:
        strcpy(eltType, "QUAD_25");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 25;
        break;
      case 55:
        strcpy(eltType, "TETRA_35");
        if (cSize == fSize && center == true) strcat(eltType, "*");
        nvpe = 35;
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
