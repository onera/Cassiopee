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
/* Build a structured array 
   IN: field: Dyn structured field
   IN: varString: variables string
   IN: ni,nj,nk: number of points in field
   OUT: PyObject created. */
//=============================================================================
PyObject* K_ARRAY::buildArray(DynArray<E_Float>& field, const char* varString, 
                              E_Int ni, E_Int nj, E_Int nk)
{
  npy_intp dim[2];
  PyArrayObject* a;
  PyObject* tpl;
  IMPORTNUMPY;

  dim[1] = field.cols();
  dim[0] = field.rows();
  a = (PyArrayObject*)PyArray_SimpleNew(2, dim, NPY_DOUBLE);
  E_Float* d = (E_Float*)PyArray_DATA(a);
  DynArray<E_Float>::iterator it = field.begin();
  E_Int dim0 = dim[0]; E_Int dim1 = dim[1];
  for (E_Int i = 0; i < dim1; i++)
    for (E_Int n = 0; n < dim0; n++)
    {
      d[i + n*dim1] = *it; it++;
    }
  tpl = Py_BuildValue("[sOlll]", varString, a,  
                      (long)ni, (long)nj, (long)nk);
  Py_DECREF(a);
  return tpl;
}

//=============================================================================
/* Build a unstructured array 
   IN: field: Dyn unstructured field
   IN: varString : variable string
   IN: c: Dyn connectivity
   IN: et: elements type: 0 (NODE), 1 (BAR), 2 (TRI), 3 (QUAD),
                          4 (TETRA), 5 (PYRA), 6 (PENTA), 7 (HEXA),
                          8 (NGON).
   IN: etString: if et=-1, etString is used instead.
   IN: center: set to true if field is localised in the centers of
   elements, otherwise let it to false.
   
   OUT: PyObject created. */
//=============================================================================
PyObject* K_ARRAY::buildArray(DynArray<E_Float>& field, const char* varString,
                              DynArray<E_Int>& c, E_Int et, const char* etString,
			                  E_Boolean center)
{
  npy_intp dim[2];
  PyArrayObject* a;
  PyArrayObject* ac;
  PyObject* tpl;
  char eltType[8];
  E_Int shift=1;// temporary hack to preserve exisitng code : need a rethink of indexation start (0 and/or 1)
  E_Int fSize = field.getSize();
  E_Int cSize = c.getSize(); // taille de la connectivite
  E_Int cEltSize = cSize; // nb d'elements dans la connectivite
  if (et == 8 || (et == -1 && K_STRING::cmp(etString, "NGON") == 0))
  {
    E_Int sizeFN = c[1]; cEltSize = c[sizeFN+2];
    shift=0;
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
  dim[1] = field.cols();
  dim[0] = field.rows();
  E_Int dim1 = dim[1];
  a = (PyArrayObject*)PyArray_SimpleNew(2, dim, NPY_DOUBLE);
  E_Float* d = (E_Float*)PyArray_DATA(a);
  DynArray<E_Float>::iterator it = field.begin();
  for (E_Int i = 0; i < dim1; i++)
    for (E_Int n = 0; n < dim[0]; n++)
    {
      d[i + n*dim1] = *it; it++;
    }

  // Build array for connectivity   
  dim[1] = c.cols();
  dim[0] = c.rows();
  dim1 = dim[1];
  ac = (PyArrayObject*)PyArray_SimpleNew(2, dim, NPY_INT);
  E_Int* di = (E_Int*)PyArray_DATA(ac);
  DynArray<E_Int>::iterator iti = c.begin();
  for (E_Int i = 0; i < dim1; i++)
    for (E_Int n = 0; n < dim[0]; n++)
    {
      di[i + n*dim1] = *iti+shift; iti++;
    }

  tpl = Py_BuildValue("[sOOs]", varString, a, ac, eltType);
  Py_DECREF(a); Py_DECREF(ac);

  return tpl;
}
