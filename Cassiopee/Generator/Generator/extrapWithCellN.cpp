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

# include "generator.h"

using namespace std;
    
//=============================================================================
/* 
   extrapWithCellN
   IN: basic elt mesh + cellN in center + field in center
   OUT: field modifie in-place
 */
//=============================================================================
PyObject* K_GENERATOR::extrapWithCellN(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* cellNObject;
  if (!PYPARSETUPLE_(args, OO_, &array, &cellNObject))
  {
    return NULL;
  }

  // Check array: must be HEXA mesh
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res != 2 || strcmp(eltType, "NGON") == 0) 
  {
    if (res != 0) RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "extrapWithCellN: a must be basic element.");
    return NULL;
  }
  
  // Check cellN (like array but in center)
  E_Int ni1, nj1, nk1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray3(cellNObject, varString1, f1, ni1, nj1, nk1, 
                                      cn1, eltType1);
  
  E_Int posCellN = K_ARRAY::isCellNatureField1Present(varString1);
  if (posCellN == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    RELEASESHAREDB(res1, cellNObject, f1, cn1);
    PyErr_SetString(PyExc_TypeError,
                    "extrapWithCellN: cellN must be present in second array.");
    return NULL;
  }
  E_Int nfld = f1->getNfld();
  E_Int nv = f->getSize();
  E_Int ne = cn->getSize();

  E_Int cellNi, cellNv, ind;

  // Construction de la connectivite elt/elt voisins
  vector< vector<E_Int> > cEEN(ne);
  
  //printf("nv=%d, ne=%d, nfld=%d\n", nv, ne, nfld);
  K_CONNECT::connectEV2EENbrs(eltType, nv, *cn, cEEN);
  
  // Extrapolation du champ
  E_Float* cellN = f1->begin(posCellN+1);
  
  for (E_Int i = 0; i < ne; i++)
  {
    cellNi = cellN[i];
    for (size_t n = 0; n < cEEN[i].size(); n++)
    {
      ind = cEEN[i][n];
      cellNv = cellN[ind];
      if (cellNi == 1 && cellNv == 0)
      {
        //printf("extrapolation en %d\n",i);
        for (E_Int n = 1; n <= nfld; n++)
          if (n-1 != posCellN) (*f1)(ind, n) = (*f1)(i, n);
      }
      else if (cellNi == 0 && cellNv == 1)
      {
        //printf("extrapolation en %d\n",i);
        for (E_Int n = 1; n <= nfld; n++)
          if (n-1 != posCellN) (*f1)(i, n) = (*f1)(ind, n);
      }
      else if (cellNi == 0)
      {
        for (E_Int n = 1; n <= nfld; n++)
          if (n-1 != posCellN) (*f1)(i, n) = 0.;
      }
    }
  }
  
  RELEASESHAREDB(res, array, f, cn);
  RELEASESHAREDB(res1, cellNObject, f1, cn1);
  
  Py_INCREF(Py_None);
  return Py_None;

}
