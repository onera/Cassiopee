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

#include "kcore.h"
#include "Nuga/include/KdTree.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
/* Transforme une liste d'indices de faces (indElt*nface + noface) en
   une connectivite.
   IN: array: array global
   IN: listFaces: liste des indices de faces
   OUT: numpy array de la connectivite correspondant.
 */
//=============================================================================
PyObject* K_KCORE::indiceFace2Connect(PyObject* self, PyObject* args)
{
  PyObject  *array, *listFaces;
  if (!PyArg_ParseTuple(args, "OO", &array, &listFaces))  return NULL;
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType, true); 

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "indiceFace2Connect: unknown type of array.");
    return NULL;
  }

  if (res == 1)
  {
    delete f;
    PyErr_SetString(PyExc_TypeError,
                    "indiceFace2Connect: can not be used on a structured array.");
    return NULL;
  }

  E_Int n = PyList_Size(listFaces);
  E_Int ind, indElt, nof, nf, nt2;

  char eltTypeFaces[10];
  int fBAR[] = { 
    1, 
    2 };
  int fTRI[] = {
    1, 2,
    2, 3,
    3, 1 };
  int fQUAD[] = {
    1, 2,
    2, 3,
    3, 4,
    4, 1 };
  int fTETRA[] = {
    1, 3, 2,
    1, 2, 4,
    2, 3, 4,
    3, 1, 4 };
  int fHEXA[] = {
    1, 4, 3, 2,
    1, 2, 6, 5,
    2, 3, 7, 6,
    3, 4, 8, 7,
    1, 5, 8, 4,
    5, 6, 7, 8 };
  
  int* fPtr = NULL;
  // nf: nbre de face de l'element d'entree
  // eltTypeFaces: type de l'element de sortie (face de l'elt d'entree)
  // nt2: nbre de noeuds de l'element de sortie
  if (K_STRING::cmp(eltType, "BAR") == 0)
  {
    nf = 2; strcpy(eltTypeFaces, "NODE"); nt2 = 1; fPtr = fBAR;
  }
  else if (K_STRING::cmp(eltType, "TRI") == 0)
  {
    nf = 3; strcpy(eltTypeFaces, "BAR"); nt2 = 2; fPtr = fTRI;
  }
  else if (K_STRING::cmp(eltType, "QUAD") == 0)
  {
    nf = 4; strcpy(eltTypeFaces, "BAR"); nt2 = 2; fPtr = fQUAD;
  }
  else if (K_STRING::cmp(eltType, "TETRA") == 0)
  {
    nf = 4; strcpy(eltTypeFaces, "TRI"); nt2 = 3; fPtr = fTETRA;
  }
  else if (K_STRING::cmp(eltType, "HEXA") == 0)
  {
    nf = 6; strcpy(eltTypeFaces, "QUAD"); nt2 = 4; fPtr = fHEXA;
  }
  else
  {
    delete f; delete cn;
    PyErr_SetString(PyExc_TypeError,
                    "indiceFace2Connect: this element type is not taken into account.");
    return NULL;
  }

  PyObject* a = K_NUMPY::buildNumpyArray(nt2, n, 1);
  E_Int* ptr = K_NUMPY::getNumpyPtrI(a);

  // Selectionne les faces subzones
  for (E_Int i = 0; i < n; i++)
  {
    ind = PyLong_AsLong(PyList_GetItem(listFaces, i));
    indElt = ind / nf;
    nof = ind - indElt*nf; 
    for (E_Int v = 0; v < nt2; v++) ptr[i+v*nt2] = (*cn)(indElt, fPtr[nof*nt2+v]);
  }

  RELEASESHAREDU(array, f, cn);
  
  return a;
}
