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
 
#include "compressor.h"

using namespace K_FLD;

//=============================================================================
// Compression delta de deux numpy d'index
// index1: a comparer
// index2: reference
//=============================================================================
PyObject* K_COMPRESSOR::deltaIndex(PyObject* self, PyObject* args)
{
  PyObject* array1; PyObject* array2;
  PyObject* del1; PyObject* del2;
  if (!PYPARSETUPLE_(args, OOO_ O_, &array1, &array2,
                     &del1, &del2)) return NULL;

  // Array d'indexes a comparer
  FldArrayI index;
  E_Int ret = K_ARRAY::getFromList(array1, index);
  if (ret == 0) 
  { 
    PyErr_SetString(
      PyExc_ValueError,
        "deltaIndex: first arg must be an integer list or a numpy.");
    return NULL;
  }
  
  // Array d'indexes de reference
  FldArrayI ref;
  ret = K_ARRAY::getFromList(array2, ref);
  if (ret == 0) 
  {
    PyErr_SetString(
      PyExc_ValueError,
        "deltaIndex: second arg must be an integer list or a numpy.");
    return NULL;
  }

  // in1d de index dans ref
  FldArrayI r1;
  ret = K_ARRAY::getFromList(del1, r1);
  if (ret == 0)
  { 
    PyErr_SetString(
      PyExc_ValueError,
        "deltaIndex: third arg must be an integer list or a numpy.");
    return NULL;
  }
  
  // ind1d de ref dans index
  FldArrayI r2;
  ret = K_ARRAY::getFromList(del2, r2);
  if (ret == 0)
  { 
    PyErr_SetString(
      PyExc_ValueError,
        "deltaIndex: fourth arg must be an integer list or a numpy.");
    return NULL;
  }

  // build deltas
  E_Int* r1p = r1.begin();
  E_Int* r2p = r2.begin();

  E_Int s1 = r1.getSize(); // taille de index
  E_Int s2 = r2.getSize(); // taille de ref
  E_Int Nadd = 0;
  for (E_Int i = 0; i < s1; i++)
  {
    if (r1p[i] == false) Nadd++;
  }
  E_Int Nsupp = 0;
  for (E_Int i = 0; i < s2; i++)
  {
    if (r2p[i] == false) Nsupp++;
  }
  // build le numpy de sortie
  PyObject* a = K_NUMPY::buildNumpyArray(Nadd + Nsupp + 2, 1, 1);
  E_Int* ptra = K_NUMPY::getNumpyPtrI(a);

  ptra[0] = Nadd;
  E_Int c = 1;
  for (E_Int i = 0; i < s1; i++)
  {
    if (r1p[i] == false) { ptra[c] = index[i]; c++; }
  }
  ptra[c] = Nsupp; c++;
  for (E_Int i = 0; i < s2; i++)
  {
    if (r2p[i] == false) { ptra[c] = ref[i]; c++; }
  }
  return a;
}
