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
# include "converter.h"
using namespace K_FLD;
#include <algorithm> 

//=============================================================================
/* Make the difference between two index numpy arrays
  On suppose que B est un sous ensemble de A
  Retourne les indices de A qui ne sont pas dans B
  Attention: l'ordre des indices dans A et B est modifie
*/
//=============================================================================
PyObject* K_CONVERTER::diffIndex(PyObject* self, PyObject* args)
{
  PyObject* arrayA; PyObject* arrayB;
  if (!PYPARSETUPLEI(args, "OO", "OO", &arrayA, &arrayB)) return NULL;

  // Check numpy (indexA)
  FldArrayI* indexA;
  E_Int res = K_NUMPY::getFromNumpyArray(arrayA, indexA, true);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "diffIndex: arrayA numpy is invalid.");
    return NULL;
  }

  // Check numpy (indexB)
  FldArrayI* indexB;
  res = K_NUMPY::getFromNumpyArray(arrayB, indexB, true);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "diffIndex: arrayB numpy is invalid.");
    return NULL;
  }

  // Construit la difference
  E_Int sizeA = indexA->getSize();
  E_Int sizeB = indexB->getSize();
  E_Int nd = sizeA-sizeB;
  E_Int* posA = indexA->begin();
  E_Int* posB = indexB->begin();
  //printf("diffindex: %d %d %d\n", sizeA, sizeB, nd);
  
  if (sizeA == sizeB) // retourne vide?
  {
    PyObject* tpl = K_NUMPY::buildNumpyArray(0, 1, 1, 1);
    return tpl;
  }

  if (sizeB > sizeA)
  {
    printf("Warning: diffIndex: listB is not included in listA.");
    PyObject* tpl = K_NUMPY::buildNumpyArray(0, 1, 1, 1);
    return tpl;
  }

  if (sizeB == 0) // retourne une copie de sizeA
  {
    PyObject* tpl = K_NUMPY::buildNumpyArray(sizeA, 1, 1, 1);
    E_Int* pos = K_NUMPY::getNumpyPtrI(tpl);
    for (E_Int i = 0; i < sizeA; i++) pos[i] = posA[i];
    return tpl;
  }    

  // Il faut deja trier A et B
  std::sort(posA, posA+sizeA);
  std::sort(posB, posB+sizeB);

  PyObject* tpl = K_NUMPY::buildNumpyArray(nd, 1, 1, 1);
  E_Int* pos = K_NUMPY::getNumpyPtrI(tpl);
  
  /*
  bool isIncluded = std::includes(indexA->begin(), indexA->begin()+sizeA,
                      indexB->begin(), indexB->begin()+sizeB);
  if (not isIncluded) printf("not included\n"); */
  std::set_difference(posA, posA+sizeA,
                      posB, posB+sizeB,
                      pos);

  /*
  for (E_Int i = 0; i < nd; i++) printf("%d\n", pos[i]);
  printf("done\n");
  */
  RELEASESHAREDN(arrayA, indexA);
  RELEASESHAREDN(arrayB, indexB);
  
  return tpl;
}
