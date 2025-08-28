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
# include "converter.h"
# include <algorithm> 

using namespace K_FLD;


//=============================================================================
/* Make the difference between two index numpy arrays
  On suppose que B est un sous ensemble de A
  Retourne les indices de A qui ne sont pas dans B
  NB: Si B contient des indices qui ne sont pas dans A, affiche un warning mais
  calcule quand meme la difference
  Attention: l'ordre des indices dans A et B est modifie
*/
//=============================================================================
PyObject* K_CONVERTER::diffIndex(PyObject* self, PyObject* args)
{
  PyObject* arrayA; PyObject* arrayB;
  if (!PYPARSETUPLE_(args, OO_, &arrayA, &arrayB)) return NULL;

  // Check numpy (indexA)
  FldArrayI* indexA;
  E_Int res = K_NUMPY::getFromNumpyArray(arrayA, indexA);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "diffIndex: arrayA numpy is invalid.");
    return NULL;
  }

  // Check numpy (indexB)
  FldArrayI* indexB;
  res = K_NUMPY::getFromNumpyArray(arrayB, indexB);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "diffIndex: arrayB numpy is invalid.");
    RELEASESHAREDN(arrayA, indexA);
    return NULL;
  }

  // Construit la difference
  E_Int sizeA = indexA->getSize();
  E_Int sizeB = indexB->getSize();
  E_Int nd = sizeA-sizeB;
  E_Int* posA = indexA->begin();
  E_Int* posB = indexB->begin();
  
  if (sizeA == sizeB) // retourne vide
  {
    PyObject* tpl = K_NUMPY::buildNumpyArray(E_Int(0), 1, 1, 1);
    RELEASESHAREDN(arrayA, indexA);
    RELEASESHAREDN(arrayB, indexB);
    return tpl;
  }

  if (sizeB == 0) // retourne une copie de sizeA
  {
    PyObject* tpl = K_NUMPY::buildNumpyArray(sizeA, 1, 1, 1);
    E_Int* pos = K_NUMPY::getNumpyPtrI(tpl);
    #pragma omp parallel for if (sizeA > __MIN_SIZE_MEAN__)
    for (E_Int i = 0; i < sizeA; i++) pos[i] = posA[i];
    RELEASESHAREDN(arrayA, indexA);
    RELEASESHAREDN(arrayB, indexB);
    return tpl;
  }
  
  if (sizeB > sizeA)
  {
    printf("Warning: diffIndex: listB is not included in listA.");
    PyObject* tpl = K_NUMPY::buildNumpyArray(E_Int(0), 1, 1, 1);
    RELEASESHAREDN(arrayA, indexA);
    RELEASESHAREDN(arrayB, indexB);
    return tpl;
  }
  
  // Il faut deja trier A et B
  std::sort(posA, posA+sizeA);
  std::sort(posB, posB+sizeB);
  
  E_Bool isSubset = std::includes(posA, posA+sizeA, posB, posB+sizeB);
  if (not isSubset)
  {
    printf("Warning: diffIndex: listB is not a subset of listA - proceed anyway");
    // On retourne les indices de A qui ne sont pas dans B meme si B
    // contient des indices qui ne sont pas dans A
    /*std::cout << "sizes: " << sizeA << ", " << sizeB << std::endl;
    for (E_Int i = 0; i < sizeA; i++) std::cout << " " << posA[i];
    std::cout << "\ndone A\n" << std::endl;
    for (E_Int i = 0; i < sizeB; i++) std::cout << " " << posB[i];
    std::cout << "\ndone B\n" << std::endl;*/
    
    std::vector<E_Int> subset(sizeA);  // alloue initialement a la taille max, ie, sizeA
    auto iter = std::set_difference(posA, posA+sizeA, posB, posB+sizeB, subset.begin());
    nd = iter - subset.begin();  // corrige la taille du subset
    //std::cout << "nd: " << nd << std::endl;
    
    PyObject* tpl = K_NUMPY::buildNumpyArray(nd, 1, 1, 1);
    E_Int* pos = K_NUMPY::getNumpyPtrI(tpl);
    #pragma omp parallel for if (nd > __MIN_SIZE_MEAN__)
    for (E_Int i = 0; i < nd; i++) pos[i] = subset[i];
    
    /*for (E_Int i = 0; i < nd; i++) std::cout << " " << pos[i];
    std::cout << "\ndone subset\n" << std::endl;*/

    RELEASESHAREDN(arrayA, indexA);
    RELEASESHAREDN(arrayB, indexB);
    return tpl;
  } 

  PyObject* tpl = K_NUMPY::buildNumpyArray(nd, 1, 1, 1);
  E_Int* pos = K_NUMPY::getNumpyPtrI(tpl);
  
  std::set_difference(posA, posA+sizeA, posB, posB+sizeB, pos);
  
  RELEASESHAREDN(arrayA, indexA);
  RELEASESHAREDN(arrayB, indexB);
  return tpl;
}
