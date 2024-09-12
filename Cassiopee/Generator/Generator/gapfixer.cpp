/*    
    Copyright 2013-2024 Onera.

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
# include "Nuga/include/GapFixer.h"
using namespace std;

//=============================================================================
/* Maillage d'un trou delimite par une bar a l'aide d'une nappe de points 
   approximant la surface contenant le trou. */
//=============================================================================
PyObject* K_GENERATOR::gapfixer(PyObject* self, PyObject* args)
{
  PyObject* arrC, *arrB0, *arrHP(0);
  E_Int refine=0;
  if (!PYPARSETUPLE_(args, OOO_ I_, &arrB0, &arrC, &arrHP, &refine)) return NULL;

  E_Int ni, nj, nk;
  K_FLD::FloatArray *fC, *fB0;
  char* varString; char* eltType;
  K_FLD::IntArray* cn1, *cn2;

  // Check the contour
  E_Int res1 = K_ARRAY::getFromArray(arrB0, varString, fB0, ni, nj, nk,
                                     cn1, eltType);

  if (res1 != 1 && res1 != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "gapfixer: invalid array.");
    return NULL;
  }
  
  if (res1 == 1) 
  {
    if (ni != 1 && nj != 1 && nk != 1)
    {
      delete fB0; delete cn1; // always delete cn1 because DynArray
      PyErr_SetString(PyExc_TypeError,
                      "gapfixer: array must define a plane.");
      return NULL;
    }
  }
  else if (res1 == 2) 
  {
    if (strcmp(eltType, "TRI") != 0 && 
        strcmp(eltType, "QUAD") != 0 &&
        strcmp(eltType, "NODE") != 0 &&
        strcmp(eltType, "BAR") != 0)
    {
      delete fB0; delete cn1;
      PyErr_SetString(PyExc_TypeError,
                      "gapfixer: array must define a plane.");
      return NULL;
    }
  }

  // Coordinates
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  
  if ((posx == -1) || (posy == -1) || (posz == -1))
  {
    delete fB0; delete cn1; // always delete cn1 because DynArray
    PyErr_SetString(PyExc_TypeError,
                    "gapfixer: can't find coordinates in array.");
    return NULL;
  }

  // Check the cloud.
  E_Int res2 = K_ARRAY::getFromArray(arrC, varString, fC, ni, nj, nk,
                                     cn2, eltType);

  if (res2 != 1) 
  {
    delete fB0; delete cn1;
    if (res2 == 2) { delete fC; delete cn2; }
    PyErr_SetString(PyExc_TypeError,
                    "gapfixer: invalid array. Must be structured.");
    return NULL;
  }

  E_Bool not_a_surface = ((ni != 1) && (nj != 1) && (nk != 1));
  not_a_surface &= ((ni*nj == 1) || (ni*nk == 1) || (nj*nk == 1));
  
  if (not_a_surface)
  {
    delete fB0; delete cn1; delete fC; delete cn2;
    PyErr_SetString(PyExc_TypeError,
                    "gapfixer: array must define a surface.");
    return NULL;
  }

  ni = (ni == 1) ? nj : ni;
  
  // Coordinates
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  
  if ((posx == -1) || (posy == -1) || (posz == -1))
  {
    delete fB0; delete cn1; delete fC; delete cn2;
    PyErr_SetString(PyExc_TypeError,
                    "gapfixer: can't find coordinates in array.");
    return NULL;
  }
  
  // Check the hard points
  K_FLD::FloatArray *fHP = 0;
  if (arrHP != Py_None) // If hard points are specified
  {
    K_FLD::IntArray* cndum;
    E_Int nidum, njdum, nkdum;
    E_Int res = K_ARRAY::getFromArray(arrHP, varString, fHP, nidum, njdum, nkdum, cndum, eltType);
    if (res != -1)
    {
      posx = K_ARRAY::isCoordinateXPresent(varString);
      posy = K_ARRAY::isCoordinateYPresent(varString);
      posz = K_ARRAY::isCoordinateZPresent(varString);
    
      if ((posx == -1) || (posy == -1) || (posz == -1))
      {
        delete fB0; delete cn1; delete fC; delete cn2;
        delete fHP; delete cndum;
        PyErr_SetString(PyExc_TypeError,
                        "gapfixer: can't find coordinates in the hard points array.");
        return NULL;
      }
    }
    delete cndum;  // always delete because DynArray
  }
  
  K_FLD::FloatArray &posB0 = *fB0;
  K_FLD::FloatArray &posC = *fC;
  K_FLD::FloatArray posG;
  K_FLD::IntArray   &connectB0 = *cn1, connectG;

  E_Int err = GapFixer::run(posC, ni, posB0, connectB0, posG, connectG, refine, fHP);

  if (err)
  {
    delete fB0; delete cn1; delete fC; delete cn2;
    delete fHP;
    PyErr_SetString(PyExc_TypeError,
                    "gapfixer: failed to proceed.");
    return NULL;
  }

  PyObject* tpl = K_ARRAY::buildArray(posG, varString, connectG, -1, "TRI",
                                      false);
  delete fB0; delete fC; delete fHP;
  //if (res1 == 2) delete cn1;  // only for FldArray
  delete cn1; delete cn2;  // always delete because DynArray (see getFromArrayDyn)
  return tpl;
}

//========================  Generator/gapfixer.cpp ========================
