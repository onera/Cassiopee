/*    
    Copyright 2013-2026 ONERA.

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

// Creates a structured points cloud over a BAR.

# include "generator.h"
# include "Nuga/include/Plaster.h"

using namespace std;

//=============================================================================
/* Generation d'une nappe de points approximant une surface s'appuyant sur la 
   BAR en input. */
//=============================================================================
PyObject* K_GENERATOR::fittingPlaster(PyObject* self, PyObject* args)
{
  PyObject *arrB0;
  E_Float pump_factor;
  if (!PYPARSETUPLE_(args, O_ R_, 
                    &arrB0, &pump_factor))
    return NULL;
  E_Int ni, nj, nk;
  K_FLD::FloatArray *fB0;
  char* varString; char* eltType;
  K_FLD::IntArray* cn;

  // Check the contour
  E_Int res = K_ARRAY::getFromArray(arrB0, varString, fB0, ni, nj, nk,
                                    cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "fittingPlaster: invalid array.");
    return NULL;
  }
  
  if (res == 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "fittingPlaster: array must define a BAR.");
    delete fB0; return NULL;
  }
  if (res == 2) 
  {
    if (strcmp(eltType, "BAR") != 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "fittingPlaster: array must define a BAR.");
      delete fB0; delete cn; return NULL;
    }
  }

  // Coordinates
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  
  if ((posx == -1) || (posy == -1) || (posz == -1))
  {
    PyErr_SetString(PyExc_TypeError,
                    "fittingPlaster: can't find coordinates in array.");
    delete fB0; delete cn; return NULL;
  }
  
  K_FLD::FloatArray &posB0(*fB0), plaster; 
  K_FLD::IntArray &connectB0(*cn);
  E_Int ni1;
  Plaster p;
  E_Int err = p.make(posB0, connectB0, plaster, ni1, pump_factor);
  
  delete fB0; delete cn;
  PyObject* tpl = NULL;
  if (!err)
    tpl = K_ARRAY::buildArray(plaster, varString, ni1, plaster.cols()/ni1, 1);
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "fittingPlaster: creation failed.");
  }
  return tpl;
}

//======================  Generator/fittingPlaster.cpp ========================
