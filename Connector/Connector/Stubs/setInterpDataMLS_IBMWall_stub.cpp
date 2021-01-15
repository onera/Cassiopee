/*
    Copyright 2013-2021 Onera.

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

// Computes MLS coefficients for the projection of IBM solution on a triangulated surface

# include "connector.h"
# include "stub.h"


// ============================================================================
/*  IN: nuage de points donneur defini par une zone NODE
    IN : zone surfacique triangulaire receptrice 
    OUT:  donnees d interpolation du nuage sur la surface */
// ============================================================================
PyObject* K_CONNECTOR::setInterpData_IBMWall(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}
