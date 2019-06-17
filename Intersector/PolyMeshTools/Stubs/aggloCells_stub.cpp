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


# include "intersector.h"
# include "stub.h"

//=============================================================================
/* Agglomerate cells with a too high aspect ratio */
//=============================================================================
PyObject* K_INTERSECTOR::agglomerateSmallCells(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* Agglomerate cells with a too high aspect ratio */
//=============================================================================
PyObject* K_INTERSECTOR::agglomerateNonStarCells(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* Agglomerate cells where polygons are specified */
//=============================================================================
PyObject* K_INTERSECTOR::agglomerateCellsWithSpecifiedFaces(PyObject* self, PyObject* args)
{
    PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=======================  Intersector/PolyMeshTools/aggloFaces.cpp ====================
