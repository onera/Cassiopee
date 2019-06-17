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

E_Int K_INTERSECTOR::check_is_NGON(PyObject* arr, K_FLD::FloatArray*& f1, K_FLD::IntArray*& cn1, char*& varString, char*& eltType)
{
  return 0;
}

PyObject* K_INTERSECTOR::updatePointLists(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL; 
}

//=============================================================================
/* Triangulates exterior faces (any Polygon). */
//=============================================================================
PyObject* K_INTERSECTOR::triangulateExteriorFaces(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* Triangulates BC. */
//=============================================================================
PyObject* K_INTERSECTOR::triangulateSpecifiedFaces(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=============================================================================
/* Convexify any concave polygon found in the input polyhedral mesh. */
//=============================================================================
PyObject* K_INTERSECTOR::convexifyFaces(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//============================================================================================================
/* Split (convexify, starify) some targeted polygons on targeted cells 
 * (typically bad polyhedra -concaves, non-centroid-star-shaped-)
 * to prepare the split of those bad cells.*/
//============================================================================================================
PyObject* K_INTERSECTOR::prepareCellsSplit(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return NULL;
}

//=======================  Intersector/PolyMeshTools/split_faces.cpp ====================
