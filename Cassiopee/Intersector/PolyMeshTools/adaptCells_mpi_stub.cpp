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

# include "intersector.h"
# include "stub.h"

//=============================================================================
/* Hierarchical Mesh Adaptation : MPI version (has MPI calls) */
//=============================================================================
PyObject* K_INTERSECTOR::adaptCells_mpi(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return nullptr;
}

//=============================================================================
/* Exchange the owned PointLists : to update PointListDonor (has MPI calls) */
//=============================================================================
PyObject* K_INTERSECTOR::exchangePointLists(PyObject* self, PyObject* args)
{
  PyErr_SetString(PyExc_NotImplementedError, STUBMSG);
  return nullptr;
}
