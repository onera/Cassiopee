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
#include "mesh.h"

PyObject *K_XCORE::IntersectMesh_ExtractMesh(PyObject *self, PyObject *args)
{
    PyObject *MESH;
    E_Int remove_periodic;
  
    if (!PYPARSETUPLE_(args, O_ I_, &MESH, &remove_periodic)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(MESH, "IntersectMesh")) {
        RAISE("Bad mesh hook.");
        return NULL;
    }

    IMesh *M = (IMesh *)PyCapsule_GetPointer(MESH, "IntersectMesh");

    puts("Extracting mesh...");

    PyObject *master_out = M->export_karray(remove_periodic);

    return master_out;
}