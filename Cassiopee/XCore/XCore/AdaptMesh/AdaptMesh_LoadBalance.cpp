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
#include "Mesh.h"

PyObject *K_XCORE::AdaptMesh_LoadBalance(PyObject *self, PyObject *args)
{
    PyObject *MESH;

    if (!PYPARSETUPLE_(args, O_, &MESH)) {
        RAISE("Wrong input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(MESH, "AdaptMesh")) {
        RAISE("Bad mesh hook.");
        return NULL;
    }

    Mesh *M = (Mesh *)PyCapsule_GetPointer(MESH, "AdaptMesh");

    //if (M->npc == 1) return Py_None;

    // Make global cell-cell connectivity
    Mesh_make_cell_cells(M);

    // Load balance
    E_Int ret = Mesh_load_balance(M);

    if (ret != 0) {
        RAISE("Failed to load-balance the mesh.");
        return NULL;
    }

    if (M->pid == 0) puts("OK LoadBalance.");

    return Py_None;
}
