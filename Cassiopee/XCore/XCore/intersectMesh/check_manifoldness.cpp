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
#include "mesh.h"
#include "common/Karray.h"

PyObject *K_XCORE::check_manifoldness(PyObject *self, PyObject *args)
{
    PyObject *MESH;
    if (!PYPARSETUPLE_(args, O_, &MESH)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray array;
    if (Karray_parse_ngon(MESH, array) != 0) {
        RAISE("Input not NGon.");
        return NULL;
    }

    IMesh M(array);

    bool is_manifold = true;

    M.make_skin();
    Smesh Mf(M, M.skin, false);

    Karray_free_ngon(array);

    return PyBool_FromLong(is_manifold);
}