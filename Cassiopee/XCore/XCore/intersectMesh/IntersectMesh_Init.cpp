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

PyObject *K_XCORE::IntersectMesh_Init(PyObject *self, PyObject *args)
{
    PyObject *ARRAY, *TAGS;

    if (!PYPARSETUPLE_(args, OO_, &ARRAY, &TAGS)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray karray;

    E_Int ret;

    ret = Karray_parse_ngon(ARRAY, karray);

    if (ret != 0) return NULL;

    // Init mesh

    IMesh *M = new IMesh(karray);

    if (TAGS != Py_None) {
        E_Int size, nfld;

        E_Float *tags = NULL;

        ret = K_NUMPY::getFromNumpyArray(TAGS, tags, size, nfld);

        if (ret != 1 || size != M->nc || nfld != 1) {
            RAISE("Bad cell tags.")
            delete M;
            Karray_free_ngon(karray);
            return NULL;
        }
        
        M->ctag.resize(M->nc);

        for (E_Int i = 0; i < M->nc; i++) M->ctag[i] = (int)tags[i];
    }

    // Init surface mesh data
    M->make_skin();
    M->Mf = Smesh::Smesh_from_mesh_skin(*M, M->skin, false);
    printf("Mf: %d tris\n", M->Mf.nf);
    M->Mf.make_fcenters();
    M->Mf.make_BVH();
    printf("Mesh bounds: ");
    M->Mf.box.print();

    // Clean-up

    Karray_free_ngon(karray);

    // TODO(Imad): Python2
    PyObject *hook = PyCapsule_New((void *)M, "IntersectMesh", NULL);

    return hook;
}
