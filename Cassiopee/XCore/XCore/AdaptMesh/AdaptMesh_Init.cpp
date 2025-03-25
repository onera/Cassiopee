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
#include "Mesh.h"
#include "common/Karray.h"
#include "common/mem.h"

PyObject *K_XCORE::AdaptMesh_Init(PyObject *self, PyObject *args)
{
    PyObject *ARRAY, *MODE2D, *BCS;
    PyObject *COMM, *CGIDS, *FGIDS;

    if (!PYPARSETUPLE_(args, OOO_ OOO_, &ARRAY, &MODE2D, &BCS, &COMM,
                       &CGIDS, &FGIDS)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray karray;

    E_Int ret;

    ret = Karray_parse_ngon(ARRAY, karray);

    if (ret != 0) return NULL;

    // Parse cell/face/point global ids
    E_Int *cgids, *fgids;

    E_Int size, nfld;
    E_Int nc = karray.ncells();
    E_Int nf = karray.nfaces();

    if (CGIDS != Py_None) {
        ret = K_NUMPY::getFromNumpyArray(CGIDS, cgids, size, nfld, true);
        if (ret != 1 || size != nc || nfld != 1) {
            RAISE("Bad cell global ids.");
            Karray_free_ngon(karray);
            return NULL;
        }
    } else {
        cgids = IntArray(nc);
        for (E_Int i = 0; i < nc; i++) cgids[i] = i;
    }

    if (FGIDS != Py_None) {
        ret = K_NUMPY::getFromNumpyArray(FGIDS, fgids, size, nfld, true);
        if (ret != 1 || size != nf || nfld != 1) {
            RAISE("Bad face global ids.");
            Karray_free_ngon(karray);
            return NULL;
        }
    } else {
        fgids = IntArray(nf);
        for (E_Int i = 0; i < nf; i++) fgids[i] = i+1;
    }

    // Init Mesh

    Mesh *M = Mesh_from_Karray(&karray);
    //printf("Initial cells: %d\n", M->nc);
    //printf("Initial faces: %d\n", M->nf);
    //printf("Initial points: %d\n", M->np);
 
    // Init global-to-local ids maps
    for (E_Int i = 0; i < M->nf; i++) M->g2lf[fgids[i]-1] = i;
    for (E_Int i = 0; i < M->nc; i++) M->g2lc[cgids[i]] = i;

    // Init local-to-global ids maps
    M->l2gf = IntArray(M->nf);
    M->l2gc = IntArray(M->nc);

    for (E_Int i = 0; i < M->nf; i++) M->l2gf[i] = fgids[i]-1;
    for (E_Int i = 0; i < M->nc; i++) M->l2gc[i] = cgids[i];

    // Parse boundary patches
    // TODO(Imad): error instead of assert
    assert(PyList_Check(BCS));
    M->nbp = PyList_Size(BCS);
    M->bps = (BPatch *)XCALLOC(M->nbp, sizeof(BPatch));

    for (E_Int i = 0; i < M->nbp; i++) {
        PyObject *PATCH = PyList_GetItem(BCS, i);

        // TODO(Imad): error instead of assert
        assert(PyList_Size(PATCH) == 4);

        // TODO(Imad): error check
        E_Int *ptr = NULL;
        K_NUMPY::getFromNumpyArray(PyList_GetItem(PATCH, 0), ptr,
            M->bps[i].nf, true);
        
        M->bps[i].pf = IntArray(M->bps[i].nf);
        for (E_Int j = 0; j < M->bps[i].nf; j++) M->bps[i].pf[j] = ptr[j];
        
        M->bps[i].gid = PyLong_AsLong(PyList_GetItem(PATCH, 1));

#if PY_VERSION_HEX >= 0x03000000
        const char *bcname = PyUnicode_AsUTF8(PyList_GetItem(PATCH, 2));
        const char *bctype = PyUnicode_AsUTF8(PyList_GetItem(PATCH, 3));
#else
        const char *bcname = PyString_AsString(PyList_GetItem(PATCH, 2));
        const char *bctype = PyString_AsString(PyList_GetItem(PATCH, 3));
#endif

        M->bps[i].name = (char *)XMALLOC(strlen(bcname)+1);
        strcpy(M->bps[i].name, bcname);

        M->bps[i].type = (char *)XMALLOC(strlen(bctype)+1);
        strcpy(M->bps[i].type, bctype);

        // Convert to local ids
        BPatch *P = &M->bps[i];

        // Make faces zero-based
        for (E_Int j = 0; j < P->nf; j++) {
            assert(P->pf[j] > 0 && P->pf[j] <= M->nf);
            P->pf[j]--;
        }
    }

    // Face to bpatch map
    for (E_Int i = 0; i < M->nbp; i++) {
        BPatch *P = &M->bps[i];

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int lfid = P->pf[j];
            M->face_to_bpatch[lfid] = i;
        }
    }

    // Parse comm patches
    // TODO(Imad): error instead of assert
    assert(PyList_Check(COMM));
    M->npp = PyList_Size(COMM);
    M->pps = (PPatch *)XCALLOC(M->npp, sizeof(PPatch));

    for (E_Int i = 0; i < M->npp; i++) {
        PyObject *PATCH = PyList_GetItem(COMM, i);

        // TODO(Imad): error instead of assert
        assert(PyList_Size(PATCH) == 3);

        // TODO(Imad): error check
        M->pps[i].nei = PyLong_AsLong(PyList_GetItem(PATCH, 0));

        // TODO(Imad): error check
        E_Int *ptr = NULL;
        K_NUMPY::getFromNumpyArray(PyList_GetItem(PATCH, 1),
            ptr, M->pps[i].nf, true);
        
        M->pps[i].pf = IntArray(M->pps[i].nf);
        memcpy(M->pps[i].pf, ptr, M->pps[i].nf * sizeof(E_Int));
        
        // TODO(Imad): error check
        ptr = NULL;
        K_NUMPY::getFromNumpyArray(PyList_GetItem(PATCH, 2),
            ptr, M->pps[i].nf, true);

        M->pps[i].pn = IntArray(M->pps[i].nf);
        memcpy(M->pps[i].pn, ptr, M->pps[i].nf * sizeof(E_Int));

        // Zero-based
        for (E_Int j = 0; j < M->pps[i].nf; j++) {
            M->pps[i].pf[j]--;
        }

        // Send/receive buffers allocated on-the-fly
        M->pps[i].sbuf_i = NULL;
        M->pps[i].rbuf_i = NULL;
    }

    // Face to ppatch map
    for (E_Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        for (E_Int j = 0; j < P->nf; j++) {
            E_Int lfid = P->pf[j];
            M->face_to_ppatch[lfid] = i;
        }
    }

    // Set face types/levels
    Mesh_set_face_types(M);
    M->flevel = IntArray(M->nf);

    // Set cell types
    Mesh_set_cell_types(M);
    M->clevel = IntArray(M->nc);

    // Prepare cells for 2D
    if (MODE2D != Py_None) {
        K_NUMPY::getFromNumpyArray(MODE2D, M->mode_2D, size, nfld, false);
        if (size != 3 || nfld != 1) {
            RAISE("Bad 2D direction array.");
            Karray_free_ngon(karray);
            Mesh_free(M);
            return NULL;
        }

        ret = Mesh_set_cells_for_2D(M);
        if (ret != 0) {
            RAISE("Failed to set cells for 2D.");
            Karray_free_ngon(karray);
            Mesh_free(M);
            return NULL;
        }
    }

    // Set orientation
    ret = Mesh_set_orientation(M);

    if (ret != 0) {
        RAISE("Failed to orient the mesh.");
        Karray_free_ngon(karray);
        Mesh_free(M);
        return NULL;
    }

    // Reorder cells, useful?


    M->fparent = IntArray(M->nf);
    for (E_Int i = 0; i < M->nf; i++) M->fparent[i] = i;

    M->cref = NULL;
    M->fref = NULL;
    M->fpattern = (E_Int *)XMALLOC(M->nf * sizeof(E_Int));
    memset(M->fpattern, -1, M->nf * sizeof(E_Int));

    M->ctag = IntArray(M->nc);
    M->ftag = IntArray(M->nf);
    M->ptag = IntArray(M->np);

    // Clean-up

    Karray_free_ngon(karray);

    // TODO(Imad): Python2
    PyObject *hook = PyCapsule_New((void *)M, "AdaptMesh", NULL);

    return hook;
}
