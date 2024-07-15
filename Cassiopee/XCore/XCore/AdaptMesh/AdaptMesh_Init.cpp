#include "Mesh.h"
#include "Karray.h"
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

    Int ret;

    ret = Karray_parse_ngon(ARRAY, karray);

    if (ret != 0) return NULL;

    // Parse cell/face/point global ids
    Int *cgids, *fgids, *pgids;

    Int size, nfld;
    Int nc = karray.ncells();
    Int nf = karray.nfaces();

    if (CGIDS != Py_None) {
        ret = K_NUMPY::getFromNumpyArray(CGIDS, cgids, size, nfld, true);
        if (ret != 1 || size != nc || nfld != 1) {
            RAISE("Bad cell global ids.");
            Karray_free_ngon(karray);
            return NULL;
        }
    } else {
        cgids = IntArray(nc);
        for (Int i = 0; i < nc; i++) cgids[i] = i;
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
        for (Int i = 0; i < nf; i++) fgids[i] = i+1;
    }

    // Init Mesh

    Mesh *M = Mesh_from_Karray(&karray);
 
    // Init global-to-local ids maps
    for (Int i = 0; i < M->nf; i++) M->g2lf[fgids[i]-1] = i;
    for (Int i = 0; i < M->nc; i++) M->g2lc[cgids[i]] = i;

    // Init local-to-global ids maps
    M->l2gf = IntArray(M->nf);
    M->l2gc = IntArray(M->nc);

    for (Int i = 0; i < M->nf; i++) M->l2gf[i] = fgids[i]-1;
    for (Int i = 0; i < M->nc; i++) M->l2gc[i] = cgids[i];

    // Parse boundary patches
    // TODO(Imad): error instead of assert
    assert(PyList_Check(BCS));
    M->nbp = PyList_Size(BCS);
    M->bps = (BPatch *)XCALLOC(M->nbp, sizeof(BPatch));

    for (Int i = 0; i < M->nbp; i++) {
        PyObject *PATCH = PyList_GetItem(BCS, i);

        // TODO(Imad): error instead of assert
        assert(PyList_Size(PATCH) == 4);

        // TODO(Imad): error check
        K_NUMPY::getFromNumpyArray(PyList_GetItem(PATCH, 0), M->bps[i].pf,
            M->bps[i].nf, false);
        
        M->bps[i].gid = PyLong_AsLong(PyList_GetItem(PATCH, 1));

        const char *bcname = PyUnicode_AsUTF8(PyList_GetItem(PATCH, 2));

        M->bps[i].name = (char *)XMALLOC(strlen(bcname)+1);
        strcpy(M->bps[i].name, bcname);
        Int len = strlen(M->bps[i].name);
        assert(M->bps[i].name[len] == '\0');

        const char *bctype = PyUnicode_AsUTF8(PyList_GetItem(PATCH, 3));

        M->bps[i].type = (char *)XMALLOC(strlen(bctype)+1);
        strcpy(M->bps[i].type, bctype);
        len = strlen(M->bps[i].type);
        assert(M->bps[i].type[len] == '\0');

        // Convert to local ids
        BPatch *P = &M->bps[i];

        // Make faces zero-based
        for (Int j = 0; j < P->nf; j++) {
            assert(P->pf[j] > 0 && P->pf[j] <= M->nf);
            P->pf[j]--;
        }
    }

    // Face to bpatch map
    for (Int i = 0; i < M->nbp; i++) {
        BPatch *P = &M->bps[i];

        for (Int j = 0; j < P->nf; j++) {
            Int lfid = P->pf[j];
            M->face_to_bpatch[lfid] = i;
        }
    }

    // Parse comm patches
    // TODO(Imad): error instead of assert
    assert(PyList_Check(COMM));
    M->npp = PyList_Size(COMM);
    M->pps = (PPatch *)XCALLOC(M->npp, sizeof(PPatch));

    for (Int i = 0; i < M->npp; i++) {
        PyObject *PATCH = PyList_GetItem(COMM, i);

        // TODO(Imad): error instead of assert
        assert(PyList_Size(PATCH) == 3);

        // TODO(Imad): error check
        M->pps[i].nei = PyLong_AsLong(PyList_GetItem(PATCH, 0));

        // TODO(Imad): error check
        K_NUMPY::getFromNumpyArray(PyList_GetItem(PATCH, 1),
            M->pps[i].pf, M->pps[i].nf, false);
        
        // TODO(Imad): error check
        K_NUMPY::getFromNumpyArray(PyList_GetItem(PATCH, 2),
            M->pps[i].pn, M->pps[i].nf, false);

        // Zero-based
        for (Int j = 0; j < M->pps[i].nf; j++) {
            M->pps[i].pf[j]--;
        }

        // Send/receive buffers allocated on-the-fly
        M->pps[i].sbuf_i = NULL;
        M->pps[i].rbuf_i = NULL;
    }

    // Face to ppatch map
    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        for (Int j = 0; j < P->nf; j++) {
            Int lfid = P->pf[j];
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
    }

    if (M->mode_2D) {
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

    M->fparent = IntArray(M->nf);
    for (Int i = 0; i < M->nf; i++) M->fparent[i] = i;

    // Clean-up

    Karray_free_ngon(karray);

    // TODO(Imad): Python2
    PyObject *hook = PyCapsule_New((void *)M, "AdaptMesh", NULL);

    return hook;
}
