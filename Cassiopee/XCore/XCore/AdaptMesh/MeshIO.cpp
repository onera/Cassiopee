#include "Mesh.h"
#include "xcore.h"

PyObject *Mesh_export_karray(Mesh *M)
{
    Int nfld = 3;
    const char *var_string = "CoordinateX,CoordinateY,CoordinateZ";
    const char *elt_type = "NGON";
    Int ngon_type = 3;
    bool center = false;
    Int api = 3;

    Int sizeNFace = Mesh_get_sizeNFace(M);
    Int sizeNGon = Mesh_get_sizeNGon(M);

    PyObject *karray = K_ARRAY::buildArray3(nfld, var_string, M->np, M->nc,
        M->nf, elt_type, sizeNGon, sizeNFace, ngon_type, center,
        api);

    K_FLD::FldArrayF *f;
    K_FLD::FldArrayI *cn;
    K_ARRAY::getFromArray3(karray, f, cn);

    E_Float *px = f->begin(1);
    E_Float *py = f->begin(2);
    E_Float *pz = f->begin(3);

    memcpy(px, M->X, M->np * sizeof(Float));
    memcpy(py, M->Y, M->np * sizeof(Float));
    memcpy(pz, M->Z, M->np * sizeof(Float));

    Int *indpg = cn->getIndPG();
    Int *ngon = cn->getNGon();
    Int *indph = cn->getIndPH();
    Int *nface = cn->getNFace();

    // NGON

    memset(indpg, 0, (M->nf + 1) * sizeof(Int));

    for (Int i = 0; i < M->nf; i++) {
        Int *frange = Mesh_get_frange(M, i);
        for (Int j = 0; j < M->fstride[i]; j++) {
            indpg[i+1] += frange[j];
        }
    }

    indpg[0] = 0;
    for (Int i = 0; i < M->nf; i++) indpg[i+1] += indpg[i];

    assert(indpg[M->nf] == sizeNGon);
    
    Int *ptr = ngon;

    for (Int i = 0; i < M->nf; i++) {
        Int *face = Mesh_get_face(M, i);
        Int *frange = Mesh_get_frange(M, i);

        for (Int j = 0; j < M->fstride[i]; j++) {
            Int *pn = face + 2*j;

            for (Int k = 0; k < frange[j]; k++) *ptr++ = pn[k] + 1;
        }
    }


    // NFACE

    memset(indph, 0, (M->nc + 1) * sizeof(Int));

    for (Int i = 0; i < M->nc; i++) {
        Int *crange = Mesh_get_crange(M, i);
        for (Int j = 0; j < M->cstride[i]; j++) {
            indph[i+1] += crange[j];
        }
    }

    indph[0] = 0;
    for (Int i = 0; i < M->nc; i++) indph[i+1] += indph[i];

    assert(indph[M->nc] == sizeNFace);
    
    ptr = nface;

    for (Int i = 0; i < M->nc; i++) {
        Int *cell = Mesh_get_cell(M, i);
        Int *crange = Mesh_get_crange(M, i);

        for (Int j = 0; j < M->cstride[i]; j++) {
            Int *pf = cell + 4*j;

            for (Int k = 0; k < crange[j]; k++) *ptr++ = pf[k] + 1;
        }
    }

    PyObject *out = PyList_New(0);

    PyList_Append(out, karray);
    Py_DECREF(karray);

    // BCs

    PyObject *bcs = PyList_New(0);

    npy_intp dims[2];
    dims[1] = 1;

    for (Int i = 0; i < M->nbp; i++) {
        BPatch *P = &M->bps[i];

        dims[0] = P->nf;

        PyArrayObject *facelist = (PyArrayObject *)
            PyArray_SimpleNew(1, dims, E_NPY_INT);
        Int *pf = (Int *)PyArray_DATA(facelist);
        
        for (Int j = 0; j < P->nf; j++) pf[j] = P->pf[j] + 1;
        
        PyObject *bc = Py_BuildValue("[Olss]", facelist, P->gid, P->name, P->type);
        Py_DECREF(facelist);
        
        PyList_Append(bcs, bc);
        Py_DECREF(bc);
    }

    PyList_Append(out, bcs);
    Py_DECREF(bcs);

    // Comm patches

    PyObject *pps = PyList_New(0);

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        dims[0] = P->nf;

        PyArrayObject *facelist = (PyArrayObject *)
            PyArray_SimpleNew(1, dims, E_NPY_INT);
        Int *pf = (Int *)PyArray_DATA(facelist); 
    
        PyArrayObject *neilist = (PyArrayObject *)
            PyArray_SimpleNew(1, dims, E_NPY_INT);
        Int *pn = (Int *)PyArray_DATA(neilist);

        for (Int j = 0; j < P->nf; j++) pf[j] = P->pf[j] + 1;
        
        for (Int j = 0; j < P->nf; j++) pn[j] = P->pn[j];

        PyObject *patch = Py_BuildValue("[lOO]", P->nei, facelist, neilist);
        Py_DECREF(facelist);
        Py_DECREF(neilist);

        PyList_Append(pps, patch);
        Py_DECREF(patch);
    }

    PyList_Append(out, pps);
    Py_DECREF(pps);

    return out;
}
