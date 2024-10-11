#include "common/Karray.h"

PyObject *K_XCORE::extractFacesFromPointTag(PyObject *self, PyObject *args)
{
    PyObject *ARR, *TAG;
    if (!PYPARSETUPLE_(args, OO_, &ARR, &TAG)) {
        RAISE("Bad input.");
        return NULL;
    }

    Karray karray;
    E_Int ret = Karray_parse_ngon(ARR, karray);
    if (ret != 0) return NULL;

    E_Float *tag = NULL;
    E_Int nfld, size;
    ret = K_NUMPY::getFromNumpyArray(TAG, tag, size, nfld, true);
    if (ret != 1 || nfld != 1 || size != karray.npoints()) {
        RAISE("Bad input tag.");
        Karray_free_ngon(karray);
        return NULL;
    }

    // Count the faces to extract
    E_Int nf = 0;

    for (E_Int fid = 0; fid < karray.nfaces(); fid++) {
        E_Int np;
        E_Int *pn = karray.get_face(fid, np);
        E_Int i;
        for (i = 0; i < np; i++) {
            if (tag[pn[i]-1] != 1.0) {
                break;
            }
        }
        if (i == np)
            nf++;
    }

    // Allocate
    npy_intp dims[2];
    dims[1] = 1;
    dims[0] = (npy_intp)nf;
    
    PyArrayObject *out = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
    E_Int *pf = (E_Int *)PyArray_DATA(out);
    E_Int *ptr = pf;
    
    for (E_Int fid = 0; fid < karray.nfaces(); fid++) {
        E_Int np;
        E_Int *pn = karray.get_face(fid, np);
        E_Int i;
        for (i = 0; i < np; i++) {
            if (tag[pn[i]-1] != 1.0) {
                break;
            }
        }
        if (i == np)
            *ptr++ = fid;
    }

    return (PyObject *)out;
}
