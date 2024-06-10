#pragma once

#include "xcore.h"

struct Karray {
    // Reference to the python array
    PyObject *pyobject;

    K_FLD::FldArrayF *f;
    E_Int ni, nj, nk;

    K_FLD::FldArrayI *cn;

    E_Float *X, *Y, *Z;
    E_Int npts;
};

void Karray_free_ngon(Karray &karray);

E_Int Karray_parse_ngon(PyObject *pyobject, Karray &karray);

void Karray_free_structured(Karray &karray);

E_Int Karray_parse_structured(PyObject *pyobject, Karray &karray);