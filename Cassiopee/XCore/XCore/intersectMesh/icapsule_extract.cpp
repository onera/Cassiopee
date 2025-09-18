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
#include "icapsule.h"

PyObject *K_XCORE::icapsule_extract_master(PyObject *self, PyObject *args)
{
    PyObject *ICAPSULE;
    if (!PYPARSETUPLE_(args, O_, &ICAPSULE)) 
    {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(ICAPSULE, "ICAPSULE"))
    {
        RAISE("Bad ICapsule hook.");
        return NULL;
    }

    ICapsule *icap = (ICapsule *)PyCapsule_GetPointer(ICAPSULE, "ICAPSULE");
    const auto &M = icap->M;

    PyObject *out = PyList_New(0);

    PyObject *Mout = M.export_karray();
    PyList_Append(out, Mout);
    Py_DECREF(Mout);

    // Extract cell tags
    npy_intp dims[2];
    assert(M.ctag.size() == (size_t)M.nc);
    dims[0] = (npy_intp)M.nc;
    dims[1] = 1;
    PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    E_Float *ptr = (E_Float *)PyArray_DATA(arr);
    for (E_Int i = 0; i < M.nc; i++) ptr[i] = M.ctag[i];
    PyList_Append(out, (PyObject *)arr);

    return out;
}

PyObject *K_XCORE::icapsule_extract_slaves(PyObject *self, PyObject *args)
{
    PyObject *ICAPSULE;

    if (!PYPARSETUPLE_(args, O_, &ICAPSULE)) 
    {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(ICAPSULE, "ICAPSULE")) 
    {
        RAISE("Bad ICapsule hook.");
        return NULL;
    }

    ICapsule *icap = (ICapsule *)PyCapsule_GetPointer(ICAPSULE, "ICAPSULE");

    PyObject *out = PyList_New(0);

    PyObject *arrays = PyList_New(0);
    PyObject *ctags = PyList_New(0);

    for (size_t i = 0; i < icap->Ss.size(); i++) {
        const auto &S = icap->Ss[i];
        PyObject *sarray = S.export_karray();
        PyList_Append(arrays, sarray);
        Py_DECREF(sarray);

        npy_intp dims[2];
        assert(S.ctag.size() == (size_t)S.nc);
        dims[0] = (npy_intp)S.nc;
        dims[1] = 1;
        PyArrayObject *arr = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
        E_Float *ptr = (E_Float *)PyArray_DATA(arr);
        for (E_Int i = 0; i < S.nc; i++) ptr[i] = S.ctag[i];
        PyList_Append(ctags, (PyObject *)arr);
        Py_DECREF(arr);
    }

    PyList_Append(out, arrays);
    Py_DECREF(arrays);
    PyList_Append(out, ctags);
    Py_DECREF(ctags);

    return out;
}

PyObject *K_XCORE::icapsule_extract_slave(PyObject *self, PyObject *args)
{
    assert(0 && "Unimplemented");
    return Py_None;
    /*
    PyObject *ICAPSULE;
    E_Int INDEX;
    if (!PYPARSETUPLE_(args, O_ I_, &ICAPSULE, &INDEX)) {
        RAISE("Bad input.");
        return NULL;
    }

    if (!PyCapsule_IsValid(ICAPSULE, "ICAPSULE")) {
        RAISE("Bad ICapsule hook.");
        return NULL;
    }

    ICapsule *icap = (ICapsule *)PyCapsule_GetPointer(ICAPSULE, "ICAPSULE");

    if (INDEX >= (E_Int)icap->Ss.size()) {
        RAISE("Bad slave index.");
        return NULL;
    }

    auto Sout = icap->Ss[INDEX].export_karray();

    return Sout;
    */
}
