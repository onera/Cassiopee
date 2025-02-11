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
#include <array>
#include <iostream>
#define K_ARRAY_UNIQUE_SYMBOL
#include "Python.h"
#include "kcore.h"
#include "zfp.h"

namespace zfp {
const char *compress_doc = R"RAW(
Compress double precision float array (or a list of arrays) with zfp
--------------------------------------------------------------------

Compress a muti-dimensionnal array (number of dimensions must be between one and four ).

Usage :
    zfp.pack(array[,reversible = True or False, rate = <float>, accuracy = <float>])
  or
    zfp.pack(list of arrays, [,reversible = True or False, rate = <float>, accuracy = <float>])

- reversible : If true, the compressor try to compress with a near loss-less compression
- Rate gives a number of bits wanted per 4**d elements (where d is the number of dimensions of the array)
on average. The rate at the return of the function can be slightly different than the provided rate. Useful if you
have to bound the compressed size or if you need random access to blocks.
- Accuracy : Fix the maximum absolute error on each value of the array.
    )RAW";
PyObject *
py_compress(PyObject *self, PyObject *args, PyObject *kwd)
{
    //int         status = 0; /* return value: 0 = success */
    zfp_type    type;       /* array scalar type */
    zfp_field * field=NULL; /* array meta data */
    zfp_stream *zfp;        /* compressed stream */
    void *      buffer;     /* storage for compressed stream */
    size_t      bufsize;    /* byte size of compressed buffer */
    bitstream * stream;     /* bit stream to write to or read from */
    size_t      zfpsize;    /* byte size of compressed stream */

    PyObject * arrays;
    int        reversible = 0;
    double     rate       = -1.;
    double     accuracy   = -1.;
    static const char *kwlist[]  = {"data", "reversible", "rate", "accuracy", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwd, "O|idd", (char**)kwlist, &arrays, &reversible, &rate, &accuracy)) {
        PyErr_SetString(PyExc_SyntaxError,
                        "Wrong syntax. Right syntax : pack(array or list of arrays,[reversible = boolean, rate = "
                        "float, accuracy = float]");
        return NULL;
    }
    bool                         is_list = false;
    std::vector<PyArrayObject *> np_arrays;
    if (PyList_Check(arrays)) {
        is_list = true;
        np_arrays.reserve(PyList_Size(arrays));
        for (Py_ssize_t i = 0; i < PyList_Size(arrays); ++i)
            np_arrays.push_back((PyArrayObject *)PyList_GetItem(arrays, i));
    }
    else if (PyArray_Check(arrays))
        np_arrays.push_back((PyArrayObject *)arrays);
    else {
        PyErr_SetString(PyExc_TypeError, "First argument must be an array or a list of array");
        return NULL;
    }
    // Specify the type to compress :
    type = zfp_type_double;

    PyObject *compressed_list = PyList_New(np_arrays.size());
    int       ind_list        = 0;
    for (const auto &an_array : np_arrays) {
        int       ndims = PyArray_NDIM(an_array);
        npy_intp *dims  = PyArray_DIMS(an_array);
        bool is_c_order = false;
        if (PyArray_CHKFLAGS(an_array, NPY_ARRAY_C_CONTIGUOUS)) is_c_order = true;

        void* data = PyArray_DATA(an_array);
        if (ndims > 4)
        {
            PyErr_SetString(PyExc_ValueError, "Array must have less than five entries");
            return NULL;
        }
        if      (ndims == 1) field = zfp_field_1d(data, type, dims[0]); 
        else 
        {
            if (is_c_order)//_ C order _
            {
                if (ndims == 2) field = zfp_field_2d(data, type, dims[1], dims[0]);
                else if (ndims == 3) field = zfp_field_3d(data, type, dims[2], dims[1], dims[0]);
                else if (ndims == 4) field = zfp_field_4d(data, type, dims[3], dims[2], dims[1], dims[0]);
            }
            else //_ Fortran order _
            {
                if (ndims == 2) field = zfp_field_2d(data, type, dims[0], dims[1]);
                else if (ndims == 3) field = zfp_field_3d(data, type, dims[0], dims[1], dims[2]);
                else if (ndims == 4) field = zfp_field_4d(data, type, dims[0], dims[1], dims[2], dims[3]);
            }
        }
        /* allocate meta data for a compressed stream */
        zfp = zfp_stream_open(NULL);
        if (reversible==1) zfp_stream_set_reversible(zfp);
        if (rate >= 0)
            zfp_stream_set_rate(zfp, rate, type, ndims, 0);
        if (accuracy >= 0.)
        {
            zfp_stream_set_accuracy(zfp, accuracy);
        }
        /* allocate buffer for compressed data */
        bufsize = zfp_stream_maximum_size(zfp, field);
        buffer = malloc(bufsize);
        /* associate bit stream with allocated buffer */
        stream = stream_open(buffer, bufsize);
        zfp_stream_set_bit_stream(zfp, stream);
        zfp_stream_rewind(zfp);

         zfp_stream_set_omp_threads(zfp,1);
        /* compress array*/
        zfpsize = zfp_compress(zfp, field);
        if (!zfpsize) {
            PyErr_SetString(PyExc_RuntimeError, "Compression failed !");
            free(buffer);
            return NULL;
        }
        npy_intp cprsize = zfpsize;
        PyObject *shape = PyTuple_New(ndims);
        for (int i = 0; i < ndims; ++i) PyTuple_SET_ITEM(shape, i, PyLong_FromLong(long(dims[i])));
        PyObject *obj = PyTuple_New(3);
        PyTuple_SET_ITEM(obj, 0, shape);
        PyArrayObject* cdata = (PyArrayObject*)PyArray_SimpleNewFromData(1, &cprsize, NPY_BYTE, buffer);
        PyArray_ENABLEFLAGS(cdata, NPY_ARRAY_OWNDATA);
        PyTuple_SET_ITEM(obj, 1, (PyObject*)cdata);
        if (is_c_order)
        {
            Py_IncRef(Py_True);
            PyTuple_SetItem(obj, 2, Py_True); 
        }
        else
        {
            Py_IncRef(Py_False);
            PyTuple_SetItem(obj, 2, Py_False); 
        }
        PyList_SetItem(compressed_list, ind_list, obj);
        ind_list++;
        zfp_field_free(field);
        zfp_stream_close(zfp);
        stream_close(stream);
    }
    if (!is_list) {
        PyObject *array = PyList_GetItem(compressed_list, 0);
        Py_INCREF(array);
        Py_DECREF(compressed_list);
        return array;
    }
    return compressed_list;
}

static const char *decompress_doc = R"RAW(
Decompress a compressed array or a list of compressed arrays.
-------------------------------------------------------------
Usage : array = zfp.unpack(compressed_array[,options]);
      or
        lst_of_arrays = zfp.unpack([comp_arr1, comp_arr2,..., comp_arrb][,options])
      where :
        compressed_array (or comp_arr1 and so.) are compressed data (with sz),
      and options (optional arguments with keywords) as same signification than options in compress function.
)RAW";
PyObject *
py_decompress(PyObject *self, PyObject *args, PyObject* kwd)
{
    //int         status = 0; /* return value: 0 = success */
    zfp_type    type;       /* array scalar type */
    zfp_field * field=NULL; /* array meta data */
    zfp_stream *zfp;        /* compressed stream */
    void *      buffer;     /* storage for compressed stream */
    size_t      bufsize;    /* byte size of compressed buffer */
    bitstream * stream;     /* bit stream to write to or read from */
    //size_t      zfpsize;    /* byte size of compressed stream */

    PyObject * arrays;
    int        reversible = 0;
    double     rate       = -1.;
    double     accuracy   = -1.;
    static const char *kwlist[6]  = {"data", "reversible", "rate", "accuracy", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwd, "O|idd", (char **)kwlist, &arrays, &reversible, &rate, &accuracy)) {
        PyErr_SetString(PyExc_SyntaxError,
                        "Wrong syntax. Right syntax : pack(array or list of arrays,[reversible = boolean, rate = "
                        "(int,int), accuracy = float]");
        return NULL;
    }
    bool                         is_list = false;
    std::vector<PyArrayObject *> cpr_arrays;
    std::vector<PyObject*>       shapes;
    std::vector<bool>            is_c_order;
    if (PyList_Check(arrays)) {
        is_list              = true;
        Py_ssize_t list_size = PyList_Size(arrays);
        cpr_arrays.reserve(list_size);
        is_c_order.reserve(list_size);
        for (Py_ssize_t i = 0; i < list_size; ++i) {
            PyObject *tuple = PyList_GetItem(arrays, i);
            if (!PyTuple_Check(tuple)) {
                PyErr_SetString(PyExc_TypeError, "The values of the list must be tuple as (shape, compressed data)");
                return NULL;
            }
            PyObject *array = PyTuple_GetItem(tuple, 1);
            if (!PyArray_Check(array)) {
                PyErr_SetString(PyExc_TypeError, "Second value of the tuple must be an array with compressed data");
                return NULL;
            }
            cpr_arrays.push_back((PyArrayObject *)array);
            PyObject *shp = PyTuple_GetItem(tuple, 0);
            if (!PyTuple_Check(shp)) {
                PyErr_SetString(PyExc_TypeError, "A shape must be a tuple of integers");
                return NULL;
            }
            shapes.push_back(shp);
            PyObject *py_is_c_order = PyTuple_GetItem(tuple, 2);
            if (py_is_c_order == Py_True)
                is_c_order.push_back(1);
            else//_ Fortran order _
            {
                assert(py_is_c_order == Py_False);
                is_c_order.push_back(0);
            }

        }
    }
    else if (PyTuple_Check(arrays)) {
        PyObject *shape = PyTuple_GetItem(arrays, 0);
        PyObject *array = PyTuple_GetItem(arrays, 1);
        if (!PyArray_Check(array)) {
            PyErr_SetString(PyExc_TypeError, "Second value of tuple must be an array with compressed data");
            return NULL;
        }
        cpr_arrays.push_back((PyArrayObject *)array);
        if (!PyTuple_Check(shape)) {
            PyErr_SetString(PyExc_TypeError,
                            "First value of tuple must be a tuple containing the shape of the original array");
            return NULL;
        }
        shapes.push_back(shape);
        PyObject *py_is_c_order = PyTuple_GetItem(arrays, 2);
        if (py_is_c_order == Py_True)
            is_c_order.push_back(1);
        else//_ Fortran order _
        {
            assert(py_is_c_order == Py_False);
            is_c_order.push_back(0);
        }
    }
    else {
        PyErr_SetString(PyExc_TypeError, "First argument must be an compressed array or a list of compressed arrays");
        return NULL;
    }
    // Specify the type to compress :
    type = zfp_type_double;
    PyObject *lst_out_arrays;
    lst_out_arrays = PyList_New(cpr_arrays.size());
    for (size_t i = 0; i < cpr_arrays.size(); ++i) {
        PyArrayObject* array = cpr_arrays[i];
        PyObject* shape = shapes[i];
        npy_intp  dims[4];
        Py_ssize_t ndim = PyTuple_Size(shape);
        if (ndim > 4) {
            PyErr_SetString(PyExc_ValueError, "Shape must have only four or less elements");
            return NULL;
        }
        for ( Py_ssize_t j = 0; j < ndim; ++j )
        {
            PyObject *py_dim = PyTuple_GetItem(shape, j);
            if (PyLong_Check(py_dim) == false && PyInt_Check(py_dim) == false) 
            {
                PyErr_SetString(PyExc_TypeError, "Values in shape must be integers");
                return NULL;
            }
            long dim = PyLong_AsLong(py_dim);
            dims[j] = dim;
        }
        bool is_c = (is_c_order[i]==1);
        PyArrayObject *py_array;
        if (is_c)  py_array = (PyArrayObject *)PyArray_EMPTY(ndim, dims, NPY_DOUBLE,0);
        else py_array = (PyArrayObject *)PyArray_EMPTY(ndim, dims, NPY_DOUBLE,1);
        unsigned char *py_array_data = (unsigned char *)PyArray_DATA(py_array);

        if (ndim == 1) field = zfp_field_1d(py_array_data, type, dims[0]);
        else 
        if (is_c)
        {
            if (ndim == 2) field = zfp_field_2d(py_array_data, type, dims[1], dims[0]);
            else if (ndim == 3) field = zfp_field_3d(py_array_data, type, dims[2], dims[1], dims[0]);
            else if (ndim == 4) field = zfp_field_4d(py_array_data, type, dims[3], dims[2], dims[1], dims[0]);
        }
        else//_ Fortran order _
        {
            if (ndim == 2) field = zfp_field_2d(py_array_data, type, dims[1], dims[0]);
            else if (ndim == 3) field = zfp_field_3d(py_array_data, type, dims[0], dims[1], dims[2]);
            else if (ndim == 4) field = zfp_field_4d(py_array_data, type, dims[0], dims[1], dims[2], dims[3]);            
        }
        /* allocate meta data for a compressed stream */
        zfp = zfp_stream_open(NULL);
        if (reversible==1) zfp_stream_set_reversible(zfp);
        if (rate >= 0.)
        {
            zfp_stream_set_rate(zfp, rate, type, ndim, 0);
        }
        if (accuracy >= 0.)
        {
            zfp_stream_set_accuracy(zfp, accuracy);
        }
        /* allocate buffer for compressed data */
        bufsize = PyArray_Size((PyObject*)array);
        buffer = PyArray_DATA(array);
        /* associate bit stream with allocated buffer */
        stream = stream_open(buffer, bufsize);
        zfp_stream_set_bit_stream(zfp, stream);
        zfp_stream_rewind(zfp);
        if (!zfp_decompress(zfp, field))
        {
            PyErr_SetString(PyExc_RuntimeError, "Failed to decompress data for an array !");
            return NULL;
        }
        zfp_field_free(field);
        zfp_stream_close(zfp);
        stream_close(stream);
        if (!is_list) {
            Py_DecRef(lst_out_arrays);
            /* clean up */
            return (PyObject *)py_array;
        }
        /* clean up */
        PyList_SetItem(lst_out_arrays, i, (PyObject *)py_array);
    }
    return lst_out_arrays;
}
} // namespace zfp

static PyMethodDef Pycompressor_zfp[] = {{"pack", (PyCFunction)zfp::py_compress, METH_VARARGS | METH_KEYWORDS, zfp::compress_doc},
                                        {"unpack", (PyCFunction)zfp::py_decompress, METH_VARARGS | METH_KEYWORDS, zfp::decompress_doc},
                                        {NULL, NULL,0,NULL}};


static const char *       module_doc = R"RAW(
Python interface for the scientific dataset compressor zfp
----------------------------------------------------------
Two functions to use :
  - pack   : compress an array or a list of arrays containing double values
  - unpack : decompress an array or a list of arrays, providing shapes for decompressed arrays
)RAW";

#if PY_MAJOR_VERSION >= 3
// =====================================================================
static struct PyModuleDef moduledef  = {PyModuleDef_HEAD_INIT,
                                       "czfp",
                                       module_doc,
                                       -1, // sizeof(struct module_state),
                                       Pycompressor_zfp,
                                       NULL,
                                       NULL, // myextension_traverse,
                                       NULL, // myextension_clear,
                                       NULL};

PyMODINIT_FUNC
PyInit_czfp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (m == NULL) return NULL;
    /* Tres important : initialise numpy afin de pouvoir l'utiliser ici !!!! */
    import_array();

    return m;
}
#else
PyMODINIT_FUNC
initczfp(void)
{
    PyObject* m = Py_InitModule3("czfp", Pycompressor_zfp, module_doc);
    if (m == NULL) return;
    /* Tres important : initialise numpy afin de pouvoir l'utiliser ici !!!! */
    import_array();

}
#endif
