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
#include <array>
#include <cstdint>
#include <iostream>
#include "fpc.h"
#include "compressor.h"

namespace K_COMPRESSOR
{
PyObject* py_fpc_compress(PyObject *self, PyObject *args)
{
    PyObject *arrays;
    if (!PYPARSETUPLE_(args, O_, &arrays)) 
    {
        PyErr_SetString(PyExc_SyntaxError,
                        "pack: wrong syntax. Right syntax: packFpc(array or list of arrays");
        return NULL;
    }
    bool is_list = false;
    std::vector<PyArrayObject *> np_arrays;
    if (PyList_Check(arrays)) 
    {
        is_list = true;
        np_arrays.reserve(PyList_Size(arrays));
        for (Py_ssize_t i = 0; i < PyList_Size(arrays); ++i)
            np_arrays.push_back((PyArrayObject *)PyList_GetItem(arrays, i));
    }
    else if (PyArray_Check(arrays))
        np_arrays.push_back((PyArrayObject *)arrays);
    else 
    {
        PyErr_SetString(PyExc_TypeError, "pack: first argument must be an array or a list of array");
        return NULL;
    }

    PyObject *compressed_list = PyList_New(np_arrays.size());
    for ( size_t i = 0; i < np_arrays.size(); ++i )
    {
        //= On recupere un tableau numpy et ses dimensions et sa taille
        PyArrayObject* an_array = np_arrays[i];
        int       ndims = PyArray_NDIM(an_array);
        npy_intp *dims  = PyArray_DIMS(an_array);
        std::size_t an_array_length = PyArray_SIZE(an_array);
        double* array_data = (double*)PyArray_DATA(an_array);
        bool is_c_order = false;
        if (PyArray_CHKFLAGS(an_array, NPY_ARRAY_C_CONTIGUOUS)) is_c_order = true;
        
        // encode with fpc
        E_Int outSize = FPC_UPPER_BOUND(an_array_length);
        uint8_t* out_compressed = new uint8_t [outSize];
        fpc_context_t ctx;
        uint64_t fcm[FPC_TABLE_SIZE_DEFAULT];
        uint64_t dfcm[FPC_TABLE_SIZE_DEFAULT];
        ctx.fcm_size = FPC_TABLE_SIZE_DEFAULT;
        ctx.fcm = fcm;
        ctx.dfcm_size = FPC_TABLE_SIZE_DEFAULT;
        ctx.dfcm = dfcm;
        ctx.hash_args = FPC_DEFAULT_HASH_ARGS;
        ctx.seed = 0.0;
        memset(fcm, 0, sizeof(fcm));
        memset(dfcm, 0, sizeof(dfcm));
        //for (E_Int i = 0; i < an_array_length; i++) printf("%g ", array_data[i]);
        
        E_Int size = fpc_encode(&ctx, array_data, an_array_length, out_compressed);
        //printf("compression: init=%d, compressed=%d, reserved=%d\n", an_array_length*8, size, outSize);
        //printf("outcompress %d\n", size);
        //for (E_Int i = 0; i < size; i++) printf("%u ", out_compressed[i]);
        
        // Decompression:
        //memset(fcm, 0, sizeof(fcm));
        //memset(dfcm, 0, sizeof(dfcm));
        //double* decompressed = new double [an_array_length];
        //fpc_decode(&ctx, out_compressed, decompressed, an_array_length);
        //delete [] decompressed;
        //for (E_Int i = 0; i < an_array_length; i++) printf("%f ", decompressed[i]);

        //= On prepare l'objet compresse
        PyObject *shape = PyTuple_New(ndims);
        for (E_Int j = 0; j < ndims; ++j) PyTuple_SET_ITEM(shape, j, PyLong_FromLong(long(dims[j])));
        PyObject *obj = PyTuple_New(3);
        PyTuple_SET_ITEM(obj, 0, shape);
    
        //= Reservation memoire pour le buffer compresse
        npy_intp sz = size;
        PyArrayObject* cpr_arr = (PyArrayObject*)PyArray_SimpleNew(1, &sz, NPY_BYTE);
        std::uint8_t* buffer = (std::uint8_t*)PyArray_DATA(cpr_arr);
        # pragma omp parallel for
        for (npy_intp ibyte = 0; ibyte < size; ibyte++)
        {
            buffer[ibyte] = out_compressed[ibyte];
        }
        delete [] out_compressed;

        //= On rajoute le tableau au tuple (shape,buffer)
        PyTuple_SET_ITEM(obj, 1, (PyObject*)cpr_arr);
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
        PyList_SetItem(compressed_list, i, obj);
    }
    if (!is_list) 
    {
        //= Si ce n'etait pas une liste au depart, on retourne un tableau
        PyObject* array = PyList_GetItem(compressed_list, 0);
        Py_INCREF(array);
        Py_DECREF(compressed_list);
        return array;
    }
    //= Sinon on retourne une liste de tableaux
    return compressed_list;
}

PyObject* py_fpc_uncompress(PyObject *self, PyObject *args)
{
    PyObject *cpr_arrays;
    if (!PYPARSETUPLE_(args, O_, &cpr_arrays)) 
    {
        PyErr_SetString(PyExc_SyntaxError, "Wrong syntax. Right syntax : unpackCellN(array or list of compressed arrays");
        return NULL;
    }
    bool is_list = false;
    std::vector<PyArrayObject *> np_cpr_arrays;
    std::vector<std::vector<npy_intp>> shape_arrays;
    std::vector<bool> is_c_order;
    if (PyList_Check(cpr_arrays)) 
    {
        is_list = true;
        Py_ssize_t list_size = PyList_Size(cpr_arrays);
        np_cpr_arrays.reserve(list_size);
        shape_arrays.resize(list_size);
        is_c_order.reserve(list_size);
        for (Py_ssize_t i = 0; i < list_size; ++i) 
        {
            PyObject *tuple = PyList_GetItem(cpr_arrays, i);
            if (!PyTuple_Check(tuple)) {
                PyErr_SetString(PyExc_TypeError, "The values of the list must be tuple as (shape, compressed data)");
                return NULL;
            }
            PyObject *array = PyTuple_GetItem(tuple, 1);
            if (!PyArray_Check(array)) 
            {
                PyErr_SetString(PyExc_TypeError, "Second value of the tuple must be an array with compressed data");
                return NULL;
            }
            np_cpr_arrays.push_back((PyArrayObject *)array);
            PyObject *shp = PyTuple_GetItem(tuple, 0);
            if (!PyTuple_Check(shp)) 
            {
                PyErr_SetString(PyExc_TypeError, "A shape must be a tuple of integers");
                return NULL;
            }
            PyObject* py_is_c_order = PyTuple_GetItem(tuple,2);
            if (py_is_c_order == Py_True)
            {
                is_c_order.push_back(true);
            }
            else
            {
                is_c_order.push_back(false);
            }
            Py_ssize_t dimshape = PyTuple_Size(shp);
            shape_arrays[i].reserve(dimshape);
            for (Py_ssize_t j = 0; j < dimshape; ++j) 
            {
                PyObject *py_dim = PyTuple_GetItem(shp, j);
                if (PyLong_Check(py_dim) == false && PyInt_Check(py_dim) == false) 
                {
                    PyErr_SetString(PyExc_TypeError, "Values in shape must be integers");
                    return NULL;
                }
                long dim = PyLong_AsLong(py_dim);
                shape_arrays[i][j] = npy_intp(dim);
            }
        } // End for (i= 0; i < list_size; ...)
    }     // Fin du cas si c'est une liste
    else if (PyTuple_Check(cpr_arrays)) 
    {
        PyObject *shape = PyTuple_GetItem(cpr_arrays, 0);
        PyObject *array = PyTuple_GetItem(cpr_arrays, 1);
        PyObject *py_is_c_order = PyTuple_GetItem(cpr_arrays, 2);
        if (!PyArray_Check(array)) 
        {
            PyErr_SetString(PyExc_TypeError, "Second value of tuple must be an array with compressed data");
            return NULL;
        }
        np_cpr_arrays.push_back((PyArrayObject*)array);
        shape_arrays.resize(1);
        if (!PyTuple_Check(shape)) 
        {
            PyErr_SetString(PyExc_TypeError,
                            "First value of tuple must be a tuple containing the shape of the original array");
            return NULL;
        }
        if (py_is_c_order == Py_True)
        {
            is_c_order.push_back(true);
        }
        else
        {
            is_c_order.push_back(false);
        }

        Py_ssize_t dimshape = PyTuple_Size(shape);
        shape_arrays[0].resize(dimshape);
        for (Py_ssize_t j = 0; j < dimshape; ++j) 
        {
            PyObject *py_dim = PyTuple_GetItem(shape, j);
            if (!PyLong_Check(py_dim)) 
            {
                PyErr_SetString(PyExc_TypeError, "Values in shape must be integers");
                return NULL;
            }
            long dim = PyLong_AsLong(py_dim);
            shape_arrays[0][j] = npy_intp(dim);
        }
    } // Fin du cas avec un seul tableau
    else 
    {
        PyErr_SetString(PyExc_TypeError, "First argument must be an compressed array or a list of compressed arrays");
        return NULL;
    }

    PyObject *lst_out_arrays;
    lst_out_arrays = PyList_New(np_cpr_arrays.size());
    for (size_t i = 0; i < np_cpr_arrays.size(); ++i) 
    {
        npy_intp  dims[5];
        E_Int     ndim;
        ndim = shape_arrays[i].size();
        for (E_Int j = 0; j < ndim; ++j) 
        { 
            dims[j] = shape_arrays[i][j]; 
        }
        PyArrayObject *py_array;
        bool is_c_ord = is_c_order[i];
        if (is_c_ord) py_array = (PyArrayObject *)PyArray_EMPTY(ndim, dims, NPY_DOUBLE, 0);
        else
        {
            py_array = (PyArrayObject *)PyArray_EMPTY(ndim, dims, NPY_DOUBLE, 1);
        }
        double* py_array_data = (double*)PyArray_DATA(py_array);
        //std::size_t cpr_length = PyArray_SIZE(np_cpr_arrays[i]);
        std::size_t array_length = PyArray_SIZE(py_array);
        std::uint8_t* cpr_data = (std::uint8_t*)PyArray_DATA(np_cpr_arrays[i]);

        // decode with fpc
        fpc_context_t ctx;
        uint64_t fcm[FPC_TABLE_SIZE_DEFAULT];
        uint64_t dfcm[FPC_TABLE_SIZE_DEFAULT];
        ctx.fcm_size = FPC_TABLE_SIZE_DEFAULT;
        ctx.fcm = fcm;
        ctx.dfcm_size = FPC_TABLE_SIZE_DEFAULT;
        ctx.dfcm = dfcm;
        ctx.hash_args = FPC_DEFAULT_HASH_ARGS;
        ctx.seed = 0.0;
        memset(fcm, 0, sizeof(fcm));
        memset(dfcm, 0, sizeof(dfcm));
        //printf("decompression: compressed size=%d uncomp=%d\n", cpr_length, array_length*8);
        //printf("cprdata %d\n", cpr_length);
        //for (E_Int i = 0; i < cpr_length; i++) printf("%u ", cpr_data[i]);
        
        //for (E_Int i = 0; i < array_length; i++) py_array_data[i] = 0.;
        fpc_decode(&ctx, cpr_data, py_array_data, array_length);
        //double* decompressed = new double [array_length];
        //fpc_decode(&ctx, cpr_data, decompressed, array_length);
        //for (E_Int j = 0; j < array_length; j++) printf("%f ", py_array_data[j]);
                
        if (!is_list)
        {
            Py_DECREF(lst_out_arrays);
            return (PyObject *)py_array;
        }
        PyList_SetItem(lst_out_arrays, i, (PyObject *)py_array);
    }
    
    return lst_out_arrays;
}
}// Fin namespace K_COMPRESSOR
