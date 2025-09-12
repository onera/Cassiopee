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
#include <cstdint>
#include <iostream>
#include "compressor.h"

namespace K_COMPRESSOR
{
PyObject* py_cellN_compress(PyObject *self, PyObject *args)
{
    PyObject *arrays;
    if (!PYPARSETUPLE_(args, O_, &arrays)) 
    {
        PyErr_SetString(PyExc_SyntaxError,
                        "pack: wrong syntax. Right syntax: packCellN(array or list of arrays");
        return NULL;
    }
    bool  is_list = false;
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
        E_Int ndims = PyArray_NDIM(an_array);
        npy_intp *dims  = PyArray_DIMS(an_array);
        std::size_t  an_array_length = PyArray_SIZE(an_array);
        double* array_data = (double*)PyArray_DATA(an_array);
        bool is_c_order = false;
        if (PyArray_CHKFLAGS(an_array, NPY_ARRAY_C_CONTIGUOUS)) is_c_order = true;
        //= On prepare l'objet decrivant la compression du cellN
        PyObject *shape = PyTuple_New(ndims);
        for (E_Int j = 0; j < ndims; ++j) PyTuple_SET_ITEM(shape, j, PyLong_FromLong(long(dims[j])));
        PyObject *obj = PyTuple_New(3);
        PyTuple_SET_ITEM(obj, 0, shape);
        //= Reservation mémoire pour le buffer compressé de cellN
        npy_intp sz = npy_intp(an_array_length+3)/4;
        PyArrayObject* cpr_arr = (PyArrayObject*)PyArray_SimpleNew(1, &sz, NPY_BYTE);
        std::uint8_t* buffer = (std::uint8_t*)PyArray_DATA(cpr_arr);
#       pragma omp parallel for        
        for (std::size_t ibyte = 0; ibyte < an_array_length/4; ++ibyte)
        {
            std::size_t ind = 4*ibyte;
            std::uint8_t c1 = std::uint8_t(array_data[ind+0])&3, 
                         c2 = std::uint8_t(array_data[ind+1])&3,
                         c3 = std::uint8_t(array_data[ind+2])&3, 
                         c4 = std::uint8_t(array_data[ind+3])&3;
            buffer[ibyte] =  c1 + (c2<<2) + (c3<<4) + (c4<<6);
        }
        //= Il faut traiter le cas où le tableau a une longueur non divisible
        //- par quatre :
        std::size_t remainder = an_array_length&3;
        if (remainder > 0)
        {
            std::size_t ind = an_array_length - remainder;
            std::uint8_t c1 = std::uint8_t(array_data[ind])&3;
            std::uint8_t c2 = 3, c3 = 3, c4 = 3;
            if (remainder > 1) c2 = std::uint8_t(array_data[ind+1])&3;
            if (remainder > 2) c3 = std::uint8_t(array_data[ind+2])&3;
            buffer[sz-1] = c1 + (c2<<2) + (c3<<4) + (c4<<6);
        }
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
        //= Si ce n'était pas une liste au départ, on retourne un tableau
        PyObject* array = PyList_GetItem(compressed_list, 0);
        Py_INCREF(array);
        Py_DECREF(compressed_list);
        return array;
    }
    //= Sinon on retourne une liste de tableaux
    return compressed_list;
}

PyObject* py_cellN_uncompress(PyObject *self, PyObject *args)
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
            for (Py_ssize_t j = 0; j < dimshape; ++j) {
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
        if (!PyArray_Check(array)) {
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
        E_Int ndim;
        ndim = shape_arrays[i].size();
        for (E_Int j = 0; j < ndim; ++j) 
        { 
            dims[j] = shape_arrays[i][j]; 
        }
        PyArrayObject *py_array;
        bool is_c_ord = is_c_order[i];
        if (is_c_ord)  py_array = (PyArrayObject *)PyArray_EMPTY(ndim, dims, NPY_DOUBLE,0);
        else 
        {
            py_array = (PyArrayObject *)PyArray_EMPTY(ndim, dims, NPY_DOUBLE,1);
        }
        double* py_array_data = (double*)PyArray_DATA(py_array);
        std::size_t cpr_length = PyArray_SIZE(np_cpr_arrays[i]);
        //std::size_t   array_length  = PyArray_SIZE(py_array);
        std::uint8_t* cpr_data = (std::uint8_t*)PyArray_DATA(np_cpr_arrays[i]);
#       pragma omp parallel for        
        for (std::size_t ibyte = 0; ibyte < cpr_length-1; ++ibyte)
        {
            std::size_t ind = 4*ibyte;
            std::int8_t byte = cpr_data[ibyte];
            py_array_data[ind + 0] = (byte   )&3;
            py_array_data[ind + 1] = (byte>>2)&3;
            py_array_data[ind + 2] = (byte>>4)&3;
            py_array_data[ind + 3] = (byte>>6)&3;
        }
        std::size_t ind = 4*(cpr_length-1);
        std::int8_t byte = cpr_data[cpr_length-1];
        py_array_data[ind+0] = (byte   )&3;
        if (((byte>>2)&3) != 3) py_array_data[ind+1] = (byte>>2)&3;
        if (((byte>>4)&3) != 3) py_array_data[ind+2] = (byte>>4)&3;
        if (((byte>>6)&3) != 3) py_array_data[ind+3] = (byte>>6)&3;
        if (!is_list) 
        {
            Py_DecRef(lst_out_arrays);
            return (PyObject *)py_array;
        }
        PyList_SetItem(lst_out_arrays, i, (PyObject *)py_array);
    }
    return lst_out_arrays;
}
}// Fin namespace K_COMPRESSOR
