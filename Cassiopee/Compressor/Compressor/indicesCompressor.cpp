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
#include "zstd/zstd.h"
#include "compressor.h"

namespace K_COMPRESSOR
{
PyObject* py_indices_compress(PyObject* self, PyObject* args)
{
    PyObject *arrays;
    if (!PYPARSETUPLE_(args, O_, &arrays)) 
    {
        PyErr_SetString(PyExc_SyntaxError,
                        "Wrong syntax. Right syntax: packIndices((n,array) or list of (n,array)s");
        return NULL;
    }
    bool is_list = false;
    std::vector<PyArrayObject *> np_arrays;
    std::vector<long> nb_verts_per_elt_arr;
    std::vector<bool> is_c_order;
    if (PyList_Check(arrays)) 
    {
        is_list = true;
        np_arrays.reserve(PyList_Size(arrays));
        nb_verts_per_elt_arr.reserve(PyList_Size(arrays));
        is_c_order.reserve(PyList_Size(arrays));
        for (Py_ssize_t i = 0; i < PyList_Size(arrays); ++i)
        {
            PyObject* tuple_indices = PyList_GetItem(arrays, i);
            if (!PyTuple_Check(tuple_indices))
            {
                PyErr_SetString(PyExc_TypeError, "Wrong syntax : each element of the list must be a tuple as (nb_verts_per_elt, indices)");
                return NULL;
            }
            PyObject* py_nb_v_p_e = PyTuple_GetItem(tuple_indices, 0);
            if (!PyLong_Check(py_nb_v_p_e))
            {
                PyErr_SetString(PyExc_TypeError, "Wrong syntax : First element of each tuple must be an integer for number of vertices per element");
                return NULL;                
            }
            long nb_v_per_elt = PyLong_AsLong(py_nb_v_p_e);
            PyObject* py_indices = PyTuple_GetItem(tuple_indices, 1);
            if (!PyArray_Check(py_indices))
            {
                PyErr_SetString(PyExc_TypeError, "Wront syntax: The second element of each tuple must be an array of indices");
                return NULL;
            }
            if (PyArray_CHKFLAGS((PyArrayObject*)py_indices, NPY_ARRAY_C_CONTIGUOUS)) 
                is_c_order.push_back(true);
            else 
                is_c_order.push_back(false);

            np_arrays.push_back((PyArrayObject *)py_indices);
            nb_verts_per_elt_arr.push_back(nb_v_per_elt);
        }
    }
    else if (PyTuple_Check(arrays))
    {
        PyObject* py_nb_v_p_e = PyTuple_GetItem(arrays, 0);
        if (!PyLong_Check(py_nb_v_p_e))
        {
            PyErr_SetString(PyExc_TypeError, "Wrong syntax: First element of each tuple must be an integer for number of vertices per element");
            return NULL;                
        }
        long nb_v_per_elt = PyLong_AsLong(py_nb_v_p_e);
        PyObject* py_indices = PyTuple_GetItem(arrays, 1);
        if (!PyArray_Check(py_indices))
        {
            PyErr_SetString(PyExc_TypeError, "Wrong syntax: Second element of each tuple must be an array of indices");
            return NULL;                            
        }
        if (PyArray_CHKFLAGS((PyArrayObject*)py_indices, NPY_ARRAY_C_CONTIGUOUS)) 
            is_c_order.push_back(true);
        else 
            is_c_order.push_back(false);
        np_arrays.push_back((PyArrayObject *)py_indices);
        nb_verts_per_elt_arr.push_back(nb_v_per_elt);
    }
    else {
        PyErr_SetString(PyExc_TypeError, "First argument must be a tuple or a list of tuples");
        return NULL;
    }

    PyObject *compressed_list = PyList_New(np_arrays.size());
    for ( size_t i = 0; i < np_arrays.size(); ++i )
    {
        //= On récupère le nombre de sommets par élément :
        long nb_v_per_elt = nb_verts_per_elt_arr[i];
        //= On récupère un tableau numpy et ses dimensions et sa taille
        PyArrayObject* an_array = np_arrays[i];
        int       ndims = PyArray_NDIM(an_array);
        if (ndims != 1)
        {
            PyErr_SetString(PyExc_TypeError, "Wrong data: The indices arrays must have one dimension only");
            return NULL;                
        }
        std::size_t  an_array_length = PyArray_SIZE(an_array);
        if ( an_array_length%nb_v_per_elt != 0 )
        {
            PyErr_SetString(PyExc_TypeError, "Wrong data: The length of the indice array is not a multiple of number of vertices per element !");
            return NULL;                
        }

        E_Int* array_data = (E_Int*)PyArray_DATA(an_array);
        //= On prépare l'objet décrivant la compression du cellN
        PyObject *py_nb_v_per_e = PyLong_FromLong(nb_v_per_elt);
        PyObject *obj = PyTuple_New(4);
        PyTuple_SET_ITEM(obj, 0, py_nb_v_per_e);
        //= Rajout de la taille du tableau au tuple
        PyObject *py_length = PyLong_FromLong(long(an_array_length));
        PyTuple_SET_ITEM(obj, 1, py_length);
        //= On transforme maintenant le tableau pour avoir des indices relatifs à un indice de référence
        std::vector<E_Int> rel_indices(an_array_length);
        for (E_Int j = 0; j < E_Int(an_array_length); j += E_Int(nb_v_per_elt) )
        {
            rel_indices[j] = array_data[j];
            for ( E_Int k = 1; k < nb_v_per_elt; ++k )
                rel_indices[k+j] = array_data[k+j] - array_data[j];
        }
        //= Compression de rel_indices with zstd :
        std::size_t maxSize = ZSTD_compressBound(an_array_length*sizeof(E_Int));
        std::vector<std::uint8_t> compressed_data(maxSize);
        std::size_t szCompress = ZSTD_compress(compressed_data.data(), maxSize, rel_indices.data(), an_array_length*sizeof(E_Int), 
                                               ZSTD_maxCLevel());
        npy_intp sz = npy_intp(szCompress);
        PyArrayObject* cpr_arr = (PyArrayObject*)PyArray_SimpleNew(1, &sz, NPY_BYTE);
        std::uint8_t* buffer = (std::uint8_t*)PyArray_DATA(cpr_arr);
        std::copy(compressed_data.data(), compressed_data.data()+szCompress, buffer);
        //= On rajoute le tableau au tuple (nb_v_p_e,size,buffer)
        PyTuple_SET_ITEM(obj, 2, (PyObject*)cpr_arr);
        //$ Rajout du flag pour savoir si on est en fortran order ou pas
        if (is_c_order[i])
        {
            Py_IncRef(Py_True);
            PyTuple_SetItem(obj, 3, Py_True);
        }
        else
        {
            Py_IncRef(Py_False);
            PyTuple_SetItem(obj, 3, Py_False);
        }
        //= Rajout du tuple à la liste éventuelle
        PyList_SetItem(compressed_list, i, obj);
    }
    if (!is_list) {
        //= Si ce n'était pas une liste au départ, on retourne un tableau
        PyObject* array = PyList_GetItem(compressed_list, 0);
        Py_INCREF(array);
        Py_DECREF(compressed_list);
        return array;
    }
    //= Sinon on retourne une liste de tableaux
    return compressed_list;
}
//@ ________________________________________________________________________________________________
PyObject* py_indices_uncompress(PyObject *self, PyObject *args)
{
    PyObject *cpr_arrays;
    if (!PYPARSETUPLE_(args, O_, &cpr_arrays)) 
    {
        PyErr_SetString(PyExc_SyntaxError, "Wrong syntax. Right syntax: unpackIndices(array or list of compressed arrays");
        return NULL;
    }
    std::vector<PyArrayObject *> np_cpr_arrays;
    std::vector<npy_intp> length_arrays;
    std::vector<E_Int>    nb_verts_per_elt;
    std::vector<bool>     is_c_order;
    if (PyList_Check(cpr_arrays)) {
        Py_ssize_t list_size = PyList_Size(cpr_arrays);
        np_cpr_arrays.reserve(list_size);
        length_arrays.reserve(list_size);
        nb_verts_per_elt.reserve(list_size);
        is_c_order.reserve(list_size);
        for (Py_ssize_t i = 0; i < list_size; ++i) {
            PyObject *tuple = PyList_GetItem(cpr_arrays, i);
            if (!PyTuple_Check(tuple)) {
                PyErr_SetString(PyExc_TypeError, "The values of the list must be tuple as (nb_v_per_e, length, compressed data)");
                return NULL;
            }
            PyObject* flag = PyTuple_GetItem(tuple,3);
            if (flag == Py_True)
                is_c_order.push_back(true);
            else
                is_c_order.push_back(false);
            PyObject *array = PyTuple_GetItem(tuple, 2);
            if (!PyArray_Check(array)) {
                PyErr_SetString(PyExc_TypeError, "Third value of the tuple must be an array with compressed data");
                return NULL;
            }
            np_cpr_arrays.push_back((PyArrayObject *)array);
            PyObject *lgth = PyTuple_GetItem(tuple, 1);
            if (!PyLong_Check(lgth)) {
                PyErr_SetString(PyExc_TypeError, "Second element of the tuple must be an integer");
                return NULL;
            }
            length_arrays.push_back(npy_intp(PyLong_AsLong(lgth)));            
            PyObject* py_nb_v_p_e = PyTuple_GetItem(tuple,0);
            if (!PyLong_Check(py_nb_v_p_e)) {
                PyErr_SetString(PyExc_TypeError, "First element of the tuple must be an integer");
                return NULL;
            }
            nb_verts_per_elt.push_back(E_Int(PyLong_AsLong(py_nb_v_p_e)));
        } // End for (i= 0; i < list_size; ...)
    }     // Fin du cas si c'est une liste
    else if (PyTuple_Check(cpr_arrays)) {
        PyObject *py_nb_v_p_e = PyTuple_GetItem(cpr_arrays, 0);
        PyObject *lgth        = PyTuple_GetItem(cpr_arrays, 1);
        PyObject *array       = PyTuple_GetItem(cpr_arrays, 2);
        PyObject *flag        = PyTuple_GetItem(cpr_arrays, 3);
        if (!PyArray_Check(array)) {
            PyErr_SetString(PyExc_TypeError, "Second value of tuple must be an array with compressed data");
            return NULL;
        }
        np_cpr_arrays.push_back((PyArrayObject *)array);
        if (!PyLong_Check(lgth)) {
            PyErr_SetString(PyExc_TypeError, "Second element of the tuple must be an integer");
            return NULL;
        }
        length_arrays.push_back(npy_intp(PyLong_AsLong(lgth)));            
        if (!PyLong_Check(py_nb_v_p_e)) {
            PyErr_SetString(PyExc_TypeError, "First element of the tuple must be an integer");
            return NULL;
        }
        nb_verts_per_elt.push_back(E_Int(PyLong_AsLong(py_nb_v_p_e)));
        if (flag == Py_True)
            is_c_order.push_back(true);
        else
            is_c_order.push_back(false);
    } // Fin du cas avec un seul tableau
    else {
        PyErr_SetString(PyExc_TypeError, "First argument must be an compressed array or a list of compressed arrays");
        return NULL;
    }

    PyObject *lst_out_arrays;
    lst_out_arrays = PyList_New(np_cpr_arrays.size());
    for (size_t i = 0; i < np_cpr_arrays.size(); ++i) {
        npy_intp  length = length_arrays[i];
        E_Int ordering = (is_c_order[i] ? 0 : 1);
        PyArrayObject *py_array = (PyArrayObject *)PyArray_EMPTY(1, &length, E_NPY_INT, ordering);
        E_Int *py_array_data = (E_Int *)PyArray_DATA(py_array);
        std::uint8_t *cpr_data = (std::uint8_t *)PyArray_DATA(np_cpr_arrays[i]);
        std::size_t compressed_size = PyArray_SIZE((PyArrayObject*)np_cpr_arrays[i]);
        std::vector<E_Int> uncompressed_data(length);
        std::size_t lgth_uncompress = ZSTD_decompress(uncompressed_data.data(), length*sizeof(E_Int), 
                                                      cpr_data, compressed_size);
        if (lgth_uncompress != length*sizeof(E_Int))
        {
            PyErr_SetString(PyExc_RuntimeError, "Uncompatible size between indice array and decompressed data size !");
            return NULL;            
        }
        //= Reconstitution des indices :
        for ( E_Int j = 0; j < E_Int(length); j += nb_verts_per_elt[i] )
        {
            py_array_data[j] = uncompressed_data[j];
            for ( E_Int k = 1; k < nb_verts_per_elt[i]; ++k )
                py_array_data[j+k] = uncompressed_data[j+k] + py_array_data[j];
        }
        PyList_SetItem(lst_out_arrays, i, (PyObject *)py_array);
    }
    return lst_out_arrays;
}
}
