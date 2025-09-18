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
#include "zstd/zstd.h"
#include "compressor.h"

namespace K_COMPRESSOR
{
PyObject* py_ngon_indices_compress(PyObject* self, PyObject* args)
{
    PyObject *arrays, *indices = nullptr;
    if (!PyArg_ParseTuple(args, "O|O", &arrays, &indices)) 
    {
        PyErr_SetString(PyExc_SyntaxError,
                        "Wrong syntax. Right syntax : packIndices((n,array) or list of (n,array)s");
        return NULL;
    }
    bool is_list = false;
    std::vector<PyArrayObject *> np_beg_indices;
    std::vector<PyArrayObject *> np_indices;
    if (PyList_Check(arrays)) {
        is_list = true;
        Py_ssize_t nb_set_of_indices = PyList_Size(arrays);
        np_beg_indices.reserve(nb_set_of_indices);
        np_indices.reserve(nb_set_of_indices);
        for (Py_ssize_t i = 0; i < PyList_Size(arrays); ++i)
        {
            PyObject* tuple_indices = PyList_GetItem(arrays, i); 
            if (!PyTuple_Check(tuple_indices))
            {
                // Si ce n'est pas un tuple, c'est un simple tableau :
                if (!PyArray_Check(tuple_indices))
                {
                    PyErr_SetString(PyExc_TypeError, 
                                    "Wrong syntax : each element of the list must be a tuple as (beg_indices, indices) or indices array");
                    return NULL;
                }
                np_beg_indices.push_back(nullptr);
                np_indices.push_back((PyArrayObject*)tuple_indices);
            }
            else { // OK, c'est un tuple...
                PyObject* arr1 = PyTuple_GetItem(tuple_indices, 0);
                PyObject* arr2 = PyTuple_GetItem(tuple_indices, 1);
                if ((!PyArray_Check(arr1)) || (!PyArray_Check(arr2)))
                {
                    PyErr_SetString(PyExc_TypeError, "Both elements of the tuple must be numpy arrays");
                    return NULL;
                }
                np_beg_indices.push_back((PyArrayObject*)arr1);
                np_indices.push_back((PyArrayObject*)arr2);
            }
        }
    }// Fin de if (PyList_Check(arrays))
    else if (PyTuple_Check(arrays))
    {
        PyObject* arr1 = PyTuple_GetItem(arrays, 0);
        PyObject* arr2 = PyTuple_GetItem(arrays, 1);
        if ((!PyArray_Check(arr1)) || (!PyArray_Check(arr2)))
        {
            PyErr_SetString(PyExc_TypeError, "Both elements of the tuple must be numpy arrays");
            return NULL;
        }
        np_beg_indices.push_back((PyArrayObject*)arr1);
        np_indices.push_back((PyArrayObject*)arr2);
    }
    else if (PyArray_Check(arrays))
    {
        np_beg_indices.push_back(nullptr);
        np_indices.push_back((PyArrayObject*)arrays);
    }
    else {
        PyErr_SetString(PyExc_TypeError, "First argument must be an array or a list of arrays or tuples with arrays");
        return NULL;
    }

    Py_ssize_t sz_list = np_indices.size();
    PyObject *compressed_list = PyList_New(sz_list);
    for ( E_Int i = 0; i < sz_list; ++i )
    {
        //= Dans tous les cas, on va mettre au même format, un tableau avec :
        //= [1 (un simple tableau à l'origine) ou 2 (beg + indices), nbre indices 1er élt/face, indice 1er indice sommet,
        //   position relative des autres sommets, nbre indices 2e élt/face, indice 1er indice sommet, ...]
        bool is_c_order = false;
        if (PyArray_CHKFLAGS(np_indices[i], NPY_ARRAY_C_CONTIGUOUS)) is_c_order = true;
        std::vector<E_Int> rel_indices;
        if (np_beg_indices[i] == nullptr)
        {   
            // On est dans le cas d'un format [nbre sommets, s1, s2, ..., nbre sommets, ...]
            Py_ssize_t sz_array = PyArray_SIZE(np_indices[i]);
            E_Int* data = (E_Int*)PyArray_DATA(np_indices[i]);
            // Transformation des data pour améliorer la compression
            int pteur = 0;
            std::vector<E_Int>(sz_array+1).swap(rel_indices);
            rel_indices[0] = 1;
            while (pteur < sz_array)
            {
                rel_indices[pteur+1] = data[pteur    ];// Le N
                assert(pteur+2 < sz_array);
                rel_indices[pteur+2] = data[pteur + 1];// Le 1er index
                for ( int j = 2; j <= data[pteur]; ++j )
                {
                    assert(pteur+j < sz_array);
                    rel_indices[pteur+j+1] = data[pteur+j] - data[pteur+1];
                }
                pteur += data[pteur]+1;
            }
        }
        else
        {// Cas où on est dans le format direct avec deux tableaux :
            Py_ssize_t sz_beg_ind = PyArray_SIZE(np_beg_indices[i]);
            E_Int* beg_indices    = (E_Int*)PyArray_DATA(np_beg_indices[i]);
            Py_ssize_t sz_ind     = PyArray_SIZE(np_indices[i]);
            E_Int* indices        = (E_Int*)PyArray_DATA(np_indices[i]);

            int pteur = 0;
            std::vector<E_Int>(sz_beg_ind+sz_ind).swap(rel_indices);
            rel_indices[0] = 2;
            E_Int* data_rel_ind = rel_indices.data() + 1;
            for ( E_Int ibeg = 0; ibeg < sz_beg_ind-1; ++ibeg)
            {
                E_Int nb_verts = beg_indices[ibeg+1] - beg_indices[ibeg];
                data_rel_ind[pteur++] = nb_verts;
                data_rel_ind[pteur++] = indices[beg_indices[ibeg]];
                for ( E_Int j = 1; j < nb_verts; ++j)
                {
                    assert(pteur < sz_beg_ind + sz_ind - 1);
                    data_rel_ind[pteur++] = indices[beg_indices[ibeg]+j] - indices[beg_indices[ibeg]];
                }
            }
        }
        //= On prépare l'objet décrivant la compression du NGon
        PyObject *obj = PyTuple_New(3);
        //= Rajout de la taille du tableau au tuple
        PyObject *py_length = PyLong_FromLong(rel_indices.size());
        PyTuple_SET_ITEM(obj, 0, py_length);
        //= Compression de rel_indices with zstd :
        std::size_t maxSize = ZSTD_compressBound(rel_indices.size()*sizeof(E_Int));
        std::vector<std::uint8_t> compressed_data(maxSize);
        std::size_t szCompress = ZSTD_compress(compressed_data.data(), maxSize, rel_indices.data(), 
                                               rel_indices.size()*sizeof(E_Int), ZSTD_maxCLevel());
        npy_intp sz = npy_intp(szCompress);
        PyArrayObject* cpr_arr = (PyArrayObject*)PyArray_SimpleNew(1, &sz, NPY_BYTE);
        std::uint8_t* buffer = (std::uint8_t*)PyArray_DATA(cpr_arr);
        std::copy(compressed_data.data(), compressed_data.data()+szCompress, buffer);
        //= On rajoute le tableau au tuple (nb_v_p_e,size,buffer)
        PyTuple_SET_ITEM(obj, 1, (PyObject*)cpr_arr);
        //$ On rajoute un flag pour savoir si on est en C order ou F order                          
        if (is_c_order)
        {
            Py_INCREF(Py_True);
            PyTuple_SetItem(obj, 2, Py_True);
        }
        else
        {
            Py_IncRef(Py_False);
            PyTuple_SetItem(obj, 2, Py_False);
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
//@ ____________________________________________________________________
PyObject *
py_ngon_indices_uncompress(PyObject *self, PyObject *args)
{
    PyObject *cpr_arrays;
    if (!PYPARSETUPLE_(args, O_, &cpr_arrays)) 
    {
        PyErr_SetString(PyExc_SyntaxError, "Wrong syntax. Right syntax : unpackIndices(array or list of compressed arrays");
        return NULL;
    }
    std::vector<PyArrayObject *> np_cpr_arrays;
    std::vector<npy_intp> length_arrays;
    std::vector<bool> is_c_order;
    if (PyList_Check(cpr_arrays)) {
        Py_ssize_t list_size = PyList_Size(cpr_arrays);
        np_cpr_arrays.reserve(list_size);
        length_arrays.reserve(list_size);
        is_c_order.reserve(list_size);
        for (Py_ssize_t i = 0; i < list_size; ++i) {
            PyObject *tuple = PyList_GetItem(cpr_arrays, i);
            if (!PyTuple_Check(tuple)) {
                PyErr_SetString(PyExc_TypeError, "The values of the list must be tuple as (length, compressed data)");
                return NULL;
            }
            PyObject *flag  = PyTuple_GetItem(tuple, 2);
            if (flag == Py_True)
                is_c_order.push_back(true);
            else
                is_c_order.push_back(false);
            PyObject *array = PyTuple_GetItem(tuple, 1);
            if (!PyArray_Check(array)) {
                PyErr_SetString(PyExc_TypeError, "Second value of the tuple must be an array with compressed data");
                return NULL;
            }
            np_cpr_arrays.push_back((PyArrayObject *)array);
            PyObject *lgth = PyTuple_GetItem(tuple, 0);
            if (!PyLong_Check(lgth)) {
                PyErr_SetString(PyExc_TypeError, "First element of the tuple must be an integer");
                return NULL;
            }
            length_arrays.push_back(npy_intp(PyLong_AsLong(lgth)));            
        } // End for (i= 0; i < list_size; ...)
    }     // Fin du cas si c'est une liste
    else if (PyTuple_Check(cpr_arrays)) {
        PyObject *lgth        = PyTuple_GetItem(cpr_arrays, 0);
        PyObject *array       = PyTuple_GetItem(cpr_arrays, 1);
        PyObject *flag        = PyTuple_GetItem(cpr_arrays, 2);
        if (!PyArray_Check(array)) {
            PyErr_SetString(PyExc_TypeError, "Second value of tuple must be an array with compressed data");
            return NULL;
        }
        np_cpr_arrays.push_back((PyArrayObject *)array);
        if (!PyLong_Check(lgth)) {
            PyErr_SetString(PyExc_TypeError, "First element of the tuple must be an integer");
            return NULL;
        }
        length_arrays.push_back(npy_intp(PyLong_AsLong(lgth)));
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
        E_Int kind_of_array = uncompressed_data[0];
        PyArrayObject *py_array, *py_index = nullptr;
        E_Int ordering = (is_c_order[i] == true ? 0 : 1);
        if (kind_of_array == 1)
        {
            npy_intp size_of_array = length-1;
            py_array = (PyArrayObject *)PyArray_EMPTY(1, &size_of_array, E_NPY_INT, ordering);
            E_Int *py_array_data = (E_Int *)PyArray_DATA(py_array);
            // Reconstitution des indices :
            size_t cpteur = 1;
            while (cpteur < uncompressed_data.size())
            {
                E_Int N = uncompressed_data[cpteur];
                assert(N>=1);
                py_array_data[cpteur-1] = N;
                py_array_data[cpteur+0] = uncompressed_data[cpteur+1];
                for ( E_Int j = 1; j < N; ++j)
                {
                    assert(cpteur+j+1 < uncompressed_data.size());
                    py_array_data[cpteur+j] = uncompressed_data[cpteur+j+1] + uncompressed_data[cpteur+1];
                }
                cpteur += N+1;
            }
        }
        else
        {
            //= On va compter le nombre de faces/éléments à traiter :
            size_t cpteur = 1; E_Int nb_elts = 0;
            while ( cpteur < uncompressed_data.size() )
            {
                nb_elts += 1;
                cpteur += uncompressed_data[cpteur];
            }
            //= Réservation (et initialisation) des tableaux numpy
            npy_intp sz_index    = nb_elts+1;
            py_index = (PyArrayObject*)PyArray_EMPTY(1, &sz_index, NPY_INT32, ordering);
            E_Int *py_index_data = (E_Int *)PyArray_DATA(py_index);
            npy_intp sz_py_array = length - nb_elts - 1;
            py_array = (PyArrayObject *)PyArray_EMPTY(1, &sz_py_array, NPY_INT32, ordering);
            E_Int *py_array_data = (E_Int *)PyArray_DATA(py_array);
            E_Int beg_elt = 0; cpteur = 1;
            for (E_Int ielt = 0; ielt < nb_elts; ++ielt )
            {
                py_index_data[ielt] = beg_elt;
                E_Int N = uncompressed_data[cpteur];
                assert(N>=1);
                py_array_data[beg_elt+0] = uncompressed_data[cpteur+1];
                for ( int j = 1; j < N; ++j )
                {
                    assert(cpteur+j+1 < uncompressed_data.size());
                    py_array_data[beg_elt+j] = uncompressed_data[cpteur+j+1] + uncompressed_data[cpteur+1];
                }
                beg_elt += N;
                cpteur  += N+1;
            }
        }// Fin des deux types de représentation des indices NGons
        if (py_index == nullptr)
            PyList_SetItem(lst_out_arrays, i, (PyObject *)py_array);
        else
        {
            PyObject* tuple = PyTuple_New(2);
            PyTuple_SetItem(tuple, 0, (PyObject *)py_index);
            PyTuple_SetItem(tuple, 1, (PyObject *)py_array);
        }
    }// End for ( E_int i = 0; i < np_cpr_arrays.size(); ++i )
    return lst_out_arrays;
}// End uncompress_NGon_Indices function
}
