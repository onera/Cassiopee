/*    
    Copyright 2013-2023 Onera.

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
#include "Array/Array.h"
#include "String/kstring.h"
#include <stdio.h>
#include <string.h>

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Build an empty structured array3
   IN: nfld: nbre de champs
   IN: varString: variables string
   IN: ni,nj,nk: number of points in field
   IN: api (1: array, 2: array2, 3: array3)
   OUT: PyObject created. */
//=============================================================================
PyObject* K_ARRAY::buildArray3(E_Int nfld, const char* varString, 
                               E_Int ni, E_Int nj, E_Int nk, E_Int api)
{
    PyObject* tpl;
    IMPORTNUMPY;

    if (api == 1) // Array1
    {
        npy_intp dim[2];
        dim[1] = ni*nj*nk;
        dim[0] = nfld;
        PyArrayObject* a = (PyArrayObject*)PyArray_SimpleNew(2, dim, NPY_DOUBLE);
        tpl = Py_BuildValue("[sOlll]", varString, a, (long)ni, (long)nj, (long)nk);
        Py_DECREF(a);
    }
    else // Array2 ou Array3
    {
        npy_intp dim[3]; int ndim=3;
        dim[0] = ni; dim[1] = nj; dim[2] = nk;
        if (nk == 1) ndim--;
        if (nj == 1) ndim--;
        PyObject* rake = PyList_New(0);
        for (E_Int n=0; n < nfld; n++)
        {
            PyArrayObject* a = (PyArrayObject*)PyArray_EMPTY(3, dim, NPY_DOUBLE, 1);
            PyList_Append(rake, (PyObject*)a); Py_DECREF(a);
        }
        tpl = Py_BuildValue("[sOlll]", varString, rake, (long)ni, (long)nj, (long)nk);
        Py_DECREF(rake);
    }
    return tpl;
}

//=============================================================================
/* Build an empty NGON array 
   IN: nfld: number of fields
   IN: varString: variable string
   IN: nvertex: number of vertex

   IN: nelt: number total of elements
   IN: etString: NGON ou NGON*
   IN: center: set to true if field is localised in the centers of
   elements, otherwise let it to false.
   IN: sizeNGon, sizeNFace, nface: connectivity size.
   if sizeNFace == -1, NFACE is not created 

   OUT: PyObject created. */
//=============================================================================
// build pour les NGONS
PyObject* K_ARRAY::buildArray3(E_Int nfld, const char* varString,
                               E_Int nvertex, E_Int nelt, E_Int nface,
                               const char* etString,
                               E_Int sizeNGon, E_Int sizeNFace, 
                               E_Boolean center, E_Int api)
{
    npy_intp dim[2];
    PyObject* a; PyObject* ac; PyObject* tpl;
    char eltType[12];

    // taille de f
    E_Int fSize;
    if (center == true) fSize = nelt;
    else fSize = nvertex;

    IMPORTNUMPY;

    // element string - ajoute * pour les centres
    strcpy(eltType, etString);
    E_Int pos = strlen(eltType)-1;
    pos = 0;
    if (eltType[pos] != '*' && center == true) strcat(eltType, "*");
    else if (eltType[pos] == '*') eltType[pos] = '\0';

    // Build array of fields
    if (api == 1) // Array1
    { 
        dim[1] = fSize; dim[0] = nfld;
        a = PyArray_SimpleNew(2, dim, NPY_DOUBLE);
    }
    else // Array2 or Array3
    {
        dim[0] = fSize;
        a = PyList_New(0);
        for (E_Int n=0; n < nfld; n++)
        {
            PyArrayObject* ar = (PyArrayObject*)PyArray_EMPTY(1, dim, NPY_DOUBLE, 1);
            PyList_Append(a, (PyObject*)ar); Py_DECREF(ar);
        }
    } 

    // Build array for connectivity
    if (api == 1) // Array 1
    {
        dim[1] = 4+sizeNGon+sizeNFace; dim[0] = 1;
        ac = PyArray_SimpleNew(2, dim, E_NPY_INT);
        E_Int* data = (E_Int*)PyArray_DATA((PyArrayObject*)ac);
        data[0] = nface;
        data[1] = sizeNGon;
        data[sizeNGon+2] = nelt;
        data[sizeNGon+3] = sizeNFace;
    }
    else if (api == 2) // Array2
    {
        ac = PyList_New(0);
        // ngons - NGON - sizeNGon
        dim[0] = sizeNGon;
        PyObject* ar = PyArray_EMPTY(1, dim, E_NPY_INT, 0);
        PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
        // ngons - NFACE - sizeNFace
        dim[0] = sizeNFace;
        ar = PyArray_EMPTY(1, dim, E_NPY_INT, 0);
        PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
        // ngons - indPG - nfaces
        dim[0] = nface;
        ar = PyArray_EMPTY(1, dim, E_NPY_INT, 0);
        PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
        // ngons - indPH - nelts
        dim[0] = nelt;
        ar = PyArray_EMPTY(1, dim, E_NPY_INT, 0);
        PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
        // Eventuellement PE - 2*nface
        //dim[0] = nface; dim[1] = 2;
        //PyObject* ar = PyArray_EMPTY(2, dim, E_NPY_INT, 0);
        //PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
    }
    else
    {
        ac = PyList_New(0);
        // NGON - sizeNGon
        dim[0] = sizeNGon;
        PyObject* ar = PyArray_EMPTY(1, dim, E_NPY_INT, 0);
        PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
        // NFACE - sizeNFace
        dim[0] = sizeNFace;
        ar = PyArray_EMPTY(1, dim, E_NPY_INT, 0);
        PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
        // NGON - StartOffset
        dim[0] = nface+1;
        ar = PyArray_EMPTY(1, dim, E_NPY_INT, 0);
        PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
        // NFACE - startOffset
        dim[0] = nelt+1;
        ar = PyArray_EMPTY(1, dim, E_NPY_INT, 0);
        PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
        // Eventuellement PE - 2*nface
        //dim[0] = nface; dim[1] = 2;
        //PyObject* ar = PyArray_EMPTY(2, dim, E_NPY_INT, 0);
        //PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
    }
  
    tpl = Py_BuildValue("[sOOs]", varString, a, ac, eltType);
    Py_DECREF(a); Py_DECREF(ac);

    return tpl;
}


//=============================================================================
/* Build an empty BE array 
   IN: nfld: number of fields
   IN: varString: variable string
   IN: nvertex: number of vertex
   IN: nelt: number of elements
   IN: etString: "TRI" ou avec *
   IN: center: set to true if field is localised in the centers of
   elements, otherwise let it to false.

   OUT: PyObject created. */
//=============================================================================
// build pour les single Element (BE)
PyObject* K_ARRAY::buildArray3(E_Int nfld, const char* varString,
                               E_Int nvertex,
                               E_Int nelts,
                               const char* etString,
                               E_Boolean center, E_Int api)
{
    npy_intp dim[2];
    PyObject* a; PyObject* ac; PyObject* tpl;
    char eltType[256];
    strcpy(eltType, etString);

    // taille de f
    E_Int fSize;
    if (center == true) fSize = nelts;
    else fSize = nvertex;

    IMPORTNUMPY;

    // Build array of fields
    if (api == 1) // Array1
    { 
        dim[1] = fSize; dim[0] = nfld;
        a = PyArray_SimpleNew(2, dim, NPY_DOUBLE);
    }
    else // Array2 or Array3
    {
        dim[0] = fSize;
        a = PyList_New(0);
        for (E_Int n=0; n < nfld; n++)
        {
            PyArrayObject* ar = (PyArrayObject*)PyArray_EMPTY(1, dim, NPY_DOUBLE, 1);
            PyList_Append(a, (PyObject*)ar); Py_DECREF(ar);
        }
    } 

    // Connectivite
    if (api == 1) // Array1
    {
        E_Int cSize = nelts;
        char st[256]; E_Int dummy; E_Int nvpe;
        eltString2TypeId(eltType, st, nvpe, dummy, dummy);
        dim[1] = cSize; dim[0] = nvpe;
        ac = PyArray_SimpleNew(2, dim, E_NPY_INT);
    }
    else if (api == 2 || api == 3) // Array2 ou 3
    {
        E_Int cSize = nelts;
        char st[256]; E_Int dummy; E_Int nvpe;
        eltString2TypeId(eltType, st, nvpe, dummy, dummy);
        ac = PyList_New(0);
        dim[0] = cSize; dim[1] = nvpe;
        PyObject* ar = PyArray_EMPTY(2, dim, E_NPY_INT, 0);
        PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "buildArray: unkown api.");
        return NULL;
    }

    tpl = Py_BuildValue("[sOOs]", varString, a, ac, eltType);
    Py_DECREF(a); Py_DECREF(ac);

    return tpl;
}

//=============================================================================
/* Build an empty ME array 
   IN: nfld: number of fields
   IN: varString: variable string
   IN: nvertex: number of vertex
   IN: neltsPerType: number of elements per type
   IN: etString: "TRI,QUAD" ou avec *
   IN: center: set to true if field is localised in the centers of
   elements, otherwise let it to false.

   OUT: PyObject created. */
//=============================================================================
// build pour les Multiple Element (ME)
PyObject* K_ARRAY::buildArray3(E_Int nfld, const char* varString,
                               E_Int nvertex,
                               std::vector<E_Int>& neltsPerType,
                               const char* etString,
                               E_Boolean center, E_Int api)
{
    npy_intp dim[2];
    PyObject* a; PyObject* ac; PyObject* tpl;
    char eltType[256];
    strcpy(eltType, etString);

    // taille de f
    E_Int nelt = 0;
    for (size_t i = 0; i < neltsPerType.size(); i++) nelt += neltsPerType[i];
    E_Int fSize;
    if (center == true) fSize = nelt;
    else fSize = nvertex;

    IMPORTNUMPY;

    // Build array of fields
    if (api == 1) // Array1
    { 
        dim[1] = fSize; dim[0] = nfld;
        a = PyArray_SimpleNew(2, dim, NPY_DOUBLE);
    }
    else // Array2 or Array3
    {
        dim[0] = fSize;
        a = PyList_New(0);
        for (E_Int n=0; n < nfld; n++)
        {
            PyArrayObject* ar = (PyArrayObject*)PyArray_EMPTY(1, dim, NPY_DOUBLE, 1);
            PyList_Append(a, (PyObject*)ar); Py_DECREF(ar);
        }
    } 

    // Connectivite
    if (api == 1) // Array1
    {
        E_Int cSize = nelt;
        char st[256]; E_Int dummy; E_Int nvpe;
        eltString2TypeId(eltType, st, nvpe, dummy, dummy);
        dim[1] = cSize; dim[0] = nvpe;
        ac = PyArray_SimpleNew(2, dim, E_NPY_INT);
    }
    else if (api == 2) // Array2
    {
        E_Int cSize = nelt;
        char st[256]; E_Int dummy; E_Int nvpe;
        eltString2TypeId(eltType, st, nvpe, dummy, dummy);
        ac = PyList_New(0);
        dim[0] = cSize; dim[1] = nvpe;
        PyObject* ar = PyArray_EMPTY(2, dim, E_NPY_INT, 0);
        PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
    }
    else // Array3
    {
        std::vector<char*> eltTypes;
        K_ARRAY::extractVars(eltType, eltTypes);
        char st[256]; E_Int dummy; E_Int nvpe;
        ac = PyList_New(0);
        //printf("size=%d %s\n", eltTypes.size(), eltType);
        for (size_t i = 0; i < eltTypes.size(); i++)
        {
            E_Int cSize = neltsPerType[i];
            eltString2TypeId(eltTypes[i], st, nvpe, dummy, dummy);
            dim[0] = cSize; dim[1] = nvpe;
            PyObject* ar = PyArray_EMPTY(2, dim, E_NPY_INT, 0);
            PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
        }
        for (size_t i = 0; i < eltTypes.size(); i++) delete [] eltTypes[i];
    }

    tpl = Py_BuildValue("[sOOs]", varString, a, ac, eltType);
    Py_DECREF(a); Py_DECREF(ac);

    return tpl;
}
