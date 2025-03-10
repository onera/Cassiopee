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

// Routines pour extraire les donnees d'un objet python BCFaces

# include "GenIO.h"
# include <stdio.h>
# include <string.h>
# include "Array/Array.h"

using namespace K_FLD;

//=============================================================================
// Retourne le nombre de BCs stockees dans l'objet python BCFaces
// -1: BCFaces n'est pas valide
//=============================================================================
E_Int K_IO::GenIO::getSizeOfBCFaces(PyObject* BCFaces)
{
  if (PyList_Check(BCFaces) == false) return -1;
  E_Int l = PyList_Size(BCFaces);
  return l/2;
}

//=============================================================================
// Retourne le nom de la BC no i dans BCFaces,
// une copie des indices des faces de la BC (commencant a 1)
// a deleter ensuite.
// IN: BCFaces: objet python
// IN: i: no de la BC voulue
// OUT: name: nom de la BC, doit etre deja alloue
// OUT: faces: indices des faces de la BC i
// Retourne 1: OK
// Retourne -1: FAILED
//=============================================================================
E_Int K_IO::GenIO::getBCFaces(PyObject* BCFaces, E_Int i, char* name,
                              FldArrayI& faces)
{
  // BC name
  PyObject* o = PyList_GetItem(BCFaces, 2*i);
  if (PyString_Check(o)) strcpy(name, PyString_AsString(o));
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(o)) strcpy(name, PyBytes_AsString(PyUnicode_AsUTF8String(o))); 
#endif
  else return -1;
  
  // BC data
  o = PyList_GetItem(BCFaces, 2*i+1);
  PyArrayObject* a;
  //a = (PyArrayObject*)PyArray_ContiguousFromObject(o, NPY_INT,
  //                                                 1, 10000000);
  a = (PyArrayObject*)o; Py_INCREF(a);
  if (a == NULL) return -1;
  
  E_Int s;
  if (PyArray_NDIM(a) == 2) s = PyArray_DIMS(a)[1]*PyArray_DIMS(a)[0];
  else if (PyArray_NDIM(a) == 1) s = PyArray_DIMS(a)[0];
  else return -1;
  
  faces.malloc(s);
  E_Int* facesp = faces.begin();
  E_Int* ptr = (E_Int*)PyArray_DATA(a);
  for (E_Int i = 0; i < s; i++) facesp[i] = ptr[i];
  Py_DECREF(a);
  return 1;
}
