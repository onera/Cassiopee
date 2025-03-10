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
#include "PyTree/PyTree.h"
#include "Numpy/Numpy.h"

#define CREATENODE  \
PyObject* p1 = PyString_FromString(name);       \
PyObject* p3 = PyList_New(0);                   \
PyObject* p4 = PyString_FromString(type);       \
PyObject* p = PyList_New(0);                    \
PyList_Append(p, p1); PyList_Append(p, p2);     \
PyList_Append(p, p3); PyList_Append(p, p4);     \
PyObject* children = PyList_GetItem(o, 2);      \
if (pos == -1) PyList_Append(children, p);      \
else PyList_Insert(children, pos, p);

//=============================================================================
// Ajoute un noeud a o.
// IN: o: objet representant un noeud de pyTree
// IN: name: le nom du noeud a creer
// IN: type: le type du noeud a creer
// IN: value: la valeur a lui donner
// OUT: modifie o
//=============================================================================
PyObject* K_PYTREE::createChild(PyObject* o, const char* name, 
                                const char* type, 
                                E_Float* value, E_Int size, E_Int nfld,
                                E_Int pos)
{
  PyObject* p2 = K_NUMPY::buildNumpyArray(value, size, nfld, 1);
  CREATENODE;
  return p;
}
PyObject* K_PYTREE::createChild(PyObject* o, const char* name, 
                                const char* type, 
                                PyObject* value,
                                E_Int pos)
{
  PyObject* p2 = value;
  if (value == NULL) { Py_INCREF(Py_None); p2 = Py_None; }
  CREATENODE;
  return p;
}
PyObject* K_PYTREE::createChild(PyObject* o, const char* name, 
                                const char* type, 
                                E_Int* value, E_Int size, E_Int nfld,
                                E_Int pos)
{
  PyObject* p2 = K_NUMPY::buildNumpyArray(value, size, nfld, 1);
  CREATENODE;
  return p;
}
