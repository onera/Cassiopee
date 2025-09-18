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

# include "transform.h"
# include "Nuga/include/TSSplitter.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* 
   splitTRI: 
   Splits a TRI into several TRIs delimited by the input polyline.
*/
//=============================================================================
PyObject* K_TRANSFORM::splitTRI(PyObject* self, PyObject* args)
{
  PyObject *array, *polyLineList;
  if (!PYPARSETUPLE_(args, OO_, &array, &polyLineList))
  {
    return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FloatArray* f; IntArray* cn;
  char* varString; char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType); 

  if (res == 1)
  {
    delete f;
    PyErr_SetString(PyExc_TypeError,
                    "splitTRI: can not be used on a structured array.");
    return NULL;
  }
  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitTRI: unknown type of array.");
    return NULL;
  }
  
  if (strcmp(eltType, "TRI") != 0)
  {
    delete f; delete cn;
    PyErr_SetString(PyExc_TypeError,
                    "splitTRI: must be used on a TRI-array.");
    return NULL;
  }
  
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
   
  if (posx == -1 || posy == -1)
  {
    delete f; delete cn;
    PyErr_SetString(PyExc_TypeError,
                    "splitTRI: can't find coordinates in array.");
    return NULL;
  }

  // Check and get polyLines.
  PyObject *polyLine, *pN0, *pN1;
  IntArray connectP; 
  if (PyList_Check(polyLineList) == 0)
  {
    PyErr_Warn(PyExc_Warning,
               "splitTRI: an array must be a list of list of indices. Check list.");
    return NULL;
  }

  E_Int lsize = PyList_Size(polyLineList), sz, Edge[2];
  for (E_Int i = 0; i < lsize; ++i)
  {
    polyLine = PyList_GetItem(polyLineList, i);
    
    if (PyList_Check(polyLineList) == 0)
    {
      PyErr_Warn(PyExc_Warning,
               "splitTRI: an array must be a list of list of indices. Check list.");
      return NULL;
    }

    sz = PyList_Size(polyLine);

    for (E_Int j = 0; j < sz-1; ++j)
    {
      pN0 = PyList_GetItem(polyLine, j);
      pN1 = PyList_GetItem(polyLine, j+1);
      
      if (PyLong_Check(pN0) == 0 && PyInt_Check(pN0) == 0)
      {
        PyErr_Warn(PyExc_Warning,
                 "splitTRI: indices must be integers.");
        return NULL;
      }

      if (PyLong_Check(pN1) == 0 && PyInt_Check(pN1) == 0)
      {
        PyErr_Warn(PyExc_Warning,
                 "splitTRI: indices must be integers.");
        return NULL;
      }
      Edge[0] = PyLong_AsLong(pN0);
      Edge[1] = PyLong_AsLong(pN1);
      connectP.pushBack(Edge, Edge+2);
    }
  }
  
  vector<IntArray>   cs;
  vector<FloatArray> fs;
  TSSplitter::split(*f, *cn, connectP, fs, cs);
  delete f; delete cn;

  // Formation des array de sortie
  PyObject* tpl;
  PyObject* l = PyList_New(0);
  
  for (size_t i = 0; i < cs.size(); ++i)
  {
    tpl = K_ARRAY::buildArray(fs[i], varString, cs[i], -1, eltType);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  return l;
}
