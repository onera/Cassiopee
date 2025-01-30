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
#include "geom.h"

using namespace K_FLD;

//===========================================================================
/* Polyline */
//===========================================================================
PyObject* K_GEOM::polyline( PyObject* self, PyObject* args )
{
  PyObject* listPts;

  /* 1: verifications */
  if (!PyArg_ParseTuple(args, "O", &listPts)) return NULL;

  // verifie que l'objet est une liste de tuples
  if (PyList_Check(listPts) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "polyline: argument must be a list of tuples.");
    return NULL;
  }
  
  /* Taille de la liste de listes et conversion en C */
  E_Int npts = PyList_Size(listPts); // nb d'elements ds la liste
  
  for (E_Int i = 0; i < npts; i++)
  {
    PyObject* tpli = PyList_GetItem(listPts,i);
    if (PyTuple_Check(tpli) == 0 && PyList_Check(tpli) == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "polyline: each element of the list must be (x,y,z).");
      
      return NULL;
    }
    
    E_Int dim;  
    if (PyTuple_Check(tpli) != 0)
    {
      dim = PyTuple_Size(tpli); 
    }
    else
    {
      dim = PyList_Size(tpli); 
    }
    
    // verifie que chq element de la liste est un triplet (x,y,z)
    if (dim != 3 && dim != 2)
    {        
      PyErr_SetString(PyExc_TypeError,
                      "polyline: 2 or 3 coordinates must be found in each element of the list.");
      return NULL;  
    }
  }
 
  /* 2: lecture des coordonnees et creation du tableau de pts */
  PyObject* tpl = K_ARRAY::buildArray(3, "x,y,z", npts, 1, 1);

  E_Float* x = K_ARRAY::getFieldPtr(tpl);
  E_Float* y = x + npts;
  E_Float* z = y + npts;
  
  PyObject *o, *c1, *c2, *c3;
 
  for (E_Int j = 0; j < npts; j++)
  {
    o = PyList_GetItem(listPts, j); // on recupere les listes des elements
    if (PyTuple_Check(o) == true && PyTuple_Size(o) == 2)
    {
      c1 = PyTuple_GetItem(o, 0);
      c2 = PyTuple_GetItem(o, 1);
      x[j] = PyFloat_AsDouble(c1);
      y[j] = PyFloat_AsDouble(c2);
      z[j] = 0.;
    }
    else if (PyList_Check(o) == true && PyList_Size(o) == 2)
    {
      c1 = PyList_GetItem(o, 0);
      c2 = PyList_GetItem(o, 1);
      x[j] = PyFloat_AsDouble(c1);
      y[j] = PyFloat_AsDouble(c2);
      z[j] = 0.;
    }
    else if (PyTuple_Check(o) == true && PyTuple_Size(o) == 3)
    {
      c1 = PyTuple_GetItem(o, 0);
      c2 = PyTuple_GetItem(o, 1);
      c3 = PyTuple_GetItem(o, 2);
      x[j] = PyFloat_AsDouble(c1);
      y[j] = PyFloat_AsDouble(c2);
      z[j] = PyFloat_AsDouble(c3);
    }
    else if (PyList_Check(o) == true && PyList_Size(o) == 3)
    {
      c1 = PyList_GetItem(o, 0);
      c2 = PyList_GetItem(o, 1);
      c3 = PyList_GetItem(o, 2);
      x[j] = PyFloat_AsDouble(c1);
      y[j] = PyFloat_AsDouble(c2);
      z[j] = PyFloat_AsDouble(c3);
    }
    else
    {
      x[j] = 0.;
      y[j] = 0.;
      z[j] = 0.;
    }
  }
  
  return tpl;
}
