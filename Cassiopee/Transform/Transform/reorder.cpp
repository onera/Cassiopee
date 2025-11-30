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

#include "transform.h"

using namespace K_FLD;

// ============================================================================
// Reorder the vertex indices of the mesh elements
// ============================================================================
PyObject* K_TRANSFORM::reorder(PyObject* self, PyObject* args)
{
  E_Int oi = 1, oj = 1, ok = 1;
  PyObject* array; PyObject* order;
  if (!PYPARSETUPLE_(args, OO_, &array, &order))
  {
    return NULL;
  }

  // Check tuple
  if (PyTuple_Check(order) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "reorder: order argument must be a tuple.");
    return NULL;
  }
  E_Int size = PyTuple_Size(order);
  PyObject *tpl, *tplo;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  E_Int api = f->getApi();
  FldArrayF* f2;

  if (res == 1)
  {
    if (size == 3)
    {
      tplo = PyTuple_GetItem(order, 0); oi = PyLong_AsLong(tplo);
      tplo = PyTuple_GetItem(order, 1); oj = PyLong_AsLong(tplo);
      tplo = PyTuple_GetItem(order, 2); ok = PyLong_AsLong(tplo);
    }
    else
    {
      PyErr_SetString(PyExc_TypeError,
                      "reorder: order tuple must be of size 3, eg, (1,2,3).");
      RELEASESHAREDS(array, f);
      return NULL;
    }
  
    tpl = K_ARRAY::buildArray3(*f, varString, im, jm, km, api);
    K_ARRAY::getFromArray3(tpl, f2);
    K_CONNECT::reorderStructField(im, jm, km, *f, *f2,
                                  E_Int(oi), E_Int(oj), E_Int(ok));
    PyList_SetItem(tpl, 2, PyInt_FromLong(im));
    PyList_SetItem(tpl, 3, PyInt_FromLong(jm));
    PyList_SetItem(tpl, 4, PyInt_FromLong(km));
    RELEASESHAREDS(tpl, f2);
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else if (res == 2)
  {
    FldArrayI* cn2;
    if (strncmp(eltType, "NGON", 4) == 0)
    {
      if (size == 1)
      {
        tplo = PyTuple_GetItem(order, 0); oi = PyLong_AsLong(tplo);
      }
      else
      {
        PyErr_SetString(PyExc_TypeError,
                        "reorder: order tuple must be (1,) or (-1,) for NGONs.");
        RELEASESHAREDU(array, f, cn);
        return NULL;
      }
  
      // check si le NGON est surfacique
      E_Int dim = cn->getDim();

      if (dim == 3)
      {
        PyErr_SetString(PyExc_TypeError,
                        "reorderNGON: NGON array must be a surface.");
        RELEASESHAREDU(array, f, cn);
        return NULL;
      }
      else if (dim == 1)
      {
        tpl = K_ARRAY::buildArray3(*f, varString, *cn, eltType, api);
        RELEASESHAREDU(array, f, cn);
        return tpl;
      }

      tpl = K_ARRAY::buildArray3(*f, varString, *cn, eltType, api);
      K_ARRAY::getFromArray3(tpl, f2, cn2);
      K_CONNECT::reorderNGON(*f2, *cn2, E_Int(oi));
      RELEASESHAREDU(tpl, f2, cn2);
      RELEASESHAREDU(array, f, cn);
      return tpl;
    }
    else  // BE/ME
    {
      E_Int dim = K_CONNECT::getDimME(eltType);
      tpl = K_ARRAY::buildArray3(*f, varString, *cn, eltType, api);
      K_ARRAY::getFromArray3(tpl, f2, cn2);
      if (dim == 2)
      {
        if (size == 1)
        {
          tplo = PyTuple_GetItem(order, 0); oi = PyLong_AsLong(tplo);
        }
        else
        {
          PyErr_SetString(PyExc_TypeError,
                          "reorder: order tuple must be (1,) or (-1,) for TRI or QUAD.");
          RELEASESHAREDU(tpl, f2, cn2);
          RELEASESHAREDU(array, f, cn);
          return NULL;
        }
        K_CONNECT::reorderUnstruct2D(*f2, *cn2, E_Int(oi));
      }
      else if (dim == 3)
      {
        K_CONNECT::reorderUnstruct3D(varString, *f2, *cn2, eltType);
      }
      RELEASESHAREDU(tpl, f2, cn2);
      RELEASESHAREDU(array, f, cn);
      return tpl;
    }
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "reorder: unknown type of array.");
    return NULL;
  }
}
