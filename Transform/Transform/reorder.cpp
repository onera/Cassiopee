/*    
    Copyright 2013-2019 Onera.

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

using namespace std;
using namespace K_FLD;

// ============================================================================
/* Reorder the numerotation of mesh */
// ============================================================================
PyObject* K_TRANSFORM::reorder(PyObject* self, PyObject* args)
{
  E_Int oi=1, oj=1, ok=1;
  PyObject* array; PyObject* order;
  if (!PyArg_ParseTuple(args, "OO",
                        &array, &order))
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
  PyObject* tpl;
  if (size == 3)
  {
    tpl = PyTuple_GetItem(order, 0);
    oi = PyLong_AsLong(tpl);
    tpl = PyTuple_GetItem(order, 1);
    oj = PyLong_AsLong(tpl);
    tpl = PyTuple_GetItem(order, 2);
    ok = PyLong_AsLong(tpl);
  }
  else if (size == 1)
  {
    tpl = PyTuple_GetItem(order, 0);
    oi = PyLong_AsLong(tpl);
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, 
                    "reorder: order must be like (1,2,3) or (1,).");
    return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType,
                          true); 
  
  E_Int nfld = f->getNfld();
  E_Int npts = f->getSize();

  if (res == 1)
  {
    PyObject* tpl = K_ARRAY::buildArray(nfld, varString, im, jm, km);
    E_Float* foutp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF fout(npts, nfld, foutp, true);
    K_CONNECT::reorderStructField(im, jm, km, *f, fout, 
                                  E_Int(oi), E_Int(oj), E_Int(ok));
    PyList_SetItem(tpl,2, PyInt_FromLong(im));
    PyList_SetItem(tpl,3, PyInt_FromLong(jm));
    PyList_SetItem(tpl,4, PyInt_FromLong(km));
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else if (res == 2)
  {
    if (strcmp(eltType, "QUAD") == 0 || strcmp(eltType, "TRI") == 0 ||
        strcmp(eltType, "QUAD*") == 0 || strcmp(eltType, "TRI*") == 0)
    {
      PyObject* tpl = K_ARRAY::buildArray(nfld, varString,
                                          npts, cn->getSize(),
                                          -1, eltType);
      E_Float* fp = K_ARRAY::getFieldPtr(tpl);
      FldArrayF fn(npts, nfld, fp, true); fn = *f;
      E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
      FldArrayI cnn(cn->getSize(), cn->getNfld(), cnnp, true); cnn = *cn;
      K_CONNECT::reorderQuadTriField(fn, cnn, E_Int(oi));
      RELEASESHAREDU(array, f, cn);
      return tpl;
    }
    else if (strcmp(eltType, "NGON") == 0) 
    {
      // check si le NGON est surfacique
      E_Int* cnp = cn->begin();

      if (cnp[2] > 2) // la face a plus de 2 sommets ce n'est pas une arete
      {
        PyErr_SetString(PyExc_TypeError,
                        "reorderNGON: NGON array must be a surface.");
        RELEASESHAREDU(array, f, cn); return NULL;
      }
      if (cnp[2] == 1) // la face a 1 seul sommet
      {
        tpl = K_ARRAY::buildArray(*f, varString, *cn, -1, eltType);
        RELEASESHAREDU(array, f, cn);
        return tpl;
      }
      
      E_Int csize = cn->getSize()*cn->getNfld();
      PyObject* tpl = K_ARRAY::buildArray(nfld, varString, npts, cn->getSize(),
                                          -1, eltType, false, csize);
      E_Float* fp = K_ARRAY::getFieldPtr(tpl);
      FldArrayF fn(npts, nfld, fp, true); fn = *f;
      E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
      FldArrayI cnn(cn->getSize(), cn->getNfld(), cnnp, true); cnn = *cn;
      K_CONNECT::reorderNGON(fn, cnn, E_Int(oi));
      RELEASESHAREDU(array, f, cn);
      return tpl;
    }
    else
    {
      RELEASESHAREDB(res, array, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "reorder: only for TRI-array or QUAD-array.");
      return NULL;
    }
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "reorder: unknown type of array.");
    return NULL;
  }
}
