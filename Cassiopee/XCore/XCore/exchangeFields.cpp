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
#include "xcore.h"
#include <mpi.h>

PyObject *K_XCORE::exchangeFields(PyObject *self, PyObject *args)
{
  PyObject *arr, *pe, *flds, *comm_list;
  if (!PYPARSETUPLE_(args, OOOO_, &arr, &pe, &flds, &comm_list)) 
  {
    PyErr_SetString(PyExc_ValueError, "Bad input.");
    return NULL;
  }

  // Process comm data
  E_Int psize = PyList_Size(comm_list);
  if (psize == 0) { Py_INCREF(Py_None); return Py_None; }

  // Check array
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  char* varString; char* eltType;
  E_Int ret = K_ARRAY::getFromArray3(arr, varString, f, ni, nj, nk, cn, eltType);

  if (ret <= 0) {
    PyErr_SetString(PyExc_TypeError, "Bad mesh array.");
    return NULL;
  }

  if (ret == 1) {
    PyErr_SetString(PyExc_TypeError, "Only for NGons.");
    RELEASESHAREDS(arr, f);
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(ret, arr, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "Can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  // Process solution fields
  E_Int fsize = PyList_Size(flds);
  if (fsize == 0) {
    RELEASESHAREDU(arr, f, cn);
    Py_INCREF(Py_None);
    return Py_None;
  }

  E_Float **fields = (E_Float **)calloc(fsize, sizeof(E_Float *));
  for (E_Int i = 0; i < fsize; i++) {
    PyObject *fld = PyList_GetItem(flds, i);
    E_Int nc = -1;
    ret = K_NUMPY::getFromNumpyArray(fld, fields[i], nc);
    assert(nc == cn->getNElts());
  }

  // Parent Elements
  E_Int *PE;
  E_Int pe_size;
  K_NUMPY::getFromNumpyArray(pe, PE, pe_size);
  assert(pe_size == 2*cn->getNFaces());

  // Comm data
  E_Int *procs = (E_Int *)malloc(psize * sizeof(E_Int));
  E_Int **ptlists = (E_Int **)calloc(psize, sizeof(E_Int *));
  E_Int *npfaces = (E_Int *)malloc(psize * sizeof(E_Int));
  for (E_Int i = 0; i < psize; i++) 
  {
    PyObject *proc_and_list = PyList_GetItem(comm_list, i);
    PyObject *proc = PyList_GetItem(proc_and_list, 0);
    PyObject *list = PyList_GetItem(proc_and_list, 1);
    procs[i] = PyLong_AsLong(proc);
    ret = K_NUMPY::getFromNumpyArray(list, ptlists[i], npfaces[i]);
  }

  // Exchange fields one by one
  E_Float *send_buf = NULL;
  E_Int *owner = PE;
  E_Int req_size = fsize*psize*2;
  MPI_Request *reqs = (MPI_Request *)malloc(req_size * sizeof(MPI_Request));
  E_Int nreq = 0;

  PyObject *out = PyList_New(0);
  
  for (E_Int i = 0; i < psize; i++) {
    E_Int npf = npfaces[i];
    E_Int *pfaces = ptlists[i];
    E_Int dest = procs[i];
    npy_intp dims[2];
    dims[0] = npf;
    dims[1] = 1;

    PyObject *proc_rdata = PyList_New(0);
    
    send_buf = (E_Float *)realloc(send_buf, npf * sizeof(E_Float));

    for (E_Int j = 0; j < fsize; j++) {
      const auto &data = fields[j];
      PyArrayObject *recv = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
      
      for (E_Int k = 0; k < npf; k++) {
        E_Int face = pfaces[k]-1;
        E_Int own = owner[face]-1;
        E_Float val = data[own];
        send_buf[k] = val;
      }

      MPI_Irecv(PyArray_DATA(recv), npf, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD,
        &reqs[nreq++]); 
      MPI_Isend((void *)send_buf, npf, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD,
        &reqs[nreq++]);
      
      PyList_Append(proc_rdata, (PyObject *)recv);
      Py_DECREF(recv);
    }

    PyList_Append(out, proc_rdata);
    Py_DECREF(proc_rdata);
  }

  // WARNING(Imad): Where do we put the wait_all?
  MPI_Waitall(nreq, reqs, MPI_STATUSES_IGNORE);
  nreq = 0;

  // FREE MPI STUFF
  free(reqs);
  free(send_buf);

  // FREE COMM DATA
  free(procs);
  for (E_Int i = 0; i < psize; i++) {
    PyObject *o = PyList_GetItem(comm_list, i);
    Py_DECREF(PyList_GetItem(o, 1));
  }
  free(ptlists);
  free(npfaces);

  // FREE FIELDS
  for (E_Int i = 0; i < fsize; i++) Py_DECREF(PyList_GetItem(flds, i));
  free(fields);

  // FREE PE
  Py_DECREF(pe);

  return out;
}
