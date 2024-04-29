/*
    Copyright 2013-2024 Onera.

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
#ifndef _XCORE_XCORE_H_
#define _XCORE_XCORE_H_

#ifdef E_DOUBLEINT
#define XMPI_INT MPI_INT64_T
#else
#define XMPI_INT MPI_INT
#endif

# include "kcore.h"
namespace K_XCORE
{
  //PyObject* setBCDataInGhostCellsStruct(PyObject* self, PyObject* args);
  PyObject *zoltan1(PyObject *self, PyObject *args);

  PyObject *chunk2partNGon(PyObject *self, PyObject *args);
  PyObject *chunk2partElt(PyObject *self, PyObject *args);

  PyObject *exchangeFields(PyObject *self, PyObject *args);

  PyObject *adaptMesh(PyObject *self, PyObject *args);
  
  PyObject *adaptMeshSeq(PyObject *self, PyObject *args);
  PyObject *extractLeafMesh(PyObject *self, PyObject *args);
  PyObject *createAdaptMesh(PyObject *self, PyObject *args);

  PyObject *adaptMeshDir(PyObject *self, PyObject *args);

  PyObject *AdaptMesh(PyObject *self, PyObject *args);
  PyObject *CreateAdaptMesh(PyObject *self, PyObject *args);
  PyObject *computeHessian(PyObject *self, PyObject *args);
  PyObject *computeGradient(PyObject *self, PyObject *args);
  PyObject *hessianToMetric(PyObject *self, PyObject *args);
  PyObject *_makeRefDataFromGradAndHess(PyObject *self, PyObject *args);
  PyObject *_prepareMeshForAdaptation(PyObject *self, PyObject *args);
  PyObject *ExtractLeafMesh(PyObject *self, PyObject *args);
  PyObject *_assignRefDataToAM(PyObject *self, PyObject *args);
  PyObject *extractBoundaryMesh(PyObject *self, PyObject *args);

  PyObject *intersectSurf(PyObject *self, PyObject *args);
  PyObject *removeIntersectingKPlanes(PyObject *self, PyObject *args);
}
#endif
