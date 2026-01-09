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
#ifndef _XCORE_XCORE_H_
#define _XCORE_XCORE_H_

#ifdef E_DOUBLEINT
#define XMPI_INT MPI_INT64_T
#else
#define XMPI_INT MPI_INT
#endif

#include "kcore.h"
#include <limits>

constexpr E_Float EFLOATMIN = -std::numeric_limits<E_Float>::max();
constexpr E_Float EFLOATMAX = std::numeric_limits<E_Float>::max();

constexpr E_Int EINTMIN = std::numeric_limits<E_Int>::min();
constexpr E_Int EINTMAX = std::numeric_limits<E_Int>::max();

namespace K_XCORE
{
    //PyObject* setBCDataInGhostCellsStruct(PyObject* self, PyObject* args);
    PyObject *zoltan1(PyObject *self, PyObject *args);

    PyObject *chunk2partNGon(PyObject *self, PyObject *args);
    PyObject *chunk2partElt(PyObject *self, PyObject *args);

    PyObject *exchangeFields(PyObject *self, PyObject *args);

    PyObject *AdaptMesh_Init(PyObject *self, PyObject *args);
    PyObject *AdaptMesh_AssignRefData(PyObject *self, PyObject *args);
    PyObject *AdaptMesh_LoadBalance(PyObject *self, PyObject *args);
    PyObject *AdaptMesh_Adapt(PyObject *self, PyObject *args);
    PyObject *AdaptMesh_ExtractMesh(PyObject *self, PyObject *args);
    PyObject *AdaptMesh_Exit(PyObject *self, PyObject *args);
    PyObject *AdaptMesh_ExtractOwners(PyObject *self, PyObject *args);
    PyObject *AdaptMesh_ExtractNeighbours(PyObject *self, PyObject *args);
    PyObject *AdaptMesh_ExtractCellLevels(PyObject *self, PyObject *args);
    PyObject *AdaptMesh_ExtractCellRanges(PyObject *self, PyObject *args);
    PyObject *AdaptMesh_ExtractHaloCellLevels(PyObject *self, PyObject *args);

    PyObject *removeIntersectingKPlanes(PyObject *self, PyObject *args);

    PyObject *IntersectMesh_Init(PyObject *self, PyObject *args);
    PyObject *IntersectMesh_ExtractMesh(PyObject *self, PyObject *args);
    PyObject *IntersectMesh_TriangulateFaceSet(PyObject *self, PyObject *args);
    PyObject *IntersectMesh_Exit(PyObject *self, PyObject *args);
    PyObject *IntersectMesh_ExtractFaceSet(PyObject *self, PyObject *args);
    
    PyObject *icapsule_init2(PyObject *self, PyObject *args);
    PyObject *icapsule_set_master(PyObject *self, PyObject *args);
    PyObject *icapsule_set_slaves(PyObject *self, PyObject *args);
    PyObject *icapsule_adapt2(PyObject *self, PyObject *args);
    PyObject *icapsule_intersect2(PyObject *self, PyObject *args);

    PyObject *split_connex(PyObject *self, PyObject *args);

    PyObject *icapsule_init(PyObject *self, PyObject *args);
    PyObject *icapsule_adapt(PyObject *self, PyObject *args);
    PyObject *icapsule_intersect(PyObject *self, PyObject *args);
    PyObject *icapsule_extract_master(PyObject *self, PyObject *args);
    PyObject *icapsule_extract_slave(PyObject *self, PyObject *args);
    PyObject *icapsule_extract_slaves(PyObject *self, PyObject *args);

    PyObject *triangulate_skin(PyObject *self, PyObject *args);

    PyObject *extractCell(PyObject *self, PyObject *args);

    PyObject *extractFacesFromPointTag(PyObject *self, PyObject *args);
}

#endif
