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

#include "xcore.h"

PyObject* K_XCORE::CreateAdaptMesh(PyObject *self, PyObject *args)
{
    PyErr_SetString(PyExc_TypeError,
        "CreateAdaptMesh: not available (no mpi).");
    return NULL;
}

PyObject* K_XCORE::AdaptMesh(PyObject *self, PyObject *args)
{
    PyErr_SetString(PyExc_TypeError,
        "AdaptMesh: not available (no mpi).");
    return NULL;
}

PyObject *K_XCORE::computeHessian(PyObject *self, PyObject *args)
{
    PyErr_SetString(PyExc_TypeError,
        "computeHessian: not available (no mpi).");
    return NULL;
}

PyObject *K_XCORE::computeGradient(PyObject *self, PyObject *args)
{
    PyErr_SetString(PyExc_TypeError,
        "computeGradient: not available (no mpi).");
    return NULL;
}

PyObject *K_XCORE::hessianToMetric(PyObject *self, PyObject *args)
{
    PyErr_SetString(PyExc_TypeError,
        "hessianToMetric: not available (no mpi).");
    return NULL;
}

PyObject *K_XCORE::_makeRefDataFromGradAndHess(PyObject *self, PyObject *args)
{
    PyErr_SetString(PyExc_TypeError,
        "_makeRefDataFromGradAndHess: not available (no mpi).");
    return NULL;
}

PyObject *K_XCORE::_prepareMeshForAdaptation(PyObject *self, PyObject *args)
{
    PyErr_SetString(PyExc_TypeError,
        "_prepareMeshForAdaptation: not available (no mpi).");
    return NULL;
}

PyObject *K_XCORE::ExtractLeafMesh(PyObject *self, PyObject *args)
{
    PyErr_SetString(PyExc_TypeError,
        "ExtractLeafMesh: not available (no mpi).");
    return NULL;
}

PyObject *K_XCORE::_assignRefDataToAM(PyObject *self, PyObject *args)
{
    PyErr_SetString(PyExc_TypeError,
        "_assignRefDataToAM: not available (no mpi).");
    return NULL;
}

PyObject *K_XCORE::extractBoundaryMesh(PyObject *self, PyObject *args)
{
    PyErr_SetString(PyExc_TypeError,
        "extractBoundaryMesh: not available (no mpi).");
    return NULL;
}