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
#ifndef _KCORE_H_
#define _KCORE_H_

#include "Array/Array.h"
#include "Def/DefTypes.h"
#include "Def/DefFunction.h"
#include "Numpy/Numpy.h"
#include "PyTree/PyTree.h"
#include "String/kstring.h"
#include "Interp/Interp.h"
#include "CompGeom/compGeom.h"
#include "Connect/connect.h"
#include "Sort/sort.h"
#include "Noise/noise.h"
#include "Loc/loc.h"
#include "Linear/linear.h"
#include "Metric/metric.h"
#include "Math/math.h"
#include "Nuga/include/KdTree.h"
#include "parallel.h"
#include "kPython.h"
#include <vector>

namespace K_KCORE
{
  void testFooKCore();
  PyObject* isCoordinateXPresent(PyObject* self, PyObject* args);
  PyObject* isCoordinateYPresent(PyObject* self, PyObject* args);
  PyObject* isCoordinateZPresent(PyObject* self, PyObject* args);
  PyObject* isNamePresent(PyObject* self, PyObject* args);
  PyObject* indiceStruct2Unstr(PyObject* self, PyObject* args);
  PyObject* indiceStruct2Unstr2(PyObject* self, PyObject* args);
  PyObject* indiceFace2Connect(PyObject* self, PyObject* args);
  PyObject* setOmpMaxThreads(PyObject* self, PyObject* args);
  PyObject* getOmpMaxThreads(PyObject* self, PyObject* args);
  PyObject* empty(PyObject* self, PyObject* args);
  PyObject* tester(PyObject* self, PyObject* args);
  PyObject* testerAcc(PyObject* self, PyObject* args);
  PyObject* copyto(PyObject* self, PyObject* args);
  PyObject* copyfrom(PyObject* self, PyObject* args);
  void memcpy__(E_Int* a, E_Int* b, E_Int s);
  void memcpy__(E_Float* a, E_Float* b, E_Int s);
} 

#endif
