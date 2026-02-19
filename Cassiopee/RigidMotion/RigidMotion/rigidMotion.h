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

# ifndef _RIGIDMOTION_RIGIDMOTION_H_
# define _RIGIDMOTION_RIGIDMOTION_H_

# include "kcore.h"

namespace K_RIGIDMOTION
{
  PyObject* move(PyObject* self, PyObject* args);
  PyObject* moveN(PyObject* self, PyObject* args);
  PyObject* evalGridMotionN(PyObject* self, PyObject* args);

  PyObject* _computeRotorMotionZ(PyObject* self, PyObject* args);
  //Return a list of numpys: [r0,x0,rotMat,s0]
  PyObject* _computeRotorMotionInfo(PyObject* self, PyObject* args);
  PyObject* evalSpeed3(PyObject* self, PyObject* args);
  PyObject* copyCoords(PyObject* self, PyObject* args);
}
#endif
