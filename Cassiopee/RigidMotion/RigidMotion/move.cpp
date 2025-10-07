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

# include "rigidMotion.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
// Move a mesh with a defined motion
// Xp = R* (X-C) + d
//=============================================================================
PyObject* K_RIGIDMOTION::move(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float dx, dy, dz;
  E_Float cx, cy, cz;
  E_Float r11, r12, r13, r21, r22, r23, r31, r32, r33;

  if (!PYPARSETUPLE_(args,
                     O_ RRRR_ RRRR_ RRRR_ RRR_,
                     &array,
                     &dx, &dy, &dz,
                     &cx, &cy, &cz,
                     &r11, &r12, &r13,
                     &r21, &r22, &r23,
                     &r31, &r32, &r33)) return NULL;
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);
  if (res < 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "move: 1st arg is an invalid array.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f,cn);
    PyErr_SetString(PyExc_TypeError,
                    "move: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Int npts = f->getSize();
  E_Float* fx = f->begin(posx);
  E_Float* fy = f->begin(posy);
  E_Float* fz = f->begin(posz);

#pragma omp parallel shared (npts, fx, fy, fz, cx, cy, cz, dx, dy, dz, r11, r12, r13, r21, r22, r23, r31, r32, r33) if (npts > 100)
  {
    E_Float x, y, z;
  #pragma omp for
    for (E_Int ind = 0; ind < npts; ind++)
    {
      x = fx[ind]; y = fy[ind]; z = fz[ind];

      fx[ind] = r11*(x-cx) + r12*(y-cy) + r13*(z-cz) + dx;
      fy[ind] = r21*(x-cx) + r22*(y-cy) + r23*(z-cz) + dy;
      fz[ind] = r31*(x-cx) + r32*(y-cy) + r33*(z-cz) + dz;
    }
  }
  RELEASESHAREDB(res, array, f, cn);
  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
// Move coordinates defined by 3 numpys
// Xp = R* (X-C) + d
//=============================================================================
PyObject* K_RIGIDMOTION::moveN(PyObject* self, PyObject* args)
{
  PyObject *pyCoordsN;
  E_Float dx, dy, dz;
  E_Float cx, cy, cz;
  E_Float r11, r12, r13, r21, r22, r23, r31, r32, r33;

  if (!PYPARSETUPLE_(args, O_ TRRR_ TRRR_ RRRR_ RRRR_ R_,
                     &pyCoordsN,
                     &dx, &dy, &dz,
                     &cx, &cy, &cz,
                     &r11, &r12, &r13,
                     &r21, &r22, &r23,
                     &r31, &r32, &r33)) return NULL;
  if (PyList_Check(pyCoordsN) == 0)
  {
    PyErr_SetString(PyExc_TypeError,"moveN: 1st arg must be a list.");
    return NULL;
  }
  int size=PyList_Size(pyCoordsN);
  if (size != 3)
  {
    PyErr_SetString(PyExc_TypeError,"moveN: 1st arg must be a list of 3 elements.");
    return NULL;
  }
  vector<FldArrayF*> coords(size);
  vector<PyObject*> l(size);
  for (int i = 0; i < size; i++)
  {
    PyObject* tpl = PyList_GetItem(pyCoordsN,i);
    K_NUMPY::getFromNumpyArray(tpl, coords[i]);
    l[i]=tpl;
  }
  E_Float* xt = coords[0]->begin();
  E_Float* yt = coords[1]->begin();
  E_Float* zt = coords[2]->begin();

  E_Int npts = coords[0]->getSize() * coords[0]->getNfld();

#pragma omp parallel default(shared)
  {
    E_Float x, y, z;
  #pragma omp for
    for (E_Int ind = 0; ind < npts; ind++)
    {
      x = xt[ind]; y = yt[ind]; z = zt[ind];
      xt[ind] = r11*(x-cx) + r12*(y-cy) + r13*(z-cz) + dx;
      yt[ind] = r21*(x-cx) + r22*(y-cy) + r23*(z-cz) + dy;
      zt[ind] = r31*(x-cx) + r32*(y-cy) + r33*(z-cz) + dz;
    }
  }

  for (int i = 0; i < 3; i++) RELEASESHAREDN(l[i], coords[i]);
  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
// Compute the grid velocity :
// se = s0-omgp ^ c + omgp ^ x
//=============================================================================
PyObject* K_RIGIDMOTION::evalGridMotionN(PyObject* self, PyObject* args)
{
  E_Float s01, s02, s03;
  E_Float cx, cy, cz;
  E_Float omg1, omg2, omg3; // vitesse de rotation
  PyObject *pyCoordsN;
  PyObject *pySeN;
  if (!PYPARSETUPLE_(args, OO_ TRRR_ TRRR_ TRRR_,
                     &pyCoordsN, &pySeN,
                     &s01, &s02, &s03,
                     &cx, &cy, &cz,
                     &omg1, &omg2, &omg3)) return NULL;
  if (PyList_Check(pyCoordsN) == 0)
  {
    PyErr_SetString(PyExc_TypeError,"evalGridMotionN: 1st arg must be a list.");
    return NULL;
  }
  int size = PyList_Size(pyCoordsN);
  if (size != 3)
  {
    PyErr_SetString(PyExc_TypeError,"evalGridMotionN: 1st arg must be a list of 3 elements.");
    return NULL;
  }
  if (PyList_Check(pySeN) == 0)
  {
    PyErr_SetString(PyExc_TypeError,"evalGridMotionN: 2nd arg must be a list.");
    return NULL;
  }
  int size2 = PyList_Size(pySeN);
  if (size2 != 3)
  {
    PyErr_SetString(PyExc_TypeError,"evalGridMotionN: 2nd arg must be a list of 3 elements.");
    return NULL;
  }
  vector<FldArrayF*> coords(size);
  vector<PyObject*> l(size);
  vector<PyObject*> l2(size);
  vector<FldArrayF*> se(size);
  for (int i = 0; i < size; i++)
  {
    PyObject* tpl = PyList_GetItem(pyCoordsN,i);
    K_NUMPY::getFromNumpyArray(tpl, coords[i]);
    l[i] = tpl;
    tpl = PyList_GetItem(pySeN, i);
    K_NUMPY::getFromNumpyArray(tpl, se[i]);
    l2[i] = tpl;
  }

  E_Float* xt = coords[0]->begin();
  E_Float* yt = coords[1]->begin();
  E_Float* zt = coords[2]->begin();
  E_Int npts = coords[0]->getSize() * coords[0]->getNfld();
  E_Float* se1 = se[0]->begin();
  E_Float* se2 = se[1]->begin();
  E_Float* se3 = se[2]->begin();

  E_Float tx = s01 - (omg2 * cz - omg3 * cy);
  E_Float ty = s02 - (omg3 * cx - omg1 * cz);
  E_Float tz = s03 - (omg1 * cy - omg2 * cx);

#pragma omp parallel default(shared)
  {
  #pragma omp for
    for (E_Int ind = 0; ind < npts; ind++)
    {
      se1[ind] = tx + (omg2 * zt[ind] - omg3 * yt[ind]);
      se2[ind] = ty + (omg3 * xt[ind] - omg1 * zt[ind]);
      se3[ind] = tz + (omg1 * yt[ind] - omg2 * xt[ind]);
    }
  }

  for (int i = 0; i < 3; i++)
  {
    RELEASESHAREDN(l[i], coords[i]);
    RELEASESHAREDN(l2[i], se[i]);
  }
  Py_INCREF(Py_None);
  return Py_None;
}
