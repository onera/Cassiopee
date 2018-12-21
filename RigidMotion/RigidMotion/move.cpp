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
 
# include "rigidMotion.h"

using namespace K_FLD;
using namespace std;
//=============================================================================
// Move a mesh with a defined motion
//=============================================================================
PyObject* K_RIGIDMOTION::move(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float dx, dy, dz;
  E_Float cx, cy, cz;
  E_Float r11, r12, r13, r21, r22, r23, r31, r32, r33;

  if (!PYPARSETUPLEF(args, 
                     "Oddddddddddddddd", 
                     "Offfffffffffffff",
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
  E_Int res = K_ARRAY::getFromArray2(array, varString, f, im, jm, km, cn, eltType); 
  if ( res < 0 )
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
  #pragma omp for nowait
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
//=============================================================================
PyObject* K_RIGIDMOTION::moveN(PyObject* self, PyObject* args)
{
  PyObject *pyCoordsN;
  E_Float dx, dy, dz;
  E_Float cx, cy, cz;
  E_Float r11, r12, r13, r21, r22, r23, r31, r32, r33;

  if (!PYPARSETUPLEF(args, 
                     "O(ddd)(ddd)ddddddddd", 
                     "O(fff)(fff)fffffffff",
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
  if ( size != 3)
  {
    PyErr_SetString(PyExc_TypeError,"moveN: 1st arg must be a list of 3 elements.");    
    return NULL;
  }
  vector<FldArrayF*> coords(size);
  vector<PyObject*> l(size);
  for (int i = 0; i < size; i++)
  {
    PyObject* tpl = PyList_GetItem(pyCoordsN,i);
    K_NUMPY::getFromNumpyArray(tpl, coords[i], true);
    l[i]=tpl;
  }
  E_Float* xt = coords[0]->begin();
  E_Float* yt = coords[1]->begin();
  E_Float* zt = coords[2]->begin();

  E_Int npts = coords[0]->getSize();

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

  for (int i = 0; i < 3; i++)
    RELEASESHAREDN(l[i], coords[i]);
  Py_INCREF(Py_None);
  return Py_None;
}
