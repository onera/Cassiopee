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

# include "transform.h"
using namespace std;
using namespace K_FLD;
using namespace K_SEARCH;

// ============================================================================
/* perform a cylindric radial projection of a vector field */
// ============================================================================
PyObject* K_TRANSFORM::_alignVectorFieldWithRadialCylindricProjection(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float CenterX, CenterY, CenterZ;
  E_Float AxisX, AxisY, AxisZ; // axis unitary vector
  PyObject* varList;
  if (!PYPARSETUPLE_(args, O_ TRRR_ TRRR_ O_,
                    &array, &CenterX, &CenterY, &CenterZ, &AxisX, &AxisY, &AxisZ, &varList))
  {
    return NULL;
  }


  // Check varList
  if (PyList_Check(varList) == false)
  {
    PyErr_SetString(PyExc_TypeError,
                    "alignVectorFieldFollowingRadialCylindricProjection: varList must be a list of variables.");
    return NULL;
  }

  // nbre de variables a lisser
  E_Int nvars = PyList_Size(varList);

  if (nvars != 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "alignVectorFieldFollowingRadialCylindricProjection: must provide exactly 3 fields.");
    return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString,
                                     f, im, jm, km, cn, eltType);

  // position des variables
  vector<E_Int> posVars(nvars);
  for (E_Int i = 0; i < nvars; i++)
  {
      PyObject* varname = PyList_GetItem(varList, i);
      char* var = NULL;
      if (PyString_Check(varname)) var = PyString_AsString(varname);
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(varname)) var = (char*)PyUnicode_AsUTF8(varname);
#endif
      E_Int pos = K_ARRAY::isNamePresent(var, varString);
      if (pos == -1)
      {
          PyErr_SetString(PyExc_TypeError,
                          "alignVectorFieldFollowingRadialCylindricProjection: variable doesn't exist.");
        return NULL;
      }
      posVars[i] = pos;
  }
  E_Float* vx = NULL; E_Float* vy = NULL; E_Float* vz = NULL;
  vx = f->begin(posVars[0]+1);
  vy = f->begin(posVars[1]+1);
  vz = f->begin(posVars[2]+1);


  // Get coordinates
  E_Float* cx = NULL; E_Float* cy = NULL; E_Float* cz = NULL;
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "alignVectorFieldFollowingRadialCylindricProjection: requires coordinates.");
    return NULL;
  }
  cx = f->begin(posx+1);
  cy = f->begin(posy+1);
  cz = f->begin(posz+1);


  E_Int npts = f->getSize();

  E_Float qx, qy, qz;
  qx = CenterX + AxisX;
  qy = CenterY + AxisY;
  qz = CenterZ + AxisZ;

  E_Float qcx, qcy, qcz;
  qcx = CenterX - qx;
  qcy = CenterY - qy;
  qcz = CenterZ - qz;

  #pragma omp parallel
  {
  E_Float qpx, qpy, qpz;
  E_Float qpnorm;
  E_Float alignment;
  E_Float alignmentTol(1.0-1.e-10);
  E_Float nx, ny, nz;
  E_Float bx, by, bz;
  E_Float tx, ty, tz;
  E_Float vnorm, tnorm;
  E_Float vxi, vyi, vzi;

  #pragma omp for
  for (E_Int i = 0; i < npts; i++)
  {
    vxi = vx[i];
    vyi = vy[i];
    vzi = vz[i];
    vnorm = sqrt(vxi*vxi + vyi*vyi + vzi*vzi);

    qpx = cx[i] - qx;
    qpy = cy[i] - qy;
    qpz = cz[i] - qz;
    qpnorm = sqrt(qpx*qpx + qpy*qpy + qpz*qpz);
    qpx /= qpnorm;
    qpy /= qpnorm;
    qpz /= qpnorm;
    alignment = abs(qpx*AxisX + qpy*AxisY + qpz*AxisZ);
    if (alignment >= alignmentTol) continue;

    nx = qcy*qpz - qcz*qpy; //
    ny = qcz*qpx - qcx*qpz; // qc x qp
    nz = qcx*qpy - qcy*qpx; //

    bx = ny*vzi - nz*vyi; //
    by = nz*vxi - nx*vzi; // n x v
    bz = nx*vyi - ny*vxi; //

    tx = by*nz - bz*ny; //
    ty = bz*nx - bx*nz; // b x n
    tz = bx*ny - by*nx; //

    tnorm = sqrt(tx*tx + ty*ty + tz*tz);
    tx /= tnorm;
    ty /= tnorm;
    tz /= tnorm;

    vx[i] = tx * vnorm;
    vy[i] = ty * vnorm;
    vz[i] = tz * vnorm;

  }
  }

  RELEASESHAREDB(res, array, f, cn);
  Py_INCREF(Py_None);
  return Py_None;
}
