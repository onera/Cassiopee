/*
    Copyright 2013-2018 Onera.

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
using namespace K_FLD;
using namespace std;

// ============================================================================
/* cyl2Cart in place, array version */
// ============================================================================
PyObject* K_TRANSFORM::_cyl2CartA(PyObject* self, PyObject* args)
{
    PyObject *array;
    E_Float X0, Y0, Z0;
    E_Float ex, ey, ez;
    if (!PYPARSETUPLEF(args,"O(ddd)(ddd)", "O(fff)(fff)",
                       &array, &X0, &Y0, &Z0, &ex, &ey, &ez))
      return NULL;

    // check array
    E_Int im, jm, km;
    FldArrayF* f; FldArrayI* cn;
    char* varString; char* eltType;
    E_Int res = K_ARRAY::getFromArray2(array, varString, f, im, jm, km,
                                       cn, eltType);
    if (res != 1 && res != 2)
    {
      PyErr_SetString(PyExc_TypeError,
                      "cyl2CartA: 1st arg is an invalid array.");
      return NULL;
    }

    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "cyl2Cart: cannot find coordinates in zone.");
      RELEASESHAREDB(res, array, f, cn); return NULL;
    }
    posx++; posy++; posz++;
    E_Int npts = f->getSize();

    E_Float eps = K_CONST::E_GEOM_CUTOFF;
    E_Int posxc, posyc, poszc;

    if (ex > eps && ey < eps && ez < eps)
    {
      posxc = posz; posyc = posy; poszc = posx;
    }
    else if (ey > eps && ex < eps && ez < eps)
    {
      posxc = posx; posyc = posz; poszc = posy;
    }
    else if (ez > eps && ey < eps && ex < eps)
    {
      posxc = posx; posyc = posy; poszc = posz;
    }
    else
    {
      PyErr_SetString(PyExc_TypeError,
                      "cyl2Cart: axis must be x,y or z.");
      RELEASESHAREDB(res, array, f, cn); return NULL;
    }
    E_Float* xt = f->begin(posxc);
    E_Float* yt = f->begin(posyc);
    //E_Float* zt = f->begin(poszc);

#pragma omp parallel default(shared)
    {
#pragma omp for
      for (E_Int ind = 0; ind < npts; ind++)
      {
        E_Float r     = xt[ind];
        E_Float theta = yt[ind];
        E_Float x = r*cos(theta);
        E_Float y = r*sin(theta);
        xt[ind] = x;
        yt[ind] = y;
      }
    }
    RELEASESHAREDB(res, array, f, cn);

    Py_INCREF(Py_None);
    return Py_None;
}

// ============================================================================
/*  cyl2Cart: in place, pyTree version */
// ============================================================================
PyObject* K_TRANSFORM::_cyl2CartZ(PyObject* self, PyObject* args)
{
    PyObject *zone;
    char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
    E_Float X0, Y0, Z0;
    E_Float ex, ey, ez;
    if (!PYPARSETUPLEF(args,"O(ddd)(ddd)sss", "O(fff)(fff)sss",
                       &zone, &X0, &Y0, &Z0, &ex, &ey, &ez, &GridCoordinates, &FlowSolutionNodes, &FlowSolutionCenters))
        return NULL;

    vector<PyArrayObject*> hook;
    E_Int im, jm, km, cnSize, cnNfld;
    char* varString; char* eltType;
    vector<E_Float*> fields; vector<E_Int> locs;
    vector<E_Int*> cn;
    E_Int res = K_PYTREE::getFromZone(zone, 1, 0, varString, fields, locs,
                                      im, jm, km,
                                      cn, cnSize, cnNfld, eltType, hook,
                                      GridCoordinates,
                                      FlowSolutionNodes, FlowSolutionCenters);
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
      PyErr_SetString(PyExc_TypeError,
                      "cyl2Cart: cannot find coordinates in zone.");
      return NULL;
    }

    E_Int npts;
    if (res == 1) npts = im*jm*km;
    else npts = im;

    E_Float eps = K_CONST::E_GEOM_CUTOFF;
    E_Int posR, posTHETA, posZ;
    if (ex > eps && ey < eps && ez < eps)// AXE (OX)
    {
      posR = posy; posTHETA = posz; posZ = posx;
    }
    else if (ey > eps && ex < eps && ez < eps)// AXE (OY)
    {
      posR = posx; posTHETA = posz; posZ = posy;
    }
    else if (ez > eps && ey < eps && ex < eps)// AXE (OZ)
    {
      posR = posx; posTHETA = posy; posZ = posz;
    }
    else
    {
      RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
      PyErr_SetString(PyExc_TypeError,
                      "cyl2Cart: axis must be x,y or z.");
      return NULL;
    }
    E_Float* xt = fields[posR];
    E_Float* yt = fields[posTHETA];
    //E_Float* zt = fields[posZ];
#pragma omp parallel default(shared)
    {
#pragma omp for
      for (E_Int ind = 0; ind < npts; ind++)
      {
        E_Float r     = xt[ind];
        E_Float theta = yt[ind];
        E_Float x = r*cos(theta);
        E_Float y = r*sin(theta);
        xt[ind] = x;
        yt[ind] = y;
      }
    }

    delete [] eltType;
    RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
    Py_INCREF(Py_None);
    return Py_None;
}
