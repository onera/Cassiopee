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

# include "transform.h"
using namespace K_FLD;
using namespace std;

// ============================================================================
/* cart2cyl in place, array version */
// ============================================================================
PyObject* K_TRANSFORM::_cart2CylA(PyObject* self, PyObject* args)
{
    PyObject *array;
    E_Float X0, Y0, Z0;
    E_Float ex, ey, ez;
    E_Int depth;
    E_Float thetaShift;
    if (!PYPARSETUPLE_(args, O_ TRRR_ TRRR_ I_ R_,
                      &array, &X0, &Y0, &Z0, &ex, &ey, &ez, &depth, &thetaShift))
        return NULL;

    // check array
    E_Int im, jm, km;
    FldArrayF* f; FldArrayI* cn;
    char* varString; char* eltType;
    E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km,
                                       cn, eltType);
    if (res != 1 && res != 2)
    {
      PyErr_SetString(PyExc_TypeError,
                      "cart2Cyl: 1st arg is an invalid array.");
      return NULL;
    }

    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "cart2Cyl: cannot find coordinates in zone.");
      RELEASESHAREDB(res, array, f, cn); return NULL;
    }
    posx++; posy++; posz++;
    E_Int npts = f->getSize();

    E_Float eps = 1.e-12;
    E_Float *rt, *thetat;
    if (ex > eps && ey < eps && ez < eps) // axe X
    {
      rt = f->begin(posy); thetat = f->begin(posz);
    }
    else if (ey > eps && ex < eps && ez < eps) // axe Y
    {
      rt = f->begin(posx); thetat = f->begin(posz);
    }
    else if (ez > eps && ey < eps && ex < eps) // axe Z
    {
      rt = f->begin(posx); thetat = f->begin(posy);
    }
    else
    {
      PyErr_SetString(PyExc_TypeError,
                      "cart2Cyl: axis must be canonical.");
      RELEASESHAREDB(res, array, f, cn); return NULL;
    }
    E_Int ret = K_LOC::cart2Cyl(npts, f->begin(posx), f->begin(posy), f->begin(posz),
                                X0, Y0, Z0, ex, ey, ez, rt, thetat,
                                im, jm, km, depth, thetaShift);
    if (ret == 1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "cart2Cyl: axis must be canonical.");
      RELEASESHAREDB(res, array, f, cn); return NULL;
    }

    RELEASESHAREDB(res, array, f, cn);

    Py_INCREF(Py_None);
    return Py_None;
}

// ============================================================================
/*  cart2Cyl: in place, pyTree version */
// ============================================================================
PyObject* K_TRANSFORM::_cart2CylZ(PyObject* self, PyObject* args)
{
    PyObject *zone;
    char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
    E_Float X0, Y0, Z0;
    E_Float ex, ey, ez;
    E_Int depth; E_Float thetaShift;
    if (!PYPARSETUPLE_(args, O_ TRRR_ TRRR_ I_ R_ SSS_,
                      &zone, &X0, &Y0, &Z0, &ex, &ey, &ez,
                      &depth, &thetaShift, &GridCoordinates, &FlowSolutionNodes, &FlowSolutionCenters))
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
      if (res == 2) delete [] eltType;
      RELEASESHAREDZ(hook, varString, (char*)NULL);
      PyErr_SetString(PyExc_TypeError,
                      "cart2Cyl: cannot find coordinates in zone.");
      return NULL;
    }

    E_Int npts;
    if (res == 1) npts = im*jm*km;
    else npts = im;
    E_Float eps = 1.e-12;

    E_Float *rt, *thetat;
    thetat = NULL; rt = NULL;
    if (ex > eps && ey < eps && ez < eps) // axe X
    {
      rt = fields[posy]; thetat = fields[posz];
    }
    else if (ey > eps && ex < eps && ez < eps) // axe Y
    {
      rt = fields[posz]; thetat = fields[posx];
    }
    else if (ez > eps && ey < eps && ex < eps) // axe Z
    {
      rt = fields[posx]; thetat = fields[posy];
    }

    E_Int ret = K_LOC::cart2Cyl(npts, fields[posx], fields[posy], fields[posz],
                                X0, Y0, Z0, ex, ey, ez, rt, thetat,  im, jm, km, depth, thetaShift);
    if (ret == 1)
    {
      if (res == 2) delete [] eltType;
      RELEASESHAREDZ(hook, varString, (char*)NULL);
      PyErr_SetString(PyExc_TypeError,
                      "cart2Cyl: axis must be canonical.");
      return NULL;
    }
    if (res == 2) delete [] eltType;
    RELEASESHAREDZ(hook, varString, (char*)NULL);
    Py_INCREF(Py_None);
    return Py_None;
}
