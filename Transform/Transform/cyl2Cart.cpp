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
    E_Int posr, postheta;

    if (ex > eps && ey < eps && ez < eps)
    {
      posr = posy; postheta = posz; 
    }
    else if (ey > eps && ex < eps && ez < eps)
    {
      posr = posz; postheta = posx;
    }
    else if (ez > eps && ey < eps && ex < eps)
    {
      posr = posx; postheta = posy;
    }
    else
    {
      PyErr_SetString(PyExc_TypeError,
                      "cyl2Cart: axis must be x,y or z.");
      RELEASESHAREDB(res, array, f, cn); return NULL;
    }
    E_Float* rt = f->begin(posr);
    E_Float* thetat = f->begin(postheta);
    
#pragma omp parallel default(shared)
    {
#pragma omp for
      for (E_Int ind = 0; ind < npts; ind++)
      {
        E_Float r     = rt[ind];
        E_Float theta = thetat[ind];
        rt[ind] = r*cos(theta);
        thetat[ind] = r*sin(theta);//(y,z) dans le cas (X;R;Theta), (x,y) dans le cas (R,Theta,Z), (z,x) dans le cas (Theta,Y,R)
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
    E_Int posR, posTHETA;
    if (ex > eps && ey < eps && ez < eps)// AXE (OX)
    {
      posR = posy; posTHETA = posz;
    }
    else if (ey > eps && ex < eps && ez < eps)// AXE (OY)
    {
      posR = posz; posTHETA = posx;
    }
    else if (ez > eps && ey < eps && ex < eps)// AXE (OZ)
    {
      posR = posx; posTHETA = posy;
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

#pragma omp parallel default(shared)
    {
#pragma omp for
      for (E_Int ind = 0; ind < npts; ind++)
      {
        E_Float r     = xt[ind];
        E_Float theta = yt[ind];
        xt[ind] = r*cos(theta);
        yt[ind] = r*sin(theta);
      }//(y,z) dans le cas (X;R;Theta), (x,y) dans le cas (R,Theta,Z), (z,x) dans le cas (Theta,Y,R)
    }

    delete [] eltType;
    RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
    Py_INCREF(Py_None);
    return Py_None;
}
