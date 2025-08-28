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
// Copy les 3 numpys en entree dans les 3 numpys en sortie
//=============================================================================
PyObject* K_RIGIDMOTION::copyCoords(PyObject* self, PyObject* args)
{
    PyObject *xin, *yin, *zin;
    PyObject *xout, *yout, *zout;
    
    if (!PyArg_ParseTuple(args, "OOOOOO", &xin, &yin, &zin, &xout, &yout, &zout)) return NULL;
    
    E_Int size;
    E_Float *xi, *yi, *zi;
    K_NUMPY::getFromNumpyArray(xin, xi, size);
    K_NUMPY::getFromNumpyArray(yin, yi, size);
    K_NUMPY::getFromNumpyArray(zin, zi, size);
    E_Float *xo, *yo, *zo;
    K_NUMPY::getFromNumpyArray(xout, xo, size);
    K_NUMPY::getFromNumpyArray(yout, yo, size);
    K_NUMPY::getFromNumpyArray(zout, zo, size);

#pragma omp parallel
    {
#pragma omp for
        for (E_Int i = 0; i < size; i++)
        {
            xo[i] = xi[i];
        }
#pragma omp for
        for (E_Int i = 0; i < size; i++)
        {
            yo[i] = yi[i];
        }
#pragma omp for
        for (E_Int i = 0; i < size; i++)
        {
            zo[i] = zi[i];
        }
    }
    Py_DECREF(xin); Py_DECREF(yin); Py_DECREF(zin);
    Py_DECREF(xout); Py_DECREF(yout); Py_DECREF(zout);
    Py_INCREF(Py_None);
    return Py_None;
}