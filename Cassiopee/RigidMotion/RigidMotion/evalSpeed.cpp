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
// Eval part of speed (motion3)
//=============================================================================
PyObject* K_RIGIDMOTION::evalSpeed3(PyObject* self, PyObject* args)
{
  PyObject *xo, *yo, *zo;
  PyObject *sxo, *syo, *szo;
  E_Float omega, teta;
  E_Float tx, ty, tz;
  E_Float cx, cy, cz;
  E_Float kx, ky, kz;

  if (!PYPARSETUPLE_(args, OOOO_ OO_ RRRR_ RRRR_ RRR_,
                     &xo, &yo, &zo,
                     &sxo, &syo, &szo,
                     &omega, &teta,
                     &tx, &ty, &tz,
                     &cx, &cy, &cz,
                     &kx, &ky, &kz)) return NULL;
  // Check numpys
  E_Float* x; E_Int size;
  K_NUMPY::getFromNumpyArray(xo, x, size);
  E_Float* y;
  K_NUMPY::getFromNumpyArray(yo, y, size);
  E_Float* z;
  K_NUMPY::getFromNumpyArray(zo, z, size);
  E_Float* sx;
  K_NUMPY::getFromNumpyArray(sxo, sx, size);
  E_Float* sy;
  K_NUMPY::getFromNumpyArray(syo, sy, size);
  E_Float* sz;
  K_NUMPY::getFromNumpyArray(szo, sz, size);
  
#pragma omp parallel
  {
    //E_Float sin_teta, cos_teta, kcm, 
    E_Float cmx, cmy, cmz;
    E_Float kvcmx, kvcmy, kvcmz;

    //sin_teta = sin(teta);
    //cos_teta = cos(teta);

#pragma omp for
    for (E_Int i = 0; i < size; i++)
    {
      // cm
      cmx = x[i] - cx;
      cmy = y[i] - cy;
      cmz = z[i] - cz;
    
      // k.cm
      //kcm = kx*cmx + ky*cmy + kz*cmz;
    
      // k x cm
      kvcmx = ky*cmz-kz*cmy;
      kvcmy = kz*cmx-kx*cmz;
      kvcmz = kx*cmy-ky*cmx;

      // grid speed
      //sx[i] = tx-omega*sin_teta*(cmx-kcm*kx)+omega*cos_teta*kvcmx;
      //sy[i] = ty-omega*sin_teta*(cmy-kcm*ky)+omega*cos_teta*kvcmy;
      //sz[i] = tz-omega*sin_teta*(cmz-kcm*kz)+omega*cos_teta*kvcmz;
      sx[i] = tx + omega*kvcmx;
      sy[i] = ty + omega*kvcmy;
      sz[i] = tz + omega*kvcmz;
      //printf("%f %f %f\n",sx[i],sy[i],sz[i]);
    }
  }
  Py_DECREF(xo); Py_DECREF(yo); Py_DECREF(zo);
  Py_DECREF(sxo); Py_DECREF(syo); Py_DECREF(szo);
  Py_INCREF(Py_None);
  return Py_None; 
}
