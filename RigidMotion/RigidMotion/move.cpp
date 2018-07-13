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
 
# include "rigidMotion.h"

using namespace K_FLD;

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
  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType); 

  E_Int posx, posy, posz;
      
  if (res == 1 || res == 2)
  {
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      if (res == 2) delete cn;
      PyErr_SetString(PyExc_TypeError,
                      "move: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;
    
    FldArrayF* fn = new FldArrayF(*f);
    E_Int npts = f->getSize();
    E_Float* fx = f->begin(posx);
    E_Float* fy = f->begin(posy);
    E_Float* fz = f->begin(posz);
    E_Float* fnx = fn->begin(posx);
    E_Float* fny = fn->begin(posy);
    E_Float* fnz = fn->begin(posz);

#pragma omp parallel shared (npts, fnx, fny, fnz, fx, fy, fz, cx, cy, cz, dx, dy, dz, r11, r12, r13, r21, r22, r23, r31, r32, r33) if (npts > 100)
    {
      E_Float x, y, z;
#pragma omp for nowait
      for (E_Int ind = 0; ind < npts; ind++)
      {
        x = fx[ind]; y = fy[ind]; z = fz[ind];
        
        fnx[ind] = r11*(x-cx) + r12*(y-cy) + r13*(z-cz) + dx;
        fny[ind] = r21*(x-cx) + r22*(y-cy) + r23*(z-cz) + dy;
        fnz[ind] = r31*(x-cx) + r32*(y-cy) + r33*(z-cz) + dz;
      }
    }
    delete f;

    // Build array
    PyObject* tpl;
    if (res == 1) //structured
    {
      tpl = K_ARRAY::buildArray(*fn, varString, 
                                im, jm, km);
      delete fn;
    } 
    else //unstructured 
    {  
      tpl = K_ARRAY::buildArray(*fn, varString,
                                *cn, -1, eltType);
      delete fn; delete cn;
    }
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "move: not a valid array.");
    return NULL;
  }             
}
