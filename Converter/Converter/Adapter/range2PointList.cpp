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
# include "converter.h"

using namespace K_FLD;

//=============================================================================
/* Convert a Range of a structured block into a face indices numpy. 
   The return indices of faces start 0 and contains i faces, then j
   faces and then k faces */
//=============================================================================
PyObject* K_CONVERTER::range2PointList(PyObject* self, PyObject* args)
{
  E_Int imin, imax, jmin, jmax, kmin, kmax, ni, nj, nk;
  if (!PYPARSETUPLEI(args, "lllllllll", "iiiiiiiii", &imin, &imax, &jmin, &jmax, &kmin, &kmax,
                     &ni, &nj, &nk)) return NULL;
  
  E_Int ni1 = std::max(ni-1,1);
  E_Int nj1 = std::max(nj-1,1);
  E_Int nk1 = std::max(nk-1,1);
  E_Int shift, size, ind;
  E_Int ii = 0;
  if ( kmin == kmax && nk == 1) kmax = 2;// 2D

  // fenetre en i
  if (imin == imax)
  {
    size = (kmax-kmin)*(jmax-jmin);
    PyObject* o = K_NUMPY::buildNumpyArray(size, 1, 1);
    E_Int* p = K_NUMPY::getNumpyPtrI(o);
    for (E_Int k = kmin-1; k < kmax-1; k++)
    {
      for (E_Int j = jmin-1; j < jmax-1; j++)
      {
        ind = imin-1 + j*ni + k*ni*nj1;
        p[ii] = ind; ii++;
      }
    }
    return o;
  }
  // fenetre en j
  else if (jmin == jmax)
  {
    shift = ni*nj1*nk1;
    size = (kmax-kmin)*(imax-imin);
    PyObject* o = K_NUMPY::buildNumpyArray(size, 1, 1);
    E_Int* p = K_NUMPY::getNumpyPtrI(o);
    for (E_Int k = kmin-1; k < kmax-1; k++)
    {
      for (E_Int i = imin-1; i < imax-1; i++)
      {
        ind = shift + i + (jmin-1)*ni1 + k*ni1*nj;
        p[ii] = ind; ii++;
      }
    }
    return o;
  }
  // fenetre en k
  else if (kmin == kmax)
  {
    shift = ni*nj1*nk1 + ni1*nj*nk1;
    size = (jmax-jmin)*(imax-imin);
    PyObject* o = K_NUMPY::buildNumpyArray(size, 1, 1);
    E_Int* p = K_NUMPY::getNumpyPtrI(o);
    for (E_Int j = jmin-1; j < jmax-1; j++)
    {
      for (E_Int i = imin-1; i < imax-1; i++)
      {
        ind = shift + i + j*ni1 + (kmin-1)*ni1*nj1;
        p[ii] = ind; ii++;
      }
    }
    return o;
  }
  else
  {
    printf("Warning: range2PointList: requires a 2D range.\n");
    return NULL;
  }
}
