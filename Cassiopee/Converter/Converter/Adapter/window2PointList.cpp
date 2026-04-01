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
# include "converter.h"

using namespace K_FLD;

//=============================================================================
// Convert the range of a structured block into a face indices numpy. 
// Face indices start at 0 and contain all i faces, followed by all j faces and
// k faces.
//=============================================================================
PyObject* K_CONVERTER::window2FacePointList(PyObject* self, PyObject* args)
{
  E_Int imin, imax, jmin, jmax, kmin, kmax, ni, nj, nk;
  if (!PYPARSETUPLE_(args, IIII_ IIII_ I_, &imin, &imax, &jmin, &jmax,
                     &kmin, &kmax, &ni, &nj, &nk)) return NULL;
  
  if (kmin == kmax && nk == 1) kmax = 2; // 2D
  
  E_Int shift, size;
  E_Int ni1 = K_FUNC::E_max(ni - 1, 1);
  E_Int nj1 = K_FUNC::E_max(nj - 1, 1);
  E_Int nk1 = K_FUNC::E_max(nk - 1, 1);
  E_Int imin1 = imin - 1, jmin1 = jmin - 1, kmin1 = kmin - 1;
  E_Int imax1 = imax - 1, jmax1 = jmax - 1, kmax1 = kmax - 1;
  E_Int di = imax - imin;
  E_Int dj = jmax - jmin;
  E_Int dk = kmax - kmin;
  
  if (imin == imax)  // window in i
  {
    size = dj*dk;
    PyObject* o = K_NUMPY::buildNumpyArray(size, 1, 1);
    E_Int* p = K_NUMPY::getNumpyPtrI(o);
    #pragma omp parallel
    {
      E_Int ind;
      #pragma omp for collapse(2)
      for (E_Int k = kmin1; k < kmax1; k++)
      for (E_Int j = jmin1; j < jmax1; j++)
      {
        ind = (k - kmin + 1)*dj + (j - jmin + 1);
        p[ind] = imin - 1 + j*ni + k*ni*nj1;
      }
    }
    return o;
  }
  else if (jmin == jmax)  // window in j
  {
    shift = ni*nj1*nk1;
    size = di*dk;
    PyObject* o = K_NUMPY::buildNumpyArray(size, 1, 1);
    E_Int* p = K_NUMPY::getNumpyPtrI(o);
    #pragma omp parallel
    {
      E_Int ind;
      #pragma omp for collapse(2)
      for (E_Int k = kmin1; k < kmax1; k++)
      for (E_Int i = imin1; i < imax1; i++)
      {
        ind = (k - kmin + 1)*di + (i - imin + 1);
        p[ind] = shift + i + (jmin - 1)*ni1 + k*ni1*nj;
      }
    }
    return o;
  }
  else if (kmin == kmax)  // window in k
  {
    shift = ni*nj1*nk1 + ni1*nj*nk1;
    size = di*dj;
    PyObject* o = K_NUMPY::buildNumpyArray(size, 1, 1);
    E_Int* p = K_NUMPY::getNumpyPtrI(o);
    #pragma omp parallel
    {
      E_Int ind;
      #pragma omp for collapse(2)
      for (E_Int j = jmin1; j < jmax1; j++)
      for (E_Int i = imin1; i < imax1; i++)
      {
        ind = (j - jmin + 1)*di + (i - imin + 1);
        p[ind] = shift + i + j*ni1 + (kmin - 1)*ni1*nj1;
      }
    }
    return o;
  }
  else
  {
    PyErr_SetString(PyExc_ValueError,
                    "window2FacePointList: requires a 2D range.\n");
    return NULL;
  }
}

//=============================================================================
// Convert the range of a structured block into a vertex indices numpy. 
// Vertex indices start at 0 and contain all i vertices, followed by all
// j vertices and k vertices.
//=============================================================================
PyObject* K_CONVERTER::window2VertexPointList(PyObject* self, PyObject* args)
{
  E_Int imin, imax, jmin, jmax, kmin, kmax, ni, nj, nk;
  if (!PYPARSETUPLE_(args, IIII_ IIII_ I_, &imin, &imax, &jmin, &jmax,
                     &kmin, &kmax, &ni, &nj, &nk)) return NULL;

  E_Int size;
  E_Int di = imax - imin + 1;
  E_Int dj = jmax - jmin + 1;
  E_Int dk = kmax - kmin + 1;
  E_Int imin1 = imin - 1, jmin1 = jmin - 1, kmin1 = kmin - 1;

  if (imin == imax)  // window in i
  {
    size = dj*dk;
    PyObject* o = K_NUMPY::buildNumpyArray(size, 1, 1);
    E_Int* p = K_NUMPY::getNumpyPtrI(o);
    #pragma omp parallel
    {
        #pragma omp for collapse(2)
        for (E_Int k = kmin1; k < kmax; k++)
        for (E_Int j = jmin1; j < jmax; j++)
        {
          E_Int ind = (k - kmin1)*dj + (j - jmin1);
          p[ind] = imin1 + j*ni + k*ni*nj;
        }
    }
    return o;
  }
  else if (jmin == jmax)
  {
    size = di*dk;
    PyObject* o = K_NUMPY::buildNumpyArray(size, 1, 1);
    E_Int* p = K_NUMPY::getNumpyPtrI(o);

    #pragma omp parallel
    {
      #pragma omp for collapse(2)
      for (E_Int k = kmin1; k < kmax; k++)
      for (E_Int i = imin1; i < imax; i++)
      {
        E_Int ind = (k - kmin1)*di + (i - imin1);
        p[ind] = i + jmin1*ni + k*ni*nj;
      }
    }
    return o;
  }
  else if (kmin == kmax)
  {
    size = di*dj;
    PyObject* o = K_NUMPY::buildNumpyArray(size, 1, 1);
    E_Int* p = K_NUMPY::getNumpyPtrI(o);

    #pragma omp parallel
    {
      #pragma omp for collapse(2)
      for (E_Int j = jmin1; j < jmax; j++)
      for (E_Int i = imin1; i < imax; i++)
      {
        E_Int ind = (j - jmin1)*di + (i - imin1);
        p[ind] = i + j*ni + kmin1*ni*nj;
      }
    }
    return o;
  }
  else
  {
    PyErr_SetString(PyExc_ValueError,
                    "window2VertexPointList: requires a 2D range.");
    return NULL;
  }
}