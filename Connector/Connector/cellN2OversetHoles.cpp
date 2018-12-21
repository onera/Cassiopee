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
# include "connector.h"

// ============================================================================
/* put blanked cells in array PointList for structured grids */
// ============================================================================
PyObject* K_CONNECTOR::cellN2OversetHolesStruct(PyObject* self, PyObject* args)
{
  IMPORTNUMPY;
  PyObject* PointList;
  PyObject* CellN;
  E_Int Im, Jm, CellDim;
  E_Int SizeCellN;
  if (!PYPARSETUPLEI(args,
                    "OOllll", "OOiiii",
                    &PointList, &CellN, &Im, &Jm, &CellDim, &SizeCellN))
  {
      return NULL;
  }
  // Get E_Int values
  E_Int im            = Im;
  E_Int jm            = Jm;
  E_Int cellDim       = CellDim;
  E_Int sizeCellN     = SizeCellN;

  // Get CellN array
  PyArrayObject* cellNatureField = (PyArrayObject*)CellN;
  double* cellN = (double*)PyArray_DATA(cellNatureField);

  // Get PointList array
  PyArrayObject* pointListArray = (PyArrayObject*)PointList;
  int* pointList = (int*)PyArray_DATA(pointListArray);

  // Local values
  E_Int nb = 0;
  E_Int i  = 0;
  E_Int j  = 0;
  E_Int k  = 0;
  E_Int imjm;
  E_Float eps = K_CONST::E_GEOM_CUTOFF;

  // fill pointList array with (i,j,k)-indices from blanked cells
  nb = 0;
  imjm = im*jm;

  switch (cellDim)
  {
    case 3: // 3D
      for (E_Int ind = 0;ind < sizeCellN; ind++)
      {
        k = ind/imjm;
        j = (ind - k*imjm)/im;
        i = ind - k*imjm - j*im;
        if (K_FUNC::fEqualZero(cellN[ind], eps) == true)
        {
          pointList[nb*cellDim] = i+1;
          pointList[nb*cellDim+1] = j+1;
          pointList[nb*cellDim+2] = k+1;
          nb++;
        }
      }
      break;
    case 2: // 2D
      for (E_Int ind = 0;ind < sizeCellN; ind++)
      {
        j = (ind - k*imjm)/im;
        i = ind - k*imjm - j*im;
        if (K_FUNC::fEqualZero(cellN[ind], eps) == true)
        {
          pointList[nb*cellDim] = i+1;
          pointList[nb*cellDim+1] = j+1;
          nb++;
        }
      }
      break;
    case 1: // 1D
      for (E_Int ind = 0;ind < sizeCellN; ind++)
        if (K_FUNC::fEqualZero(cellN[ind], eps) == true)
        {
          pointList[nb] = ind+1;
          nb++;
        }        
      break;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

// ============================================================================
/* put blanked cells in array PointList for unstructured grids */
// ============================================================================
PyObject* K_CONNECTOR::cellN2OversetHolesUnStruct(PyObject* self, 
						  PyObject* args)
{
  IMPORTNUMPY;
  PyObject* PointList;
  PyObject* CellN;
  E_Int SizeCellN;
  if (!PYPARSETUPLEI(args,
                    "OOl", "OOi",
                    &PointList, &CellN, &SizeCellN))
  {
      return NULL;
  }
  
  // Get E_Int values
  E_Int sizeCellN = SizeCellN;

  // Get CellN array
  PyArrayObject* cellNatureField = (PyArrayObject*)CellN;
  double* cellN = (double*)PyArray_DATA(cellNatureField);

  // Get PointList array
  PyArrayObject* pointListArray = (PyArrayObject*)PointList;
  int* pointList = (int*)PyArray_DATA(pointListArray);

  // Local values
  E_Int nb = 0;

  E_Float eps = K_CONST::E_GEOM_CUTOFF;

  // fill pointList array with (i,j,k)-indices from blanked cells
  nb = 0;
  
  for (E_Int ind = 0;ind < sizeCellN; ind++)
  {
    if (K_FUNC::fEqualZero(cellN[ind], eps) == true)
    {
      pointList[nb] = ind+1;
      nb++;
    }
  }
  Py_INCREF(Py_None);
  return Py_None;
}
