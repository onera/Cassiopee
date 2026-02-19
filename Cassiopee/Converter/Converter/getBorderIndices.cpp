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

// get indices from join border
# include "converter.h"

#define isPositive(a) (a < 0 ? 0 : 1)
#define signDir(a) (a < 0 ? -1 : 1)

using namespace K_FUNC;

//=============================================================================
/* Get indices from join border */
//=============================================================================
PyObject* K_CONVERTER::getJoinBorderIndices(PyObject* self, PyObject* args)
{
  IMPORTNUMPY;
  PyObject* arrayI;
  E_Int Wimin, Wimax, Wjmin, Wjmax, Wkmin, Wkmax;
  E_Int Im, Jm, Km;
  E_Int Direction, D, Dim, Dim1;
  E_Int Shift;
  if (!PYPARSETUPLE_(args, O_ IIII_ IIII_ IIII_ II_, 
                     &arrayI, &Dim1, &Im, &Jm, &Km,
                     &Wimin, &Wimax, &Wjmin, &Wjmax, &Wkmin, &Wkmax,
                     &Direction, &Dim, &D, &Shift)) return NULL;

  // Get array
  PyArrayObject* arI = (PyArrayObject*)arrayI;
  E_Int* arrayIndices = (E_Int*)PyArray_DATA(arI);

  E_Int im = Im;
  E_Int jm = Jm;
  E_Int km = Km;
  E_Int dim1 = Dim1;
  E_Int wimin= Wimin;
  E_Int wimax= Wimax;
  E_Int wjmin= Wjmin;
  E_Int wjmax= Wjmax;
  E_Int wkmin= Wkmin;
  E_Int wkmax= Wkmax;
  E_Int dim = Dim;
  E_Int direction = Direction;
  E_Int absdir = E_abs(direction);
  E_Int posdir = isPositive(direction);
  E_Int d = D;
  
  E_Int shiftDir = 0;
  if ( Shift != 0 ) shiftDir = -Shift*signDir(direction);
  // local values
  E_Int ind;
  E_Int im2,im2jm2;
  switch (dim)
  {
    case 3: // 3D Treatment
      im2 = im+2*d;
      im2jm2 = im2*(jm+2*d);
      // initialization of array of indices and build list
      if (absdir == 1)
      {
        for (E_Int k = 0; k < wkmax-wkmin+1; k++)
          for (E_Int j = 0; j < wjmax-wjmin+1; j++)
          {
            ind = posdir*(im-1)+d+shiftDir + (j+wjmin-1+d)*im2+(k+wkmin-1+d)*im2jm2;
            arrayIndices[k*dim1+j]=ind;
          }
      }
      else if (absdir == 2)
      {
        for (E_Int k=0; k<wkmax-wkmin+1;k++)
          for (E_Int i=0; i<wimax-wimin+1;i++)
          {
            ind = i+wimin-1+d  + (posdir*(jm-1)+d+shiftDir)*im2+(k+wkmin-1+d)*im2jm2;
            arrayIndices[k*dim1+i]=ind;
          }
      }
      else if (absdir == 3)
      {
        for (E_Int j=0; j<wjmax-wjmin+1;j++)
          for (E_Int i=0; i<wimax-wimin+1;i++)
          {
            ind = i+wimin-1+d  + (j+wjmin-1+d)*im2+(posdir*(km-1)+d+shiftDir)*im2jm2;
            arrayIndices[j*dim1+i]=ind;
          }
      }
      else
      {
        PyErr_SetString(PyExc_ValueError, "getJoinBorderIndices: invalid direction for join border.");
        return NULL;
      }
      break;
      
    case 2: // 2D Treatment
      im2 = im+2*d;
      if (absdir == 1)
      {
	for (E_Int j=0; j<wjmax-wjmin+1;j++)
	{
	  ind = posdir*(im-1)+d+shiftDir + (j+wjmin-1+d)*im2;
	  arrayIndices[j] = ind;
	}
      }
      else if (absdir == 2)
      {
	for (E_Int i=0; i<wimax-wimin+1;i++)
	{
	  ind = i+wimin-1+d  + (posdir*(jm-1)+d+shiftDir)*im2;
	  arrayIndices[i] = ind;
	}
      }
      break;
  
    case 1:
      PyErr_SetString(PyExc_ValueError, "getJoinBorderIndices: not implemented for 1D problems.");
      return NULL;
  
  }

  Py_INCREF(Py_None);
  return Py_None;
}
