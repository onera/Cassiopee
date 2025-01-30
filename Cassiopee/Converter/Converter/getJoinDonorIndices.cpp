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

// get a list of donor indices for a matching join
# include "converter.h"

#define signVar(a) (a < 0 ? -1 : 1)

using namespace std;
using namespace K_FUNC;

//=============================================================================
/* Get indices from join border */
//=============================================================================
PyObject* K_CONVERTER::getJoinDonorIndices(PyObject* self, PyObject* args)
{
  IMPORTNUMPY;
  PyObject* arrayI;
  PyObject* arrayBorderI;
  E_Int T1, T2, T3;
  E_Int Direction, Dim, DimI, DimBI;
  if (!PYPARSETUPLE_(args, OO_ IIII_ III_,
			              &arrayI, &arrayBorderI, &T1, &T2, &T3,
			              &Direction, &Dim, &DimI, &DimBI)) return NULL;

  // Get arrays
  PyArrayObject* ar = (PyArrayObject*)arrayI;
  E_Int* array = (E_Int*)PyArray_DATA(ar);
  PyArrayObject* arB = (PyArrayObject*)arrayBorderI;
  E_Int* arrayborder = (E_Int*)PyArray_DATA(arB);
  
  // Get integer values
  E_Int t1 = T1;
  E_Int t2 = T2;
  E_Int t3 = T3;
  E_Int dim = Dim;
  E_Int dim1 = DimI;
  E_Int dimb1 = DimBI;
  E_Int direction = Direction;

  // local values
  E_Int abst1 = E_abs(t1);
  E_Int abst2 = E_abs(t2);
  E_Int abst3 = E_abs(t3);
  E_Int absdir = E_abs(direction);
  vector<E_Int> listjoindonor;

  switch (dim)
  {
    case 3: // 3D Treatment
      if (absdir == 1)
      {
        if (t2 < 0)
        {
          array[abst2-1] = array[dim1 + abst2-1]-1;
          array[dim1 + abst2-1] = -1;
        }
        if (t3 < 0)
        {
          array[abst3-1] = array[dim1 + abst3-1]-1;
          array[dim1 + abst3-1] = -1;
        }
        if ((signVar(t3) == 1)&&(signVar(t2) == 1 ))
        {
          for (E_Int j=array[abst3-1]; j < array[dim1 + abst3-1]; j++)
            for (E_Int i=array[abst2-1]; i < array[dim1 + abst2-1]; i++)
            {
              if (abst3 > abst2)
                listjoindonor.push_back(arrayborder[j*dimb1+i]);
              else
                listjoindonor.push_back(arrayborder[i*dimb1+j]);
            }
        }
        else if ((signVar(t3) == 1)&&(signVar(t2) == -1))
        {
          for (E_Int j=array[abst3-1]; j < array[dim1 + abst3-1]; j++)
            for (E_Int i=array[abst2-1]; i > array[dim1 + abst2-1]; i--)
            {
              if (abst3 > abst2)
                listjoindonor.push_back(arrayborder[j*dimb1+i]);
              else
                listjoindonor.push_back(arrayborder[i*dimb1+j]);
            }
        }
        else if ((signVar(t3) == -1)&&(signVar(t2) == 1))
        {
          for (E_Int j=array[abst3-1]; j > array[dim1 + abst3-1]; j--)
            for (E_Int i=array[abst2-1]; i < array[dim1 + abst2-1]; i++)
            {
              if (abst3 > abst2)
                listjoindonor.push_back(arrayborder[j*dimb1+i]);
              else
                listjoindonor.push_back(arrayborder[i*dimb1+j]);
            }
        }
        else
        {
          for (E_Int j=array[abst3-1]; j > array[dim1 + abst3-1]; j--)
            for (E_Int i=array[abst2-1]; i > array[dim1 + abst2-1]; i--)
            {
              if (abst3 > abst2)
                listjoindonor.push_back(arrayborder[j*dimb1+i]);
              else
                listjoindonor.push_back(arrayborder[i*dimb1+j]);
            }
        }
      }
      else if (absdir == 2)
      {
        if (t1 < 0)
        {
          array[abst1-1]=array[dim1 + abst1-1]-1;
          array[dim1 + abst1-1]=-1;
        }
        if (t3 < 0)
        {
          array[abst3-1]=array[dim1 + abst3-1]-1;
          array[dim1 + abst3-1]=-1;
        }
        if ((signVar(t3) == 1)&&(signVar(t1) == 1 ))
        {
          for (E_Int j=array[abst3-1]; j < array[dim1 + abst3-1]; j++)
            for (E_Int i=array[abst1-1]; i < array[dim1 + abst1-1]; i++)
            {
              if (abst3 > abst1)
                listjoindonor.push_back(arrayborder[j*dimb1+i]);
              else
                listjoindonor.push_back(arrayborder[i*dimb1+j]);
            }
        }
        else if ((signVar(t3) == 1)&&(signVar(t1) == -1 ))
        {
          for (E_Int j=array[abst3-1]; j < array[dim1 + abst3-1]; j++)
            for (E_Int i=array[abst1-1]; i > array[dim1 + abst1-1]; i--)
            {
              if (abst3 > abst1)
                listjoindonor.push_back(arrayborder[j*dimb1+i]);
              else
                listjoindonor.push_back(arrayborder[i*dimb1+j]);
            }
        }
        else if ((signVar(t3) == -1)&&(signVar(t1) == 1 ))
        {
          for (E_Int j=array[abst3-1]; j > array[dim1 + abst3-1]; j--)
            for (E_Int i=array[abst1-1]; i < array[dim1 + abst1-1]; i++)
            {
              if (abst3 > abst1)
                listjoindonor.push_back(arrayborder[j*dimb1+i]);
              else
                listjoindonor.push_back(arrayborder[i*dimb1+j]);
            }
        }
        else
        {
          for (E_Int j=array[abst3-1]; j > array[dim1 + abst3-1]; j--)
            for (E_Int i=array[abst1-1]; i > array[dim1 + abst1-1]; i--)
            {
              if (abst3 > abst1)
                listjoindonor.push_back(arrayborder[j*dimb1+i]);
              else
                listjoindonor.push_back(arrayborder[i*dimb1+j]);
            }
        }
      }
      else if (absdir == 3)
      {
        if (t1 < 0)
        {
          array[abst1-1]=array[dim1 + abst1-1]-1;
          array[dim1 + abst1-1]=-1;
        }
        if (t2 < 0)
        {
          array[abst2-1]=array[dim1 + abst2-1]-1;
          array[dim1 + abst2-1]=-1;
        }
        if ((signVar(t2) == 1)&&(signVar(t1) == 1 ))
        {
          for (E_Int j=array[abst2-1]; j < array[dim1 + abst2-1]; j++)
            for (E_Int i=array[abst1-1]; i < array[dim1 + abst1-1]; i++)
            {
              if (abst2 > abst1)
                listjoindonor.push_back(arrayborder[j*dimb1+i]);
              else
                listjoindonor.push_back(arrayborder[i*dimb1+j]);
            }
        }
        else if ((signVar(t2) == 1)&&(signVar(t1) == -1 ))
        {
          for (E_Int j=array[abst2-1]; j < array[dim1 + abst2-1]; j++)
            for (E_Int i=array[abst1-1]; i > array[dim1 + abst1-1]; i--)
            {
              if (abst2 > abst1)
                listjoindonor.push_back(arrayborder[j*dimb1+i]);
              else
                listjoindonor.push_back(arrayborder[i*dimb1+j]);
            }
        }
        else if ((signVar(t2) == -1)&&(signVar(t1) == 1 ))
        {
          for (E_Int j=array[abst2-1]; j > array[dim1 + abst2-1]; j--)
            for (E_Int i=array[abst1-1]; i < array[dim1 + abst1-1];i++)
            {
              if (abst2 > abst1)
                listjoindonor.push_back(arrayborder[j*dimb1+i]);
              else
                listjoindonor.push_back(arrayborder[i*dimb1+j]);
            }
        }
        else
        {
          for (E_Int j=array[abst2-1]; j > array[dim1 + abst2-1]; j--)
            for (E_Int i=array[abst1-1]; i > array[dim1 + abst1-1]; i--)
            {
              if (abst2 > abst1)
                listjoindonor.push_back(arrayborder[j*dimb1+i]);
              else
                listjoindonor.push_back(arrayborder[i*dimb1+j]);
            }
        }
      }
      else
      {
        PyErr_SetString(PyExc_ValueError, "getJoinBorderIndices: unvalid direction for join border.");
        return NULL;
      }
      break;
      
    case 2: // 2D Treatment
      if (absdir == 1)
      {
        if (t2 < 0)
        {
         array[abst2-1]=array[dim1 + abst2-1]-1;
         array[dim1 + abst2-1]=-1;
         for (E_Int i=array[abst2-1]; i > array[dim1 + abst2-1]; i--)
           listjoindonor.push_back(arrayborder[i]);
        }
        else
        {
          for (E_Int i=array[abst2-1]; i < array[dim1 + abst2-1]; i++)
            listjoindonor.push_back(arrayborder[i]);
        }
      }
      else if (absdir == 2)
      {
        if (t1 <0)
        {
          array[abst1-1]=array[dim1 + abst1-1]-1;
          array[dim1 + abst1-1]=-1;
          for (E_Int i=array[abst1-1]; i > array[dim1 + abst1-1]; i--)
            listjoindonor.push_back(arrayborder[i]);
        }
        else
        {
          for (E_Int i=array[abst1-1]; i < array[dim1 + abst1-1]; i++)
            listjoindonor.push_back(arrayborder[i]);
        }
      }
      break;
  
    case 1:
      PyErr_SetString(PyExc_ValueError, "getJoinBorderIndices: not implemented for 1D problems.");
      return NULL;
  }

  // Convert vector<E_Int> listjoindonor in PyObject tpljoindonor
  E_Int size = listjoindonor.size();
  PyObject* tpljoindonor;
  PyObject* l = PyList_New(0);
  for (E_Int i = 0; i < size; i++)
  {
    tpljoindonor = Py_BuildValue(I_, listjoindonor[i]);
    PyList_Append(l, tpljoindonor);
    Py_DECREF(tpljoindonor);
  }

  return l;
}
