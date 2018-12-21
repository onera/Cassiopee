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

// fill join border values of an array

# include "converter.h"
using namespace K_FUNC;

//=============================================================================
/* Fill join border values of an array with values of donor block */
//=============================================================================
PyObject* K_CONVERTER::fillJoin(PyObject* self, PyObject* args)
{
  IMPORTNUMPY;
  PyObject* arrayR; PyObject* arrayD; PyObject* arrayIR;
  PyObject* listID;
  char* loc;
  E_Int arrayRdim1, arrayRdim2;
  E_Int Direction, DirDonor, IncrR, IncrD, DimZone, D;
  E_Int Im=0;
  E_Int ImDonor=0;
  if (!PYPARSETUPLEI(args, "OOOOsllllllllll", "OOOOsiiiiiiiiii",
                     &arrayR, &arrayD, 
                     &arrayIR, &listID, &loc, 
                     &arrayRdim1, &arrayRdim2,
                     &Direction, &DirDonor, &IncrR, &IncrD, &DimZone, &D , 
                     &Im, &ImDonor)) return NULL;

  // Get E_Int values
  E_Int dim = DimZone;
  E_Int dim1 = arrayRdim1;
  E_Int dim2 = arrayRdim2;
  E_Int im = Im;
  E_Int imdonor = ImDonor;
  E_Int direction = Direction;
  E_Int dirdonor = DirDonor;
  E_Int d = D;
  E_Int incrrecv = IncrR;    // increment between to cells in the join border direction for reciever block
  E_Int incrdonor = IncrD;   // increment between to cells in the join border direction for donor block

  // check if listID is a list
  if (PyList_Check (listID) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "fillJoin: donor indices argument must be a list.");
    return NULL;
  }
  E_Int listIDSize = PyList_Size(listID);

  // Convert Python listID in array arrayDonor: 
  // array of indices of donor block cells
  E_Int* arrayDonor = new E_Int[listIDSize]; 
  for (E_Int i = 0; i < listIDSize; i++)
  {
    PyObject* tpl = PyList_GetItem(listID, i);
    arrayDonor[i] = PyInt_AsLong(tpl);
  }

  // dataR: data filled for the join (receiver block)
  // dataD: data from donor block
  // arrayBorder: array of indices of join border (reciever block)
  PyArrayObject* arR = (PyArrayObject*)arrayR;
  PyArrayObject* arD = (PyArrayObject*)arrayD;
  PyArrayObject* arI = (PyArrayObject*)arrayIR;
  double* dataR = (double*)PyArray_DATA(arR);
  double* dataD = (double*)PyArray_DATA(arD);
  int* arrayBorder = (int*)PyArray_DATA(arI);

  // local values
  E_Int indadj;               // cell index of reciever block cell adjacent to the join border
  E_Int indadjdonor;          // cell index of donor block cell
  E_Int inddonor, indghost;
  E_Int shift_dir      = direction/E_abs(direction)*incrrecv;
  E_Int shift_dirdonor = dirdonor/E_abs(dirdonor)*incrdonor;
  E_Int ind, ghost;
  E_Int localisation = 0;
  if (!strcmp(loc, "Vertex")) localisation = 1;

  switch (dim)
  {
    case 3:
      ind = 0;
      for (E_Int j = 0; j < dim2; j++)
        for (E_Int i = 0; i < dim1; i++)
        {
         indadj = arrayBorder[j*dim1+i];
         indadjdonor = arrayDonor[ind];
         for (E_Int index = 0; index < d; index++)
         {
          indghost = indadj+shift_dir*(index+1);
          inddonor = indadjdonor-shift_dirdonor*(index+localisation);
          dataR[indghost] = dataD[inddonor];
        }
        ind++;
      }
      break;

    case 2:
      for (E_Int i = 0; i < dim1; i++)
      {
        indadj = arrayBorder[i];
        indadjdonor = arrayDonor[i];
        for (E_Int index = 0; index < d; index++)
        {
          indghost = indadj+shift_dir*(index+1);
          inddonor = indadjdonor-shift_dirdonor*(index+localisation);
          dataR[indghost] = dataD[inddonor];
        }
      }
      break;
  
    case 1:
      if (direction == -1) // min
      {
        for (E_Int i = -d; i < 0; i++)
        {
          ghost = i+d;
          if (dirdonor == -1) inddonor = -i+localisation;
          else inddonor = imdonor+i-localisation;
          dataR[ghost] = dataD[inddonor]; 
        }
      }       
      else // max
      {
        for (E_Int i = im; i < im+d; i++)
        {
          ghost = i+d;
          if (dirdonor == -1) inddonor = i-im+localisation;
          else inddonor = imdonor - (i-im) -1 - localisation;
          dataR[ghost] = dataD[inddonor]; 
        }
      }
      break;
  }

  delete [] arrayDonor;

  Py_INCREF(Py_None);
  return Py_None;
}
