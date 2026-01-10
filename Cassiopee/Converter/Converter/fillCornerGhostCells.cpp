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

// fill join border values of an array
# include "converter.h"

using namespace K_FUNC;

//=============================================================================
/* Fill join border values of an array with values of donor block */
//=============================================================================
// Current join has a local frame (I,J) with increments incrRcvI and incrRcvJ 
// to shift from a cell to the next one in direction I or J
// Opposite join has increments incrDonorI and incrDonorJ with correspond to 
// the same local frame (I,J):
// if cell with address RcvCell has for donor opposit cell with address 
// DonorCell, cell with address (RcvCell + incrRcvI) has for donor opposite 
// cell with address (DonorCell + incrDonorI)
//=============================================================================
PyObject* K_CONVERTER::fillCornerGhostCells(PyObject* self, PyObject* args)
{
  IMPORTNUMPY;
  PyObject* arrayR; PyObject* arrayD;
  PyObject* arrayIR; PyObject* listID;
  char* loc;
  E_Int IncrRcvI, IncrRcvJ, IncrDonorI, IncrDonorJ;
  E_Int arrayRdim1, arrayRdim2;
  E_Int Direction, DirDonor, IncrR, IncrD, DimZone, D;
  E_Int Passage;
  if (!PYPARSETUPLE_(args, OOOO_ S_ IIII_ IIII_ IIII_ I_,
                     &arrayR, &arrayD,
                     &arrayIR, &listID, &loc, 
                     &arrayRdim1, &arrayRdim2, 
                     &IncrRcvI, &IncrRcvJ, &IncrDonorI, &IncrDonorJ,
                     &Direction, &DirDonor, &IncrR, &IncrD, &DimZone, &D, 
                     &Passage)) return NULL;

  // Get E_Int values
  E_Int dim = DimZone;
  E_Int dim1 = arrayRdim1;
  E_Int dim2 = arrayRdim2;
  E_Int direction  = Direction;
  E_Int dirdonor  = DirDonor;
  E_Int d = D;
  E_Int passage = Passage;
  E_Int incrrecv = IncrR;    // increment between to cells in the join border direction for reciever block
  E_Int incrdonor = IncrD;   // increment between to cells in the join border direction for donor block
  // join windows of reciever and donor blocks have a same local (I,J) frame 
  // with following corresponding increment :
  E_Int incrRcvI   = IncrRcvI;
  E_Int incrRcvJ   = IncrRcvJ;
  E_Int incrDonorI = IncrDonorI;
  E_Int incrDonorJ = IncrDonorJ;

  // check if listID is a list
  if (PyList_Check(listID) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "fillCornerGostCells: donor indices argument must be a list.");
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

  // dataR: data filled for the join (reciever block)
  // dataD: data from donor block
  // arrayBorder: array of indices of join border (receiver block)
  PyArrayObject* arR = (PyArrayObject*)arrayR;
  PyArrayObject* arD = (PyArrayObject*)arrayD;
  PyArrayObject* arI = (PyArrayObject*)arrayIR;
  E_Float* dataR = (E_Float*)PyArray_DATA(arR);
  E_Float* dataD = (E_Float*)PyArray_DATA(arD);
  E_Int* arrayBorder = (E_Int*)PyArray_DATA(arI);

  // local values
  E_Int indadj;                 // cell index of reciever block cell adjacent to the join border
  E_Int indadjdonor;            // cell index of donor block cell
  E_Int inddonor, indghost;
  E_Int shift_dir = direction/E_abs(direction)*incrrecv;
  E_Int shift_dirdonor = dirdonor/E_abs(dirdonor)*incrdonor;
  E_Int indadjmin,indadjmax;
  E_Int indadjdonormin,indadjdonormax;
  E_Int indRghostmin, indRghostmax;
  E_Int indDghostmin, indDghostmax;
  E_Int indghostmin, indghostmax;
  E_Int indDmin, indDmax;
  E_Int incrR2, indRg1, indRg2;
  E_Int incrD2, indDg1, indDg2;
  E_Int localisation = 0;
  if (!strcmp(loc, "Vertex")) localisation = 1;

  switch (dim)
  {
    case 3:
      if (passage == 1)
      {
        for (E_Int i = 0; i < dim1; i++)
        {
          indadjmin = arrayBorder[i];
          indadjmax = arrayBorder[(dim2-1)*dim1+i];
          indadjdonormin = arrayDonor[i];
          indadjdonormax = arrayDonor[(dim2-1)*dim1+i];
          for (E_Int index = 0; index < d; index++)
          {
            indghostmin = indadjmin+shift_dir*(index+1);
            indghostmax = indadjmax+shift_dir*(index+1);
            indDmin = indadjdonormin-shift_dirdonor*(index+localisation);
            indDmax = indadjdonormax-shift_dirdonor*(index+localisation);
            for (E_Int indcorner = 0; indcorner < d; indcorner++)
            {
              // index of corner ghost cells
              indRghostmin = indghostmin-(indcorner+1)*incrRcvJ;
              indRghostmax = indghostmax+(indcorner+1)*incrRcvJ;
              // index of donor cells for corner ghost cells
              indDghostmin = indDmin-(indcorner+1)*incrDonorJ;
              indDghostmax = indDmax+(indcorner+1)*incrDonorJ;
              dataR[indRghostmin] = dataD[indDghostmin];
              dataR[indRghostmax] = dataD[indDghostmax];
            }
          }
        }
      }
      for (E_Int j = 0; j < dim2; j++)
      {
	indadjmin = arrayBorder[j*dim1];
	indadjmax = arrayBorder[j*dim1+dim1-1];
        indadjdonormin = arrayDonor[j*dim1];
        indadjdonormax = arrayDonor[j*dim1+dim1-1];
	for (E_Int index = 0; index < d; index++)
	{
	  indghostmin = indadjmin+shift_dir*(index+1);
	  indghostmax = indadjmax+shift_dir*(index+1);
          indDmin = indadjdonormin-shift_dirdonor*(index+localisation);
          indDmax = indadjdonormax-shift_dirdonor*(index+localisation);
	  for (E_Int indc1 = 0; indc1 < d; indc1++)
          {
            indRg1 = indghostmin-(indc1+1)*incrRcvI;
            indDg1 = indDmin - (indc1+1)*incrDonorI;
             if (passage == 1)
             {
               dataR[indRg1] = dataD[indDg1];
             }
             else if (passage == 2)
             {
               if (j == 0)
               {
                 incrR2 = -incrRcvJ; incrD2 = -incrDonorJ;                 
               }
               else if (j == dim2-1)
               {
                 incrR2 = incrRcvJ;
                 incrD2 = incrDonorJ;
               }
               if (j == 0 || j == dim2-1)
               {
                 for (E_Int indc2 = 0; indc2 < d; indc2++)
                 {
                   dataR[indRg1+(indc2+1)*incrR2] = dataD[indDg1+(indc2+1)*incrD2];
                 }
               }
             }
             indRg2 = indghostmax+(indc1+1)*incrRcvI;
             indDg2 = indDmax+(indc1+1)*incrDonorI;
             if (passage == 1)
             {
               dataR[indRg2]= dataD[indDg2];                   
             }
             else if (passage == 2)
             {
               if (j == 0)
               {
                 incrR2=-incrRcvJ;incrD2=-incrDonorJ;
               }
               else if (j == dim2-1)
               {
                 incrR2=+incrRcvJ;incrD2=+incrDonorJ;
               }
               if (j == 0 || j == dim2-1)
               {
                 for (E_Int indc2 = 0; indc2 < d; indc2++)
                 {
                   dataR[indRg2+(indc2+1)*incrR2] = dataD[indDg2+(indc2+1)*incrD2];
                 }
               }
             }
          }
	}
      }
      break;

    case 2:
      for (E_Int i = 0; i < dim1; i++)
      {
	indadj = arrayBorder[i];
        indadjdonor = arrayDonor[i];
	if (i == 0 || i == dim1-1)
	{
	  incrRcvI = -incrRcvI;
	  incrDonorI = -incrDonorI;
	  for (E_Int index = 0; index < d; index++)
	  {
	    indghost = indadj+shift_dir*(index+1);
            inddonor = indadjdonor-shift_dirdonor*(index+localisation);
	    for (E_Int index2 = 0; index2 < d; index2++)
	      dataR[indghost+(index2+1)*incrRcvI] = dataD[inddonor+(index2+1)*incrDonorI];
	  }
	}
      }
      break;
  }
  delete [] arrayDonor;
  Py_INCREF(Py_None);
  return Py_None;
}
