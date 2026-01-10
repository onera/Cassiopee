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

// fill ghost cell values for nearmatch connectivities

# include "converter.h"
using namespace K_FUNC;

//=============================================================================
/* Coordinates and fields located at nodes */
//=============================================================================
PyObject* K_CONVERTER::fillJoinNMNodes(PyObject* self, PyObject* args)
{
  IMPORTNUMPY;
  PyObject* arrayR; PyObject* arrayD; PyObject* arrayIR;
  PyObject* listID;
  E_Int arrayRdim1, arrayRdim2;
  E_Int Direction, DirDonor, IncrR, IncrD, DimZone, D;
  E_Int Im=0;
  E_Int ImDonor=0;
  E_Int ShiftDir1, ShiftDir2, IsFine;
  if (!PYPARSETUPLE_(args, OOOO_ IIII_ IIII_ IIII_ I_, 
                     &arrayR, &arrayD, &arrayIR, &listID,
                     &arrayRdim1, &arrayRdim2, &Direction, &DirDonor, 
                     &IncrR, &IncrD, &DimZone, &D, &Im, &ImDonor,
                     &ShiftDir1, &ShiftDir2, &IsFine)) return NULL;

  // Get E_Int values
  E_Int dim = DimZone;
  E_Int dim1 = arrayRdim1;
  E_Int dim2 = arrayRdim2;
  E_Int direction = Direction;
  E_Int dirdonor = DirDonor;
  E_Int d = D;
  E_Int incrrecv = IncrR;    
  E_Int incrdonor = IncrD; 
  E_Int shiftDir1 = ShiftDir1;
  E_Int shiftDir2 = ShiftDir2;
  E_Int isFine=IsFine;// 1: rcv is finest, -1: rcv is coarsest
  // check if listID is a list
  if (PyList_Check (listID) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "fillJoinNM: donor indices argument must be a list.");
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
  // dataR: data filled for the join (receptor block)
  // dataD: data from donor block
  // arrayBorder: array of indices of join border (receptor block)
  PyArrayObject* arR = (PyArrayObject*)arrayR;
  PyArrayObject* arD = (PyArrayObject*)arrayD;
  PyArrayObject* arI = (PyArrayObject*)arrayIR;
  E_Float* dataR = (E_Float*)PyArray_DATA(arR);
  E_Float* dataD = (E_Float*)PyArray_DATA(arD);
  E_Int* arrayBorder = (E_Int*)PyArray_DATA(arI);

  // local values
  E_Int indadj, indadjdonor, indadjprev, indadjnext;             
  E_Int inddonor, indghost, indghostprev, indghostnext;
  E_Int shift_dir       = direction/E_abs(direction)*incrrecv;
  E_Int shift_dirdonor  = dirdonor/E_abs(dirdonor)*incrdonor;
  E_Int ind, i, j, iprev, inext;
  E_Float invShiftDir1 = 1./shiftDir1;
  E_Float invShiftDir2 = 1./shiftDir2;
  E_Int dim1D = (dim1-1)*shiftDir1+1;

  if (isFine == 1) // finest zone
  {
    switch (dim)
    {
      case 3:
        // First: fill oneovern ghost cells
        ind = 0;
        for (j=0; j < dim2; j+= shiftDir2)
          for (i=0; i< dim1; i+= shiftDir1)
          {
            indadj = arrayBorder[j*dim1+i];
            indadjdonor = arrayDonor[ind];
            for (E_Int index = 0; index < d; index++)
            {
              indghost = indadj+shift_dir*(index+1);
              inddonor = indadjdonor-shift_dirdonor*(index+1);
              dataR[indghost] = dataD[inddonor];
            }
            ind++;
          }

        if ( shiftDir1 > 1 ) 
        {   
          // Then fill points between oneovern pts
          for (j=0; j < dim2; j++)
            for (i=0; i< dim1-shiftDir1; i+= shiftDir1)
            {
              indadjprev = arrayBorder[j*dim1+i];
              indadjnext = arrayBorder[j*dim1+i+shiftDir1];
              for (E_Int index = 0; index < d; index++)
              {
                indghostprev = indadjprev+ shift_dir*(index+1);
                indghostnext = indadjnext+ shift_dir*(index+1);
                  
                for (E_Int il = i+1; il < i+shiftDir1; il++)
                {
                  indadj = arrayBorder[j*dim1+il];
                  indghost = indadj+shift_dir*(index+1);
                  dataR[indghost]=dataR[indghostprev]+(il-i)*invShiftDir1*(dataR[indghostnext]-dataR[indghostprev]);
                }
              }            
            }
        }      
        if ( shiftDir2 > 1 ) 
        {
          // Then fill points between oneovern pts
          for (i=0; i< dim1; i++)
            for (j=0; j < dim2-shiftDir2; j+= shiftDir2)
            {
              indadjprev = arrayBorder[j*dim1+i];
              indadjnext = arrayBorder[(j+shiftDir2)*dim1+i];
              for (E_Int index = 0; index < d; index++)
              {
                indghostprev = indadjprev+ shift_dir*(index+1);
                indghostnext = indadjnext+ shift_dir*(index+1);
                
                for (E_Int jl = j+1; jl < j+shiftDir2; jl++)
                {
                  indadj = arrayBorder[jl*dim1+i];
                  indghost = indadj+shift_dir*(index+1);
                  dataR[indghost]=dataR[indghostprev]+(jl-j)*invShiftDir2*(dataR[indghostnext]-dataR[indghostprev]);
                }
              }            
            }
        }

        break;
        
      case 2:
        ind = 0;
        // First: fill oneovern ghost cells
        for (i=0; i < dim1; i += shiftDir1)
        {
          indadj = arrayBorder[i];
          indadjdonor = arrayDonor[ind];
          for (E_Int index = 0; index < d; index++)
          {
            indghost = indadj+shift_dir*(index+1);
            inddonor = indadjdonor-shift_dirdonor*(index+1);
            dataR[indghost] = dataD[inddonor];
          }
          ind++;// increment donneur
        }

        // Then fill points between oneovern pts
        for (inext = shiftDir1; inext < dim1; inext+= shiftDir1)
        {
          iprev = inext-shiftDir1;
          indadjprev = arrayBorder[iprev];
          indadjnext = arrayBorder[inext];
          for (i = iprev+1; i < inext; i++)
          {
            indadj = arrayBorder[i];       
            for (E_Int index = 0; index < d; index++)
            {
              indghost = indadj+shift_dir*(index+1);
              indghostprev = indadjprev+shift_dir*(index+1);
              indghostnext = indadjnext+shift_dir*(index+1);
              dataR[indghost]=dataR[indghostprev]+
                (i-iprev)*invShiftDir1*(dataR[indghostnext]-dataR[indghostprev]);
              
            }          
          }
        }        
        break;
    }
  }
  else if (isFine == -1) //coarsest
  {
    switch (dim)
    {
      case 3:
        for (j=0; j < dim2; j++)
        {
          ind=shiftDir2*dim1D*j;
          for (i=0; i< dim1; i++)
          {
            indadj = arrayBorder[j*dim1+i];
            indadjdonor = arrayDonor[ind];
            for (E_Int index = 0; index < d; index++)
            {
              indghost = indadj+shift_dir*(index+1);
              inddonor = indadjdonor-shift_dirdonor*(index+1);
              dataR[indghost] = dataD[inddonor];
            }
            ind+=shiftDir1;
          }      
        }
        break;
        
      case 2:
        ind = 0;
        for (i=0; i < dim1; i++)
        {
          indadj = arrayBorder[i];
          indadjdonor = arrayDonor[ind];
          for (E_Int index = 0; index < d; index++)
          {
            indghost = indadj+shift_dir*(index+1);
            inddonor = indadjdonor-shift_dirdonor*(index+1);
            dataR[indghost] = dataD[inddonor];
          }
          ind+=shiftDir1;// increment donneur
        }        
      break;
    }
  }

  delete [] arrayDonor;
  Py_INCREF(Py_None);
  return Py_None;
}
//=============================================================================
/* Fields located at centers */
//=============================================================================
PyObject* K_CONVERTER::fillJoinNMCenters(PyObject* self, PyObject* args)
{
  IMPORTNUMPY;
  PyObject* arrayR; PyObject* arrayD; PyObject* arrayIR;
  PyObject* listID;
  E_Int arrayRdim1, arrayRdim2;
  E_Int Direction, DirDonor, IncrR, IncrD, DimZone, D;
  E_Int Im=0;
  E_Int ImDonor=0;
  E_Int ShiftDir1, ShiftDir2, IsFine;
  if (!PYPARSETUPLE_(args, OOOO_ IIII_ IIII_ IIII_ I_, 
                        &arrayR, &arrayD, &arrayIR, &listID,
                        &arrayRdim1, &arrayRdim2, &Direction, &DirDonor, 
                        &IncrR, &IncrD, &DimZone, &D, &Im, &ImDonor,
                        &ShiftDir1, &ShiftDir2, &IsFine)) return NULL;
    
  // Get E_Int values
  E_Int dim = DimZone;
  E_Int dim1 = arrayRdim1;
  E_Int dim2 = arrayRdim2;
  E_Int direction = Direction;
  E_Int dirdonor = DirDonor;
  E_Int d = D;
  E_Int incrrecv = IncrR;    
  E_Int incrdonor = IncrD; 
  E_Int shiftDir1 = ShiftDir1; 
  E_Int shiftDir2 = ShiftDir2;
  E_Int isFine=IsFine;// 1: rcv is finest, -1: rcv is coarsest
  // check if listID is a list
  if (PyList_Check (listID) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "fillJoinNM: donor indices argument must be a list.");
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
  // dataR: data filled for the join (receptor block)
  // dataD: data from donor block
  // arrayBorder: array of indices of join border (receptor block)
  PyArrayObject* arR = (PyArrayObject*)arrayR;
  PyArrayObject* arD = (PyArrayObject*)arrayD;
  PyArrayObject* arI = (PyArrayObject*)arrayIR;
  E_Float* dataR = (E_Float*)PyArray_DATA(arR);
  E_Float* dataD = (E_Float*)PyArray_DATA(arD);
  E_Int* arrayBorder = (E_Int*)PyArray_DATA(arI);

  // local values
  E_Int indadj, indadjdonor;             
  E_Int inddonor, indghost;
  E_Int shift_dir       = direction/E_abs(direction)*incrrecv;
  E_Int shift_dirdonor  = dirdonor/E_abs(dirdonor)*incrdonor;
  E_Int ind, i, j;
  E_Float invShiftDir1 = 1./shiftDir1;
  E_Float invShiftDir2 = 1./shiftDir2;
  //E_Int dim1D = (dim1-1)*shiftDir1+1;
  E_Float inv12 = invShiftDir1*invShiftDir2;
  E_Int indsav;
  E_Int sizeDir1Dnr = listIDSize/shiftDir2;
  if (isFine == 1) //finest zone
  {
    switch (dim)
    {
      case 3:
        ind = 0;
        for (j=0; j < dim2; j+=shiftDir2)
        {
          // ind=dim1D*j;
          for (i=0; i < dim1; i+=shiftDir1)
          {
            indadjdonor = arrayDonor[ind];
            for (E_Int jl = j; jl < j+shiftDir2;jl++)
              for (E_Int il = i; il < i+shiftDir1;il++)
              {
                indadj = arrayBorder[il+jl*dim1];              
                for (E_Int index = 0; index < d; index++)
                {
                  indghost = indadj+shift_dir*(index+1);
                  inddonor = indadjdonor-shift_dirdonor*index;              
                  dataR[indghost] = dataD[inddonor];
                }
              }
            ind++;// increment donneur
          }
          
        }
        break;
        
      case 2:
        ind = 0;
        for (i=0; i < dim1; i+=shiftDir1)
        {
          indadjdonor = arrayDonor[ind];
          for (E_Int il = i; il < i+shiftDir1;il++)
          {
            indadj = arrayBorder[il];
            for (E_Int index = 0; index < d; index++)
            {
              indghost = indadj+shift_dir*(index+1);
              inddonor = indadjdonor-shift_dirdonor*index;
              dataR[indghost] = dataD[inddonor];
            }
          }
          ind++;// increment donneur
        }        
        break;
    }
  }

  else if ( isFine == -1 ) // coarsest zone
  {
    // a mean of the corresponding centers on the finest zone is computed
    switch (dim)
    {
      case 3:
        ind = 0;
        for (j=0; j < dim2; j++)
        {
          for (i=0; i < dim1; i++)
          {
            indadj = arrayBorder[i+j*dim1];
            indsav = ind;
            for (E_Int index = 0; index < d; index++)
            {
              ind = indsav;
              indghost = indadj+shift_dir*(index+1);
              dataR[indghost] = 0.;
              for (E_Int il = 0; il < shiftDir1; il++)
              {
                for (E_Int jl = 0; jl < shiftDir2; jl++)           
                {
                  indadjdonor = arrayDonor[ind+jl*sizeDir1Dnr];            
                  inddonor = indadjdonor-shift_dirdonor*index;
                  dataR[indghost] += dataD[inddonor]; 
                }
                ind++;
              }                      
              dataR[indghost] *= inv12;
            }
          }
        }
        break;

      case 2:
        ind = 0;
        for (i=0; i < dim1; i++)
        {
          indadj = arrayBorder[i];
          for (E_Int index = 0; index < d; index++)
          {
            indghost = indadj+shift_dir*(index+1);
            dataR[indghost] = 0.;
            for (E_Int il = 0; il < shiftDir1; il++)
            {
              indadjdonor = arrayDonor[ind+il];            
              inddonor = indadjdonor-shift_dirdonor*index;
              dataR[indghost] += dataD[inddonor];   
            }
            dataR[indghost] *= invShiftDir1;
          }
          ind += shiftDir1;
        } 
        break;
    }
  }
  delete [] arrayDonor;
  Py_INCREF(Py_None);
  return Py_None;
}
