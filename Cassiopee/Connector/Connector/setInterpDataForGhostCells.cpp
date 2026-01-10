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
# include "connector.h"
using namespace K_FLD;


//=============================================================================
/* Returns the indices of ghost points (cells or nodes) and the indices of 
   donor points according to the 1-to-1 grid connectivity for structured grids */
//=============================================================================
PyObject* K_CONNECTOR::setInterpDataForGC(PyObject* self, PyObject* args)
{
  E_Int dim; // pb dimension (1, 2 or 3 for 1D, 2D or 3D)
  E_Int loc; // 0 : nodes, 1 : centers
  E_Int incrR, incrD;
  E_Int depth;
  PyObject* pyArrayBorderIndices;
  PyObject* listOfDonorIndices;

  if (!PYPARSETUPLE_(args, OO_ IIII_ I_,
                    &pyArrayBorderIndices, &listOfDonorIndices,
                    &dim, &loc, &depth, &incrR, &incrD))
  {
    return NULL;
  }
  
  FldArrayI* arrayBorderI;
  E_Int res = K_NUMPY::getFromNumpyArray(pyArrayBorderIndices, arrayBorderI);
  if (res == 0) 
  {    
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpDataForGhostCells: 1st arg is not a valid numpy array.");
    return NULL;
  }
  E_Int* arrayBorder = arrayBorderI->begin();

  // check if listOfDonorIndices is a list
  if (PyList_Check (listOfDonorIndices) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpDataForGhostCells: 2nd arg must be a list.");
    RELEASESHAREDN(pyArrayBorderIndices, arrayBorderI);
    return NULL;
  }
  E_Int listDnrIndSize = PyList_Size(listOfDonorIndices);

  // Convert Python listOfDonorIndices: 
  // array of indices of donor block cells
  E_Int* arrayDonor = new E_Int[listDnrIndSize]; 
  for (E_Int i = 0; i < listDnrIndSize; i++)
  {
    PyObject* tpl = PyList_GetItem(listOfDonorIndices, i);
    arrayDonor[i] = PyInt_AsLong(tpl);
  }
  //
  E_Int arrayBorderSize1 = arrayBorderI->getSize();
  E_Int arrayBorderSize2 = arrayBorderI->getNfld();

  E_Int ngc = arrayBorderSize1*arrayBorderSize2*depth;
  
  FldArrayI* donorInd1D = new FldArrayI(ngc); E_Int* donorPtr = donorInd1D->begin(); 
  FldArrayI* rcvInd1D = new FldArrayI(ngc); E_Int* rcvPtr = rcvInd1D->begin();
  FldArrayF* coefs = new FldArrayF(ngc); coefs->setAllValuesAt(1.);
  FldArrayI* donorTypes = new FldArrayI(ngc); donorTypes->setAllValuesAt(1); // coincident

  // local values
  E_Int indadjr;         // index of receptor block cell
  E_Int indadjd;         // index of donor block cell
  E_Int inddonor, indghost;
  E_Int ind;
  //E_Int loc0 = 0; // voir si en Vertex c est toujours valable

  E_Int indg = 0;
  E_Int noind = 0;
  switch (dim)
  {
    case 3:
      for (E_Int j = 0; j < arrayBorderSize2; j++)
        for (E_Int i = 0; i < arrayBorderSize1; i++)
        {
          ind = i + j*arrayBorderSize1;
          indadjr = arrayBorder[ind];// index of the 1st interieur pt near the border
          indadjd = arrayDonor[noind];// corresponding opposite index  

          for (E_Int index = 0; index < depth; index++)
          {
            indghost = indadjr+incrR*index; inddonor=indadjd+incrD*index;
            donorPtr[indg] = inddonor; rcvPtr[indg] = indghost; indg++;
          }
          noind++;
        }
      break;

    case 2:
      for (E_Int i = 0; i < arrayBorderSize1; i++)
      {
        indadjr = arrayBorder[i];// index of the 1st interieur pt near the border
        indadjd = arrayDonor[i];// corresponding opposite index  

        for (E_Int index = 0; index < depth; index++)
        {
          indghost = indadjr+incrR*index; inddonor=indadjd+incrD*index;
          donorPtr[indg] = inddonor; rcvPtr[indg] = indghost; indg++;
        }
      }
      break;
      
      break;
    default:
      PyErr_SetString(PyExc_TypeError, 
                      "setInterpDataForGhostCells: dimension is not valid.");
      RELEASESHAREDN(pyArrayBorderIndices, arrayBorderI);
      delete [] arrayDonor;
      return NULL;
      
  }  
  RELEASESHAREDN(pyArrayBorderIndices, arrayBorderI);
  delete [] arrayDonor;

  //buildNumpy with Fortran storage (1)
  // donorIndices1D
  PyObject* PyDonorInd = K_NUMPY::buildNumpyArray(*donorInd1D,1);
  delete donorInd1D;

  // receptorIndices1D
  PyObject* PyRcvInd = K_NUMPY::buildNumpyArray(*rcvInd1D,1);
  delete rcvInd1D;
 
  // donorType
  PyObject* PyDonorTypes = K_NUMPY::buildNumpyArray(*donorTypes,1);
  delete donorTypes;

  // coefs
  PyObject* PyCoefs = K_NUMPY::buildNumpyArray(*coefs,1);
  delete coefs;

  PyObject* tpl = Py_BuildValue("[OOOO]", PyRcvInd, PyDonorInd, PyDonorTypes, PyCoefs);
  Py_DECREF(PyRcvInd); Py_DECREF(PyDonorInd);  Py_DECREF(PyDonorTypes); Py_DECREF(PyCoefs);
  return tpl;
}
