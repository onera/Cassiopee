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
// extrapInterior2BCFaceStruct

# include "converter.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
/* set values in BCFace by extrapolation from interior */
//=============================================================================
PyObject* K_CONVERTER::extrapInterior2BCFaceStruct(PyObject* self, PyObject* args)
{
  E_Int imin, imax, jmin, jmax, kmin, kmax;
  PyObject *zone, *dataBC;
  E_Int Loc;
  char* varname; char *GridCoordinates, *FlowSolutionNodes, *FlowSolutionCenters;
  if (!PYPARSETUPLE_(args, OO_ IIII_ III_ SSSS_, 
                     &zone, &dataBC, &imin, &imax, &jmin, &jmax, &kmin, &kmax, 
                     &Loc, &varname, &GridCoordinates, &FlowSolutionNodes, &FlowSolutionCenters)) return NULL;

  vector<PyArrayObject*> hook;
  E_Int im, jm, km, cnSize, cnNfld;
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  K_PYTREE::getFromZone(zone, 0, Loc, varString, fields, locs, im, jm, km, 
                        cn, cnSize, cnNfld, eltType, hook, GridCoordinates, 
                        FlowSolutionNodes, FlowSolutionCenters);
  E_Int posf = K_ARRAY::isNamePresent(varname,varString);
  if (posf == -1) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "extrapInterior2BCFaceStruct: variable not found in zone.");
    RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
    return NULL; 
  }
  E_Float* ptrField = fields[posf];

  // flow field at BC (Vertex or Faces)
  FldArrayF* fieldInt;
  K_NUMPY::getFromNumpyArray(dataBC, fieldInt);
  E_Float* ptrBCField = fieldInt->begin();

  E_Int imjm = im*jm;
  E_Int noindint = 0;
  E_Int indv, shift;
  if ( Loc == 0)// vertex
  {
    for (E_Int k = kmin-1; k < kmax; k++)
      for (E_Int j = jmin-1; j < jmax; j++)
        for (E_Int i = imin-1; i < imax; i++)
        {
          indv = i+j*im+k*imjm;
          ptrBCField[noindint] = ptrField[indv];
          noindint++;
        }
  }
  else //Face centers
  {
    E_Int im1 = K_FUNC::E_max(1, im-1); 
    E_Int jm1 = K_FUNC::E_max(1, jm-1); 
    E_Int km1 = K_FUNC::E_max(1, km-1); 
    E_Int im1jm1 = im1*jm1;
    if (imin == imax)
    {
      if (imin == 1) shift = 0;
      else shift= im1-1;  
      for (E_Int k = kmin-1; k < kmax-1; k++)
        for (E_Int j = jmin-1; j < jmax-1; j++)
        {
          indv=shift+j*im1+k*im1jm1;
          ptrBCField[noindint]=ptrField[indv];
          noindint++;
        }
    }
    
    else if ( jmin == jmax)
    {
      if (jmin == 1) shift = 0;
      else shift= (jm1-1)*im1;  
      for (E_Int k = kmin-1; k < kmax-1; k++)
        for (E_Int i = imin-1; i < imax-1; i++)
        {
          indv=i+shift+k*im1jm1;
          ptrBCField[noindint]=ptrField[indv];
          noindint++;
        }
    }      
    
    else if ( kmin == kmax)
    {
      if (kmin == 1) shift = 0;
      else shift= (km1-1)*im1jm1;  
      for (E_Int j = jmin-1; j < jmax-1; j++)
        for (E_Int i = imin-1; i < imax-1; i++)
        {
          indv=i+j*im1+shift;
          ptrBCField[noindint]=ptrField[indv];
          noindint++;
        }      
    }
  }
  RELEASESHAREDN(dataBC, fieldInt);
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
  Py_INCREF(Py_None);
  return Py_None;  
}
