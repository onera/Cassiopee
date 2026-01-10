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
#include "converter.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/**/
//=============================================================================
PyObject* K_CONVERTER::extractBCFields(PyObject* self, PyObject* args)
{
  PyObject *zone, *pyIndices, *pyVariables;
  E_Int locI; // 0 = nodes, 1 centers
  char *GridCoordinates, *FlowSolutionNodes, *FlowSolutionCenters;
  if (!PYPARSETUPLE_(args, OOO_ I_ SSS_, &zone, &pyIndices, &pyVariables, 
                     &locI, 
                     &GridCoordinates, &FlowSolutionNodes, &FlowSolutionCenters)) 
    return NULL; 
  if (locI != 1)
  {
    PyErr_SetString(PyExc_TypeError, "extractBCFields: not yet implemented for fields located at nodes.");
    return NULL;
  }
  E_Int ni, nj, nk, cnSize, cnNfld;
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  vector<PyArrayObject*> hook;
  E_Int zoneType = K_PYTREE::getFromZone(zone, 0, locI, varString, fields, locs, ni, nj, nk, 
                                         cn, cnSize, cnNfld, eltType, hook, GridCoordinates, 
                                         FlowSolutionNodes, FlowSolutionCenters);
  if (zoneType == 0) 
  {
    PyErr_SetString(PyExc_TypeError, "extractBCFields: not a valid zone.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }
  E_Int* PE = NULL;
  if (zoneType == 2)
  {
    if (cn.size() < 3) //PE does not exist
    {
      PyErr_SetString(PyExc_TypeError, "extractBCFields: ParentElements node must be defined in zone.");
      RELEASESHAREDZ(hook, varString, eltType);
      return NULL;  
    }
    else PE = cn[2];
  }
  (void)PE;
  FldArrayI* indicesBC;
  E_Int res = K_NUMPY::getFromNumpyArray(pyIndices, indicesBC);

  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, "extractBCFields: not a valid numpy for indices of BC.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;   
  }
  /*--------------------------------------*/
  /*  Positions des variables a extraire  */
  /*--------------------------------------*/
  vector<E_Int> posvars;
  E_Int posvar;
  if (PyList_Check(pyVariables) != 0)
  {
    int nvariables = PyList_Size(pyVariables);
    if (nvariables > 0)
    {
      for (int i = 0; i < nvariables; i++)
      {
        PyObject* tpl0 = PyList_GetItem(pyVariables, i);
        if (PyString_Check(tpl0)) 
        {
          char* varname = PyString_AsString(tpl0);
          posvar = K_ARRAY::isNamePresent(varname, varString);      
          if (posvar != -1 ) posvars.push_back(posvar);  
        }
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(tpl0))
        {
           const char* varname = PyUnicode_AsUTF8(tpl0);
           posvar = K_ARRAY::isNamePresent(varname, varString);
           if (posvar != -1 ) posvars.push_back(posvar);  
        }
#endif
        else
        {
          PyErr_Warn(PyExc_Warning, "extractBCFields: variable must be a string. Skipped.");
        }
      }
    }
  }  
  E_Int nIndices = indicesBC->getSize();
  E_Int* ptrIndicesBC = indicesBC->begin();
  PyObject* PyListOfFieldsBC = PyList_New(0);

  E_Int nvars = posvars.size();
  if (zoneType == 1)//structure
  {
    E_Int dimZ = 3;
    if (nk==1) dimZ = 2;
        
    E_Int nic = max(E_Int(1),ni-1); E_Int njc = max(E_Int(1),nj-1); E_Int nkc = max(E_Int(1),nk-1);
    E_Int nbIntI = ni*njc*nkc; E_Int nbIntJ = nic*nj*nkc;
    E_Int nbIntIJ = nbIntI+nbIntJ;
    for (E_Int novar = 0; novar < nvars; novar++)     
    {
      PyObject* pyFieldBC = K_NUMPY::buildNumpyArray(nIndices, 1, 0, 1);  
      E_Float* ptrF = K_NUMPY::getNumpyPtrF(pyFieldBC);

      E_Int posv = posvars[novar];
      E_Float* fieldV = fields[posv];
      for (E_Int noint = 0; noint < nIndices; noint++)
      {
        E_Int indint = ptrIndicesBC[noint];
        E_Int kcell = 0;
        E_Int jcell = 0;
        E_Int icell = 0;

        if ( indint <nbIntI)// I-interface
        {
          if ( dimZ == 3) kcell = indint/(ni*njc);
          jcell = (indint-kcell*ni*njc)/ni;
          icell = indint -jcell*ni-kcell*ni*njc;                    
        }
        else if ( indint < nbIntI+nbIntJ)
        {
          //indint = icell + jcell * nic + kcell*nic*nj + nbIntI
          if ( dimZ == 3) kcell = (indint-nbIntI)/(nic*nj);
          jcell = ( indint-nbIntI-kcell*nic*nj)/nic;
          icell = indint-nbIntI -jcell*nic-kcell*nic*nj;                    
        }
        else 
        {
          //indint = icell + jcell * nic + kcell*nic*njc + nbIntIJ
          if ( dimZ == 3) kcell = (indint-nbIntIJ)/(nic*njc);
          jcell = ( indint-nbIntIJ-kcell*nic*njc)/nic;
          icell = indint-nbIntIJ -jcell*nic-kcell*nic*njc;   
        }
        E_Int indcell = icell + jcell*nic + kcell*nic*njc;
        ptrF[noint] = fieldV[indcell]; // extrap from adjacent cell
      }
      PyList_Append(PyListOfFieldsBC, pyFieldBC); Py_DECREF(pyFieldBC);                                                              
    }
  }
  else //non structure
  {
    for (E_Int novar = 0; novar < nvars; novar++)     
    {
      PyObject* pyFieldBC = K_NUMPY::buildNumpyArray(nIndices, 1, 0, 1);  
      E_Float* ptrF = K_NUMPY::getNumpyPtrF(pyFieldBC);
      E_Int* PE = cn[2];
      E_Int posv = posvars[novar];
      E_Float* fieldV = fields[posv];
      for (E_Int noint = 0; noint < nIndices; noint++)
      {
        E_Int indint = ptrIndicesBC[noint]-1;        
        E_Int indcell = PE[indint]-1;
        ptrF[noint] = fieldV[indcell];
      }
      PyList_Append(PyListOfFieldsBC, pyFieldBC); Py_DECREF(pyFieldBC);                                                              
    }
  }
  RELEASESHAREDN(pyIndices, indicesBC);
  RELEASESHAREDZ(hook, varString, eltType);
  return PyListOfFieldsBC;  
}
