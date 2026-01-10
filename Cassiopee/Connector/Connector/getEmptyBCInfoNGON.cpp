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
using namespace std;
using namespace K_FLD;

//=============================================================================
/* getEmptyBCInfoNGON: */
//=============================================================================
PyObject* K_CONNECTOR::_getEmptyBCInfoNGON(PyObject* self, PyObject* args)
{
  // INPUT 
  PyObject* ExteriorFacesZone;
  PyObject* ExteriorFaceIndices;
  PyObject* ExteriorDefinedFaceIndices;
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  if (!PYPARSETUPLE_(args, OOO_ SSS_,
                        &ExteriorFacesZone, &ExteriorFaceIndices, 
                        &ExteriorDefinedFaceIndices, &GridCoordinates, 
                        &FlowSolutionNodes, &FlowSolutionCenters))
  {
    return NULL;
  }
  // Recup du numpy des faces externes
  FldArrayI* exteriorFaceIndices;
  K_NUMPY::getFromNumpyArray(ExteriorFaceIndices, exteriorFaceIndices);
  // numpy des faces externes definies
  FldArrayI* exteriorDefinedFaceIndices;
  K_NUMPY::getFromNumpyArray(ExteriorDefinedFaceIndices, exteriorDefinedFaceIndices);
      
  vector<PyArrayObject*> hook;
  E_Int im, jm, km, cnSize, cnNfld;
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  E_Int res = K_PYTREE::getFromZone(ExteriorFacesZone, 0, 1, varString,
                                    fields, locs, im, jm, km,
                                    cn, cnSize, cnNfld, 
                                    eltType, hook, GridCoordinates, 
                                    FlowSolutionNodes, FlowSolutionCenters);
  E_Int err = 0;
  E_Int postag1 = K_ARRAY::isNamePresent("tag1", varString);      
  E_Int nextfaces = exteriorFaceIndices->getSize();
  E_Int ndefextfaces = exteriorDefinedFaceIndices->getSize();
  E_Int* exteriorFacesIndp = exteriorFaceIndices->begin();
  E_Int* exteriorDefinedFacesIndp = exteriorDefinedFaceIndices->begin();
  E_Float* tag1p = fields[postag1];
  if (res != 2) {err = 1; goto end;}
  if (postag1 == -1) {err = 2; goto end;}

  for (E_Int nofe = 0; nofe < nextfaces; nofe++)
  {
    E_Int indfaceExt = exteriorFacesIndp[nofe];
    for (E_Int nofed = 0; nofed < ndefextfaces; nofed++)
    {
      if (indfaceExt == exteriorDefinedFacesIndp[nofed])
      {
        tag1p[nofe] = -1000.;//definie
        break;
      }
    }
  }

  end:;
  delete [] eltType;
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL); 
  RELEASESHAREDN(ExteriorFaceIndices, exteriorFaceIndices);
  RELEASESHAREDN(ExteriorDefinedFaceIndices, exteriorDefinedFaceIndices);
  if (err == 0 )
  {
    Py_INCREF(Py_None);
    return Py_None;
  }
  
  if (err == 1)
    PyErr_SetString(PyExc_TypeError, 
                    "getEmptyBCInfoNGON: not valid for structured zones.");
  else if (err==2)
    PyErr_SetString(PyExc_TypeError, 
                    "getEmptyBCInfoNGON: tag1 or tag2 variable not found in zone.");

  
  return NULL;
}
