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

# include "connector.h"

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* update cellN for IBM (zone, in place) */
// ============================================================================
PyObject* K_CONNECTOR::_updateNatureForIBM(PyObject* self, PyObject* args)
{
  PyObject* zone;
  char* GridCoordinates; char* FlowSolutionNodes;
  char* FlowSolutionCenters;
  E_Int ibctype;
  if (!PYPARSETUPLE_(args, O_ I_ SSS_,
                    &zone, &ibctype, &GridCoordinates, &FlowSolutionNodes, &FlowSolutionCenters))
  {
    return NULL;
  }
  vector<PyArrayObject*> hook;
  E_Int im, jm, km, cnSize, cnNfld;
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  E_Int res = K_PYTREE::getFromZone(
    zone, 0, 1, varString, fields, locs, im, jm, km, 
    cn, cnSize, cnNfld, eltType, hook, GridCoordinates, 
    FlowSolutionNodes, FlowSolutionCenters);

  if (res != 1)
  {
    if (res == 2)
    { 
      delete [] eltType; delete [] varString;
      RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
      PyErr_SetString(PyExc_TypeError,
                      "updateNatureForIBM: not valid for unstructured grids.");
      return NULL;
    }
    else 
    {
      RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
      PyErr_SetString(PyExc_TypeError,
                      "updateNatureForIBM: invalid zone.");
      return NULL;
    }
  }  
  E_Int poscellni = K_ARRAY::isNamePresent("cellNIBC", varString);      
  E_Int poscellnc = K_ARRAY::isNamePresent("cellNChim", varString);      
  E_Int poscellnf = K_ARRAY::isNamePresent("cellNFront", varString);      
  if (poscellni == -1 || poscellnf == -1 || poscellnc == -1)
  {
    RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
    PyErr_SetString(PyExc_TypeError,
                    "updateNatureForIBM: zone must contain cellNIBC, cellNChim and cellNFront variables.");
    return NULL;
  }
  
  E_Float* ptrCellNIBC = fields[poscellni];
  E_Float* ptrCellNChim = fields[poscellnc];
  //E_Float* ptrCellNFront = fields[poscellnf];
  E_Int imc = K_FUNC::E_max(1, im-1);
  E_Int jmc = K_FUNC::E_max(1, jm-1);
  E_Int kmc = K_FUNC::E_max(1, km-1);
  E_Int ncells = imc*jmc*kmc;
#pragma omp parallel default(shared)
  {
#pragma omp for 
    for (E_Int ind = 0; ind < ncells; ind++)
    {
      E_Float& cellNChim = ptrCellNChim[ind];
      E_Float& cellNIBC = ptrCellNIBC[ind];
      //E_Float& cellNFront = ptrCellNFront[ind];

      if (K_FUNC::fEqualZero(cellNChim - 1.))
      {
        if (K_FUNC::fEqualZero(cellNIBC)) cellNChim = -3.;//~ blanked
      }
      else if (K_FUNC::fEqualZero(cellNChim - 2.))
      {
        if (K_FUNC::fEqualZero(cellNIBC - 1.)) cellNIBC = 3.;//not a donor
        else if (K_FUNC::fEqualZero(cellNIBC - 2.)) cellNIBC = -3.;
        else if (K_FUNC::fEqualZero(cellNIBC)) cellNChim = -3.;// ~ blanked
      }
      else if (K_FUNC::fEqualZero(cellNChim))
      {
        if (K_FUNC::fEqualZero(cellNIBC - 1.)) cellNIBC = -3.; //~ blanked
      }

      // c'est commente : on suppose que les corps IBM n'intersectent jamais les corps Chimere
      // if ( cellNFront != 0.)
      // {
      //   if (cellNIBC == -3.) cellNFront = 0.;
      // }
    }
  }
  delete [] eltType; delete [] varString;
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
  Py_INCREF(Py_None); 
  return Py_None;
}
