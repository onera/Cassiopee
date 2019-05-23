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

# include "converter.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* SetPartialFields + filter 
  Modify fields for points of indices in listIndices
  IN: z: zone to be updated
  IN: listIndices : numpy of indices of points to be updated
  IN: fArrays: fields in array2. 
  On suppose que tous les arrays de fArrays ont les champs dans la meme position 
  dans le tableau 
*/
//=============================================================================
PyObject* K_CONVERTER::filterPartialFields(PyObject* self, PyObject* args)
{
  E_Float ZEROVOL=1.e-16;
  E_Float penaltyExtrap = 1e6;
  E_Float penaltyOrphan = 1e12;
  E_Int countExtrap= 0; E_Int countOrphan = 0;

  char *GridCoordinates, *FlowSolutionNodes, *FlowSolutionCenters;
  PyObject* zone;
  PyObject* fArrays; // list of arrays to be filtered 
  PyObject* listIndicesO;
  E_Int loc;
  E_Int startFrom;
  char *filterName;// nom du filtre dans fArrays
  E_Int verbose;
  if (!PYPARSETUPLEI(args, "OOOllssssl", "OOOiissssi", &zone, &fArrays, &listIndicesO, &loc, &startFrom, 
                     &filterName, &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters, &verbose))
    return NULL; 

  /* zone a modifier */
  vector<PyArrayObject*> hook;
  E_Int imZ, jmZ, kmZ, cnZSize, cnZNfld;
  char* varStringZ; char* eltTypeZ;
  vector<E_Float*> fieldsZ; vector<E_Int> locsZ;
  vector<E_Int*> cnZ;

  K_PYTREE::getFromZone(zone, 0, loc, varStringZ,
                        fieldsZ, locsZ, imZ, jmZ, kmZ,
                        cnZ, cnZSize, cnZNfld, eltTypeZ, hook,
                        GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);

  E_Int nfldZ = fieldsZ.size();
  if (nfldZ == 0) 
  {
    RELEASESHAREDZ(hook, varStringZ, eltTypeZ);
    PyErr_SetString(PyExc_TypeError,
                    "filterPartialFields: no field to set.");
    return NULL;
  }

  // Check fArrays: champs a inserer
  vector<E_Int> resD;  vector<char*> varStringD;
  vector<FldArrayF*> fieldsD;
  vector<void*> a2; //ni,nj,nk ou cnt en NS
  vector<void*> a3; //eltType en NS
  vector<void*> a4;
  vector<PyObject*> objsD;
  E_Boolean skipNoCoord = false;  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = true;  E_Boolean skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(fArrays, resD, varStringD, fieldsD, a2, a3, a4, objsD,  
                                      skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nzonesD = objsD.size();
  if (isOk == -1)
  {
    RELEASESHAREDZ(hook, varStringZ, eltTypeZ);
    for (E_Int no = 0; no < nzonesD; no++)
      RELEASESHAREDA(resD[no],objsD[no],fieldsD[no],a2[no],a3[no],a4[no]); 
    PyErr_SetString(PyExc_TypeError,
                    "filterPartialFields: 2nd arg is not a valid list of arrays.");
    return NULL; 
  }

  /*--------------------------------------------*/
  /* Extraction des indices des pts a modifier */
  /*-------------------------------------------*/
  FldArrayI* listIndices;
  E_Int resi = K_NUMPY::getFromNumpyArray(listIndicesO, listIndices, true);
  if (resi == 0)
  {
    RELEASESHAREDZ(hook, varStringZ, eltTypeZ);
    for (E_Int no = 0; no < nzonesD; no++)
      RELEASESHAREDA(resD[no],objsD[no],fieldsD[no],a2[no],a3[no],a4[no]); 
    PyErr_SetString(PyExc_TypeError, 
                    "filterPartialFields: 3rd arg must be a numpy of integers.");
    return NULL;
  }

  E_Int nPts = listIndices->getSize();
  E_Int* indices = listIndices->begin();


  /* Variables communes entre la zone et les champs a inserer */
  char* varStringC; // chaine de caractere commune
  E_Int lenD = strlen(varStringD[0]);
  E_Int lenZ = strlen(varStringZ);
  E_Int l = min(lenD,lenZ);
  varStringC = new char [l+1];
  vector<E_Int> posvZ;//pos demarrent a 1
  vector<E_Int> posvD; 
  K_ARRAY::getPosition(varStringD[0], varStringZ, posvD, posvZ, varStringC);
  delete [] varStringC;
  
  E_Int ncommonfields = posvZ.size(); // nb de variables communes

  vector<E_Int> posfD;
  for (E_Int i = 0; i < nzonesD; i++)
  {
    E_Int posf = K_ARRAY::isNamePresent(filterName, varStringD[i]);
    if ( posf == -1) 
    {
        RELEASESHAREDZ(hook, varStringZ, eltTypeZ);
        RELEASESHAREDN(listIndicesO, listIndices);
        for (E_Int no = 0; no < nzonesD; no++)
          RELEASESHAREDA(resD[no],objsD[no],fieldsD[no],a2[no],a3[no],a4[no]); 
        PyErr_SetString(PyExc_TypeError, 
                        "filterPartialFields: filter variable not found in 2nd arg.");
        return NULL;        
    }
    posfD.push_back(posf+1);
  }
#pragma omp for
  for (E_Int i = 0; i < nPts; i++)
  {
    E_Float filterMax = K_CONST::E_MAX_FLOAT;
    E_Int bestDnr=-1;
    for (E_Int nozD=0; nozD < nzonesD; nozD++)
    {
      E_Int posf = posfD[nozD];//demarre a 1
      E_Float* ptrFilter = fieldsD[nozD]->begin(posf);
      E_Float filterVal = ptrFilter[i];
      if (filterVal < filterMax && filterVal> ZEROVOL)
      {
        filterMax = filterVal;
        bestDnr=nozD;
      }
    }
    if ( filterMax >= penaltyOrphan ) countOrphan++;
    else
    {
      if ( filterMax >= penaltyExtrap) {
        countExtrap++;}
      if ( bestDnr>-1)
      {
        E_Int ind = indices[i]-startFrom;

        for (E_Int eq = 0; eq < ncommonfields; eq++)
        {
          E_Int posZ = posvZ[eq]; 
          E_Float* fZ = fieldsZ[posZ-1];
          E_Int posD = posvD[eq]; 
          E_Float* ptrFieldD = fieldsD[bestDnr]->begin(posD);          
          fZ[ind] = ptrFieldD[i];
        }
      }
    }
  }

  if (countExtrap>0 || countOrphan>0)
  {
    E_Int countInterp = nPts-countOrphan-countExtrap;
    PyObject* v = PyList_GetItem(zone, 0);
    char* zname = NULL;
    if (PyString_Check(v)) zname = PyString_AsString(v);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(v)) zname = PyBytes_AsString(PyUnicode_AsUTF8String(v)); 
#endif
    
    if (verbose == 1) 
    {
      printf("Zone %s : interpolated=%d; extrapolated=%d; orphans=%d.\n",zname, countInterp, countExtrap, countOrphan);
      if (countOrphan > 0)
        printf("WARNING: Zone %s has %d orphan points.\n",zname,countOrphan);
    }
  }
  for (E_Int no = 0; no < nzonesD; no++)
    RELEASESHAREDA(resD[no],objsD[no],fieldsD[no],a2[no],a3[no],a4[no]); 
  RELEASESHAREDN(listIndicesO, listIndices);
  RELEASESHAREDZ(hook, varStringZ, eltTypeZ);
  Py_INCREF(Py_None);
  return Py_None;
}
