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

# include "converter.h"
# include "Nuga/include/KdTree.h"

using namespace K_FLD;
using namespace std;
using namespace K_SEARCH;

// ============================================================================
/* Identify points from donor zones with points of receptor zones
   and set the solution from donor to receptor
   IN: coordsRcv: ones receveuses 
   IN: hook: stockage des coordonnees des donneurs
   IN: soln: solution des zones donneuses localisees comme coordsDnr
   IN: match avec une tolerance tol */
// ============================================================================
PyObject* K_CONVERTER::identifySolutions(PyObject* self, PyObject* args)
{
  PyObject *hook, *coordsRcv, *solDnr;
  E_Float tol;
  
  if (!PYPARSETUPLE_(args, OOO_ R_, &hook, &solDnr, &coordsRcv, &tol))
    return NULL;

  E_Float tol2 = tol*tol;

// recupere le hook
  void** packet = NULL;
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
  packet = (void**) PyCObject_AsVoidPtr(hook);
#else
  packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
  E_Int* type = (E_Int*)packet[0];
  if (*type != 100 && *type != 102 &&  *type != 103)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifySolutions: this function requires a hook on a global kdtree.");
    return NULL;
  }
  //FldArrayF* coordsDnr = (FldArrayF*)packet[1];
  K_SEARCH::KdTree<FldArrayF>* globalKdt = (K_SEARCH::KdTree<FldArrayF>*) packet[3];
  
  // Receptor zones
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = false;
  E_Bool skipUnstructured = false;
  E_Bool skipDiffVars = false;
  
  vector<E_Int> resR;  vector<char*> varStringR;
  vector<FldArrayF*> coordsR;
  vector<void*> aR2; //ni,nj,nk ou cnt en NS
  vector<void*> aR3; //eltType en NS
  vector<void*> aR4;
  vector<PyObject*> objsR;
  E_Int isOk = K_ARRAY::getFromArrays(
    coordsRcv, resR, varStringR, coordsR, aR2, aR3, aR4, objsR,  
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nRcvZones = objsR.size();
  if (isOk == -1)
  {
    for (E_Int no = 0; no < nRcvZones; no++)
      RELEASESHAREDA(resR[no],objsR[no],coordsR[no],aR2[no],aR3[no],aR4[no]);
    PyErr_SetString(PyExc_TypeError,
                    "identifySolutions: 3rd argument is not valid.");
    return NULL;
  }   
  if (nRcvZones == 0) 
  { 
    PyErr_SetString(PyExc_TypeError,
                    "identifySolutions: 3rd arg: no valid zone found.");
    return NULL;
  }

  // Donor zones : fields
  skipNoCoord = false;
  skipStructured = false;
  skipUnstructured = false;
  skipDiffVars = true;
  
  vector<E_Int> res;  vector<char*> varString;
  vector<FldArrayF*> fields;
  vector<void*> a2; //ni,nj,nk ou cnt en NS
  vector<void*> a3; //eltType en NS
  vector<void*> a4;
  vector<PyObject*> objs;
  isOk = K_ARRAY::getFromArrays(
    solDnr, res, varString, fields, a2, a3, a4, objs,  
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nDnrFields = objs.size();
  if (isOk == -1)
  {
    for (E_Int no = 0; no < nRcvZones; no++)
      RELEASESHAREDA(resR[no],objsR[no],coordsR[no],aR2[no],aR3[no],aR4[no]);
    for (E_Int no = 0; no < nDnrFields; no++)
      RELEASESHAREDA(res[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
    PyErr_SetString(PyExc_TypeError,
                    "identifySolutions: 2nd argument is not valid.");
    return NULL;
  }
  // Build arrays
  PyObject* l = PyList_New(0);
  PyObject* tpl;  
  E_Int nfldout = fields[0]->getNfld();
  char* varStringOut;
  if (varString.size() > 0)
  {
    varStringOut = new char [strlen(varString[0])+2];
    strcpy(varStringOut, varString[0]);
  }
   else
  {
    varStringOut = new char [2]; varStringOut[0] = '\0';
  }

  FldArrayI indirZones(nDnrFields);
  E_Int sizeP = 0;
  for (E_Int nod = 0; nod < nDnrFields; nod++)
  { 
    sizeP += fields[nod]->getSize();
    indirZones[nod] = sizeP;
  }
  E_Int posx1, posy1, posz1;

  for (E_Int nor = 0; nor < nRcvZones; nor++)
  {
    posx1 = K_ARRAY::isCoordinateXPresent(varStringR[nor]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(varStringR[nor]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(varStringR[nor]); posz1++;
    E_Int nptsR = coordsR[nor]->getSize();
    E_Float* xR = coordsR[nor]->begin(posx1);
    E_Float* yR = coordsR[nor]->begin(posy1);
    E_Float* zR = coordsR[nor]->begin(posz1);
    FldArrayF* fout = new FldArrayF(nptsR, nfldout); fout->setAllValuesAtNull();
    E_Int api = fields[0]->getApi();
    if (resR[nor] == 1)
    {
      E_Int nir = *(E_Int*)aR2[nor];
      E_Int njr = *(E_Int*)aR3[nor];
      E_Int nkr = *(E_Int*)aR4[nor];
      tpl = K_ARRAY::buildArray3(nfldout, varStringOut, nir, njr, nkr, api);
    }
    else // unstr: NGON or BE/ME
    {
      FldArrayI* cnR = (FldArrayI*)aR2[nor];
      char* eltTypeR = (char*)aR3[nor];
      tpl = K_ARRAY::buildArray3(*fout, varStringOut, *cnR, eltTypeR, api);
    }
    FldArrayI* cnout;
    K_ARRAY::getFromArray3(tpl, fout, cnout);
    
    // identification
#pragma omp parallel default(shared)
    {
      E_Float pt[3];
      E_Int indkdt, indD, noblkD;
      E_Float ds2;
#pragma omp for schedule(dynamic)
      for (E_Int indR = 0; indR < nptsR; indR++)
      {
        pt[0] = xR[indR]; pt[1] = yR[indR]; pt[2] = zR[indR];
        indkdt = globalKdt->getClosest(pt, ds2);
        if (ds2 < tol2 && indkdt > -1)
        {
          indD = -1; noblkD = -1;
          // recherche du numero du bloc correspondant
          // algo dichotomique a brancher eventuellement       
          for (E_Int nod = 0; nod < nDnrFields; nod++)
          {
            if (indkdt < indirZones[nod]) 
            {
              noblkD = nod;
              if (nod == 0) indD = indkdt;
              else indD = indkdt-indirZones[nod-1];
              break;
            }
          }
          FldArrayF& fieldD = *fields[noblkD];
          for (E_Int eq = 1; eq <= nfldout; eq++)
            (*fout)(indR,eq) = fieldD(indD,eq); 
        }
      }
    }
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  delete [] varStringOut;

  for (E_Int no = 0; no < nRcvZones; no++)
    RELEASESHAREDA(resR[no],objsR[no],coordsR[no],aR2[no],aR3[no],aR4[no]); 
  for (E_Int no = 0; no < nDnrFields; no++)
    RELEASESHAREDA(res[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
  return l;
}
