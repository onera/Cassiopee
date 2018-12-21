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
// Tag degenerated cells with tag1=-3.

# include "connector.h"

using namespace K_FUNC;
using namespace std;
using namespace K_FLD;
//=============================================================================
// Identify matching points in windows
//=============================================================================
PyObject* K_CONNECTOR::identifyDegenerated(PyObject* self, PyObject* args)
{
  PyObject *listOfAllWins;
  E_Float tol;
  if (!PYPARSETUPLEF(args,
                    "Od", "Of",
                    &listOfAllWins, &tol))
  {
      return NULL;
  }

  if (PyList_Check(listOfAllWins) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyDegenerated: 1st argument must be a list.");
    return NULL;
  }
  E_Float minvol = tol*tol;
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = true;
  E_Boolean skipDiffVars = true;

  E_Int isOk = K_ARRAY::getFromArrays(
    listOfAllWins, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltTypet, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nzones = structF.size();
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError, "identifyDegenerated: 1st list of arrays is not valid.");
    for (E_Int is = 0; is < nzones; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    return NULL;
  } 
  if (nzones == 0) 
  {
    PyErr_SetString(PyExc_TypeError, "identifyDegenerated: 1st list of arrays does not contain valid structured zones.");
    for (E_Int is = 0; is < nzones; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    return NULL;
  }
  vector<E_Int> post1; vector<E_Int> posvt;
  E_Int posvol, postag1;
  for (E_Int i = 0; i < nzones; i++)
  {   
    postag1 = K_ARRAY::isNamePresent("tag1", structVarString[i]); postag1++;
    posvol  = K_ARRAY::isNamePresent("vol", structVarString[i]); posvol++;
    if (postag1 < 1 || posvol < 1) 
    {
      PyErr_SetString(PyExc_TypeError, "identifyDegenerated: tag1 and vol must be defined in listOfAllWins.");
      for (E_Int is = 0; is < nzones; is++)
        RELEASESHAREDS(objst[is], structF[is]);
      return NULL;
    } 
    post1.push_back(postag1); posvt.push_back(posvol); 
  }
  /*-------------------- Fin des verifs -------------------------------------*/
  PyObject* l = PyList_New(0);
  vector<E_Float*> fields;
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    PyObject* tpl = K_ARRAY::buildArray(
      structF[noz]->getNfld(), structVarString[noz],
      nit[noz], njt[noz], nkt[noz]);
    E_Float* fp = K_ARRAY::getFieldPtr(tpl);
    E_Int npts = structF[noz]->getSize();
    FldArrayF ftemp0(npts,structF[noz]->getNfld(), fp, true); ftemp0 = *structF[noz];
    posvol = posvt[noz]; postag1 = post1[noz];
    E_Float* volp = ftemp0.begin(posvol);
    E_Float* tagp = ftemp0.begin(postag1);
    for (E_Int ind = 0; ind < npts; ind++)
    {
      if (K_FUNC::E_abs(volp[ind]) < minvol) tagp[ind] = -3.;
    }
    fields.push_back(fp);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
 
  for (E_Int is = 0; is < nzones; is++) RELEASESHAREDS(objst[is], structF[is]);
  return l;
}
