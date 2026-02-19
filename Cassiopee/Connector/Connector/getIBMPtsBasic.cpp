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
# include "Nuga/include/BbTree.h"
# include "Nuga/include/KdTree.h"
# include "Nuga/include/ArrayAccessor.h"
using namespace K_FLD;
using namespace std;
using namespace K_SEARCH;

#define RELEASEZONES \
  for (E_Int nos = 0; nos < nzonesS; nos++)\
    RELEASESHAREDS(objst[nos], structF[nos]);\
  for (E_Int nos = 0; nos < nzones; nos++)\
    RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);

// ============================================================================
/* Get the wall and image points if front is not provided 
    wall pts are obtained by projDir (projOrtho, closest pt if impossible)
    image pts are obtained by translation of corrected pts of hi/he*/
// ============================================================================
PyObject* K_CONNECTOR::getIBMPtsBasic(PyObject* self, PyObject* args)
{
  PyObject *allCorrectedPts, *distName, *normalNames;
  if (!PYPARSETUPLE_(args, OOO_, &allCorrectedPts, &normalNames, &distName))
    return NULL;
  // check distname
  char* distname;
  if (PyString_Check(distName)) distname = PyString_AsString(distName);
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(distName)) distname = (char*)PyUnicode_AsUTF8(distName);
#endif
  else
  {    
    PyErr_SetString(PyExc_TypeError, 
                    "getIBMPts: distName must be a string.");
    return NULL;
  }

  // Check normal components
  if (PyList_Check(normalNames) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "getIBMPts: normal vars must be a list of strings.");
    return NULL;
  }
  E_Int nvars = PyList_Size(normalNames);
  if (nvars != 3)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "getIBMPts: normal vars must be a 3-component vector.");
    return NULL;
  }
  char* var;
  vector<char*> varsn;// normal components names
  for (E_Int v = 0; v < 3; v++)
  {
    PyObject* l = PyList_GetItem(normalNames, v);
    if (PyString_Check(l))
    {
      var = PyString_AsString(l);
      varsn.push_back(var);
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(l)) 
    {
      var = (char*)PyUnicode_AsUTF8(l);
      varsn.push_back(var);
    } 
#endif
    else
    {
      PyErr_SetString(PyExc_TypeError,
                      "getIBMPts: invalid string for normal component.");
      return NULL;
    }
  }
  // Extract correctedPts
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt; vector<char*> eltType;
  vector<PyObject*> objst, objut;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = true;
  E_Bool skipUnstructured = false;
  E_Bool skipDiffVars = false;
  E_Int isOk = K_ARRAY::getFromArrays(allCorrectedPts, resl, structVarString, unstrVarString,
                                      structF, unstrF, nit, njt, nkt, cnt, eltType, objst, objut, 
                                      skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nzones = objut.size(); E_Int nzonesS = objst.size();
  if (isOk == -1)   
  {
    PyErr_SetString(PyExc_TypeError,"getIBMPts: 1st arg is not valid.");
    RELEASEZONES;
    return NULL;
  }
  for (E_Int no = 0; no < nzones; no++)
  {
    if (K_STRING::cmp(eltType[no], "NODE") != 0) 
    {
      PyErr_SetString(PyExc_TypeError,"getIBMPts: 1st arg must be a NODE type zone.");
      RELEASEZONES;
      return NULL;
    }
  }
  // corrected points coordinates and normal components: position 
  E_Int posx1, posy1, posz1, poshi, poshe;
  vector<E_Int> posxt; vector<E_Int> posyt; vector<E_Int> poszt;
  vector<E_Int> posnxt; vector<E_Int> posnyt; vector<E_Int> posnzt; vector<E_Int> posdist;
  vector<E_Int> poshit; vector<E_Int> poshet;
  for (E_Int no = 0; no < nzones; no++)        
  {
    posx1 = K_ARRAY::isCoordinateXPresent(unstrVarString[no]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(unstrVarString[no]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(unstrVarString[no]); posz1++;
    posxt.push_back(posx1); posyt.push_back(posy1); poszt.push_back(posz1);
    poshi = K_ARRAY::isNamePresent("hi",unstrVarString[no]); poshi++;
    poshe = K_ARRAY::isNamePresent("he",unstrVarString[no]); poshe++;
    poshit.push_back(poshi); poshet.push_back(poshe);
    posx1 = K_ARRAY::isNamePresent(varsn[0], unstrVarString[no]);      
    posy1 = K_ARRAY::isNamePresent(varsn[1], unstrVarString[no]);      
    posz1 = K_ARRAY::isNamePresent(varsn[2], unstrVarString[no]);           
    if (posx1==-1 || posy1==-1 || posz1==-1) 
    {
      PyErr_SetString(PyExc_TypeError,"getIBMPts: one of normal components not found in 1st arg.");
      RELEASEZONES; 
      return NULL;
    }
                
    posx1++; posy1++; posz1++;
    posnxt.push_back(posx1); posnyt.push_back(posy1); posnzt.push_back(posz1);
        
    posx1 = K_ARRAY::isNamePresent(distname, unstrVarString[no]);   
    if ( posx1 == -1)
    {
      PyErr_SetString(PyExc_TypeError,"getIBMPts: dist variable is not present in 1st arg.");
      RELEASEZONES; return NULL;
    }
    posx1++; posdist.push_back(posx1);
  }

  // Checks completed...build arrays
  PyObject* PyListOfImagePts = PyList_New(0);  
  PyObject* PyListOfWallPts  = PyList_New(0);

  E_Int nvarOut = 3;
  char varStringOut[K_ARRAY::VARSTRINGLENGTH]; 
  strcpy(varStringOut,"CoordinateX,CoordinateY,CoordinateZ");
  char eltTypeOut[8]; strcpy(eltTypeOut,"NODE");
  FldArrayI cnl(0);
  E_Int nelts = 0;
  vector<E_Float*> xit; vector<E_Float*> yit; vector<E_Float*> zit;
  vector<E_Float*> xwt; vector<E_Float*> ywt; vector<E_Float*> zwt;

  for (E_Int no = 0; no < nzones; no++)
  {
    E_Int npts = unstrF[no]->getSize();
    PyObject* tpli = K_ARRAY::buildArray(nvarOut, varStringOut, npts, nelts, -1, eltTypeOut, false);
    E_Float* coordip = K_ARRAY::getFieldPtr(tpli);
    FldArrayF coordi(npts, 3, coordip, true); // non initialized yet
    coordi.setOneField(*unstrF[no],posxt[no],1);
    coordi.setOneField(*unstrF[no],posyt[no],2);
    coordi.setOneField(*unstrF[no],poszt[no],3);
    xit.push_back(coordi.begin(1));
    yit.push_back(coordi.begin(2));
    zit.push_back(coordi.begin(3));
    PyList_Append(PyListOfImagePts, tpli); Py_DECREF(tpli);

    PyObject* tplw = K_ARRAY::buildArray(nvarOut, varStringOut, npts, nelts, -1, eltTypeOut, false);
    E_Float* coordwp = K_ARRAY::getFieldPtr(tplw);
    FldArrayF coordw(npts, 3, coordwp, true); // non initialized yet
    coordw.setOneField(*unstrF[no],posxt[no],1);
    coordw.setOneField(*unstrF[no],posyt[no],2);
    coordw.setOneField(*unstrF[no],poszt[no],3);
    xwt.push_back(coordw.begin(1));
    ywt.push_back(coordw.begin(2));
    zwt.push_back(coordw.begin(3));
    PyList_Append(PyListOfWallPts, tplw); Py_DECREF(tplw);        
  }

  /* Algorithm */
  E_Float dist0, delta;
  E_Float dirx, diry, dirz ,dirn;
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    FldArrayF* correctedPts = unstrF[noz];
    E_Int npts = correctedPts->getSize();
    E_Float* ptrNX = correctedPts->begin(posnxt[noz]);
    E_Float* ptrNY = correctedPts->begin(posnyt[noz]);
    E_Float* ptrNZ = correctedPts->begin(posnzt[noz]);
    E_Float* ptrXC = correctedPts->begin(posxt[noz]);
    E_Float* ptrYC = correctedPts->begin(posyt[noz]);
    E_Float* ptrZC = correctedPts->begin(poszt[noz]);
    E_Float* ptrXW = xwt[noz];
    E_Float* ptrYW = ywt[noz];
    E_Float* ptrZW = zwt[noz];
    E_Float* ptrXI = xit[noz];
    E_Float* ptrYI = yit[noz];
    E_Float* ptrZI = zit[noz];
    E_Float* distp  = correctedPts->begin(posdist[noz]);
    E_Float* hit = correctedPts->begin(poshit[noz]);
    E_Float* het = correctedPts->begin(poshet[noz]);

    for (E_Int ind = 0; ind < npts; ind++)
    {
      //xc0 = ptrXC[ind]; yc0 = ptrYC[ind]; zc0 = ptrZC[ind];
      dist0 = distp[ind];
      dirx = ptrNX[ind]; diry = ptrNY[ind]; dirz = ptrNZ[ind];
      dirn = sqrt(dirx*dirx+diry*diry+dirz*dirz);
      dirx = dirx/dirn; // normalized
      diry = diry/dirn; 
      dirz = dirz/dirn;

      // wall pts
      ptrXW[ind] = ptrXC[ind] - dist0 * dirx;
      ptrYW[ind] = ptrYC[ind] - dist0 * diry;
      ptrZW[ind] = ptrZC[ind] - dist0 * dirz;

      // image pts
      if (dist0 > K_CONST::E_GEOM_CUTOFF) delta = het[ind];
      else
      {
        if (hit[ind] <= K_CONST::E_GEOM_CUTOFF) delta = -dist0;//symmetrical if hi=0
        else delta = hit[ind];
      }

      ptrXI[ind] = ptrXW[ind] + delta * dirx;
      ptrYI[ind] = ptrYW[ind] + delta * diry;
      ptrZI[ind] = ptrZW[ind] + delta * dirz;
    }
  }

  // Sortie
  PyObject* tpl = Py_BuildValue("[OO]", PyListOfWallPts, PyListOfImagePts);
  Py_DECREF(PyListOfWallPts);
  Py_DECREF(PyListOfImagePts);
  RELEASEZONES;
  return tpl;
}
