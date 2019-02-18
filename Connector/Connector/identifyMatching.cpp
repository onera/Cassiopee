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

// Identify matching points in windows

# include "connector.h"
# include "Connect/merge.h"
# include "Search/KdTree.h"
# include "Fld/ArrayAccessor.h"

using namespace K_FUNC;
using namespace std;
using namespace K_FLD;

//=============================================================================
/* Identify matching points in windows
   arrays listOfAllWins must contain first structured arrays, then NGON arrays
*/
//=============================================================================
PyObject* K_CONNECTOR::identifyMatchingP(PyObject* self, PyObject* args)
{
  PyObject *listOfAllWins, *listOfAllWinsP;
  E_Float tolmatch;
  if (!PYPARSETUPLEF(args,
                    "OOd", "OOf",
                    &listOfAllWins, &listOfAllWinsP, &tolmatch))
  {
      return NULL;
  }

  if (PyList_Check(listOfAllWins) == 0 )
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyMatching: 1st argument must be a list.");
    return NULL;
  }
  if (PyList_Check(listOfAllWinsP) == 0 )
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyMatching: 2nd argument must be a list.");
    return NULL;
  }
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false;
  E_Boolean skipDiffVars = false;//true;
  E_Int isOk = K_ARRAY::getFromArrays(
    listOfAllWins, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltTypet, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int ns = structF.size(); E_Int nu = unstrF.size();

  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError, "identifyMatching: 1st list of arrays is not valid.");
    for (E_Int is = 0; is < ns; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    return NULL;
  } 
  if (ns+nu == 0) 
  {
    PyErr_SetString(PyExc_TypeError, "identifyMatching: 1st list of arrays does not contain valid zones.");
    for (E_Int is = 0; is < ns; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    return NULL;
  }
  vector<E_Int> posxt; vector<E_Int> posyt; vector<E_Int> poszt; 
  vector<E_Int> posvt; vector<E_Int> post1; vector<E_Int> post2;
  E_Int posxi, posyi, poszi, posvi, postag1, postag2;
  for (E_Int i = 0; i < ns; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(structVarString[i]);
    posyi = K_ARRAY::isCoordinateYPresent(structVarString[i]);
    poszi = K_ARRAY::isCoordinateZPresent(structVarString[i]);
    posvi = K_ARRAY::isNamePresent("vol", structVarString[i]);
    postag1 = K_ARRAY::isNamePresent("tag1", structVarString[i]);
    postag2 = K_ARRAY::isNamePresent("tag2", structVarString[i]);
    if (postag1 < 0 || postag2 < 0 || posvi <0) 
    {
      PyErr_SetString(PyExc_TypeError, "identifyMatching: vol, tag1 and tag2 must be defined in listOfAllWins.");
      for (E_Int is = 0; is < ns; is++)
        RELEASESHAREDS(objst[is], structF[is]);
      for (E_Int nos = 0; nos < nu; nos++)
        RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
      return NULL;
    } 
    posxt.push_back(posxi); posyt.push_back(posyi); poszt.push_back(poszi);
    posvt.push_back(posvi); post1.push_back(postag1); post2.push_back(postag2); 
  }
  for (E_Int i = 0; i < nu; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(unstrVarString[i]);
    posyi = K_ARRAY::isCoordinateYPresent(unstrVarString[i]);
    poszi = K_ARRAY::isCoordinateZPresent(unstrVarString[i]);
    posvi = K_ARRAY::isNamePresent("vol", unstrVarString[i]);
    postag1 = K_ARRAY::isNamePresent("tag1", unstrVarString[i]);
    postag2 = K_ARRAY::isNamePresent("tag2", unstrVarString[i]);
    if (postag1 < 0 || postag2 < 0 || posvi <0) 
    {
      PyErr_SetString(PyExc_TypeError, "identifyMatching: vol, tag1 and tag2 must be defined in listOfAllWins.");
      for (E_Int is = 0; is < ns; is++)
        RELEASESHAREDS(objst[is], structF[is]);
      for (E_Int nos = 0; nos < nu; nos++)
        RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
      return NULL;
    } 
    if (strcmp(eltTypet[i],"NGON") != 0 )
    {
      PyErr_SetString(PyExc_TypeError, "identifyMatching: unstructured arrays must be NGON.");
      for (E_Int is = 0; is < ns; is++)
        RELEASESHAREDS(objst[is], structF[is]);
      for (E_Int nos = 0; nos < nu; nos++)
        RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
      return NULL;
    } 
    posxt.push_back(posxi); posyt.push_back(posyi); poszt.push_back(poszi);
    posvt.push_back(posvi); post1.push_back(postag1); post2.push_back(postag2); 
  }
  // periodiques
  vector<E_Int> reslp;
  vector<char*> structVarStringp; vector<char*> unstrVarStringp;
  vector<FldArrayF*> structFp; vector<FldArrayF*> unstrFp;
  vector<E_Int> nitp; vector<E_Int> njtp; vector<E_Int> nktp;
  vector<FldArrayI*> cntp;
  vector<char*> eltTypetp;
  vector<PyObject*> objstp, objutp;
  isOk = K_ARRAY::getFromArrays(listOfAllWinsP, reslp, structVarStringp, unstrVarStringp,
                                structFp, unstrFp, nitp, njtp, nktp, cntp, eltTypetp, objstp, objutp, 
                                skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  vector<E_Int> posxp; vector<E_Int> posyp; vector<E_Int> poszp;
  for (E_Int i = 0; i < ns; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(structVarStringp[i]);
    posyi = K_ARRAY::isCoordinateYPresent(structVarStringp[i]);
    poszi = K_ARRAY::isCoordinateZPresent(structVarStringp[i]);
    posxp.push_back(posxi); posyp.push_back(posyi); poszp.push_back(poszi);
  }
  for (E_Int i = 0; i < nu; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(unstrVarStringp[i]);
    posyi = K_ARRAY::isCoordinateYPresent(unstrVarStringp[i]);
    poszi = K_ARRAY::isCoordinateZPresent(unstrVarStringp[i]);
    if (strcmp(eltTypetp[i],"NGON") != 0 )
    {
      PyErr_SetString(PyExc_TypeError, "identifyMatching: unstructured arrays must be NGON.");
      for (E_Int is = 0; is < ns; is++)
      {
        RELEASESHAREDS(objst[is], structF[is]);
        RELEASESHAREDS(objstp[is], structFp[is]);
      }
      for (E_Int nos = 0; nos < nu; nos++)
      {
        RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
        RELEASESHAREDU(objutp[nos], unstrFp[nos], cntp[nos]);
      }
      return NULL;
    } 
    posxp.push_back(posxi); posyp.push_back(posyi); poszp.push_back(poszi);
  }
  
  /*-------------------- Fin des verifs --------------------------------------*/
  E_Int sizemaxS = 0;
  E_Int sizemaxU = 0;
  for (E_Int noz = 0; noz < ns; noz++)
  { sizemaxS += structF[noz]->getSize(); }
  for (E_Int noz = 0; noz < nu; noz++)
  { sizemaxU += unstrF[noz]->getSize(); }
  E_Int sizemax = sizemaxS+sizemaxU;

  // Creation du kdtree des frontieres periodiques et des tableaux d indirection
  FldArrayF ftemp(sizemax,3);
  E_Float* xp = ftemp.begin(1);
  E_Float* yp = ftemp.begin(2);
  E_Float* zp = ftemp.begin(3);
  E_Int indt = 0;
  FldArrayI indirB(sizemax);
  FldArrayI indirI(sizemax);

  for (E_Int noz = 0; noz < ns; noz++)
  {
    E_Float* xt = structFp[noz]->begin(posxp[noz]+1);
    E_Float* yt = structFp[noz]->begin(posyp[noz]+1);
    E_Float* zt = structFp[noz]->begin(poszp[noz]+1);
    for (E_Int ind = 0; ind < structFp[noz]->getSize(); ind++)
    {
      xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
      indirB[indt] = noz; indirI[indt] = ind;
      indt++;
    }        
  }
  for (E_Int noz = 0; noz < nu; noz++)
  {
    E_Float* xt = unstrFp[noz]->begin(posxp[noz+ns]+1);
    E_Float* yt = unstrFp[noz]->begin(posyp[noz+ns]+1);
    E_Float* zt = unstrFp[noz]->begin(poszp[noz+ns]+1);
    for (E_Int ind = 0; ind < unstrFp[noz]->getSize(); ind++)
    {
      xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
      indirB[indt] = ns+noz; indirI[indt] = ind;
      indt++;
    }        
  }

  PyObject* l = PyList_New(0);
  vector<E_Float*> fields; vector<E_Int> sizeOfFields;
  for (E_Int noz = 0; noz < ns; noz++)
  {
    E_Int npts = structF[noz]->getSize();
    E_Int nfld = structF[noz]->getNfld();
    PyObject* tpl = K_ARRAY::buildArray(nfld, structVarString[noz],
                                        nit[noz], njt[noz], nkt[noz]);
    E_Float* fp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF ftemp0(npts,nfld, fp, true); ftemp0 = *structF[noz];
    fields.push_back(fp); sizeOfFields.push_back(npts);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  const char* eltType = "NGON";
  for (E_Int noz = 0; noz < nu; noz++)
  {
    E_Int* cnp = cnt[noz]->begin(); // pointeur sur la connectivite NGon
    E_Int sizeFN = cnp[1]; //  taille de la connectivite Face/Noeuds
    E_Int nelts = cnp[sizeFN+2];  // nombre total d elements
    E_Int nfld = unstrF[noz]->getNfld();
    E_Int npts = unstrF[noz]->getSize();
    E_Int sizecn = cnt[noz]->getSize();
    PyObject* tpl = K_ARRAY::buildArray(nfld, unstrVarString[noz],
                                        npts, nelts,-1, eltType, false, sizecn);
    E_Float* fp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF ftemp0(npts, nfld, fp, true); ftemp0 = *unstrF[noz];
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    FldArrayI cnn(sizecn, 1, cnnp, true); cnn = *cnt[noz];
    fields.push_back(fp); sizeOfFields.push_back(npts);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  ArrayAccessor<FldArrayF> posAcc(ftemp,1,2,3);
  K_SEARCH::KdTree<FldArrayF> globalKdtPeriodic(posAcc);
  E_Float pt[3];
  E_Float dist2;
  FldArrayIS tagDnr(sizemax); tagDnr.setAllValuesAtNull(); short* tagDnrp = tagDnr.begin();
  E_Float tolmatch2 = tolmatch*tolmatch;
  E_Float tolmatchg = 200.*tolmatch;
  E_Float tolmatchg2 = tolmatchg*tolmatchg;
  E_Float tolsurf = K_CONST::E_GEOM_CUTOFF; tolsurf = tolsurf*tolsurf;
  //1ere passe : exacte
  vector<E_Int> indicesR; vector<E_Int> indicesD;
  vector<E_Int> zonesR; vector<E_Int> zonesD;
 
  E_Int nzones = ns+nu;
  for (E_Int nozR = 0; nozR < nzones; nozR++)
  {
    E_Int nptsR = sizeOfFields[nozR];
    E_Float* xtR = fields[nozR]+posxt[nozR]*nptsR;
    E_Float* ytR = fields[nozR]+posyt[nozR]*nptsR;
    E_Float* ztR = fields[nozR]+poszt[nozR]*nptsR;
    E_Float* voltR = fields[nozR]+posvt[nozR]*nptsR;
    E_Float* tagR2 = fields[nozR]+post2[nozR]*nptsR;
    E_Float* tagR1 = fields[nozR]+post1[nozR]*nptsR;
    for (E_Int indR = 0; indR < nptsR; indR++)
    {
      if (tagR1[indR] == -1.)
      {
        pt[0] = xtR[indR]; pt[1] = ytR[indR]; pt[2] = ztR[indR];
        E_Int indgopp = globalKdtPeriodic.getClosest(pt, dist2); // closest pt
        E_Int nozD = indirB[indgopp]; 
        E_Int indD = indirI[indgopp];
        E_Int nptsD = sizeOfFields[nozD];
        if (tagDnrp[indgopp] == 0)
        {          
          E_Float* tagD1 = fields[nozD]+post1[nozD]*nptsD;
          E_Float* tagD2 = fields[nozD]+post2[nozD]*nptsD;
          E_Float* voltD = fields[nozD]+posvt[nozD]*nptsD;
          E_Float ecart_abs_vol = K_FUNC::E_abs(voltR[indR]-voltD[indD]);
          E_Float ecart_rel_vol = ecart_abs_vol/voltR[indR];
          ecart_rel_vol = K_FUNC::E_max(ecart_rel_vol, ecart_abs_vol/voltD[indD]);

          if (dist2 < tolmatch2)
          {
            if (ecart_rel_vol<0.05)
            {
              tagR1[indR] = nozD; tagR2[indR] = indD;        
              tagD1[indD] = nozR; tagD2[indD] = indR;        
              tagDnrp[indgopp] = 1;
            } 
          }
          else if (dist2 < tolmatchg2)
          {
            if (ecart_rel_vol<0.15)
            {
              indicesR.push_back(indR); indicesD.push_back(indD);
              zonesR.push_back(nozR); zonesD.push_back(nozD);
            }
          }
        }
      }
    }
  }
  tagDnr.malloc(0); indirB.malloc(0); indirI.malloc(0);
  E_Int sizeI = indicesR.size();
  if (sizeI > 0 )
  {
    //2eme passe : on est plus tolerant 
    for (E_Int i = 0; i < sizeI; i++)
    {
      E_Int indR = indicesR[i];//indice local a la zone    
      E_Int indD = indicesD[i];
      E_Int nozR = zonesR[i];
      E_Int nozD = zonesD[i]; 
      E_Int nptsR = sizeOfFields[nozR];
      E_Int nptsD = sizeOfFields[nozD];
      E_Float* tagR1 = fields[nozR]+post1[nozR]*nptsR;
      E_Float* tagR2 = fields[nozR]+post2[nozR]*nptsR;
      E_Float* tagD1 = fields[nozD]+post1[nozD]*nptsD;
      E_Float* tagD2 = fields[nozD]+post2[nozD]*nptsD;
      tagR1[indR] = nozD; tagR2[indR] = indD;        
      tagD1[indD] = nozR; tagD2[indD] = indR;        
    }
  }

  for (E_Int is = 0; is < ns; is++)
    RELEASESHAREDS(objst[is], structF[is]);
  for (E_Int nos = 0; nos < nu; nos++)
    RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
  return l;
}
//=============================================================================
/* Identify matching points in windows
   arrays listOfAllWins must contain first structured arrays, then NGON arrays
*/
//=============================================================================
PyObject* K_CONNECTOR::identifyMatching(PyObject* self, PyObject* args)
{
  PyObject *listOfAllWins;
  E_Float tolmatch;
  if (!PYPARSETUPLEF(args,
                    "Od", "Of",
                    &listOfAllWins, &tolmatch))
  {
      return NULL;
  }

  if (PyList_Check(listOfAllWins) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyMatching: 1st argument must be a list.");
    return NULL;
  }
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false;
  E_Boolean skipDiffVars = false;//true;

  E_Int isOk = K_ARRAY::getFromArrays(
    listOfAllWins, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltTypet, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int ns = structF.size(); E_Int nu = unstrF.size();

  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError, "identifyMatching: 1st list of arrays is not valid.");
    for (E_Int is = 0; is < ns; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    return NULL;
  } 
  if (ns+nu == 0) 
  {
    PyErr_SetString(PyExc_TypeError, "identifyMatching: 1st list of arrays does not contain valid zones.");
    for (E_Int is = 0; is < ns; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    return NULL;
  }
  vector<E_Int> posxs; vector<E_Int> posys; vector<E_Int> poszs; 
  vector<E_Int> posxu; vector<E_Int> posyu; vector<E_Int> poszu; 
  vector<E_Int> post1; vector<E_Int> post2;
  E_Int posxi, posyi, poszi, postag1, postag2;
  for (E_Int i = 0; i < ns; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(structVarString[i]); posxi++; 
    posyi = K_ARRAY::isCoordinateYPresent(structVarString[i]); posyi++;
    poszi = K_ARRAY::isCoordinateZPresent(structVarString[i]); poszi++;
    postag1 = K_ARRAY::isNamePresent("tag1", structVarString[i]); postag1++;
    postag2 = K_ARRAY::isNamePresent("tag2", structVarString[i]); postag2++;
    if (postag1 < 1 || postag2 < 1) 
    {
      PyErr_SetString(PyExc_TypeError, "identifyMatching: tag1 and tag2 must be defined in listOfAllWins.");
      for (E_Int is = 0; is < ns; is++)
        RELEASESHAREDS(objst[is], structF[is]);
      for (E_Int nos = 0; nos < nu; nos++)
        RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
      return NULL;
    } 
    posxs.push_back(posxi); posys.push_back(posyi); poszs.push_back(poszi);
    post1.push_back(postag1); post2.push_back(postag2); 
  }
  for (E_Int i = 0; i < nu; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(unstrVarString[i]); posxi++; 
    posyi = K_ARRAY::isCoordinateYPresent(unstrVarString[i]); posyi++;
    poszi = K_ARRAY::isCoordinateZPresent(unstrVarString[i]); poszi++;
    postag1 = K_ARRAY::isNamePresent("tag1", unstrVarString[i]); postag1++;
    postag2 = K_ARRAY::isNamePresent("tag2", unstrVarString[i]); postag2++;
    if (postag1 < 1 || postag2 < 1) 
    {
      PyErr_SetString(PyExc_TypeError, "identifyMatching: tag1 and tag2 must be defined in listOfAllWins.");
      for (E_Int is = 0; is < ns; is++)
        RELEASESHAREDS(objst[is], structF[is]);
      for (E_Int nos = 0; nos < nu; nos++)
        RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
      return NULL;
    } 
    if (strcmp(eltTypet[i],"NGON") != 0 )
    {
      PyErr_SetString(PyExc_TypeError, "identifyMatching: unstructured arrays must be NGON.");
      for (E_Int is = 0; is < ns; is++)
        RELEASESHAREDS(objst[is], structF[is]);
      for (E_Int nos = 0; nos < nu; nos++)
        RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
      return NULL;
    } 
    posxu.push_back(posxi); posyu.push_back(posyi); poszu.push_back(poszi);
    post1.push_back(postag1); post2.push_back(postag2); 
  }
  /*-------------------- Fin des verifs --------------------------------------*/
  E_Int sizemax = 0;
  for (E_Int noz = 0; noz < ns; noz++)
  { sizemax += structF[noz]->getSize(); }
  for (E_Int noz = 0; noz < nu; noz++)
  { sizemax += unstrF[noz]->getSize(); }

  // Creation du kdtree et des tableaux d indirection
  vector<E_Int> newID(sizemax);
  FldArrayF ftemp(sizemax,3);
  E_Float* xp = ftemp.begin(1);
  E_Float* yp = ftemp.begin(2);
  E_Float* zp = ftemp.begin(3);
  E_Int indt = 0;
  FldArrayI indirB(sizemax);
  FldArrayI indirI(sizemax);
  
  for (E_Int noz = 0; noz < ns; noz++)
  {
    E_Float* xt = structF[noz]->begin(posxs[noz]);
    E_Float* yt = structF[noz]->begin(posys[noz]);
    E_Float* zt = structF[noz]->begin(poszs[noz]);
    for (E_Int ind = 0; ind < structF[noz]->getSize(); ind++)
    {
      xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
      indirB[indt] = noz; indirI[indt] = ind;
      indt++;
    }        
  }
  for (E_Int noz = 0; noz < nu; noz++)
  {
    E_Float* xt = unstrF[noz]->begin(posxu[noz]);
    E_Float* yt = unstrF[noz]->begin(posyu[noz]);
    E_Float* zt = unstrF[noz]->begin(poszu[noz]);
    for (E_Int ind = 0; ind < unstrF[noz]->getSize(); ind++)
    {
      xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
      indirB[indt] = ns+noz; indirI[indt] = ind;
      indt++;
    }        
  }

  PyObject* l = PyList_New(0);
  vector<E_Float*> fields; vector<E_Int> sizeOfFields;
  for (E_Int noz = 0; noz < ns; noz++)
  {
    E_Int npts = structF[noz]->getSize();
    E_Int nfld = structF[noz]->getNfld();
    PyObject* tpl = K_ARRAY::buildArray(nfld, structVarString[noz],
                                        nit[noz], njt[noz], nkt[noz]);
    E_Float* fp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF ftemp0(npts, nfld, fp, true); ftemp0 = *structF[noz];
    fields.push_back(fp); sizeOfFields.push_back(npts);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  const char* eltType = "NGON";

  for (E_Int noz = 0; noz < nu; noz++)
  {
    E_Int* cnp = cnt[noz]->begin(); // pointeur sur la connectivite NGon
    E_Int sizeFN = cnp[1]; //  taille de la connectivite Face/Noeuds
    E_Int nelts = cnp[sizeFN+2];  // nombre total d elements
    E_Int nfld = unstrF[noz]->getNfld();
    E_Int npts = unstrF[noz]->getSize();
    E_Int sizecn = cnt[noz]->getSize();
    PyObject* tpl = K_ARRAY::buildArray(nfld, unstrVarString[noz],
                                        npts, nelts,-1, eltType, false, sizecn);
    E_Float* fp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF ftemp0(npts, nfld, fp, true); ftemp0 = *unstrF[noz];
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    FldArrayI cnn(sizecn, 1, cnnp, true); cnn = *cnt[noz];
    fields.push_back(fp); sizeOfFields.push_back(npts);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }

  ArrayAccessor<FldArrayF> posAcc(ftemp,1,2,3);
  E_Int nbMerges = merge(posAcc, tolmatch, newID);
  if (nbMerges > 0)
  {
    for (E_Int indt = 0; indt < sizemax; indt++)
    {
      E_Int indtopp = newID[indt];
      if (indtopp != indt)
      {
        E_Int noz = indirB[indt]; E_Int nozopp = indirB[indtopp];
        E_Int ind = indirI[indt]; E_Int indopp = indirI[indtopp];   
        // update receiver tags
        E_Int nptsz = sizeOfFields[noz];
        E_Float* tag1 = fields[noz]+(post1[noz]-1)*nptsz;
        E_Float* tag2 = fields[noz]+(post2[noz]-1)*nptsz;
        E_Int nptszopp = sizeOfFields[nozopp];
        E_Float* tago1 = fields[nozopp]+(post1[nozopp]-1)*nptszopp;
        E_Float* tago2 = fields[nozopp]+(post2[nozopp]-1)*nptszopp;
        if (tag1[ind] == -2. && tago1[indopp] == -2.)
          continue;
        
        // if (tag1[ind] == -1. && tago1[indopp] == -1.)
        if (tag1[ind] <0. && tago1[indopp] <0.)
        {
          tag1[ind] = nozopp; tag2[ind] = indopp;
          // update donor tags
          tago1[indopp] = noz; tago2[indopp] = ind;     
        }
      }      
    }
  }
  for (E_Int is = 0; is < ns; is++)
    RELEASESHAREDS(objst[is], structF[is]);
  for (E_Int nos = 0; nos < nu; nos++)
    RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
  return l;
}

//=============================================================================
/* Identify matching points in windows - used in connectNearMatch
   listOfAllWinsR must define oneovern windows
   listOfAllWinsD must define real windows  */
//=============================================================================
PyObject* K_CONNECTOR::identifyMatchingNM(PyObject* self, PyObject* args)
{
  PyObject *listOfAllWinsR, *listOfAllWinsD;
  E_Float tolmatch;
  if (!PYPARSETUPLEF(args,
                    "OOd", "OOf",
                    &listOfAllWinsR, &listOfAllWinsD, &tolmatch))
  {
      return NULL;
  }

  if (PyList_Check(listOfAllWinsR) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyMatching: 1st argument must be a list.");
    return NULL;
  }
  if (PyList_Check(listOfAllWinsD) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "identifyMatching: 2nd argument must be a list.");
    return NULL;
  }
  vector<E_Int> reslR;
  vector<char*> structVarStringR; vector<char*> unstrVarStringR;
  vector<FldArrayF*> structFR; vector<FldArrayF*> unstrFR;
  vector<E_Int> nitR; vector<E_Int> njtR; vector<E_Int> nktR;
  vector<FldArrayI*> cntR;
  vector<char*> eltTypetR;
  vector<PyObject*> objstR, objutR;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = true;
  E_Boolean skipDiffVars = false;//true;

  E_Int isOk = K_ARRAY::getFromArrays(
    listOfAllWinsR, reslR, structVarStringR, unstrVarStringR,
    structFR, unstrFR, nitR, njtR, nktR, cntR, eltTypetR, objstR, objutR, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nzonesR = structFR.size();
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError, "identifyMatching: 1st list of arrays is not valid.");
    for (E_Int is = 0; is < nzonesR; is++)
      RELEASESHAREDS(objstR[is], structFR[is]);
    return NULL;
  } 
  if (nzonesR == 0) 
  {
    PyErr_SetString(PyExc_TypeError, "identifyMatching: 1st list of arrays does not contain valid structured zones.");
    for (E_Int is = 0; is < nzonesR; is++)
      RELEASESHAREDS(objstR[is], structFR[is]);
    return NULL;
  }
  vector<E_Int> posxtR; vector<E_Int> posytR; vector<E_Int> posztR; 
  vector<E_Int> postR1; vector<E_Int> postR2;
  E_Int posxi, posyi, poszi, postag1, postag2;
  for (E_Int i = 0; i < nzonesR; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(structVarStringR[i]); posxi++; 
    posyi = K_ARRAY::isCoordinateYPresent(structVarStringR[i]); posyi++;
    poszi = K_ARRAY::isCoordinateZPresent(structVarStringR[i]); poszi++;
    postag1 = K_ARRAY::isNamePresent("tag1", structVarStringR[i]); postag1++;
    postag2 = K_ARRAY::isNamePresent("tag2", structVarStringR[i]); postag2++;
    if (postag1 < 1 || postag2 < 1) 
    {
      PyErr_SetString(PyExc_TypeError, "identifyMatching: tag1 and tag2 must be defined in listOfAllWins.");
      for (E_Int is = 0; is < nzonesR; is++)
        RELEASESHAREDS(objstR[is], structFR[is]);
      return NULL;
    } 
    posxtR.push_back(posxi); posytR.push_back(posyi); posztR.push_back(poszi);
    postR1.push_back(postag1); postR2.push_back(postag2); 
  }

  vector<E_Int> reslD;
  vector<char*> structVarStringD; vector<char*> unstrVarStringD;
  vector<FldArrayF*> structFD; vector<FldArrayF*> unstrFD;
  vector<E_Int> nitD; vector<E_Int> njtD; vector<E_Int> nktD;
  vector<FldArrayI*> cntD;
  vector<char*> eltTypetD;
  vector<PyObject*> objstD, objutD;
  
  isOk = K_ARRAY::getFromArrays(
    listOfAllWinsD, reslD, structVarStringD, unstrVarStringD,
    structFD, unstrFD, nitD, njtD, nktD, cntD, eltTypetD, objstD, objutD, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nzonesD = structFD.size();
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError, "identifyMatching: 2nd list of arrays is not valid.");
    for (E_Int is = 0; is < nzonesR; is++)
      RELEASESHAREDS(objstR[is], structFR[is]);
    for (E_Int is = 0; is < nzonesD; is++)
      RELEASESHAREDS(objstD[is], structFD[is]);
    return NULL;
  } 
  if (nzonesD == 0) 
  {
    PyErr_SetString(PyExc_TypeError, "identifyMatching: 2nd list of arrays does not contain valid structured zones.");
    for (E_Int is = 0; is < nzonesR; is++)
      RELEASESHAREDS(objstR[is], structFR[is]);
    for (E_Int is = 0; is < nzonesD; is++)
      RELEASESHAREDS(objstD[is], structFD[is]);
    return NULL;
  }
  if (nzonesD != nzonesR)
  {
    PyErr_SetString(PyExc_TypeError, "identifyMatching: 1st and 2nd list must be of same size.");
    for (E_Int is = 0; is < nzonesR; is++)
      RELEASESHAREDS(objstR[is], structFR[is]);
    for (E_Int is = 0; is < nzonesD; is++)
      RELEASESHAREDS(objstD[is], structFD[is]);
    return NULL;
  }
  vector<E_Int> posxtD; vector<E_Int> posytD; vector<E_Int> posztD; 
  vector<E_Int> postD1; vector<E_Int> postD2;
  for (E_Int i = 0; i < nzonesD; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(structVarStringD[i]); posxi++; 
    posyi = K_ARRAY::isCoordinateYPresent(structVarStringD[i]); posyi++;
    poszi = K_ARRAY::isCoordinateZPresent(structVarStringD[i]); poszi++;
    postag1 = K_ARRAY::isNamePresent("tag1", structVarStringD[i]); postag1++;
    postag2 = K_ARRAY::isNamePresent("tag2", structVarStringD[i]); postag2++;
    if (postag1 < 1 || postag2 < 1) 
    {
      PyErr_SetString(PyExc_TypeError, "identifyMatching: tag1 and tag2 must be defined in listOfAllWinsD.");
      for (E_Int is = 0; is < nzonesR; is++)
      {
        RELEASESHAREDS(objstR[is], structFR[is]);
        RELEASESHAREDS(objstD[is], structFD[is]);
      }
      return NULL;
    } 
    posxtD.push_back(posxi); posytD.push_back(posyi); posztD.push_back(poszi);
    postD1.push_back(postag1); postD2.push_back(postag2); 
  }
  /*-------------------- Fin des verifs --------------------------------------*/
  E_Int sizemaxR = 0;
  for (E_Int noz = 0; noz < nzonesR; noz++)
  { sizemaxR += structFR[noz]->getSize(); }
  
  // Creation du kdtree et des tableaux d indirection
  FldArrayF ftempR(sizemaxR,3);
  E_Float* xp = ftempR.begin(1);
  E_Float* yp = ftempR.begin(2);
  E_Float* zp = ftempR.begin(3);
  FldArrayI indirBR(sizemaxR);// original block number 
  FldArrayI indirIR(sizemaxR);// index in original block
  E_Int indt = 0;
  for (E_Int noz = 0; noz < nzonesR; noz++)
  {
    E_Float* xt = structFR[noz]->begin(posxtR[noz]);
    E_Float* yt = structFR[noz]->begin(posytR[noz]);
    E_Float* zt = structFR[noz]->begin(posztR[noz]);
    for (E_Int ind = 0; ind < structFR[noz]->getSize(); ind++)
    {
      xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
      indirBR[indt] = noz; indirIR[indt] = ind;
      indt++;
    }        
  }
  ArrayAccessor<FldArrayF> coordAcc(ftempR, 1, 2, 3);
  K_SEARCH::KdTree<FldArrayF> globalKdt(coordAcc);
  PyObject* l = PyList_New(0);
  vector<E_Float*> fields;
  for (E_Int noz = 0; noz < nzonesR; noz++)
  {
    PyObject* tpl = K_ARRAY::buildArray(structFR[noz]->getNfld(), structVarStringR[noz],
                                        nitR[noz], njtR[noz], nktR[noz]);
    E_Float* fp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF ftemp0(structFR[noz]->getSize(),structFR[noz]->getNfld(), fp, true); ftemp0 = *structFR[noz];
    fields.push_back(fp);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  for (E_Int noz = 0; noz < nzonesD; noz++)
  {
    PyObject* tpl = K_ARRAY::buildArray(structFD[noz]->getNfld(), structVarStringD[noz],
                                        nitD[noz], njtD[noz], nktD[noz]);
    E_Float* fp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF ftemp0(structFD[noz]->getSize(),structFD[noz]->getNfld(), fp, true); ftemp0 = *structFD[noz];
    fields.push_back(fp);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  E_Float eps = 10*tolmatch;
  E_Float bbmin[3], bbmax[3];
  vector<E_Int> candidates;
  for (E_Int nozd = 0; nozd < nzonesD; nozd++)
  {
    E_Int nptszD = structFD[nozd]->getSize();
    E_Float* xdt = fields[nozd+nzonesR]+(posxtD[nozd]-1)*nptszD;
    E_Float* ydt = fields[nozd+nzonesR]+(posytD[nozd]-1)*nptszD;
    E_Float* zdt = fields[nozd+nzonesR]+(posztD[nozd]-1)*nptszD;
    E_Float* tagD1 = fields[nozd+nzonesR]+(postD1[nozd]-1)*nptszD;
    E_Float* tagD2 = fields[nozd+nzonesR]+(postD2[nozd]-1)*nptszD;

    for (E_Int indd = 0; indd < structFD[nozd]->getSize(); indd++)
    {
      bbmin[0] = xdt[indd]-eps; bbmin[1] = ydt[indd]-eps; bbmin[2] = zdt[indd]-eps;
      bbmax[0] = xdt[indd]+eps; bbmax[1] = ydt[indd]+eps; bbmax[2] = zdt[indd]+eps;
      candidates.clear();
      globalKdt.getInBox(bbmin,bbmax,candidates);
      for (size_t cr = 0; cr < candidates.size(); cr++)
      {
        E_Int indkr = candidates[cr];
        E_Int nozr = indirBR[indkr]; E_Int indr = indirIR[indkr];
        if (nozr != nozd) 
        {
          E_Int nptszR = structFR[nozr]->getSize();
          E_Float* tagR1 = fields[nozr]+(postR1[nozr]-1)*nptszR;
          E_Float* tagR2 = fields[nozr]+(postR2[nozr]-1)*nptszR;
          E_Float* xrt = fields[nozr]+(posxtR[nozr]-1)*nptszR;
          E_Float* yrt = fields[nozr]+(posytR[nozr]-1)*nptszR;
          E_Float* zrt = fields[nozr]+(posztR[nozr]-1)*nptszR;
          E_Float dx = xrt[indr]-xdt[indd];
          E_Float dy = yrt[indr]-ydt[indd];
          E_Float dz = zrt[indr]-zdt[indd];
          if (K_FUNC::fEqualZero(dx,tolmatch) == true &&
              K_FUNC::fEqualZero(dy,tolmatch) == true &&
              K_FUNC::fEqualZero(dz,tolmatch) == true && tagD1[indd] == -1.)// && tagR1[indr] == -1.)
          {
            tagR1[indr] = nozd+nzonesR; tagR2[indr] = indd;
            // update donor tags
            tagD1[indd] = nozr; tagD2[indd] = indr;
          }
        }
      }
    }
  }
  for (E_Int is = 0; is < nzonesR; is++)
    RELEASESHAREDS(objstR[is], structFR[is]);
  for (E_Int is = 0; is < nzonesD; is++)
    RELEASESHAREDS(objstD[is], structFD[is]);
  return l;
}
