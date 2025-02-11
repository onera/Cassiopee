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

# include "dist2walls.h"

using namespace std;
using namespace K_FLD;
using namespace K_SEARCH;

//=============================================================================
/* Calcul de la distance a la paroi par interface de distance minimum signee
   les maillages  sont deja en centres ici */
//=============================================================================
PyObject* K_DIST2WALLS::distance2WallsSigned(PyObject* self, PyObject* args)
{
  PyObject *blks, *bodiesC;
  if (!PyArg_ParseTuple(args, "OO", &blks, &bodiesC)) return NULL;
  if (PyList_Check(blks) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "distance2WallsSigned: 1st argument must be a list.");
    return NULL;
  }

  if (PyList_Check(bodiesC) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "distance2WallsSigned: 2nd argument must be a list.");
    return NULL;
  }
  
  // Maillage en centres
  // les coordonnees doivent etre en premier
  vector<PyObject*> objst, objut;
  vector<E_Int> res0;
  vector<char*> structVarString0;
  vector<char*> unstrVarString0;
  vector<FldArrayF*> structF0;
  vector<FldArrayF*> unstrF0;
  vector<E_Int> nit0;
  vector<E_Int> njt0; 
  vector<E_Int> nkt0;
  vector<FldArrayI*> cnt0;
  vector<char*> eltType0;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false;
  E_Boolean skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    blks, res0, structVarString0, unstrVarString0,
    structF0, unstrF0, nit0, njt0, nkt0, cnt0, eltType0, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int ns0 = structF0.size(); E_Int nu0 = unstrF0.size();
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "distance2WallsSigned: invalid list of mesh arrays.");
    for (E_Int nos = 0; nos < ns0; nos++)
      RELEASESHAREDS(objst[nos], structF0[nos]);
    for (E_Int nos = 0; nos < nu0; nos++)
      RELEASESHAREDU(objut[nos], unstrF0[nos], cnt0[nos]);
    return NULL;
  }
  E_Int posx=-1, posy=-1, posz=-1; 
  vector<E_Int> ncellss; vector<E_Int> ncellsu; 
  // Verification de posxi, posyi, poszi dans listFields: vars 1,2,3 imposees
  for (E_Int i = 0; i < ns0; i++)
  {
    E_Int posxi = K_ARRAY::isCoordinateXPresent(structVarString0[i]);
    E_Int posyi = K_ARRAY::isCoordinateYPresent(structVarString0[i]);
    E_Int poszi = K_ARRAY::isCoordinateZPresent(structVarString0[i]);
    posxi++; posyi++; poszi++; 
    if (posx == -1) posx = posxi;
    if (posy == -1) posy = posyi;
    if (posz == -1) posz = poszi;

    if (posxi != posx || posyi != posy || poszi != posz) 
    {
      PyErr_SetString( PyExc_TypeError,
        "distance2WallsSigned: coordinates must be located at same position for all zones.");
      for (E_Int nos = 0; nos < ns0; nos++)
        RELEASESHAREDS(objst[nos], structF0[nos]);
      for (E_Int nos = 0; nos < nu0; nos++)
        RELEASESHAREDU(objut[nos], unstrF0[nos], cnt0[nos]);
      return NULL;
    }
    ncellss.push_back(structF0[i]->getSize());
  }
  for (E_Int i = 0; i < nu0; i++)
  {
    E_Int posxi = K_ARRAY::isCoordinateXPresent(unstrVarString0[i]);
    E_Int posyi = K_ARRAY::isCoordinateYPresent(unstrVarString0[i]);
    E_Int poszi = K_ARRAY::isCoordinateZPresent(unstrVarString0[i]);
    posxi++; posyi++; poszi++;
    if (posx == -1) posx = posxi;
    if (posy == -1) posy = posyi;
    if (posz == -1) posz = poszi;

    if (posxi != posx || posyi != posy || poszi != posz) 
    {
      PyErr_SetString( PyExc_TypeError,
        "distance2WallsSigned: coordinates must be located at same position for all zones.");
      for (E_Int nos = 0; nos < ns0; nos++)
        RELEASESHAREDS(objst[nos], structF0[nos]);
      for (E_Int nos = 0; nos < nu0; nos++)
        RELEASESHAREDU(objut[nos], unstrF0[nos], cnt0[nos]);
      return NULL;
    }
    ncellsu.push_back(unstrF0[i]->getSize());
  }
  // Extract infos from body surfaces + cellN + sx, sy, sz
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypeb;
  vector<PyObject*> objs, obju;
  skipNoCoord = true;
  skipStructured = true;
  skipUnstructured = false; 
  skipDiffVars = true;
  K_ARRAY::getFromArrays(
    bodiesC, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltTypeb, objs, obju, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nwalls = unstrF.size();

  if (nwalls == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "distance2WallsSigned: invalid list of blocks.");
    for (E_Int nos = 0; nos < ns0; nos++)
      RELEASESHAREDS(objst[nos], structF0[nos]);
    for (E_Int nos = 0; nos < nu0; nos++)
      RELEASESHAREDU(objut[nos], unstrF0[nos], cnt0[nos]);
    return NULL;
  }
  
  // verification de posxv, posyv, poszv dans bodiesC  
  vector<E_Int> posxv;  vector<E_Int> posyv;  vector<E_Int> poszv;
  vector<E_Int> poscv;
  E_Int possx = K_ARRAY::isNamePresent("sx",unstrVarString[0]); possx++;
  E_Int possy = K_ARRAY::isNamePresent("sy",unstrVarString[0]); possy++;
  E_Int possz = K_ARRAY::isNamePresent("sz",unstrVarString[0]); possz++;
  for (E_Int v = 0; v < nwalls; v++)
  {
    E_Int posxv0 = K_ARRAY::isCoordinateXPresent(unstrVarString[v]);
    E_Int posyv0 = K_ARRAY::isCoordinateYPresent(unstrVarString[v]);
    E_Int poszv0 = K_ARRAY::isCoordinateZPresent(unstrVarString[v]);
    E_Int poscv0 = K_ARRAY::isCellNatureField2Present(unstrVarString[v]);
    posxv0++; posxv.push_back(posxv0);
    posyv0++; posyv.push_back(posyv0);
    poszv0++; poszv.push_back(poszv0);
    poscv0++; poscv.push_back(poscv0);
    if (poscv0 == 0) 
    {
      PyErr_SetString(PyExc_TypeError,"distance2WallsSigned: cellN must be defined for bodies.");
      for (E_Int nos = 0; nos < ns0; nos++)
        RELEASESHAREDS(objst[nos], structF0[nos]);
      for (E_Int nos = 0; nos < nu0; nos++)
        RELEASESHAREDU(objut[nos], unstrF0[nos], cnt0[nos]);
      for (E_Int nos = 0; nos < nwalls; nos++)
        RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
      return NULL;
    }
    if (K_ARRAY::isNamePresent("sx",unstrVarString[v])+1 != possx ||
        K_ARRAY::isNamePresent("sy",unstrVarString[v])+1 != possy ||
        K_ARRAY::isNamePresent("sz",unstrVarString[v])+1 != possz)
    {
      PyErr_SetString(PyExc_TypeError,"distance2WallsSigned: normals vectors must be located at the same place.");
      for (E_Int nos = 0; nos < ns0; nos++)
        RELEASESHAREDS(objst[nos], structF0[nos]);
      for (E_Int nos = 0; nos < nu0; nos++)
        RELEASESHAREDU(objut[nos], unstrF0[nos], cnt0[nos]);
      for (E_Int nos = 0; nos < nwalls; nos++)
        RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
      return NULL;
    }
  }
  
  // Calcul de la distance a la paroi
  vector<FldArrayF*> distances;
  for (E_Int i = 0; i < ns0; i++)
  {
    E_Int ncells = ncellss[i];
    FldArrayF* distance = new FldArrayF(ncells); 
    distance->setAllValuesAt(K_CONST::E_INFINITE);
    distances.push_back(distance);
  }
  if (structF0.size() > 0) 
    computeMininterfSigned(ncellss, posx, posy, posz, 
                           structF0, posxv, posyv, poszv, poscv, unstrF, 
                           possx, possy, possz, 
                           distances);
  vector<FldArrayF*> distancesu;
  for (E_Int i = 0; i < nu0; i++)
  {
    E_Int ncells = ncellsu[i];
    FldArrayF* distance = new FldArrayF(ncells); 
    distance->setAllValuesAt(K_CONST::E_INFINITE);
    distancesu.push_back(distance);
  }
  if (unstrF0.size() > 0) 
    computeMininterfSigned(ncellsu, posx, posy, posz,
                           unstrF0, posxv, posyv, poszv, poscv, unstrF, 
                           possx, possy, possz, 
                           distancesu);

  for (E_Int nos = 0; nos < nwalls; nos++)
    RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
  // Build arrays
  PyObject* l = PyList_New(0);
  PyObject* tpl;    
  for (E_Int nos = 0; nos < ns0; nos++)
  {
    tpl = K_ARRAY::buildArray(*distances[nos], "TurbulentDistance", 
                              nit0[nos], njt0[nos], nkt0[nos]);
    PyList_Append(l, tpl); Py_DECREF(tpl);
    delete distances[nos];
    RELEASESHAREDS(objst[nos], structF0[nos]);
  }
  for (E_Int nou = 0; nou < nu0; nou++)
  {
    K_FLD::FldArrayI* cnout = new K_FLD::FldArrayI(*cnt0[nou]);
    tpl = K_ARRAY::buildArray(*distancesu[nou], "TurbulentDistance", *cnout,
                              -1, eltType0[nou]);
    PyList_Append(l, tpl); Py_DECREF(tpl);
    RELEASESHAREDU(objut[nou], unstrF0[nou],cnt0[nou]);
    delete distancesu[nou]; delete cnout;
  }
  return l;
}

//=============================================================================
/* Calcul de la distance a la paroi */
//=============================================================================
void K_DIST2WALLS::computeMininterfSigned(
  vector<E_Int>& ncellst,
  E_Int posx, E_Int posy, E_Int posz, 
  vector<FldArrayF*>& fields, 
  vector<E_Int>& posxv, vector<E_Int>& posyv, vector<E_Int>& poszv, 
  vector<E_Int>& poscv, vector<FldArrayF*>& fieldsw, 
  E_Int possx, E_Int possy, E_Int possz,
  vector<FldArrayF*>& distances)
{
  /*1 - creation du kdtree*/
  E_Int nwalls = fieldsw.size();
  E_Int nptsmax = 0;
  for (E_Int v = 0; v < nwalls; v++)
  {
    FldArrayF* fieldv = fieldsw[v];
    nptsmax = nptsmax + fieldv->getSize();
  }
  FldArrayF* wallpts = new FldArrayF(nptsmax,6);
  E_Float* xw2 = wallpts->begin(1);
  E_Float* yw2 = wallpts->begin(2);
  E_Float* zw2 = wallpts->begin(3);
  E_Float* sxw2 = wallpts->begin(4);
  E_Float* syw2 = wallpts->begin(5);
  E_Float* szw2 = wallpts->begin(6);
  
  E_Int c = 0; E_Int nzones = fields.size();
  for (E_Int v = 0; v < nwalls; v++)
  {
    FldArrayF* fieldv = fieldsw[v];
    E_Int posxw = posxv[v]; E_Int posyw = posyv[v]; E_Int poszw = poszv[v];
    E_Float* xw = fieldv->begin(posxw);
    E_Float* yw = fieldv->begin(posyw);
    E_Float* zw = fieldv->begin(poszw);
    E_Float* sxw = fieldv->begin(possx);
    E_Float* syw = fieldv->begin(possy);
    E_Float* szw = fieldv->begin(possz);
    E_Int ncellsw = fieldv->getSize();
    /* recuperation des points calcules uniquement 
       pas de pts masques et interpoles dans kdtree */
    E_Int poscw = poscv[v]; E_Float* cellnw0 = fieldv->begin(poscw);
    for (E_Int i = 0; i < ncellsw; i++)
    {
      if (cellnw0[i] == 1.) 
      {
        xw2[c] = xw[i]; yw2[c] = yw[i]; zw2[c] = zw[i]; 
        sxw2[c] = sxw[i]; syw2[c] = syw[i]; szw2[c] = szw[i]; 
        c++;
      }
    }
  }//fin kdtree
  if (c == 0) 
  {
    delete wallpts;
    for (E_Int v = 0; v < nzones; v++) 
    {
      E_Int ncells = ncellst[v];
      E_Float* distancep = distances[v]->begin();
      for (E_Int ind = 0; ind < ncells; ind++)
        distancep[ind] = sqrt(distancep[ind]); 
    }
    return;
  }
  wallpts->reAllocMat(c, 6);
  xw2 = wallpts->begin(1); yw2 = wallpts->begin(2); zw2 = wallpts->begin(3);
  sxw2 = wallpts->begin(4); syw2 = wallpts->begin(5); szw2 = wallpts->begin(6);
  
  ArrayAccessor<FldArrayF> coordAcc(*wallpts, 1,2,3);
  KdTree<FldArrayF> kdt(coordAcc, E_EPSILON);
    
  /* detection de la paroi la plus proche */
  E_Float pt[3];
  for (E_Int v = 0; v < nzones; v++)
  {
    E_Int ncells = ncellst[v];
    E_Float* xt = fields[v]->begin(posx);
    E_Float* yt = fields[v]->begin(posy);
    E_Float* zt = fields[v]->begin(posz);
    E_Float* distancep = distances[v]->begin();

#pragma omp parallel for default(shared) private(pt) schedule(dynamic)
    for (E_Int ind = 0; ind < ncells; ind++)
    {
      E_Float sgn; E_Int indw2;
      E_Float dist, rx, ry, rz, rad, sx, sy, sz;
      pt[0] = xt[ind]; pt[1] = yt[ind]; pt[2] = zt[ind];
      indw2 = kdt.getClosest(pt);
      rx = xw2[indw2]-pt[0]; ry = yw2[indw2]-pt[1]; rz = zw2[indw2]-pt[2];
      sx = sxw2[indw2]; sy = syw2[indw2]; sz = szw2[indw2]; 
      dist = rx*rx + ry*ry + rz*rz; rad = sqrt(dist);
      sgn = rx*sx + ry*sy + rz*sz;
      if (sgn < 0.) distancep[ind] = rad;
      else distancep[ind] = -rad;
    } // fin boucle

  }
  delete wallpts;
}
//=============================================================================
