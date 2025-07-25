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
/* Calcul de la distance a la paroi par interface de distance minimum
   les maillages  sont deja en centres ici */
//=============================================================================
PyObject* K_DIST2WALLS::distance2Walls(PyObject* self, PyObject* args)
{
  PyObject *blks, *bodiesC;
  E_Int isminortho;
  if (!PYPARSETUPLE_(args, OO_ I_, &blks, &bodiesC, &isminortho)) return NULL;
  
  if (PyList_Check(blks) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "dist2Walls: 1st argument must be a list.");
    return NULL;
  }

  if (PyList_Check(bodiesC) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "dist2Walls: 2nd argument must be a list.");
    return NULL;
  }
  
  // Maillage: les coordonnees doivent etre en premier
  vector<PyObject*> objst, objut;
  vector<char*> structVarString0; vector<char*> unstrVarString0;
  vector<FldArrayF*> structF0; vector<FldArrayF*> unstrF0;
  vector<E_Int> nit0; vector<E_Int> njt0; vector<E_Int> nkt0;
  vector<FldArrayI*> cnt0;
  vector<char*> eltType0;
  char* varStringl; char* eltTypel;
  E_Int nil, njl, nkl;
  FldArrayF* fl; FldArrayI* cnl;
  PyObject* o; E_Int res;
  E_Int nz = PyList_Size(blks);
  for (E_Int i = 0; i < nz; i++)
  {
    o = PyList_GetItem(blks, i);
    res = K_ARRAY::getFromArray3(o, varStringl, fl, 
                                 nil, njl, nkl, cnl, eltTypel);
    if (res == 1)
    {
      structVarString0.push_back(varStringl);
      structF0.push_back(fl);
      nit0.push_back(nil); njt0.push_back(njl), nkt0.push_back(nkl);
      objst.push_back(o);
    }
    else if (res == 2)
    {
      unstrVarString0.push_back(varStringl);
      unstrF0.push_back(fl);
      eltType0.push_back(eltTypel); cnt0.push_back(cnl);
      objut.push_back(o);
    }
    else printf("Warning: dist2Walls: array " SF_D_ " is invalid. Discarded.\n", i);
  }

  E_Int ns0 = structF0.size(); E_Int nu0 = unstrF0.size();
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
        "dist2Walls: coordinates must be located at same position for all zones.");
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
        "dist2Walls: coordinates must be located at same position for all zones.");
      for (E_Int nos = 0; nos < ns0; nos++)
        RELEASESHAREDS(objst[nos], structF0[nos]);
      for (E_Int nos = 0; nos < nu0; nos++)
        RELEASESHAREDU(objut[nos], unstrF0[nos], cnt0[nos]);
      return NULL;
    }
    ncellsu.push_back(unstrF0[i]->getSize());
  }
  // Extract infos from body surfaces + cellN
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypeb;
  vector<PyObject*> objs, obju;  
  nz = PyList_Size(bodiesC);
  for (E_Int i = 0; i < nz; i++)
  {
    o = PyList_GetItem(bodiesC, i);
    res = K_ARRAY::getFromArray3(o, varStringl, fl, 
                                 nil, njl, nkl, cnl, eltTypel);
    if (res == 1)
    {
      structVarString.push_back(varStringl);
      structF.push_back(fl);
      nit.push_back(nil); njt.push_back(njl), nkt.push_back(nkl);
      objs.push_back(o);
    }
    else if (res == 2)
    {
      unstrVarString.push_back(varStringl);
      unstrF.push_back(fl);
      eltTypeb.push_back(eltTypel); cnt.push_back(cnl);
      obju.push_back(o);
    }
    else printf("Warning: dist2Walls: array " SF_D_ " is invalid. Discarded.\n", i);
  }

  E_Int nwalls = unstrF.size();
  if (nwalls == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "dist2Walls: invalid list of blocks.");
    for (E_Int nos = 0; nos < ns0; nos++)
      RELEASESHAREDS(objst[nos], structF0[nos]);
    for (E_Int nos = 0; nos < nu0; nos++)
      RELEASESHAREDU(objut[nos], unstrF0[nos], cnt0[nos]);
    for (E_Int nos = 0; nos < nwalls; nos++)
      RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
    return NULL;
  }

  // verification de posxv, posyv, poszv dans bodiesC  
  vector<E_Int> posxv; vector<E_Int> posyv; vector<E_Int> poszv;
  vector<E_Int> poscv;
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
    computeMininterf(ncellss, posx, posy, posz, 
                     structF0, posxv, posyv, poszv, poscv, unstrF, 
                     distances,cnt,isminortho);
  vector<FldArrayF*> distancesu;
  for (E_Int i = 0; i < nu0; i++)
  {
    E_Int ncells = ncellsu[i];
    FldArrayF* distance = new FldArrayF(ncells); 
    distance->setAllValuesAt(K_CONST::E_INFINITE);
    distancesu.push_back(distance);
  }
  if (unstrF0.size() > 0)
    computeMininterf(ncellsu, posx, posy, posz,
                     unstrF0, posxv, posyv, poszv, poscv, unstrF, 
                     distancesu, cnt, isminortho);

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
void K_DIST2WALLS::computeMininterf(
  vector<E_Int>& ncellst,
  E_Int posx, E_Int posy, E_Int posz, 
  vector<FldArrayF*>& fields, 
  vector<E_Int>& posxv, vector<E_Int>& posyv, vector<E_Int>& poszv, 
  vector<E_Int>& poscv, vector<FldArrayF*>& fieldsw, 
  vector<FldArrayF*>& distances,vector<FldArrayI*>& cntw,
  E_Int isminortho)
{

  /* 1 - creation du kdtree */
  E_Int nwalls = fieldsw.size(); // number of zones
  E_Int nptsmax = 0;
  E_Int npts_local = 0;
  vector<vector< vector<E_Int>  > >cVE_all;
  vector<E_Int> npts_walls_limit;
  for (E_Int v = 0; v < nwalls; v++)
  {
    nptsmax += fieldsw[v]->getSize();
    #include "mininterf_ortho_npts_limit.h"
  }
  FldArrayF* wallpts = new FldArrayF(nptsmax,3);
  E_Float* xw2 = wallpts->begin(1);
  E_Float* yw2 = wallpts->begin(2);
  E_Float* zw2 = wallpts->begin(3);
  E_Int c = 0;
  E_Int nzones = fields.size();
  
  // concatenate walls in a single array
  for (E_Int v = 0; v < nwalls; v++)
  {
    FldArrayF* fieldv = fieldsw[v];
    E_Int posxw = posxv[v]; E_Int posyw = posyv[v]; E_Int poszw = poszv[v];
    E_Float* xw = fieldv->begin(posxw);
    E_Float* yw = fieldv->begin(posyw);
    E_Float* zw = fieldv->begin(poszw);
    E_Int ncellsw = fieldv->getSize();
    /* recuperation des points calcules uniquement 
       pas de pts masques et interpoles dans kdtree */
    E_Int poscw = poscv[v]; E_Float* cellnw0 = fieldv->begin(poscw);
    
    for (E_Int i = 0; i < ncellsw; i++)
    {
      if (cellnw0[i] == 1.)
      { xw2[c] = xw[i]; yw2[c] = yw[i]; zw2[c] = zw[i]; c++; }
    }
  }

  if (c != wallpts->getSize()) wallpts->reAllocMat(c, 3); // si cellN
  
  if (c == 0) // no wall
  {
    for (E_Int v = 0; v < nzones; v++) 
    {
      E_Int ncells = ncellst[v];
      E_Float* distancep = distances[v]->begin();
      for (E_Int ind = 0; ind < ncells; ind++)
        distancep[ind] = sqrt(distancep[ind]);
    }
    return;
  }
  
  xw2 = wallpts->begin(1);
  yw2 = wallpts->begin(2);
  zw2 = wallpts->begin(3);
  
  ArrayAccessor<FldArrayF> coordAcc(*wallpts, 1,2,3);
  KdTree<FldArrayF> kdt(coordAcc, E_EPSILON);
  
  /* Detection de la paroi la plus proche */
  #pragma omp parallel
  {
    E_Float pt[3];
    E_Float p0[3]; E_Float p1[3]; E_Float p2[3]; E_Float p[3];
    E_Int ret, vw, indw2;
	  E_Float dist, distmin, dx, dy, dz;

    for (E_Int v = 0; v < nzones; v++)
    {
      E_Int ncells = ncellst[v];
      E_Float* xt = fields[v]->begin(posx);
      E_Float* yt = fields[v]->begin(posy);
      E_Float* zt = fields[v]->begin(posz);
      E_Float* distancep = distances[v]->begin();
      if (isminortho == 1) // mininterf_ortho
      {
        #pragma omp for schedule(dynamic)
        for (E_Int ind = 0; ind < ncells; ind++)
	      {
	        pt[0] = xt[ind]; pt[1] = yt[ind]; pt[2] = zt[ind];
	        indw2 = kdt.getClosest(pt);
          #include "mininterf.h"
          #include "mininterf_ortho.h"
	        if (ret != -1)
          {
	          dx = xp_local-pt[0]; dy = yp_local-pt[1]; dz = zp_local-pt[2];
	          dist = dx*dx + dy*dy + dz*dz;
	          if (dist < distmin) { distancep[ind] = sqrt(dist); distmin = dist; }
	        }
	      } // fin boucle
      }
      else // mininterf
      {
        #pragma omp for schedule(dynamic)
        for (E_Int ind = 0; ind < ncells; ind++)
	      {
	        pt[0] = xt[ind]; pt[1] = yt[ind]; pt[2] = zt[ind];
	        indw2 = kdt.getClosest(pt);
          #include "mininterf.h"
	      } // fin boucle
      }
    }
  }
  delete wallpts;
}
//=============================================================================
