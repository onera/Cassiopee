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

# include <stdio.h>
# include "transform.h"
using namespace std;
using namespace K_FLD;
using namespace K_SEARCH;

// ============================================================================
/* Computes a deformation vector for each point of a list of arrays */
// ============================================================================
PyObject* K_TRANSFORM::computeDeformationVector(PyObject* self,
                                                PyObject* args)
{
  PyObject* arrays; PyObject* deltas;
  E_Float beta;
  if (!PYPARSETUPLE_(args, OO_ R_,
                    &arrays, &deltas, &beta))
  {
      return NULL;
  }

  // Extract infos from mesh arrays for which delta must be computed
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt; vector<char*> eltType;
  vector<PyObject*> objst, objut;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = false;
  E_Bool skipUnstructured = false;
  E_Bool skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltType, objst, objut,
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nu = objut.size(); E_Int ns = objst.size();
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDeformationVector: invalid list of arrays.");
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    return NULL;
  }
  E_Int posx1, posy1, posz1;
  vector<E_Int> posxs; vector<E_Int> posys; vector<E_Int> poszs;
  vector<E_Int> posxu; vector<E_Int> posyu; vector<E_Int> poszu;
  for (E_Int nos = 0; nos < ns; nos++)
  {
    posx1 = K_ARRAY::isCoordinateXPresent(structVarString[nos]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(structVarString[nos]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(structVarString[nos]); posz1++;
    posxs.push_back(posx1); posys.push_back(posy1); poszs.push_back(posz1);
  }

  for (E_Int nou = 0; nou < nu; nou++)
  {
    posx1 = K_ARRAY::isCoordinateXPresent(unstrVarString[nou]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(unstrVarString[nou]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(unstrVarString[nou]); posz1++;
    posxu.push_back(posx1); posyu.push_back(posy1); poszu.push_back(posz1);
  }

  // Extract infos from surface arrays for which delta is known
  vector<E_Int> resls;
  vector<char*> structVarStrings; vector<char*> unstrVarStrings;
  vector<FldArrayF*> structFs; vector<FldArrayF*> unstrFs;
  vector<E_Int> nits; vector<E_Int> njts; vector<E_Int> nkts;
  vector<FldArrayI*> cnts;
  vector<char*> eltTypes;
  vector<PyObject*> objsts, objuts;
  skipNoCoord = true;
  skipStructured = true;
  skipUnstructured = false;
  skipDiffVars = true;
  isOk = K_ARRAY::getFromArrays(
    deltas, resls, structVarStrings, unstrVarStrings,
    structFs, unstrFs, nits, njts, nkts, cnts, eltTypes, objsts, objuts,
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);

  E_Int nwalls = objuts.size();
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDeformationVector: invalid list of arrays.");
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    for (E_Int nos = 0; nos < nwalls; nos++)
      RELEASESHAREDU(objuts[nos], unstrFs[nos], cnts[nos]);
    return NULL;
  }
  if (nwalls == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDeformationVector: no valid surface provided.");
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    for (E_Int nos = 0; nos < nwalls; nos++)
      RELEASESHAREDU(objuts[nos], unstrFs[nos], cnts[nos]);
    return NULL;
  }
  vector<E_Int> posxw; vector<E_Int> posyw; vector<E_Int> poszw;
  vector<E_Int> posdxw; vector<E_Int> posdyw; vector<E_Int> posdzw;
  E_Int posdx, posdy, posdz;
  E_Int found = 1;
  for (E_Int nou = 0; nou < nwalls; nou++)
  {
    posx1 = K_ARRAY::isCoordinateXPresent(unstrVarStrings[nou]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(unstrVarStrings[nou]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(unstrVarStrings[nou]); posz1++;
    posdx = K_ARRAY::isNamePresent("dx",unstrVarStrings[nou]); posdx++;
    if (posdx < 1) { found = 0; break; }
    posdy = K_ARRAY::isNamePresent("dy",unstrVarStrings[nou]); posdy++;
    if (posdx < 1) { found = 0; break; }
    posdz = K_ARRAY::isNamePresent("dz",unstrVarStrings[nou]); posdz++;
    if (posdx < 1) { found = 0; break; }
    posxw.push_back(posx1); posyw.push_back(posy1); poszw.push_back(posz1);
    posdxw.push_back(posdx); posdyw.push_back(posdy); posdzw.push_back(posdz);
  }
  if (found == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "computeDeformationVector: surface arrays must all contain dx, dy, dz variables.");
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    for (E_Int nos = 0; nos < nwalls; nos++)
      RELEASESHAREDU(objuts[nos], unstrFs[nos], cnts[nos]);
    return NULL;
  }

  // Creation du kdtree
  E_Int nptsmax = 0;
  for (E_Int v = 0; v < nwalls; v++) nptsmax += unstrFs[v]->getSize();
  FldArrayF* surfaces = new FldArrayF(nptsmax, 6);
  E_Float* xw =  surfaces->begin(1);
  E_Float* yw =  surfaces->begin(2);
  E_Float* zw =  surfaces->begin(3);
  E_Float* dxw = surfaces->begin(4);
  E_Float* dyw = surfaces->begin(5);
  E_Float* dzw = surfaces->begin(6);
  E_Int c = 0;
  for (E_Int v = 0; v < nwalls; v++)
  {
    FldArrayF* fieldv = unstrFs[v];
    E_Int posxv = posxw[v]; E_Int posyv = posyw[v]; E_Int poszv = poszw[v];
    E_Int posdxv = posdxw[v]; E_Int posdyv = posdyw[v]; E_Int posdzv = posdzw[v];
    E_Float* xw0 = fieldv->begin(posxv); E_Float* dxw0 = fieldv->begin(posdxv);
    E_Float* yw0 = fieldv->begin(posyv); E_Float* dyw0 = fieldv->begin(posdyv);
    E_Float* zw0 = fieldv->begin(poszv); E_Float* dzw0 = fieldv->begin(posdzv);
    E_Int nptsw = fieldv->getSize();
    for (E_Int i = 0; i < nptsw; i++)
    {
      xw[c] = xw0[i]; yw[c] = yw0[i]; zw[c] = zw0[i];
      dxw[c] = dxw0[i]; dyw[c] = dyw0[i]; dzw[c] = dzw0[i];
      c++;
    }
  } //fin kdtree

  ArrayAccessor<FldArrayF> coordAcc(*surfaces, 1, 2, 3);
  KdTree<FldArrayF> kdt(coordAcc);

  // Build arrays
  PyObject* l = PyList_New(0);
  vector<E_Float*> coordx; vector<E_Float*> dxs;
  vector<E_Float*> coordy; vector<E_Float*> dys;
  vector<E_Float*> coordz; vector<E_Float*> dzs;
  vector<E_Int> sizes;
  PyObject* tpl;
  FldArrayF* f2;
  E_Int nfld = 3;
  for (E_Int nos = 0; nos < ns; nos++)
  {
    E_Int api = structF[nos]->getApi();
    E_Int npts = structF[nos]->getSize();
    tpl = K_ARRAY::buildArray3(
      nfld, "dx,dy,dz", nit[nos], njt[nos], nkt[nos], api
    );
    K_ARRAY::getFromArray3(tpl, f2);
    f2->setAllValuesAtNull();
    coordx.push_back(structF[nos]->begin(posxs[nos])); dxs.push_back(f2->begin(1));
    coordy.push_back(structF[nos]->begin(posys[nos])); dys.push_back(f2->begin(2));
    coordz.push_back(structF[nos]->begin(poszs[nos])); dzs.push_back(f2->begin(3));
    sizes.push_back(npts);
    RELEASESHAREDS(tpl, f2);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }

  for (E_Int nou = 0; nou < nu; nou++)
  {
    E_Int api = unstrF[nou]->getApi();
    E_Int npts = unstrF[nou]->getSize();
    tpl = K_ARRAY::buildArray3(
      nfld, "dx,dy,dz", npts,
      *cnt[nou], eltType[nou], false, api, true
    );
    K_ARRAY::getFromArray3(tpl, f2);
    f2->setAllValuesAtNull();
    coordx.push_back(unstrF[nou]->begin(posxu[nou])); dxs.push_back(f2->begin(1));
    coordy.push_back(unstrF[nou]->begin(posyu[nou])); dys.push_back(f2->begin(2));
    coordz.push_back(unstrF[nou]->begin(poszu[nou])); dzs.push_back(f2->begin(3));
    sizes.push_back(npts);
    RELEASESHAREDS(tpl, f2);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }

  // deformation vector computation for each node
  // detection of the nearest surface point
  // On utilise une exponentielle en puissance de 2
  // On attend 70% d'amortissement a une distance de beta*deformation locale
  //E_Float beta = 7;
  //E_Int factor = 8;
  E_Float pt[3];
  beta = 1./(beta*beta);
  //beta = 1./pow(beta, factor);
  E_Float eps = 1.e-10;
  E_Float theta; E_Float beta0; E_Float delta;
  E_Float x0, y0, z0, dx0, dy0, dz0, dist;
  E_Int nzones = sizes.size();
  for (E_Int v = 0; v < nzones; v++)
  {
    E_Float* xt = coordx[v]; E_Float* yt = coordy[v]; E_Float* zt = coordz[v];
    E_Float* dxt = dxs[v]; E_Float* dyt = dys[v]; E_Float* dzt = dzs[v];
    E_Int npts = sizes[v];
    for (E_Int ind = 0; ind < npts; ind++)
    {
      pt[0] = xt[ind]; pt[1] = yt[ind]; pt[2] = zt[ind];
      E_Int indw = kdt.getClosest(pt);
      x0 = xw[indw]; y0 = yw[indw]; z0 = zw[indw];//nearest surface pt
      dx0 = dxw[indw]; dy0 = dyw[indw]; dz0 = dzw[indw];//delta vector of nearest surface pt
      delta = dx0*dx0+dy0*dy0+dz0*dz0;//norm of the delta vector
      dist = (xt[ind]-x0)*(xt[ind]-x0)+(yt[ind]-y0)*(yt[ind]-y0)+(zt[ind]-z0)*(zt[ind]-z0);
      beta0 = dist/(eps+delta);
      //beta0 = pow(beta0, factor/2);
      //beta0 = K_FUNC::E_min(100.,beta0);
      theta = exp(-beta*beta0);
      dxt[ind] = theta*dx0; dyt[ind] = theta*dy0; dzt[ind] = theta*dz0;
    }
  }

  delete surfaces;
  // Release memory
  for (E_Int nos = 0; nos < ns; nos++)
    RELEASESHAREDS(objst[nos], structF[nos]);
  for (E_Int nou = 0; nou < nu; nou++)
    RELEASESHAREDU(objut[nou], unstrF[nou], cnt[nou]);
  for (E_Int nou = 0; nou < nwalls; nou++)
    RELEASESHAREDU(objuts[nou], unstrFs[nou], cnts[nou]);
  return l;
}
