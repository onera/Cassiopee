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

# include "generator.h"
# include "Nuga/include/ArrayAccessor.h"
# include "Nuga/include/merge.h"
# include "Nuga/include/KdTree.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace K_CONST;
using namespace std;

// ============================================================================
/* Close borders of a set of zones if external vertices are distant from eps */
// ============================================================================
PyObject* K_GENERATOR::closeBorders(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  PyObject* arraysEF; //exterior faces for unstructured meshes
  E_Float eps;
  
  if (!PYPARSETUPLE_(args, O_ O_ R_ , &arrays, &arraysEF, &eps))   
    return NULL;

  // Extract infos from arrays
  vector<E_Int> res;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt; vector<char*> eltType;
  vector<PyObject*> objs, obju;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = false;
  E_Bool skipUnstructured = false;
  E_Bool skipDiffVars = false;

  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, res, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltType, objs, obju, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);

  E_Int size;
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "closeBorders: invalid list of arrays.");
    for (size_t v = 0; v < structF.size(); v++) RELEASESHAREDS(objs[v], structF[v]);
    for (size_t v = 0; v < unstrF.size(); v++) RELEASESHAREDU(obju[v], unstrF[v], cnt[v]);
    return NULL;
  }
  
  // Position des coordonnees
  vector<E_Int> posxs; vector<E_Int> posys; vector<E_Int> poszs;
  E_Int posx, posy, posz;
  size = structF.size();
  for (E_Int v = 0; v < size; v++)
  {
    posx = K_ARRAY::isCoordinateXPresent(structVarString[v]);
    posy = K_ARRAY::isCoordinateYPresent(structVarString[v]);
    posz = K_ARRAY::isCoordinateZPresent(structVarString[v]);
    posxs.push_back(posx+1);
    posys.push_back(posy+1);
    poszs.push_back(posz+1);
  }
  vector<E_Int> posxu; vector<E_Int> posyu; vector<E_Int> poszu;
  size = unstrF.size();
  for (E_Int v = 0; v < size; v++)
  {
    posx = K_ARRAY::isCoordinateXPresent(unstrVarString[v]);
    posy = K_ARRAY::isCoordinateYPresent(unstrVarString[v]);
    posz = K_ARRAY::isCoordinateZPresent(unstrVarString[v]);
    posxu.push_back(posx+1);
    posyu.push_back(posy+1);
    poszu.push_back(posz+1);
  }

  // output
  PyObject* out = PyList_New(0);
  PyObject* tpl;

  // Structure
  if (structF.size() > 0)
  {
    std::vector<FldArrayF*> cstructF;
    for (size_t v = 0; v < structF.size(); v++)
    {
      E_Int api = structF[v]->getApi();
      tpl = K_ARRAY::buildArray3(*structF[v], structVarString[v],
                                 nit[v], njt[v], nkt[v], api);
      RELEASESHAREDS(objs[v], structF[v]);
      FldArrayF* f2;
      K_ARRAY::getFromArray3(tpl, f2);
      cstructF.push_back(f2);
      PyList_Append(out, tpl);
      Py_DECREF(tpl);
    }
    closeAllStructuredMeshes(cstructF, nit, njt, nkt, posxs, posys, poszs, eps);
    for (size_t v = 0; v < cstructF.size(); v++) RELEASESHAREDS(PyList_GetItem(out,v), cstructF[v]);
  }

  // Non structure
  if (unstrF.size() > 0)
  { 
    // Extract infos from extFaces
    vector<E_Int> rese;
    vector<char*> structVarStringe; vector<char*> unstrVarStringe;
    vector<FldArrayF*> structEF; vector<FldArrayF*> unstrEF;
    vector<E_Int> nie; vector<E_Int> nje; vector<E_Int> nke;
    vector<FldArrayI*> cne; vector<char*> eltTypeE;
    vector<PyObject*> objse, objue;
    skipNoCoord = true;
    skipStructured = true;
    skipUnstructured = false;
    skipDiffVars = false;

    isOk = K_ARRAY::getFromArrays(
      arraysEF, rese, structVarStringe, unstrVarStringe,
      structEF, unstrEF, nie, nje, nke, cne, eltTypeE, objse, objue, 
      skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);

    if (isOk == -1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "closeBorders: invalid list of arrays for exterior faces.");

      for (size_t v = 0; v < structF.size(); v++) RELEASESHAREDS(objs[v], structF[v]);
      for (size_t v = 0; v < unstrF.size(); v++) RELEASESHAREDU(obju[v], unstrF[v], cnt[v]);
  
      for (size_t v = 0; v < structEF.size(); v++) RELEASESHAREDS(objse[v], structEF[v]);
      for (size_t v = 0; v < unstrEF.size(); v++) RELEASESHAREDU(objue[v], unstrEF[v], cne[v]);

      return NULL;
    }
    
    size = unstrEF.size();
    if (size > 0)
    {
      // Position des coordonnees
      vector<E_Int> posxe; vector<E_Int> posye; vector<E_Int> posze;

      for (E_Int v = 0; v < size; v++)
      {
        posx = K_ARRAY::isCoordinateXPresent(unstrVarStringe[v]);
        posy = K_ARRAY::isCoordinateYPresent(unstrVarStringe[v]);
        posz = K_ARRAY::isCoordinateZPresent(unstrVarStringe[v]);
        posxe.push_back(posx+1);
        posye.push_back(posy+1);
        posze.push_back(posz+1);
      }
      std::vector<FldArrayF*> cunstrF; std::vector<FldArrayI*> ccnt;
      for (size_t v = 0; v < unstrF.size(); v++)
      {
        E_Int api = unstrF[v]->getApi();
        tpl = K_ARRAY::buildArray3(*unstrF[v], unstrVarString[v],
                                   *cnt[v], eltType[v], api);
        RELEASESHAREDU(obju[v], unstrF[v], cnt[v]);
        FldArrayF* f2; FldArrayI* cn2;
        K_ARRAY::getFromArray3(tpl, f2, cn2);
        cunstrF.push_back(f2); ccnt.push_back(cn2);
        PyList_Append(out, tpl);
        Py_DECREF(tpl);
      }
      closeAllUnstructuredMeshes(cunstrF, ccnt, posxu, posyu, poszu, 
                                 unstrEF, cne, posxe, posye, posze, 
                                 eps);
      for (size_t v = 0; v < cunstrF.size(); v++) RELEASESHAREDU(PyList_GetItem(out,v), cunstrF[v], ccnt[v]);
    }
    for (size_t v = 0; v < unstrEF.size(); v++) RELEASESHAREDU(objue[v], unstrEF[v], cne[v]);
  }// unstructured

  /*
  PyObject* out = PyList_New(0);
  PyObject* tpl;
  for (size_t v = 0; v < structF.size(); v++)
  {
    FldArrayF& f0 = *structF[v];
    tpl = K_ARRAY::buildArray(f0, structVarString[v], nit[v], njt[v], nkt[v]);
    RELEASESHAREDS(objs[v], structF[v]);
    PyList_Append(out, tpl);
    Py_DECREF(tpl);
  }
  for (size_t v = 0; v < unstrF.size(); v++)
  {
    FldArrayF& f0 = *unstrF[v];
    tpl = K_ARRAY::buildArray(f0, unstrVarString[v], *cnt[v], -1, eltType[v],
                              false);
    RELEASESHAREDU(obju[v], unstrF[v], cnt[v]);
    PyList_Append(out, tpl);
    Py_DECREF(tpl);
  }
  */
  return out;
}

//=============================================================================
/* closeAllUnstructuredMeshes */
//=============================================================================
void K_GENERATOR::closeAllUnstructuredMeshes(
  vector<FldArrayF*>& unstructF, vector<FldArrayI*>& connectsEV,
  vector<E_Int>& posxt, vector<E_Int>& posyt, vector<E_Int>& poszt,
  vector<FldArrayF*>& exteriorFacesF, vector<FldArrayI*>& connectsEF,
  std::vector<E_Int>& posxe, std::vector<E_Int>& posye, std::vector<E_Int>& posze,
  E_Float eps)
{
  // Create a kdtree with all the external windows of all zones
  E_Int nzones = unstructF.size();

  FldArrayF dmin(nzones);
  FldArrayI indmin(nzones,2); // stockage indtmin, indmin par zone
  E_Int sizemax = 0; E_Int sizemaxloc = 0;
  for (E_Int v1 = 0; v1 < nzones; v1++)
  {    
    sizemaxloc = exteriorFacesF[v1]->getSize();
    sizemax += sizemaxloc;
  }

  // Creation du kdtree et des tableaux d'indirection
  FldArrayF ftemp(sizemax,3);
  E_Float* xp = ftemp.begin(1);
  E_Float* yp = ftemp.begin(2);
  E_Float* zp = ftemp.begin(3); 
  FldArrayI indirB(sizemax); FldArrayI indirI(sizemax);
  E_Int indt = 0; // index in kdtree
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    E_Float* xt = exteriorFacesF[noz]->begin(posxe[noz]);
    E_Float* yt = exteriorFacesF[noz]->begin(posye[noz]);
    E_Float* zt = exteriorFacesF[noz]->begin(posze[noz]);   
    E_Int nptsExt = exteriorFacesF[noz]->getSize();

    for (E_Int ind = 0; ind < nptsExt; ind++)
    {
      xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
      indirB[indt] = noz; indirI[indt] = ind;
      indt++;
    }
  }
  ArrayAccessor<FldArrayF> coordAcc(ftemp, 1, 2, 3);
  K_SEARCH::KdTree<FldArrayF> globalKdt(coordAcc);

  // indices of mesh points to be modified
  nzones = unstructF.size();
  vector<FldArrayI*> indicesOrig(nzones);
  E_Float tol_match = 1e-8;
  for (E_Int nov = 0; nov < nzones; nov++)
  {
    E_Int nPts = unstructF[nov]->getSize();
    FldArrayI* indicesOrigLoc = new FldArrayI(nPts);
    indicesOrigLoc->setAllValuesAt(-1);
    indicesOrig[nov] = indicesOrigLoc;
  }

  for (E_Int nov = 0; nov <nzones; nov++)
  {
    E_Int* indicesOrigLoc = indicesOrig[nov]->begin();
#pragma omp parallel default(shared)
  {    
    E_Float pt[3];
    E_Float xf,yf,zf,dx,dy,dz;
    E_Int ind;
    E_Int nPts = unstructF[nov]->getSize();
    E_Float* xt = unstructF[nov]->begin(posxt[nov]);
    E_Float* yt = unstructF[nov]->begin(posyt[nov]);
    E_Float* zt = unstructF[nov]->begin(poszt[nov]);   
#pragma omp for schedule(dynamic)
    for (E_Int i = 0; i < nPts; i++)
    {
      xf = xt[i]; yf = yt[i]; zf = zt[i];
      pt[0] = xf; pt[1] = yf; pt[2] = zf;
      ind = globalKdt.getClosest(pt); // closest pt      
      dx = xp[ind]-xf; dy = yp[ind]-yf; dz = zp[ind]-zf;
      if (K_FUNC::E_abs(dx) < tol_match && K_FUNC::E_abs(dy) < tol_match && K_FUNC::E_abs(dz) < tol_match) 
        indicesOrigLoc[i] = ind;
    }
    }
  }
  FldArrayI tag(sizemax); tag.setAllValuesAtNull();
  E_Float eps10 = 10*eps;
  E_Float bbmin[3], bbmax[3];
  vector<E_Int> candidates;
  E_Float dx, dy, dz;
  for (E_Int indt = 0; indt < sizemax; indt++)
  {
    E_Int noz = indirB[indt]; E_Int ind = indirI[indt];
    E_Float* xt = exteriorFacesF[noz]->begin(posxe[noz]);
    E_Float* yt = exteriorFacesF[noz]->begin(posye[noz]);
    E_Float* zt = exteriorFacesF[noz]->begin(posze[noz]);
    bbmin[0] = xt[ind]-eps10; bbmin[1] = yt[ind]-eps10; bbmin[2] = zt[ind]-eps10;
    bbmax[0] = xt[ind]+eps10; bbmax[1] = yt[ind]+eps10; bbmax[2] = zt[ind]+eps10;
    candidates.clear();
    globalKdt.getInBox(bbmin, bbmax, candidates);

    //dmin = K_CONST::E_MAX_FLOAT;
    //indmin = -1; indtmin = -1;
    dmin.setAllValuesAt(K_CONST::E_MAX_FLOAT);
    indmin.setAllValuesAt(-1);
    for (size_t cr = 0; cr < candidates.size(); cr++)
    {
      E_Int indkr = candidates[cr];
      E_Int nozr = indirB[indkr]; E_Int indr = indirI[indkr];
      if (nozr != noz)
      {
        E_Float* xr = exteriorFacesF[nozr]->begin(posxe[nozr]);
        E_Float* yr = exteriorFacesF[nozr]->begin(posye[nozr]);
        E_Float* zr = exteriorFacesF[nozr]->begin(posze[nozr]);
        dx = xr[indr]-xt[ind];
        dy = yr[indr]-yt[ind];
        dz = zr[indr]-zt[ind];
        if (K_FUNC::fEqualZero(dx,eps) == true &&
            K_FUNC::fEqualZero(dy,eps) == true &&
            K_FUNC::fEqualZero(dz,eps) == true && 
            tag[indt] == 0 && tag[indkr] == 0)
        {
          if (dx*dx+dy*dy+dz*dz < dmin[nozr])
          {
            dmin[nozr] = dx*dx+dy*dy+dz*dz;
            indmin(nozr,1) = indkr; indmin(nozr,2) = indr;
          }
        }
      }
    }
    E_Int accu = 1;
    E_Float xmean = xt[ind];
    E_Float ymean = yt[ind];
    E_Float zmean = zt[ind];
    for (E_Int nz = 0; nz < nzones; nz++)
    {
      E_Int im = indmin(nz,2);
      if (im != -1)
      {
        accu++;
        E_Float* xm = exteriorFacesF[nz]->begin(posxe[nz]);
        E_Float* ym = exteriorFacesF[nz]->begin(posye[nz]);
        E_Float* zm = exteriorFacesF[nz]->begin(posze[nz]);
        xmean += xm[im];
        ymean += ym[im];
        zmean += zm[im];
      }
    }
    if (accu > 1)
    {
      E_Float w = 1./accu;
      xmean = xmean*w; ymean = ymean*w; zmean = zmean*w;
      xt[ind] = xmean; yt[ind] = ymean; zt[ind] = zmean; tag[indt] = 1;
      for (E_Int nz = 0; nz < nzones; nz++)
      {
        E_Int im = indmin(nz,2);
        if (im != -1)
        {
          E_Float* xm = exteriorFacesF[nz]->begin(posxe[nz]);
          E_Float* ym = exteriorFacesF[nz]->begin(posye[nz]);
          E_Float* zm = exteriorFacesF[nz]->begin(posze[nz]);
          xm[im] = xmean; ym[im] = ymean; zm[im] = zmean;
          tag[indmin(nz,1)] = 1;
        }
      }
    }
  }  

  for (E_Int nov = 0; nov <nzones; nov++)
  {
    E_Int* indicesOrigLoc = indicesOrig[nov]->begin();
    E_Int nPts = unstructF[nov]->getSize();
    E_Float* xt = unstructF[nov]->begin(posxt[nov]);
    E_Float* yt = unstructF[nov]->begin(posyt[nov]);
    E_Float* zt = unstructF[nov]->begin(poszt[nov]);

    for (E_Int ind = 0; ind < nPts; ind++)
    {
      E_Int indt = indicesOrigLoc[ind];
      if (indt != -1)
      {
        E_Int noe = indirB[indt]; E_Int inde = indirI[indt];   
        E_Float* xm = exteriorFacesF[noe]->begin(posxe[noe]); 
        E_Float* ym = exteriorFacesF[noe]->begin(posye[noe]); 
        E_Float* zm = exteriorFacesF[noe]->begin(posze[noe]);
        xt[ind] = xm[inde]; 
        yt[ind] = ym[inde]; 
        zt[ind] = zm[inde]; 
      }
    }
  }
  
  for (E_Int nov = 0; nov < nzones; nov++) delete indicesOrig[nov];
}

//=============================================================================
/* closeAllStructuredMeshes */
//=============================================================================
void K_GENERATOR::closeAllStructuredMeshes(
  vector<FldArrayF*>& structF,
  vector<E_Int>& nit, vector<E_Int>& njt, vector<E_Int>& nkt,
  vector<E_Int>& posxt, vector<E_Int>& posyt, vector<E_Int>& poszt,
  E_Float eps)
{
  // Create a kdtree with all the external windows of all zones
  E_Int nzones = structF.size();

  FldArrayF dmin(nzones);
  FldArrayI indmin(nzones,2); // stockage indtmin, indmin par zone

  E_Int sizemax = 0;
  E_Int sizemaxloc = 0;
  for (E_Int v1 = 0; v1 < nzones; v1++)
  {    
    if (nkt[v1] > 1) sizemaxloc = 2*njt[v1]*nkt[v1]+2*nit[v1]*nkt[v1]+2*nit[v1]*njt[v1];
    else if (njt[v1] > 1) sizemaxloc = 2*njt[v1]+2*nit[v1];
    else sizemaxloc = 2;
    sizemax += sizemaxloc;
  }

  // Creation du kdtree et des tableaux d'indirection
  FldArrayF ftemp(sizemax,3);
  E_Float* xp = ftemp.begin(1);
  E_Float* yp = ftemp.begin(2);
  E_Float* zp = ftemp.begin(3);
  FldArrayI indirB(sizemax); FldArrayI indirI(sizemax);
  
  E_Int ind, shifti, shiftj, shiftk;
  E_Int indt = 0; // index in kdtree
  for (E_Int noz = 0; noz < nzones; noz++)
  {    
    E_Int ni = nit[noz]; E_Int nj = njt[noz]; E_Int nk = nkt[noz];
    E_Int ninj = ni*nj;
    E_Float* xt = structF[noz]->begin(posxt[noz]);
    E_Float* yt = structF[noz]->begin(posyt[noz]);
    E_Float* zt = structF[noz]->begin(poszt[noz]);

    // i = 1
    shifti = 0;
    for (E_Int k = 0; k < nk; k++)
      for (E_Int j = 0; j < nj; j++)
      {
        ind = shifti + j*ni + k*ninj;
        xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
        indirB[indt] = noz; indirI[indt] = ind;
        indt++;
      }
    // i = ni
    shifti = ni-1;
    for (E_Int k = 0; k < nk; k++)
      for (E_Int j = 0; j < nj; j++)
      {
        ind = shifti + j*ni + k*ninj;
        xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
        indirB[indt] = noz; indirI[indt] = ind;
        indt++;
      }
    if (nj > 1)
    {
      // j = 1
      shiftj = 0;
      for (E_Int k = 0; k < nk; k++)
        for (E_Int i = 0; i < ni; i++)
        {
          ind = i + shiftj*ni + k*ninj;
          xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
          indirB[indt] = noz; indirI[indt] = ind;
          indt++;
        }
      // j = nj
      shiftj = nj-1;
      for (E_Int k = 0; k < nk; k++)
        for (E_Int i = 0; i < ni; i++)
        {
          ind = i + shiftj*ni + k*ninj;
          xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
          indirB[indt] = noz; indirI[indt] = ind;
          indt++;
        }
    }
    if (nk > 1)
    {
      // k = 1
      shiftk = 0;
      for (E_Int j = 0; j < nj; j++)
        for (E_Int i = 0; i < ni; i++)
        {
          ind = i + j*ni + shiftk*ninj;
          xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
          indirB[indt] = noz; indirI[indt] = ind;
          indt++;
        }
      // k = nk
      shiftk = nk-1;
      for (E_Int j = 0; j < nj; j++)
        for (E_Int i = 0; i < ni; i++)
        {
          ind = i + j*ni + shiftk*ninj;
          xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
          indirB[indt] = noz; indirI[indt] = ind;
          indt++;
        }    
    }
  }
  ArrayAccessor<FldArrayF> coordAcc(ftemp, 1, 2, 3);
  K_SEARCH::KdTree<FldArrayF> globalKdt(coordAcc);
  FldArrayI tag(sizemax); tag.setAllValuesAtNull();
  E_Float eps10 = 10*eps;
  E_Float bbmin[3], bbmax[3];
  vector<E_Int> candidates;
  E_Float dx, dy, dz;
  for (E_Int indt = 0; indt < sizemax; indt++)
  {
    E_Int noz = indirB[indt]; E_Int ind = indirI[indt];
    
    E_Float* xt = structF[noz]->begin(posxt[noz]);
    E_Float* yt = structF[noz]->begin(posyt[noz]);
    E_Float* zt = structF[noz]->begin(poszt[noz]);
    bbmin[0] = xt[ind]-eps10; bbmin[1] = yt[ind]-eps10; bbmin[2] = zt[ind]-eps10;
    bbmax[0] = xt[ind]+eps10; bbmax[1] = yt[ind]+eps10; bbmax[2] = zt[ind]+eps10;
    candidates.clear();
    globalKdt.getInBox(bbmin, bbmax, candidates);
    //dmin = K_CONST::E_MAX_FLOAT;
    //indmin = -1; indtmin = -1;
    dmin.setAllValuesAt(K_CONST::E_MAX_FLOAT);
    indmin.setAllValuesAt(-1);
    for (size_t cr = 0; cr < candidates.size(); cr++)
    {
      E_Int indkr = candidates[cr];
      E_Int nozr = indirB[indkr]; E_Int indr = indirI[indkr];
      if (nozr != noz)
      {
        E_Float* xr = structF[nozr]->begin(posxt[nozr]);
        E_Float* yr = structF[nozr]->begin(posyt[nozr]);
        E_Float* zr = structF[nozr]->begin(poszt[nozr]);
        dx = xr[indr]-xt[ind];
        dy = yr[indr]-yt[ind];
        dz = zr[indr]-zt[ind];
        if (K_FUNC::fEqualZero(dx,eps) == true &&
            K_FUNC::fEqualZero(dy,eps) == true &&
            K_FUNC::fEqualZero(dz,eps) == true && 
            tag[indt] == 0 && tag[indkr] == 0)
        {
          if (dx*dx+dy*dy+dz*dz < dmin[nozr])
          {
            dmin[nozr] = dx*dx+dy*dy+dz*dz;
            indmin(nozr,1) = indkr; indmin(nozr,2) = indr;
          }
        }
      }
    }
    
    E_Int accu = 1;
    E_Float xmean = xt[ind];
    E_Float ymean = yt[ind];
    E_Float zmean = zt[ind];
    for (E_Int nz = 0; nz < nzones; nz++)
    {
      E_Int im = indmin(nz,2);
      if (im != -1)
      {
        accu++;
        E_Float* xm = structF[nz]->begin(posxt[nz]);
        E_Float* ym = structF[nz]->begin(posyt[nz]);
        E_Float* zm = structF[nz]->begin(poszt[nz]);
        xmean += xm[im];
        ymean += ym[im];
        zmean += zm[im];
      }
    }

    if (accu > 1)
    {
      E_Float w = 1./accu;
      xmean = xmean*w; ymean = ymean*w; zmean = zmean*w;
      xt[ind] = xmean; yt[ind] = ymean; zt[ind] = zmean; tag[indt] = 1;
      for (E_Int nz = 0; nz < nzones; nz++)
      {
        E_Int im = indmin(nz,2);
        if (im != -1)
        {
          E_Float* xm = structF[nz]->begin(posxt[nz]);
          E_Float* ym = structF[nz]->begin(posyt[nz]);
          E_Float* zm = structF[nz]->begin(poszt[nz]);
          xm[im] = xmean; ym[im] = ymean; zm[im] = zmean;
          tag[indmin(nz,1)] = 1;
        }
      }
    }
  } // tous les pts
}
