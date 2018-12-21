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
# include <string.h>
# include <stdio.h>
# include "generator.h"
#include "Def/DefTypes.h"

using namespace K_FLD;
using namespace std;

// ============================================================================
/* Effectue une projection d'un front defini dans arrays par un cellN
   sur la surface surface. */
// ============================================================================
PyObject* K_GENERATOR::snapFront(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  PyObject* surface;
  E_Int optimized;
#ifdef E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "OOl", &arrays, &surface, &optimized)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "OOi", &arrays, &surface, &optimized)) return NULL;
#endif
  // Extract infos from arrays
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt; vector<char*> eltType;
  vector<PyObject*> objst, objut;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false;
  E_Boolean skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltType, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nu = objut.size(); E_Int ns = objst.size();
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "snapFront: invalid list of arrays.");
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nou = 0; nou < nu; nou++)
      RELEASESHAREDU(objut[nou], unstrF[nou], cnt[nou]);
    return NULL;
  }

  E_Int posx1, posy1, posz1, posc1;
  vector<E_Int> posxs; vector<E_Int> posys; vector<E_Int> poszs; 
  vector<E_Int> poscs;
  vector<E_Int> posxu; vector<E_Int> posyu; vector<E_Int> poszu;
  vector<E_Int> poscu;
  for (E_Int nos = 0; nos < ns; nos++)
  {
    posx1 = K_ARRAY::isCoordinateXPresent(structVarString[nos]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(structVarString[nos]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(structVarString[nos]); posz1++;
    posc1 = K_ARRAY::isCellNatureField1Present(structVarString[nos]); posc1++;
    if (posx1 == 0 || posy1 == 0 || posz1 == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "snapFront: arrays must contain coordinates.");
      for (E_Int nos = 0; nos < ns; nos++)
        RELEASESHAREDS(objst[nos], structF[nos]);
      for (E_Int nou = 0; nou < nu; nou++)
        RELEASESHAREDU(objut[nou], unstrF[nou], cnt[nou]);
      return NULL;
    }
    if (posc1 == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "snapFront: arrays must contain cellN.");
      for (E_Int nos = 0; nos < ns; nos++)
        RELEASESHAREDS(objst[nos], structF[nos]);
      for (E_Int nou = 0; nou < nu; nou++)
        RELEASESHAREDU(objut[nou], unstrF[nou], cnt[nou]);
      return NULL;
    }
    posxs.push_back(posx1); posys.push_back(posy1); poszs.push_back(posz1);
    poscs.push_back(posc1);
  }
  for (E_Int nou = 0; nou < nu; nou++)
  {
    posx1 = K_ARRAY::isCoordinateXPresent(unstrVarString[nou]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(unstrVarString[nou]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(unstrVarString[nou]); posz1++;
    posc1 = K_ARRAY::isCellNatureField1Present(unstrVarString[nou]); posc1++;
    if (posx1 == 0 || posy1 == 0 || posz1 == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "snapFront: arrays must contain coordinates.");
      for (E_Int nos = 0; nos < ns; nos++)
        RELEASESHAREDS(objst[nos], structF[nos]);
      for (E_Int nou = 0; nou < nu; nou++)
        RELEASESHAREDU(objut[nou], unstrF[nou], cnt[nou]);
      return NULL;
    }
    if (posc1 == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "snapFront: arrays must contain cellN.");
      for (E_Int nos = 0; nos < ns; nos++)
        RELEASESHAREDS(objst[nos], structF[nos]);
      for (E_Int nou = 0; nou < nu; nou++)
        RELEASESHAREDU(objut[nou], unstrF[nou], cnt[nou]);
      return NULL;
    }
    posxu.push_back(posx1); posyu.push_back(posy1); poszu.push_back(posz1);
    poscu.push_back(posc1);
  }

  // Surface de projection
  E_Int im2, jm2, km2;
  FldArrayF* f2;
  FldArrayI* cn2;
  char* varString2;
  char* eltType2;
  E_Int res2 = K_ARRAY::getFromArray(surface, varString2, 
                                     f2, im2, jm2, km2, cn2, eltType2, true); 
  if (res2 != 2)
  {
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    RELEASESHAREDB(res2, surface, f2, cn2);
    PyErr_SetString(PyExc_TypeError,
                    "adapt2Front: surface must be unstructured.");
    return NULL;
  }
  if (strcmp(eltType2, "TRI") != 0 && strcmp(eltType2, "BAR") != 0)
  {
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    RELEASESHAREDU(surface, f2, cn2);
    PyErr_SetString(PyExc_TypeError,
                    "adapt2Front: surface must be a TRI or BAR array.");
    return NULL;
  }

  E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
   
  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    RELEASESHAREDU(surface, f2, cn2);
    PyErr_SetString(PyExc_TypeError,
                    "adapt2Front: coordinates not found in surface.");
    return NULL;
  }
  posx2++; posy2++; posz2++;
  
  // Build arrays
  PyObject* l = PyList_New(0);
  vector<E_Float*> fs;
  vector<E_Float*> fu;
  for (E_Int nos = 0; nos < ns; nos++)
  {
    E_Int nfld = structF[nos]->getNfld(); E_Int npts = structF[nos]->getSize();
    PyObject* tpl = K_ARRAY::buildArray(nfld, structVarString[nos], nit[nos], njt[nos], nkt[nos]);
    E_Float* fp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f(npts, nfld, fp, true); f = *structF[nos];
    fs.push_back(fp);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }

  for (E_Int nou = 0; nou < nu; nou++)
  {
    E_Int nfld = unstrF[nou]->getNfld(); E_Int npts = unstrF[nou]->getSize();
    E_Int nelts = cnt[nou]->getSize(); E_Int nvert = cnt[nou]->getNfld();
    PyObject* tpl = K_ARRAY::buildArray(nfld, unstrVarString[nou], 
                                        npts, nelts, -1, eltType[nou], false, nelts);
    E_Float* fp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f(npts, nfld, fp, true); f = *unstrF[nou];
    E_Int* cnpo = K_ARRAY::getConnectPtr(tpl);
    FldArrayI cno(nelts, nvert, cnpo, true); cno = *cnt[nou];
    fu.push_back(fp);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }

  // liste des points du front
  E_Int ind, indv;
  E_Float x1,y1,z1,xo1,yo1,zo1,d1,h;
  vector<E_Float*> coordx;
  vector<E_Float*> coordy;
  vector<E_Float*> coordz;
  vector<E_Int> sizet;
  vector<E_Int*> indices;

  // optimized=1 : Front optimise avec node-in
  // optimized=2 : Front optimise avec cell intersect (le front 0 passe a 1)

  if (optimized == 1 || optimized == 2)
  {
    for (E_Int nos = 0; nos < ns; nos++)
    {
      E_Int npts = structF[nos]->getSize();
      FldArrayF& f = *structF[nos];
      E_Float* xp = f.begin(posxs[nos]);
      E_Float* yp = f.begin(posys[nos]);
      E_Float* zp = f.begin(poszs[nos]);
      E_Float* cellN = f.begin(poscs[nos]);
      E_Int ni = nit[nos]; E_Int nj = njt[nos]; E_Int nk = nkt[nos];
      E_Int* indi = new E_Int [npts]; // pts du front

      E_Int np = 0;
      E_Int corners = 0;
      if (optimized == 2) corners = 1;
      computeStructFrontIndi(cellN, ni, nj, nk, 0, 1, corners, indi, np);

      if (optimized == 2)
      {
        for (E_Int i = 0; i < np; i++) cellN[indi[i]] = 1;
      }
     
      // realloc
      E_Int* indin = new E_Int [np];
      for (E_Int i = 0; i < np; i++) indin[i] = indi[i];
      delete [] indi; 
      indices.push_back(indin);
      sizet.push_back(np);
        
      // calcul des coords
      E_Float* xt = new E_Float [np];
      E_Float* yt = new E_Float [np];
      E_Float* zt = new E_Float [np];
      for (E_Int i = 0; i < np; i++)
      { 
        xt[i] = xp[indin[i]];
        yt[i] = yp[indin[i]];
        zt[i] = zp[indin[i]];
      }
      coordx.push_back(xt);
      coordy.push_back(yt);
      coordz.push_back(zt);
    }

    for (E_Int nou = 0; nou < nu; nou++)
    {
      E_Int npts = unstrF[nou]->getSize();
      FldArrayF& f = *unstrF[nou];
      E_Float* xp = f.begin(posxu[nou]);
      E_Float* yp = f.begin(posyu[nou]);
      E_Float* zp = f.begin(poszu[nou]);
      E_Float* cellN = f.begin(poscu[nou]);
      E_Int* indi = new E_Int [npts]; // pts du front

      E_Int np = 0;
      E_Int corners = 0;
      if (optimized == 2) corners = 1;
      computeUnstrFrontIndi(cellN, npts, 0, 1, corners, *cnt[nou], indi, np);

      if (optimized == 2)
      {
        for (E_Int i = 0; i < np; i++) cellN[indi[i]] = 1;
      }     
      // realloc
      E_Int* indin = new E_Int [np];
      for (E_Int i = 0; i < np; i++) indin[i] = indi[i];
      delete [] indi; 
      indices.push_back(indin);
      sizet.push_back(np);
        
      // calcul des coords
      E_Float* xt = new E_Float [np];
      E_Float* yt = new E_Float [np];
      E_Float* zt = new E_Float [np];
      for (E_Int i = 0; i < np; i++)
      { 
        xt[i] = xp[indin[i]];
        yt[i] = yp[indin[i]];
        zt[i] = zp[indin[i]];
      }
      coordx.push_back(xt);
      coordy.push_back(yt);
      coordz.push_back(zt);
    }

    if (optimized == 2) goto nextf;
    // Calcul la projection sur les surfaces
    K_COMPGEOM::projectOrthoWithPrecond(posx2, posy2, posz2, *cn2, *f2, sizet, 
                                        coordx, coordy, coordz);

    // Met le cellN a 1 pour les pts proches de la paroi
    for (E_Int nos = 0; nos < ns; nos++)
    {
      FldArrayF& f = *structF[nos];
      E_Int npts = f.getSize();
      E_Float* cellN = f.begin(poscs[nos]);
      E_Float* xp = f.begin(posxs[nos]);
      E_Float* yp = f.begin(posys[nos]);
      E_Float* zp = f.begin(poszs[nos]);
      E_Float* xpt = fs[nos] + (posxs[nos]-1)*npts;
      E_Float* ypt = fs[nos] + (posys[nos]-1)*npts;
      E_Float* zpt = fs[nos] + (poszs[nos]-1)*npts;
      E_Float* cellNpt = fs[nos] + (poscs[nos]-1)*npts;
      E_Int* indi = indices[nos];
      E_Float* xt = coordx[nos];
      E_Float* yt = coordy[nos];
      E_Float* zt = coordz[nos];
      E_Int np = sizet[nos];
      E_Int ni = nit[nos]; E_Int nj = njt[nos];
      h = K_FUNC::E_abs(xp[ni/2+nj/2*ni+1]-xp[ni/2+nj/2*ni]); // pas h valide seulement pour un bloc cartesien regulier
     
      for (E_Int i = 0; i < np; i++)
      { 
        ind = indi[i];
        xo1 = xp[ind]; yo1 = yp[ind]; zo1 = zp[ind]; // point du front non projete
        x1 = xt[i]; y1 = yt[i]; z1 = zt[i]; // point du front projete
        d1 = (x1-xo1)*(x1-xo1)+(y1-yo1)*(y1-yo1)+(z1-zo1)*(z1-zo1);
        if (d1 < 0.25*h*h) 
        { 
          cellN[ind] = 1;
          /* force la projection: peut creer des cellules plates */
          cellNpt[ind] = 1;
          xpt[ind] = xt[i]; ypt[ind] = yt[i]; zpt[ind] = zt[i];
        }
      }
    }
  
    // Met le cellN a 1 pour les pts proches
    for (E_Int nou = 0; nou < nu; nou++)
    {
      FldArrayF& f = *unstrF[nou];
      E_Int npts = f.getSize();
      E_Float* cellN = f.begin(poscu[nou]);
      E_Float* xp = f.begin(posxu[nou]);
      E_Float* yp = f.begin(posyu[nou]);
      E_Float* zp = f.begin(poszu[nou]);
      E_Float* xpt = fu[nou] + (posxu[nou]-1)*npts;
      E_Float* ypt = fu[nou] + (posyu[nou]-1)*npts;
      E_Float* zpt = fu[nou] + (poszu[nou]-1)*npts;
      E_Float* cellNpt = fu[nou] + (poscu[nou]-1)*npts;
      E_Int* indi = indices[ns+nou];
      E_Float* xt = coordx[ns+nou];
      E_Float* yt = coordy[ns+nou];
      E_Float* zt = coordz[ns+nou];
      E_Int np = sizet[ns+nou];
      FldArrayI& cn = *cnt[nou];
      
      E_Float x0 = xp[cn(0,1)-1];
      h = 0.;
      for (E_Int nv = 2; nv <= cn.getNfld(); nv++)
        h = K_FUNC::E_max(h, K_FUNC::E_abs(xp[cn(0,nv)-1]-x0));
     
      for (E_Int i = 0; i < np; i++)
      { 
        ind = indi[i];
        xo1 = xp[ind]; yo1 = yp[ind]; zo1 = zp[ind];
        x1 = xt[i]; y1 = yt[i]; z1 = zt[i];
        d1 = (x1-xo1)*(x1-xo1)+(y1-yo1)*(y1-yo1)+(z1-zo1)*(z1-zo1);
        if (d1 < 0.25*h*h) 
        { 
          cellN[ind] = 1;
          /* force la projection: peut creer des cellules plates */
          cellNpt[ind] = 1;
          xpt[ind] = xt[i]; ypt[ind] = yt[i]; zpt[ind] = zt[i];
        }
      }
    }

    // delete
    for (E_Int nos = 0; nos < ns; nos++)
    {
      delete [] indices[nos]; // for output
      delete [] coordx[nos];
      delete [] coordy[nos];
      delete [] coordz[nos];
    }
    for (E_Int nou = 0; nou < nu; nou++)
    {
      delete [] indices[ns+nou]; // for output
      delete [] coordx[ns+nou];
      delete [] coordy[ns+nou];
      delete [] coordz[ns+nou];
    } 
    indices.clear();
    coordx.clear(); coordy.clear(); coordz.clear(); sizet.clear();

  } // Fin du front optimise --

  // Calcul le front cellN=1
  for (E_Int nos = 0; nos < ns; nos++)
  {
    E_Int npts = structF[nos]->getSize();
    FldArrayF& f = *structF[nos];
    E_Float* xp = f.begin(posxs[nos]);
    E_Float* yp = f.begin(posys[nos]);
    E_Float* zp = f.begin(poszs[nos]); 
    E_Float* cellN = f.begin(poscs[nos]);
    E_Int ni = nit[nos]; E_Int nj = njt[nos]; E_Int nk = nkt[nos];
    E_Int* indi = new E_Int [npts]; // pts du front

      E_Int np = 0;
      computeStructFrontIndi(cellN, ni, nj, nk, 1, 0, 1, indi, np);
    
    // realloc
    E_Int* indin = new E_Int [np];
    for (E_Int i = 0; i < np; i++) indin[i] = indi[i];
    delete [] indi;
    indices.push_back(indin);
    sizet.push_back(np);

    // calcul des coords
    E_Float* xt = new E_Float [np];
    E_Float* yt = new E_Float [np];
    E_Float* zt = new E_Float [np];
    for (E_Int i = 0; i < np; i++)
    { 
      xt[i] = xp[indin[i]];
      yt[i] = yp[indin[i]];
      zt[i] = zp[indin[i]];
    }
    coordx.push_back(xt);
    coordy.push_back(yt);
    coordz.push_back(zt);
  }

  for (E_Int nou = 0; nou < nu; nou++)
  {
    E_Int npts = unstrF[nou]->getSize();
    FldArrayF& f = *unstrF[nou];
    E_Float* xp = f.begin(posxu[nou]);
    E_Float* yp = f.begin(posyu[nou]);
    E_Float* zp = f.begin(poszu[nou]); 
    E_Float* cellN = f.begin(poscu[nou]);
    FldArrayI& cn = *cnt[nou];

    // tagg indi
    E_Int ne = cn.getSize();
    E_Int nfld = cn.getNfld();
    E_Int zero;
    E_Int* tag = new E_Int [npts]; // pts du front
    for (E_Int i = 0; i < npts; i++) tag[i] = 0;

    for (E_Int ie = 0; ie < ne; ie++) // elts
    {
      zero = 0;
      for (E_Int nv = 1; nv <= nfld; nv++)
      {
        indv = cn(ie, nv)-1;
        if (cellN[indv] == 0) zero++;
      }
      if (zero > 0 && zero < nfld)
      {
        for (E_Int nv = 1; nv <= nfld; nv++)
        {
          indv = cn(ie, nv)-1;
          if (cellN[indv] == 1) tag[indv] = 1;
        }
      }
    }

    // gathers indi
    E_Int np = 0;
    E_Int* indi = new E_Int [npts]; // pts du front
    for (E_Int ind = 0; ind < npts; ind++)
    {
      if (tag[ind] == 1) { indi[np] = ind; np++; }
    }
    delete [] tag;

    // realloc
    E_Int* indin = new E_Int [np];
    for (E_Int i = 0; i < np; i++) indin[i] = indi[i];
    delete [] indi;
    indices.push_back(indin);
    sizet.push_back(np);

    // calcul des coords
    E_Float* xt = new E_Float [np];
    E_Float* yt = new E_Float [np];
    E_Float* zt = new E_Float [np];
    for (E_Int i = 0; i < np; i++)
    { 
      xt[i] = xp[indin[i]];
      yt[i] = yp[indin[i]];
      zt[i] = zp[indin[i]];
    }
    coordx.push_back(xt);
    coordy.push_back(yt);
    coordz.push_back(zt);
  }

  nextf:;
  // Projete les points du front
  K_COMPGEOM::projectOrthoWithPrecond(posx2, posy2, posz2, *cn2, *f2, sizet, 
                                      coordx, coordy, coordz);
 // 
  // Patch back pour le structure
  for (E_Int nos = 0; nos < ns; nos++)
  {
    FldArrayF& f = *structF[nos];
    E_Int npts = f.getSize();
    E_Float* xp = fs[nos] + (posxs[nos]-1)*npts;
    E_Float* yp = fs[nos] + (posys[nos]-1)*npts;
    E_Float* zp = fs[nos] + (poszs[nos]-1)*npts;
    E_Float* cellN = fs[nos] + (poscs[nos]-1)*npts;
    E_Int* indi = indices[nos];
    E_Float* xt = coordx[nos];
    E_Float* yt = coordy[nos];
    E_Float* zt = coordz[nos];
    E_Int np = sizet[nos];
    for (E_Int i = 0; i < np; i++)
    { 
      ind = indi[i];
      xp[ind] = xt[i];
      yp[ind] = yt[i];
      zp[ind] = zt[i];
      cellN[ind] = 1;
    }
  }

  // Patch back pour le non structure
  for (E_Int nou = 0; nou < nu; nou++)
  {
    FldArrayF& f = *unstrF[nou];
    E_Int npts = f.getSize();
    E_Float* xp = fu[nou] + (posxu[nou]-1)*npts;
    E_Float* yp = fu[nou] + (posyu[nou]-1)*npts;
    E_Float* zp = fu[nou] + (poszu[nou]-1)*npts;
    E_Float* cellN = fu[nou] + (poscu[nou]-1)*npts;
    E_Int* indi = indices[ns+nou];
    E_Float* xt = coordx[ns+nou];
    E_Float* yt = coordy[ns+nou];
    E_Float* zt = coordz[ns+nou];
    E_Int np = sizet[ns+nou];
    for (E_Int i = 0; i < np; i++)
    { 
      ind = indi[i];
      xp[ind] = xt[i];
      yp[ind] = yt[i];
      zp[ind] = zt[i];
      cellN[ind] = 1;
    }
  }
  
  // delete
  for (E_Int nos = 0; nos < ns; nos++)
  {
    delete [] indices[nos]; // for output
    delete [] coordx[nos];
    delete [] coordy[nos];
    delete [] coordz[nos];
  }
  for (E_Int nou = 0; nou < nu; nou++)
  {
    delete [] indices[ns+nou]; // for output
    delete [] coordx[ns+nou];
    delete [] coordy[ns+nou];
    delete [] coordz[ns+nou];
  }

  RELEASESHAREDU(surface, f2, cn2);
  for (E_Int nos = 0; nos < ns; nos++)
    RELEASESHAREDS(objst[nos], structF[nos]);
  for (E_Int nou = 0; nou < nu; nou++)
    RELEASESHAREDU(objut[nou], unstrF[nou], cnt[nou]);
  return l;
}
// ============================================================================
/* Construit un tableau d indirection pour determiner un front par un cellN pour un 
   maillage structure:
   Le front construit correspond aux cellules cellN=var1 ayant au moins 
   un voisin avec cellN=var2.
   IN: cellN: champ cellNatureField du maillage
   IN: var1, var2: valeurs de cellN determinant le front
   IN: meshType: 0: maillage structure. 1: maillage non structure
   IN: corners: prise en compte des voisins par les coins dans le front (valeur 0 ou 1)
   OUT: indi: tableau d indirection entre le maillage et le front
   OUT: np : nombre de points du front
 */
// ============================================================================
void K_GENERATOR::computeStructFrontIndi(E_Float* cellN, E_Int ni, E_Int nj, E_Int nk, E_Int var1, E_Int var2,
                                         E_Int corners,
                                         E_Int* indi, E_Int& np)
{
  // indices de maillage pour calculer le front
  E_Int ind, indv, ivm, ivp, jvm, jvp, kvm, kvp;
  E_Int nij = ni*nj;


#define BRANCH indi[np] = ind; np++; goto next2;
  // Calcul le front cellN=var1
  for (E_Int k = 0; k < nk; k++)
    for (E_Int j = 0; j < nj; j++)
      for (E_Int i = 0; i < ni; i++)
      {
        ind = i + j*ni + k*nij;
          
        if (cellN[ind] == var1)
        {
          ivm = max(i-1, 0);
          indv = ivm + j*ni + k*nij;
          if (cellN[indv] == var2) { BRANCH; }
          ivp = min(i+1, ni-1);
          indv = ivp + j*ni + k*nij;
          if (cellN[indv] == var2) { BRANCH; }
          jvm = max(j-1, 0);
          indv = i + jvm*ni + k*nij;
          if (cellN[indv] == var2) { BRANCH; }
          jvp = min(j+1, nj-1);
          indv = i + jvp*ni + k*nij;
          if (cellN[indv] == var2) { BRANCH; }
          kvm = max(k-1, 0);
          indv = i + j*ni + kvm*nij;
          if (cellN[indv] == var2) { BRANCH; }
          kvp = min(k+1, nk-1);
          indv = i + j*ni + kvp*nij;
          if (cellN[indv] == var2) { BRANCH; }
              
          /* Les pts de coins peuvent nuire a la connexite de l'enveloppe */
          if (corners == 1)
          {
            indv = ivm + jvm*ni + k*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = ivm + jvp*ni + k*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = ivp + jvm*ni + k*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = ivp + jvp*ni + k*nij;
            if (cellN[indv] == var2) { BRANCH; }
                
            indv = ivm + jvm*ni + kvm*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = ivm + jvp*ni + kvm*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = ivp + jvm*ni + kvm*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = ivp + jvp*ni + kvm*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = ivm + j*ni + kvm*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = ivp + j*ni + kvm*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = i + jvm*ni + kvm*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = i + jvp*ni + kvm*nij;
            if (cellN[indv] == var2) { BRANCH; }
                
            indv = ivm + jvm*ni + kvp*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = ivm + jvp*ni + kvp*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = ivp + jvm*ni + kvp*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = ivp + jvp*ni + kvp*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = ivm + j*ni + kvp*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = ivp + j*ni + kvp*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = i + jvm*ni + kvp*nij;
            if (cellN[indv] == var2) { BRANCH; }
            indv = i + jvp*ni + kvp*nij;
            if (cellN[indv] == var2) { BRANCH; }
          }
        }
        next2:;
      }
#undef BRANCH
  return;
}
// ============================================================================
/* Construit un tableau d indirection pour determiner un front par un cellN pour un 
   maillage non structure:
   Le front construit correspond aux cellules cellN=var1 ayant au moins 
   un voisin avec cellN=var2.
   IN: cellN: champ cellNatureField du maillage
   IN: npts: taille du maillage
   IN: var1, var2: valeurs de cellN determinant le front
   IN: corners: prise en compte des voisins par les coins dans le front (valeur 0 ou 1)
   IN: cn: connectivite
   OUT: indi: tableau d indirection entre le maillage et le front
   OUT: np : nombre de points du front
 */
// ============================================================================
void K_GENERATOR::computeUnstrFrontIndi(E_Float* cellN, E_Int npts, E_Int var1, E_Int var2,
                                        E_Int corners, FldArrayI& cn,
                                        E_Int* indi, E_Int& np)
{
  // indices du voisin
  E_Int indv;

#define BRANCH indi[np] = ind; np++; goto next2u;
  vector< vector<E_Int> > cVN(npts); 
  K_CONNECT::connectEV2VNbrs(cn, cVN, corners);
  // Calcul le front cellN=var1
  for (E_Int ind = 0; ind < npts; ind++)
  {          
    if (cellN[ind] == var1)
    {
      vector<E_Int>& voisin = cVN[ind];
      for (unsigned int nv = 0; nv < voisin.size(); nv++)
      {
        indv = voisin[nv]-1;
        if (cellN[indv] == var2) { BRANCH; }
      }
    }
    next2u:;
  }
#undef BRANCH
  return;
}
