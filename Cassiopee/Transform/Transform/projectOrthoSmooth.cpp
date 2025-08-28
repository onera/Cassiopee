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
# include <string.h>
# include <stdio.h>
# include "transform.h"

using namespace K_FLD;
using namespace std;

// ============================================================================
/* Projete un surface array1 sur un surface array2 (TRI).
   La normale est lissee si le projete n'est pas regulier. */
// ============================================================================
PyObject* K_TRANSFORM::projectOrthoSmooth(PyObject* self, PyObject* args)
{
  PyObject* arrays; PyObject* array2;
  if (!PyArg_ParseTuple(args, "OO", &arrays, &array2))
  {
    return NULL;
  }

  // Extract infos from arrays
  vector<E_Int> resl;
  vector<char*> structVarString;
  vector<char*> unstrVarString;
  vector<FldArrayF*> structF;
  vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltType;
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
                    "projectOrthoSmooth: invalid list of arrays.");
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
  // Projection surface array
  E_Int im2, jm2, km2;
  FldArrayF* f2; FldArrayI* cn2;
  char* varString2; char* eltType2;
  E_Int res2 = K_ARRAY::getFromArray3(array2, varString2, 
                                      f2, im2, jm2, km2, cn2, eltType2); 
  if (res2 != 2)
  {
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    RELEASESHAREDB(res2, array2, f2, cn2);
    PyErr_SetString(PyExc_TypeError,
                    "projectOrthoSmooth: array2 must be unstructured.");
    return NULL;
  }
  if (strcmp(eltType2, "TRI") != 0 && strcmp(eltType2, "BAR") != 0)
  {
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    RELEASESHAREDU(array2, f2, cn2);
    PyErr_SetString(PyExc_TypeError,
                    "projectOrthoSmooth: array2 must be a TRI or BAR array.");
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
    RELEASESHAREDU(array2, f2, cn2);
    PyErr_SetString(PyExc_TypeError,
                    "projectOrthoSmooth: coordinates not found in array2.");
    return NULL;
  }
  posx2++; posy2++; posz2++;
  
  // Build arrays
  PyObject* l = PyList_New(0);
  vector<E_Float*> coordx; vector<E_Float*> coordy; vector<E_Float*> coordz;
  vector<E_Int> sizet;
  for (E_Int nos = 0; nos < ns; nos++)
  {
    E_Int nfld = structF[nos]->getNfld(); E_Int npts = structF[nos]->getSize();
    PyObject* tpl = K_ARRAY::buildArray(nfld, structVarString[nos], 
                                        nit[nos], njt[nos], nkt[nos]);
    E_Float* fp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f(npts, nfld, fp, true); f = *structF[nos];
    coordx.push_back(f.begin(posxs[nos]));
    coordy.push_back(f.begin(posys[nos]));
    coordz.push_back(f.begin(poszs[nos]));
    sizet.push_back(f.getSize());
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }

  for (E_Int nou = 0; nou < nu; nou++)
  {
    E_Int nfld = unstrF[nou]->getNfld(); E_Int npts = unstrF[nou]->getSize();
    E_Int nelts = cnt[nou]->getSize(); E_Int nvert = cnt[nou]->getNfld();
    PyObject* tpl = K_ARRAY::buildArray(nfld, unstrVarString[nou], npts, 
                                        nelts, -1, eltType[nou], false, nelts);
    E_Float* fp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f(npts, nfld, fp, true); f = *unstrF[nou];
    E_Int* cnpo = K_ARRAY::getConnectPtr(tpl);
    FldArrayI cno(nelts, nvert, cnpo, true); cno = *cnt[nou];    
    coordx.push_back(f.begin(posxu[nou]));
    coordy.push_back(f.begin(posyu[nou]));
    coordz.push_back(f.begin(poszu[nou]));
    sizet.push_back(f.getSize());
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  // Projete
  K_COMPGEOM::projectOrthoWithPrecond(posx2, posy2, posz2, *cn2, *f2, 
                                      sizet, coordx, coordy, coordz);

  // On cree les directions de projection et la qualite
  E_Int no = 0;
  vector<E_Float*> nxt;
  vector<E_Float*> nyt;
  vector<E_Float*> nzt;
  vector<E_Float*> qual;
  for (E_Int nos = 0; nos < ns; nos++)
  {
    E_Int npts = structF[nos]->getSize();
    E_Float* nx = new E_Float[npts]; nxt.push_back(nx);
    E_Float* ny = new E_Float[npts]; nyt.push_back(ny);
    E_Float* nz = new E_Float[npts]; nzt.push_back(nz);
    E_Float* q = new E_Float[npts]; qual.push_back(q);
    no++;
  }
  for (E_Int nou = 0; nou < nu; nou++)
  {
    E_Int npts = unstrF[nou]->getSize();
    E_Float* nx = new E_Float[npts]; nxt.push_back(nx);
    E_Float* ny = new E_Float[npts]; nyt.push_back(ny);
    E_Float* nz = new E_Float[npts]; nzt.push_back(nz);
    E_Float* q = new E_Float[npts]; qual.push_back(q);
    no++;
  }

  // Calcul des directions et de la qualite
  no = 0;
  for (E_Int nos = 0; nos < ns; nos++)
  {
    FldArrayF& orig = *structF[nos];
    E_Float* xo = orig.begin(posxs[nos]); // original surface
    E_Float* yo = orig.begin(posys[nos]);
    E_Float* zo = orig.begin(poszs[nos]);
    E_Int npts = structF[nos]->getSize();
    E_Float* x = coordx[no]; // projete 
    E_Float* y = coordy[no];
    E_Float* z = coordz[no];
    E_Float* nx = nxt[no];
    E_Float* ny = nyt[no];
    E_Float* nz = nzt[no];
    E_Float* q = qual[no];
    
    for (E_Int i = 0; i < npts; i++)
    {
      nx[i] = x[i]-xo[i];
      ny[i] = y[i]-yo[i];
      nz[i] = z[i]-zo[i];
    }
    E_Int im = nit[nos]; E_Int jm = njt[nos]; E_Int km = nkt[nos];
    K_METRIC::compStructSurft(im, jm, km, x, y, z, q);
    no++;
  }
  for (E_Int nou = 0; nou < nu; nou++)
  {
    FldArrayF& orig = *unstrF[nou];
    E_Float* xo = orig.begin(posxu[nou]); // original surface
    E_Float* yo = orig.begin(posyu[nou]);
    E_Float* zo = orig.begin(poszu[nou]);
    E_Int npts = unstrF[nou]->getSize();
    E_Float* x = coordx[nou]; // projete 
    E_Float* y = coordy[nou];
    E_Float* z = coordz[nou];
    E_Float* nx = nxt[nou];
    E_Float* ny = nyt[nou];
    E_Float* nz = nzt[nou];
    
    for (E_Int i = 0; i < npts; i++)
    {
      nx[i] = x[i]-xo[i];
      ny[i] = y[i]-yo[i];
      nz[i] = z[i]-zo[i];
    }
    
    no++;
  }

  // delete direction et qualite
   for (E_Int no = 0; no < ns+nu; no++)
  {
    delete [] nxt[no];
    delete [] nyt[no];
    delete [] nzt[no];
    delete [] qual[no];
  }
                                   
  RELEASESHAREDU(array2, f2, cn2);
  for (E_Int nos = 0; nos < ns; nos++)
    RELEASESHAREDS(objst[nos], structF[nos]);
  for (E_Int nou = 0; nou < nu; nou++)
    RELEASESHAREDU(objut[nou], unstrF[nou], cnt[nou]);
  return l;
}
