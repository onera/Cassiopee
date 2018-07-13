/*    
    Copyright 2013-2018 Onera.

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

#include "generator.h"
using namespace K_FLD;
using namespace std;

//=============================================================================
/* Balance the octree */
//=============================================================================
PyObject* K_GENERATOR::balanceOctree(PyObject* self, PyObject* args)
{
  E_Int ratio;
  PyObject *octree;
#ifdef E_DOUBLEINT
  if (!PyArg_ParseTuple(args, "Ol", &octree, &ratio)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "Oi", &octree, &ratio)) return NULL;
#endif
  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(octree, varString, f, ni, nj, nk, cn, 
                                    eltType, true);
  if (res != 2) 
  {
    if (res == 1) RELEASESHAREDS(octree, f); 
    PyErr_SetString(PyExc_TypeError,
                    "balanceOctree: array must be unstructured.");
    return NULL;
  }
  if (strcmp(eltType, "HEXA") != 0 && strcmp(eltType, "QUAD") != 0)
  {
    RELEASESHAREDU(octree, f, cn); 
    PyErr_SetString(PyExc_TypeError,
                    "balanceOctree: the octree must be HEXA or QUAD.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(octree, f, cn); 
    PyErr_SetString(PyExc_TypeError,
                    "balanceOctree: the octree must contain coordinates.");
    return NULL;
  }
  posx++; posy++; posz++;
 
  if (ratio == 2) checkBalancing2(*cn, *f);
  else checkBalancing3(*cn, *f);
  PyObject* tpl = K_ARRAY::buildArray(*f, "x,y,z", *cn, -1, eltType);
  RELEASESHAREDU(octree, f, cn); 
  return tpl;
}
//=============================================================================
/* Verifie l'equilibrage du maillage engendre. Cas octree 2. */
//=============================================================================
void K_GENERATOR::checkBalancing2(FldArrayI& cn, FldArrayF& coords)
{
  //E_Float eps = 1.e-10;
  E_Int ind, ind1, ind2;
  FldArrayI cno = cn; FldArrayF fo = coords; // copies
  E_Int nvert = cn.getNfld();
  E_Int npts = coords.getSize();
  E_Int nelts = cn.getSize();
  E_Int no = 0; E_Int eto = 0;
  E_Int et2, nvoisins;
  E_Float dh1, dh2, dh1s2;
  E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2);
  E_Float* xt = coords.begin(1); E_Float* yt = coords.begin(2); E_Float* zt = coords.begin(3);
  E_Int splitl = 1;
  const char* eltType;
  if (nvert == 4) eltType = "QUAD";
  else eltType = "HEXA";
  while (splitl == 1)
  {
    npts = coords.getSize(); nelts = cn.getSize();
    xt = coords.begin(1); yt = coords.begin(2); zt = coords.begin(3);
    FldArrayF indict(nelts); indict.setAllValuesAtNull();
    FldArrayI indir(npts); indir.setAllValuesAt(-1);
    E_Int* indirp = indir.begin();
    //calcul de dh
    FldArrayF dht(nelts);
    E_Float* dhtp = dht.begin();
    for (E_Int et = 0; et < nelts; et++)
    {
      ind1 = cn1[et]-1; ind2 = cn2[et]-1;
      dhtp[et] = xt[ind2]-xt[ind1];
    }  
    // determination des elements voisins
    vector< vector<E_Int> > cEEN(nelts);
    E_Int ok = getNeighbourElts(npts, xt, yt, zt, cn, cEEN); 
    if (ok == -1) return;

    eto = 0; no = 0; splitl = 0;
    for (E_Int et1 = 0; et1 < nelts; et1++)
    {
      dh1 = dhtp[et1]; dh1s2 = 0.5*dh1;
      vector<E_Int>& voisins = cEEN[et1]; nvoisins = voisins.size();
      for (E_Int noet2 = 0; noet2 < nvoisins; noet2++)
      {
        et2 = voisins[noet2]; dh2 = dhtp[et2];
        if (dh1s2 > 1.1*dh2)
        {
          splitElement(et1, npts, nelts, indict[et1], indirp, cn, xt, yt, zt, 
                       cno, fo, indict, no, eto);
          splitl = 1; goto nextet1;
        }
      }
      // construction dans cno, fo
      for (E_Int nov = 1; nov <= nvert; nov++)
      {
        ind = cn(et1, nov)-1;
        if (indirp[ind] == -1) 
        {
          fo(no,1) = xt[ind]; fo(no,2) = yt[ind]; fo(no,3) = zt[ind];
          no++; 
          cno(eto, nov) = no; indirp[ind] = no;
          if (no+10 > fo.getSize()) fo.reAllocMat(fo.getSize()+npts, 3);
        }
        else cno(eto,nov) = indirp[ind];
      }
      eto++;
      if (eto+10 > indict.getSize()) indict.reAlloc(indict.getSize()+nelts);
      if (eto+10 > cno.getSize()) cno.reAllocMat(cno.getSize()+nelts, nvert);

      nextet1:;
    }//fin parcours et1

    // Realloc
    if (splitl == 1)
    {
      fo.reAllocMat(no,3); cno.reAllocMat(eto,nvert);
      cn = cno; coords = fo; cn1 = cn.begin(1); cn2 = cn.begin(2);
    }
  }
  K_CONNECT::cleanConnectivity(1, 2, 3, 1.e-10, eltType, coords, cn);
  return;
}
//=============================================================================
/* Verifie l'equilibrage du maillage engendre. Cas d'un 9/27-tree */
//=============================================================================
void K_GENERATOR::checkBalancing3(FldArrayI& cn, FldArrayF& coords)
{
  E_Int ind, ind1, ind2;
  E_Float eps = 1.e-10;
  E_Int split = 1;

  while (split == 1) 
  {
    FldArrayI cno = cn; FldArrayF fo = coords; // copies
    E_Int npts = coords.getSize(); E_Int nelts = cn.getSize(); E_Int nvert = cn.getNfld();
    FldArrayF indict(nelts); indict.setAllValuesAtNull();
    //calcul de dh
    FldArrayF dht(nelts);
    E_Float* dhtp = dht.begin();
    E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2);
    E_Float* xt = coords.begin(1); E_Float* yt = coords.begin(2); E_Float* zt = coords.begin(3);
    for (E_Int et = 0; et < nelts; et++)
    {
      ind1 = cn1[et]-1; ind2 = cn2[et]-1;
      dhtp[et] = xt[ind2]-xt[ind1];
    }  
    
    // determination des elements voisins
    vector< vector<E_Int> > cEEN(nelts);
    E_Int ok = getNeighbourElts(npts, xt, yt, zt, cn, cEEN); 
    if (ok == -1) return;
    E_Int no = 0; E_Int eto = 0;
    E_Int et2, nvoisins;
    E_Float dh1, dh2;
    split = 0;

    for (E_Int et1 = 0; et1 < nelts; et1++)
    {
      dh1 = dhtp[et1];
      vector<E_Int>& voisins = cEEN[et1]; nvoisins = voisins.size();
      for (E_Int noet2 = 0; noet2 < nvoisins; noet2++)
      {
        et2 = voisins[noet2]; dh2 = dhtp[et2];
        if (dh1 > 3.*dh2 + eps)
        {
          splitElement27(et1, npts, indict[et1], cn, xt, yt, zt, cno, fo, 
                         indict, no, eto);
          split = 1; goto nextet1;
        }
      }
      // construction dans cno, fo
       for (E_Int nov = 1; nov <= nvert; nov++)
       {
         ind = cn(et1, nov)-1;
         fo(no,1) = xt[ind]; fo(no,2) = yt[ind]; fo(no,3) = zt[ind];
         no++; cno(eto, nov) = no; 
         if (no+10 > fo.getSize()) fo.reAllocMat(fo.getSize()+npts,3);
       }
       eto++;
       if (eto + 28 > indict.getSize()) indict.reAllocMat(indict.getSize()+nelts,1);
       if (eto + 28 > cno.getSize()) cno.reAllocMat(cno.getSize()+nelts,nvert); 
      
      nextet1:;
    }
    // Nettoyage
    for (E_Int v = 0; v < nelts; v++) cEEN[v].clear();
    cEEN.clear();
    // Realloc
    fo.reAllocMat(no,3); cno.reAllocMat(eto,nvert); indict.reAllocMat(eto,1);
    if (nvert == 4) K_CONNECT::cleanConnectivity(1, 2, 3, 1.e-10, "QUAD", 
                                                 fo, cno);
    else K_CONNECT::cleanConnectivity(1, 2, 3, 1.e-10, "HEXA", fo, cno);
    coords = fo; cn = cno;
  }
  return;
}
