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
# include "generator.h"
#include "Def/DefTypes.h"

using namespace K_FLD;
using namespace std;

// ============================================================================
/* Deplace les points des maillages arrays sur les points 
   de la surface discretisee surface. */
// ============================================================================
PyObject* K_GENERATOR::snapSharpEdges(PyObject* self, PyObject* args)
{
  PyObject* arrays; PyObject* surface;
  if (!PYPARSETUPLE_(args, OO_, &arrays, &surface)) return NULL;

  // Extract infos from arrays
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltType;
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
                    "snapSharpEdges: invalid list of arrays.");
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nou = 0; nou < nu; nou++)
      RELEASESHAREDU(objut[nou], unstrF[nou], cnt[nou]);
    return NULL;
  }

  E_Int posx1, posy1, posz1, posi1;
  vector<E_Int> posxs; vector<E_Int> posys; vector<E_Int> poszs; 
  vector<E_Int> posis;
  vector<E_Int> posxu; vector<E_Int> posyu; vector<E_Int> poszu;
  vector<E_Int> posiu;
  for (E_Int nos = 0; nos < ns; nos++)
  {
    posx1 = K_ARRAY::isCoordinateXPresent(structVarString[nos]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(structVarString[nos]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(structVarString[nos]); posz1++;
    posi1 = K_ARRAY::isNamePresent("indic", structVarString[nos]); posi1++;
    if (posx1 == 0 || posy1 == 0 || posz1 == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "snapSharpEdges: arrays must contain coordinates.");
      for (E_Int nos = 0; nos < ns; nos++)
        RELEASESHAREDS(objst[nos], structF[nos]);
      for (E_Int nou = 0; nou < nu; nou++)
        RELEASESHAREDU(objut[nou], unstrF[nou], cnt[nou]);
      return NULL;
    }
    posxs.push_back(posx1); posys.push_back(posy1); poszs.push_back(posz1);
    posis.push_back(posi1);
  }
  for (E_Int nou = 0; nou < nu; nou++)
  {
    posx1 = K_ARRAY::isCoordinateXPresent(unstrVarString[nou]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(unstrVarString[nou]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(unstrVarString[nou]); posz1++;
    posi1 = K_ARRAY::isNamePresent("indic", unstrVarString[nou]); posi1++;
    if (posx1 == 0 || posy1 == 0 || posz1 == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "snapSharpEdges: arrays must contain coordinates.");
      for (E_Int nos = 0; nos < ns; nos++)
        RELEASESHAREDS(objst[nos], structF[nos]);
      for (E_Int nou = 0; nou < nu; nou++)
        RELEASESHAREDU(objut[nou], unstrF[nou], cnt[nou]);
      return NULL;
    }
    posxu.push_back(posx1); posyu.push_back(posy1); poszu.push_back(posz1);
    posiu.push_back(posi1);
  }
 
  // verification pour la surface
  E_Int im2, jm2, km2;
  FldArrayF* f2; FldArrayI* cn2;
  char* varString2; char* eltType2;
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
                    "snapSharpEdges: surface must be unstructured.");
    return NULL;
  }
  if (strcmp(eltType2, "TRI") != 0 && strcmp(eltType2, "BAR") != 0 
      && strcmp(eltType2, "NODE") != 0)
  {
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    RELEASESHAREDU(surface, f2, cn2);
    PyErr_SetString(PyExc_TypeError,
                    "snapSharpEdges: surface must be a TRI or BAR array.");
    return NULL;
  }

  E_Int posxu2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posyu2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int poszu2 = K_ARRAY::isCoordinateZPresent(varString2);
   
  if (posxu2 == -1 || posyu2 == -1 || poszu2 == -1)
  {
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    RELEASESHAREDU(surface, f2, cn2);
    PyErr_SetString(PyExc_TypeError,
                    "snapSharpEdges: coordinates not found in surface.");
    return NULL;
  }
  posxu2++; posyu2++; poszu2++;
  
  // Build arrays
  PyObject* l = PyList_New(0);
  vector<E_Float*> fs; vector<E_Float*> fu;
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

  // liste des coordonnees du maillage
  vector<E_Float*> coordx;
  vector<E_Float*> coordy;
  vector<E_Float*> coordz;
  vector<E_Float*> indic;
  vector<FldArrayI*> connect;
  // vecteur des tailles des zones
  vector<E_Int> sizet;

  // coordonnees du maillage pour les zones structures
  for (E_Int nos = 0; nos < ns; nos++)
  {
    E_Int npts = structF[nos]->getSize();
    E_Float* xp = fs[nos] + (posxs[nos]-1)*npts;
    E_Float* yp = fs[nos] + (posys[nos]-1)*npts;
    E_Float* zp = fs[nos] + (poszs[nos]-1)*npts;
    sizet.push_back(npts);
    
    // calcul des coords
    E_Float* xt = new E_Float [npts];
    E_Float* yt = new E_Float [npts];
    E_Float* zt = new E_Float [npts];
    for (E_Int i = 0; i < npts; i++)
    { 
      xt[i] = xp[i]; yt[i] = yp[i]; zt[i] = zp[i];
    }
    coordx.push_back(xt);
    coordy.push_back(yt);
    coordz.push_back(zt); 
    if (posis[nos] > 0) indic.push_back(structF[nos]->begin(posis[nos]));
    else indic.push_back(NULL);
    connect.push_back(NULL);
  }

  // coordonnees du maillage pour les zones non structures
  for (E_Int nou = 0; nou < nu; nou++)
  {
    E_Int npts = unstrF[nou]->getSize();
    E_Float* xp = fu[nou] + (posxu[nou]-1)*npts;
    E_Float* yp = fu[nou] + (posyu[nou]-1)*npts;
    E_Float* zp = fu[nou] + (poszu[nou]-1)*npts;
    sizet.push_back(npts);
      
    // calcul des coords
    E_Float* xt = new E_Float [npts];
    E_Float* yt = new E_Float [npts];
    E_Float* zt = new E_Float [npts];
    for (E_Int i = 0; i < npts; i++)
    { 
      xt[i] = xp[i]; yt[i] = yp[i]; zt[i] = zp[i];
    }
    coordx.push_back(xt);
    coordy.push_back(yt);
    coordz.push_back(zt);
    if (posiu[nou] > 0) indic.push_back(unstrF[nou]->begin(posiu[nou]));
    else indic.push_back(NULL);
    if (strcmp(eltType[nou], "NGON") == 0) connect.push_back(NULL);
    else connect.push_back(cnt[nou]);
  }

  // snap des maillages (structure et non structure)
  K_GENERATOR::snapMesh(*f2, sizet, posxu2, posyu2, poszu2,
                        coordx, coordy, coordz, indic, connect);

  // Patch back pour le structure
  for (E_Int nos = 0; nos < ns; nos++)
  {
    FldArrayF& f = *structF[nos];
    E_Int npts = f.getSize();
    E_Float* xp = fs[nos] + (posxs[nos]-1)*npts;
    E_Float* yp = fs[nos] + (posys[nos]-1)*npts;
    E_Float* zp = fs[nos] + (poszs[nos]-1)*npts;
    E_Float* indicp = fs[nos] + (posis[nos]-1)*npts;
    for (E_Int i = 0; i < npts; i++)
    {
      xp[i] = coordx[nos][i]; yp[i] = coordy[nos][i]; zp[i] = coordz[nos][i];
    }
    if (indic[nos] != NULL)
    {
      for (E_Int i = 0; i < npts; i++)
      {
        indicp[i] = indic[nos][i];
      }
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
    E_Float* indicp = fu[nou] + (posiu[nou]-1)*npts;
    for (E_Int i = 0; i < npts; i++)
    {
      xp[i] = coordx[ns+nou][i];
      yp[i] = coordy[ns+nou][i];
      zp[i] = coordz[ns+nou][i];
    }
    if (indic[ns+nou] != NULL)
    {
      for (E_Int i = 0; i < npts; i++)
      {
        indicp[i] = indic[ns+nou][i];
      }
    } 
  }

  // delete
  for (E_Int nos = 0; nos < ns; nos++)
  {
    delete [] coordx[nos];
    delete [] coordy[nos];
    delete [] coordz[nos];
  }
  for (E_Int nou = 0; nou < nu; nou++)
  {
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
/* snap un maillage par rapport a une surface.
   IN: surface: surface sur laquelle le maillage est snappe
   IN: sizet: vecteur des tailles des tableaux de (coordx,coordy,coordz) 
   de chaque zone
   IN: (posx,posy,posz): positions des coordonnees dans surface
   IN/OUT: (coordx,coordy,coordz): maillage a snapper 
   IN/OUT: indic: =1 sur les noeuds snappes */
// ============================================================================
void K_GENERATOR::snapMesh(
  FldArrayF& surface, vector<E_Int> sizet,
  E_Int posx, E_Int posy, E_Int posz,
  vector<E_Float*>& coordx, vector<E_Float*>& coordy, 
  vector<E_Float*>& coordz, vector<E_Float*>& indic,
  vector<FldArrayI*>& connect)
{
  E_Int n = sizet.size();  // nombre de zones a snapper
  E_Int totalSize = 0;     // taille globale pour les maillages
  E_Int s;                 // taille d'une zone donnee
  E_Int indp;              // indice du pt du front le plus proche du pt du contour
  E_Float pt[3];           // coordonnees d'un point du front

  for (E_Int no = 0; no < n; no++) {totalSize += sizet[no];}
  // Creation d un nouveau kdtree contenant les points du maillage
  // pour prendre en compte les modifications dues aux contours
  FldArrayF coords(totalSize, 4);
  E_Float* coords1 = coords.begin(1);
  E_Float* coords2 = coords.begin(2);
  E_Float* coords3 = coords.begin(3);
  E_Float* indic1 = coords.begin(4);
  FldArrayI index(totalSize, 2);
  E_Int* index1 = index.begin(1);
  E_Int* index2 = index.begin(2);

  // Ajoute pour controler le nbre de points projetes par element
  vector< vector< vector<E_Int> >* > connectVE;
  for (size_t i = 0; i < connect.size(); i++)
  {
    if (connect[i] != NULL)
    {
      //connectVE.push_back(NULL);
      vector< vector<E_Int> >* cVE = new vector< vector<E_Int> > (sizet[i]);
      K_CONNECT::connectEV2VE(*connect[i], *cVE);
      connectVE.push_back(cVE);
    }
    else connectVE.push_back(NULL);
  }
  E_Int c = 0;
  for (E_Int no = 0; no < n; no++)
  {
    s = sizet[no];
    if (indic[no] == NULL)
    {
      for (E_Int i = 0; i < s; i++)
      {
        coords1[c] = coordx[no][i];
        coords2[c] = coordy[no][i];
        coords3[c] = coordz[no][i];
        index1[c] = no; index2[c] = i;
        indic1[c] = 0;
        c++;
      }
    }
    else
    {
      for (E_Int i = 0; i < s; i++)
      {
        coords1[c] = coordx[no][i];
        coords2[c] = coordy[no][i];
        coords3[c] = coordz[no][i];
        indic1[c] = indic[no][i];
        index1[c] = no; index2[c] = i;
        c++;
      } 
    }
  }
  
  ArrayAccessor<FldArrayF> coordAccs(coords, 1, 2, 3);
  K_SEARCH::KdTree<FldArrayF> kdt(coordAccs);
        
  E_Float* fx = surface.begin(posx);
  E_Float* fy = surface.begin(posy);
  E_Float* fz = surface.begin(posz);
 
  E_Float dist1, dist2, xi, yi, zi;
  E_Int no, i, ind;
  for (E_Int indpt = 0; indpt < surface.getSize(); indpt++) 
  {
    pt[0] = fx[indpt]; pt[1] = fy[indpt]; pt[2] = fz[indpt];
    indp = kdt.getClosest(pt);
    if (indp != E_IDX_NONE)
    { 
      no = index1[indp]; i = index2[indp];
      xi = coordx[no][i]; yi = coordy[no][i]; zi = coordz[no][i];
      // check le nbre de pts projete par elements
      E_Int cptg = 0;
      vector< vector<E_Int> >* cVE = connectVE[no];
      if (cVE != NULL && indic[no] != NULL)
      {
        vector<E_Int>& elts = (*cVE)[i]; E_Int cpt = 0;
        FldArrayI& cEV = *connect[no];
        for (size_t e = 0; e < elts.size(); e++)
        {
          cpt = 0;
          for (E_Int j = 0; j < cEV.getNfld(); j++)
          {
            ind = cEV(elts[e],j+1)-1;
            if (ind != i) cpt += indic[no][ind];
          }
          cptg = K_FUNC::E_max(cpt, cptg);
        }
        //printf("pt " SF_D_ " -> cptg:" SF_D_ "\n",indp, cptg);
      }
      // distance si je snap ce point
      dist1 = (xi-fx[indpt])*(xi-fx[indpt]) +
        (yi-fy[indpt])*(yi-fy[indpt]) +
        (zi-fx[indpt])*(zi-fz[indpt]);
      // distance si le point a deja ete snappe
      dist2 = (coords1[indp]-xi)*(coords1[indp]-xi) +
        (coords2[indp]-yi)*(coords2[indp]-yi) +
        (coords3[indp]-zi)*(coords3[indp]-zi);
      //printf("indpt " SF_D_ " - dist1 %g, dist2=%g, cptg=" SF_D_ "\n", indpt, dist1, dist2, cptg); 
    
      // snap si le pt n'a jamais ete snappe ou si un pt plus proche existe
      if ((dist2 < 1.e-12 || dist1 <= dist2 || dist1 < 1.e-12) && cptg < 2)
      {
        coords1[indp] = fx[indpt];
        coords2[indp] = fy[indpt];
        coords3[indp] = fz[indpt];
        if (indic[no] != NULL) { indic1[indp] = 1.; indic[no][i] = 1.; }
        //printf("Le pt " SF_D_ " de la contrainte attire " SF_D_ " (" SF_F_ ")\n", indpt, indp, dist1);
      }
    }
    else printf("Can not snap point " SF_F3_ "\n", pt[0], pt[1], pt[2]);
  }
  c = 0;
  for (E_Int no = 0; no < n; no++)
  {
    s = sizet[no];
    if (indic[no] == NULL)
    {
      for (E_Int i = 0; i < s; i++)
      {
        coordx[no][i] = coords1[c];
        coordy[no][i] = coords2[c];
        coordz[no][i] = coords3[c]; c++;
      }
    }
    else
    {
      for (E_Int i = 0; i < s; i++)
      {
        coordx[no][i] = coords1[c];
        coordy[no][i] = coords2[c];
        coordz[no][i] = coords3[c];
        indic[no][i] = K_FUNC::E_max(indic1[c], indic[no][i]); c++;
      }
    }
  }

  // Compte les pts snappes par element
  for (E_Int no = 0; no < n; no++)
  {
    if (indic[no] != NULL)
    {
      E_Int cpt;
      FldArrayI& c = *connect[no];
      E_Int nfld = c.getNfld();
      for (E_Int i = 0; i < c.getSize(); i++)
      {
        cpt = 0;
        for (E_Int q = 0; q < nfld; q++) cpt += indic[no][c(i,q+1)-1]; 
        if (cpt >= 2) printf("warning: elt " SF_D_ " snaped more than 2 times\n", i);
      }
    }
  }

  // free cVE
  for (size_t i = 0; i < connectVE.size(); i++)
  {
    if (connectVE[i] != NULL) { delete connectVE[i]; }
  }
}
