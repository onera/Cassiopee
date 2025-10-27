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

# include "connector.h"
# include "Nuga/include/BbTree.h"
# include "Nuga/include/ArrayAccessor.h"
# include "Connect/connect.h"

using namespace std;
using namespace K_FLD;

struct facette 
{
    // facette structuree decoupee en 4 triangles
    // facette non structuree QUAD (cas prismatique/hexaedrique) decoupee en 4 triangles
    E_Int indA;
    E_Int indB;
    E_Int indC;
    E_Int indD;
    E_Float xo;
    E_Float yo;
    E_Float zo;
    E_Int indcell1;//cellule adjacente 1
    E_Int indcell2;//cellule adjacente 2
    E_Int nozone;// numero de la zone contenant la facette dans la liste des zones
};
struct trifacette 
{
    // facette triangulaire
    E_Float x1, y1, z1;
    E_Float x2, y2, z2;
    E_Float x3, y3, z3;
};
 

//=============================================================================
/* blankIntersectingCells: en structure, on teste uniquement les facettes 
   (i,k) et (j,k) 
   Attention: dans le cas non structure, ne prend pas en compte les elements
   non conformes. */
//=============================================================================
PyObject* K_CONNECTOR::blankIntersectingCells(PyObject* self, PyObject* args)
{
  PyObject *arrays, *cellnArrays;
  E_Float eps;
  if (!PYPARSETUPLE_(args, OO_ R_,
                    &arrays, &cellnArrays, &eps))
  {
      return NULL;
  }
  // Check every arrays
  if (PyList_Check(arrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "blankIntersectingCells: arrays argument must be a list.");
    return NULL;
  }
  // Check every arrays
  if (PyList_Check(arrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "blankIntersectingCells: arrays argument must be a list.");
    return NULL;
  }
  
  if (PyList_Check(cellnArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "blankIntersectingCells: cellnArrays argument must be a list.");
    return NULL;
  }
  
  E_Bool skipNoCoord = true; E_Bool skipDiffVars = true;
  E_Bool skipStructured = false; E_Bool skipUnstructured = false;
  E_Bool skipNoCoordc = false;
  // Extract infos from arrays
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt; vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltTypet, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int ns = structF.size(); E_Int nu = unstrF.size();
  if (isOk == -1) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "blankIntersectingCells: invalid list of arrays.");
    for (E_Int is = 0; is < ns; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    return NULL;   
  }
  // Extract infos from cellnArrays
  vector<E_Int> rescl;
  vector<char*> structVarStringc; vector<char*> unstrVarStringc; 
  vector<FldArrayF*> structFc; vector<FldArrayF*> unstrFc;
  vector<E_Int> nict; vector<E_Int> njct; vector<E_Int> nkct;
  vector<FldArrayI*> cnct; vector<char*> eltTypect;
  vector<PyObject*> objsct, objuct;
  isOk = K_ARRAY::getFromArrays(
    cellnArrays, rescl, structVarStringc, unstrVarStringc,
    structFc, unstrFc, nict, njct, nkct, cnct, eltTypect, objsct, objuct, 
    skipDiffVars, skipNoCoordc, skipStructured, skipUnstructured, true);
  E_Int nsc = structFc.size(); E_Int nuc = unstrFc.size();
  if (isOk == -1) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "blankIntersectingCells: invalid list of celln arrays.");
    for (E_Int is = 0; is < ns; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    for (E_Int is = 0; is < nsc; is++)
      RELEASESHAREDS(objsct[is], structFc[is]);
 
    for (E_Int is = 0; is < nu; is++)
      RELEASESHAREDU(objut[is], unstrF[is], cnt[is]);
    for (E_Int is = 0; is < nuc; is++)
      RELEASESHAREDU(objuct[is], unstrFc[is], cnct[is]);
    return NULL;   
  }
  
  E_Int nzones, ncelln;
  E_Int res = 0;
  if (structF.size() > 0 && unstrF.size() == 0) // structure uniquement  
  {
    nzones = structF.size(); ncelln = structFc.size();

    if (nzones != ncelln) 
    {
      PyErr_SetString(PyExc_TypeError, 
                      "blankIntersectingCells: cellnArrays and arrays must be of same size.");
      for (E_Int is = 0; is < nsc; is++)
        RELEASESHAREDS(objsct[is], structFc[is]);
      for (E_Int is = 0; is < ns; is++)
        RELEASESHAREDS(objst[is], structF[is]); 
      return NULL;
    }
    /* Verification des positions de x,y,z*/
    E_Int posx = K_ARRAY::isCoordinateXPresent(structVarString[0]);
    E_Int posy = K_ARRAY::isCoordinateYPresent(structVarString[0]);
    E_Int posz = K_ARRAY::isCoordinateZPresent(structVarString[0]);
    if (posx == -1 || posy == -1 || posz == -1) 
    {
      for (E_Int is = 0; is < ns; is++)
        RELEASESHAREDS(objst[is], structF[is]);
      for (E_Int is = 0; is < nsc; is++)
        RELEASESHAREDS(objsct[is], structFc[is]);
      PyErr_SetString(PyExc_TypeError,
                      "blankIntersectingCells: coordinates not found in array.");
      return NULL;
    }
    posx++; posy++; posz++;
    E_Int posc = K_ARRAY::isCellNatureField2Present(structVarStringc[0]);
    if (posc == -1)
    {
      for (E_Int is = 0; is < ns; is++)
        RELEASESHAREDS(objst[is], structF[is]);
      for (E_Int is = 0; is < nsc; is++)
        RELEASESHAREDS(objsct[is], structFc[is]);
      PyErr_SetString(PyExc_TypeError,
                      "blankIntersectingCells: celln variable not found.");
      return NULL;
    }
    posc++;
    res = blankIntersectingCellsStruct(eps, posx, posy, posz, posc, 
                                       nit, njt, nkt, structF, nict, njct, nkct, structFc);
    if (res == 0) 
    {
      for (E_Int is = 0; is < ns; is++)
        RELEASESHAREDS(objst[is], structF[is]);
      for (E_Int is = 0; is < nsc; is++)
        RELEASESHAREDS(objsct[is], structFc[is]);
      return NULL;
    }
    PyObject* l = PyList_New(0); 
    for (E_Int v = 0; v < nzones; v++)
    {
      E_Int api = structFc[v]->getApi();
      PyObject* tpl = K_ARRAY::buildArray3(*structFc[v], structVarStringc[v],
                                           nict[v], njct[v], nkct[v], api);
      PyList_Append(l, tpl); Py_DECREF(tpl);
    }
    for (E_Int is = 0; is < ns; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    for (E_Int is = 0; is < nsc; is++)
      RELEASESHAREDS(objsct[is], structFc[is]);
    
    return l;
  }// structure uniquement  
  else if (structF.size() == 0 && unstrF.size() > 0) // non structure uniqt
  {
    nzones = unstrF.size(); ncelln = unstrFc.size();
    if (nzones != ncelln) 
    {
      PyErr_SetString(PyExc_TypeError, 
                      "blankIntersectingCells: cellnArrays and arrays must be of same size.");
      for (E_Int is = 0; is < nu; is++)
        RELEASESHAREDU(objut[is], unstrF[is], cnt[is]);
      for (E_Int is = 0; is < nuc; is++)
        RELEASESHAREDU(objuct[is], unstrFc[is], cnct[is]);
      return NULL;
    }
    else
    { 
      /* Verification des positions de x,y,z */
      E_Int posx = K_ARRAY::isCoordinateXPresent(unstrVarString[0]);
      E_Int posy = K_ARRAY::isCoordinateYPresent(unstrVarString[0]);
      E_Int posz = K_ARRAY::isCoordinateZPresent(unstrVarString[0]);
      if (posx == -1 || posy == -1 || posz == -1) 
      {
        PyErr_SetString(PyExc_TypeError,
                        "blankIntersectingCells: coordinates not found in array.");
        for (E_Int is = 0; is < nu; is++)
          RELEASESHAREDU(objut[is], unstrF[is], cnt[is]);
        for (E_Int is = 0; is < nuc; is++)
          RELEASESHAREDU(objuct[is], unstrFc[is], cnct[is]);
        return NULL;
      }
      posx++; posy++; posz++;
      E_Int posc = K_ARRAY::isCellNatureField2Present(unstrVarStringc[0]);
      if (posc == -1)
      {
        PyErr_SetString(PyExc_TypeError,
                        "blankIntersectingCells: celln variable not found.");
        for (E_Int is = 0; is < nu; is++)
          RELEASESHAREDU(objut[is], unstrF[is], cnt[is]);
        for (E_Int is = 0; is < nuc; is++)
          RELEASESHAREDU(objuct[is], unstrFc[is], cnct[is]);
        return NULL;
      }
      posc++;
      E_Int type = 0; // 6 PENTA - 7 HEXA
      if (strcmp(eltTypect[0], "PENTA*") == 0) type = 6;
      else if (strcmp(eltTypect[0], "HEXA*") == 0) type = 7;
      else 
      {
        PyErr_SetString(PyExc_TypeError,
                        "blankIntersectingCells: elt type must be PENTA* or HEXA*.");
        for (E_Int is = 0; is < nu; is++)
          RELEASESHAREDU(objut[is], unstrF[is], cnt[is]);
        for (E_Int is = 0; is < nuc; is++)
          RELEASESHAREDU(objuct[is], unstrFc[is], cnct[is]);
        return NULL;
      }
      E_Int neltsc = eltTypect.size();
      for (E_Int v = 1; v < neltsc; v++)
      {
        if (type == 6 && strcmp(eltTypect[0], "PENTA*") != 0) 
        {
          PyErr_SetString(PyExc_TypeError,
                          "blankIntersectingCells: elt type must be PENTA*.");
          for (E_Int is = 0; is < nu; is++)
            RELEASESHAREDU(objut[is], unstrF[is], cnt[is]);
          for (E_Int is = 0; is < nuc; is++)
            RELEASESHAREDU(objuct[is], unstrFc[is], cnct[is]);
          return NULL;
        }
        else if (type == 7 && strcmp(eltTypect[0],"HEXA*") != 0) 
        {
          PyErr_SetString(PyExc_TypeError,
                          "blankIntersectingCells: elt type must be HEXA*.");
          for (E_Int is = 0; is < nu; is++)
            RELEASESHAREDU(objut[is], unstrF[is], cnt[is]);
          for (E_Int is = 0; is < nuc; is++)
            RELEASESHAREDU(objuct[is], unstrFc[is], cnct[is]);
          return NULL;
        }
      }
      if (type == 6) res = blankIntersectingCellsPenta(eps, posx, posy, posz, posc, 
                                                       cnt, unstrF, cnct, unstrFc);
      else res = blankIntersectingCellsHexa(eps, posx, posy, posz, posc, 
                                            cnt, unstrF, cnct, unstrFc);
      if (res == 0) 
      {
        for (E_Int is = 0; is < nu; is++)
          RELEASESHAREDU(objut[is], unstrF[is], cnt[is]);
        for (E_Int is = 0; is < nuc; is++)
          RELEASESHAREDU(objuct[is], unstrFc[is], cnct[is]);
        return NULL;
      }
      PyObject* l = PyList_New(0); 
      for (E_Int v = 0; v < nzones; v++)
      {
        E_Int api = unstrFc[v]->getApi();
        PyObject* tpl = K_ARRAY::buildArray3(*unstrFc[v], unstrVarStringc[v],
                                             *cnct[v], eltTypect[v], api);
        PyList_Append(l, tpl); Py_DECREF(tpl);
      }
      for (E_Int is = 0; is < nu; is++)
        RELEASESHAREDU(objut[is], unstrF[is], cnt[is]);
      for (E_Int is = 0; is < nuc; is++)
        RELEASESHAREDU(objuct[is], unstrFc[is], cnct[is]);
      return l;
    }
  }//fin  non structure
  else 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "blankIntersectingCells: cannot manage both structured and unstructured grids.");
    for (E_Int is = 0; is < ns; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    for (E_Int is = 0; is < nsc; is++)
      RELEASESHAREDS(objsct[is], structFc[is]);
    for (E_Int is = 0; is < nu; is++)
      RELEASESHAREDU(objut[is], unstrF[is], cnt[is]);
    for (E_Int is = 0; is < nuc; is++)
      RELEASESHAREDU(objuct[is], unstrFc[is], cnct[is]);
    return NULL;
  }
  return cellnArrays;
}

//=============================================================================
/* cas non structure prismatique : les facettes sont les facettes quad et 
   non triangulaires */
//=============================================================================
E_Int K_CONNECTOR::blankIntersectingCellsPenta( 
  E_Float eps, E_Int posx, E_Int posy, E_Int posz, E_Int posc,
  vector<FldArrayI*>& cnt, vector<FldArrayF*>& unstrF, 
  vector<FldArrayI*>& cnct, vector<FldArrayF*>& unstrFc)
{
  E_Int nzones = unstrF.size();
  for (E_Int v = 0; v < nzones; v++)
  {
    if (cnt[v]->getSize() != cnct[v]->getSize()) 
    {
      PyErr_SetString(PyExc_TypeError,
                      "blankIntersectingCells: coords and celln arrays are not consistent.");
      return 0;
    }
  }
  E_Int nk = getNbOfPentaLayers(unstrF[0]->getSize(), *cnt[0]);
  if (nk == 0) 
  { 
    PyErr_SetString(PyExc_TypeError,
                    "blankIntersectingCells: PENTA mesh must be structured in the k direction.");
    return 0;
  }
  E_Int ok = blankInvalidCellsPenta(eps, posx, posy, posz, posc, cnt, unstrF, cnct, unstrFc);
  if (ok == 0) return 0;

  /* construction de la liste des facettes non structurees */
  vector < vector< vector<E_Int> > > allCEEN(nzones);
  for (E_Int v = 0; v < nzones; v++)
  {
    FldArrayI* cn = cnt[v]; FldArrayF* f = unstrF[v];
    E_Int nelts = cn->getSize(); E_Int nvert = f->getSize();
    vector< vector<E_Int> > cEEN(nelts);
    E_Int ok = K_CONNECT::connectEV2EENbrs("PENTA", nvert, *cn, cEEN); 
    if (ok == 0) {printf("Error: connectivity elt/neigbour elts not valid.\n"); return 0;}
    allCEEN[v]=cEEN;
  }

  vector<facette*> listFacettes;
  E_Int indA, indB, indC, indD, etvoisin;
  for (E_Int k = 0; k < nk-1; k++)
  {
    for (E_Int v = 0; v < nzones; v++)
    {
      FldArrayF* f = unstrF[v]; FldArrayI* cn = cnt[v];
      E_Float* xt = f->begin(posx); E_Float* yt = f->begin(posy); E_Float* zt = f->begin(posz);
      E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2); E_Int* cn3 = cn->begin(3);
      E_Int* cn4 = cn->begin(4); E_Int* cn5 = cn->begin(5); E_Int* cn6 = cn->begin(6);
      E_Int nelts = cn->getSize();
      vector< vector<E_Int> >& cEEN = allCEEN[v];
      E_Int nelts0 = nelts/(nk-1);
      for (E_Int noet = 0; noet < nelts0; noet++)
      {
        E_Int et = noet + k * nelts0;
        //facettes ABED
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn5[et]-1; indD = cn4[et]-1;
        facette* face = new facette;
        face->indA = indA; face->indB = indB; face->indC = indC; face->indD = indD; 
        face->xo = 0.25*(xt[indA]+xt[indB]+xt[indC]+xt[indD]);
        face->yo = 0.25*(yt[indA]+yt[indB]+yt[indC]+yt[indD]);
        face->zo = 0.25*(zt[indA]+zt[indB]+zt[indC]+zt[indD]);
        etvoisin = K_CONNECT::getNbrForQuadFace(indA, indB, indC, indD, *cn, cEEN[et]);
        if ( etvoisin == -1 ) etvoisin = et;
        face->indcell1 = et; face->indcell2 = etvoisin;
        face->nozone = v;
        listFacettes.push_back(face);
        //facettes BCFE
        indA = cn2[et]-1; indB = cn3[et]-1; indC = cn6[et]-1; indD = cn5[et]-1;
        face = new facette;
        face->indA = indA; face->indB = indB; face->indC = indC; face->indD = indD; 
        face->xo = 0.25*(xt[indA]+xt[indB]+xt[indC]+xt[indD]);
        face->yo = 0.25*(yt[indA]+yt[indB]+yt[indC]+yt[indD]);
        face->zo = 0.25*(zt[indA]+zt[indB]+zt[indC]+zt[indD]);
        etvoisin = K_CONNECT::getNbrForQuadFace(indA, indB, indC, indD, *cn, cEEN[et]);
        if ( etvoisin == -1 ) etvoisin = et;
        face->indcell1 = et; face->indcell2 = etvoisin;
        face->nozone = v;
        listFacettes.push_back(face);      

        //facettes CADF
        indA = cn3[et]-1; indB = cn1[et]-1; indC = cn4[et]-1; indD = cn6[et]-1;
        face = new facette;
        face->indA = indA; face->indB = indB; face->indC = indC; face->indD = indD; 
        face->xo = 0.25*(xt[indA]+xt[indB]+xt[indC]+xt[indD]);
        face->yo = 0.25*(yt[indA]+yt[indB]+yt[indC]+yt[indD]);
        face->zo = 0.25*(zt[indA]+zt[indB]+zt[indC]+zt[indD]);
        etvoisin = K_CONNECT::getNbrForQuadFace(indA, indB, indC, indD, *cn, cEEN[et]);
        if ( etvoisin == -1 ) etvoisin = et;
        face->indcell1 = et; face->indcell2 = etvoisin;
        face->nozone = v;
        listFacettes.push_back(face);
      } 
    }
  }
  // preconditionnement par BBtree
  E_Int nfacettes = listFacettes.size();
  typedef K_SEARCH::BoundingBox<3>  BBox3DType; 
  E_Float minB[3];  E_Float maxB[3];
  FldArrayF bbox(nfacettes,6);
  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(2);
  E_Float* yminp = bbox.begin(3); E_Float* ymaxp = bbox.begin(4);
  E_Float* zminp = bbox.begin(5); E_Float* zmaxp = bbox.begin(6);
  vector<BBox3DType*> boxes(nfacettes); 
  E_Float xmin, xmax, ymin, ymax, zmin, zmax;
  E_Int noz = 0;
  for (E_Int v = 0; v < nfacettes; v++)
  {
    facette* face = listFacettes[v]; noz = face->nozone;
    indA = face->indA; indB = face->indB; indC = face->indC; indD = face->indD;
    E_Float* xt = unstrF[noz]->begin(posx); E_Float* yt = unstrF[noz]->begin(posy); E_Float* zt = unstrF[noz]->begin(posz);
    xmin = K_FUNC::E_min(xt[indA],xt[indB]); xmin = K_FUNC::E_min(xmin,xt[indC]); xmin = K_FUNC::E_min(xmin,xt[indD]);
    ymin = K_FUNC::E_min(yt[indA],yt[indB]); ymin = K_FUNC::E_min(ymin,yt[indC]); ymin = K_FUNC::E_min(ymin,yt[indD]);
    zmin = K_FUNC::E_min(zt[indA],zt[indB]); zmin = K_FUNC::E_min(zmin,zt[indC]); zmin = K_FUNC::E_min(zmin,zt[indD]);
    xmax = K_FUNC::E_max(xt[indA],xt[indB]); xmax = K_FUNC::E_max(xmax,xt[indC]); xmax = K_FUNC::E_max(xmax,xt[indD]);
    ymax = K_FUNC::E_max(yt[indA],yt[indB]); ymax = K_FUNC::E_max(ymax,yt[indC]); ymax = K_FUNC::E_max(ymax,yt[indD]);
    zmax = K_FUNC::E_max(zt[indA],zt[indB]); zmax = K_FUNC::E_max(zmax,zt[indC]); zmax = K_FUNC::E_max(zmax,zt[indD]);
    xminp[v] = xmin; yminp[v] = ymin; zminp[v] = zmin; xmaxp[v] = xmax; ymaxp[v] = ymax; zmaxp[v] = zmax;
    minB[0] = xmin-eps; minB[1] = ymin-eps; minB[2] = zmin-eps; 
    maxB[0] = xmax+eps; maxB[1] = ymax+eps; maxB[2] = zmax+eps; 
    boxes[v] = new BBox3DType(minB, maxB);
  }
  //======================
  // Build the box tree.
  //======================
  K_SEARCH::BbTree3D* bbtree = new K_SEARCH::BbTree3D(boxes);
  E_Int intersect = 0;
  E_Int etg, etd, etg2, etd2;
  E_Float ptA1[3]; E_Float ptB1[3]; E_Float ptC1[3]; E_Float ptD1[3]; E_Float ptO1[3];
  E_Float ptA2[3]; E_Float ptB2[3]; E_Float ptC2[3]; E_Float ptD2[3]; E_Float ptO2[3];
  vector<E_Int> indicesBB;
  E_Int noz1, noz2;
  for (E_Int v1 = 0; v1 < nfacettes; v1++)
  {
    minB[0] = xminp[v1]; minB[1] = yminp[v1]; minB[2] = zminp[v1];
    maxB[0] = xmaxp[v1]; maxB[1] = ymaxp[v1]; maxB[2] = zmaxp[v1]; 
    facette* face1 = listFacettes[v1]; noz1 = face1->nozone;
    etg = face1->indcell1; etd = face1->indcell2;
    E_Float* xt1 = unstrF[noz1]->begin(posx); 
    E_Float* yt1 = unstrF[noz1]->begin(posy); 
    E_Float* zt1 = unstrF[noz1]->begin(posz);
    E_Float* cellN = unstrFc[noz1]->begin(posc);
    ptA1[0] = xt1[face1->indA]; ptA1[1] = yt1[face1->indA]; ptA1[2] = zt1[face1->indA];
    ptB1[0] = xt1[face1->indB]; ptB1[1] = yt1[face1->indB]; ptB1[2] = zt1[face1->indB];
    ptC1[0] = xt1[face1->indC]; ptC1[1] = yt1[face1->indC]; ptC1[2] = zt1[face1->indC];
    ptD1[0] = xt1[face1->indD]; ptD1[1] = yt1[face1->indD]; ptD1[2] = zt1[face1->indD];
    ptO1[0] = face1->xo; ptO1[1] = face1->yo; ptO1[2] = face1->zo;

    indicesBB.clear();
    if ( cellN[etg] != 0. || cellN[etd] != 0. ) 
    {
      bbtree->getOverlappingBoxes(minB, maxB, indicesBB);
      E_Int nindicesBB = indicesBB.size();
      for (E_Int nov2 = 0; nov2 < nindicesBB; nov2++)
      {
        E_Int v2 = indicesBB[nov2];
        if ( v2 != v1 ) 
        {
          facette* face2 = listFacettes[v2]; noz2 = face2->nozone;
          E_Float* cellN2 = unstrFc[noz2]->begin(posc); 
          etg2 = face2->indcell1; etd2 = face2->indcell2;
          E_Float* xt2 = unstrF[noz2]->begin(posx); 
          E_Float* yt2 = unstrF[noz2]->begin(posy); 
          E_Float* zt2 = unstrF[noz2]->begin(posz);
          ptA2[0] = xt2[face2->indA]; ptA2[1] = yt2[face2->indA]; ptA2[2] = zt2[face2->indA];
          ptB2[0] = xt2[face2->indB]; ptB2[1] = yt2[face2->indB]; ptB2[2] = zt2[face2->indB];
          ptC2[0] = xt2[face2->indC]; ptC2[1] = yt2[face2->indC]; ptC2[2] = zt2[face2->indC];
          ptD2[0] = xt2[face2->indD]; ptD2[1] = yt2[face2->indD]; ptD2[2] = zt2[face2->indD];
          ptO2[0] = face2->xo; ptO2[1] = face2->yo; ptO2[2] = face2->zo;
          
          //A2B2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsPenta(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsPenta(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsPenta(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsPenta(etg2, etd2, *cnt[noz2], cellN2);

          //B2C2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptB2, ptC2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsPenta(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptB2, ptC2, ptO2, eps);        
          if ( intersect < 0 ) blankAboveCellsPenta(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptB2, ptC2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsPenta(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptB2, ptC2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsPenta(etg2, etd2, *cnt[noz2], cellN2);

          //C2D2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsPenta(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsPenta(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsPenta(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsPenta(etg2, etd2, *cnt[noz2], cellN2);

          //A2D2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsPenta(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsPenta(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsPenta(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsPenta(etg2, etd2, *cnt[noz2], cellN2);
        }
      }
    }
  }
  for (E_Int v = 0; v < nfacettes; v++) { delete listFacettes[v]; delete boxes[v]; }
  delete bbtree;
  return 1;
}

//=============================================================================
/* cas non structure hexa : les facettes sont les facettes normales aux 
   facettes 1234 et 5678  */
//=============================================================================
E_Int K_CONNECTOR::blankIntersectingCellsHexa( 
  E_Float eps, E_Int posx, E_Int posy, E_Int posz, E_Int posc,
  vector<FldArrayI*>& cnt, vector<FldArrayF*>& unstrF, 
  vector<FldArrayI*>& cnct, vector<FldArrayF*>& unstrFc)
{
  E_Int nzones = unstrF.size();
  for (E_Int v = 0; v < nzones; v++)
  {
    if ( cnt[v]->getSize() != cnct[v]->getSize() ) 
    {
      PyErr_SetString(PyExc_TypeError,
                      "blankIntersectingCells: coords and celln arrays are not consistent.");
      return 0;
    }
  }
  E_Int nk = getNbOfHexaLayers(unstrF[0]->getSize(), *cnt[0]);
  if ( nk == 0 ) 
  { 
    PyErr_SetString(PyExc_TypeError,
                    "blankIntersectingCells: HEXA mesh must be structured in the k direction.");
    return 0;
  }
  E_Int ok = blankInvalidCellsHexa(eps, posx, posy, posz, posc, cnt, unstrF, cnct, unstrFc);
  if ( ok == 0 ) return 0;

  /* construction de la liste des facettes non structurees */
  vector<facette*> listFacettes;
  E_Int indA, indB, indC, indD, etvoisin;
  vector < vector< vector<E_Int> > > allCEEN(nzones);
  for (E_Int v = 0; v < nzones; v++)
  {
    FldArrayI* cn = cnt[v]; FldArrayF* f = unstrF[v];
    E_Int nelts = cn->getSize(); E_Int nvert = f->getSize();
    vector< vector<E_Int> > cEEN(nelts);
    E_Int ok = K_CONNECT::connectEV2EENbrs("HEXA", nvert, *cn, cEEN); 
    if ( ok == 0 ) {printf("Error: connectivity elt/neigbour elts not valid.\n"); return 0;}
    allCEEN[v]=cEEN;
  }
  
  for (E_Int k = 0; k < nk-1; k++)
  {
    for (E_Int v = 0; v < nzones; v++)
    {
      FldArrayF* f = unstrF[v]; FldArrayI* cn = cnt[v];
      E_Float* xt = f->begin(posx); E_Float* yt = f->begin(posy); E_Float* zt = f->begin(posz);
      E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2); E_Int* cn3 = cn->begin(3); E_Int* cn4 = cn->begin(4); 
      E_Int* cn5 = cn->begin(5); E_Int* cn6 = cn->begin(6); E_Int* cn7 = cn->begin(7); E_Int* cn8 = cn->begin(8); 
      E_Int nelts = cn->getSize();
      vector< vector<E_Int> >& cEEN = allCEEN[v];
      E_Int nelts0 = nelts/(nk-1);
      for (E_Int noet = 0; noet < nelts0; noet++)
      {
        E_Int et = noet + k * nelts0;
        //facettes 1485
        indA = cn1[et]-1; indB = cn4[et]-1; indC = cn8[et]-1; indD = cn5[et]-1;
        facette* face = new facette;
        face->indA = indA; face->indB = indB; face->indC = indC; face->indD = indD; 
        face->xo = 0.25*(xt[indA]+xt[indB]+xt[indC]+xt[indD]);
        face->yo = 0.25*(yt[indA]+yt[indB]+yt[indC]+yt[indD]);
        face->zo = 0.25*(zt[indA]+zt[indB]+zt[indC]+zt[indD]);
        etvoisin = K_CONNECT::getNbrForQuadFace(indA, indB, indC, indD, *cn, cEEN[et]);
        if ( etvoisin == -1 ) etvoisin = et;
        face->indcell1 = et; face->indcell2 = etvoisin;
        face->nozone = v;
        listFacettes.push_back(face);
        //facettes 2376
        indA = cn2[et]-1; indB = cn3[et]-1; indC = cn7[et]-1; indD = cn6[et]-1;
        face = new facette;
        face->indA = indA; face->indB = indB; face->indC = indC; face->indD = indD; 
        face->xo = 0.25*(xt[indA]+xt[indB]+xt[indC]+xt[indD]);
        face->yo = 0.25*(yt[indA]+yt[indB]+yt[indC]+yt[indD]);
        face->zo = 0.25*(zt[indA]+zt[indB]+zt[indC]+zt[indD]);
        etvoisin = K_CONNECT::getNbrForQuadFace(indA, indB, indC, indD, *cn, cEEN[et]);
        if ( etvoisin == -1 ) etvoisin = et;
        face->indcell1 = et; face->indcell2 = etvoisin;
        face->nozone = v;
        listFacettes.push_back(face);      

        //facettes 1265
        indA = cn1[et]-1; indB = cn2[et]-1; indC = cn6[et]-1; indD = cn5[et]-1;
        face = new facette;
        face->indA = indA; face->indB = indB; face->indC = indC; face->indD = indD; 
        face->xo = 0.25*(xt[indA]+xt[indB]+xt[indC]+xt[indD]);
        face->yo = 0.25*(yt[indA]+yt[indB]+yt[indC]+yt[indD]);
        face->zo = 0.25*(zt[indA]+zt[indB]+zt[indC]+zt[indD]);
        etvoisin = K_CONNECT::getNbrForQuadFace(indA, indB, indC, indD, *cn, cEEN[et]);
        if ( etvoisin == -1 ) etvoisin = et;
        face->indcell1 = et; face->indcell2 = etvoisin;
        face->nozone = v;
        listFacettes.push_back(face);

        //facettes 4378
        indA = cn4[et]-1; indB = cn3[et]-1; indC = cn7[et]-1; indD = cn8[et]-1;
        face = new facette;
        face->indA = indA; face->indB = indB; face->indC = indC; face->indD = indD; 
        face->xo = 0.25*(xt[indA]+xt[indB]+xt[indC]+xt[indD]);
        face->yo = 0.25*(yt[indA]+yt[indB]+yt[indC]+yt[indD]);
        face->zo = 0.25*(zt[indA]+zt[indB]+zt[indC]+zt[indD]);
        etvoisin = K_CONNECT::getNbrForQuadFace(indA, indB, indC, indD, *cn, cEEN[et]);
        if ( etvoisin == -1 ) etvoisin = et;
        face->indcell1 = et; face->indcell2 = etvoisin;
        face->nozone = v;
        listFacettes.push_back(face);
      }
    }
  }
  
  // preconditionnement par BBtree
  E_Int nfacettes = listFacettes.size();
  typedef K_SEARCH::BoundingBox<3>  BBox3DType; 
  E_Float minB[3];  E_Float maxB[3];
  FldArrayF bbox(nfacettes,6);
  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(2);
  E_Float* yminp = bbox.begin(3); E_Float* ymaxp = bbox.begin(4);
  E_Float* zminp = bbox.begin(5); E_Float* zmaxp = bbox.begin(6);
  vector<BBox3DType*> boxes(nfacettes); 
  E_Float xmin, xmax, ymin, ymax, zmin, zmax;
  E_Int noz = 0;
  for (E_Int v = 0; v < nfacettes; v++)
  {
    facette* face = listFacettes[v]; noz = face->nozone;
    indA = face->indA; indB = face->indB; indC = face->indC; indD = face->indD;
    E_Float* xt = unstrF[noz]->begin(posx); E_Float* yt = unstrF[noz]->begin(posy); E_Float* zt = unstrF[noz]->begin(posz);
    xmin = K_FUNC::E_min(xt[indA],xt[indB]); xmin = K_FUNC::E_min(xmin,xt[indC]); xmin = K_FUNC::E_min(xmin,xt[indD]);
    ymin = K_FUNC::E_min(yt[indA],yt[indB]); ymin = K_FUNC::E_min(ymin,yt[indC]); ymin = K_FUNC::E_min(ymin,yt[indD]);
    zmin = K_FUNC::E_min(zt[indA],zt[indB]); zmin = K_FUNC::E_min(zmin,zt[indC]); zmin = K_FUNC::E_min(zmin,zt[indD]);
    xmax = K_FUNC::E_max(xt[indA],xt[indB]); xmax = K_FUNC::E_max(xmax,xt[indC]); xmax = K_FUNC::E_max(xmax,xt[indD]);
    ymax = K_FUNC::E_max(yt[indA],yt[indB]); ymax = K_FUNC::E_max(ymax,yt[indC]); ymax = K_FUNC::E_max(ymax,yt[indD]);
    zmax = K_FUNC::E_max(zt[indA],zt[indB]); zmax = K_FUNC::E_max(zmax,zt[indC]); zmax = K_FUNC::E_max(zmax,zt[indD]);
    xminp[v] = xmin; yminp[v] = ymin; zminp[v] = zmin; xmaxp[v] = xmax; ymaxp[v] = ymax; zmaxp[v] = zmax;
    minB[0] = xmin-eps; minB[1] = ymin-eps; minB[2] = zmin-eps; 
    maxB[0] = xmax+eps; maxB[1] = ymax+eps; maxB[2] = zmax+eps; 
    boxes[v] = new BBox3DType(minB, maxB);
  }
  //======================
  // Build the box tree.
  //======================
  K_SEARCH::BbTree3D* bbtree = new K_SEARCH::BbTree3D(boxes);
  E_Int intersect = 0;
  E_Int etg, etd, etg2, etd2;
  E_Float ptA1[3]; E_Float ptB1[3]; E_Float ptC1[3]; E_Float ptD1[3]; E_Float ptO1[3];
  E_Float ptA2[3]; E_Float ptB2[3]; E_Float ptC2[3]; E_Float ptD2[3]; E_Float ptO2[3];
  vector<E_Int> indicesBB;
  E_Int noz1, noz2;
  for (E_Int v1 = 0; v1 < nfacettes; v1++)
  {
    minB[0] = xminp[v1]; minB[1] = yminp[v1]; minB[2] = zminp[v1];
    maxB[0] = xmaxp[v1]; maxB[1] = ymaxp[v1]; maxB[2] = zmaxp[v1]; 
    facette* face1 = listFacettes[v1]; noz1 = face1->nozone;
    etg = face1->indcell1; etd = face1->indcell2;
    E_Float* xt1 = unstrF[noz1]->begin(posx); 
    E_Float* yt1 = unstrF[noz1]->begin(posy); 
    E_Float* zt1 = unstrF[noz1]->begin(posz);
    E_Float* cellN = unstrFc[noz1]->begin(posc);
    ptA1[0] = xt1[face1->indA]; ptA1[1] = yt1[face1->indA]; ptA1[2] = zt1[face1->indA];
    ptB1[0] = xt1[face1->indB]; ptB1[1] = yt1[face1->indB]; ptB1[2] = zt1[face1->indB];
    ptC1[0] = xt1[face1->indC]; ptC1[1] = yt1[face1->indC]; ptC1[2] = zt1[face1->indC];
    ptD1[0] = xt1[face1->indD]; ptD1[1] = yt1[face1->indD]; ptD1[2] = zt1[face1->indD];
    ptO1[0] = face1->xo; ptO1[1] = face1->yo; ptO1[2] = face1->zo;

    indicesBB.clear();
    if ( cellN[etg] != 0. || cellN[etd] != 0. ) 
    {
      bbtree->getOverlappingBoxes(minB, maxB, indicesBB);
      E_Int nindicesBB = indicesBB.size();
      for (E_Int nov2 = 0; nov2 < nindicesBB; nov2++)
      {
        E_Int v2 = indicesBB[nov2];
        if ( v2 != v1 ) 
        {
          facette* face2 = listFacettes[v2]; noz2 = face2->nozone;
          E_Float* cellN2 = unstrFc[noz2]->begin(posc); 
          etg2 = face2->indcell1; etd2 = face2->indcell2;
          E_Float* xt2 = unstrF[noz2]->begin(posx); 
          E_Float* yt2 = unstrF[noz2]->begin(posy); 
          E_Float* zt2 = unstrF[noz2]->begin(posz);
          ptA2[0] = xt2[face2->indA]; ptA2[1] = yt2[face2->indA]; ptA2[2] = zt2[face2->indA];
          ptB2[0] = xt2[face2->indB]; ptB2[1] = yt2[face2->indB]; ptB2[2] = zt2[face2->indB];
          ptC2[0] = xt2[face2->indC]; ptC2[1] = yt2[face2->indC]; ptC2[2] = zt2[face2->indC];
          ptD2[0] = xt2[face2->indD]; ptD2[1] = yt2[face2->indD]; ptD2[2] = zt2[face2->indD];
          ptO2[0] = face2->xo; ptO2[1] = face2->yo; ptO2[2] = face2->zo;

          // A2B2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsHexa(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsHexa(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsHexa(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsHexa(etg2, etd2, *cnt[noz2], cellN2);

          //B2C2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptB2, ptC2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsHexa(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptB2, ptC2, ptO2, eps);        
          if ( intersect < 0 ) blankAboveCellsHexa(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptB2, ptC2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsHexa(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptB2, ptC2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsHexa(etg2, etd2, *cnt[noz2], cellN2);

          //C2D2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsHexa(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsHexa(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsHexa(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsHexa(etg2, etd2, *cnt[noz2], cellN2);

          //A2D2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsHexa(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsHexa(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsHexa(etg2, etd2, *cnt[noz2], cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsHexa(etg2, etd2, *cnt[noz2], cellN2);
          }
      }
    }
  }
  for (E_Int v = 0; v < nfacettes; v++) { delete listFacettes[v]; delete boxes[v]; }
  delete bbtree;
  return 1;
}
//=============================================================================
/* cas structure : les facettes sont les facettes normales au plan (i,j)*/
//=============================================================================
E_Int K_CONNECTOR::blankIntersectingCellsStruct( 
  E_Float eps, E_Int posx, E_Int posy, E_Int posz, E_Int posc,
  vector<E_Int>& nit, vector<E_Int>& njt, vector<E_Int>& nkt, vector<FldArrayF*>& structF, 
  vector<E_Int>& nict, vector<E_Int>& njct, vector<E_Int>& nkct, vector<FldArrayF*>& structFc)
{
  E_Int nzones = structF.size();
  for (E_Int v = 0; v < nzones; v++)
  {
    if (K_FUNC::E_max(1,nit[v]-1) != nict[v] || 
        K_FUNC::E_max(1,njt[v]-1) != njct[v] || 
        K_FUNC::E_max(1,nkt[v]-1) != nkct[v])
    {
      PyErr_SetString(PyExc_TypeError,
                      "blankIntersectingCells: coords and celln arrays are not consistent.");
      return 0;
    }
  }
  E_Int ok = blankInvalidCellsStruct(eps, posx, posy, posz, posc, nit, njt, nkt, structF, 
                                     nict, njct, nkct, structFc);
  if ( ok == 0 ) return 0;

  /* construction de la liste des facettes structurees */
  vector<facette*> listFacettes;
  E_Int indA, indB, indC, indD;
  E_Int indcellg, indcelld, indcell;
  E_Int incd = 0; E_Int incg = 0;  
  E_Int nk = nkt[0];
  for (E_Int k = 0; k < nk-1; k++)
  {
    for (E_Int v = 0; v < nzones; v++)
    {
      FldArrayF* f = structF[v];
      E_Int ni = nit[v]; E_Int nj = njt[v]; E_Int ninj = ni*nj; 
      E_Int nic = nict[v]; E_Int njc = njct[v]; E_Int nicnjc = nic*njc;
      E_Float* xt = f->begin(posx); E_Float* yt = f->begin(posy); E_Float* zt = f->begin(posz);
      //facettes en i
      for (E_Int i = 0; i < ni; i++)
      {
        incg = -1; incd = 0; if (i == 0) incg = 0; if (i == ni-1) incd = -1;  
        for (E_Int j = 0; j < nj-1; j++)
        {
          indA = i + j * ni + k * ninj; indB = indA+ni; indC = indB + ninj; indD = indA + ninj;
          facette* face = new facette;
          indcell = i + j * nic + k * nicnjc;
          indcelld = indcell + incd; indcellg = indcell + incg; 
          face->indA = indA; face->indB = indB; face->indC = indC; face->indD = indD; 
          face->xo = 0.25*(xt[indA]+xt[indB]+xt[indC]+xt[indD]);
          face->yo = 0.25*(yt[indA]+yt[indB]+yt[indC]+yt[indD]);
          face->zo = 0.25*(zt[indA]+zt[indB]+zt[indC]+zt[indD]);
          face->indcell1 = indcellg; face->indcell2 = indcelld;
          face->nozone = v;
          listFacettes.push_back(face);
        }
      }
      
      //facettes en j 
      for (E_Int j = 0; j < nj; j++)
      {
        incg = -nic; incd = 0; if (j == 0) incg = 0; if (j == nj-1) incd = -nic;  
        for (E_Int i = 0; i < ni-1; i++)
        {
          indA = i + j * ni + k * ninj; indB = indA+1; indC = indB + ninj; indD = indA + ninj;
          facette* face = new facette;
          indcell = i + j * nic + k * nicnjc; 
          indcelld = indcell+incd; indcellg = indcell + incg; 
          face->indA = indA; face->indB = indB; face->indC = indC; face->indD = indD; 
          face->xo = 0.25*(xt[indA]+xt[indB]+xt[indC]+xt[indD]);
          face->yo = 0.25*(yt[indA]+yt[indB]+yt[indC]+yt[indD]);
          face->zo = 0.25*(zt[indA]+zt[indB]+zt[indC]+zt[indD]);
          face->indcell1 = indcellg; face->indcell2 = indcelld;
          face->nozone = v; listFacettes.push_back(face);
        }
      }
    }// boucle en k 
  }
  // preconditionnement par BBtree
  E_Int nfacettes = listFacettes.size();
  typedef K_SEARCH::BoundingBox<3> BBox3DType; 
  E_Float minB[3];  E_Float maxB[3];
  FldArrayF bbox(nfacettes,6);
  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(2);
  E_Float* yminp = bbox.begin(3); E_Float* ymaxp = bbox.begin(4);
  E_Float* zminp = bbox.begin(5); E_Float* zmaxp = bbox.begin(6);
  vector<BBox3DType*> boxes(nfacettes); 
  E_Float xmin, xmax, ymin, ymax, zmin, zmax;
  E_Int noz = 0;
  for (E_Int v = 0; v < nfacettes; v++)
  {
    facette* face = listFacettes[v]; noz = face->nozone;
    indA = face->indA; indB = face->indB; indC = face->indC; indD = face->indD;
    E_Float* xt = structF[noz]->begin(posx); E_Float* yt = structF[noz]->begin(posy); E_Float* zt = structF[noz]->begin(posz);
    xmin = K_FUNC::E_min(xt[indA],xt[indB]); xmin = K_FUNC::E_min(xmin,xt[indC]); xmin = K_FUNC::E_min(xmin,xt[indD]);
    ymin = K_FUNC::E_min(yt[indA],yt[indB]); ymin = K_FUNC::E_min(ymin,yt[indC]); ymin = K_FUNC::E_min(ymin,yt[indD]);
    zmin = K_FUNC::E_min(zt[indA],zt[indB]); zmin = K_FUNC::E_min(zmin,zt[indC]); zmin = K_FUNC::E_min(zmin,zt[indD]);
    xmax = K_FUNC::E_max(xt[indA],xt[indB]); xmax = K_FUNC::E_max(xmax,xt[indC]); xmax = K_FUNC::E_max(xmax,xt[indD]);
    ymax = K_FUNC::E_max(yt[indA],yt[indB]); ymax = K_FUNC::E_max(ymax,yt[indC]); ymax = K_FUNC::E_max(ymax,yt[indD]);
    zmax = K_FUNC::E_max(zt[indA],zt[indB]); zmax = K_FUNC::E_max(zmax,zt[indC]); zmax = K_FUNC::E_max(zmax,zt[indD]);
    xminp[v] = xmin; yminp[v] = ymin; zminp[v] = zmin; xmaxp[v] = xmax; ymaxp[v] = ymax; zmaxp[v] = zmax;
    minB[0] = xmin-eps; minB[1] = ymin-eps; minB[2] = zmin-eps; 
    maxB[0] = xmax+eps; maxB[1] = ymax+eps; maxB[2] = zmax+eps; 
    boxes[v] = new BBox3DType(minB, maxB);
  }
  //======================
  // Build the box tree.
  //======================
  K_SEARCH::BbTree3D* bbtree = new K_SEARCH::BbTree3D(boxes);
  E_Int intersect = 0;
  E_Int indcellg2, indcelld2;
  E_Float ptA1[3]; E_Float ptB1[3]; E_Float ptC1[3]; E_Float ptD1[3]; E_Float ptO1[3];
  E_Float ptA2[3]; E_Float ptB2[3]; E_Float ptC2[3]; E_Float ptD2[3]; E_Float ptO2[3];
  vector<E_Int> indicesBB;
  E_Int noz1, noz2;
  for (E_Int v1 = 0; v1 < nfacettes; v1++)
  {
    minB[0] = xminp[v1]; minB[1] = yminp[v1]; minB[2] = zminp[v1];
    maxB[0] = xmaxp[v1]; maxB[1] = ymaxp[v1]; maxB[2] = zmaxp[v1]; 

    facette* face1 = listFacettes[v1]; noz1 = face1->nozone;
    indcellg = face1->indcell1; indcelld = face1->indcell2;

    E_Float* xt1 = structF[noz1]->begin(posx); 
    E_Float* yt1 = structF[noz1]->begin(posy); 
    E_Float* zt1 = structF[noz1]->begin(posz);
    E_Float* cellN = structFc[noz1]->begin(posc);
    ptA1[0] = xt1[face1->indA]; ptA1[1] = yt1[face1->indA]; ptA1[2] = zt1[face1->indA];
    ptB1[0] = xt1[face1->indB]; ptB1[1] = yt1[face1->indB]; ptB1[2] = zt1[face1->indB];
    ptC1[0] = xt1[face1->indC]; ptC1[1] = yt1[face1->indC]; ptC1[2] = zt1[face1->indC];
    ptD1[0] = xt1[face1->indD]; ptD1[1] = yt1[face1->indD]; ptD1[2] = zt1[face1->indD];
    ptO1[0] = face1->xo; ptO1[1] = face1->yo; ptO1[2] = face1->zo;     
    indicesBB.clear();

    if ( cellN[indcellg] != 0. || cellN[indcelld] != 0. ) 
    {
      bbtree->getOverlappingBoxes(minB, maxB, indicesBB);
      E_Int nindicesBB = indicesBB.size();
      for (E_Int nov2 = 0; nov2 < nindicesBB; nov2++)
      {
        E_Int v2 = indicesBB[nov2];
        if ( v2 > v1)
        {
          facette* face2 = listFacettes[v2]; noz2 = face2->nozone;
          E_Float* cellN2 = structFc[noz2]->begin(posc);
          indcellg2 = face2->indcell1; indcelld2 = face2->indcell2;
          E_Float* xt2 = structF[noz2]->begin(posx); 
          E_Float* yt2 = structF[noz2]->begin(posy); 
          E_Float* zt2 = structF[noz2]->begin(posz);
          ptA2[0] = xt2[face2->indA]; ptA2[1] = yt2[face2->indA]; ptA2[2] = zt2[face2->indA];
          ptB2[0] = xt2[face2->indB]; ptB2[1] = yt2[face2->indB]; ptB2[2] = zt2[face2->indB];
          ptC2[0] = xt2[face2->indC]; ptC2[1] = yt2[face2->indC]; ptC2[2] = zt2[face2->indC];
          ptD2[0] = xt2[face2->indD]; ptD2[1] = yt2[face2->indD]; ptD2[2] = zt2[face2->indD];
          ptO2[0] = face2->xo; ptO2[1] = face2->yo; ptO2[2] = face2->zo;
          E_Int nic2 = nict[noz2]; E_Int njc2 = njct[noz2]; E_Int nkc2 = nkct[noz2]; 
          E_Int nic2njc2 = nic2*njc2;
   
          //A2B2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsStruct(indcellg2, indcelld2, nic2, njc2, nkc2, nic2njc2, cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsStruct(indcellg2, indcelld2, nic2, njc2, nkc2, nic2njc2, cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsStruct(indcellg2, indcelld2, nic2, njc2, nkc2, nic2njc2, cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsStruct(indcellg2, indcelld2, nic2, njc2, nkc2, nic2njc2, cellN2);

          //B2C2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptB2, ptC2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsStruct(indcellg2, indcelld2, nic2, njc2, nkc2, nic2njc2, cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptB2, ptC2, ptO2, eps);       
          if ( intersect < 0 ) blankAboveCellsStruct(indcellg2, indcelld2, nic2, njc2, nkc2, nic2njc2, cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptB2, ptC2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsStruct(indcellg2, indcelld2, nic2, njc2, nkc2, nic2njc2, cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptB2, ptC2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsStruct(indcellg2, indcelld2, nic2, njc2, nkc2, nic2njc2, cellN2);

          //C2D2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsStruct(indcellg2, indcelld2, nic2, njc2, nkc2, nic2njc2, cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsStruct(indcellg2, indcelld2, nic2, njc2, nkc2, nic2njc2, cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsStruct(indcellg2, indcelld2, nic2, njc2, nkc2, nic2njc2, cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsStruct(indcellg2, indcelld2, nic2, njc2, nkc2, nic2njc2, cellN2);

          //A2D2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsStruct(indcellg2, indcelld2, nic2, njc2, nkc2, nic2njc2, cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsStruct(indcellg2, indcelld2, nic2, njc2, nkc2, nic2njc2, cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsStruct(indcellg2, indcelld2, nic2, njc2, nkc2, nic2njc2, cellN2);

          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) blankAboveCellsStruct(indcellg2, indcelld2, nic2, njc2, nkc2, nic2njc2, cellN2);
        }
      }
    }
  }
  for (E_Int v = 0; v < nfacettes; v++) { delete listFacettes[v]; delete boxes[v]; }
  delete bbtree;
  return 1;
}
//=====================================================================================
/* Passer le cellN a 0 pour les points en k au dessus des points indcell1 et indcell2*/
//=====================================================================================
E_Int K_CONNECTOR::blankAboveCellsStruct(E_Int indcell1, E_Int indcell2, 
                                         E_Int nic, E_Int njc, E_Int nkc, E_Int nicnjc, 
                                         E_Float* cellN)
{
  E_Int k1 = indcell1/nicnjc;  E_Int k2 = indcell2/nicnjc; 
  if ( k1 != k2 ) 
  {printf("Error: k1 must be equal to k2.\n"); return 0;}
  E_Int indp1, indp2;
  E_Int c = 0;

  for (E_Int k = k1; k < nkc; k++)
  {
    indp1 = indcell1 + c*nicnjc; indp2 = indcell2 + c*nicnjc;
    cellN[indp1] = 0.; cellN[indp2] = 0.; c++;
  }
  return 1;
}
//=====================================================================================
/* Passer le cellN a 0 pour les elements situes au dessus des elts etg et etd */
//=====================================================================================
E_Int K_CONNECTOR::blankAboveCellsPenta(E_Int etg, E_Int etd, FldArrayI& cn,
                                        E_Float* cellN)
{
  E_Int nelts = cn.getSize();
  E_Int ind4, ind5, ind6;
  E_Int indp1, indp2, indp3;
  E_Int* cn1 = cn.begin(1);  E_Int* cn2 = cn.begin(2); E_Int* cn3 = cn.begin(3);
  E_Int* cn4 = cn.begin(4);  E_Int* cn5 = cn.begin(5); E_Int* cn6 = cn.begin(6);

  FldArrayIS dejaVu(nelts);
  // elt de gauche + voisins au dessus
  ind4 = cn4[etg]-1; ind5 = cn5[etg]-1; ind6 = cn6[etg]-1;
  dejaVu.setAllValuesAtNull(); dejaVu[etg] = 1;dejaVu[etd] = 1;
  cellN[etg] = 0.; cellN[etd] = 0.;
  E_Int et0 = etg;

  restart:;
  for (E_Int et = et0; et < nelts; et++)
  {
    if ( dejaVu[et] == 0 ) 
    {
      indp1 = cn1[et]-1; indp2 = cn2[et]-1; indp3 = cn3[et]-1;
      if ( indp1 == ind4 || indp1 == ind5 || indp1 == ind6 ) 
      {
        if ( indp2 == ind4 || indp2 == ind5 || indp2 == ind6 ) 
        {
          if ( indp3 == ind4 || indp3 == ind5 || indp3 == ind6 )
          {
            cellN[et] = 0.; dejaVu[et] = 1;  et0 = et;
            ind4 = cn4[et]-1; ind5 = cn5[et]-1; ind6 = cn6[et]-1;
            goto restart;
          }
        }
      }
    }
  }
  if ( etd == etg ) return 1;

  // elt de droite + voisins au dessus
  ind4 = cn4[etd]-1; ind5 = cn5[etd]-1; ind6 = cn6[etd]-1;
  et0 = etd;
  restart2:;
  for (E_Int et = 0; et < nelts; et++)
  {
    if ( dejaVu[et] == 0 ) 
    {
      indp1 = cn1[et]-1; indp2 = cn2[et]-1; indp3 = cn3[et]-1;
      if ( indp1 == ind4 || indp1 == ind5 || indp1 == ind6 ) 
      {
        if ( indp2 == ind4 || indp2 == ind5 || indp2 == ind6 ) 
        {
          if ( indp3 == ind4 || indp3 == ind5 || indp3 == ind6 )
          {
            cellN[et] = 0.; dejaVu[et] = 1; et0 = et;
            ind4 = cn4[et]-1; ind5 = cn5[et]-1; ind6 = cn6[et]-1;
            goto restart2;
          }
        }
      }
    }
  }

  return 1;
}
//=====================================================================================
/* Passer le cellN a 0 pour les elements situes au dessus des elts etg et etd */
//=====================================================================================
E_Int K_CONNECTOR::blankAboveCellsHexa(E_Int etg, E_Int etd, FldArrayI& cn,
                                       E_Float* cellN)
{
  E_Int nelts = cn.getSize();
  E_Int ind5, ind6, ind7, ind8;
  E_Int indp1, indp2, indp3, indp4;
  FldArrayIS dejaVu(nelts);
  E_Int* cn1 = cn.begin(1); E_Int* cn2 = cn.begin(2); E_Int* cn3 = cn.begin(3); E_Int* cn4 = cn.begin(4); 
  E_Int* cn5 = cn.begin(5); E_Int* cn6 = cn.begin(6); E_Int* cn7 = cn.begin(7); E_Int* cn8 = cn.begin(8); 

  // elt de gauche + voisins au dessus
  ind5 = cn5[etg]-1; ind6 = cn6[etg]-1; ind7 = cn7[etg]-1; ind8 = cn8[etg]-1;
  dejaVu.setAllValuesAtNull(); dejaVu[etg] = 1; dejaVu[etd] = 1;
  cellN[etg] = 0.; cellN[etd] = 0.;
  E_Int et0 = etg;
  restart:;
  for (E_Int et = et0; et < nelts; et++)
  {
    if ( dejaVu[et] == 0 ) 
    {
      indp1 = cn1[et]-1; indp2 = cn2[et]-1; indp3 = cn3[et]-1; indp4 = cn4[et]-1;
      if ( indp1 == ind5 || indp1 == ind6 || indp1 == ind7 || indp1 == ind8) 
      {
        if ( indp2 == ind5 || indp2 == ind6 || indp2 == ind7 || indp2 == ind8) 
        {
          if ( indp3 == ind5 || indp3 == ind6 || indp3 == ind7 || indp3 == ind8)
          { 
            if ( indp4 == ind5 || indp4 == ind6 || indp4 == ind7 || indp4 == ind8) 
            {
              cellN[et] = 0; dejaVu[et] = 1; et0 = et;
              ind5 = cn5[et]-1; ind6 = cn6[et]-1; ind7 = cn7[et]-1; ind8 = cn8[et]-1;
              goto restart;
            }
          }
        }
      }
    }
  }
  if ( etd == etg ) return 1;
  // elt de droite + voisins au dessus
  ind5 = cn5[etd]-1; ind6 = cn6[etd]-1; ind7 = cn7[etd]-1; ind8 = cn8[etd]-1;
  et0 = etd;
  restart2:;
  for (E_Int et = et0; et < nelts; et++)
  {
    if ( dejaVu[et] == 0 ) 
    {
      indp1 = cn1[et]-1; indp2 = cn2[et]-1; indp3 = cn3[et]-1; indp4 = cn4[et]-1;
      if ( indp1 == ind5 || indp1 == ind6 || indp1 == ind7 || indp1 == ind8) 
      {
        if ( indp2 == ind5 || indp2 == ind6 || indp2 == ind7 || indp2 == ind8) 
        {
          if ( indp3 == ind5 || indp3 == ind6 || indp3 == ind7 || indp3 == ind8)
          { 
            if ( indp4 == ind5 || indp4 == ind6 || indp4 == ind7 || indp4 == ind8) 
            {
              cellN[et] = 0; dejaVu[et] = 1; et0 = et;
              ind5 = cn5[et]-1; ind6 = cn6[et]-1; ind7 = cn7[et]-1; ind8 = cn8[et]-1;
              goto restart2;
            }
          }
        }
      }
    }
  }

  return 1;
}
//=============================================================================
/* blanking cells with self-intersecting faces and negative volume */
//=============================================================================
E_Int K_CONNECTOR::blankInvalidCellsPenta( 
  E_Float eps, E_Int posx, E_Int posy, E_Int posz, E_Int posc,
  vector<FldArrayI*>& cnt, vector<FldArrayF*>& unstrF, 
  vector<FldArrayI*>& cnct, vector<FldArrayF*>& unstrFc)
{
  E_Int nzones = unstrF.size();
  // blanking des cellules de volume negatif
  for (E_Int v = 0; v < nzones; v++)
  {
    E_Int npts = unstrF[v]->getSize();
    FldArrayF coord(npts, 3);
    coord.setOneField(*unstrF[v], posx, 1);
    coord.setOneField(*unstrF[v], posy, 2);
    coord.setOneField(*unstrF[v], posz, 3);

    FldArrayI& cm = *(cnt[v]->getConnect(0));
    E_Int nelts = cm.getSize();
    E_Int nfacets = 5*nelts;
    FldArrayF snx(nfacets), sny(nfacets), snz(nfacets), surf(nfacets);
    FldArrayF vol(nelts);
    K_METRIC::compMetricUnstruct(
      *cnt[v], "PENTA", coord.begin(1), coord.begin(2), coord.begin(3),
      snx.begin(), sny.begin(), snz.begin(), surf.begin(), vol.begin()
    );

    E_Float* cellN = unstrFc[v]->begin(posc);
    for (E_Int ind = 0; ind < nelts; ind++)
    {
      if (vol[ind] < 0.) cellN[ind] = 0.;
    }
  }
  //blanking des cellules dont leur propres facettes s intersectent
  vector<trifacette*> listFacettes;
  E_Int indA, indB, indC, indD, indE, indF;
  E_Int intersect = 0;
  E_Float ptA1[3]; E_Float ptB1[3]; E_Float ptC1[3];
  E_Float ptA2[3]; E_Float ptB2[3]; E_Float ptC2[3]; 
  E_Float xo, yo, zo;
  for (E_Int v = 0; v < nzones; v++)
  {
    FldArrayF* f = unstrF[v]; FldArrayI* cn = cnt[v];
    E_Float* xt = f->begin(posx); E_Float* yt = f->begin(posy); E_Float* zt = f->begin(posz);
    E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2); E_Int* cn3 = cn->begin(3); 
    E_Int* cn4 = cn->begin(4); E_Int* cn5 = cn->begin(5); E_Int* cn6 = cn->begin(6); 
    E_Int nelts = cn->getSize(); 
    vector< vector<E_Int> > cEEN(nelts);
    E_Float* cellN = unstrFc[v]->begin(posc);

    for (E_Int et = 0; et < nelts; et++)
    {
      indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1;
      indD = cn4[et]-1; indE = cn5[et]-1; indF = cn6[et]-1;
      //facette ABC
      trifacette* face = new trifacette;
      face->x1 = xt[indA]; face->y1 = yt[indA]; face->z1 = zt[indA]; 
      face->x2 = xt[indB]; face->y2 = yt[indB]; face->z2 = zt[indB]; 
      face->x3 = xt[indC]; face->y3 = yt[indC]; face->z3 = zt[indC]; 
      listFacettes.push_back(face);
      //facette DEF
      face = new trifacette;
      face->x1 = xt[indD]; face->y1 = yt[indD]; face->z1 = zt[indD]; 
      face->x2 = xt[indE]; face->y2 = yt[indE]; face->z2 = zt[indE]; 
      face->x3 = xt[indF]; face->y3 = yt[indF]; face->z3 = zt[indF]; 
      listFacettes.push_back(face);

      //facette ACFD decoupee en 4 triangles
      xo = 0.25*(xt[indA]+xt[indC]+xt[indF]+xt[indD]);
      yo = 0.25*(yt[indA]+yt[indC]+yt[indF]+yt[indD]);
      zo = 0.25*(zt[indA]+zt[indC]+zt[indF]+zt[indD]);
      face = new trifacette;//ACO
      face->x1 = xt[indA]; face->y1 = yt[indA]; face->z1 = zt[indA]; 
      face->x2 = xt[indC]; face->y2 = yt[indC]; face->z2 = zt[indC]; 
      face->x3 = xo; face->y3 = yo; face->z3 = zo; 
      listFacettes.push_back(face);
      face = new trifacette;//ADO
      face->x1 = xt[indA]; face->y1 = yt[indA]; face->z1 = zt[indA]; 
      face->x2 = xt[indD]; face->y2 = yt[indD]; face->z2 = zt[indD]; 
      face->x3 = xo; face->y3 = yo; face->z3 = zo; 
      listFacettes.push_back(face);
      face = new trifacette;//CFO
      face->x1 = xt[indF]; face->y1 = yt[indF]; face->z1 = zt[indF]; 
      face->x2 = xt[indC]; face->y2 = yt[indC]; face->z2 = zt[indC]; 
      face->x3 = xo; face->y3 = yo; face->z3 = zo; 
      listFacettes.push_back(face);
      face = new trifacette;//DFO
      face->x1 = xt[indF]; face->y1 = yt[indF]; face->z1 = zt[indF]; 
      face->x2 = xt[indD]; face->y2 = yt[indD]; face->z2 = zt[indD]; 
      face->x3 = xo; face->y3 = yo; face->z3 = zo; 
      listFacettes.push_back(face);

      //facette ABED decoupee en 4 triangles
      xo = 0.25*(xt[indA]+xt[indB]+xt[indE]+xt[indD]);
      yo = 0.25*(yt[indA]+yt[indB]+yt[indE]+yt[indD]);
      zo = 0.25*(zt[indA]+zt[indB]+zt[indE]+zt[indD]);
      face = new trifacette;//ABO
      face->x1 = xt[indA]; face->y1 = yt[indA]; face->z1 = zt[indA]; 
      face->x2 = xt[indB]; face->y2 = yt[indB]; face->z2 = zt[indB]; 
      face->x3 = xo; face->y3 = yo; face->z3 = zo; 
      listFacettes.push_back(face);
      face = new trifacette;//ADO
      face->x1 = xt[indA]; face->y1 = yt[indA]; face->z1 = zt[indA]; 
      face->x2 = xt[indD]; face->y2 = yt[indD]; face->z2 = zt[indD]; 
      face->x3 = xo; face->y3 = yo; face->z3 = zo; 
      listFacettes.push_back(face);
      face = new trifacette;//BEO
      face->x1 = xt[indB]; face->y1 = yt[indB]; face->z1 = zt[indB]; 
      face->x2 = xt[indE]; face->y2 = yt[indE]; face->z2 = zt[indE]; 
      face->x3 = xo; face->y3 = yo; face->z3 = zo; 
      listFacettes.push_back(face);
      face = new trifacette;//EDO
      face->x1 = xt[indE]; face->y1 = yt[indE]; face->z1 = zt[indE]; 
      face->x2 = xt[indD]; face->y2 = yt[indD]; face->z2 = zt[indD]; 
      face->x3 = xo; face->y3 = yo; face->z3 = zo; 
      listFacettes.push_back(face);

      //facette BCFE decoupee en 4 triangles
      xo = 0.25*(xt[indB]+xt[indC]+xt[indF]+xt[indE]);
      yo = 0.25*(yt[indB]+yt[indC]+yt[indF]+yt[indE]);
      zo = 0.25*(zt[indB]+zt[indC]+zt[indF]+zt[indE]);
      face = new trifacette;//BCO
      face->x1 = xt[indC]; face->y1 = yt[indC]; face->z1 = zt[indC]; 
      face->x2 = xt[indB]; face->y2 = yt[indB]; face->z2 = zt[indB]; 
      face->x3 = xo; face->y3 = yo; face->z3 = zo; 
      listFacettes.push_back(face);
      face = new trifacette;//CFO
      face->x1 = xt[indC]; face->y1 = yt[indC]; face->z1 = zt[indC]; 
      face->x2 = xt[indF]; face->y2 = yt[indF]; face->z2 = zt[indF]; 
      face->x3 = xo; face->y3 = yo; face->z3 = zo; 
      listFacettes.push_back(face);
      face = new trifacette;//EFO
      face->x1 = xt[indF]; face->y1 = yt[indF]; face->z1 = zt[indF]; 
      face->x2 = xt[indE]; face->y2 = yt[indE]; face->z2 = zt[indE]; 
      face->x3 = xo; face->y3 = yo; face->z3 = zo; 
      listFacettes.push_back(face);
      face = new trifacette;//BEO
      face->x1 = xt[indE]; face->y1 = yt[indE]; face->z1 = zt[indE]; 
      face->x2 = xt[indB]; face->y2 = yt[indB]; face->z2 = zt[indB]; 
      face->x3 = xo; face->y3 = yo; face->z3 = zo; 
      listFacettes.push_back(face);

      E_Int nfacettes = listFacettes.size();
      for (E_Int v1 = 0; v1 < nfacettes; v1++)
      {
        trifacette* face1 = listFacettes[v1];
        ptA1[0] = face1->x1; ptA1[1] = face1->y1; ptA1[2] = face1->z1;
        ptB1[0] = face1->x2; ptB1[1] = face1->y2; ptB1[2] = face1->z2;
        ptC1[0] = face1->x3; ptC1[1] = face1->y3; ptC1[2] = face1->z3;
       
        for (E_Int v2 = v1+1; v2 < nfacettes; v2++)
        {
          trifacette* face2 = listFacettes[v2];
          ptA2[0] = face2->x1; ptA2[1] = face2->y1; ptA2[2] = face2->z1;
          ptB2[0] = face2->x2; ptB2[1] = face2->y2; ptB2[2] = face2->z2;
          ptC2[0] = face2->x3; ptC2[1] = face2->y3; ptC2[2] = face2->z3;
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptC1, ptA2, ptB2, ptC2, eps);
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}
        }
      }
      end:;
      for (E_Int v0 = 0; v0 < nfacettes; v0++) delete listFacettes[v0];
      listFacettes.clear();
    }
  }
  return 1;
}
//=============================================================================
/* blanking cells with self-intersecting faces and negative volume */
//=============================================================================
E_Int K_CONNECTOR::blankInvalidCellsHexa( 
  E_Float eps, E_Int posx, E_Int posy, E_Int posz, E_Int posc,
  vector<FldArrayI*>& cnt, vector<FldArrayF*>& unstrF, 
  vector<FldArrayI*>& cnct, vector<FldArrayF*>& unstrFc)
{
  E_Int nzones = unstrF.size();
  // blanking des cellules de volume negatif
  for (E_Int v = 0; v < nzones; v++)
  {
    E_Int npts = unstrF[v]->getSize();
    FldArrayF coord(npts, 3);
    coord.setOneField(*unstrF[v], posx, 1);
    coord.setOneField(*unstrF[v], posy, 2);
    coord.setOneField(*unstrF[v], posz, 3);

    FldArrayI& cm = *(cnt[v]->getConnect(0));
    E_Int nelts = cm.getSize();
    E_Int nfacets = 6*nelts;
    FldArrayF snx(nfacets), sny(nfacets), snz(nfacets), surf(nfacets);
    FldArrayF vol(nelts);
    K_METRIC::compMetricUnstruct(
      *cnt[v], "HEXA", coord.begin(1), coord.begin(2), coord.begin(3),
      snx.begin(), sny.begin(), snz.begin(), surf.begin(), vol.begin()
    );

    E_Float* cellN = unstrFc[v]->begin(posc);
    for (E_Int ind = 0; ind < nelts; ind++)
    {
      if (vol[ind] < 0.) cellN[ind] = 0.;
    }
  }
  //blanking des cellules dont leur propres facettes s intersectent
  vector<facette*> listFacettes;
  E_Int indA, indB, indC, indD, indE, indF, indG, indH;
  E_Int intersect = 0;
  E_Float ptA1[3]; E_Float ptB1[3]; E_Float ptC1[3]; E_Float ptD1[3]; E_Float ptO1[3];
  E_Float ptA2[3]; E_Float ptB2[3]; E_Float ptC2[3]; E_Float ptD2[3]; E_Float ptO2[3];
  for (E_Int v = 0; v < nzones; v++)
  {
    FldArrayF* f = unstrF[v]; FldArrayI* cn = cnt[v];
    E_Float* xt = f->begin(posx); E_Float* yt = f->begin(posy); E_Float* zt = f->begin(posz);
    E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2); E_Int* cn3 = cn->begin(3); E_Int* cn4 = cn->begin(4); 
    E_Int* cn5 = cn->begin(5); E_Int* cn6 = cn->begin(6); E_Int* cn7 = cn->begin(7); E_Int* cn8 = cn->begin(8); 
    E_Int nelts = cn->getSize();
    E_Float* cellN = unstrFc[v]->begin(posc);

    for (E_Int et = 0; et < nelts; et++)
    {
      indA = cn1[et]-1; indB = cn2[et]-1; indC = cn3[et]-1; indD = cn4[et]-1;
      indE = cn5[et]-1; indF = cn6[et]-1; indG = cn7[et]-1; indH = cn8[et]-1;
      facette* face = new facette;  // facette ABCD
      face->indA = indA; face->indB = indB; face->indC = indC; face->indD = indD; 
      face->xo = 0.25*(xt[indA]+xt[indB]+xt[indC]+xt[indD]);
      face->yo = 0.25*(yt[indA]+yt[indB]+yt[indC]+yt[indD]);
      face->zo = 0.25*(zt[indA]+zt[indB]+zt[indC]+zt[indD]);
      face->nozone = v; listFacettes.push_back(face);

      face = new facette;//facette EFGH
      face->indA = indE; face->indB = indF; face->indC = indG; face->indD = indH; 
      face->xo = 0.25*(xt[indE]+xt[indF]+xt[indG]+xt[indH]);
      face->yo = 0.25*(yt[indE]+yt[indF]+yt[indG]+yt[indH]);
      face->zo = 0.25*(zt[indE]+zt[indF]+zt[indG]+zt[indH]);
      face->nozone = v; listFacettes.push_back(face);

      face = new facette; //face ABFE
      face->indA = indA; face->indB = indB; face->indC = indF; face->indD = indE; 
      face->xo = 0.25*(xt[indA]+xt[indB]+xt[indF]+xt[indE]);
      face->yo = 0.25*(yt[indA]+yt[indB]+yt[indF]+yt[indE]);
      face->zo = 0.25*(zt[indA]+zt[indB]+zt[indF]+zt[indE]);
      face->nozone = v; listFacettes.push_back(face);

      face = new facette; //face DCGH
      face->indA = indD; face->indB = indC; face->indC = indG; face->indD = indH; 
      face->xo = 0.25*(xt[indD]+xt[indC]+xt[indG]+xt[indH]);
      face->yo = 0.25*(yt[indD]+yt[indC]+yt[indG]+yt[indH]);
      face->zo = 0.25*(zt[indD]+zt[indC]+zt[indG]+zt[indH]);
      face->nozone = v; listFacettes.push_back(face);
  
      face = new facette; //face ADHE
      face->indA = indA; face->indB = indD; face->indC = indH; face->indD = indE; 
      face->xo = 0.25*(xt[indA]+xt[indD]+xt[indH]+xt[indE]);
      face->yo = 0.25*(yt[indA]+yt[indD]+yt[indH]+yt[indE]);
      face->zo = 0.25*(zt[indA]+zt[indD]+zt[indH]+zt[indE]);
      face->nozone = v; listFacettes.push_back(face);
  
      face = new facette; //face BCGF
      face->indA = indB; face->indB = indC; face->indC = indG; face->indD = indF; 
      face->xo = 0.25*(xt[indB]+xt[indC]+xt[indG]+xt[indF]);
      face->yo = 0.25*(yt[indB]+yt[indC]+yt[indG]+yt[indF]);
      face->zo = 0.25*(zt[indB]+zt[indC]+zt[indG]+zt[indF]);
      face->nozone = v; listFacettes.push_back(face);
   
      E_Int nfacettes = listFacettes.size();
      for (E_Int v1 = 0; v1 < nfacettes; v1++)
      {
        facette* face1 = listFacettes[v1];
        ptA1[0] = xt[face1->indA]; ptA1[1] = yt[face1->indA]; ptA1[2] = zt[face1->indA];
        ptB1[0] = xt[face1->indB]; ptB1[1] = yt[face1->indB]; ptB1[2] = zt[face1->indB];
        ptC1[0] = xt[face1->indC]; ptC1[1] = yt[face1->indC]; ptC1[2] = zt[face1->indC];
        ptD1[0] = xt[face1->indD]; ptD1[1] = yt[face1->indD]; ptD1[2] = zt[face1->indD];
        ptO1[0] = face1->xo; ptO1[1] = face1->yo; ptO1[2] = face1->zo;     
        for (E_Int v2 = v1+1; v2 < nfacettes; v2++)
        {
          facette* face2 = listFacettes[v2];
          ptA2[0] = xt[face2->indA]; ptA2[1] = yt[face2->indA]; ptA2[2] = zt[face2->indA];
          ptB2[0] = xt[face2->indB]; ptB2[1] = yt[face2->indB]; ptB2[2] = zt[face2->indB];
          ptC2[0] = xt[face2->indC]; ptC2[1] = yt[face2->indC]; ptC2[2] = zt[face2->indC];
          ptD2[0] = xt[face2->indD]; ptD2[1] = yt[face2->indD]; ptD2[2] = zt[face2->indD];
          ptO2[0] = face2->xo; ptO2[1] = face2->yo; ptO2[2] = face2->zo;
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}

          //B2C2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptB2, ptC2, ptO2, eps);
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptB2, ptC2, ptO2, eps);        
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptB2, ptC2, ptO2, eps);
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptB2, ptC2, ptO2, eps);
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}
          //C2D2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}
          //A2D2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) {cellN[et] = 0.; goto end;}
        }
      }
      end:;
      for (E_Int v0 = 0; v0 < nfacettes; v0++) delete listFacettes[v0];
      listFacettes.clear();
    }
  }
  return 1;
}
//=============================================================================
/* blanking cells with self-intersecting faces and negative volume */
//=============================================================================
E_Int K_CONNECTOR::blankInvalidCellsStruct( 
  E_Float eps, E_Int posx, E_Int posy, E_Int posz, E_Int posc,
  vector<E_Int>& nit, vector<E_Int>& njt, vector<E_Int>& nkt, vector<FldArrayF*>& structF, 
  vector<E_Int>& nict, vector<E_Int>& njct, vector<E_Int>& nkct, vector<FldArrayF*>& structFc)
{
  E_Int nzones = structF.size();

  /* blanking des cellules de volume negatif */
  for (E_Int v = 0; v < nzones; v++)
  {
    E_Int ni = nit[v]; E_Int nj = njt[v]; E_Int nk = nkt[v]; 
    E_Int nic = nict[v]; E_Int njc = njct[v]; E_Int nkc = nkct[v];
    E_Int ncells = nic*njc*nkc;
    E_Int ninti = ni*njc*nkc; E_Int nintj = nic*nj*nkc; E_Int nintk = nic*njc*nk;
    E_Int nint =  ninti + nintj + nintk;
    FldArrayF vol(ncells); FldArrayF surf(nint,3);
    FldArrayF snorm(nint); FldArrayF centerInt(nint, 3);
    FldArrayF coord(ni*nj*nk, 3);
    coord.setOneField(*structF[v], posx, 1);
    coord.setOneField(*structF[v], posy, 2);
    coord.setOneField(*structF[v], posz, 3);
    K_METRIC::compMetricStruct(
      ni, nj, nk, ninti, nintj, nintk,
      coord.begin(1), coord.begin(2), coord.begin(3), 
      vol.begin(), 
      surf.begin(1), surf.begin(2), surf.begin(3), 
      snorm.begin(), 
      centerInt.begin(1), centerInt.begin(2), centerInt.begin(3));
    E_Float* cellN = structFc[v]->begin(posc);
    for (E_Int ind = 0; ind < ncells; ind++)
    {
      if ( vol[ind] < 0. ) cellN[ind] = 0.;
    }
  }
  //blanking des cellules dont leur propres facettes s intersectent
  vector<facette*> listFacettes;
  E_Int indA, indB, indC, indD, indE, indF, indG, indH;
  E_Int intersect = 0;
  E_Float ptA1[3]; E_Float ptB1[3]; E_Float ptC1[3]; E_Float ptD1[3]; E_Float ptO1[3];
  E_Float ptA2[3]; E_Float ptB2[3]; E_Float ptC2[3]; E_Float ptD2[3]; E_Float ptO2[3];
  for (E_Int v = 0; v < nzones; v++)
  {
    FldArrayF* f = structF[v];
    E_Int ni = nit[v]; E_Int nj = njt[v]; E_Int ninj = ni*nj; 
    E_Int nic = nict[v]; E_Int njc = njct[v]; E_Int nicnjc = nic*njc;
    E_Float* xt = f->begin(posx); E_Float* yt = f->begin(posy); E_Float* zt = f->begin(posz);
    E_Int ncells = structFc[v]->getSize();
    E_Float* cellN = structFc[v]->begin(posc);
    for (E_Int ind = 0; ind < ncells; ind++)
    {
      E_Int k = ind/nicnjc; E_Int j = (ind-k*nicnjc)/nic; E_Int i = ind-j*nic-k*nicnjc;
      indA = i+j*ni+k*ninj; indB = indA+1; indC = indB+ni; indD = indA+ni;
      indE = indA+ninj; indF = indB+ninj; indG = indC+ninj; indH = indD+ninj;
      facette* face = new facette; //face ABCD
      face->indA = indA; face->indB = indB; face->indC = indC; face->indD = indD; 
      face->xo = 0.25*(xt[indA]+xt[indB]+xt[indC]+xt[indD]);
      face->yo = 0.25*(yt[indA]+yt[indB]+yt[indC]+yt[indD]);
      face->zo = 0.25*(zt[indA]+zt[indB]+zt[indC]+zt[indD]);
      face->indcell1 = ind; face->indcell2 = -1;
      face->nozone = v; listFacettes.push_back(face);
      face = new facette; //face EFGH
      face->indA = indE; face->indB = indF; face->indC = indG; face->indD = indH; 
      face->xo = 0.25*(xt[indE]+xt[indF]+xt[indG]+xt[indH]);
      face->yo = 0.25*(yt[indE]+yt[indF]+yt[indG]+yt[indH]);
      face->zo = 0.25*(zt[indE]+zt[indF]+zt[indG]+zt[indH]);
      face->indcell1 = ind; face->indcell2 = -1;
      face->nozone = v; listFacettes.push_back(face);
      face = new facette; //face ABFE
      face->indA = indA; face->indB = indB; face->indC = indF; face->indD = indE; 
      face->xo = 0.25*(xt[indA]+xt[indB]+xt[indF]+xt[indE]);
      face->yo = 0.25*(yt[indA]+yt[indB]+yt[indF]+yt[indE]);
      face->zo = 0.25*(zt[indA]+zt[indB]+zt[indF]+zt[indE]);
      face->indcell1 = ind; face->indcell2 = -1;
      face->nozone = v; listFacettes.push_back(face);      
      face = new facette; //face DCGH
      face->indA = indD; face->indB = indC; face->indC = indG; face->indD = indH; 
      face->xo = 0.25*(xt[indD]+xt[indC]+xt[indG]+xt[indH]);
      face->yo = 0.25*(yt[indD]+yt[indC]+yt[indG]+yt[indH]);
      face->zo = 0.25*(zt[indD]+zt[indC]+zt[indG]+zt[indH]);
      face->indcell1 = ind; face->indcell2 = -1;
      face->nozone = v; listFacettes.push_back(face);       
      face = new facette; //face ADHE
      face->indA = indA; face->indB = indD; face->indC = indH; face->indD = indE; 
      face->xo = 0.25*(xt[indA]+xt[indD]+xt[indH]+xt[indE]);
      face->yo = 0.25*(yt[indA]+yt[indD]+yt[indH]+yt[indE]);
      face->zo = 0.25*(zt[indA]+zt[indD]+zt[indH]+zt[indE]);
      face->indcell1 = ind; face->indcell2 = -1;
      face->nozone = v; listFacettes.push_back(face);    
      face = new facette; //face BCGF
      face->indA = indB; face->indB = indC; face->indC = indG; face->indD = indF; 
      face->xo = 0.25*(xt[indB]+xt[indC]+xt[indG]+xt[indF]);
      face->yo = 0.25*(yt[indB]+yt[indC]+yt[indG]+yt[indF]);
      face->zo = 0.25*(zt[indB]+zt[indC]+zt[indG]+zt[indF]);
      face->indcell1 = ind; face->indcell2 = -1;
      face->nozone = v; listFacettes.push_back(face);    
      E_Int nfacettes = listFacettes.size();
      for (E_Int v1 = 0; v1 < nfacettes; v1++)
      {
        facette* face1 = listFacettes[v1];
        ptA1[0] = xt[face1->indA]; ptA1[1] = yt[face1->indA]; ptA1[2] = zt[face1->indA];
        ptB1[0] = xt[face1->indB]; ptB1[1] = yt[face1->indB]; ptB1[2] = zt[face1->indB];
        ptC1[0] = xt[face1->indC]; ptC1[1] = yt[face1->indC]; ptC1[2] = zt[face1->indC];
        ptD1[0] = xt[face1->indD]; ptD1[1] = yt[face1->indD]; ptD1[2] = zt[face1->indD];
        ptO1[0] = face1->xo; ptO1[1] = face1->yo; ptO1[2] = face1->zo;     
        for (E_Int v2 = v1+1; v2 < nfacettes; v2++)
        {
          facette* face2 = listFacettes[v2];
          ptA2[0] = xt[face2->indA]; ptA2[1] = yt[face2->indA]; ptA2[2] = zt[face2->indA];
          ptB2[0] = xt[face2->indB]; ptB2[1] = yt[face2->indB]; ptB2[2] = zt[face2->indB];
          ptC2[0] = xt[face2->indC]; ptC2[1] = yt[face2->indC]; ptC2[2] = zt[face2->indC];
          ptD2[0] = xt[face2->indD]; ptD2[1] = yt[face2->indD]; ptD2[2] = zt[face2->indD];
          ptO2[0] = face2->xo; ptO2[1] = face2->yo; ptO2[2] = face2->zo;
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) {cellN[ind] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) {cellN[ind] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) {cellN[ind] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptA2, ptB2, ptO2, eps);
          if ( intersect < 0 ) {cellN[ind] = 0.; goto end;}

          //B2C2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptB2, ptC2, ptO2, eps);
          if ( intersect < 0 ) {cellN[ind] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptB2, ptC2, ptO2, eps);        
          if ( intersect < 0 ) {cellN[ind] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptB2, ptC2, ptO2, eps);
          if ( intersect < 0 ) {cellN[ind] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptB2, ptC2, ptO2, eps);
          if ( intersect < 0 ) {cellN[ind] = 0.; goto end;}
          //C2D2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) {cellN[ind] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) {cellN[ind] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) {cellN[ind] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptC2, ptD2, ptO2, eps);
          if ( intersect < 0 ) {cellN[ind] = 0.; goto end;}
          //A2D2O2
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptA1, ptB1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) {cellN[ind] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptB1, ptC1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) {cellN[ind] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptC1, ptD1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) {cellN[ind] = 0.; goto end;}
          intersect = K_COMPGEOM::crossIntersectionOfTriangles(ptD1, ptA1, ptO1, ptA2, ptD2, ptO2, eps);
          if ( intersect < 0 ) {cellN[ind] = 0.; goto end;}
        }
      }
      end:;
      for (E_Int v0 = 0; v0 < nfacettes; v0++) delete listFacettes[v0];
      listFacettes.clear();
    }
  }
  return 1;
}
//=============================================================================
/* For HEXA mesh generated by extrusion, determines the number of nk layers in
   the normal direction to the initial surface 
   IN : npts : nb de noeuds dans le maillage
   IN : cEV : connectivite elts->noeuds
   retourne 0 si pb */
//=============================================================================
E_Int K_CONNECTOR::getNbOfHexaLayers(E_Int npts, FldArrayI& cEV)
{
  E_Int nelts = cEV.getSize();
  vector< vector<E_Int> > cEEN(nelts);
  E_Int ok = K_CONNECT::connectEV2EENbrs("HEXA", npts, cEV, cEEN); 
  if ( ok == 0) return 0;
  E_Int ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;
  E_Int* cn1 = cEV.begin(1); E_Int* cn2 = cEV.begin(2);
  E_Int* cn3 = cEV.begin(3); E_Int* cn4 = cEV.begin(4);
  E_Int* cn5 = cEV.begin(5); E_Int* cn6 = cEV.begin(6);
  E_Int* cn7 = cEV.begin(7); E_Int* cn8 = cEV.begin(8);
  E_Int nfaces = 0;
  E_Int found = 0;
  for (E_Int et = 0; et < nelts; et++)
  {
    ind5 = cn5[et]; ind6 = cn6[et]; ind7 = cn7[et]; ind8 = cn8[et];
    vector<E_Int>& voisins = cEEN[et]; E_Int nvoisins = voisins.size();
    found = 0;
    for (E_Int noet2 = 0; noet2 < nvoisins; noet2++)
    {
      E_Int et2 = voisins[noet2];
      ind1 = cn1[et2]; ind2 = cn2[et2]; ind3 = cn3[et2]; ind4 = cn4[et2]; 
      if ( ind1 == ind5 && ind2 == ind6 && ind3 == ind7 && ind4 == ind8 ) found = 1;
    }
    if ( found == 0 ) nfaces++;
  }
  if ( nfaces == 0 ) return 0;
  E_Int nk = nelts/nfaces+1;
  return nk;

}
//=============================================================================
/* For PENTA mesh generated by extrusion, determines the number of nk layers in
   the normal direction to the initial surface 
   IN : npts : nb de noeuds dans le maillage
   IN : cEV : connectivite elts->noeuds
   retourne 0 si pb */
//=============================================================================
E_Int K_CONNECTOR::getNbOfPentaLayers(E_Int npts, FldArrayI& cEV)
{
  E_Int nelts = cEV.getSize();
  vector< vector<E_Int> > cEEN(nelts);
  E_Int ok = K_CONNECT::connectEV2EENbrs("PENTA", npts, cEV, cEEN); 
  if ( ok == 0 ) return 0;
  E_Int ind1, ind2, ind3, ind4, ind5, ind6;
  E_Int* cn1 = cEV.begin(1); E_Int* cn2 = cEV.begin(2);
  E_Int* cn3 = cEV.begin(3); E_Int* cn4 = cEV.begin(4);
  E_Int* cn5 = cEV.begin(5); E_Int* cn6 = cEV.begin(6);
  E_Int nfaces = 0; E_Int found = 0;
  for (E_Int et = 0; et < nelts; et++)
  {
    ind4 = cn4[et]; ind5 = cn5[et]; ind6 = cn6[et];
    vector<E_Int>& voisins = cEEN[et]; E_Int nvoisins = voisins.size();
    found = 0;
    for (E_Int noet2 = 0; noet2 < nvoisins; noet2++)
    {
      E_Int et2 = voisins[noet2];
      ind1 = cn1[et2]; ind2 = cn2[et2]; ind3 = cn3[et2];
      if ( ind1 == ind4 && ind2 == ind5 && ind3 == ind6 ) found = 1;
    }
    if ( found == 0 ) nfaces++;
  }
  if ( nfaces == 0 ) return 0;
  E_Int nk = nelts/nfaces+1;
  return nk;

}
