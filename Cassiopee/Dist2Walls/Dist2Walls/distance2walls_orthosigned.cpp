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
# include "Nuga/include/KdTree.h"
# include "Nuga/include/BbTree.h"
# include "Nuga/include/ArrayAccessor.h"

using namespace std;
using namespace K_FLD; 
using namespace K_SEARCH;

//=============================================================================
/* Calcul du distancefield signe par projection orthogonale
   Les corps doivent etre en TRI, avec le cellN localise aux noeuds du 
   maillage TRI */
//=============================================================================
PyObject* K_DIST2WALLS::distance2WallsOrthoSigned(PyObject* self, 
                                                  PyObject* args)
{
  PyObject *centers, *bodiesC;
  if (!PyArg_ParseTuple(args, "OO", &centers, &bodiesC)) return NULL;
  
  if (PyList_Check(centers) == 0)
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

  // Maillage en noeuds
  // les coordonnees doivent etre en premier
  vector<PyObject*> objsn, objun;
  vector<char*> structVarStringn; vector<char*> unstrVarStringn;
  vector<FldArrayF*> structFn; vector<FldArrayF*> unstrFn;
  vector<E_Int> nitn; vector<E_Int> njtn; vector<E_Int> nktn;
  vector<FldArrayI*> cntn;
  vector<char*> eltTypen;
  char* varStringl; char* eltTypel;
  E_Int nil, njl, nkl;
  FldArrayF* fl; FldArrayI* cnl;
  PyObject* o; E_Int res;
  E_Int nz = PyList_Size(centers);
  for (E_Int i = 0; i < nz; i++)
  {
    o = PyList_GetItem(centers, i);
    res = K_ARRAY::getFromArray3(o, varStringl, fl, 
                                 nil, njl, nkl, cnl, eltTypel);
    if (res == 1)
    {
      structVarStringn.push_back(varStringl);
      structFn.push_back(fl);
      nitn.push_back(nil); njtn.push_back(njl), nktn.push_back(nkl);
      objsn.push_back(o);
    }
    else if (res == 2)
    {
      unstrVarStringn.push_back(varStringl);
      unstrFn.push_back(fl);
      eltTypen.push_back(eltTypel); cntn.push_back(cnl);
      objun.push_back(o);
    }
    else printf("Warning: dist2Walls: array " SF_D_ " is invalid. Discarded.\n", i);
  }

  E_Int nsn = structFn.size(); E_Int nun = unstrFn.size();
  E_Int posx=-1, posy=-1, posz=-1;
  vector<E_Int> ncellss; vector<E_Int> ncellsu;
  // Verification de posxi, posyi, poszi dans listFields: vars 1,2,3 imposees
  for (E_Int i = 0; i < nsn; i++)
  {
    E_Int posxi = K_ARRAY::isCoordinateXPresent(structVarStringn[i]);
    E_Int posyi = K_ARRAY::isCoordinateYPresent(structVarStringn[i]);
    E_Int poszi = K_ARRAY::isCoordinateZPresent(structVarStringn[i]);
    posxi++; posyi++; poszi++; 
    if (posx == -1) posx = posxi;
    if (posy == -1) posy = posyi;
    if (posz == -1) posz = poszi;
    if (posxi != posx || posyi != posy || poszi != posz) 
    {
      PyErr_SetString(PyExc_TypeError,
                      "dist2Walls: coordinates must be located at same position for all zones.");
      for (E_Int nos = 0; nos < nsn; nos++)
        RELEASESHAREDS(objsn[nos], structFn[nos]);
      for (E_Int nos = 0; nos < nun; nos++)
        RELEASESHAREDU(objun[nos], unstrFn[nos], cntn[nos]);
      return NULL;
    }
    ncellss.push_back(structFn[i]->getSize());
  }
  for (E_Int i = 0; i < nun; i++)
  {
    E_Int posxi = K_ARRAY::isCoordinateXPresent(unstrVarStringn[i]);
    E_Int posyi = K_ARRAY::isCoordinateYPresent(unstrVarStringn[i]);
    E_Int poszi = K_ARRAY::isCoordinateZPresent(unstrVarStringn[i]);
    posxi++; posyi++; poszi++;
    if (posx == -1) posx = posxi;
    if (posy == -1) posy = posyi;
    if (posz == -1) posz = poszi;

    if (posxi != posx || posyi != posy || poszi != posz) 
    {
      PyErr_SetString( PyExc_TypeError,
        "dist2Walls: coordinates must be located at same position for all zones.");
      for (E_Int nos = 0; nos < nsn; nos++)
        RELEASESHAREDS(objsn[nos], structFn[nos]);
      for (E_Int nos = 0; nos < nun; nos++)
        RELEASESHAREDU(objun[nos], unstrFn[nos], cntn[nos]);
      return NULL;
    }
    ncellsu.push_back(unstrFn[i]->getSize());
  }
  // Extract infos from body surfaces + sx, sy, sz
  vector<E_Int> resl;
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
    PyErr_SetString(PyExc_TypeError,"distance2Walls: invalid list of surfaces.");
    for (E_Int nos = 0; nos < nsn; nos++)
      RELEASESHAREDS(objsn[nos], structFn[nos]);
    for (E_Int nos = 0; nos < nun; nos++)
      RELEASESHAREDU(objun[nos], unstrFn[nos], cntn[nos]);
    return NULL;
  }
  
  // Get posxv, posyv, poszv in bodiesC, already checked  
  vector<E_Int> posxv;  vector<E_Int> posyv;  vector<E_Int> poszv;
  vector<E_Int> poscv;
  E_Int possx = K_ARRAY::isNamePresent("sx", unstrVarString[0]); possx++;
  E_Int possy = K_ARRAY::isNamePresent("sy", unstrVarString[0]); possy++;
  E_Int possz = K_ARRAY::isNamePresent("sz", unstrVarString[0]); possz++;
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
      PyErr_SetString(PyExc_TypeError,"distance2Walls: cellN must be defined for bodies.");
      for (E_Int nos = 0; nos < nsn; nos++)
        RELEASESHAREDS(objsn[nos], structFn[nos]);
      for (E_Int nos = 0; nos < nun; nos++)
        RELEASESHAREDU(objun[nos], unstrFn[nos], cntn[nos]);
      for (E_Int nos = 0; nos < nwalls; nos++)
        RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
      return NULL;
    }
    if (K_ARRAY::isNamePresent("sx", unstrVarString[v])+1 != possx ||
        K_ARRAY::isNamePresent("sy", unstrVarString[v])+1 != possy ||
        K_ARRAY::isNamePresent("sz", unstrVarString[v])+1 != possz)
    {
      PyErr_SetString(PyExc_TypeError,"distance2Walls: normals vectors must have the same localization.");
      for (E_Int nos = 0; nos < nsn; nos++)
        RELEASESHAREDS(objsn[nos], structFn[nos]);
      for (E_Int nos = 0; nos < nun; nos++)
        RELEASESHAREDU(objun[nos], unstrFn[nos], cntn[nos]);
      for (E_Int nos = 0; nos < nwalls; nos++)
        RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);
      return NULL;
    }
  }

  // Initialisation de la distance pour les maillages structures
  vector<FldArrayF*> distances;
  for (E_Int i = 0; i < nsn; i++)
  {
    E_Int ncells = ncellss[i];
    FldArrayF* distance = new FldArrayF(ncells); 
    distance->setAllValuesAt(K_CONST::E_INFINITE);
    distances.push_back(distance);
  }
  // Initialisation de la distance pour les maillages non structures
  vector<FldArrayF*> distancesu;
  for (E_Int i = 0; i < nun; i++)
  {
    E_Int ncells = ncellsu[i];
    FldArrayF* distance = new FldArrayF(ncells); 
    distance->setAllValuesAt(K_CONST::E_INFINITE);
    distancesu.push_back(distance);
  }
  if (structFn.size() > 0)
    computeSignedOrthoDist(ncellss, posx, posy, posz, 
                           structFn, posxv, posyv, poszv, poscv, unstrF, cnt,
                           possx, possy, possz, 
                           distances);
  if (unstrFn.size() > 0)
    computeSignedOrthoDist(ncellsu, posx, posy, posz,
                           unstrFn, posxv, posyv, poszv, poscv, unstrF, cnt, 
                           possx, possy, possz, 
                           distancesu);
  
  for (E_Int nos = 0; nos < nwalls; nos++)
    RELEASESHAREDU(obju[nos], unstrF[nos], cnt[nos]);

  // Build arrays
  PyObject* l = PyList_New(0);
  PyObject* tpl;    
  for (E_Int nos = 0; nos < nsn; nos++)
  {
    tpl = K_ARRAY::buildArray(*distances[nos], "TurbulentDistance", 
                              nitn[nos], njtn[nos], nktn[nos]);
    PyList_Append(l, tpl); Py_DECREF(tpl);
    delete distances[nos];
    RELEASESHAREDS(objsn[nos], structFn[nos]);
  }
  for (E_Int nou = 0; nou < nun; nou++)
  {
    K_FLD::FldArrayI* cnout = new K_FLD::FldArrayI(*cntn[nou]);
    tpl = K_ARRAY::buildArray(*distancesu[nou], "TurbulentDistance", *cnout,
                              -1, eltTypen[nou]);
    PyList_Append(l, tpl); Py_DECREF(tpl);
    RELEASESHAREDU(objun[nou], unstrFn[nou], cntn[nou]);
    delete distancesu[nou]; delete cnout;
  }
  return l;
}

//=============================================================================
/* Calcul de la distance a la paroi */
//=============================================================================
void K_DIST2WALLS::computeSignedOrthoDist(
  vector<E_Int>& ncellst,
  E_Int posx, E_Int posy, E_Int posz, 
  vector<FldArrayF*>& fields, 
  vector<E_Int>& posxv, vector<E_Int>& posyv, vector<E_Int>& poszv, 
  vector<E_Int>& poscv, vector<FldArrayF*>& fieldsw, vector<FldArrayI*>& cntw,
  E_Int possx, E_Int possy, E_Int possz, 
  vector<FldArrayF*>& distances)
{
  E_Int nzones = fields.size();
  /* 1 - creation du kdtree et du bbtree */
  typedef K_SEARCH::BoundingBox<3>  BBox3DType; 
  E_Float minB[3];  E_Float maxB[3];
  vector< vector<BBox3DType*> > vectOfBoxes;// a detruire a la fin
  
  // allocate kdtree array : kdtree points are cell vertices
  E_Int nwalls = cntw.size(); E_Int nptsmax = 0; 
  for (E_Int v = 0; v < nwalls; v++) nptsmax += fieldsw[v]->getSize();
  FldArrayF* wallpts = new FldArrayF(nptsmax,6);
  FldArrayF* lmax = new FldArrayF(nptsmax);
  E_Float* xw2 = wallpts->begin(1);
  E_Float* yw2 = wallpts->begin(2);
  E_Float* zw2 = wallpts->begin(3);
  E_Float* lmaxp = lmax->begin();
  E_Float* sxw2 = wallpts->begin(4);
  E_Float* syw2 = wallpts->begin(5);
  E_Float* szw2 = wallpts->begin(6);
  // create kdtree elements and bbtree information
  E_Int nop = 0; 
  vector<FldArrayF> bboxes;
  vector<E_Int> indirW; vector<E_Int> indirZ;
  for (E_Int now = 0; now < nwalls; now++)
  {
    FldArrayF* fieldv = fieldsw[now];
    E_Int posxw = posxv[now]; E_Int posyw = posyv[now]; E_Int poszw = poszv[now];  
    E_Float* xw = fieldv->begin(posxw);
    E_Float* yw = fieldv->begin(posyw);
    E_Float* zw = fieldv->begin(poszw);
    E_Float* sxw = fieldv->begin(possx);
    E_Float* syw = fieldv->begin(possy);
    E_Float* szw = fieldv->begin(possz);
    E_Int poscw = poscv[now]; E_Float* cellnw = fieldv->begin(poscw);
    E_Int npts = fieldv->getSize(); 
    FldArrayI& cnloc = *cntw[now];
    E_Int nelts = cnloc.getSize(); E_Int nvert = cnloc.getNfld();
    vector<BBox3DType*> boxes(nelts);// liste des bbox de ts les elements de la paroi courante
    FldArrayF bbox(nelts,6);// xmin, ymin, zmin, xmax, ymax, zmax
    K_COMPGEOM::boundingBoxOfUnstrCells(cnloc, xw, yw, zw, bbox); bboxes.push_back(bbox);
    E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
    E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
    E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
    FldArrayIS dejavu(npts); dejavu.setAllValuesAtNull();
    short* dejavup = dejavu.begin();
    for (E_Int et = 0; et < nelts; et++)
    {
      minB[0] = xminp[et]; minB[1] = yminp[et]; minB[2] = zminp[et];
      maxB[0] = xmaxp[et]; maxB[1] = ymaxp[et]; maxB[2] = zmaxp[et]; 
      boxes[et] = new BBox3DType(minB, maxB);
      for (E_Int nov = 1; nov <= nvert; nov++)
      {
        E_Int ind = cnloc(et,nov)-1;
        if (cellnw[ind] == 1. && dejavup[ind] == 0)
        {
          xw2[nop] = xw[ind]; yw2[nop] = yw[ind]; zw2[nop] = zw[ind]; 
          sxw2[nop] = sxw[ind]; syw2[nop] = syw[ind]; szw2[nop] = szw[ind];
          lmaxp[nop] = K_FUNC::E_max(maxB[0]-minB[0], maxB[1]-minB[1], maxB[2]-minB[2]);
          nop++; dejavup[ind] = 1; 
          indirW.push_back(et); indirZ.push_back(now); 
        }
      }
    }
    vectOfBoxes.push_back(boxes);
  }
  wallpts->reAllocMat(nop, 6); lmax->resize(nop);
  xw2 = wallpts->begin(1); yw2 = wallpts->begin(2); zw2 = wallpts->begin(3);
  sxw2 = wallpts->begin(4); syw2 = wallpts->begin(5); szw2 = wallpts->begin(6);
  lmaxp = lmax->begin();

  if (nop == 0)
  {
    E_Int nboxes = vectOfBoxes.size();
    for (E_Int v0 = 0; v0 < nboxes; v0++)
    {
      vector<BBox3DType*>& boxes = vectOfBoxes[v0];
      E_Int size = boxes.size();
      for (E_Int v = 0; v < size; v++) delete boxes[v];
    }
    for (E_Int v = 0; v < nzones; v++) 
    {
      E_Int ncells = ncellst[v];
      E_Float* distancep = distances[v]->begin();
      for (E_Int ind = 0; ind < ncells; ind++)
        distancep[ind] = sqrt(distancep[ind]); 
    }
    return;
  }

  // Build the kdtree
  ArrayAccessor<FldArrayF> coordAcc(*wallpts, 1,2,3);
  KdTree<FldArrayF> kdt(coordAcc, E_EPSILON);
  // Build the bbtrees
  vector<K_SEARCH::BbTree3D*> vectOfBBTrees;
  for (E_Int v = 0; v < nwalls; v++)
  {
    K_SEARCH::BbTree3D* bbtree = new K_SEARCH::BbTree3D(vectOfBoxes[v]);
    vectOfBBTrees.push_back(bbtree);
  }

  /* Compute the distance with orthogonal projection */
  E_Float pt[3];
  E_Float p0[3]; E_Float p1[3]; E_Float p2[3]; E_Float p[3];
  vector<E_Int> indicesBB; vector<E_Int> candidates;

  for (E_Int v = 0; v < nzones; v++)
  {
    E_Float* xt = fields[v]->begin(posx);
    E_Float* yt = fields[v]->begin(posy);
    E_Float* zt = fields[v]->begin(posz);
    E_Float* distancep = distances[v]->begin();
    E_Int npts = distances[v]->getSize();

#pragma omp parallel for default(shared) private(minB, maxB, pt, p0, p1, p2, p, indicesBB, candidates) schedule(dynamic)
    for (E_Int ind = 0; ind < npts; ind++)
    {
      E_Int ret, ets, nows; E_Float rxsav, rysav, rzsav;
      E_Float dist, xp, yp, zp, rx, ry, rz, rad, sx, sy, sz;
      E_Float distmin;
      E_Int indw2;
      ets = -1; nows = -1; rxsav = 0.; rysav = 0.; rzsav = 0.; 
      sx = 0.; sy = 0.; sz = 0.;
      indicesBB.clear(); candidates.clear();
      distmin = distancep[ind];
      pt[0] = xt[ind]; pt[1] = yt[ind]; pt[2] = zt[ind];
      // recherche du sommet P' des parois le plus proche de P
      indw2 = kdt.getClosest(pt);
      // calcul de la bounding box de la sphere de rayon PP'
      rx = xw2[indw2]-pt[0]; ry = yw2[indw2]-pt[1]; rz = zw2[indw2]-pt[2];
      dist = rx*rx + ry*ry + rz*rz; rad = sqrt(dist);

      rad = sqrt(dist);
      E_Float A = 1./(10.*lmaxp[indw2]);
      E_Float rad2 = exp(-A*rad);
      E_Float alpha = 1.-rad2;
      E_Float R = rad*rad2;
      E_Float xQ = pt[0] + alpha*(xw2[indw2]-pt[0]);
      E_Float yQ = pt[1] + alpha*(yw2[indw2]-pt[1]);
      E_Float zQ = pt[2] + alpha*(zw2[indw2]-pt[2]);
      minB[0] = xQ-R; minB[1] = yQ-R; minB[2] = zQ-R;
      maxB[0] = xQ+R; maxB[1] = yQ+R; maxB[2] = zQ+R;

      if (dist < distmin) 
      {
        distmin = dist; rxsav = rx; rysav = ry; rzsav = rz; //radsav = rad;
        sx = sxw2[indw2]; sy = syw2[indw2]; sz = szw2[indw2]; 
      }
      // calcul des cellules intersectantes
      for (E_Int now = 0; now < nwalls; now++)
      {
        indicesBB.clear(); candidates.clear();
        K_SEARCH::BbTree3D* bbtree = vectOfBBTrees[now];
        bbtree->getOverlappingBoxes(minB, maxB, indicesBB);
        FldArrayF* fieldv = fieldsw[now];
        E_Int posxw = posxv[now]; E_Int posyw = posyv[now]; E_Int poszw = poszv[now]; E_Int poscw = poscv[now];
        E_Float* xw = fieldv->begin(posxw);
        E_Float* yw = fieldv->begin(posyw);
        E_Float* zw = fieldv->begin(poszw);
        E_Float* cellnw = fieldv->begin(poscw);
        FldArrayI& cnloc = *cntw[now];
        E_Int nbb = indicesBB.size();
        for (E_Int i = 0; i < nbb; i++)
        {
          E_Int et = indicesBB[i];
          E_Int ind10 = cnloc(et,1)-1;
          E_Int ind20 = cnloc(et,2)-1;
          E_Int ind30 = cnloc(et,3)-1;
          E_Float prod = cellnw[ind10]*cellnw[ind20]*cellnw[ind30];
          if (prod == 1.) candidates.push_back(et);
        }

        ret = K_COMPGEOM::projectOrthoPrecond(pt[0], pt[1], pt[2], xw, yw, zw, 
                                              candidates, cnloc, xp, yp, zp,
                                              p0, p1, p2, p);
        if (ret != -1)
        {
          rx = xp-pt[0]; ry = yp-pt[1]; rz = zp-pt[2];
          dist = rx*rx + ry*ry + rz*rz;    
          if (dist < distmin) 
          {
            distmin = dist; rxsav = rx; rysav = ry; rzsav = rz; 
            ets = ret; nows = now;
          }
        }
      } // boucle sur les parois de projection
      if (ets != -1)
      {
        FldArrayF* fieldv = fieldsw[nows];
        E_Float* sxw = fieldv->begin(possx);
        E_Float* syw = fieldv->begin(possy);
        E_Float* szw = fieldv->begin(possz);
        FldArrayI& cnloc = *cntw[nows];
        E_Int ind10 = cnloc(ets,1)-1; 
        E_Int ind20 = cnloc(ets,2)-1;
        E_Int ind30 = cnloc(ets,3)-1;
        sx = (sxw[ind10]+sxw[ind20]+sxw[ind30]);
        sy = (syw[ind10]+syw[ind20]+syw[ind30]);
        sz = (szw[ind10]+szw[ind20]+szw[ind30]);
      }
      E_Float sgn = rxsav*sx + rysav*sy + rzsav*sz;
      E_Float dp = sqrt(rxsav*rxsav+rysav*rysav+rzsav*rzsav);
      if (sgn < 0.) distancep[ind] = dp;
      else distancep[ind] = -dp;
    } // fin boucle sur les pts sur lesquels la distance est calculee
  }// fin boucle sur les zones ou la distance est a calculer

  // Cleaning
  E_Int nboxes = vectOfBoxes.size();
  for (E_Int v0 = 0; v0 < nboxes; v0++)
  {
    vector<BBox3DType*>& boxes = vectOfBoxes[v0];
    E_Int size = boxes.size();
    for (E_Int v = 0; v < size; v++) delete boxes[v];
    delete vectOfBBTrees[v0];
  }
  vectOfBoxes.clear(); vectOfBBTrees.clear();
  delete wallpts; delete lmax;
  return;
}
