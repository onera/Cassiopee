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
/* Calcul de la distance a la paroi par projection orthogonale
   Les corps doivent etre en TRI, avec le cellN localise aux noeuds du 
   maillage TRI */
//=============================================================================
PyObject* K_DIST2WALLS::distance2WallsOrtho(PyObject* self, PyObject* args)
{
  PyObject *centers, *bodiesC;
  E_Int isminortho;
  E_Int isIBM_F1;
  E_Float dTarget;
  if (!PyArg_ParseTuple(args, "OOiid", &centers, &bodiesC, &isminortho, &isIBM_F1, &dTarget)) return NULL;
  
  if (PyList_Check(centers) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "distance2Walls: 1st argument must be a list.");
    return NULL;
  }

  if (PyList_Check(bodiesC) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "distance2Walls: 2nd argument must be a list.");
    return NULL;
  }
  
  // Maillage en noeuds
  // les coordonnees doivent etre en premier
  vector<PyObject*> objsn, objun;
  vector<E_Int> resn;
  vector<char*> structVarStringn; vector<char*> unstrVarStringn;
  vector<FldArrayF*> structFn; vector<FldArrayF*> unstrFn;
  vector<E_Int> nitn; vector<E_Int> njtn; vector<E_Int> nktn;
  vector<FldArrayI*> cntn;
  vector<char*> eltTypen;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false;
  E_Boolean skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    centers, resn, structVarStringn, unstrVarStringn,
    structFn, unstrFn, nitn, njtn, nktn, cntn, eltTypen, objsn, objun, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nsn = structFn.size(); E_Int nun = unstrFn.size();
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "distance2Walls: invalid list of mesh arrays.");
    for (E_Int nos = 0; nos < nsn; nos++)
      RELEASESHAREDS(objsn[nos], structFn[nos]);
    for (E_Int nos = 0; nos < nun; nos++)
      RELEASESHAREDU(objun[nos], unstrFn[nos], cntn[nos]);
    return NULL;
  }
  E_Int posx=-1, posy=-1, posz=-1;
  vector<E_Int> ncellss; vector<E_Int> ncellsu;
  vector<E_Int> posflags; vector<E_Int> posflagu; 

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
      PyErr_SetString( PyExc_TypeError,
        "distance2Walls: coordinates must be located at same position for all zones.");
      for (E_Int nos = 0; nos < nsn; nos++)
        RELEASESHAREDS(objsn[nos], structFn[nos]);
      for (E_Int nos = 0; nos < nun; nos++)
        RELEASESHAREDU(objun[nos], unstrFn[nos], cntn[nos]);
      return NULL;
    }
    E_Int posfi = K_ARRAY::isNamePresent("flag",structVarStringn[i]);
    posfi++;
    posflags.push_back(posfi);
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
        "distance2Walls: coordinates must be located at same position for all zones.");
      for (E_Int nos = 0; nos < nsn; nos++)
        RELEASESHAREDS(objsn[nos], structFn[nos]);
      for (E_Int nos = 0; nos < nun; nos++)
        RELEASESHAREDU(objun[nos], unstrFn[nos], cntn[nos]);
      return NULL;
    }
    E_Int posfi = K_ARRAY::isNamePresent("flag",unstrVarStringn[i]);
    posfi++;
    posflagu.push_back(posfi);
    ncellsu.push_back(unstrFn[i]->getSize());
  }
  // Extract infos from body surfaces
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
      PyErr_SetString(PyExc_TypeError, "distance2Walls: cellN must be defined for bodies.");
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
    computeOrthoDist(ncellss, posx, posy, posz, posflags, 
                     structFn, posxv, posyv, poszv, poscv, unstrF, cnt,
                     distances, isminortho, isIBM_F1, dTarget);
  if (unstrFn.size() > 0)
    computeOrthoDist(ncellsu, posx, posy, posz, posflagu,
                     unstrFn, posxv, posyv, poszv, poscv, unstrF, cnt, 
                     distancesu, isminortho, isIBM_F1, dTarget);

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
    RELEASESHAREDU(objun[nou], unstrFn[nou],cntn[nou]);
    delete distancesu[nou]; delete cnout;
  }
  return l;
}

//=============================================================================
/* Calcul de la distance a la paroi */
//=============================================================================
void K_DIST2WALLS::computeOrthoDist(
  vector<E_Int>& ncellst,
  E_Int posx, E_Int posy, E_Int posz,  vector<E_Int>& posflag,
  vector<FldArrayF*>& fields, 
  vector<E_Int>& posxv, vector<E_Int>& posyv, vector<E_Int>& poszv, 
  vector<E_Int>& poscv,
  vector<FldArrayF*>& fieldsw, vector<FldArrayI*>& cntw,
  vector<FldArrayF*>& distances,E_Int isminortho, E_Int isIBM_F1, E_Float dTarget)
{
  E_Int nzones = fields.size();
  /* 1 - creation du kdtree et du bbtree */
  typedef K_SEARCH::BoundingBox<3> BBox3DType; 
  vector< vector<BBox3DType*> > vectOfBoxes; // a detruire a la fin
  // allocate kdtree array: kdtree points are cell vertices
  E_Int nwalls = cntw.size();
  E_Int nptsmax = 0;
  E_Int npts_local = 0;
  vector<vector< vector<E_Int>  > >cVE_all;
  vector<E_Int> npts_walls_limit;
  for (E_Int v = 0; v < nwalls; v++)
  {
    nptsmax += fieldsw[v]->getSize();
    #include "mininterf_ortho_npts_limit.h"
  }
  FldArrayF* wallpts = new FldArrayF(nptsmax, 3);
  FldArrayF* lmax = new FldArrayF(nptsmax);

  E_Float* xw2 = wallpts->begin(1);
  E_Float* yw2 = wallpts->begin(2);
  E_Float* zw2 = wallpts->begin(3);
  E_Float* lmaxp = lmax->begin();

  // create kdtree elements and bbtree information
  E_Int nop = 0;
  E_Int ind;
  vector<FldArrayF> bboxes;
  E_Float minB[3]; E_Float maxB[3];
  
  for (E_Int now = 0; now < nwalls; now++)
  {
    FldArrayF* fieldv = fieldsw[now];
    E_Int posxw = posxv[now]; E_Int posyw = posyv[now]; E_Int poszw = poszv[now];  
    E_Float* xw = fieldv->begin(posxw);
    E_Float* yw = fieldv->begin(posyw);
    E_Float* zw = fieldv->begin(poszw);
    E_Int poscw = poscv[now]; E_Float* cellnw = fieldv->begin(poscw);
    E_Int npts = fieldv->getSize();
    FldArrayI& cnloc = *cntw[now];
    E_Int nelts = cnloc.getSize(); E_Int nvert = cnloc.getNfld();
    vector<BBox3DType*> boxes(nelts); // liste des bbox de ts les elements de la paroi courante
    FldArrayF bbox(nelts, 6); // xmin, ymin, zmin, xmax, ymax, zmax
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
        ind = cnloc(et, nov)-1;
        if (cellnw[ind] == 1. && dejavup[ind] == 0) 
        {
          xw2[nop] = xw[ind]; yw2[nop] = yw[ind]; zw2[nop] = zw[ind];
          lmaxp[nop] = K_FUNC::E_max(maxB[0]-minB[0], maxB[1]-minB[1], maxB[2]-minB[2]);
          nop++; dejavup[ind] = 1;
        }
      }
    }
    vectOfBoxes.push_back(boxes);
  }
  wallpts->reAllocMat(nop, 3); lmax->resize(nop);
  xw2 = wallpts->begin(1); yw2 = wallpts->begin(2); zw2 = wallpts->begin(3);
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
      for (ind = 0; ind < ncells; ind++)
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
  #pragma omp parallel
  {
    E_Float pt[3]; 
    vector<E_Int> indicesBB; vector<E_Int> candidates;
    E_Float minB[3]; E_Float maxB[3];
    E_Int ret,vw;
    E_Float dist, dx, dy, dz, xp, yp, zp, rx, ry, rz, rad;
    E_Float distmin, prod;
    E_Int et, ind10, indw2, nbb, nvert;
    E_Int posxw, posyw, poszw, poscw;
    E_Float A, rad2, alpha, R, xQ, yQ, zQ, rmax;
    E_Float* xw; E_Float* yw; E_Float* zw; E_Float* cellnw;
    E_Float p0[3]; E_Float p1[3]; E_Float p2[3]; E_Float p[3];

    for (E_Int v = 0; v < nzones; v++)
    {
      E_Float* xt = fields[v]->begin(posx);
      E_Float* yt = fields[v]->begin(posy);
      E_Float* zt = fields[v]->begin(posz);
      E_Float* distancep = distances[v]->begin(); 
      E_Int npts = distances[v]->getSize();
      E_Int isFlagged=false;
      E_Float* flagp = NULL;
      if (posflag[v] > 0 )
      { 
        flagp = fields[v]->begin(posflag[v]);
        isFlagged = true;
      }
      if (isFlagged == true)
      {
        #pragma omp for schedule(dynamic)
        for (E_Int ind = 0; ind < npts; ind++)
        {
          if (flagp[ind] == 0.) { ; }
          else
          {
            #include "algoOrtho.h"
          }
        }   
      }
      else
      {
        #pragma omp for schedule(dynamic)
        for (E_Int ind = 0; ind < npts; ind++)
        {   
          #include "algoOrtho.h"
        }  
      }
    }
  }
  
  // Computes the distance (sqrt)
  #pragma omp parallel
  {
    for (E_Int v = 0; v < nzones; v++)
    {
        E_Float* distancep = distances[v]->begin();
        E_Int npts = distances[v]->getSize();
        #pragma omp for
        for (E_Int ind = 0; ind < npts; ind++)
          distancep[ind] = sqrt(distancep[ind]);
    }
  }

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
