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

// Pour conformisation de l'octree: 27 branches dans l'arbre au lieu de 8

# include "generator.h"
# include <stack>
# include "Nuga/include/BbTree.h"
# include "Search/OctreeNode.h"
# include "Nuga/include/ArrayAccessor.h"

using namespace std;
using namespace K_FLD;
using namespace K_SEARCH;

namespace K_GENERATOR 
{
struct stackData3
{
    OctreeNode* current;
};

OctreeNode* splitNode27(OctreeNode* current, E_Int levelMax, E_Int dim, E_Int& split);
OctreeNode* splitNode9(OctreeNode* current, E_Int levelMax, E_Int& split);

OctreeNode* splitVoisinBoth(OctreeNode* node, E_Int levelMax, E_Int dim, stack<stackData3>& stack);
OctreeNode* splitVoisin9(OctreeNode* node, E_Int levelMax, stack<stackData3>& stack);
OctreeNode* splitVoisin27(OctreeNode* node, E_Int levelMax, 
                          stack<stackData3>& stack);

OctreeNode* addSplitVoisin9(OctreeNode* voisin, E_Int l0, E_Int levelMax, 
                             stack<stackData3>& stack);
OctreeNode* addSplitVoisin27(OctreeNode* voisin, E_Int l0, E_Int levelMax, 
                             stack<stackData3>& stack);

OctreeNode* updateVoisin1_9(OctreeNode* node);
OctreeNode* updateVoisin2_9(OctreeNode* node);
OctreeNode* updateVoisin3_9(OctreeNode* node);
OctreeNode* updateVoisin4_9(OctreeNode* node);
OctreeNode* updateVoisin5_9(OctreeNode* node);
OctreeNode* updateVoisin6_9(OctreeNode* node);
OctreeNode* updateVoisin1_27(OctreeNode* node);
OctreeNode* updateVoisin2_27(OctreeNode* node);
OctreeNode* updateVoisin3_27(OctreeNode* node);
OctreeNode* updateVoisin4_27(OctreeNode* node);
OctreeNode* updateVoisin5_27(OctreeNode* node);
OctreeNode* updateVoisin6_27(OctreeNode* node);

//============================================================================
/* Generation d'un octree � 27 branches 
   a partir d'une liste de surfaces et de snear */
//============================================================================
PyObject* octree3(PyObject* self, PyObject* args)
{
  PyObject *stlArrays, *listOfSnears;
  E_Float dfar; E_Int levelMax;
  PyObject* octant;
  if (!PYPARSETUPLE_(args, OO_ R_ I_ O_,
                    &stlArrays, &listOfSnears, 
                    &dfar, &levelMax, &octant)) return NULL;
  
  if (PyList_Size(stlArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "octree3: 1st argument is an empty list.");
    return NULL;
  }
  if (PyList_Size(listOfSnears) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "octree3: 2nd argument is an empty list.");
    return NULL;
  }
  if (PyList_Size(stlArrays) != PyList_Size(listOfSnears)) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "octree3: 1st and 2nd args must be of same length.");
    return NULL;
  }

  //recuperations des stl
  vector<E_Int> resl;
  vector<char*> structVarString;
  vector<char*> unstrVarString;
  vector<FldArrayF*> structF;
  vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = true;
  E_Bool skipUnstructured = false;
  E_Bool skipDiffVars = true;

  E_Int res = K_ARRAY::getFromArrays(
    stlArrays, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltTypet, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  if (res == -1) 
  {
    for (size_t v = 0; v < structF.size(); v++) RELEASESHAREDS(objst[v], structF[v]);
    for (size_t v = 0; v < unstrF.size(); v++) RELEASESHAREDU(objut[v], unstrF[v], cnt[v]);            
    PyErr_SetString(PyExc_TypeError, 
                    "octree3: 1st arg is not valid.");
    return NULL;
  }
  E_Int nzones = unstrF.size();

  E_Int dim = -1;
  for (E_Int i = 0; i < nzones; i++)
  {
    if ( strcmp(eltTypet[i],"TRI") == 0 ) 
    {
      if ( dim == -1 ) dim = 3;
      else if ( dim != 3) 
      {
        for (size_t v = 0; v < structF.size(); v++) RELEASESHAREDS(objst[v], structF[v]);
        for (size_t v = 0; v < unstrF.size(); v++) RELEASESHAREDU(objut[v], unstrF[v], cnt[v]);                
        PyErr_SetString(PyExc_TypeError, 
                        "octree3: 1st arg must be a list of TRI zones.");
        return NULL;
      }
    }
    else if (strcmp(eltTypet[i],"BAR") == 0 ) 
    {
      if ( dim == -1 ) dim = 2;
      else if ( dim != 2) 
      {
        for (size_t v = 0; v < structF.size(); v++) RELEASESHAREDS(objst[v], structF[v]);
        for (size_t v = 0; v < unstrF.size(); v++) RELEASESHAREDU(objut[v], unstrF[v], cnt[v]);                
        PyErr_SetString(PyExc_TypeError, 
                        "octree3: 1st arg must be a list of BAR zones.");
        return NULL;
      }
    }
    else 
    {
      for (size_t v = 0; v < structF.size(); v++) RELEASESHAREDS(objst[v], structF[v]);
      for (size_t v = 0; v < unstrF.size(); v++) RELEASESHAREDU(objut[v], unstrF[v], cnt[v]);              
      PyErr_SetString(PyExc_TypeError, 
                      "octree3: 1st arg must be a list of TRI or BAR zones.");
      return NULL; 
    }
  }
  E_Int posxi, posyi, poszi;
  vector<E_Int> posxt; vector<E_Int> posyt; vector<E_Int> poszt;
  for (E_Int i = 0; i < nzones; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(unstrVarString[i]);
    posyi = K_ARRAY::isCoordinateYPresent(unstrVarString[i]);
    poszi = K_ARRAY::isCoordinateZPresent(unstrVarString[i]);
    posxi++; posyi++; poszi++;
    posxt.push_back(posxi); posyt.push_back(posyi); poszt.push_back(poszi);
  }

  // recuperation des snears 
  E_Int nsnear = PyList_Size(listOfSnears);
  if ( nzones != nsnear )
  {
    for (size_t v = 0; v < structF.size(); v++) RELEASESHAREDS(objst[v], structF[v]);
    for (size_t v = 0; v < unstrF.size(); v++) RELEASESHAREDU(objut[v], unstrF[v], cnt[v]);            
    PyErr_SetString(PyExc_TypeError, 
                    "octree3: 1st and 2nd args must be consistent.");
    return NULL;
  }    
  PyObject* tpl = NULL;
  vector<E_Float> snears(nzones);
  for (int i = 0; i < nzones; i++)
  {
    tpl = PyList_GetItem(listOfSnears, i);
    if (PyFloat_Check(tpl) == 0 && PyInt_Check(tpl) == 0)
    {
      for (size_t v = 0; v < structF.size(); v++) RELEASESHAREDS(objst[v], structF[v]);
      for (size_t v = 0; v < unstrF.size(); v++) RELEASESHAREDU(objut[v], unstrF[v], cnt[v]);              
      PyErr_SetString(PyExc_TypeError, 
                      "octree3: not a valid value for snear.");
      return NULL;
    }
    else snears[i] = PyFloat_AsDouble(tpl);
  }
  E_Float dfxm = dfar; E_Float dfym = dfar; E_Float dfzm = dfar;
  E_Float dfxp = dfar; E_Float dfyp = dfar; E_Float dfzp = dfar;  
  if ( dim == 2 ) {dfzp = 0.;dfzm = 0.;}

  // creation des bbox trees pour ttes les surfaces
  typedef K_SEARCH::BoundingBox<3>  BBox3DType;
  vector<K_SEARCH::BbTree3D*> bboxtrees(nzones);
  vector< vector<BBox3DType*> > vectOfBBoxes;// pour etre detruit a la fin
  E_Float minB[3];  E_Float maxB[3];
  E_Int posx2, posy2, posz2, nelts2;
  E_Float xmino =  K_CONST::E_MAX_FLOAT;
  E_Float ymino =  K_CONST::E_MAX_FLOAT;
  E_Float zmino =  K_CONST::E_MAX_FLOAT;
  E_Float xmaxo = -K_CONST::E_MAX_FLOAT;
  E_Float ymaxo = -K_CONST::E_MAX_FLOAT;
  E_Float zmaxo = -K_CONST::E_MAX_FLOAT;
  if (dim == 2) {zmino = 0.; zmaxo= 0.;}
  E_Float xminloc, yminloc, zminloc, xmaxloc, ymaxloc, zmaxloc;

  for (E_Int v = 0; v < nzones; v++)
  {
    FldArrayF& f2 = *unstrF[v]; FldArrayI& cn2 = *cnt[v];
    posx2 = posxt[v]; posy2 = posyt[v]; posz2 = poszt[v];
    //bounding box globale ? 
    K_COMPGEOM::boundingBoxUnstruct(f2.getSize(),
                                    f2.begin(posx2), f2.begin(posy2), f2.begin(posz2),
                                    xminloc, yminloc, zminloc, xmaxloc, ymaxloc, zmaxloc);
    xmino = K_FUNC::E_min(xminloc,xmino); xmaxo = K_FUNC::E_max(xmaxloc,xmaxo);
    ymino = K_FUNC::E_min(yminloc,ymino); ymaxo = K_FUNC::E_max(ymaxloc,ymaxo);
    if (dim == 2) {zmino = 0.; zmaxo = 0.;}
    else {zmino = K_FUNC::E_min(zminloc,zmino); zmaxo = K_FUNC::E_max(zmaxloc,zmaxo);}
    // Creation de la bboxtree
    nelts2 = cn2.getSize();
    vector<BBox3DType*> boxes(nelts2);// liste des bbox de ts les elements de a2
    K_FLD::FldArrayF bbox(nelts2,6);// xmin, ymin, zmin, xmax, ymax, zmax
    K_COMPGEOM::boundingBoxOfUnstrCells(
      cn2, f2.begin(posx2), f2.begin(posy2),f2.begin(posz2),bbox);
    E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
    E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
    E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
    for (E_Int et = 0; et < nelts2; et++)
    {
      minB[0] = xminp[et]; minB[1] = yminp[et]; minB[2] = zminp[et];
      maxB[0] = xmaxp[et]; maxB[1] = ymaxp[et]; maxB[2] = zmaxp[et]; 
      boxes[et] = new BBox3DType(minB, maxB);
    }
    // Build the box tree
    bboxtrees[v] = new K_SEARCH::BbTree3D(boxes);
    vectOfBBoxes.push_back(boxes);
  }
  K_ARRAY::cleanUnstrFields(unstrF, cnt, eltTypet);

  // construction de l'octree (octant)
  E_Float Deltax = xmaxo+dfxp-xmino+dfxm;
  E_Float Deltay = ymaxo+dfyp-ymino+dfym;
  E_Float Deltaz = zmaxo+dfzp-zmino+dfzm;
  E_Float Delta = K_FUNC::E_max(Deltax,Deltay); 
  if (dim == 3) Delta = K_FUNC::E_max(Delta,Deltaz);
  Delta = 0.5*Delta;
  E_Float xc = 0.5*(xmaxo+dfxp+xmino-dfxm);
  E_Float yc = 0.5*(ymaxo+dfyp+ymino-dfym);
  E_Float zc = 0.5*(zmaxo+dfzp+zmino-dfzm);
  zmino = 0.; zmaxo = 0.;
  xmino = xc-Delta; ymino = yc-Delta; if (dim == 3) zmino = zc-Delta;
  xmaxo = xc+Delta; ymaxo = yc+Delta; if (dim == 3) zmaxo = zc+Delta;

  // si octant est present, ecrase xmino,...
  if (octant != Py_None && PyList_Check(octant) == true)
  {
    E_Int s = PyList_Size(octant);
    if (s == 2) 
    { 
      xmino = PyFloat_AsDouble(PyList_GetItem(octant,0)); 
      xmaxo = PyFloat_AsDouble(PyList_GetItem(octant,1)); 
    }
    else if (s == 4)
    { 
      xmino = PyFloat_AsDouble(PyList_GetItem(octant,0)); 
      ymino = PyFloat_AsDouble(PyList_GetItem(octant,1)); 
      xmaxo = PyFloat_AsDouble(PyList_GetItem(octant,2)); 
      ymaxo = PyFloat_AsDouble(PyList_GetItem(octant,3)); 
    } 
    else if (s == 6)
    {
      xmino = PyFloat_AsDouble(PyList_GetItem(octant,0)); 
      ymino = PyFloat_AsDouble(PyList_GetItem(octant,1)); 
      zmino = PyFloat_AsDouble(PyList_GetItem(octant,2));
      xmaxo = PyFloat_AsDouble(PyList_GetItem(octant,3)); 
      ymaxo = PyFloat_AsDouble(PyList_GetItem(octant,4));
      zmaxo = PyFloat_AsDouble(PyList_GetItem(octant,5));
    }
  }

  E_Int l0 = 0;// niveau le plus grossier
  OctreeNode* toptree = new OctreeNode(xmino, ymino, zmino, xmaxo-xmino, l0);
  vector<E_Int> indicesBB;
  E_Float snear, dh0, dh;
  E_Float tol = 1.e-6;
  E_Int found;
  E_Int npts = 8; E_Int nelts = 1;
  OctreeNode* current = toptree;
  stack<stackData3> stack;
  stackData3 dataForStack;
  dataForStack.current = current;
  stack.push(dataForStack);
  E_Int incn = 56; if ( dim == 2 ) incn = 12;
  E_Int ince = 26; if ( dim == 2 ) ince = 8;
  E_Int split = 0;
  while (stack.size() != 0) 
  {
    dataForStack = stack.top(); stack.pop();
    current = dataForStack.current;
    dh0 = current->getDh();
    minB[0] = current->getXmin(); maxB[0] = minB[0]+dh0;
    minB[1] = current->getYmin(); maxB[1] = minB[1]+dh0;
    minB[2] = current->getZmin(); maxB[2] = minB[2]+dh0;
    dh = K_FUNC::E_abs(maxB[0]-minB[0]); //dxs3 = dh/3.; dys3 = dxs3; dzs3 = dxs3; if (dim == 2) dzs3 = 0.;
    snear = K_CONST::E_MAX_FLOAT; found = 0;
    for (E_Int v = 0; v < nzones; v++)
    {
      bboxtrees[v]->getOverlappingBoxes(minB, maxB, indicesBB);
      if (indicesBB.size() != 0)
      {
        //regarder si dh > snear ? 
        snear = K_FUNC::E_min(snears[v],snear);
        if  (dh > snear - tol) found = 1;
      }
      indicesBB.clear();
    }

    l0 = current->getLevel();
    if ( found == 1 && l0 < levelMax) //decouper la cellule en 27 si dh > snear
    {
      current = splitNode27(current, levelMax, dim, split);
      if ( split != 1 ) goto next;
      // fils1->voisins
      dataForStack.current = current->getNext1(); stack.push(dataForStack);
      
      // fils2->voisins
      dataForStack.current = current->getNext2(); stack.push(dataForStack);
      
      // fils3->voisins
      dataForStack.current = current->getNext3(); stack.push(dataForStack);
      
      // fils4->voisins
      dataForStack.current = current->getNext4(); stack.push(dataForStack);
     
      // fils5->voisins
      dataForStack.current = current->getNext5(); stack.push(dataForStack);
     
      // fils6->voisins
      dataForStack.current = current->getNext6(); stack.push(dataForStack);
    
      // fils7->voisins
      dataForStack.current = current->getNext7(); stack.push(dataForStack);
     
      // fils8->voisins
      dataForStack.current = current->getNext8(); stack.push(dataForStack);
     
      // fils9->voisins
      dataForStack.current = current->getNext9(); stack.push(dataForStack);
      
      if ( dim == 3 ) 
      {
        // fils10->voisins
        dataForStack.current = current->getNext10(); stack.push(dataForStack);
        dataForStack.current = current->getNext11(); stack.push(dataForStack);
        dataForStack.current = current->getNext12(); stack.push(dataForStack);
        dataForStack.current = current->getNext13(); stack.push(dataForStack);
        dataForStack.current = current->getNext14(); stack.push(dataForStack);
        dataForStack.current = current->getNext15(); stack.push(dataForStack);
        dataForStack.current = current->getNext16(); stack.push(dataForStack);
        dataForStack.current = current->getNext17(); stack.push(dataForStack);
        dataForStack.current = current->getNext18(); stack.push(dataForStack);
        dataForStack.current = current->getNext19(); stack.push(dataForStack);
        dataForStack.current = current->getNext20(); stack.push(dataForStack);
        dataForStack.current = current->getNext21(); stack.push(dataForStack);
        dataForStack.current = current->getNext22(); stack.push(dataForStack);
        dataForStack.current = current->getNext23(); stack.push(dataForStack);
        dataForStack.current = current->getNext24(); stack.push(dataForStack);
        dataForStack.current = current->getNext25(); stack.push(dataForStack);
        dataForStack.current = current->getNext26(); stack.push(dataForStack);
        dataForStack.current = current->getNext27(); stack.push(dataForStack);
      }

      npts = npts + incn; nelts = nelts + ince;
      next:;
    }
  }
  //nettoyages...
  for (E_Int v = 0; v < nzones; v++)
  {
    delete bboxtrees[v]; 
    E_Int nbboxes =  vectOfBBoxes[v].size();
    for (E_Int v2 = 0; v2 < nbboxes; v2++)
      delete vectOfBBoxes[v][v2];
  }
  vectOfBBoxes.clear();
  // Construction du maillage octree
  E_Int size1 = npts; E_Int size2 = nelts;
  FldArrayF* coords = new FldArrayF(size1,3); coords->setAllValuesAt(K_CONST::E_MAX_FLOAT);
  E_Int nvert = 8; if ( dim == 2 ) nvert = 4;
  FldArrayI* cn = new FldArrayI(size2,nvert); // HEXA ou QUAD
  E_Float* xt = coords->begin(1);
  E_Float* yt = coords->begin(2);
  E_Float* zt = coords->begin(3);
  E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2);
  E_Int* cn3 = cn->begin(3); E_Int* cn4 = cn->begin(4);
  E_Int* cn5 = NULL; E_Int* cn6 = NULL;
  E_Int* cn7 = NULL; E_Int* cn8 = NULL;
  
  if ( dim == 3 )
  { cn5 = cn->begin(5); cn6 = cn->begin(6); cn7 = cn->begin(7); cn8 = cn->begin(8); }

  dataForStack.current = toptree; stack.push(dataForStack);
  E_Int et = 0;
  E_Int ind = 0; 
  while ( stack.size() != 0 ) 
  {
    xt = coords->begin(1); yt = coords->begin(2); zt = coords->begin(3);
    cn1 = cn->begin(1); cn2 = cn->begin(2); cn3 = cn->begin(3); cn4 = cn->begin(4);
    if ( dim == 3) {cn5 = cn->begin(5); cn6 = cn->begin(6); cn7 = cn->begin(7); cn8 = cn->begin(8);}
    dataForStack = stack.top(); stack.pop();
    current = dataForStack.current;
    dh0 = current->getDh();
    xmino = current->getXmin(); xmaxo = xmino+dh0;//current->getXmax();
    ymino = current->getYmin(); ymaxo = ymino+dh0;//current->getYmax();
    zmino = current->getZmin(); zmaxo = zmino+dh0;//current->getZmax();

    // feuille->inserer dans le maillage
    if (  current->getNext1() == NULL ) 
    { 
      //creation des 9 ou 27 sommets de la cellule 1
      xt[ind] = xmino; yt[ind] = ymino; zt[ind] = zmino; cn1[et] = ind+1; ind++; 
      xt[ind] = xmaxo; yt[ind] = ymino; zt[ind] = zmino; cn2[et] = ind+1; ind++;
      xt[ind] = xmaxo; yt[ind] = ymaxo; zt[ind] = zmino; cn3[et] = ind+1; ind++;
      xt[ind] = xmino; yt[ind] = ymaxo; zt[ind] = zmino; cn4[et] = ind+1; ind++;
      if ( dim == 3 ) 
      {
        xt[ind] = xmino; yt[ind] = ymino; zt[ind] = zmaxo; cn5[et] = ind+1; ind++;
        xt[ind] = xmaxo; yt[ind] = ymino; zt[ind] = zmaxo; cn6[et] = ind+1; ind++;
        xt[ind] = xmaxo; yt[ind] = ymaxo; zt[ind] = zmaxo; cn7[et] = ind+1; ind++;
        xt[ind] = xmino; yt[ind] = ymaxo; zt[ind] = zmaxo; cn8[et] = ind+1; ind++;
      }
      et++;

      if ( ind + nvert >= size1 ) 
      { size1 = size1+npts;  coords->reAllocMat(size1,3); }
      if ( et + 1 >= size2 ) 
      { size2 = size2 + nelts; cn->reAllocMat(size2,nvert); }
    }
    else 
    {
      dataForStack.current = current->getNext1(); stack.push(dataForStack);
      dataForStack.current = current->getNext2(); stack.push(dataForStack);
      dataForStack.current = current->getNext3(); stack.push(dataForStack);
      dataForStack.current = current->getNext4(); stack.push(dataForStack);
      dataForStack.current = current->getNext5(); stack.push(dataForStack);
      dataForStack.current = current->getNext6(); stack.push(dataForStack);
      dataForStack.current = current->getNext7(); stack.push(dataForStack);
      dataForStack.current = current->getNext8(); stack.push(dataForStack);
      dataForStack.current = current->getNext9(); stack.push(dataForStack);

      if ( dim == 3 ) 
      {
        dataForStack.current = current->getNext10(); stack.push(dataForStack);
        dataForStack.current = current->getNext11(); stack.push(dataForStack);
        dataForStack.current = current->getNext12(); stack.push(dataForStack);
        dataForStack.current = current->getNext13(); stack.push(dataForStack);
        dataForStack.current = current->getNext14(); stack.push(dataForStack);
        dataForStack.current = current->getNext15(); stack.push(dataForStack);
        dataForStack.current = current->getNext16(); stack.push(dataForStack);
        dataForStack.current = current->getNext17(); stack.push(dataForStack);
        dataForStack.current = current->getNext18(); stack.push(dataForStack);
        dataForStack.current = current->getNext19(); stack.push(dataForStack);
        dataForStack.current = current->getNext20(); stack.push(dataForStack);
        dataForStack.current = current->getNext21(); stack.push(dataForStack);
        dataForStack.current = current->getNext22(); stack.push(dataForStack);
        dataForStack.current = current->getNext23(); stack.push(dataForStack);
        dataForStack.current = current->getNext24(); stack.push(dataForStack);
        dataForStack.current = current->getNext25(); stack.push(dataForStack);
        dataForStack.current = current->getNext26(); stack.push(dataForStack);
        dataForStack.current = current->getNext27(); stack.push(dataForStack);
      }
    }
    delete current;
  }
  toptree = NULL;
  coords->reAllocMat(ind,3); cn->reAllocMat(et,nvert);
  char eltType[8];
  if (dim == 2) strcpy(eltType, "QUAD");
  else strcpy(eltType, "HEXA");
  K_CONNECT::cleanConnectivity(1, 2, 3, 1.e-6, eltType,*coords, *cn);
    
  // buildArray
  tpl = K_ARRAY::buildArray(*coords, "x,y,z", *cn, -1, eltType, false);
  // nettoyage
  delete coords; delete cn;
  for (size_t v = 0; v < structF.size(); v++) RELEASESHAREDS(objst[v], structF[v]);
  for (size_t v = 0; v < unstrF.size(); v++) RELEASESHAREDU(objut[v], unstrF[v], cnt[v]);            

  return tpl;
}
//=============================================================================
OctreeNode* splitVoisinBoth(OctreeNode* node, E_Int levelMax, E_Int dim, 
                            stack<stackData3>& stack)
{
  if (dim == 2) return splitVoisin9(node, levelMax, stack);
  else return splitVoisin27(node, levelMax, stack);
}
//=============================================================================
OctreeNode* splitVoisin27(OctreeNode* node, E_Int levelMax, 
                          stack<stackData3>& stack)
{
  OctreeNode** nodep = &node;
  E_Int l0 = (*nodep)->getLevel();
  OctreeNode* voisin = NULL; OctreeNode** voisinp = NULL;
 
  //voisin1 de v1
  *nodep = updateVoisin1_27(*nodep);
  voisin = (*nodep)->getVoisin1(); voisinp = &voisin;
  *voisinp = addSplitVoisin27(*voisinp, l0, levelMax, stack);
  (*nodep)->setVoisin1(*voisinp);
//   *nodep = updateVoisin1_27(*nodep);

  //voisin2 de v1
  *nodep = updateVoisin2_27(*nodep);
  voisin = (*nodep)->getVoisin2(); voisinp = &voisin;
  *voisinp = addSplitVoisin27(*voisinp, l0, levelMax, stack);
  (*nodep)->setVoisin2(*voisinp);
//   *nodep = updateVoisin2_27(*nodep);

  //voisin3 de v1
  *nodep = updateVoisin3_27(*nodep);
  voisin = (*nodep)->getVoisin3(); voisinp = &voisin;
  *voisinp = addSplitVoisin27(*voisinp, l0, levelMax, stack);
  (*nodep)->setVoisin3(*voisinp);
//   *nodep = updateVoisin3_27(*nodep);

  //voisin4 de v1
  *nodep = updateVoisin4_27(*nodep);
  voisin = (*nodep)->getVoisin4(); voisinp = &voisin;
  *voisinp = addSplitVoisin27(*voisinp, l0, levelMax, stack);
  (*nodep)->setVoisin4(*voisinp);
//   *nodep = updateVoisin4_27(*nodep);

  //voisin5 de v1
  *nodep = updateVoisin5_27(*nodep);
  voisin = (*nodep)->getVoisin5(); voisinp = &voisin;
  *voisinp = addSplitVoisin27(*voisinp, l0, levelMax, stack);
  (*nodep)->setVoisin5(*voisinp);
//   *nodep = updateVoisin5_27(*nodep);

  //voisin6 de v1
  *nodep = updateVoisin6_27(*nodep);
  voisin = (*nodep)->getVoisin6(); voisinp = &voisin;
  *voisinp = addSplitVoisin27(*voisinp, l0, levelMax, stack);
  (*nodep)->setVoisin6(*voisinp);
//   *nodep = updateVoisin6_27(*nodep);

  return *nodep;
}
//=============================================================================
OctreeNode* addSplitVoisin27(OctreeNode* voisin, E_Int l0, E_Int levelMax, 
                             stack<stackData3>& stack)
{
  E_Int dim = 3; 
  stackData3 dataForStack;
  OctreeNode* fils = NULL;
  OctreeNode** filsp = NULL;
  E_Int split = 0;   E_Int dl = 0;

  if ( voisin != NULL )
  {
    OctreeNode** voisinp = &voisin;
    dl = l0-(*voisinp)->getLevel();
    if ( dl > 1 )
    { 
      *voisinp = splitNode27(*voisinp, levelMax, dim, split);
      if ( split == 1 )
      {
        fils = (*voisinp)->getNext1(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack);
        (*voisinp)->setNext1(*filsp); 
        dataForStack.current = (*voisinp)->getNext1(); stack.push(dataForStack);

        fils = (*voisinp)->getNext2(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack);
        (*voisinp)->setNext2(*filsp); 
        dataForStack.current = (*voisinp)->getNext2(); stack.push(dataForStack);

        fils = (*voisinp)->getNext3(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack);
        (*voisinp)->setNext3(*filsp);
        dataForStack.current = (*voisinp)->getNext3(); stack.push(dataForStack);

        fils = (*voisinp)->getNext4(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext4(*filsp);
        dataForStack.current = (*voisinp)->getNext4(); stack.push(dataForStack);

        fils = (*voisinp)->getNext5(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack);
        (*voisinp)->setNext5(*filsp); 
        dataForStack.current = (*voisinp)->getNext5(); stack.push(dataForStack);

        fils = (*voisinp)->getNext6(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack);
        (*voisinp)->setNext6(*filsp); 
        dataForStack.current = (*voisinp)->getNext6(); stack.push(dataForStack);

        fils = (*voisinp)->getNext7(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack);
        (*voisinp)->setNext7(*filsp);
        dataForStack.current = (*voisinp)->getNext7(); stack.push(dataForStack);

        fils = (*voisinp)->getNext8(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext8(*filsp);
        dataForStack.current = (*voisinp)->getNext8(); stack.push(dataForStack);

        fils = (*voisinp)->getNext9(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext9(*filsp);
        dataForStack.current = (*voisinp)->getNext9(); stack.push(dataForStack);

        fils = (*voisinp)->getNext10(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext10(*filsp);
        dataForStack.current = (*voisinp)->getNext10(); stack.push(dataForStack);

        fils = (*voisinp)->getNext11(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext11(*filsp);
        dataForStack.current = (*voisinp)->getNext11(); stack.push(dataForStack);

        fils = (*voisinp)->getNext12(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext12(*filsp);
        dataForStack.current = (*voisinp)->getNext12(); stack.push(dataForStack);

        fils = (*voisinp)->getNext13(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext13(*filsp);
        dataForStack.current = (*voisinp)->getNext13(); stack.push(dataForStack);

        fils = (*voisinp)->getNext14(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext14(*filsp);
        dataForStack.current = (*voisinp)->getNext14(); stack.push(dataForStack);

        fils = (*voisinp)->getNext15(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext15(*filsp);
        dataForStack.current = (*voisinp)->getNext15(); stack.push(dataForStack);

        fils = (*voisinp)->getNext16(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext16(*filsp);
        dataForStack.current = (*voisinp)->getNext16(); stack.push(dataForStack);

        fils = (*voisinp)->getNext17(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext17(*filsp);
        dataForStack.current = (*voisinp)->getNext17(); stack.push(dataForStack);

        fils = (*voisinp)->getNext18(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext18(*filsp);
        dataForStack.current = (*voisinp)->getNext18(); stack.push(dataForStack);

        fils = (*voisinp)->getNext19(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext19(*filsp);
        dataForStack.current = (*voisinp)->getNext19(); stack.push(dataForStack);

        fils = (*voisinp)->getNext20(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext20(*filsp);
        dataForStack.current = (*voisinp)->getNext20(); stack.push(dataForStack);

        fils = (*voisinp)->getNext21(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext21(*filsp);
        dataForStack.current = (*voisinp)->getNext21(); stack.push(dataForStack);

        fils = (*voisinp)->getNext22(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext22(*filsp);
        dataForStack.current = (*voisinp)->getNext22(); stack.push(dataForStack);

        fils = (*voisinp)->getNext23(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext23(*filsp);
        dataForStack.current = (*voisinp)->getNext23(); stack.push(dataForStack);

        fils = (*voisinp)->getNext24(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext24(*filsp);
        dataForStack.current = (*voisinp)->getNext24(); stack.push(dataForStack);

        fils = (*voisinp)->getNext25(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext25(*filsp);
        dataForStack.current = (*voisinp)->getNext25(); stack.push(dataForStack);

        fils = (*voisinp)->getNext26(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext26(*filsp);
        dataForStack.current = (*voisinp)->getNext26(); stack.push(dataForStack);

        fils = (*voisinp)->getNext27(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext27(*filsp);
        dataForStack.current = (*voisinp)->getNext27(); stack.push(dataForStack);

        return *voisinp;
      }
    }//dl > 1
  }
  return voisin;
}
//=============================================================================
OctreeNode* addSplitVoisin9(OctreeNode* voisin, E_Int l0, E_Int levelMax, 
                             stack<stackData3>& stack)
{
  E_Int dim = 2;  
  stackData3 dataForStack;
  OctreeNode* fils = NULL;
  OctreeNode** filsp = NULL;
  E_Int split = 0;   E_Int dl = 0;
  if ( voisin != NULL )
  {
    OctreeNode** voisinp = &voisin;
    dl = l0-(*voisinp)->getLevel();
    if ( dl > 1 )
    { 
      *voisinp = splitNode27(*voisinp, levelMax, dim, split);
      split = 0;
      if ( split == 1 )
      {
        fils = (*voisinp)->getNext1(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack);
        (*voisinp)->setNext1(*filsp); 
        dataForStack.current = (*voisinp)->getNext1(); stack.push(dataForStack);

        fils = (*voisinp)->getNext2(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack);
        (*voisinp)->setNext2(*filsp); 
        dataForStack.current = (*voisinp)->getNext2(); stack.push(dataForStack);

        fils = (*voisinp)->getNext3(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack);
        (*voisinp)->setNext3(*filsp);
        dataForStack.current = (*voisinp)->getNext3(); stack.push(dataForStack);

        fils = (*voisinp)->getNext4(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext4(*filsp);
        dataForStack.current = (*voisinp)->getNext4(); stack.push(dataForStack);

        fils = (*voisinp)->getNext5(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack);
        (*voisinp)->setNext5(*filsp); 
        dataForStack.current = (*voisinp)->getNext5(); stack.push(dataForStack);

        fils = (*voisinp)->getNext6(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack);
        (*voisinp)->setNext6(*filsp); 
        dataForStack.current = (*voisinp)->getNext6(); stack.push(dataForStack);

        fils = (*voisinp)->getNext7(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack);
        (*voisinp)->setNext7(*filsp);
        dataForStack.current = (*voisinp)->getNext7(); stack.push(dataForStack);

        fils = (*voisinp)->getNext8(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext8(*filsp);
        dataForStack.current = (*voisinp)->getNext8(); stack.push(dataForStack);

        fils = (*voisinp)->getNext9(); filsp = &fils;
        *filsp = splitVoisinBoth(*filsp, levelMax, dim, stack); 
        (*voisinp)->setNext9(*filsp);
        dataForStack.current = (*voisinp)->getNext9(); stack.push(dataForStack);

        return *voisinp;
      }
    }//dl > 1
  }
  return voisin;
}

//=============================================================================
OctreeNode* splitVoisin9(OctreeNode* node, E_Int levelMax, stack<stackData3>& stack)
{
  OctreeNode** nodep = &node;
  E_Int l0 = (*nodep)->getLevel();
  OctreeNode* voisin = NULL; OctreeNode** voisinp = NULL;
  //voisin1 de v1
//   *nodep = updateVoisin1_9(*nodep);
  voisin =  (*nodep)->getVoisin1(); voisinp = &voisin;
  *voisinp = addSplitVoisin9(*voisinp, l0, levelMax, stack);
  (*nodep)->setVoisin1(*voisinp);

  //voisin2 de v1
//   *nodep = updateVoisin2_9(*nodep);
  voisin =  (*nodep)->getVoisin2(); voisinp = &voisin;
  *voisinp = addSplitVoisin9(*voisinp, l0, levelMax, stack);
  (*nodep)->setVoisin2(*voisinp);

  //voisin3 de v1
//   *nodep = updateVoisin3_9(*nodep);
  voisin =  (*nodep)->getVoisin3(); voisinp = &voisin;
  *voisinp = addSplitVoisin9(*voisinp, l0, levelMax, stack);
  (*nodep)->setVoisin3(*voisinp);

  //voisin4 de v1
//   *nodep = updateVoisin4_9(*nodep);
  voisin =  (*nodep)->getVoisin4(); voisinp = &voisin;
  *voisinp = addSplitVoisin9(*voisinp, l0, levelMax, stack);
  (*nodep)->setVoisin4(*voisinp);
  return *nodep;
}
//=============================================================================
OctreeNode* updateVoisin1_9(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData3> stack;
  OctreeNode* voisin = node->getVoisin1();
  E_Int lvl0 = node->getLevel();
  vector<OctreeNode*> candidats;
  stackData3 dataForStack;
  E_Float xmin = node->getXmin(); //E_Float xmax = xmin+node->getDh();
  E_Float ymin = node->getYmin(); E_Float ymax = ymin+node->getDh();
  E_Float xminv, yminv, xmaxv, ymaxv;
  if (voisin == NULL) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if (voisin->getNext1() == NULL ) return *nodep;

  if ( voisin->getNext3() == NULL && voisin->getNext6() == NULL && voisin->getNext9() == NULL ) return *nodep;
  if ( voisin->getNext3() != NULL ) 
  {
    dataForStack.current = voisin->getNext3(); stack.push(dataForStack);
  }
  if ( voisin->getNext6() != NULL ) 
  {
    dataForStack.current = voisin->getNext6(); stack.push(dataForStack);
  }
  if ( voisin->getNext9() != NULL ) 
  {
    dataForStack.current = voisin->getNext9(); stack.push(dataForStack);
  }
  stackData3 s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    xminv = stackp->getXmin(); xmaxv = xminv+stackp->getDh();
    yminv = stackp->getYmin(); ymaxv = yminv+stackp->getDh();
    //sous-voisin3
    if ( stackp->getNext3() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext3(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmin-xmaxv,eps) == true )
      {
        if ( yminv < ymax - eps && ymaxv > ymin + eps ) candidats.push_back(stackp);
      }
    }
    //sous-voisin6
    if ( stackp->getNext6() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext6();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmin-xmaxv,eps) == true )
      {
        if ( yminv < ymax - eps && ymaxv > ymin + eps ) candidats.push_back(stackp);
      }
    }
    //sous-voisin9
    if ( stackp->getNext9() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext9();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmin-xmaxv,eps) == true )
      {
        if ( yminv < ymax - eps && ymaxv > ymin + eps ) candidats.push_back(stackp);
      }
    }
  }
  E_Int l0 = -1;

  E_Int ncandidats = candidats.size();
  for (E_Int v = 0; v < ncandidats; v++)
  {
    if ( candidats[v]->getLevel() > l0 ) 
    {l0 = candidats[v]->getLevel(); (*nodep)->setVoisin1(candidats[v]);}
  }

  return *nodep;
}
//=============================================================================
OctreeNode* updateVoisin2_9(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData3> stack;
  OctreeNode* voisin = node->getVoisin2();
  E_Int lvl0 = node->getLevel();
  vector<OctreeNode*> candidats;
  stackData3 dataForStack;
  E_Float xmin = node->getXmin(); E_Float xmax = xmin+node->getDh();
  E_Float ymin = node->getYmin(); E_Float ymax = ymin+node->getDh();
  E_Float xminv, yminv, ymaxv;  
  if ( voisin == NULL ) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if (voisin->getNext1() == NULL ) return *nodep;

  if ( voisin->getNext1() == NULL && voisin->getNext4() == NULL && voisin->getNext7() == NULL) return *nodep;
  if ( voisin->getNext1() != NULL ) 
  {
    dataForStack.current = voisin->getNext1(); stack.push(dataForStack);
  }
  if ( voisin->getNext4() != NULL ) 
  {
    dataForStack.current = voisin->getNext4(); stack.push(dataForStack);
  }
  if ( voisin->getNext7() != NULL ) 
  {
    dataForStack.current = voisin->getNext7(); stack.push(dataForStack);
  }
  stackData3 s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    xminv = stackp->getXmin(); //xmaxv = xminv+stackp->getDh();
    yminv = stackp->getYmin(); ymaxv = yminv+stackp->getDh();
    if ( stackp->getNext1() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext1(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmax-xminv,eps) == true )
      {
        if ( yminv < ymax - eps && ymaxv > ymin + eps ) candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext4() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext4();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmax-xminv,eps) == true )
      {
        if ( yminv < ymax - eps && ymaxv > ymin + eps ) candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext7() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext7();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmax-xminv,eps) == true )
      {
        if ( yminv < ymax - eps && ymaxv > ymin + eps ) candidats.push_back(stackp); 
      }
    }
  }
  E_Int l0 = -1;
  E_Int ncandidats = candidats.size();
  for (E_Int v = 0; v < ncandidats; v++)
  {
    if ( candidats[v]->getLevel() > l0 ) 
    {l0 = candidats[v]->getLevel(); (*nodep)->setVoisin2(candidats[v]);}
  }

  return *nodep;
}
//=============================================================================
OctreeNode* updateVoisin3_9(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData3> stack;
  OctreeNode* voisin = node->getVoisin3();//voisin en ymin
  E_Int lvl0 = node->getLevel();
  vector<OctreeNode*> candidats;
  stackData3 dataForStack;
  E_Float xmin = node->getXmin(); E_Float xmax = xmin+node->getDh();
  E_Float ymin = node->getYmin(); //E_Float ymax = ymin+node->getDh();
  E_Float xminv, yminv, xmaxv, ymaxv;
  if ( voisin == NULL ) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if (voisin->getNext1() == NULL ) return *nodep;

  if ( voisin->getNext7() == NULL && voisin->getNext8() == NULL && voisin->getNext9() == NULL ) return *nodep;
  if ( voisin->getNext7() != NULL ) 
  {
    dataForStack.current = voisin->getNext7(); stack.push(dataForStack);
  }
  if ( voisin->getNext8() != NULL ) 
  {
    dataForStack.current = voisin->getNext8(); stack.push(dataForStack);
  }
  if ( voisin->getNext9() != NULL ) 
  {
    dataForStack.current = voisin->getNext9(); stack.push(dataForStack);
  }
  stackData3 s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    xminv = stackp->getXmin(); xmaxv = xminv+stackp->getDh();
    yminv = stackp->getYmin(); ymaxv = yminv+stackp->getDh();
    if ( stackp->getNext7() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext7();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymin-ymaxv,eps) == true )
      {
        if ( xminv < xmax - eps && xmaxv > xmin + eps ) candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext8() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext8(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymin-ymaxv,eps) == true )
      {
        if ( xminv < xmax - eps && xmaxv > xmin + eps ) candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext9() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext9(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymin-ymaxv,eps) == true )
      {
        if ( xminv < xmax - eps && xmaxv > xmin + eps ) candidats.push_back(stackp); 
      }
    }
  }
  E_Int l0 = -1;
  E_Int ncandidats = candidats.size();
  for (E_Int v = 0; v < ncandidats; v++)
  {
    if ( candidats[v]->getLevel() > l0 ) 
    {l0 = candidats[v]->getLevel(); (*nodep)->setVoisin3(candidats[v]);}
  }

  return *nodep;
}
//=============================================================================
OctreeNode* updateVoisin4_9(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData3> stack;
  OctreeNode* voisin = node->getVoisin4();//voisin en ymin
  E_Int lvl0 = node->getLevel();
  vector<OctreeNode*> candidats;
  stackData3 dataForStack;
  E_Float xmin = node->getXmin(); E_Float xmax = xmin+node->getDh();
  E_Float ymin = node->getYmin(); E_Float ymax = ymin+node->getDh();
  E_Float xminv, yminv, xmaxv;
  if (voisin == NULL) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if (voisin->getNext1() == NULL ) return *nodep;
  if ( voisin->getNext1() == NULL && voisin->getNext2() == NULL && voisin->getNext3() == NULL) return *nodep;
  if ( voisin->getNext1() != NULL ) 
  {
    dataForStack.current = voisin->getNext1(); stack.push(dataForStack);
  }
  if ( voisin->getNext2() != NULL ) 
  {
    dataForStack.current = voisin->getNext2(); stack.push(dataForStack);
  }
  if ( voisin->getNext3() != NULL ) 
  {
    dataForStack.current = voisin->getNext3(); stack.push(dataForStack);
  }
  stackData3 s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    xminv = stackp->getXmin(); xmaxv = xminv+stackp->getDh();
    yminv = stackp->getYmin(); //ymaxv = yminv+stackp->getDh();
    if ( stackp->getNext1() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext1();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymax-yminv,eps) == true )
      {
        if ( xminv < xmax - eps && xmaxv > xmin + eps ) candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext2() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext2(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymax-yminv,eps) == true )
      {
        if ( xminv < xmax - eps && xmaxv > xmin + eps ) candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext3() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext3(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymax-yminv,eps) == true )
      {
        if ( xminv < xmax - eps && xmaxv > xmin + eps ) candidats.push_back(stackp); 
      }
    }
  }
  E_Int l0 = -1;
  E_Int ncandidats = candidats.size();
  for (E_Int v = 0; v < ncandidats; v++)
  {
    if ( candidats[v]->getLevel() > l0 ) 
    {l0 = candidats[v]->getLevel(); (*nodep)->setVoisin4(candidats[v]);}
  }

  return *nodep;
}

//=============================================================================
/* Decoupage du noeud courant en 9 ou 27 noeuds
   Si equilibrage de l octree : rajoute les voisins du noeud courant et les
   autres fils voisins a chaque noeud fils
*/
//=============================================================================
OctreeNode* splitNode27(OctreeNode* current, E_Int levelMax, 
                        E_Int dim, E_Int& split)
{
  if ( dim == 2 ) return splitNode9(current, levelMax, split);
  else 
  {
    split = 0;
    if ( current->getLevel() == levelMax ) return current;
    if ( current == NULL ) return current;
    if ( current->getNext1() != NULL ) return current;
    split = 1;
    E_Int lvl = current->getLevel();
    OctreeNode** currentp = &current;
    E_Float xmin = (*currentp)->getXmin(); E_Float xmax = xmin+(*currentp)->getDh(); 
    E_Float ymin = (*currentp)->getYmin(); //E_Float ymax = ymin+(*currentp)->getDh();  
    E_Float zmin = (*currentp)->getZmin(); //E_Float zmax = zmin+(*currentp)->getDh();  
    E_Float dxs3 = (xmax-xmin)/3.; E_Float dys3 = dxs3; E_Float dzs3 = dxs3;
    lvl = lvl + 1;
    OctreeNode* node1 = new OctreeNode(xmin,         ymin,         zmin, dzs3, lvl);
    OctreeNode* node2 = new OctreeNode(xmin+dxs3,    ymin,         zmin, dzs3, lvl);
    OctreeNode* node3 = new OctreeNode(xmin+2.*dxs3, ymin,         zmin, dzs3, lvl);
    OctreeNode* node4 = new OctreeNode(xmin,         ymin+dys3,    zmin, dzs3, lvl);
    OctreeNode* node5 = new OctreeNode(xmin+dxs3,    ymin+dys3,    zmin, dzs3, lvl);
    OctreeNode* node6 = new OctreeNode(xmin+2.*dxs3, ymin+dys3,    zmin, dzs3, lvl);
    OctreeNode* node7 = new OctreeNode(xmin,         ymin+2.*dys3, zmin, dzs3, lvl);
    OctreeNode* node8 = new OctreeNode(xmin+dxs3,    ymin+2.*dys3, zmin, dzs3, lvl);
    OctreeNode* node9 = new OctreeNode(xmin+2.*dxs3, ymin+2.*dys3, zmin, dzs3, lvl);
   
    OctreeNode* node10 = new OctreeNode(xmin,        ymin,        zmin+dzs3, dzs3, lvl);
    OctreeNode* node11 = new OctreeNode(xmin+dxs3,   ymin,        zmin+dzs3, dzs3, lvl);
    OctreeNode* node12 = new OctreeNode(xmin+2.*dxs3,ymin,        zmin+dzs3, dzs3, lvl);
    OctreeNode* node13 = new OctreeNode(xmin,        ymin+dys3,   zmin+dzs3, dzs3, lvl);
    OctreeNode* node14 = new OctreeNode(xmin+dxs3,   ymin+dys3,   zmin+dzs3, dzs3, lvl);
    OctreeNode* node15 = new OctreeNode(xmin+2.*dxs3,ymin+dys3,   zmin+dzs3, dzs3, lvl);
    OctreeNode* node16 = new OctreeNode(xmin,        ymin+2.*dys3,zmin+dzs3, dzs3, lvl);
    OctreeNode* node17 = new OctreeNode(xmin+dxs3,   ymin+2.*dys3,zmin+dzs3, dzs3, lvl);
    OctreeNode* node18 = new OctreeNode(xmin+2.*dxs3,ymin+2.*dys3,zmin+dzs3, dzs3, lvl); 

    OctreeNode* node19 = new OctreeNode(xmin,        ymin,        zmin+2.*dzs3,dzs3, lvl);
    OctreeNode* node20 = new OctreeNode(xmin+dxs3,   ymin,        zmin+2.*dzs3,dzs3, lvl);
    OctreeNode* node21 = new OctreeNode(xmin+2.*dxs3,ymin,        zmin+2.*dzs3,dzs3, lvl);
    OctreeNode* node22 = new OctreeNode(xmin,        ymin+dys3,   zmin+2.*dzs3,dzs3, lvl);
    OctreeNode* node23 = new OctreeNode(xmin+dxs3,   ymin+dys3,   zmin+2.*dzs3,dzs3, lvl);
    OctreeNode* node24 = new OctreeNode(xmin+2.*dxs3,ymin+dys3,   zmin+2.*dzs3,dzs3, lvl);
    OctreeNode* node25 = new OctreeNode(xmin,        ymin+2.*dys3,zmin+2.*dzs3,dzs3, lvl);
    OctreeNode* node26 = new OctreeNode(xmin+dxs3,   ymin+2.*dys3,zmin+2.*dzs3,dzs3, lvl);
    OctreeNode* node27 = new OctreeNode(xmin+2.*dxs3,ymin+2.*dys3,zmin+2.*dzs3,dzs3, lvl);
    
    OctreeNode** nodep1 = &node1;
    OctreeNode** nodep2 = &node2;
    OctreeNode** nodep3 = &node3;
    OctreeNode** nodep4 = &node4;
    OctreeNode** nodep5 = &node5;
    OctreeNode** nodep6 = &node6;
    OctreeNode** nodep7 = &node7;
    OctreeNode** nodep8 = &node8;
    OctreeNode** nodep9 = &node9;
    OctreeNode** nodep10 = &node10;
    OctreeNode** nodep11 = &node11;
    OctreeNode** nodep12 = &node12;
    OctreeNode** nodep13 = &node13;
    OctreeNode** nodep14 = &node14;
    OctreeNode** nodep15 = &node15;
    OctreeNode** nodep16 = &node16;
    OctreeNode** nodep17 = &node17;
    OctreeNode** nodep18 = &node18;
    OctreeNode** nodep19 = &node19;
    OctreeNode** nodep20 = &node20;
    OctreeNode** nodep21 = &node21;
    OctreeNode** nodep22 = &node22;
    OctreeNode** nodep23 = &node23;
    OctreeNode** nodep24 = &node24;
    OctreeNode** nodep25 = &node25;
    OctreeNode** nodep26 = &node26;
    OctreeNode** nodep27 = &node27;

    (*currentp)->setNext1(*nodep1); 
    (*currentp)->setNext2(*nodep2);   
    (*currentp)->setNext3(*nodep3); 
    (*currentp)->setNext4(*nodep4); 
    (*currentp)->setNext5(*nodep5); 
    (*currentp)->setNext6(*nodep6);   
    (*currentp)->setNext7(*nodep7); 
    (*currentp)->setNext8(*nodep8);
    (*currentp)->setNext9(*nodep9); 
    (*currentp)->setNext10(*nodep10); 
    (*currentp)->setNext11(*nodep11);   
    (*currentp)->setNext12(*nodep12); 
    (*currentp)->setNext13(*nodep13); 
    (*currentp)->setNext14(*nodep14); 
    (*currentp)->setNext15(*nodep15);   
    (*currentp)->setNext16(*nodep16); 
    (*currentp)->setNext17(*nodep17);
    (*currentp)->setNext18(*nodep18); 
    (*currentp)->setNext19(*nodep19); 
    (*currentp)->setNext20(*nodep20);   
    (*currentp)->setNext21(*nodep21); 
    (*currentp)->setNext22(*nodep22); 
    (*currentp)->setNext23(*nodep23); 
    (*currentp)->setNext24(*nodep24);   
    (*currentp)->setNext25(*nodep25); 
    (*currentp)->setNext26(*nodep26);
    (*currentp)->setNext27(*nodep27); 
    return *currentp; 
  }
}
//=============================================================================
OctreeNode* 
splitNode9(OctreeNode* current, E_Int levelMax, E_Int& split)
{
  split = 0;
  if ( current->getLevel() == levelMax ) return current;
  if ( current == NULL ) return current;
  if ( current->getNext1() != NULL ) { split = -1; return current;}
  split = 1;
  E_Int lvl = current->getLevel();
  OctreeNode** currentp = &current;
  E_Float xmin = (*currentp)->getXmin(); E_Float xmax = xmin+(*currentp)->getDh(); 
  E_Float ymin = (*currentp)->getYmin(); //E_Float ymax = ymin+(*currentp)->getDh(); 
  E_Float zmin = (*currentp)->getZmin(); //E_Float zmax = zmin+(*currentp)->getDh(); 

  E_Float dxs3 = (xmax-xmin)/3.; E_Float dys3 = dxs3; 
  lvl = lvl + 1;
  OctreeNode* node1 = new OctreeNode(xmin,        ymin,        zmin, dxs3, lvl);
  OctreeNode* node2 = new OctreeNode(xmin+dxs3,   ymin,        zmin, dxs3, lvl);
  OctreeNode* node3 = new OctreeNode(xmin+2.*dxs3,ymin,        zmin, dxs3, lvl);
  OctreeNode* node4 = new OctreeNode(xmin,        ymin+dys3,   zmin, dxs3, lvl);
  OctreeNode* node5 = new OctreeNode(xmin+dxs3,   ymin+dys3,   zmin, dxs3, lvl);
  OctreeNode* node6 = new OctreeNode(xmin+2.*dxs3,ymin+dys3,   zmin, dxs3, lvl);
  OctreeNode* node7 = new OctreeNode(xmin,        ymin+2.*dys3,zmin, dxs3, lvl);
  OctreeNode* node8 = new OctreeNode(xmin+dxs3,   ymin+2.*dys3,zmin, dxs3, lvl);
  OctreeNode* node9 = new OctreeNode(xmin+2.*dxs3,ymin+2.*dys3,zmin, dxs3, lvl);
  OctreeNode** nodep1 = &node1;
  OctreeNode** nodep2 = &node2;
  OctreeNode** nodep3 = &node3;
  OctreeNode** nodep4 = &node4;
  OctreeNode** nodep5 = &node5;
  OctreeNode** nodep6 = &node6;
  OctreeNode** nodep7 = &node7;
  OctreeNode** nodep8 = &node8;
  OctreeNode** nodep9 = &node9;

  (*currentp)->setNext1(*nodep1); 
  (*currentp)->setNext2(*nodep2);   
  (*currentp)->setNext3(*nodep3); 
  (*currentp)->setNext4(*nodep4); 
  (*currentp)->setNext5(*nodep5); 
  (*currentp)->setNext6(*nodep6);   
  (*currentp)->setNext7(*nodep7); 
  (*currentp)->setNext8(*nodep8); 
  (*currentp)->setNext9(*nodep9); 
  return *currentp;  
}
//=============================================================================
OctreeNode* updateVoisin1_27(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData3> stack;
  OctreeNode* voisin = node->getVoisin1();
  E_Int lvl0 = node->getLevel();
  E_Float xmin = node->getXmin(); //E_Float xmax = xmin+node->getDh();
  E_Float ymin = node->getYmin(); E_Float ymax = ymin+node->getDh();
  E_Float zmin = node->getZmin(); E_Float zmax = zmin+node->getDh();
  E_Float xminv, yminv, zminv, xmaxv, ymaxv, zmaxv;
  vector<OctreeNode*> candidats;
  stackData3 dataForStack;

  if ( voisin == NULL ) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if (voisin->getNext1() == NULL ) return *nodep;
  // 3-6-9-12-15-18-21-24-27
  if ( voisin->getNext3() == NULL && voisin->getNext6() == NULL && voisin->getNext9() == NULL &&
       voisin->getNext12() == NULL && voisin->getNext15() == NULL && voisin->getNext18() == NULL &&
       voisin->getNext21() == NULL && voisin->getNext24() == NULL && voisin->getNext27() == NULL )
    return *nodep;

  if ( voisin->getNext3() != NULL ) 
  {dataForStack.current = voisin->getNext3(); stack.push(dataForStack);}
  if ( voisin->getNext6() != NULL ) 
  {dataForStack.current = voisin->getNext6(); stack.push(dataForStack);}
  if ( voisin->getNext9() != NULL ) 
  {dataForStack.current = voisin->getNext9(); stack.push(dataForStack);}
  if ( voisin->getNext12() != NULL ) 
  {dataForStack.current = voisin->getNext12(); stack.push(dataForStack);}
  if ( voisin->getNext15() != NULL ) 
  {dataForStack.current = voisin->getNext15(); stack.push(dataForStack);}
  if ( voisin->getNext18() != NULL ) 
  {dataForStack.current = voisin->getNext18(); stack.push(dataForStack);}
  if ( voisin->getNext21() != NULL ) 
  {dataForStack.current = voisin->getNext21(); stack.push(dataForStack);}
  if ( voisin->getNext24() != NULL ) 
  {dataForStack.current = voisin->getNext24(); stack.push(dataForStack);}
  if ( voisin->getNext27() != NULL ) 
  {dataForStack.current = voisin->getNext27(); stack.push(dataForStack);}

  stackData3 s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    xminv = stackp->getXmin(); xmaxv = xminv+stackp->getDh();
    yminv = stackp->getYmin(); ymaxv = yminv+stackp->getDh();
    zminv = stackp->getZmin(); zmaxv = zminv+stackp->getDh();

    if ( stackp->getNext3() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext3();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmin-xmaxv,eps) == true )
      {      
        if ( yminv < ymax-eps && ymaxv > ymin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext6() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext6();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmin-xmaxv,eps) == true )
      {
        if ( yminv < ymax - eps && ymaxv > ymin + eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext9() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext9();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmin-xmaxv,eps) == true )
      {
        if ( yminv < ymax - eps && ymaxv > ymin + eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext12() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext12();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmin-xmaxv,eps) == true )
      {
        if ( yminv < ymax - eps && ymaxv > ymin + eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp);
      }
    }
    if ( stackp->getNext15() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext15();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmin-xmaxv,eps) == true )
      {
        if ( yminv < ymax - eps && ymaxv > ymin + eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp);
      }
    }
    if ( stackp->getNext18() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext18();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmin-xmaxv,eps) == true )
      {
        if ( yminv < ymax - eps && ymaxv > ymin + eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp);
      }
    }
    if ( stackp->getNext21() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext21();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmin-xmaxv,eps) == true )
      {
        if ( yminv < ymax - eps && ymaxv > ymin + eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp);
      }
    }
    if ( stackp->getNext24() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext24();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmin-xmaxv,eps) == true )
      {
        if ( yminv < ymax - eps && ymaxv > ymin + eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp);
      }
    }
    if ( stackp->getNext27() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext27();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmin-xmaxv,eps) == true )
      {
        if ( yminv < ymax - eps && ymaxv > ymin + eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp);
      }
    }
  }
  E_Int l0 = -1;
  E_Int ncandidats = candidats.size();

  for (E_Int v = 0; v < ncandidats; v++)
  {
    if ( candidats[v]->getLevel() > l0 ) 
    {l0 = candidats[v]->getLevel(); (*nodep)->setVoisin1(candidats[v]);}
  }
  return *nodep;
}
//=============================================================================
OctreeNode* updateVoisin2_27(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData3> stack;
  OctreeNode* voisin = node->getVoisin2();
  E_Int lvl0 = node->getLevel();
  E_Float xmin = node->getXmin(); E_Float xmax = xmin+node->getDh();
  E_Float ymin = node->getYmin(); E_Float ymax = ymin+node->getDh();
  E_Float zmin = node->getZmin(); E_Float zmax = zmin+node->getDh();
  E_Float xminv, yminv, zminv, ymaxv, zmaxv;
  vector<OctreeNode*> candidats;
  stackData3 dataForStack;

  if ( voisin == NULL ) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if (voisin->getNext1() == NULL ) return *nodep;
  //1-4-7-10-13-16-19-22-25
  if ( voisin->getNext1() == NULL && voisin->getNext4() == NULL && voisin->getNext7() == NULL &&
       voisin->getNext10() == NULL && voisin->getNext13() == NULL && voisin->getNext16() == NULL &&
       voisin->getNext19() == NULL && voisin->getNext22() == NULL && voisin->getNext25() == NULL )
    return *nodep;
 
  if ( voisin->getNext1() != NULL ) 
  {dataForStack.current = voisin->getNext1(); stack.push(dataForStack);}
  if ( voisin->getNext4() != NULL ) 
  {dataForStack.current = voisin->getNext4(); stack.push(dataForStack);}
  if ( voisin->getNext7() != NULL ) 
  {dataForStack.current = voisin->getNext7(); stack.push(dataForStack);}
  if ( voisin->getNext10() != NULL ) 
  {dataForStack.current = voisin->getNext10(); stack.push(dataForStack);}
  if ( voisin->getNext13() != NULL ) 
  {dataForStack.current = voisin->getNext13(); stack.push(dataForStack);}
  if ( voisin->getNext16() != NULL ) 
  {dataForStack.current = voisin->getNext16(); stack.push(dataForStack);}
  if ( voisin->getNext19() != NULL ) 
  {dataForStack.current = voisin->getNext19(); stack.push(dataForStack);}
  if ( voisin->getNext22() != NULL ) 
  {dataForStack.current = voisin->getNext22(); stack.push(dataForStack);}
  if ( voisin->getNext25() != NULL ) 
  {dataForStack.current = voisin->getNext25(); stack.push(dataForStack);}
  stackData3 s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    xminv = stackp->getXmin(); //xmaxv = xminv+stackp->getDh();
    yminv = stackp->getYmin(); ymaxv = yminv+stackp->getDh();
    zminv = stackp->getZmin(); zmaxv = zminv+stackp->getDh();
    if ( stackp->getNext1() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext1();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmax-xminv,eps) == true )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext4() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext4();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmax-xminv,eps) == true )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }

    if ( stackp->getNext7() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext7();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmax-xminv,eps) == true )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }

    if ( stackp->getNext10() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext10();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmax-xminv,eps) == true )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }

    if ( stackp->getNext13() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext13();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmax-xminv,eps) == true )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }

    if ( stackp->getNext16() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext16();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmax-xminv,eps) == true )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }

    if ( stackp->getNext19() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext19();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmax-xminv,eps) == true )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }

    if ( stackp->getNext22() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext22();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmax-xminv,eps) == true )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }

    if ( stackp->getNext25() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext25();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmax-xminv,eps) == true )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }

  }
  E_Int l0 = -1;
  E_Int ncandidats = candidats.size();
  for (E_Int v = 0; v < ncandidats; v++)
  {
    if ( candidats[v]->getLevel() > l0 ) 
    {l0 = candidats[v]->getLevel(); (*nodep)->setVoisin2(candidats[v]);}
  }

  return *nodep;
}
//=============================================================================
OctreeNode* updateVoisin3_27(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData3> stack;
  OctreeNode* voisin = node->getVoisin3();
  E_Int lvl0 = node->getLevel();
  E_Float xmin = node->getXmin(); E_Float xmax = xmin+node->getDh();
  E_Float ymin = node->getYmin(); //E_Float ymax = ymin+node->getDh();
  E_Float zmin = node->getZmin(); E_Float zmax = zmin+node->getDh();
  E_Float xminv, yminv, zminv, xmaxv, ymaxv, zmaxv;
  vector<OctreeNode*> candidats;
  stackData3 dataForStack;

  if ( voisin == NULL ) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if (voisin->getNext1() == NULL ) return *nodep;
  //7-8-9-16-17-18-25-26-27
  if ( voisin->getNext7() == NULL && voisin->getNext8() == NULL && voisin->getNext9() == NULL &&
       voisin->getNext16() == NULL && voisin->getNext17() == NULL && voisin->getNext18() == NULL &&
       voisin->getNext25() == NULL && voisin->getNext26() == NULL && voisin->getNext27() == NULL )
    return *nodep;

  if ( voisin->getNext7() != NULL ) 
  {dataForStack.current = voisin->getNext3(); stack.push(dataForStack);}
  if ( voisin->getNext8() != NULL ) 
  {dataForStack.current = voisin->getNext6(); stack.push(dataForStack);}
  if ( voisin->getNext9() != NULL ) 
  {dataForStack.current = voisin->getNext9(); stack.push(dataForStack);}
  if ( voisin->getNext16() != NULL ) 
  {dataForStack.current = voisin->getNext12(); stack.push(dataForStack);}
  if ( voisin->getNext17() != NULL ) 
  {dataForStack.current = voisin->getNext15(); stack.push(dataForStack);}
  if ( voisin->getNext18() != NULL ) 
  {dataForStack.current = voisin->getNext18(); stack.push(dataForStack);}
  if ( voisin->getNext25() != NULL ) 
  {dataForStack.current = voisin->getNext21(); stack.push(dataForStack);}
  if ( voisin->getNext26() != NULL ) 
  {dataForStack.current = voisin->getNext24(); stack.push(dataForStack);}
  if ( voisin->getNext27() != NULL ) 
  {dataForStack.current = voisin->getNext27(); stack.push(dataForStack);}

  stackData3 s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    xminv = stackp->getXmin(); xmaxv = xminv+stackp->getDh();
    yminv = stackp->getYmin(); ymaxv = yminv+stackp->getDh();
    zminv = stackp->getZmin(); zmaxv = zminv+stackp->getDh();

    if ( stackp->getNext7() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext7(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymin-ymaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext8() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext8(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymin-ymaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext9() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext9(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymin-ymaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext16() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext16(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymin-ymaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext17() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext17(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymin-ymaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext18() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext18(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymin-ymaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext25() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext25(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymin-ymaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext26() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext26(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymin-ymaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext27() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext27(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymin-ymaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }
  }
  E_Int l0 = -1;
  E_Int ncandidats = candidats.size();
  for (E_Int v = 0; v < ncandidats; v++)
  {
    if ( candidats[v]->getLevel() > l0 ) 
    {l0 = candidats[v]->getLevel(); (*nodep)->setVoisin3(candidats[v]);}
  }

  return *nodep;
}
//=============================================================================
OctreeNode* updateVoisin4_27(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData3> stack;
  OctreeNode* voisin = node->getVoisin4();
  E_Int lvl0 = node->getLevel();
  E_Float xmin = node->getXmin(); E_Float xmax = xmin+node->getDh();
  E_Float ymin = node->getYmin(); E_Float ymax = ymin+node->getDh();
  E_Float zmin = node->getZmin(); E_Float zmax = zmin+node->getDh();
  E_Float xminv, yminv, zminv, xmaxv, zmaxv;
  vector<OctreeNode*> candidats;
  stackData3 dataForStack;

  if ( voisin == NULL ) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if (voisin->getNext1() == NULL ) return *nodep;
  //1-2-3-10-11-12-19-20-21
  if ( voisin->getNext1() == NULL && voisin->getNext2() == NULL && voisin->getNext3() == NULL &&
       voisin->getNext10() == NULL && voisin->getNext11() == NULL && voisin->getNext12() == NULL &&
       voisin->getNext19() == NULL && voisin->getNext20() == NULL && voisin->getNext21() == NULL )
    return *nodep;

  if ( voisin->getNext1() != NULL ) 
  {dataForStack.current = voisin->getNext1(); stack.push(dataForStack);}
  if ( voisin->getNext2() != NULL ) 
  {dataForStack.current = voisin->getNext2(); stack.push(dataForStack);}
  if ( voisin->getNext3() != NULL ) 
  {dataForStack.current = voisin->getNext3(); stack.push(dataForStack);}
  if ( voisin->getNext10() != NULL ) 
  {dataForStack.current = voisin->getNext10(); stack.push(dataForStack);}
  if ( voisin->getNext11() != NULL ) 
  {dataForStack.current = voisin->getNext11(); stack.push(dataForStack);}
  if ( voisin->getNext12() != NULL ) 
  {dataForStack.current = voisin->getNext12(); stack.push(dataForStack);}
  if ( voisin->getNext19() != NULL ) 
  {dataForStack.current = voisin->getNext19(); stack.push(dataForStack);}
  if ( voisin->getNext20() != NULL ) 
  {dataForStack.current = voisin->getNext20(); stack.push(dataForStack);}
  if ( voisin->getNext21() != NULL ) 
  {dataForStack.current = voisin->getNext21(); stack.push(dataForStack);}

  stackData3 s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;
  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    xminv = stackp->getXmin(); xmaxv = xminv+stackp->getDh();
    yminv = stackp->getYmin(); //ymaxv = yminv+stackp->getDh();
    zminv = stackp->getZmin(); zmaxv = zminv+stackp->getDh();
    if ( stackp->getNext1() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext1(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymax-yminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }     
    }
    if ( stackp->getNext2() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext2(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymax-yminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }     
    }
    if ( stackp->getNext3() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext3();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymax-yminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }    
    }
    if ( stackp->getNext10() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext10();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymax-yminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }    
    }
    if ( stackp->getNext11() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext11();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymax-yminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }    
    }
    if ( stackp->getNext12() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext12();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymax-yminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }    
    }
    if ( stackp->getNext19() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext19();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymax-yminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      { 
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }    
    }
    if ( stackp->getNext20() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext20();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymax-yminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }    
    }
    if ( stackp->getNext21() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext21();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymax-yminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }    
    }
  }

  E_Int l0 = -1;
  E_Int ncandidats = candidats.size();
  for (E_Int v = 0; v < ncandidats; v++)
  {
    if ( candidats[v]->getLevel() > l0 ) 
    {l0 = candidats[v]->getLevel(); (*nodep)->setVoisin4(candidats[v]);}
  }
  return *nodep;
}
//=============================================================================
OctreeNode* updateVoisin5_27(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData3> stack;
  OctreeNode* voisin = node->getVoisin5();
  E_Int lvl0 = node->getLevel();
  E_Float xmin = node->getXmin(); E_Float xmax = xmin+node->getDh();
  E_Float ymin = node->getYmin(); E_Float ymax = ymin+node->getDh();
  E_Float zmin = node->getZmin(); //E_Float zmax = zmin+node->getDh();
  E_Float xminv, yminv, zminv, xmaxv, ymaxv, zmaxv;
  vector<OctreeNode*> candidats;
  stackData3 dataForStack;

  if ( voisin == NULL ) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if (voisin->getNext1() == NULL ) return *nodep;
  //19 a 27
   if ( voisin->getNext19() == NULL && voisin->getNext20() == NULL && voisin->getNext21() == NULL &&
        voisin->getNext22() == NULL && voisin->getNext23() == NULL && voisin->getNext24() == NULL &&
        voisin->getNext25() == NULL && voisin->getNext26() == NULL && voisin->getNext27() == NULL )
    return *nodep;
 
  if ( voisin->getNext19() != NULL ) 
  {dataForStack.current = voisin->getNext19(); stack.push(dataForStack);}
  if ( voisin->getNext20() != NULL ) 
  {dataForStack.current = voisin->getNext20(); stack.push(dataForStack);}
  if ( voisin->getNext21() != NULL ) 
  {dataForStack.current = voisin->getNext21(); stack.push(dataForStack);}
  if ( voisin->getNext22() != NULL ) 
  {dataForStack.current = voisin->getNext22(); stack.push(dataForStack);}
  if ( voisin->getNext23() != NULL ) 
  {dataForStack.current = voisin->getNext23(); stack.push(dataForStack);}
  if ( voisin->getNext24() != NULL ) 
  {dataForStack.current = voisin->getNext24(); stack.push(dataForStack);}
  if ( voisin->getNext25() != NULL ) 
  {dataForStack.current = voisin->getNext25(); stack.push(dataForStack);}
  if ( voisin->getNext26() != NULL ) 
  {dataForStack.current = voisin->getNext26(); stack.push(dataForStack);}
  if ( voisin->getNext27() != NULL ) 
  {dataForStack.current = voisin->getNext27(); stack.push(dataForStack);}

  stackData3 s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;
  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    xminv = stackp->getXmin(); xmaxv = xminv+stackp->getDh();
    yminv = stackp->getYmin(); ymaxv = yminv+stackp->getDh();
    zminv = stackp->getZmin(); zmaxv = zminv+stackp->getDh();
    if ( stackp->getNext19() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext19(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(zmin-zmaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext20() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext20(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(zmin-zmaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      }

    }
    if ( stackp->getNext21() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext21();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(zmin-zmaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext22() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext22();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(zmin-zmaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      }
    }
   if ( stackp->getNext23() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext23(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(zmin-zmaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      }
    }
   if ( stackp->getNext24() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext24(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(zmin-zmaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      }
    }
   if ( stackp->getNext25() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext25(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(zmin-zmaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      }
    }
   if ( stackp->getNext26() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext26(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(zmin-zmaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      }
    }
   if ( stackp->getNext27() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext27(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(zmin-zmaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      }
    }
  }
 
  E_Int l0 = -1;
  E_Int ncandidats = candidats.size();
  for (E_Int v = 0; v < ncandidats; v++)
  {
    if ( candidats[v]->getLevel() > l0 ) 
    {l0 = candidats[v]->getLevel(); (*nodep)->setVoisin5(candidats[v]);}
  }

  return *nodep;
}
//=============================================================================
OctreeNode* updateVoisin6_27(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData3> stack;
  OctreeNode* voisin = node->getVoisin6();
  E_Int lvl0 = node->getLevel();
  E_Float xmin = node->getXmin(); E_Float xmax = xmin+node->getDh();
  E_Float ymin = node->getYmin(); E_Float ymax = ymin+node->getDh();
  E_Float zmin = node->getZmin(); E_Float zmax = zmin+node->getDh();
  E_Float xminv, yminv, zminv, xmaxv, ymaxv;
  vector<OctreeNode*> candidats;
  stackData3 dataForStack;

  if ( voisin == NULL ) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if (voisin->getNext1() == NULL ) return *nodep;
  //1 a 9
  if ( voisin->getNext1() == NULL && voisin->getNext2() == NULL && voisin->getNext3() == NULL &&
       voisin->getNext4() == NULL && voisin->getNext5() == NULL && voisin->getNext6() == NULL &&
       voisin->getNext7() == NULL && voisin->getNext8() == NULL && voisin->getNext9() == NULL )
    return *nodep;

  if ( voisin->getNext1() != NULL ) 
  {dataForStack.current = voisin->getNext1(); stack.push(dataForStack);}
  if ( voisin->getNext2() != NULL ) 
  {dataForStack.current = voisin->getNext2(); stack.push(dataForStack);}
  if ( voisin->getNext3() != NULL ) 
  {dataForStack.current = voisin->getNext3(); stack.push(dataForStack);}
  if ( voisin->getNext4() != NULL ) 
  {dataForStack.current = voisin->getNext4(); stack.push(dataForStack);}
  if ( voisin->getNext5() != NULL ) 
  {dataForStack.current = voisin->getNext5(); stack.push(dataForStack);}
  if ( voisin->getNext6() != NULL ) 
  {dataForStack.current = voisin->getNext6(); stack.push(dataForStack);}
  if ( voisin->getNext7() != NULL ) 
  {dataForStack.current = voisin->getNext7(); stack.push(dataForStack);}
  if ( voisin->getNext8() != NULL ) 
  {dataForStack.current = voisin->getNext8(); stack.push(dataForStack);}
  if ( voisin->getNext9() != NULL ) 
  {dataForStack.current = voisin->getNext9(); stack.push(dataForStack);}

  stackData3 s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;
  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    xminv = stackp->getXmin(); xmaxv = xminv+stackp->getDh();
    yminv = stackp->getYmin(); ymaxv = yminv+stackp->getDh();
    zminv = stackp->getZmin(); //zmaxv = zminv+stackp->getDh();
    if ( stackp->getNext1() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext1(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(zmax-zminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      } 
    }
    if ( stackp->getNext2() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext2(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(zmax-zminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      } 
    }
    if ( stackp->getNext3() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext3();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(zmax-zminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      } 
    }
    if ( stackp->getNext4() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext4();stack.push(dataForStack);}
    }
    else 
    { 
      if( K_FUNC::fEqualZero(zmax-zminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      } 
    }
    if ( stackp->getNext5() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext5();stack.push(dataForStack);}
    }
    else 
    { 
      if( K_FUNC::fEqualZero(zmax-zminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      } 
    }
    if ( stackp->getNext6() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext6();stack.push(dataForStack);}
    }
    else 
    { 
      if( K_FUNC::fEqualZero(zmax-zminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      } 
    }
    if ( stackp->getNext7() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext7();stack.push(dataForStack);}
    }
    else 
    { 
      if( K_FUNC::fEqualZero(zmax-zminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      } 
    }
    if ( stackp->getNext8() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext8();stack.push(dataForStack);}
    }
    else 
    { 
      if( K_FUNC::fEqualZero(zmax-zminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      } 
    }
    if ( stackp->getNext9() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext9();stack.push(dataForStack);}
    }
    else 
    { 
      if( K_FUNC::fEqualZero(zmax-zminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      } 
    }
  }

  E_Int l0 = -1;
  E_Int ncandidats = candidats.size();
  for (E_Int v = 0; v < ncandidats; v++)
  {
    if ( candidats[v]->getLevel() > l0 ) 
    {l0 = candidats[v]->getLevel(); (*nodep)->setVoisin6(candidats[v]);}
  }

  return *nodep;
}
}
//=============================================================================
