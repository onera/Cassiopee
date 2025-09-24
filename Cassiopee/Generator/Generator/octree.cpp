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
# include <stack>
# include "Nuga/include/BbTree.h"
# include "Search/OctreeNode.h"
# include "Nuga/include/ArrayAccessor.h"
using namespace std;
using namespace K_FLD;
using namespace K_SEARCH;

namespace K_GENERATOR 
{
struct stackData
{
    OctreeNode* current;
};

  OctreeNode* splitNode(OctreeNode* current, E_Int levelMax, E_Int dim, 
                        E_Int& split, E_Int& indmax);
  OctreeNode* splitNode2D(OctreeNode* current, E_Int levelMax, E_Int& split, 
                          E_Int& indmax);
  OctreeNode* splitVoisin(OctreeNode* node, E_Int levelMax, E_Int dim, 
                          stack<stackData>& stack, E_Int& indmax);
  OctreeNode* splitVoisin2D(OctreeNode* node, E_Int levelMax, 
                            stack<stackData>& stack, E_Int& indmax);
  OctreeNode* splitVoisin3D(OctreeNode* node, E_Int levelMax, 
                            stack<stackData>& stack, E_Int& indmax);
  OctreeNode* addSplitVoisin2D(OctreeNode* voisin, E_Int l0, E_Int levelMax, 
                               stack<stackData>& stack, E_Int& indmax);
  OctreeNode* addSplitVoisin3D(OctreeNode* voisin, E_Int l0, E_Int levelMax, 
                               stack<stackData>& stack, E_Int& indmax);
  
  OctreeNode* updateVoisin2D1(OctreeNode* node);
  OctreeNode* updateVoisin2D2(OctreeNode* node);
  OctreeNode* updateVoisin2D3(OctreeNode* node);
  OctreeNode* updateVoisin2D4(OctreeNode* node);
  OctreeNode* updateVoisin2D5(OctreeNode* node);
  OctreeNode* updateVoisin2D6(OctreeNode* node);
  OctreeNode* updateVoisin1(OctreeNode* node);
  OctreeNode* updateVoisin2(OctreeNode* node);
  OctreeNode* updateVoisin3(OctreeNode* node);
  OctreeNode* updateVoisin4(OctreeNode* node);
  OctreeNode* updateVoisin5(OctreeNode* node);
  OctreeNode* updateVoisin6(OctreeNode* node);

//============================================================================
/* Generation d un octree a partir d une liste de surfaces et de snear      */
//============================================================================
PyObject* octree(PyObject* self, PyObject* args)
{
  E_Float __DFARTOL__ = -0.5;// dfar=-1 par defaut = pas pris en compte
  PyObject* stlArrays; PyObject* listOfSnears;
  PyObject* listOfDfars;
  E_Float dfar; E_Int levelMax; E_Int dfarDir; E_Int mode;
  PyObject* octant;
  if (!PYPARSETUPLE_(args, OOO_ R_ I_ O_ I_ I_,
                    &stlArrays, &listOfSnears, &listOfDfars,
                    &dfar, &levelMax, &octant, &dfarDir, &mode)) return NULL;
  
  if (PyList_Size(stlArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "octree: 1st argument is an empty list.");
    return NULL;
  }
  if (PyList_Size(listOfSnears) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "octree: 2nd argument is an empty list.");
    return NULL;
  }
  if (PyList_Size(stlArrays) != PyList_Size(listOfSnears)) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "octree: 1st and 2nd args must be of same length.");
    return NULL;
  }
  E_Int ndfars = PyList_Size(listOfDfars);
  if (ndfars == 0 && dfar < __DFARTOL__)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "octree: you must set a global dfar as a positive value or define a list of dfars.");
    return NULL;
  }
  if (ndfars > 0 && dfar > __DFARTOL__)
  {
    printf("octree: both dfar and listOfDfars are defined; the list of dfars is taken into account.\n");
  }

  // recuperations des stl
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt; vector<char*> eltTypet;
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
    PyErr_SetString(PyExc_TypeError, 
                    "octree: 1st arg is not valid.");
    return NULL;
  }
  E_Int nzones = unstrF.size();

  E_Int dim = -1;
  for (E_Int i = 0; i < nzones; i++)
  {
    if (strcmp(eltTypet[i], "TRI") == 0) 
    {
      if (dim == -1) dim = 3;
      else if (dim != 3) 
      {
        for (E_Int no = 0; no < nzones; no++)
          RELEASESHAREDU(objut[no], unstrF[no], cnt[no]);
        PyErr_SetString(PyExc_TypeError, 
                        "octree: 1st arg must be a list of TRI zones.");
        return NULL;
      }
    }
    else if (strcmp(eltTypet[i], "BAR") == 0) 
    {
      if (dim == -1) dim = 2;
      else if (dim != 2)
      {
        for (E_Int no = 0; no < nzones; no++)
          RELEASESHAREDU(objut[no], unstrF[no], cnt[no]);
        PyErr_SetString(PyExc_TypeError, 
                        "octree: 1st arg must be a list of BAR zones.");
        return NULL;
      }
    }
    else 
    {
      for (E_Int no = 0; no < nzones; no++)
        RELEASESHAREDU(objut[no], unstrF[no], cnt[no]);
      PyErr_SetString(PyExc_TypeError, 
                      "octree: 1st arg must be a list of TRI or BAR zones.");
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
  if (nzones != nsnear)
  {
    for (E_Int no = 0; no < nzones; no++)
      RELEASESHAREDU(objut[no], unstrF[no], cnt[no]);
    PyErr_SetString(PyExc_TypeError, 
                    "octree: 1st and 2nd args must be consistent.");
    return NULL;
  }    
  // recuperation des dfars
  if (ndfars > 0 && ndfars != nzones)
  {
    for (E_Int no = 0; no < nzones; no++)
      RELEASESHAREDU(objut[no], unstrF[no], cnt[no]);
    PyErr_SetString(PyExc_TypeError, 
                    "octree: size of Dfarlist must be equal to the size of 1st arg (list of surface zones).");
    return NULL;
  }
  vector<E_Float> vectOfDfars;
  if (ndfars > 0)
  {
    PyObject* tpl0 = NULL;

    for (E_Int i = 0; i < ndfars; i++)
    {
      E_Float dfarloc;
      tpl0 = PyList_GetItem(listOfDfars,i);
      if (PyFloat_Check(tpl0) == 0)
      {
        for (E_Int no = 0; no < nzones; no++)
          RELEASESHAREDU(objut[no], unstrF[no], cnt[no]);
        PyErr_SetString(PyExc_TypeError, 
                        "octree: dfarList must be a list of floats.");
        return NULL;
      }
      else dfarloc = PyFloat_AsDouble(tpl0);
      vectOfDfars.push_back(dfarloc);
    }
  }

  PyObject* tpl = NULL;
  vector<E_Float> snears(nzones);
  for (int i = 0; i < nzones; i++)
  {
    tpl = PyList_GetItem(listOfSnears, i);
    if (PyFloat_Check(tpl) == 0)
    {
      for (E_Int no = 0; no < nzones; no++) RELEASESHAREDU(objut[no], unstrF[no], cnt[no]);
      PyErr_SetString(PyExc_TypeError,
                      "octree: not a valid value for snear.");
      return NULL;
    }
    else snears[i] = PyFloat_AsDouble(tpl);
  }
  
  typedef K_SEARCH::BoundingBox<3> BBox3DType;
  vector<K_SEARCH::BbTree3D*> bboxtrees(nzones);
  vector< vector<BBox3DType*> > vectOfBBoxes;// pour etre detruit a la fin
  E_Float minB[3]; E_Float maxB[3];
  E_Int posx2, posy2, posz2, nelts2;
  vector<E_Float> xminZ(nzones);
  vector<E_Float> yminZ(nzones);
  vector<E_Float> zminZ(nzones);
  vector<E_Float> xmaxZ(nzones);
  vector<E_Float> ymaxZ(nzones);
  vector<E_Float> zmaxZ(nzones);
  for (E_Int v = 0; v < nzones; v++)
  {
    FldArrayF& f2 = *unstrF[v]; FldArrayI& cn2 = *cnt[v];
    posx2 = posxt[v]; posy2 = posyt[v]; posz2 = poszt[v];

    // Creation de la bboxtree
    nelts2 = cn2.getSize();
    vector<BBox3DType*> boxes(nelts2);// liste des bbox de ts les elements de a2
    K_FLD::FldArrayF bbox(nelts2,6);// xmin, ymin, zmin, xmax, ymax, zmax
    K_COMPGEOM::boundingBoxOfUnstrCells(
      cn2, f2.begin(posx2), f2.begin(posy2),f2.begin(posz2),bbox);
    E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
    E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
    E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
    E_Float xminl = K_CONST::E_MAX_FLOAT; E_Float xmaxl =-K_CONST::E_MAX_FLOAT;
    E_Float yminl = K_CONST::E_MAX_FLOAT; E_Float ymaxl =-K_CONST::E_MAX_FLOAT;
    E_Float zminl = K_CONST::E_MAX_FLOAT; E_Float zmaxl =-K_CONST::E_MAX_FLOAT;
    for (E_Int et = 0; et < nelts2; et++)
    {
      minB[0] = xminp[et]; minB[1] = yminp[et]; minB[2] = zminp[et];
      maxB[0] = xmaxp[et]; maxB[1] = ymaxp[et]; maxB[2] = zmaxp[et]; 
      boxes[et] = new BBox3DType(minB, maxB);
      xminl = K_FUNC::E_min(minB[0],xminl); xmaxl = K_FUNC::E_max(maxB[0],xmaxl);
      yminl = K_FUNC::E_min(minB[1],yminl); ymaxl = K_FUNC::E_max(maxB[1],ymaxl);
      if (dim == 3) {zminl = K_FUNC::E_min(minB[2],zminl); zmaxl = K_FUNC::E_max(maxB[2],zmaxl);}
      else {zminl= 0.; zmaxl=0.;}
    }
    // Build the box tree
    bboxtrees[v] = new K_SEARCH::BbTree3D(boxes);
    vectOfBBoxes.push_back(boxes);
    // bbox de la surface noz 
    xminZ[v]=xminl; yminZ[v]=yminl; zminZ[v]=zminl;  
    xmaxZ[v]=xmaxl; ymaxZ[v]=ymaxl; zmaxZ[v]=zmaxl;
  }
  for (E_Int no = 0; no < nzones; no++)
    RELEASESHAREDU(objut[no], unstrF[no], cnt[no]);

  // calcul du dfar reel a partir de la bbox des bbox
  E_Float xmino, ymino, zmino, xmaxo, ymaxo, zmaxo;
  E_Float xc, yc, zc, Delta, Deltax, Deltay, Deltaz;
  if (ndfars == 0) // starts from the global dfar
  {
    E_Float dfxm = dfar; E_Float dfym = dfar; E_Float dfzm = dfar;
    E_Float dfxp = dfar; E_Float dfyp = dfar; E_Float dfzp = dfar;  
    if (dim == 2) { dfzp = 0.; dfzm = 0.; }

    xmino =  K_CONST::E_MAX_FLOAT; xmaxo = -K_CONST::E_MAX_FLOAT;
    ymino =  K_CONST::E_MAX_FLOAT; ymaxo = -K_CONST::E_MAX_FLOAT;
    zmino =  K_CONST::E_MAX_FLOAT; zmaxo = -K_CONST::E_MAX_FLOAT;
    if (dim == 2) {zmino = 0.; zmaxo= 0.;}
    for (E_Int v = 0; v < nzones; v++)
    {
      xmino = K_FUNC::E_min(xminZ[v],xmino); xmaxo = K_FUNC::E_max(xmaxZ[v],xmaxo);
      ymino = K_FUNC::E_min(yminZ[v],ymino); ymaxo = K_FUNC::E_max(ymaxZ[v],ymaxo);
      if (dim == 3)
      {
        zmino = K_FUNC::E_min(zminZ[v],zmino); zmaxo = K_FUNC::E_max(zmaxZ[v],zmaxo);
      }
    }
    if (dim == 2) {zmino = 0.; zmaxo = 0.;}

    Deltax = xmaxo+dfxp-xmino+dfxm;
    Deltay = ymaxo+dfyp-ymino+dfym;
    Deltaz = zmaxo+dfzp-zmino+dfzm;
    Delta = K_FUNC::E_max(Deltax, Deltay); 
    if (dim == 3) Delta = K_FUNC::E_max(Delta, Deltaz);
    Delta = 0.5*Delta;
    xc = 0.5*(xmaxo+dfxp+xmino-dfxm);
    yc = 0.5*(ymaxo+dfyp+ymino-dfym);
    zc = 0.5*(zmaxo+dfzp+zmino-dfzm);
    zmino = 0.; zmaxo = 0.;
    xmino = xc-Delta; ymino = yc-Delta; if (dim == 3) zmino = zc-Delta;
    xmaxo = xc+Delta; ymaxo = yc+Delta; if (dim == 3) zmaxo = zc+Delta;
  }
  else // if (ndfars > 0)// local dfar
  {
    // calcul de la bbox etendue de dfar local pour ttes les surfaces sauf celles tq dfarloc=-1
    xmino = K_CONST::E_MAX_FLOAT; xmaxo =-K_CONST::E_MAX_FLOAT;
    ymino = K_CONST::E_MAX_FLOAT; ymaxo =-K_CONST::E_MAX_FLOAT;
    zmino = K_CONST::E_MAX_FLOAT; zmaxo =-K_CONST::E_MAX_FLOAT;
    if (dim == 2) {zmino = 0.; zmaxo= 0.;}

    for (E_Int v = 0; v < nzones; v++)
    {
      E_Float dfarloc = vectOfDfars[v];
      if (dfarloc > __DFARTOL__) // extensions non isotropes
      {
        E_Float xminl = xminZ[v]; E_Float xmaxl = xmaxZ[v];
        E_Float yminl = yminZ[v]; E_Float ymaxl = ymaxZ[v];
        E_Float zminl = zminZ[v]; E_Float zmaxl = zmaxZ[v];

        Deltax = xmaxl-xminl+2*dfarloc;
        Deltay = ymaxl-yminl+2*dfarloc;
        Delta = K_FUNC::E_max(Deltax, Deltay); 
        if (dim == 3) 
        {
          Deltaz=zmaxl-zminl+2*dfarloc;
          Delta = K_FUNC::E_max(Delta, Deltaz);
        }
        Delta = 0.5*Delta;
        xc = 0.5*(xmaxl+xminl);
        yc = 0.5*(ymaxl+yminl);
        xminl=xc-Delta; yminl=yc-Delta;
        xmaxl=xc+Delta; ymaxl=yc+Delta; 

        if (dim == 2)
        {zc = 0.; zminl = 0.; zmaxl = 0.;}
        else 
        {
          zc = 0.5*(zmaxl+zminl);
          zminl=zc-Delta; zmaxl=zc+Delta;
        }
        xmino = K_FUNC::E_min(xmino,xminl); xmaxo = K_FUNC::E_max(xmaxo,xmaxl);
        ymino = K_FUNC::E_min(ymino,yminl); ymaxo = K_FUNC::E_max(ymaxo,ymaxl);
        if (dim == 3)
        {
          zmino = K_FUNC::E_min(zmino,zminl); zmaxo = K_FUNC::E_max(zmaxo,zmaxl);
        }        
      }//only zones with dfar > -1
    }// loop on all zones
  }// list of dfars

  // Nous rendons cubique la boite finale
  E_Float d = 0.;
  if (dfarDir == 0) // impose dfar max
  { 
    d = K_FUNC::E_max(xmaxo-xmino,ymaxo-ymino);
    if (dim == 3) d = K_FUNC::E_max(d, zmaxo-zmino);
  }
  else if (dfarDir == 1)
    d = xmaxo-xmino;
  else if (dfarDir == 2)
    d = ymaxo-ymino;
  else if (dfarDir == 3)
    d = zmaxo-zmino;

  if (mode == 1)
  {
    printf("size domain before octree size adaptation = %f\n",d);
    E_Float dd = d;
    E_Float snear_target = 1.e6;
    for (E_Int v = 0; v < nzones; v++) snear_target = K_FUNC::E_min(snears[v], snear_target);
    while (dd > snear_target) dd = dd/2.;
    d = d*snear_target/(2*dd);
    d = d*2; // better having a bigger domain than a lower one...
    printf("size domain after octree size adaptation = %f\n",d);
  }
  
  xc = 0.5*(xmino+xmaxo);
  yc = 0.5*(ymino+ymaxo);
  xmino = xc-0.5*d; ymino = yc-0.5*d;
  xmaxo = xc+0.5*d; ymaxo = yc+0.5*d;
  if (dim == 3) { zc = 0.5*(zmino+zmaxo); zmino = zc-0.5*d; zmaxo = zc+0.5*d; }
    
  // construction de l'octree (octant)
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
  E_Int l0 = 0; // niveau le plus grossier
  E_Int indmax = 0;
  
  OctreeNode* toptree;
  if (dim == 2)
  {
    toptree = new OctreeNode(xmino, ymino, zmino, xmaxo-xmino, l0, indmax, indmax+1, indmax+2, indmax+3);
    indmax += 4;
  }
  else
  {
    toptree = new OctreeNode(xmino, ymino, zmino, xmaxo-xmino, l0, indmax, indmax+1, indmax+2, indmax+3, indmax+4, indmax+5, indmax+6, indmax+7);
    indmax += 8;
  }
  vector<E_Int> indicesBB;
  E_Float snear, dh;
  E_Float tol = 1.e-6;
  E_Int found;
  E_Int nelts = 1;
  OctreeNode* current = toptree;
  stack<stackData> stack;
  stackData dataForStack;
  dataForStack.current = current;
  stack.push(dataForStack);
  //E_Int incn = 19; if (dim == 2) incn = 4;
  E_Int ince = 7; if (dim == 2) ince = 3;
  E_Int split = 0;
  while (stack.size() != 0) 
  {
    dataForStack = stack.top(); stack.pop();
    current = dataForStack.current;
    dh = current->getDh();
    minB[0] = current->getXmin(); maxB[0] = minB[0]+dh;
    minB[1] = current->getYmin(); maxB[1] = minB[1]+dh; 
    minB[2] = current->getZmin(); maxB[2] = minB[2]+dh; 
    snear = K_CONST::E_MAX_FLOAT; found = 0;
    
    for (E_Int v = 0; v < nzones; v++)
    {
      bboxtrees[v]->getOverlappingBoxes(minB, maxB, indicesBB);
      if (indicesBB.size() != 0)
      {
        //regarder si dh > snear ? 
        snear = K_FUNC::E_min(snears[v], snear);
        if (dh > snear-tol && dh != snear) found = 1;
      }
      indicesBB.clear();
    }

    l0 = current->getLevel();
    if (found == 1 && l0 < levelMax) //decouper la cellule en 8 si dh > snear
    {
      current = splitNode(current, levelMax, dim, split, indmax);
      if (split != 1) goto next;

      // fils1->voisins
      dataForStack.current = current->getNext1(); stack.push(dataForStack);
     
      // fils2->voisins
      dataForStack.current = current->getNext2(); stack.push(dataForStack);
      
      // fils3->voisins
      dataForStack.current = current->getNext3(); stack.push(dataForStack);
      
      // fils4->voisins
      dataForStack.current = current->getNext4(); stack.push(dataForStack);

      if (dim == 3) 
      {
        dataForStack.current = current->getNext5(); stack.push(dataForStack);
        dataForStack.current = current->getNext6(); stack.push(dataForStack);
        dataForStack.current = current->getNext7(); stack.push(dataForStack);
        dataForStack.current = current->getNext8(); stack.push(dataForStack);
      }
      nelts += ince;
      next:;
    }
  }
  //nettoyages...
  for (E_Int v = 0; v < nzones; v++)
  {
    delete bboxtrees[v]; 
    E_Int nbboxes =  vectOfBBoxes[v].size();
    for (E_Int v2 = 0; v2 < nbboxes; v2++) delete vectOfBBoxes[v][v2];
  }
  vectOfBBoxes.clear();
  // Construction du maillage octree
  FldArrayI indir(indmax);// nb de pts uniques dans l octree
  indir.setAllValuesAt(-1);
  E_Int* indirp = indir.begin();
  E_Int size1 = indmax; E_Int size2 = nelts;
  FldArrayF* coords = new FldArrayF(size1,3); coords->setAllValuesAt(K_CONST::E_MAX_FLOAT);
  E_Int nvert = 8; if (dim == 2) nvert = 4;
  FldArrayI* cn = new FldArrayI(size2,nvert); // HEXA ou QUAD
  E_Float* xt = coords->begin(1);
  E_Float* yt = coords->begin(2);
  E_Float* zt = coords->begin(3);
  E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2);
  E_Int* cn3 = cn->begin(3); E_Int* cn4 = cn->begin(4);
  E_Int* cn5 = NULL; E_Int* cn6 = NULL;
  E_Int* cn7 = NULL; E_Int* cn8 = NULL;
  
  if (dim == 3)
  {cn5 = cn->begin(5); cn6 = cn->begin(6); cn7 = cn->begin(7); cn8 = cn->begin(8);}

  dataForStack.current = toptree; stack.push(dataForStack);
  E_Int et = 0;
  E_Int ind = 0;
  E_Int ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8;

  while (stack.size() != 0) 
  {
    dataForStack = stack.top(); stack.pop();
    current = dataForStack.current;
    dh = current->getDh();
    xmino = current->getXmin(); xmaxo = xmino+dh;
    ymino = current->getYmin(); ymaxo = ymino+dh;
    zmino = current->getZmin(); zmaxo = zmino+dh;
    ind1 = current->getInd1(); ind2 = current->getInd2();
    ind3 = current->getInd3(); ind4 = current->getInd4();
    ind5 = current->getInd5(); ind6 = current->getInd6();
    ind7 = current->getInd7(); ind8 = current->getInd8();

    // feuille->inserer dans le maillage
    if (current->getNext1() == NULL) 
    { 
      //creation des 8 sommets de la cellule 1
      if (indirp[ind1] == -1)
      {
        xt[ind] = xmino; yt[ind] = ymino; zt[ind] = zmino; ind++;
        cn1[et] = ind; indirp[ind1] = ind; 
      }
      else cn1[et] = indirp[ind1]; 

      if (indirp[ind2] == -1)
      {
        xt[ind] = xmaxo; yt[ind] = ymino; zt[ind] = zmino; ind++;
        cn2[et] = ind; indirp[ind2] = ind; 
      }
      else cn2[et] = indirp[ind2]; 

      if (indirp[ind3] == -1)
      {
        xt[ind] = xmaxo; yt[ind] = ymaxo; zt[ind] = zmino; ind++;
        cn3[et] = ind; indirp[ind3] = ind; 
      }
      else cn3[et] = indirp[ind3];

      if (indirp[ind4] == -1)
      {
        xt[ind] = xmino; yt[ind] = ymaxo; zt[ind] = zmino; ind++;
        cn4[et] = ind; indirp[ind4] = ind; 
      }
      else cn4[et] = indirp[ind4];

      if (dim == 3) 
      {
        if (indirp[ind5] == -1)
        {
          xt[ind] = xmino; yt[ind] = ymino; zt[ind] = zmaxo; ind++;
          cn5[et] = ind; indirp[ind5] = ind; 
        }
        else cn5[et] = indirp[ind5];

        if (indirp[ind6] == -1)
        {
          xt[ind] = xmaxo; yt[ind] = ymino; zt[ind] = zmaxo; ind++;
          cn6[et] = ind; indirp[ind6] = ind;
        }
        else cn6[et] = indirp[ind6];

        if (indirp[ind7] == -1)
        {
          xt[ind] = xmaxo; yt[ind] = ymaxo; zt[ind] = zmaxo; ind++;
          cn7[et] = ind; indirp[ind7] = ind;
        }
        else cn7[et] = indirp[ind7];

        if (indirp[ind8] == -1)
        {
          xt[ind] = xmino; yt[ind] = ymaxo; zt[ind] = zmaxo; ind++;
          cn8[et] = ind; indirp[ind8] = ind;
        }
        else cn8[et] = indirp[ind8];
      }
      et++;

      if (ind + nvert >= size1) 
      { 
        size1 += indmax;  coords->reAllocMat(size1,3); 
        xt = coords->begin(1); yt = coords->begin(2); zt = coords->begin(3);
      }
      if (et + 10 >= size2) 
      { 
        size2 += nelts; cn->reAllocMat(size2,nvert); 
        cn1 = cn->begin(1); cn2 = cn->begin(2); cn3 = cn->begin(3); cn4 = cn->begin(4);
        if (dim == 3) {cn5 = cn->begin(5); cn6 = cn->begin(6); cn7 = cn->begin(7); cn8 = cn->begin(8);}
      }
    }
    else
    {
      dataForStack.current = current->getNext1(); stack.push(dataForStack);
      dataForStack.current = current->getNext2(); stack.push(dataForStack);
      dataForStack.current = current->getNext3(); stack.push(dataForStack);
      dataForStack.current = current->getNext4(); stack.push(dataForStack);
      if (dim == 3) 
      {
        dataForStack.current = current->getNext5(); stack.push(dataForStack);
        dataForStack.current = current->getNext6(); stack.push(dataForStack);
        dataForStack.current = current->getNext7(); stack.push(dataForStack);
        dataForStack.current = current->getNext8(); stack.push(dataForStack);
      }
    }
    delete current;
  }
  toptree = NULL;
  coords->reAllocMat(ind, 3); cn->reAllocMat(et, nvert);
  const char* eltType = "HEXA"; if (dim == 2) eltType = "QUAD"; 
  K_CONNECT::cleanConnectivity(1, 2, 3, 1.e-6, eltType, *coords, *cn);
  //buildArray
  tpl = K_ARRAY::buildArray(*coords, "x,y,z", *cn, -1, eltType, false);

  //nettoyage
  delete coords; delete cn;
  return tpl;
}
//=============================================================================
OctreeNode* splitVoisin(OctreeNode* node, E_Int levelMax, E_Int dim, 
                        stack<stackData>& stack, E_Int& indmax)
{
  if (dim == 2) return splitVoisin2D(node, levelMax, stack, indmax);
  else return splitVoisin3D(node, levelMax, stack, indmax);
}
//=============================================================================
OctreeNode* splitVoisin3D(OctreeNode* node, E_Int levelMax, 
                          stack<stackData>& stack, E_Int& indmax)
{
  OctreeNode** nodep = &node;
  E_Int l0 = (*nodep)->getLevel();
  OctreeNode* voisin = NULL; OctreeNode** voisinp = NULL;
 
  //voisin1 de v1
  *nodep = updateVoisin1(*nodep);
  voisin = (*nodep)->getVoisin1(); voisinp = &voisin;
  *voisinp = addSplitVoisin3D(*voisinp, l0, levelMax, stack, indmax);
  (*nodep)->setVoisin1(*voisinp);
//   *nodep = updateVoisin1(*nodep);

  //voisin2 de v1
  *nodep = updateVoisin2(*nodep);
  voisin = (*nodep)->getVoisin2(); voisinp = &voisin;
  *voisinp = addSplitVoisin3D(*voisinp, l0, levelMax, stack, indmax);
  (*nodep)->setVoisin2(*voisinp);
//   *nodep = updateVoisin2(*nodep);

  //voisin3 de v1
  *nodep = updateVoisin3(*nodep);
  voisin = (*nodep)->getVoisin3(); voisinp = &voisin;
  *voisinp = addSplitVoisin3D(*voisinp, l0, levelMax, stack, indmax);
  (*nodep)->setVoisin3(*voisinp);
//   *nodep = updateVoisin3(*nodep);

  //voisin4 de v1
  *nodep = updateVoisin4(*nodep);
  voisin = (*nodep)->getVoisin4(); voisinp = &voisin;
  *voisinp = addSplitVoisin3D(*voisinp, l0, levelMax, stack, indmax);
  (*nodep)->setVoisin4(*voisinp);
//   *nodep = updateVoisin4(*nodep);

  //voisin5 de v1
  *nodep = updateVoisin5(*nodep);
  voisin = (*nodep)->getVoisin5(); voisinp = &voisin;
  *voisinp = addSplitVoisin3D(*voisinp, l0, levelMax, stack, indmax);
  (*nodep)->setVoisin5(*voisinp);
//   *nodep = updateVoisin5(*nodep);

  //voisin6 de v1
  *nodep = updateVoisin6(*nodep);
  voisin = (*nodep)->getVoisin6(); voisinp = &voisin;
  *voisinp = addSplitVoisin3D(*voisinp, l0, levelMax, stack, indmax);
  (*nodep)->setVoisin6(*voisinp);
//   *nodep = updateVoisin6(*nodep);

  return *nodep;
}
//=============================================================================
OctreeNode* addSplitVoisin3D(OctreeNode* voisin, E_Int l0, E_Int levelMax, 
                             stack<stackData>& stack, E_Int& indmax)
{
  E_Int dim = 3; 
  stackData dataForStack;
  OctreeNode* fils = NULL;
  OctreeNode** filsp = NULL;
  E_Int split = 0; E_Int dl = 0;

  if (voisin != NULL)
  {
    OctreeNode** voisinp = &voisin;
    dl = l0-(*voisinp)->getLevel();
    if (dl > 1)
    { 
      *voisinp = splitNode(*voisinp, levelMax, dim, split, indmax);
      if (split == 1)
      {
        fils = (*voisinp)->getNext1(); filsp = &fils;
        *filsp = splitVoisin(*filsp, levelMax, dim, stack, indmax);
        (*voisinp)->setNext1(*filsp); 
        dataForStack.current = (*voisinp)->getNext1(); stack.push(dataForStack);

        fils = (*voisinp)->getNext2(); filsp = &fils;
        *filsp = splitVoisin(*filsp, levelMax, dim, stack, indmax);
        (*voisinp)->setNext2(*filsp); 
        dataForStack.current = (*voisinp)->getNext2(); stack.push(dataForStack);

        fils = (*voisinp)->getNext3(); filsp = &fils;
        *filsp = splitVoisin(*filsp, levelMax, dim, stack, indmax);
        (*voisinp)->setNext3(*filsp);
        dataForStack.current = (*voisinp)->getNext3(); stack.push(dataForStack);

        fils = (*voisinp)->getNext4(); filsp = &fils;
        *filsp = splitVoisin(*filsp, levelMax, dim, stack, indmax); 
        (*voisinp)->setNext4(*filsp);
        dataForStack.current = (*voisinp)->getNext4(); stack.push(dataForStack);

        fils = (*voisinp)->getNext5(); filsp = &fils;
        *filsp = splitVoisin(*filsp, levelMax, dim, stack, indmax);
        (*voisinp)->setNext5(*filsp); 
        dataForStack.current = (*voisinp)->getNext5(); stack.push(dataForStack);

        fils = (*voisinp)->getNext6(); filsp = &fils;
        *filsp = splitVoisin(*filsp, levelMax, dim, stack, indmax);
        (*voisinp)->setNext6(*filsp); 
        dataForStack.current = (*voisinp)->getNext6(); stack.push(dataForStack);

        fils = (*voisinp)->getNext7(); filsp = &fils;
        *filsp = splitVoisin(*filsp, levelMax, dim, stack, indmax);
        (*voisinp)->setNext7(*filsp);
        dataForStack.current = (*voisinp)->getNext7(); stack.push(dataForStack);

        fils = (*voisinp)->getNext8(); filsp = &fils;
        *filsp = splitVoisin(*filsp, levelMax, dim, stack, indmax); 
        (*voisinp)->setNext8(*filsp);
        dataForStack.current = (*voisinp)->getNext8(); stack.push(dataForStack);

        return *voisinp;
      }
    }//dl > 1
  }
  return voisin;
}
//=============================================================================
OctreeNode* addSplitVoisin2D(OctreeNode* voisin, E_Int l0, E_Int levelMax, 
                             stack<stackData>& stack, E_Int& indmax)
{
  E_Int dim = 2; 
  stackData dataForStack;
  OctreeNode* fils = NULL;
  OctreeNode** filsp = NULL;
  E_Int split = 0;  E_Int dl = 0;
  if (voisin != NULL)
  {
    OctreeNode** voisinp = &voisin;
    dl = l0-(*voisinp)->getLevel();
    if (dl > 1)
    { 
      *voisinp = splitNode(*voisinp, levelMax, dim, split, indmax);
      if ( split == 1 )
      {
        fils = (*voisinp)->getNext1(); filsp = &fils;
        *filsp = splitVoisin(*filsp, levelMax, dim, stack, indmax);
        (*voisinp)->setNext1(*filsp); 
        dataForStack.current = (*voisinp)->getNext1(); stack.push(dataForStack);

        fils = (*voisinp)->getNext2(); filsp = &fils;
        *filsp = splitVoisin(*filsp, levelMax, dim, stack, indmax);
        (*voisinp)->setNext2(*filsp); 
        dataForStack.current = (*voisinp)->getNext2(); stack.push(dataForStack);

        fils = (*voisinp)->getNext3(); filsp = &fils;
        *filsp = splitVoisin(*filsp, levelMax, dim, stack, indmax);
        (*voisinp)->setNext3(*filsp);
        dataForStack.current = (*voisinp)->getNext3(); stack.push(dataForStack);

        fils = (*voisinp)->getNext4(); filsp = &fils;
        *filsp = splitVoisin(*filsp, levelMax, dim, stack, indmax); 
        (*voisinp)->setNext4(*filsp);
        dataForStack.current = (*voisinp)->getNext4(); stack.push(dataForStack);
        return *voisinp;
      }
    }//dl > 1
  }
  return voisin;
}

//=============================================================================
OctreeNode* splitVoisin2D(OctreeNode* node, E_Int levelMax, stack<stackData>& stack, 
                          E_Int& indmax)
{
  OctreeNode** nodep = &node;
  E_Int l0 = (*nodep)->getLevel();
  OctreeNode* voisin = NULL; OctreeNode** voisinp = NULL;
  //voisin1 de v1
  *nodep = updateVoisin2D1(*nodep);
  voisin =  (*nodep)->getVoisin1(); voisinp = &voisin;
  *voisinp = addSplitVoisin2D(*voisinp, l0, levelMax, stack, indmax);
  (*nodep)->setVoisin1(*voisinp);

  //voisin2 de v1
  *nodep = updateVoisin2D2(*nodep);
  voisin =  (*nodep)->getVoisin2(); voisinp = &voisin;
  *voisinp = addSplitVoisin2D(*voisinp, l0, levelMax, stack, indmax);
  (*nodep)->setVoisin2(*voisinp);

  //voisin3 de v1
  *nodep = updateVoisin2D3(*nodep);
  voisin =  (*nodep)->getVoisin3(); voisinp = &voisin;
  *voisinp = addSplitVoisin2D(*voisinp, l0, levelMax, stack, indmax);
  (*nodep)->setVoisin3(*voisinp);

  //voisin4 de v1
  *nodep = updateVoisin2D4(*nodep);
  voisin =  (*nodep)->getVoisin4(); voisinp = &voisin;
  *voisinp = addSplitVoisin2D(*voisinp, l0, levelMax, stack, indmax);
  (*nodep)->setVoisin4(*voisinp);
  return *nodep;
}
//=============================================================================
OctreeNode* updateVoisin2D1(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData> stack;
  OctreeNode* voisin = node->getVoisin1();
  E_Int lvl0 = node->getLevel();
  vector<OctreeNode*> candidats;
  stackData dataForStack;
  E_Float xminv, yminv, xmaxv, ymaxv, dhv;
  E_Float dh = node->getDh(); 
  E_Float xmin = node->getXmin(); E_Float ymin = node->getYmin(); 
  //E_Float xmax = xmin+dh; 
  E_Float ymax = ymin+dh;

  if (voisin == NULL) return *nodep;
  if (lvl0-voisin->getLevel() < 2) return *nodep;
  if (voisin->getNext2() == NULL && voisin->getNext4() == NULL) return *nodep;
  if (voisin->getNext2() != NULL) 
  {
    dataForStack.current = voisin->getNext2(); stack.push(dataForStack);
  }
  if (voisin->getNext4() != NULL) 
  {
    dataForStack.current = voisin->getNext4(); stack.push(dataForStack);
  }
  stackData s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    dhv = stackp->getDh(); 
    xminv = stackp->getXmin(); xmaxv = xminv+dhv;
    yminv = stackp->getYmin(); ymaxv = yminv+dhv;

    //sous-voisin2
    if ( stackp->getNext2() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext2(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmin-xmaxv,eps) == true )
      {
        if ( yminv < ymax - eps && ymaxv > ymin + eps ) candidats.push_back(stackp);
      }
    }
    //sous-voisin4
    if (stackp->getNext4() != NULL)
    {
      if (stackp->getLevel() < lvl0) 
      {dataForStack.current = stackp->getNext4();stack.push(dataForStack);}
    }
    else 
    {
      if (K_FUNC::fEqualZero(xmin-xmaxv,eps) == true)
      {
        if (yminv < ymax - eps && ymaxv > ymin + eps) candidats.push_back(stackp);
      }
    }
  }
  E_Int l0 = -1;

  E_Int ncandidats = candidats.size();
  for (E_Int v = 0; v < ncandidats; v++)
  {
    if (candidats[v]->getLevel() > l0) 
    {l0 = candidats[v]->getLevel(); (*nodep)->setVoisin1(candidats[v]);}
  }

  return *nodep;
}
//=============================================================================
OctreeNode* updateVoisin2D2(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData> stack;
  OctreeNode* voisin = node->getVoisin2();
  E_Int lvl0 = node->getLevel();
  vector<OctreeNode*> candidats;
  stackData dataForStack;
  E_Float xminv, yminv, ymaxv, dhv;
  E_Float dh = node->getDh(); 
  E_Float xmin = node->getXmin(); E_Float ymin = node->getYmin(); 
  E_Float xmax = xmin+dh; E_Float ymax = ymin+dh;

  if (voisin == NULL) return *nodep;
  if (lvl0-voisin->getLevel() < 2) return *nodep;
  if (voisin->getNext1() == NULL && voisin->getNext3() == NULL) return *nodep;
  if (voisin->getNext1() != NULL) 
  {
    dataForStack.current = voisin->getNext1(); stack.push(dataForStack);
  }
  if (voisin->getNext3() != NULL) 
  {
    dataForStack.current = voisin->getNext3(); stack.push(dataForStack);
  }
  stackData s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    dhv = stackp->getDh(); 
    xminv = stackp->getXmin(); //xmaxv = xminv+dhv;
    yminv = stackp->getYmin(); ymaxv = yminv+dhv;
    if (stackp->getNext1() != NULL)
    {
      if (stackp->getLevel() < lvl0) 
      {dataForStack.current = stackp->getNext1(); stack.push(dataForStack);}
    }
    else 
    {
      if (K_FUNC::fEqualZero(xmax-xminv,eps) == true)
      {
        if (yminv < ymax - eps && ymaxv > ymin + eps) candidats.push_back(stackp); 
      }
    }
    if (stackp->getNext3() != NULL)
    {
      if (stackp->getLevel() < lvl0) 
      {dataForStack.current = stackp->getNext3();stack.push(dataForStack);}
    }
    else 
    {
      if (K_FUNC::fEqualZero(xmax-xminv,eps) == true)
      {
        if (yminv < ymax - eps && ymaxv > ymin + eps) candidats.push_back(stackp); 
      }
    }
  }
  E_Int l0 = -1;
  E_Int ncandidats = candidats.size();
  for (E_Int v = 0; v < ncandidats; v++)
  {
    if (candidats[v]->getLevel() > l0) 
    {l0 = candidats[v]->getLevel(); (*nodep)->setVoisin2(candidats[v]);}
  }

  return *nodep;
}
//=============================================================================
OctreeNode* updateVoisin2D3(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData> stack;
  OctreeNode* voisin = node->getVoisin3();//voisin en ymin
  E_Int lvl0 = node->getLevel();
  vector<OctreeNode*> candidats;
  stackData dataForStack;
  E_Float xminv, yminv, xmaxv, ymaxv, dhv;
  E_Float dh = node->getDh(); 
  E_Float xmin = node->getXmin(); E_Float ymin = node->getYmin(); 
  E_Float xmax = xmin+dh; //E_Float ymax = ymin+dh;

  if ( voisin == NULL ) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if ( voisin->getNext3() == NULL && voisin->getNext4() == NULL ) return *nodep;
  if ( voisin->getNext3() != NULL ) 
  {
    dataForStack.current = voisin->getNext3(); stack.push(dataForStack);
  }
  if ( voisin->getNext4() != NULL ) 
  {
    dataForStack.current = voisin->getNext4(); stack.push(dataForStack);
  }
  stackData s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    dhv = stackp->getDh(); 
    xminv = stackp->getXmin(); xmaxv = xminv+dhv;
    yminv = stackp->getYmin(); ymaxv = yminv+dhv;
    if ( stackp->getNext3() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext3();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymin-ymaxv,eps) == true )
      {
        if ( xminv < xmax - eps && xmaxv > xmin + eps ) candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext4() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext4(); stack.push(dataForStack);}
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
OctreeNode* updateVoisin2D4(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData> stack;
  OctreeNode* voisin = node->getVoisin4();//voisin en ymin
  E_Int lvl0 = node->getLevel();
  vector<OctreeNode*> candidats;
  stackData dataForStack;
  E_Float xminv, yminv, xmaxv, dhv;
  E_Float dh = node->getDh(); 
  E_Float xmin = node->getXmin(); E_Float ymin = node->getYmin(); 
  E_Float xmax = xmin+dh; E_Float ymax = ymin+dh;

  if ( voisin == NULL ) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if ( voisin->getNext1() == NULL && voisin->getNext2() == NULL ) return *nodep;
  if ( voisin->getNext1() != NULL ) 
  {
    dataForStack.current = voisin->getNext1(); stack.push(dataForStack);
  }
  if ( voisin->getNext2() != NULL ) 
  {
    dataForStack.current = voisin->getNext2(); stack.push(dataForStack);
  }
  stackData s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    dhv = stackp->getDh(); 
    xminv = stackp->getXmin(); xmaxv = xminv+dhv;
    yminv = stackp->getYmin(); //ymaxv = yminv+dhv;

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
/* Decoupage du noeud courant en 4 - 8 */
//=============================================================================
OctreeNode* splitNode(OctreeNode* current, E_Int levelMax, 
                      E_Int dim, E_Int& split, E_Int& indmax)
{
  if ( dim == 2 ) return splitNode2D(current, levelMax, split, indmax);
  else 
  {
    split = 0;
    if ( current->getLevel() == levelMax ) return current;
    if ( current == NULL ) return current;
    if ( current->getNext1() != NULL ) return current;
    split = 1;
    E_Int lvl = current->getLevel();
    OctreeNode** currentp = &current;
    E_Float dh0 = (*currentp)->getDh();
    E_Float xmin = (*currentp)->getXmin(); //E_Float xmax = xmin+dh0;
    E_Float ymin = (*currentp)->getYmin(); //E_Float ymax = ymin+dh0;  
    E_Float zmin = (*currentp)->getZmin(); //E_Float zmax = zmin+dh0;
    E_Float dh = 0.5*dh0; 
    E_Int ind1 = (*currentp)->getInd1(); E_Int ind2 = (*currentp)->getInd2();
    E_Int ind3 = (*currentp)->getInd3(); E_Int ind4 = (*currentp)->getInd4();
    E_Int ind5 = (*currentp)->getInd5(); E_Int ind6 = (*currentp)->getInd6();
    E_Int ind7 = (*currentp)->getInd7(); E_Int ind8 = (*currentp)->getInd8();
    lvl = lvl + 1;
    OctreeNode* node1 = new OctreeNode(xmin, ymin, zmin, dh, lvl, 
                                       ind1, indmax, indmax+4, indmax+2,
                                       indmax+6, indmax+7, indmax+10, indmax+9);
    OctreeNode* node2 = new OctreeNode(xmin+dh, ymin, zmin, dh, lvl,
                                       indmax, ind2, indmax+3, indmax+4, 
                                       indmax+7, indmax+8, indmax+11, indmax+10);
    OctreeNode* node3 = new OctreeNode(xmin, ymin+dh, zmin, dh, lvl,
                                       indmax+2,indmax+4,indmax+1,ind4,
                                       indmax+9, indmax+10,indmax+13,indmax+12);
    OctreeNode* node4 = new OctreeNode(xmin+dh, ymin+dh, zmin, dh, lvl,
                                       indmax+4,indmax+3,ind3,indmax+1,
                                       indmax+10,indmax+11,indmax+14,indmax+13);
    OctreeNode* node5 = new OctreeNode(xmin, ymin, zmin+dh, dh, lvl,
                                       indmax+6, indmax+7, indmax+10, indmax+9,
                                       ind5, indmax+15, indmax+19,indmax+17);
    OctreeNode* node6 = new OctreeNode(xmin+dh, ymin, zmin+dh, dh, lvl, 
                                       indmax+7, indmax+8, indmax+11, indmax+10,
                                       indmax+15, ind6, indmax+18, indmax+19);
    OctreeNode* node7 = new OctreeNode(xmin, ymin+dh, zmin+dh, dh, lvl, 
                                       indmax+9, indmax+10,indmax+13,indmax+12,
                                       indmax+17,indmax+19,indmax+16,ind8);
    OctreeNode* node8 = new OctreeNode(xmin+dh, ymin+dh, zmin+dh, dh, lvl,
                                       indmax+10,indmax+11,indmax+14,indmax+13,
                                       indmax+19,indmax+18,ind7,indmax+16);
    indmax+= 20;
    OctreeNode** nodep1 = &node1; OctreeNode** nodep2 = &node2;
    OctreeNode** nodep3 = &node3; OctreeNode** nodep4 = &node4;
    OctreeNode** nodep5 = &node5; OctreeNode** nodep6 = &node6;
    OctreeNode** nodep7 = &node7; OctreeNode** nodep8 = &node8;
    
    (*currentp)->setNext1(*nodep1); 
    (*currentp)->setNext2(*nodep2);   
    (*currentp)->setNext3(*nodep3); 
    (*currentp)->setNext4(*nodep4); 
    (*currentp)->setNext5(*nodep5); 
    (*currentp)->setNext6(*nodep6);   
    (*currentp)->setNext7(*nodep7); 
    (*currentp)->setNext8(*nodep8);
    return *currentp; 
  }
}
//=============================================================================
OctreeNode* 
splitNode2D(OctreeNode* current, E_Int levelMax, E_Int& split, E_Int& indmax)
{
  split = 0;
  if ( current->getLevel() == levelMax ) return current;
  if ( current == NULL ) return current;
  if ( current->getNext1() != NULL ) { split = -1; return current;}
  split = 1;
  E_Int lvl = current->getLevel();
  OctreeNode** currentp = &current;
  E_Float dh0 = (*currentp)->getDh();
  E_Float xmin = (*currentp)->getXmin(); //E_Float xmax = xmin+dh0;
  E_Float ymin = (*currentp)->getYmin(); //E_Float ymax = ymin+dh0;  
  E_Float zmin = (*currentp)->getZmin(); //E_Float zmax = zmin+dh0; 
  E_Int ind1 = (*currentp)->getInd1(); E_Int ind2 = (*currentp)->getInd2();
  E_Int ind3 = (*currentp)->getInd3(); E_Int ind4 = (*currentp)->getInd4();
  E_Float dh = 0.5*dh0; 
  lvl = lvl + 1;
  OctreeNode* node1 = new OctreeNode(xmin, ymin, zmin, dh, lvl, ind1, indmax, indmax+4, indmax+2);
  OctreeNode* node2 = new OctreeNode(xmin+dh, ymin, zmin, dh, lvl, indmax, ind2, indmax+3, indmax+4);
  OctreeNode* node3 = new OctreeNode(xmin, ymin+dh, zmin, dh, lvl, indmax+2, indmax+4, indmax+1, ind4);
  OctreeNode* node4 = new OctreeNode(xmin+dh, ymin+dh, zmin, dh, lvl, indmax+4, indmax+3, ind3, indmax+1);
  indmax+= 5;
  OctreeNode** nodep1 = &node1;
  OctreeNode** nodep2 = &node2;
  OctreeNode** nodep3 = &node3;
  OctreeNode** nodep4 = &node4;

  (*currentp)->setNext1(*nodep1); 
  (*currentp)->setNext2(*nodep2);   
  (*currentp)->setNext3(*nodep3); 
  (*currentp)->setNext4(*nodep4); 

  return *currentp;  
}
//=============================================================================
OctreeNode* updateVoisin1(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData> stack;
  OctreeNode* voisin = node->getVoisin1();
  E_Int lvl0 = node->getLevel();
  E_Float xminv, yminv, zminv, xmaxv, ymaxv, zmaxv, dhv;
  E_Float dh = node->getDh(); 
  E_Float xmin = node->getXmin(); //E_Float xmax = xmin+dh; 
  E_Float ymin = node->getYmin(); E_Float ymax = ymin+dh;
  E_Float zmin = node->getZmin(); E_Float zmax = zmin+dh;
  vector<OctreeNode*> candidats;
  stackData dataForStack;

  if ( voisin == NULL ) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if ( voisin->getNext2() == NULL && voisin->getNext4() == NULL &&
       voisin->getNext6() == NULL && voisin->getNext8() == NULL ) return *nodep;

  if ( voisin->getNext2() != NULL ) 
  {
    dataForStack.current = voisin->getNext2(); stack.push(dataForStack);
  }
  if ( voisin->getNext4() != NULL ) 
  {
    dataForStack.current = voisin->getNext4(); stack.push(dataForStack);
  }
  if ( voisin->getNext6() != NULL ) 
  {
    dataForStack.current = voisin->getNext6(); stack.push(dataForStack);
  }
  if ( voisin->getNext8() != NULL ) 
  {
    dataForStack.current = voisin->getNext8(); stack.push(dataForStack);
  }
  stackData s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    dhv = stackp->getDh(); 
    xminv = stackp->getXmin(); xmaxv = xminv+dhv;
    yminv = stackp->getYmin(); ymaxv = yminv+dhv;
    zminv = stackp->getZmin(); zmaxv = zminv+dhv;

    if ( stackp->getNext2() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext2();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmin-xmaxv,eps) == true )
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
      if( K_FUNC::fEqualZero(xmin-xmaxv,eps) == true )
      {
        if ( yminv < ymax - eps && ymaxv > ymin + eps && zminv < zmax-eps && zmaxv > zmin+eps) 
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
    if ( stackp->getNext8() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext8();stack.push(dataForStack);}
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
OctreeNode* updateVoisin2(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData> stack;
  OctreeNode* voisin = node->getVoisin2();
  E_Int lvl0 = node->getLevel();
  vector<OctreeNode*> candidats;
  stackData dataForStack;
  E_Float xminv, yminv, zminv, ymaxv, zmaxv, dhv;
  E_Float dh = node->getDh(); 
  E_Float xmin = node->getXmin(); E_Float xmax = xmin+dh; 
  E_Float ymin = node->getYmin(); E_Float ymax = ymin+dh;
  E_Float zmin = node->getZmin(); E_Float zmax = zmin+dh;

  if ( voisin == NULL ) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if ( voisin->getNext1() == NULL && voisin->getNext3() == NULL &&
       voisin->getNext5() == NULL && voisin->getNext7() == NULL ) return *nodep;
  if ( voisin->getNext1() != NULL ) 
  {
    dataForStack.current = voisin->getNext1(); stack.push(dataForStack);
  }
  if ( voisin->getNext3() != NULL ) 
  {
    dataForStack.current = voisin->getNext3(); stack.push(dataForStack);
  }
  if ( voisin->getNext5() != NULL ) 
  {
    dataForStack.current = voisin->getNext5(); stack.push(dataForStack);
  }
  if ( voisin->getNext7() != NULL ) 
  {
    dataForStack.current = voisin->getNext7(); stack.push(dataForStack);
  }
  stackData s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    dhv = stackp->getDh(); 
    xminv = stackp->getXmin(); //xmaxv = xminv+dhv;
    yminv = stackp->getYmin(); ymaxv = yminv+dhv;
    zminv = stackp->getZmin(); zmaxv = zminv+dhv;

    if ( stackp->getNext1() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext1(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(xmax-xminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if (yminv < ymax-eps && ymaxv > ymin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }      
    }
    if ( stackp->getNext3() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext3(); stack.push(dataForStack);}
    }
    else 
    {
      xminv = stackp->getXmin();
      if(K_FUNC::fEqualZero(xmax-xminv,eps) == true) //&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
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
      if( K_FUNC::fEqualZero(xmax-xminv,eps) == true )//&& stackp->getLevel() < lvl0 )
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
      xminv = stackp->getXmin();
      if( K_FUNC::fEqualZero(xmax-xminv,eps) == true )//&& stackp->getLevel() < lvl0 )
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
OctreeNode* updateVoisin3(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData> stack;
  OctreeNode* voisin = node->getVoisin3();
  E_Int lvl0 = node->getLevel();
  E_Float xminv, yminv, zminv, xmaxv, ymaxv, zmaxv, dhv;
  E_Float dh = node->getDh(); 
  E_Float xmin = node->getXmin(); E_Float xmax = xmin+dh; 
  E_Float ymin = node->getYmin(); //E_Float ymax = ymin+dh;
  E_Float zmin = node->getZmin(); E_Float zmax = zmin+dh;
  vector<OctreeNode*> candidats;
  stackData dataForStack;

  if ( voisin == NULL ) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if ( voisin->getNext3() == NULL && voisin->getNext4() == NULL &&
       voisin->getNext7() == NULL && voisin->getNext8() == NULL ) return *nodep;

  if ( voisin->getNext3() != NULL ) 
  {
    dataForStack.current = voisin->getNext3(); stack.push(dataForStack);
  }
  if ( voisin->getNext4() != NULL ) 
  {
    dataForStack.current = voisin->getNext4(); stack.push(dataForStack);
  }
  if ( voisin->getNext7() != NULL ) 
  {
    dataForStack.current = voisin->getNext7(); stack.push(dataForStack);
  }
  if ( voisin->getNext8() != NULL ) 
  {
    dataForStack.current = voisin->getNext8(); stack.push(dataForStack);
  }
  stackData s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    dhv = stackp->getDh(); 
    xminv = stackp->getXmin(); xmaxv = xminv+dhv;
    yminv = stackp->getYmin(); ymaxv = yminv+dhv;
    zminv = stackp->getZmin(); zmaxv = zminv+dhv;

    if ( stackp->getNext3() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext3(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymin-ymaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
          candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext4() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext4(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymin-ymaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
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
      if( K_FUNC::fEqualZero(ymin-ymaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
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
OctreeNode* updateVoisin4(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData> stack;
  OctreeNode* voisin = node->getVoisin4();
  E_Int lvl0 = node->getLevel();
  E_Float xminv, yminv, zminv, xmaxv, zmaxv, dhv;
  E_Float dh = node->getDh(); 
  E_Float xmin = node->getXmin(); E_Float xmax = xmin+dh; 
  E_Float ymin = node->getYmin(); E_Float ymax = ymin+dh;
  E_Float zmin = node->getZmin(); E_Float zmax = zmin+dh;
  vector<OctreeNode*> candidats;
  stackData dataForStack;

  if ( voisin == NULL ) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if ( voisin->getNext1() == NULL && voisin->getNext2() == NULL &&
       voisin->getNext5() == NULL && voisin->getNext6() == NULL ) return *nodep;
  if ( voisin->getNext1() != NULL ) 
  {
    dataForStack.current = voisin->getNext1(); stack.push(dataForStack);
  }
  if ( voisin->getNext2() != NULL ) 
  {
    dataForStack.current = voisin->getNext2(); stack.push(dataForStack);
  }
  if ( voisin->getNext5() != NULL ) 
  {
    dataForStack.current = voisin->getNext5(); stack.push(dataForStack);
  }
  if ( voisin->getNext6() != NULL ) 
  {
    dataForStack.current = voisin->getNext6(); stack.push(dataForStack);
  }
  stackData s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    dhv = stackp->getDh(); 
    xminv = stackp->getXmin(); xmaxv = xminv+dhv;
    yminv = stackp->getYmin(); //ymaxv = yminv+dhv;
    zminv = stackp->getZmin(); zmaxv = zminv+dhv;
 
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
    if ( stackp->getNext5() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext5();stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(ymax-yminv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( xminv < xmax-eps && xmaxv > xmin+eps && zminv < zmax-eps && zmaxv > zmin+eps) 
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
OctreeNode* updateVoisin5(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData> stack;
  OctreeNode* voisin = node->getVoisin5();
  E_Int lvl0 = node->getLevel();
  E_Float xminv, yminv, zminv, xmaxv, ymaxv, zmaxv, dhv;
  E_Float dh = node->getDh(); 
  E_Float xmin = node->getXmin(); E_Float xmax = xmin+dh; 
  E_Float ymin = node->getYmin(); E_Float ymax = ymin+dh;
  E_Float zmin = node->getZmin(); //E_Float zmax = zmin+dh;
  vector<OctreeNode*> candidats;
  stackData dataForStack;

  if ( voisin == NULL ) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if ( voisin->getNext5() == NULL && voisin->getNext6() == NULL &&
       voisin->getNext7() == NULL && voisin->getNext8() == NULL ) return *nodep;
  if ( voisin->getNext5() != NULL ) 
  {
    dataForStack.current = voisin->getNext5(); stack.push(dataForStack);
  }
  if ( voisin->getNext6() != NULL ) 
  {
    dataForStack.current = voisin->getNext6(); stack.push(dataForStack);
  }
  if ( voisin->getNext7() != NULL ) 
  {
    dataForStack.current = voisin->getNext7(); stack.push(dataForStack);
  }
  if ( voisin->getNext8() != NULL ) 
  {
    dataForStack.current = voisin->getNext8(); stack.push(dataForStack);
  }
  stackData s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    dhv = stackp->getDh(); 
    xminv = stackp->getXmin(); xmaxv = xminv+dhv;
    yminv = stackp->getYmin(); ymaxv = yminv+dhv;
    zminv = stackp->getZmin(); zmaxv = zminv+dhv;

    if ( stackp->getNext5() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext5(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(zmin-zmaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
      {
        if ( yminv < ymax-eps && ymaxv > ymin+eps && xminv < xmax-eps && xmaxv > xmin+eps) 
          candidats.push_back(stackp); 
      }
    }
    if ( stackp->getNext6() != NULL )
    {
      if ( stackp->getLevel() < lvl0 ) 
      {dataForStack.current = stackp->getNext6(); stack.push(dataForStack);}
    }
    else 
    {
      if( K_FUNC::fEqualZero(zmin-zmaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
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
      if( K_FUNC::fEqualZero(zmin-zmaxv,eps) == true )//&& stackp->getLevel() < lvl0 )
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
OctreeNode* updateVoisin6(OctreeNode* node)
{
  OctreeNode** nodep = &node;
  stack<stackData> stack;
  OctreeNode* voisin = node->getVoisin6();
  E_Int lvl0 = node->getLevel();
  E_Float xminv, yminv, zminv, xmaxv, ymaxv, dhv;
  E_Float dh = node->getDh(); 
  E_Float xmin = node->getXmin(); E_Float xmax = xmin+dh; 
  E_Float ymin = node->getYmin(); E_Float ymax = ymin+dh;
  E_Float zmin = node->getZmin(); E_Float zmax = zmin+dh;

  vector<OctreeNode*> candidats;
  stackData dataForStack;

  if ( voisin == NULL ) return *nodep;
  if ( lvl0-voisin->getLevel() < 2 ) return *nodep;
  if ( voisin->getNext1() == NULL && voisin->getNext2() == NULL &&
       voisin->getNext3() == NULL && voisin->getNext4() == NULL ) return *nodep;
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
  if ( voisin->getNext4() != NULL ) 
  {
    dataForStack.current = voisin->getNext4(); stack.push(dataForStack);
  }
  stackData s0;
  E_Float eps = 1.e-6;
  OctreeNode* stackp;

  while (stack.size() != 0)
  {
    s0 = stack.top(); stack.pop(); stackp = s0.current;
    dhv = stackp->getDh(); 
    xminv = stackp->getXmin(); xmaxv = xminv+dhv;
    yminv = stackp->getYmin(); ymaxv = yminv+dhv;
    zminv = stackp->getZmin(); //zmaxv = zminv+dhv;

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
