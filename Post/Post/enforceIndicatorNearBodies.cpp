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

# include "post.h"
# include "Search/BbTree.h"
# include "Search/KdTree.h"
# include "Fld/ArrayAccessor.h"

using namespace std;
using namespace K_FLD;
using namespace K_SEARCH;
extern "C"
{
  void k6boundboxunstr_(const E_Int& npts, 
                        const E_Float* x, const E_Float* y, const E_Float* z, 
                        E_Float& xmax, E_Float& ymax, E_Float& zmax, 
                        E_Float& xmin, E_Float& ymin, E_Float& zmin);
}

//=============================================================================
/* Force l'indicateur a 0 pres des corps, si le pas d'espace respecte 
   le snear */
//=============================================================================
PyObject* K_POST::enforceIndicatorNearBodies(PyObject* self, PyObject* args)
{
  PyObject *indicator, *octree, *bodies; 
  if (!PyArg_ParseTuple(args, "OOO", &indicator, &octree, &bodies)) 
    return NULL;
  
  if (PyList_Size(bodies) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "enforceIndicatorNearBodies: 3rd argument is an empty list.");
    return NULL;
  }
  // Verif octree HEXA/QUAD
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(octree, varString, f, ni, nj, nk, 
                                    cn, eltType, true);
  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "enforceIndicatorNearBodies: array must be unstructured.");
    RELEASESHAREDB(res, octree, f, cn); return NULL;   
  }
  E_Int dim = 3;
  if (strcmp(eltType,"HEXA") == 0) dim = 3;
  else if (strcmp(eltType,"QUAD") == 0) dim = 2;
  else
  {
    PyErr_SetString(PyExc_TypeError, 
                    "enforceIndicatorNearBodies: array must be HEXA or QUAD.");
    RELEASESHAREDU(octree, f, cn); return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "enforceIndicatorNearBodies: coordinates not found in array.");
    RELEASESHAREDU(octree, f, cn); return NULL;
  }
  posx++; posy++; posz++;

  // Verif indicator
  E_Int nii, nji, nki;
  FldArrayF* fi; FldArrayI* cni;
  char* varStringi; char* eltTypei;
  E_Int resi = K_ARRAY::getFromArray(indicator, varStringi, fi, 
                                     nii, nji, nki, cni, eltTypei, true);
  if (resi != 1 && resi != 2) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "enforceIndicatorNearBodies: indic array must be structured.");
    RELEASESHAREDU(octree, f, cn); return NULL;
  }
  E_Int posi = K_ARRAY::isNamePresent("indicator", varStringi);
  if (posi == -1) 
  { 
    RELEASESHAREDB(resi, indicator, fi, cni); RELEASESHAREDU(octree, f, cn); 
    printf("Warning: enforceIndicatorNearBodies: no refinement indicator given. Nothing done."); 
    return indicator;
  }
  posi++;
  if (fi->getSize() != cn->getSize()) 
  {
    RELEASESHAREDB(resi, indicator, fi, cni); RELEASESHAREDU(octree, f, cn); 
    printf("Warning: enforceIndicatorNearBodies: refinement indicator size must be equal to the number of elements. Nothing done."); 
    return indicator;
  }

  // Verif des bodies
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt; vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = true;
  E_Boolean skipUnstructured = false;
  E_Boolean skipDiffVars = true;
  res = K_ARRAY::getFromArrays(bodies, resl, structVarString, unstrVarString,
                               structF, unstrF, nit, njt, nkt, cnt, eltTypet,
                               objst, objut,
                               skipDiffVars, skipNoCoord, 
                               skipStructured, skipUnstructured, true);
  if (res == -1)
  {
    for (unsigned int nos = 0; nos < unstrF.size(); nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]); 
    RELEASESHAREDB(resi, indicator, fi, cni); RELEASESHAREDU(octree, f, cn); 
    PyErr_SetString(PyExc_TypeError, 
                    "enforceIndicatorNearBodies: 3rd arg is not valid.");
    return NULL;
  }
  E_Int nzones = unstrF.size();

  for (E_Int i = 0; i < nzones; i++)
  {
    if (strcmp(eltTypet[i],"TRI") == 0 && dim != 3) 
    {
      for (unsigned int nos = 0; nos < unstrF.size(); nos++)
        RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]); 
      RELEASESHAREDB(resi, indicator, fi, cni); RELEASESHAREDU(octree, f, cn); 
      PyErr_SetString(PyExc_TypeError, 
                      "enforceIndicatorNearBodies: 3rd arg must be a list of TRI zones.");
      return NULL;
    }
    else if (strcmp(eltTypet[i],"BAR") == 0 && dim != 2)  
    {
      for (unsigned int nos = 0; nos < unstrF.size(); nos++)
        RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]); 
      RELEASESHAREDB(resi, indicator, fi, cni); RELEASESHAREDU(octree, f, cn); 
      PyErr_SetString(PyExc_TypeError, 
                      "enforceIndicatorNearBodies: 3rd arg must be a list of BAR zones.");
      return NULL;
    }
    else if (strcmp(eltTypet[i],"BAR") != 0 && strcmp(eltTypet[i],"TRI") != 0)
    {
      for (unsigned int nos = 0; nos < unstrF.size(); nos++)
        RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
      RELEASESHAREDB(resi, indicator, fi, cni); RELEASESHAREDU(octree, f, cn); 
      PyErr_SetString(PyExc_TypeError, 
                      "enforceIndicatorNearBodies: 3rd arg must be a list of TRI or BAR zones.");
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
  /*-----------FIN DES VERIFS ------------------*/
  
  //creation des bbox trees pour ttes les surfaces
  typedef K_SEARCH::BoundingBox<3>  BBox3DType;
  vector<K_SEARCH::BbTree3D*> bboxtrees(nzones);
  vector< vector<BBox3DType*> > vectOfBBoxes;// pour etre detruit a la fin
  E_Float minB[3];  E_Float maxB[3];
  E_Int posx2, posy2, posz2, nelts2;
  E_Float xmin =  K_CONST::E_MAX_FLOAT;
  E_Float ymin =  K_CONST::E_MAX_FLOAT;
  E_Float zmin =  K_CONST::E_MAX_FLOAT;
  E_Float xmax = -K_CONST::E_MAX_FLOAT;
  E_Float ymax = -K_CONST::E_MAX_FLOAT;
  E_Float zmax = -K_CONST::E_MAX_FLOAT;
  if (dim == 2) { zmin = 0.; zmax = 0.; }
  E_Float xminloc, yminloc, zminloc, xmaxloc, ymaxloc, zmaxloc;

  for (E_Int v = 0; v < nzones; v++)
  {
    FldArrayF& f2 = *unstrF[v]; FldArrayI& cn2 = *cnt[v];
    posx2 = posxt[v]; posy2 = posyt[v]; posz2 = poszt[v];
    // bounding box globale ?
    k6boundboxunstr_(f2.getSize(), 
                     f2.begin(posx2), f2.begin(posy2), f2.begin(posz2),
                     xmaxloc, ymaxloc, zmaxloc, xminloc, yminloc, zminloc);
   
    xmin = K_FUNC::E_min(xminloc,xmin); xmax = K_FUNC::E_max(xmaxloc,xmax);
    ymin = K_FUNC::E_min(yminloc,ymin); ymax = K_FUNC::E_max(ymaxloc,ymax);
    if (dim == 2) { zmin = 0.; zmax = 0.; }
    else 
    {zmin = K_FUNC::E_min(zminloc,zmin); zmax = K_FUNC::E_max(zmaxloc,zmax);}
    
    // Creation de la bboxtree
    nelts2 = cn2.getSize();
    vector<BBox3DType*> boxes(nelts2);// liste des bbox de ts les elements de a2
    FldArrayF bbox(nelts2,6);// xmin, ymin, zmin, xmax, ymax, zmax
    K_COMPGEOM::boundingBoxOfUnstrCells(cn2, f2.begin(posx2), f2.begin(posy2),f2.begin(posz2),bbox);
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
  for (unsigned int nos = 0; nos < unstrF.size(); nos++)
    RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]); 
  E_Float* indict = fi->begin(posi);

  vector<E_Int> indicesBB; E_Int nelts = cn->getSize();
  FldArrayF bboxo(nelts,6);// xmin, ymin, zmin, xmax, ymax, zmax de l octree
  K_COMPGEOM::boundingBoxOfUnstrCells(*cn, f->begin(posx), f->begin(posy),
                                      f->begin(posz),bboxo);
  RELEASESHAREDU(octree, f, cn);
 
  E_Float* xmino = bboxo.begin(1); E_Float* xmaxo = bboxo.begin(4);
  E_Float* ymino = bboxo.begin(2); E_Float* ymaxo = bboxo.begin(5);
  E_Float* zmino = bboxo.begin(3); E_Float* zmaxo = bboxo.begin(6);

  for (E_Int et = 0; et < nelts; et++)
  { 
    minB[0] = xmino[et]; minB[1] = ymino[et]; minB[2] = zmino[et];
    maxB[0] = xmaxo[et]; maxB[1] = ymaxo[et]; maxB[2] = zmaxo[et]; 
    for (E_Int v = 0; v < nzones; v++)
      bboxtrees[v]->getOverlappingBoxes(minB, maxB, indicesBB);
    if (indicesBB.size() != 0) indict[et] = -1000.;
    indicesBB.clear();
  }
  // nettoyages...
  for (E_Int v = 0; v < nzones; v++)
  {
    delete bboxtrees[v]; 
    E_Int nbboxes =  vectOfBBoxes[v].size();
    for (E_Int v2 = 0; v2 < nbboxes; v2++) delete vectOfBBoxes[v][v2];
  }
  vectOfBBoxes.clear();

  /*-----------CONSTRUCTION ARRAY DE SORTIE ------------------*/
  //buildArray
  PyObject* tpl;
  if (resi == 1) 
    tpl = K_ARRAY::buildArray(*fi, varStringi, nii, nji, nki);
  else 
    tpl = K_ARRAY::buildArray(*fi, varStringi, *cni, -1, eltTypei, 
                              false);
  RELEASESHAREDB(resi, indicator, fi, cni);
  return tpl;
}
