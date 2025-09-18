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

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Modification des pts EX double wall */
//=============================================================================
PyObject* K_CONNECTOR::changeWallEX(PyObject* self, PyObject* args)
{
  PyObject *arrayEX, *arrayNodes, *arrayCenters,*firstWallCenters;//domaine a interpoler
  PyObject *projectSurfArrays; // liste des surfaces de projection : TRI
  E_Float planarTol;
  if (!PYPARSETUPLE_(args, OOOO_ O_ R_,
                     &arrayEX, &arrayNodes, &arrayCenters, &firstWallCenters, &projectSurfArrays, &planarTol))
    return NULL;

  if (PyList_Check(firstWallCenters) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "changeWall: 4th argument must be a list.");
    return NULL;
  }
  if (PyList_Check(projectSurfArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "changeWall: 5th argument must be a list.");
    return NULL;
  }  
     
  // Check: ptsEX
  E_Int imEX, jmEX, kmEX;
  FldArrayF* fEX; FldArrayI* cnEX;
  char* varStringEX; char* eltTypeEX;
  E_Int res = K_ARRAY::getFromArray3(arrayEX, varStringEX, 
                                     fEX, imEX, jmEX, kmEX, cnEX, eltTypeEX); 
  if (res != 2) 
  {
    if (res == 1) RELEASESHAREDS(arrayEX, fEX);
    PyErr_SetString(PyExc_TypeError, 
                    "changeWallEX: 1st arg must be unstructured.");
    return NULL;
  }
  if (fEX->getSize() == 0) 
  { 
    RELEASESHAREDU(arrayEX, fEX, cnEX); Py_INCREF(arrayEX);
    return arrayEX;
  }
  // Check: coordonnees des pts EX 
  E_Int posxEX = K_ARRAY::isCoordinateXPresent(varStringEX);
  E_Int posyEX = K_ARRAY::isCoordinateYPresent(varStringEX);
  E_Int poszEX = K_ARRAY::isCoordinateZPresent(varStringEX);
  if (posxEX == -1 || posyEX == -1 || poszEX == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "changeWallEX: 1st arg must contain coordinates.");
    RELEASESHAREDU(arrayEX, fEX, cnEX); return NULL;
  }
  //Check si les variables indcell1, nodemin et EXdir existent
  E_Int posind1 = K_ARRAY::isNamePresent("indcell1", varStringEX);
  E_Int posind2 = K_ARRAY::isNamePresent("indcell2", varStringEX);
  E_Int posnode = K_ARRAY::isNamePresent("nodemin", varStringEX);
  E_Int posdir = K_ARRAY::isNamePresent("EXdir", varStringEX);
  if (posind1 == -1 || posind2 == -1 || posnode == -1 || posdir == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "changeWallEX: 1st arg must contain indcell1, indcell1, nodemin, posdir variables.");
    RELEASESHAREDU(arrayEX, fEX, cnEX); return NULL;
  }
  posxEX++; posyEX++; poszEX++; posind1++; posind2++; posnode++; posdir++;

  // Check: coordonnees en noeuds de z
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray3(arrayNodes, varString, f, im, jm, km, cn, 
                               eltType); 
  if (res != 1) 
  {
    if (res == 2) RELEASESHAREDU(arrayNodes,f,cn);
    RELEASESHAREDU(arrayEX, fEX, cnEX);
    PyErr_SetString(PyExc_TypeError, 
                    "changeWallEX: 2nd arg must be structured.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "changeWallEX: 2nd arg must contain coordinates.");
    RELEASESHAREDU(arrayEX, fEX, cnEX);
    RELEASESHAREDS(arrayNodes, f);  
    return NULL;
  }
  posx++; posy++; posz++;

  // Check: coordonnees en centres de z
  E_Int imc, jmc, kmc;
  FldArrayF* fc;
  FldArrayI* cnc;
  char* varStringc;
  res = K_ARRAY::getFromArray3(arrayCenters, varStringc, fc, imc, jmc, kmc, 
                               cnc, eltType); 
  if (res != 1) 
  {
    RELEASESHAREDU(arrayEX,fEX,cnEX); 
    RELEASESHAREDS(arrayNodes,f);
    if (res == 2) RELEASESHAREDU(arrayCenters,fc,cnc);
    PyErr_SetString(PyExc_TypeError, 
                    "changeWallEX: 3rd arg must be structured.");
    return NULL;
  }  
  E_Int posxc = K_ARRAY::isCoordinateXPresent(varStringc);
  E_Int posyc = K_ARRAY::isCoordinateYPresent(varStringc);
  E_Int poszc = K_ARRAY::isCoordinateZPresent(varStringc);
  E_Int poscc = K_ARRAY::isCellNatureField2Present(varStringc);
  if (posxc == -1 || posyc == -1 || poszc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "changeWallEX: 3rd arg must contain coordinates.");
    RELEASESHAREDU(arrayEX, fEX, cnEX);
    RELEASESHAREDS(arrayNodes, f);  
    RELEASESHAREDS(arrayCenters, fc); 
    return NULL;
  }
  if (poscc == -1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "changeWall: 3rd arg must contain cellN variable.");
    RELEASESHAREDU(arrayEX, fEX, cnEX);
    RELEASESHAREDS(arrayNodes, f); 
    RELEASESHAREDS(arrayCenters, fc);
    return NULL; 
  }
  posxc++; posyc++; poszc++;  poscc++;

  // Check: premiers centres pres des parois
  E_Int im1, jm1, km1;
  FldArrayF* f1;
  FldArrayI* cn1;
  char* varString1;
  char* eltType1;
  res = K_ARRAY::getFromArray3(firstWallCenters, varString1, f1, 
                               im1, jm1, km1, cn1, eltType1); 
  if ( res != 2 ) 
  {
    RELEASESHAREDU(arrayEX,fEX,cnEX);
    RELEASESHAREDS(arrayNodes,f); 
    RELEASESHAREDS(arrayCenters,fc);
    if ( res == 1 ) RELEASESHAREDS(firstWallCenters,f1);
    PyErr_SetString(PyExc_TypeError, 
                    "changeWallEX: 4th arg must be unstructured.");
    return NULL;   
  }  
  if ( strcmp(eltType1,"NODE") != 0 ) 
  {
    RELEASESHAREDU(arrayEX,fEX,cnEX);
    RELEASESHAREDS(arrayNodes,f);
    RELEASESHAREDS(arrayCenters,fc); 
    RELEASESHAREDU(firstWallCenters,f1, cn1);
    PyErr_SetString(PyExc_TypeError, 
                    "changeWallEX: 4th arg must be a NODE array.");
    return NULL;   
  }
  E_Int posdir1 = K_ARRAY::isNamePresent("dir1",varString1);
  E_Int posdir2 = K_ARRAY::isNamePresent("dir2",varString1);
  E_Int posdir3 = K_ARRAY::isNamePresent("dir3",varString1);
  E_Int poshw   = K_ARRAY::isNamePresent("hmax",varString1);
  E_Int posindw = K_ARRAY::isNamePresent("indcellw",varString1);
  if ( posdir1 == -1 || posdir2 == -1 || posdir3 == -1 || poshw == -1 || posindw == -1 ) 
  {
    RELEASESHAREDU(arrayEX,fEX,cnEX);
    RELEASESHAREDS(arrayNodes,f);
    RELEASESHAREDS(arrayCenters,fc);
    RELEASESHAREDU(firstWallCenters,f1, cn1);
    PyErr_SetString(PyExc_TypeError, 
                    "changeWallEX: 4th arg must contain dir1,dir2,dir3,indcellw,hmax variables.");
    return NULL;
  }
  posdir1++; posdir2++; posdir3++; poshw++; posindw++;

  /*-------------------- surfaces de projection -----------------------------*/
  vector<E_Int> resl;
  vector<char*> structVarString;
  vector<char*> unstrVarString;
  vector<FldArrayF*> structF;
  vector<FldArrayF*> unstrF;
  vector<E_Int> nit; 
  vector<E_Int> njt; 
  vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = true;
  E_Bool skipUnstructured = false;
  E_Bool skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    projectSurfArrays, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltTypet, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nu = unstrF.size();
  if ( isOk == -1 )
  {
    PyErr_SetString(PyExc_TypeError,"changeWallEX: 5th list of arrays is not valid.");
    for (E_Int iu = 0; iu < nu; iu++)
      RELEASESHAREDU(objut[iu], unstrF[iu], cnt[iu]);
    RELEASESHAREDU(arrayEX,fEX,cnEX);
    RELEASESHAREDS(arrayNodes,f);
    RELEASESHAREDS(arrayCenters,fc);
    RELEASESHAREDU(firstWallCenters,f1, cn1);
    return NULL;
  } 
  E_Int posxi, posyi, poszi, poshi, posci;
  vector<E_Int> posxt; vector<E_Int> posyt; vector<E_Int> poszt; 
  vector<E_Int> posht; vector<E_Int> posct;
  E_Int nzones = unstrF.size();
  for (E_Int i = 0; i < nzones; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(unstrVarString[i]); posxi++; 
    posyi = K_ARRAY::isCoordinateYPresent(unstrVarString[i]); posyi++;
    poszi = K_ARRAY::isCoordinateZPresent(unstrVarString[i]); poszi++;
    poshi = K_ARRAY::isNamePresent("hmax", unstrVarString[i]); poshi++;
    posci = K_ARRAY::isCellNatureField2Present(unstrVarString[i]); posci++;
    if ( poshi == 0 ) 
    {
      PyErr_SetString(PyExc_TypeError,"changeWallEX: hmax variable missing in 5th argument.");
      for (E_Int iu = 0; iu < nu; iu++)
        RELEASESHAREDU(objut[iu], unstrF[iu], cnt[iu]);
      RELEASESHAREDU(arrayEX,fEX,cnEX);
      RELEASESHAREDS(arrayNodes,f);
      RELEASESHAREDS(arrayCenters,fc);
      RELEASESHAREDU(firstWallCenters,f1, cn1);
      return NULL;
    }
    if ( posci == 0 ) 
    {
      PyErr_SetString(PyExc_TypeError,"changeWallEX: cellN variable missing in 5th argument.");
      for (E_Int iu = 0; iu < nu; iu++)
        RELEASESHAREDU(objut[iu], unstrF[iu], cnt[iu]);
      RELEASESHAREDU(arrayEX,fEX,cnEX);
      RELEASESHAREDS(arrayNodes,f);
      RELEASESHAREDS(arrayCenters,fc);
      RELEASESHAREDU(firstWallCenters,f1, cn1);
      return NULL;
    }
    posxt.push_back(posxi); posyt.push_back(posyi); poszt.push_back(poszi); posht.push_back(poshi); posct.push_back(posci);
  }
  /*-------------------- Fin des verifs --------------------------------------*/
  E_Int nfld = fEX->getNfld(); E_Int npts = fEX->getSize();
  E_Int nelts = cnEX->getSize(); E_Int nvert = cnEX->getNfld();
  PyObject* tpl = K_ARRAY::buildArray(nfld, varStringEX, npts, nelts, -1, "NODE");
  E_Float* fcpEX = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fcEX(npts, nfld, fcpEX, true); fcEX = *fEX;
  E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
  FldArrayI cno(nelts, nvert, cnnp, true); cno = *cnEX;

  changeWallEX(fcEX.getSize(), fcEX.begin(posxEX), fcEX.begin(posyEX), fcEX.begin(poszEX),
               fcEX.begin(posind1), fcEX.begin(posind2), fcEX.begin(posdir), fcEX.begin(posnode), 
               im, jm, km, f->begin(posx), f->begin(posy), f->begin(posz),
               imc, jmc, kmc, fc->begin(posxc), fc->begin(posyc), fc->begin(poszc), fc->begin(poscc),
               f1->getSize(), f1->begin(posindw), f1->begin(posdir1), f1->begin(posdir2), f1->begin(posdir3),f1->begin(poshw),
               posxt, posyt, poszt, posht, posct, cnt, unstrF, planarTol);
            
  RELEASESHAREDU(arrayEX,fEX,cnEX);
  RELEASESHAREDS(arrayNodes,f);
  RELEASESHAREDS(arrayCenters,fc);
  RELEASESHAREDU(firstWallCenters,f1, cn1);
  for (E_Int iu = 0; iu < nu; iu++)
    RELEASESHAREDU(objut[iu], unstrF[iu], cnt[iu]);
  return tpl;
}
//==============================================================================
/* */
//==============================================================================
void K_CONNECTOR::changeWallEX(
  E_Int nEXPts, E_Float* xEX, E_Float* yEX, E_Float* zEX,
  E_Float* neighbourCellsp1, E_Float* neighbourCellsp2, E_Float* dirEXt, E_Float* nodeMinp,
  E_Int im, E_Int jm, E_Int km, E_Float* xn, E_Float* yn, E_Float* zn,
  E_Int imc, E_Int jmc, E_Int kmc, E_Float* xc, E_Float* yc, E_Float* zc, E_Float* cellnc,
  E_Int nbCentersW, E_Float* indicesw, E_Float* dirw1, E_Float* dirw2, E_Float* dirw3, E_Float* hmaxw,
  vector<E_Int> posxt, vector<E_Int> posyt, vector<E_Int> poszt, vector<E_Int> posht, vector<E_Int> posct, 
  vector<FldArrayI*>& cnt, vector<FldArrayF*>& unstrF, E_Float planartol)
{
  E_Float coefhmax = 2.; // tolerance de projection : coefhmax * hmax
  E_Float tolbb = 1.e-6;

  //Creation des bbtrees pour chaque surface
  E_Int nzones = unstrF.size();
  typedef K_SEARCH::BoundingBox<3>  BBox3DType; 
  vector<K_SEARCH::BbTree3D*> vectOfBBTrees;
  E_Float minB[3];  E_Float maxB[3];
  vector< vector<BBox3DType*> > vectOfBoxes;// a detruire a la fin
  for (E_Int v = 0; v < nzones; v++)
  {
    E_Int nelts2 = cnt[v]->getSize();
    vector<BBox3DType*> boxes(nelts2);// liste des bbox de ts les elements de a2
    K_FLD::FldArrayF bbox(nelts2,6);// xmin, ymin, zmin, xmax, ymax, zmax
    K_COMPGEOM::boundingBoxOfUnstrCells(*cnt[v], unstrF[v]->begin(posxt[v]),
                                        unstrF[v]->begin(posyt[v]), unstrF[v]->begin(poszt[v]),  
                                        bbox);
    E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
    E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
    E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
    for (E_Int et = 0; et < nelts2; et++)
    {
      minB[0] = xminp[et]; minB[1] = yminp[et]; minB[2] = zminp[et];
      maxB[0] = xmaxp[et]; maxB[1] = ymaxp[et]; maxB[2] = zmaxp[et]; 
      boxes[et] = new BBox3DType(minB, maxB);
    }
    vectOfBoxes.push_back(boxes);
    // Build the box tree.
    K_SEARCH::BbTree3D* bbtree = new K_SEARCH::BbTree3D(boxes);
    vectOfBBTrees.push_back(bbtree);
  }
  E_Int imcjmc = imc*jmc;
  E_Int ncells = imcjmc*kmc;
  E_Float delta = K_CONST::E_MAX_FLOAT;
  E_Float deltax=0.; E_Float deltay=0.; E_Float deltaz=0.;
  E_Float xp, yp, zp, alpha, alpha1, alpha2;
  E_Int noet, ind1, ind2, indil, indt1, indt2, indt3;
  E_Int nob, nov1, nov2, nov3;
  E_Int nbB, dir1, dir2, dir3, dir, indA;
  E_Bool isProjected;
  E_Float xA, yA, zA, xaEX, yaEX, zaEX, xbEX, ybEX, zbEX, dAP2;
  E_Float dxa, dya, dza, hmax1, hmax2, hmax, distnew;
  vector<E_Int> indicesElts; vector<E_Int> candidates;
  E_Float pr1[3]; E_Float pr2[3];
  E_Int dirEX, indEX, indcell1, indicew, nodemin, indEX2;
  E_Int iA, jA, kA;

  // Poids du vecteur delta pour chq pt au dessus du pt pr�s de la paroi � interpoler : si i < no1 : alpha = 1, si no1 <= i< no2 : alpha decroissant, 0 ensuite
  vector<E_Float> xbt; vector<E_Float> ybt; vector<E_Float> zbt; vector<E_Int> dirt;vector<E_Float> hmaxt;

  // Indirection : isWallEX[indEX] retourne le num�ro du premier centre paroi associ� dans indicesw
  vector<E_Int> isWallEX(nEXPts);
  FldArrayI indirEX(ncells,6); indirEX.setAllValuesAt(-1);
  E_Int* indirEX1 = indirEX.begin(1);
  E_Int* indirEX2 = indirEX.begin(2);
  E_Int* indirEX3 = indirEX.begin(3);
  E_Int* indirEX4 = indirEX.begin(4);
  E_Int* indirEX5 = indirEX.begin(5);
  E_Int* indirEX6 = indirEX.begin(6);
  E_Int* indirEXp=NULL;
  vector<E_Int>  dirEXp(nEXPts);
  E_Int no1 = 50; E_Int no2 = 80;  E_Float no1f = no1/100.; E_Float no2f = no2/100.;
  for (indEX = 0; indEX < nEXPts; indEX++)
  {
    isWallEX[indEX] = -1;
    indcell1 = E_Int(neighbourCellsp1[indEX]);
    //indcell2 = E_Int(neighbourCellsp2[indEX]);
    for (E_Int ii = 0; ii < nbCentersW; ii++)
    {
      indA = E_Int(indicesw[ii]);
      if ( indA == indcell1 ){isWallEX[indEX] = ii; break;}
    }
    dirEX = E_Int(dirEXt[indEX]);
    if (dirEX == 1)//pt EX a gauche en i
    {
      dirEXp[indEX] = 1; indirEX1[indcell1] = indEX;
    }
    else if (dirEX == 0)//pt EX a droite en i
    {
      dirEXp[indEX] =-1; indirEX2[indcell1] = indEX;
    }
    else if (dirEX == 3)//pt EX a gauche en j
    {
      dirEXp[indEX] = 2; indirEX3[indcell1] = indEX;
    }
    else if (dirEX == 2)//pt EX a droite en j
    {
      dirEXp[indEX] =-2; indirEX4[indcell1] = indEX;
    }
    else if (dirEX == 5)//pt EX a gauche en k
    {
      dirEXp[indEX] = 3; indirEX5[indcell1] = indEX;
    }
    else if (dirEX == 4)//pt EX a droite en k
    {
      dirEXp[indEX] =-3; indirEX6[indcell1] = indEX;
    }
  }
  
  for (indEX = 0; indEX < nEXPts; indEX++)
  {
    indicew = isWallEX[indEX];
    if (indicew != -1 ) // le pt EX a pour centre interpole d origine un premier pt pres de la paroi
    {
      xaEX = xEX[indEX]; yaEX = yEX[indEX]; zaEX = zEX[indEX];
      dirEX = dirEXp[indEX];
      nodemin = E_Int(nodeMinp[indEX]);

      //determination du pt B
      xbt.clear(); ybt.clear(); zbt.clear(); dirt.clear();
      dir1 = E_Int(dirw1[indicew]);
      dir2 = E_Int(dirw2[indicew]);
      dir3 = E_Int(dirw3[indicew]);
      hmax1 = hmaxw[indicew];
      indA = E_Int(indicesw[indicew]);//centre dont provient le pt EX
      xA = xc[indA]; yA = yc[indA]; zA = zc[indA];
      kA = indA/imcjmc;
      jA = (indA-kA*imcjmc)/imc;
      iA = indA-jA*imc-kA*imcjmc;
      // on recupere le pt B si le pt A est un pt paroi dans la direction ortho a i                                 
      if ( dir1 == 1 || dir1 == -1 )
      {
        dir = dir1;
        compAbovePointEX(dirEX, nodemin, indA, 1, im, jm, imc, jmc, xn, yn, zn, xbEX, ybEX, zbEX);       
        xbt.push_back(xbEX); ybt.push_back(ybEX); zbt.push_back(zbEX); dirt.push_back(dir); 
      }
      if ( dir2 == 1 || dir2 == -1 )
      {
        dir = dir2*2;
        compAbovePointEX(dirEX, nodemin, indA, 2, im, jm, imc, jmc, xn, yn, zn, xbEX, ybEX, zbEX);       
        xbt.push_back(xbEX); ybt.push_back(ybEX); zbt.push_back(zbEX); dirt.push_back(dir); 
      }
      if ( dir3 == 1 || dir3 == -1 )
      {
        dir = dir3*2;
        compAbovePointEX(dirEX, nodemin, indA, 3, im, jm, imc, jmc, xn, yn, zn, xbEX, ybEX, zbEX);       
        xbt.push_back(xbEX); ybt.push_back(ybEX); zbt.push_back(zbEX); dirt.push_back(dir); 
      }
      delta = K_CONST::E_MAX_FLOAT; deltax = 0.; deltay = 0.; deltaz = 0.; isProjected = false;
      nbB = xbt.size();
      for (nob = 0; nob < nbB; nob++)
      {
        xbEX = xbt[nob]; ybEX = ybt[nob]; zbEX = zbt[nob]; dir = dirt[nob];
        for (E_Int v = 0; v < nzones; v++)//parcours des surfaces de projection
        {
          E_Int* cn1 = cnt[v]->begin(1);
          E_Int* cn2 = cnt[v]->begin(2);
          E_Int* cn3 = cnt[v]->begin(3);
          E_Float* hmaxTri = unstrF[v]->begin(posht[v]);
          E_Float* cellNTri = unstrF[v]->begin(posct[v]);
          xp = K_CONST::E_MAX_FLOAT; yp = K_CONST::E_MAX_FLOAT; zp = K_CONST::E_MAX_FLOAT;
          // Preconditionnement : on ne prend que les indices des elts 
          K_SEARCH::BbTree3D* bbtree = vectOfBBTrees[v];
          //dirx = xaEX-xbEX; diry = yaEX-ybEX; dirz = zaEX-zbEX;// BA
          pr1[0] = xaEX; pr1[1] = yaEX; pr1[2] = zaEX;
          pr2[0] = xbEX; pr2[1] = ybEX; pr2[2] = zbEX;
          indicesElts.clear(); candidates.clear();
          bbtree->getIntersectingBoxes(pr1, pr2, indicesElts, tolbb);
          for (unsigned int noe = 0; noe < indicesElts.size(); noe++)
          {
            indt1 = cn1[indicesElts[noe]]-1;
            indt2 = cn2[indicesElts[noe]]-1;
            indt3 = cn3[indicesElts[noe]]-1;
            if ( cellNTri[indt1] != 0. ||  cellNTri[indt2] != 0. ||  cellNTri[indt3] != 0. )
              candidates.push_back(indicesElts[noe]);
          } 
          // projection suivant la direction BA sur la surface opposee
          noet = K_COMPGEOM::projectDir(xaEX, yaEX, zaEX, xaEX-xbEX, yaEX-ybEX, zaEX-zbEX, 
                                        unstrF[v]->begin(posxt[v]), unstrF[v]->begin(posyt[v]), unstrF[v]->begin(poszt[v]),
                                        candidates, *cnt[v], xp, yp, zp);
          if ( noet != -1 && xp < K_CONST::E_MAX_FLOAT && yp < K_CONST::E_MAX_FLOAT && zp < K_CONST::E_MAX_FLOAT) 
          {
            dxa = xp-xaEX; dya = yp-yaEX; dza = zp-zaEX; dAP2 = dxa*dxa + dya*dya + dza*dza;
            nov1 = cn1[noet]-1; nov2 = cn2[noet]-1;nov3 = cn3[noet]-1;
            hmax2 = (hmaxTri[nov1]+hmaxTri[nov2]+hmaxTri[nov3])/3.;//moyenne sur le triangle de hmax
            hmax = K_FUNC::E_max(hmax1, hmax2); hmax = hmax*hmax;
            if ( dAP2 < coefhmax*hmax && dAP2 < delta)  
            {delta = dAP2; deltax = dxa; deltay = dya; deltaz = dza; isProjected = true;}
            else if ( hmax < 0.1*K_CONST::E_GEOM_CUTOFF && dAP2 < planartol)
            {
              delta = dAP2; deltax = dxa; deltay = dya; deltaz = dza; isProjected = true; 
            }
          }//elt trouve
        }//parcours de ts les surfaces de projection
        if (isProjected == false) goto nextptB;
//         cerr << "Wall double definition : point of dir "<< dir << " " << xaEX <<" "<< yaEX <<" "<< zaEX << " projected to "
//              << xaEX+deltax <<" "<< yaEX+deltay <<" "<< zaEX+deltaz << endl;
        if (dirEX == 1) indirEXp = indirEX1;//imin
        else if (dirEX==-1) indirEXp = indirEX2;//imax
        else if (dirEX== 2)  indirEXp = indirEX3;//jmin
        else if (dirEX==-2) indirEXp = indirEX4;//jmax
        else if (dirEX== 3)  indirEXp = indirEX5;//kmin
        else if (dirEX==-3) indirEXp = indirEX6;//kmax

        if (dir == 1 || dir == -1)
        {
          dir1 = K_FUNC::E_sign(dir);
          if ( dir1 == 1)
          {
            ind1 = E_Int(no1f*(imc-1)) + jA * imc + kA * imcjmc;
            ind2 = E_Int(no2f*(imc-1)) + jA * imc + kA * imcjmc;
          }
          else 
          {
            ind1 = E_Int((1-no1f)*(imc-1)) + jA * imc + kA * imcjmc;
            ind2 = E_Int((1-no2f)*(imc-1)) + jA * imc + kA * imcjmc;
          }
         
          // calcul de alpha1 et alpha2
          if ( ind1 < ncells && ind1 >= 0) alpha1 = (xc[ind1]-xA)*(xc[ind1]-xA)+(yc[ind1]-yA)*(yc[ind1]-yA)+(zc[ind1]-zA)*(zc[ind1]-zA);
          else alpha1 = 0.;
          if ( ind2 < ncells && ind2 >= 0) alpha2 = (xc[ind2]-xA)*(xc[ind2]-xA)+(yc[ind2]-yA)*(yc[ind2]-yA)+(zc[ind2]-zA)*(zc[ind2]-zA);
          else alpha2 = 0.;
          
          for (E_Int il = 0; il < imc; il++)//on parcourt la ligne i et les pts EX a gauche ou a droite
          {
            indil = indA + dir1 * il; 
            indEX2 = indirEXp[indil];
            if ( indEX2 != -1 ) 
            {
              // Approx shift
              distnew=(xEX[indEX2]-xaEX)*(xEX[indEX2]-xaEX)+(yEX[indEX2]-yaEX)*(yEX[indEX2]-yaEX)+(zEX[indEX2]-zaEX)*(zEX[indEX2]-zaEX);
              if ( dir == 1 ) //i=1
              {
                if ( indil < ind1 ) alpha = 1.;
                else if ( indil < ind2 ) alpha = (alpha2-distnew)/(alpha2-alpha1);
                else alpha = 0.;
              }
              else 
              {
                if ( indil > ind1 ) alpha = 1.;
                else if ( indil > ind2 ) alpha = (alpha2-distnew)/(alpha2-alpha1);
                else alpha = 0.;
              }
              if ( il < no1 ) alpha = 1.;
              else if ( il >= no1 && il < no2 && alpha2 > 0.) alpha = (alpha2-distnew)/(alpha2-alpha1);
              else alpha = 0.;
              xEX[indEX2] = xEX[indEX2]+deltax*alpha;
              yEX[indEX2] = yEX[indEX2]+deltay*alpha;
              zEX[indEX2] = zEX[indEX2]+deltaz*alpha; 
            }
          }
        }//dir = 1 ou -1
        
        else if (dir == 2 || dir == -2)
        {
          dir2 = K_FUNC::E_sign(dir);
          if ( dir2 == 1 )
          {
            ind1 = iA + E_Int(no1f*(jmc-1))*imc + kA * imcjmc;
            ind2 = iA + E_Int(no2f*(jmc-1))*imc + kA * imcjmc;            
          }
          else 
          {
            ind1 = iA + E_Int((1-no1f)*(jmc-1))*imc + kA * imcjmc;
            ind2 = iA + E_Int((1-no2f)*(jmc-1))*imc + kA * imcjmc;   
          }


          // calcul de alpha1 et alpha2
          if ( ind1 < ncells && ind1 >= 0) alpha1 = (xc[ind1]-xA)*(xc[ind1]-xA)+(yc[ind1]-yA)*(yc[ind1]-yA)+(zc[ind1]-zA)*(zc[ind1]-zA);
          else alpha1 = 0.;
          if ( ind2 < ncells && ind2 >= 0) alpha2 = (xc[ind2]-xA)*(xc[ind2]-xA)+(yc[ind2]-yA)*(yc[ind2]-yA)+(zc[ind2]-zA)*(zc[ind2]-zA);
          else alpha2 = 0.;         
                  
          for (E_Int jl = 0; jl < jmc; jl++)//parcours de la ligne jl
          {
            indil = indA + dir2 * jl * imc;      
            indEX2 = indirEXp[indil];     
            // calcul de alpha1 et alpha2
            if ( indEX2 !=-1 )
            {
              // Approx shift
              distnew=(xEX[indEX2]-xaEX)*(xEX[indEX2]-xaEX)+(yEX[indEX2]-yaEX)*(yEX[indEX2]-yaEX)+(zEX[indEX2]-zaEX)*(zEX[indEX2]-zaEX);
              if ( dir == 2 ) //j=1
              {
                if ( indil < ind1 ) alpha = 1.;
                else if ( indil < ind2 ) alpha = (alpha2-distnew)/(alpha2-alpha1);
                else alpha = 0.;
              }
              else 
              {
                if ( indil > ind1 ) alpha = 1.;
                else if ( indil > ind2 ) alpha = (alpha2-distnew)/(alpha2-alpha1);
                else alpha = 0.;
              }
              xEX[indEX2] = xEX[indEX2]+deltax*alpha;
              yEX[indEX2] = yEX[indEX2]+deltay*alpha;
              zEX[indEX2] = zEX[indEX2]+deltaz*alpha; 
            }          
          }                    
        }//dir = 2 ou -2
        else if ( dir == 3 || dir == -3 )
        {
          dir3 = K_FUNC::E_sign(dir);
          if (dir3 == 1 )
          {
            ind1 = iA + jA * imc + E_Int(no1f * (kmc-1)) * imcjmc;
            ind2 = iA + jA * imc + E_Int(no2f * (kmc-1)) * imcjmc;       
          }
          else 
          {
            ind1 = iA + jA * imc + E_Int((1-no1f) * (kmc-1)) * imcjmc;
            ind2 = iA + jA * imc + E_Int((1-no2f) * (kmc-1)) * imcjmc;
          }

          // calcul de alpha1 et alpha2
          if ( ind1 < ncells && ind1 >= 0) alpha1 = (xc[ind1]-xA)*(xc[ind1]-xA)+(yc[ind1]-yA)*(yc[ind1]-yA)+(zc[ind1]-zA)*(zc[ind1]-zA);
          else alpha1 = 0.;
          if ( ind2 < ncells && ind2 >= 0) alpha2 = (xc[ind2]-xA)*(xc[ind2]-xA)+(yc[ind2]-yA)*(yc[ind2]-yA)+(zc[ind2]-zA)*(zc[ind2]-zA);
          else alpha2 = 0.; 
          for (E_Int kl = 0; kl < kmc; kl++)
          {
            indil = indA + dir3 * kl * imcjmc;
            indEX2 = indirEXp[indil];    
            if ( indEX2 != -1 )
            {               
              // Approx shift
              distnew=(xEX[indEX2]-xaEX)*(xEX[indEX2]-xaEX)+(yEX[indEX2]-yaEX)*(yEX[indEX2]-yaEX)+(zEX[indEX2]-zaEX)*(zEX[indEX2]-zaEX);
              if ( dir == 3 ) //k=1
              {
                if ( indil < ind1 ) alpha = 1.;
                else if ( indil < ind2 ) alpha = (alpha2-distnew)/(alpha2-alpha1);
                else alpha = 0.;
              }
              else 
              {
                if ( indil > ind1 ) alpha = 1.;
                else if ( indil > ind2 ) alpha = (alpha2-distnew)/(alpha2-alpha1);
                else alpha = 0.;
              }
              
              xEX[indEX2] = xEX[indEX2]+deltax*alpha;
              yEX[indEX2] = yEX[indEX2]+deltay*alpha;
              zEX[indEX2] = zEX[indEX2]+deltaz*alpha; 
            }
          }
        }//dir = 3 ou -3
        nextptB:;
      }//parcours de ts les pts B valides      
    }// fin pt EX paroi 
  }//parcours de ts les pts EX  
  E_Int nboxes = vectOfBoxes.size();
  for (E_Int v0 = 0; v0 < nboxes; v0++)
  {
    vector<BBox3DType*>& boxes = vectOfBoxes[v0];
    E_Int size = boxes.size();
    for (E_Int v = 0; v < size; v++) delete boxes[v];
    delete vectOfBBTrees[v0];
  }
  vectOfBoxes.clear(); vectOfBBTrees.clear();
}
//===========================================================================================
/* Calcul du centre d interface au dessus du centre d interface dont le noeud min est indn
   IN: dirEX : direction du pt EX (interface i=1 dirEX=1, i=imax dirEX=-1 etc)
   IN: indn: indice du noeud min associe au pt EX
   IN: indc: indice du centre dont provient le pt EX
   IN: dir: direction de recherche du above pt associe a indcell
   IN: im, jm: dimensions du maillage en noeuds
   IN: imc, jmc : dimensions du maillage en centres
   IN: xn,yn,zn: coordonnees des noeuds
   OUT: xb,yb,zb: coordonnees du pt EX au dessus de indEX */
//===========================================================================================
void K_CONNECTOR::compAbovePointEX(E_Int dirEX, E_Int indn, E_Int indc, E_Int dir,
                                   E_Int im, E_Int jm, E_Int imc, E_Int jmc,
                                   E_Float* xn, E_Float* yn, E_Float* zn,
                                   E_Float& xb, E_Float& yb, E_Float& zb)
{
  E_Int imjm = im*jm;  E_Int imcjmc = imc*jmc;
  E_Int ind1, ind2, ind3, ind4;
  E_Int kn = indn/imjm+1;// demarre a 1
  E_Int jn = (indn-(kn-1)*imjm)/im+1;
  E_Int in = indn-(jn-1)*im-(kn-1)*imjm+1;  
  E_Int ip=1, jp=1, kp=1;// indices du premier sommet de l'interface au dessus
  E_Int kc = indc/imcjmc+1;
  E_Int jc = (indc-(kc-1)*imcjmc)/imc+1;
  E_Int ic = indc-(jc-1)*imc-(kc-1)*imcjmc+1;  

  E_Int dirEX0 = K_FUNC::E_abs(dirEX);
  if (dirEX0 == 1 )// interface EX en i
  {
    switch (dir) 
    {
      case 3:
        if ( kc == 1 ) {ip = in; jp = jn; kp = kn+1;}
        else {ip = in; jp = jn; kp = kn-1;}
        break;
      case 2:
        if ( jc == 1 ){ip = in; jp = jn+1; kp = kn;}
        else {ip = in; jp = jn-1; kp = kn;}
        break;
      case 1:
        if ( ic == 1 ) {ip = in+1; jp = jn; kp = kn;}
        else {ip = in-1; jp = jn; kp = kn;}
        break;
    }
    ind1 = (ip-1) + (jp-1)*im +(kp-1)*imjm;
    ind2 = ind1 + im;
    ind3 = ind1 + imjm;
    ind4 = ind3 + im;   
  }
  else if ( dirEX0 == 2 ) //interface en j
  {
    switch (dir) 
    {
      case 3:
        if ( kc == 1 ) {ip = in; jp = jn; kp = kn+1;}
        else {ip = in; jp = jn; kp = kn-1;}
        break;
      case 2:
        if ( jc == 1 ){ip = in; jp = jn+1; kp = kn;}
        else {ip = in; jp = jn-1; kp = kn;}
        break;
      case 1:
        if ( ic == 1 ) {ip = in+1; jp = jn; kp = kn;}
        else {ip = in-1; jp = jn; kp = kn;}
        break;
    }
    ind1 = (ip-1) + (jp-1)*im +(kp-1)*imjm;
    ind2 = ind1 + 1;
    ind3 = ind1 + imjm;
    ind4 = ind3 + 1;    
  }
  else //interface en k 
  {
    switch (dir) 
    {
      case 3:
        if ( kc == 1 ) {ip = in; jp = jn; kp = kn+1;}
        else {ip = in; jp = jn; kp = kn-1;}
        break;
      case 2:
        if ( jc == 1 ){ip = in; jp = jn+1; kp = kn;}
        else {ip = in; jp = jn-1; kp = kn;}
        break;
      case 1:
        if ( ic == 1 ) {ip = in+1; jp = jn; kp = kn;}
        else {ip = in-1; jp = jn; kp = kn;}
        break;
    }
    ind1 = (ip-1) + (jp-1)*im +(kp-1)*imjm;
    ind2 = ind1 + 1;
    ind3 = ind1 + im;
    ind4 = ind3 + 1; 
  }
  xb = 0.25*(xn[ind1]+xn[ind2]+xn[ind3]+xn[ind4]);
  yb = 0.25*(yn[ind1]+yn[ind2]+yn[ind3]+yn[ind4]);
  zb = 0.25*(zn[ind1]+zn[ind2]+zn[ind3]+zn[ind4]);
}
