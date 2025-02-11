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
# include <iostream>
using namespace K_FLD;
using namespace std;

//=============================================================================
/*  Nouvel algo changeWall */
//=============================================================================
PyObject* K_CONNECTOR::changeWall(PyObject* self, PyObject* args)
{
  PyObject *arrayCenters, *firstWallCenters; //domaine a interpoler
  PyObject *projectSurfArrays; // liste des surfaces de projection: TRI
  E_Float planarTol; // tolerance de shift double wall dans les cas planaires
  if (!PYPARSETUPLE_(args, OOO_ R_, 
                     &arrayCenters, &firstWallCenters, &projectSurfArrays, &planarTol))
     return NULL;

  if (PyList_Check(firstWallCenters) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "changeWall: 2nd argument must be a list.");
    return NULL;
  }
  if (PyList_Check(projectSurfArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "changeWall: 3rd argument must be a list.");
    return NULL;
  }  

  // Check: coordonnees en centres de z
  E_Int imc, jmc, kmc;
  FldArrayF* fc; FldArrayI* cnc;
  char* varStringc; char* eltType;
  E_Int res = K_ARRAY::getFromArray(arrayCenters, varStringc, fc, 
                                    imc, jmc, kmc, cnc, eltType, true); 
  if (res != 1)
  {
    if (res == 2) RELEASESHAREDU(arrayCenters,fc,cnc);
    PyErr_SetString(PyExc_TypeError, 
                    "changeWall: 1st arg must be structured.");
    return NULL;   
  }  
  // Verif des coordonnees dans zc
  E_Int posxc = K_ARRAY::isCoordinateXPresent(varStringc);
  E_Int posyc = K_ARRAY::isCoordinateYPresent(varStringc);
  E_Int poszc = K_ARRAY::isCoordinateZPresent(varStringc);
  E_Int posc = K_ARRAY::isCellNatureField2Present(varStringc);

  if (posxc == -1 || posyc == -1 || poszc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "changeWall: 1st arg must contain coordinates.");
    RELEASESHAREDS(arrayCenters,fc); return NULL;
  }
  if (posc == -1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "changeWall: 1st arg must contain cellN variable.");
    RELEASESHAREDS(arrayCenters, fc); return NULL; 
  }
  posxc++; posyc++; poszc++; posc++;

  // Check: premiers centres pres des parois
  E_Int im1, jm1, km1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  res = K_ARRAY::getFromArray(firstWallCenters, varString1, f1, 
                              im1, jm1, km1, cn1, eltType1, true); 
  if (res != 2) 
  {
    if (res == 1) RELEASESHAREDS(firstWallCenters, f1);
    RELEASESHAREDS(arrayCenters,fc);
    PyErr_SetString(PyExc_TypeError, 
                    "changeWall: 2nd arg must be unstructured.");
    return NULL;   
  }  
  if (K_STRING::cmp(eltType1, "NODE") != 0) 
  {
    RELEASESHAREDS(arrayCenters,fc); RELEASESHAREDU(firstWallCenters,f1, cn1);
    PyErr_SetString(PyExc_TypeError, 
                    "changeWall: 2nd arg must be a NODE array.");
    return NULL;   
  }
  E_Int posdir1 = K_ARRAY::isNamePresent("dir1", varString1);
  E_Int posdir2 = K_ARRAY::isNamePresent("dir2", varString1);
  E_Int posdir3 = K_ARRAY::isNamePresent("dir3", varString1);
  E_Int poshw   = K_ARRAY::isNamePresent("hmax", varString1);
  E_Int posindw = K_ARRAY::isNamePresent("indcellw", varString1);
  if (posdir1 == -1 || posdir2 == -1 || posdir3 == -1 || 
      poshw == -1 || posindw == -1) 
  {
    RELEASESHAREDS(arrayCenters,fc); RELEASESHAREDU(firstWallCenters,f1, cn1);
    PyErr_SetString(PyExc_TypeError, 
                    "changeWall: 2nd arg must contain dir1,dir2,dir3,indcellw,hmax variables.");
    return NULL;
  }
  posdir1++; posdir2++; posdir3++; poshw++; posindw++;

  /*-------------------- surfaces de projection -----------------------------*/
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = true;
  E_Boolean skipUnstructured = false;
  E_Boolean skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    projectSurfArrays, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltTypet, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nu = unstrF.size();
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError, "changeWall: 3rd arg is not valid.");
    for (E_Int iu = 0; iu < nu; iu++)
      RELEASESHAREDU(objut[iu], unstrF[iu], cnt[iu]);
    RELEASESHAREDS(arrayCenters,fc); RELEASESHAREDU(firstWallCenters,f1, cn1);
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
    if (poshi == 0) 
    {
      PyErr_SetString(PyExc_TypeError,"changeWall: hmax variable missing in 3rd argument.");
      for (E_Int iu = 0; iu < nu; iu++)
        RELEASESHAREDU(objut[iu], unstrF[iu], cnt[iu]);
      RELEASESHAREDS(arrayCenters,fc); RELEASESHAREDU(firstWallCenters,f1, cn1);
      return NULL;
    }
    if (posci == 0) 
    {
      PyErr_SetString(PyExc_TypeError,"changeWall: cellN variable missing in 3rd argument.");
      for (E_Int iu = 0; iu < nu; iu++)
        RELEASESHAREDU(objut[iu], unstrF[iu], cnt[iu]);
      RELEASESHAREDS(arrayCenters,fc); RELEASESHAREDU(firstWallCenters,f1, cn1);
      return NULL;
    }
    posxt.push_back(posxi); posyt.push_back(posyi); poszt.push_back(poszi); 
    posht.push_back(poshi); posct.push_back(posci);
  }
  /*-------------------- Fin des verifs --------------------------------------*/
  PyObject* tpl = K_ARRAY::buildArray(fc->getNfld(), varStringc, imc, jmc, kmc);
  E_Float* fcp2 = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fc2(fc->getSize(), fc->getNfld(), fcp2, true); fc2.setAllValuesAt(*fc);
  changeWall(imc, jmc, kmc, fc2.begin(posc), 
             f1->getSize(), f1->begin(posindw), f1->begin(posdir1), f1->begin(posdir2),f1->begin(posdir3),
             f1->begin(poshw),
             posxt, posyt, poszt, posht, posct, cnt, unstrF, 
             fc->begin(posxc), fc->begin(posyc), fc->begin(poszc),
             fc2.begin(posxc), fc2.begin(posyc), fc2.begin(poszc), planarTol);

  // cleaning
  for (E_Int iu = 0; iu < nu; iu++)
    RELEASESHAREDU(objut[iu], unstrF[iu], cnt[iu]);
  RELEASESHAREDS(arrayCenters,fc);
  RELEASESHAREDU(firstWallCenters,f1, cn1);
  return tpl;
}
//=============================================================================
/* Nouvel algo de projection double wall */
//=============================================================================
void K_CONNECTOR::changeWall(
  E_Int imc, E_Int jmc, E_Int kmc, E_Float* cellN, 
  E_Int nbCentersW, E_Float* indicesw, 
  E_Float* dirw1, E_Float* dirw2, E_Float* dirw3, E_Float* hmaxw,
  vector<E_Int> posxt, vector<E_Int> posyt, vector<E_Int> poszt, 
  vector<E_Int> posht, vector<E_Int> posct,
  vector<FldArrayI*>& cnt, vector<FldArrayF*>& unstrF,
  E_Float* xc, E_Float* yc, E_Float* zc,
  E_Float* xc2, E_Float* yc2, E_Float* zc2, E_Float planartol)
{
  E_Float coefhmax = 10.; // tolerance de projection : coefhmax * hmax
  E_Float tolbb = 1.e-6;
  E_Int dirCoeff = 2; // si ghost cells pour trouver la premiere cellule non plate
  // Creation des bbtrees pour chaque surface de projection
  E_Int nzones = unstrF.size();
  typedef K_SEARCH::BoundingBox<3> BBox3DType;
  vector<K_SEARCH::BbTree3D*> vectOfBBTrees;
  E_Float minB[3]; E_Float maxB[3];
  vector< vector<BBox3DType*> > vectOfBoxes; // a detruire a la fin
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
  //E_Int ncells = imcjmc*kmc;
  E_Int cpt=0;
  E_Int nthreads = __NUMTHREADS__;
  E_Int* npts_proj = new E_Int [nthreads];
  
  #pragma omp parallel
  {
    E_Int ithread = __CURRENT_THREAD__;
    E_Float delta = K_CONST::E_MAX_FLOAT;
    E_Float deltax = 0.; E_Float deltay = 0.; E_Float deltaz = 0.;
    E_Float xp, yp, zp;
    E_Int noet, indt1, indt2, indt3;
    E_Int nob, nov1, nov2, nov3;
    E_Int nbB, dir1, dir2, dir3, dir, indA, indB;
    E_Int cpt_loc=0;
    E_Boolean isProjected;
    E_Float xa, ya, za, xb, yb, zb, dAP2, dirx, diry, dirz;
    E_Float dxa, dya, dza, hmax1, hmax2, hmax;
    vector<E_Int> indicesElts; vector<E_Int> candidates;
    E_Float pr1[3]; E_Float pr2[3];
    E_Int iA, jA, kA;
    vector<E_Float> xbt; vector<E_Float> ybt; vector<E_Float> zbt; vector<E_Int> dirt;vector<E_Float> hmaxt;
    // Projection des pts interpoles des first centers
    // Poids du vecteur delta pour chq pt au dessus du pt pres de la paroi a interpoler : 
    // Si i < no1 : alpha = 1, si no1 <= i< no2 : alpha decroissant, 0 ensuite
    
    #pragma omp for schedule(dynamic)
    for (E_Int ii = 0; ii < nbCentersW; ii++)
    {
      indA = E_Int(indicesw[ii]);
      if (cellN[indA] == 1.) goto nextpt;
      kA = indA/imcjmc;
      jA = (indA-kA*imcjmc)/imc;
      iA = indA-jA*imc-kA*imcjmc;
      xa = xc[indA]; ya = yc[indA]; za = zc[indA];
      //printf("xa %f %f %f\n", xa, ya, za);
      xbt.clear(); ybt.clear(); zbt.clear(); dirt.clear();
      dir1 = E_Int(dirw1[ii]);
      dir2 = E_Int(dirw2[ii]);
      dir3 = E_Int(dirw3[ii]);
      //printf("dir " SF_D3_ "\n", dir1, dir2, dir3);
      hmax1 = hmaxw[ii];

      // on recupere le pt B si le pt A est un pt paroi dans la direction ortho a i
      if (dir1 == 1 || dir1 == -1)
      {
        indB = indA + dirCoeff*dir1; dir = dir1;
        xb = xc[indB]; yb = yc[indB]; zb = zc[indB];
        xbt.push_back(xb); ybt.push_back(yb); zbt.push_back(zb); dirt.push_back(dir); 
      }
      if (dir2 == 1 || dir2 == -1)
      {
        indB = indA + dirCoeff*dir2*imc; dir = dir2*2;
        xb = xc[indB]; yb = yc[indB]; zb = zc[indB];
        xbt.push_back(xb); ybt.push_back(yb); zbt.push_back(zb); dirt.push_back(dir);
      }
      if (dir3 == 1 || dir3 == -1) 
      {
        indB = indA + dirCoeff*dir3*imcjmc; dir = dir3*3;
        xb = xc[indB]; yb = yc[indB]; zb = zc[indB];
        xbt.push_back(xb); ybt.push_back(yb); zbt.push_back(zb); dirt.push_back(dir);
      }
      //printf("xb %f %f %f\n", xb, yb, zb);
      delta = K_CONST::E_MAX_FLOAT; deltax = 0.; deltay = 0.; deltaz = 0.; isProjected = false;
      nbB = xbt.size();

      // on fait une projection suivant la direction AB
      for (nob = 0; nob < nbB; nob++)
      {
        xb = xbt[nob]; yb = ybt[nob]; zb = zbt[nob]; dir = dirt[nob];
        dirx = xa-xb; diry = ya-yb; dirz = za-zb;// BA

        for (E_Int v = 0; v < nzones; v++) // parcours des surfaces de projection
        {
          E_Int* cn1 = cnt[v]->begin(1);
          E_Int* cn2 = cnt[v]->begin(2);
          E_Int* cn3 = cnt[v]->begin(3);

          E_Float* hmaxTri = unstrF[v]->begin(posht[v]);
          E_Float* cellNTri = unstrF[v]->begin(posct[v]);
          xp = K_CONST::E_MAX_FLOAT; yp = K_CONST::E_MAX_FLOAT; zp = K_CONST::E_MAX_FLOAT;
          // Preconditionnement : on ne prend que les indices des elts 
          K_SEARCH::BbTree3D* bbtree = vectOfBBTrees[v];
          pr1[0] = xa; pr1[1] = ya; pr1[2] = za;
          pr2[0] = xa+dirx; pr2[1] = ya+diry; pr2[2] = za+dirz;
          indicesElts.clear(); candidates.clear();
          bbtree->getIntersectingBoxes(pr1, pr2, indicesElts, tolbb);
          for (size_t noe = 0; noe < indicesElts.size(); noe++)
          {
            indt1 = cn1[indicesElts[noe]]-1;
            indt2 = cn2[indicesElts[noe]]-1;
            indt3 = cn3[indicesElts[noe]]-1;
            if (cellNTri[indt1] > 0. ||  cellNTri[indt2] > 0. ||  cellNTri[indt3] > 0.)
              candidates.push_back(indicesElts[noe]);
          } 
          // projection suivant la direction BA sur la surface opposee
          noet = K_COMPGEOM::projectDir(xa, ya, za, dirx, diry, dirz, 
                                        unstrF[v]->begin(posxt[v]), unstrF[v]->begin(posyt[v]), unstrF[v]->begin(poszt[v]),
                                        candidates, *cnt[v], xp, yp, zp);
          if (noet != -1 && xp < K_CONST::E_MAX_FLOAT && yp < K_CONST::E_MAX_FLOAT && zp < K_CONST::E_MAX_FLOAT)
          {
            dxa = xp-xa; dya = yp-ya; dza = zp-za; dAP2 = dxa*dxa + dya*dya + dza*dza;
            nov1 = cn1[noet]-1; nov2 = cn2[noet]-1; nov3 = cn3[noet]-1;
            hmax2 = (hmaxTri[nov1]+hmaxTri[nov2]+hmaxTri[nov3])/3.; // moyenne sur le triangle de hmax
            hmax = K_FUNC::E_max(hmax1, hmax2); hmax = hmax*hmax;
            //printf("DAP2 %f delta=%f cmax=%f\n", dAP2, delta, coefhmax*hmax);
            if (dAP2 < coefhmax*hmax && dAP2 < delta) 
            { delta = dAP2; deltax = dxa; deltay = dya; deltaz = dza; isProjected = true; }
            else if (hmax < 0.1*K_CONST::E_GEOM_CUTOFF && dAP2 < planartol)
            { delta = dAP2; deltax = dxa; deltay = dya; deltaz = dza; isProjected = true; }
          }
        }
        if (isProjected == false) goto nextptB;      
        //   printf("Info: Wall double definition: point (%15.6f,%15.6f,%15.6f) of indices (" SF_D_ "," SF_D_ "," SF_D_ ") projected onto (%15.6f,%15.6f,%15.6f)\n",
        //                 xa, ya, za, iA+1, jA+1, kA+1, xa+deltax, ya+deltay,za+deltaz);

        cpt_loc++;
    
        // hope no race because shifted points are independant
        shiftAbovePoints(imc, jmc, kmc, dir, indA, iA, jA, kA, 
                         xa, ya, za, deltax, deltay, deltaz,
                         xc, yc, zc, cellN, xc2, yc2, zc2);   
      
        nextptB:;
      }//nob
      nextpt:;
    }// pour ts les pts interpoles
    npts_proj[ithread] = cpt_loc;
    // #pragma omp atomic update
    // cpt += cpt_loc; // ensures that race conditions are avoided
  }

  for (E_Int i = 0; i < nthreads; i++) cpt += npts_proj[i];
  printf("Info: changeWall: " SF_D_ " points have been projected\n", cpt);

  // cleanup
  E_Int nboxes = vectOfBoxes.size();
  for (E_Int v0 = 0; v0 < nboxes; v0++)
  {
    vector<BBox3DType*>& boxes = vectOfBoxes[v0];
    E_Int size = boxes.size();
    for (E_Int v = 0; v < size; v++) delete boxes[v];
    delete vectOfBBTrees[v0];
  }
  vectOfBoxes.clear(); vectOfBBTrees.clear();
  delete [] npts_proj;
}

//=============================================================================
void K_CONNECTOR::shiftAbovePoints(E_Int imc, E_Int jmc, E_Int kmc,
                                   E_Int dir,  E_Int indA, E_Int iA, E_Int jA, E_Int kA,
                                   E_Float xa, E_Float ya, E_Float za, 
                                   E_Float deltax, E_Float deltay, E_Float deltaz,
                                   E_Float* xc, E_Float* yc, E_Float* zc, E_Float* cellN,
                                   E_Float* xc2, E_Float* yc2, E_Float* zc2)
{
  E_Int no1 = 50; E_Int no2 = 80;  E_Float no1f = no1/100.; E_Float no2f = no2/100.;
  E_Int dir1, dir2, dir3, ind1, ind2, indil;
  E_Float alpha1, alpha2, alpha, distnew;

  E_Int imcjmc = imc*jmc;
  E_Int ncells = imcjmc*kmc;
  // Shifting all the points on the boundary
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
    if ( ind1 < ncells && ind1 >= 0) alpha1 = (xc[ind1]-xa)*(xc[ind1]-xa)+(yc[ind1]-ya)*(yc[ind1]-ya)+(zc[ind1]-za)*(zc[ind1]-za);
    else alpha1 = 0.;
    if ( ind2 < ncells && ind2 >= 0) alpha2 = (xc[ind2]-xa)*(xc[ind2]-xa)+(yc[ind2]-ya)*(yc[ind2]-ya)+(zc[ind2]-za)*(zc[ind2]-za);
    else alpha2 = 0.;
    for (E_Int il = 0; il < imc; il++)//on parcourt la ligne i
    {
      indil = indA + dir1 * il; 
      if (cellN[indil] > 1. )
      {
        // Approx shift
        distnew = (xc[indil]-xa)*(xc[indil]-xa)+(yc[indil]-ya)*(yc[indil]-ya)+(zc[indil]-za)*(zc[indil]-za);
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
        xc2[indil] = xc[indil]+deltax*alpha;
        yc2[indil] = yc[indil]+deltay*alpha;
        zc2[indil] = zc[indil]+deltaz*alpha; 
      }
    }
  }//dir = 1
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
    if ( ind1 < ncells && ind1 >= 0) alpha1 = (xc[ind1]-xa)*(xc[ind1]-xa)+(yc[ind1]-ya)*(yc[ind1]-ya)+(zc[ind1]-za)*(zc[ind1]-za);
    else alpha1 = 0.;
    if ( ind2 < ncells && ind2 >= 0) alpha2 = (xc[ind2]-xa)*(xc[ind2]-xa)+(yc[ind2]-ya)*(yc[ind2]-ya)+(zc[ind2]-za)*(zc[ind2]-za);
    else alpha2 = 0.;

    for (E_Int jl = 0; jl < jmc; jl++)
    {
      indil = indA + dir2 * jl * imc;           
      // calcul de alpha1 et alpha2
      if (cellN[indil] > 1. )
      {
        // Approx shift
        distnew = (xc[indil]-xa)*(xc[indil]-xa)+(yc[indil]-ya)*(yc[indil]-ya)+(zc[indil]-za)*(zc[indil]-za);
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
        xc2[indil] = xc[indil]+deltax*alpha;
        yc2[indil] = yc[indil]+deltay*alpha;
        zc2[indil] = zc[indil]+deltaz*alpha; 
      }          
    }
  }//dir = 2
  else //dir=3
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

    if ( ind1 < ncells && ind1 >= 0) alpha1 = (xc[ind1]-xa)*(xc[ind1]-xa)+(yc[ind1]-ya)*(yc[ind1]-ya)+(zc[ind1]-za)*(zc[ind1]-za);
    else alpha1 = 0.;
    if ( ind2 < ncells && ind2 >= 0) alpha2 = (xc[ind2]-xa)*(xc[ind2]-xa)+(yc[ind2]-ya)*(yc[ind2]-ya)+(zc[ind2]-za)*(zc[ind2]-za);
    else alpha2 = 0.;
     
    for (E_Int kl = 0; kl < kmc; kl++)
    {
      indil = indA + dir3 * kl * imcjmc;
      if (cellN[indil] > 1. )
      {               
        // Approx shift
        distnew = (xc[indil]-xa)*(xc[indil]-xa)+(yc[indil]-ya)*(yc[indil]-ya)+(zc[indil]-za)*(zc[indil]-za);
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
        xc2[indil] = xc[indil]+deltax*alpha;
        yc2[indil] = yc[indil]+deltay*alpha;
        zc2[indil] = zc[indil]+deltaz*alpha;  
      }
    }
  }//dir = 3
}
