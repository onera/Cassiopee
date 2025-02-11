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

// Interpolation on zones given a solution defined on a list of zones

# include "post.h"
# include "Nuga/include/KdTree.h"
# include "Nuga/include/BbTree.h"
# include "Nuga/include/ArrayAccessor.h"
using namespace K_FLD;
using namespace std;

#define CLOUDMAX 120;

// ============================================================================
/*  IN: nuage de points donneur defini par une zone NODE
    IN/OUT: zone surfacique triangulaire receptrice */
// ============================================================================
PyObject* K_POST::projectCloudSolution2Triangle(PyObject* self, PyObject* args)
{
  E_Int dimPb, ibm, old;
  PyObject *arrayR, *arrayD;
  if (!PYPARSETUPLE_(args, OO_ III_,
		     &arrayD, &arrayR, &dimPb, &ibm, &old))
  {
      return NULL;
  }

  /*-----------------------------------------------*/
  /* Extraction des infos sur le domaine recepteur */
  /*-----------------------------------------------*/
  E_Int imr, jmr, kmr;
  FldArrayF* fr; FldArrayI* cnr;
  char* varStringr; char* eltTyper;
  E_Int resr = K_ARRAY::getFromArray(arrayR, varStringr, fr,
                                     imr, jmr, kmr, cnr, eltTyper, true);
  if (resr != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "projectCloudSolution2Triangle: invalid receptor zone.");
    if (resr == 1) RELEASESHAREDB(resr, arrayR, fr, cnr);
    return NULL;
  }
  if (K_STRING::cmp(eltTyper, "TRI") != 0 && K_STRING::cmp(eltTyper, "BAR") != 0 && 
      K_STRING::cmp(eltTyper, "TRI*") != 0 && K_STRING::cmp(eltTyper, "BAR*") != 0  )
  {
    PyErr_SetString(PyExc_TypeError,
      "projectCloudSolution2Triangle: unstructured receptor zone must be TRI or BAR.");
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    return NULL;
  }

  /*---------------------------------------------*/
  /* Extraction des infos sur le domaine donneur */
  /*---------------------------------------------*/
  E_Int imd, jmd, kmd;
  FldArrayF* fd; FldArrayI* cnd;
  char* varStringd; char* eltTyped;
  E_Int resd = K_ARRAY::getFromArray(arrayD, varStringd, fd,
                                     imd, jmd, kmd, cnd, eltTyped, true);
  if (resd != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "projectCloudSolution2Triangle: 2nd arg is not a valid array.");
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    if (resd == 1) RELEASESHAREDB(resd, arrayD, fd, cnd);
    return NULL;
  }

  if (K_STRING::cmp(eltTyped, "NODE") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
      "projectCloudSolution2Triangle: unstructured donor zone must only be NODE.");
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    RELEASESHAREDB(resd, arrayD, fd, cnd);
    return NULL;
  }
  E_Int posxd = K_ARRAY::isCoordinateXPresent(varStringd);
  E_Int posyd = K_ARRAY::isCoordinateYPresent(varStringd);
  E_Int poszd = K_ARRAY::isCoordinateZPresent(varStringd);
  if (posxd == -1 || posyd == -1 || poszd == -1)
  {
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    RELEASESHAREDB(resd, arrayD, fd, cnd);
    PyErr_SetString(PyExc_TypeError,
                    "projectCloudSolution2Triangle: coordinates not found in donor zone.");
    return NULL;
  }

  E_Int posxr = K_ARRAY::isCoordinateXPresent(varStringr);
  E_Int posyr = K_ARRAY::isCoordinateYPresent(varStringr);
  E_Int poszr = K_ARRAY::isCoordinateZPresent(varStringr);
  if (posxr == -1 || posyr == -1 || poszr == -1)
  {
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    RELEASESHAREDB(resd, arrayD, fd, cnd);
    PyErr_SetString(PyExc_TypeError,
                    "projectCloudSolution2Triangle: coordinates not found in receptor zone.");
    return NULL;
  }
  //Variables a transferer
  E_Int poscd = K_ARRAY::isCellNatureField2Present(varStringd);
  E_Int poscr = K_ARRAY::isCellNatureField2Present(varStringr);
  char* varStringC; // chaine de caractere commune
  E_Int l = strlen(varStringr); varStringC = new char [l+1];
  // les positions demarrent a 1
  vector<E_Int> posvarsD; vector<E_Int> posvarsR;
  K_ARRAY::getPosition(varStringr, varStringd,
                       posvarsR, posvarsD, varStringC);
  delete [] varStringC;
  E_Int sizeVarsD = posvarsD.size();
  E_Int sizeVarsR = posvarsR.size();
  for (E_Int i = 0; i < sizeVarsD; i++) posvarsD[i] -= 1;
  for (E_Int i = 0; i < sizeVarsR; i++) posvarsR[i] -= 1;
  if (poscd != -1) posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), poscd), posvarsD.end());
  if (poscr != -1) posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), poscr), posvarsR.end());
  posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), posxd), posvarsD.end());
  posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), posyd), posvarsD.end());
  posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), poszd), posvarsD.end());
  posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), posxr), posvarsR.end());
  posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), posyr), posvarsR.end());
  posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), poszr), posvarsR.end());
  posxd++; posyd++; poszd++; posxr++; posyr++; poszr++;
  E_Int nbVars = posvarsD.size();
  E_Int order = 2; // ordre d'interpolation MLS
  const E_Int sizeBasis = K_FUNC::fact(dimPb+order-1)*1./(K_FUNC::fact(dimPb)*K_FUNC::fact(order-1));

  //1. Determination du(des) triangles contenant le point dnr
  E_Float* xtRcv = fr->begin(posxr);
  E_Float* ytRcv = fr->begin(posyr);
  E_Float* ztRcv = fr->begin(poszr);

  // Creation du kdtree des points de la triangulation
  ArrayAccessor<FldArrayF>* coordAcc = new ArrayAccessor<FldArrayF>(*fr, posxr, posyr, poszr);
  K_SEARCH::KdTree<FldArrayF>* kdt = new K_SEARCH::KdTree<FldArrayF>(*coordAcc);

  // Creation du bboxtree
  E_Int nbEltsR = cnr->getSize();
  E_Int nbPtsR = fr->getSize();

  typedef K_SEARCH::BoundingBox<3>  BBox3DType;
  vector<BBox3DType*> boxes(nbEltsR);// liste des bbox de ts les elements de a2
  K_FLD::FldArrayF bbox(nbEltsR,6);// xmin, ymin, zmin, xmax, ymax, zmax
  K_COMPGEOM::boundingBoxOfUnstrCells(*cnr, xtRcv, ytRcv, ztRcv, bbox);
  E_Float minB[3];  E_Float maxB[3];
  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
  E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
  E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
  for (E_Int et = 0; et < nbEltsR; et++)
  {
    minB[0] = xminp[et]; minB[1] = yminp[et]; minB[2] = zminp[et];
    maxB[0] = xmaxp[et]; maxB[1] = ymaxp[et]; maxB[2] = zmaxp[et];
    boxes[et] = new BBox3DType(minB, maxB);
  }
  // Build the box tree.
  K_SEARCH::BbTree3D bbtree(boxes);
  //projection des points de fd(nuage) sur fr(TRI)
  E_Float xo, yo, zo, rx, ry, rz, rad;
  E_Float pt[3];
  E_Float p0[3]; E_Float p1[3]; E_Float p2[3]; E_Float p[3];
  E_Int indev, indR;
  vector<E_Int> indicesBB; // liste des indices des facettes intersectant la bbox
  E_Int nbPtsD = fd->getSize();

  vector < vector<E_Int> > cveR(nbPtsR);//connectivite vertex/elts
  K_CONNECT::connectEV2VE(*cnr, cveR);

  vector< vector<E_Int> > listOfCloudPtsPerTri(nbEltsR);

  FldArrayI matching(nbPtsR); matching.setAllValuesAt(-1);
  E_Int* matchingP = matching.begin();
  E_Float* xtDnr = fd->begin(posxd);
  E_Float* ytDnr = fd->begin(posyd);
  E_Float* ztDnr = fd->begin(poszd);

  E_Float dist2proj = 0.;

  for (E_Int indD = 0; indD < nbPtsD; indD++)
  {
    // recherche du pt le plus proche P' de P
    pt[0] = xtDnr[indD]; pt[1] = ytDnr[indD]; pt[2] = ztDnr[indD];
    indR = kdt->getClosest(pt);

    // calcul de la bounding box de la sphere de rayon PP'
    rx = pt[0]-xtRcv[indR]; ry = pt[1]-ytRcv[indR]; rz = pt[2]-ztRcv[indR];
    rad = sqrt(rx*rx+ry*ry+rz*rz);
    if (rad < K_CONST::E_GEOM_CUTOFF)//matching
    {matchingP[indR] = indD;}
      
    else //projection 
    {
      minB[0] = pt[0]-rad; minB[1] = pt[1]-rad; minB[2] = pt[2]-rad;
      maxB[0] = pt[0]+rad; maxB[1] = pt[1]+rad; maxB[2] = pt[2]+rad;
      bbtree.getOverlappingBoxes(minB, maxB, indicesBB);
      indev = K_COMPGEOM::projectOrthoPrecond(pt[0], pt[1], pt[2],
                                              xtRcv, ytRcv, ztRcv, indicesBB, *cnr, 
                                              xo, yo, zo, p0, p1, p2, p);
      dist2proj = sqrt((pt[0]-xo)*(pt[0]-xo) + (pt[1]-yo)*(pt[1]-yo) + (pt[2]-zo)*(pt[2]-zo));
      indicesBB.clear();
      if (indev != -1 && !(dist2proj>1.e-6 && ibm==1)) // donor points should be already on the triangle or nearby for IBMs
      {
        listOfCloudPtsPerTri[indev].push_back(indD);
      }
    }
  }
  
  // cleanings
  delete kdt; delete coordAcc;
  E_Int size = boxes.size();
  for (E_Int v = 0; v < size; v++) delete boxes[v];
    
  E_Int sizeMinOfCloud = sizeBasis+1;
  vector< vector<E_Int> > cEEN(nbEltsR);
  K_CONNECT::connectEV2EENbrs(eltTyper, nbPtsR, *cnr, cEEN); 
 
  E_Int axisConst[3]; // =1 si l'axe est une direction constante, 0 sinon
  axisConst[0] = axisConst[1] = axisConst[2] = 0;
  if (dimPb == 2) axisConst[2] = 1;// pas de z
  E_Float radius[3]; // longueurs des axes de l'ellipse
  radius[0] = radius[1] = radius[2] = 1.;
  E_Float axis[9];  // Axes de l'ellipse
  axis[0] = 1.; axis[1] = 0.; axis[2] = 0.;
  axis[3] = 0.; axis[4] = 1.; axis[5] = 0.;
  axis[6] = 0.; axis[7] = 0.; axis[8] = 1.;

  E_Int crsize = cnr->getSize()*cnr->getNfld();
  PyObject* tpl = K_ARRAY::buildArray(fr->getNfld(), varStringr,
    fr->getSize(), cnr->getSize(),-1, eltTyper, false, crsize);
  E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
  K_KCORE::memcpy__(cnnp, cnr->begin(), cnr->getSize()*cnr->getNfld());
  E_Float* ptrFieldOut = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fieldROut(fr->getSize(), fr->getNfld(), ptrFieldOut, true);
  fieldROut = *fr;

  vector<E_Int> indicesExtrap;
  for (E_Int indR = 0; indR < nbPtsR; indR++)
  {
    E_Int indD = matchingP[indR];
    if (indD>-1)
    {
      for (E_Int posv = 0; posv < nbVars; posv++)
      {
        E_Float* varD = fd->begin(posvarsD[posv]+1);
        E_Float* varR = fieldROut.begin(posvarsR[posv]+1);
        varR[indR] = varD[indD];
      }
    }
    else // non matching: interpolation, then extrapolation
    {
      vector<E_Int>& eltsVoisins = cveR[indR];
      // create the cloud of donor points
      vector<E_Int> listOfCloudPtsPerVertex;
      for (size_t noe = 0; noe < eltsVoisins.size(); noe++)
      {
        E_Int indEltR = eltsVoisins[noe];// index of one elt with indR as vertex
        vector<E_Int>& projectedPts = listOfCloudPtsPerTri[indEltR];// projected points on this elt
        for (size_t nop = 0; nop < projectedPts.size(); nop++)
          listOfCloudPtsPerVertex.push_back(projectedPts[nop]);
      }
      sort(listOfCloudPtsPerVertex.begin(), listOfCloudPtsPerVertex.end());
      listOfCloudPtsPerVertex.erase(unique(listOfCloudPtsPerVertex.begin(), listOfCloudPtsPerVertex.end()), 
                                    listOfCloudPtsPerVertex.end());
      // if the MLS cloud is not big enough : increase it - then closest pt should be used
      E_Int sizeOfCloud = listOfCloudPtsPerVertex.size();
      if (sizeOfCloud < sizeMinOfCloud)
      {
        for (size_t noe = 0; noe < eltsVoisins.size(); noe++)
        {
          E_Int indEltR = eltsVoisins[noe];// index of one elt with indR as vertex
          vector<E_Int>& eltsVoisins2 = cEEN[indEltR];
          for (size_t noe2 = 0; noe2 < eltsVoisins2.size(); noe2++)
          {
            E_Int indEltR2 = eltsVoisins2[noe2];
            vector<E_Int>& projectedPts = listOfCloudPtsPerTri[indEltR2];// projected points on this elt
            for (size_t nop = 0; nop < projectedPts.size(); nop++)
              listOfCloudPtsPerVertex.push_back(projectedPts[nop]);
          }
        }
      }
      sort(listOfCloudPtsPerVertex.begin(), listOfCloudPtsPerVertex.end());
      listOfCloudPtsPerVertex.erase(unique(listOfCloudPtsPerVertex.begin(), listOfCloudPtsPerVertex.end()), 
                                    listOfCloudPtsPerVertex.end()); 
      sizeOfCloud = listOfCloudPtsPerVertex.size();

      if (sizeOfCloud < sizeMinOfCloud) indicesExtrap.push_back(indR);//extrap
      else
      {
        pt[0] = xtRcv[indR]; pt[1] = ytRcv[indR]; pt[2] = ztRcv[indR];
        vector<E_Float> cfLoc(sizeOfCloud);
        axisConst[0] = 0; axisConst[1] = 0; axisConst[2] = 0;
        E_Float x_max =-K_CONST::E_MAX_FLOAT; 
        E_Float y_max =-K_CONST::E_MAX_FLOAT; 
        E_Float z_max =-K_CONST::E_MAX_FLOAT; 
        E_Float x_min = K_CONST::E_MAX_FLOAT; 
        E_Float y_min = K_CONST::E_MAX_FLOAT; 
        E_Float z_min = K_CONST::E_MAX_FLOAT; 
        for (E_Int noindd = 0; noindd < sizeOfCloud; noindd++)
        {
          E_Int indD = listOfCloudPtsPerVertex[noindd];
          x_max = K_FUNC::E_max(xtDnr[indD],x_max);
          y_max = K_FUNC::E_max(ytDnr[indD],y_max);
          z_max = K_FUNC::E_max(ztDnr[indD],z_max);
          x_min = K_FUNC::E_min(xtDnr[indD],x_min);
          y_min = K_FUNC::E_min(ytDnr[indD],y_min);
          z_min = K_FUNC::E_min(ztDnr[indD],z_min);
        }
        if ( K_FUNC::E_abs(x_max-x_min) < 1.e-6 ) axisConst[0] = 1;
        if ( K_FUNC::E_abs(y_max-y_min) < 1.e-6 ) axisConst[1] = 1;
        if ( K_FUNC::E_abs(z_max-z_min) < 1.e-6 ) axisConst[2] = 1;

        E_Int ok = K_INTERP::getInterpCoefMLS(order, dimPb, sizeBasis, pt, xtDnr, ytDnr, ztDnr, 
                                              listOfCloudPtsPerVertex,
                                              radius, axis, axisConst, &cfLoc[0]);
        if ( ok != 1)
        {
          indicesExtrap.push_back(indR);
        }
        else 
        {
          E_Int nbNonZeroCfs = 0;
          for (E_Int nocf = 0; nocf < sizeOfCloud; nocf++)
          {
            if ( K_FUNC::E_abs(cfLoc[nocf]) > K_CONST::E_GEOM_CUTOFF) nbNonZeroCfs++;
          }
          if ( nbNonZeroCfs < 2) 
          { 
            indicesExtrap.push_back(indR);           
          }
          else
          {
            E_Float sumCf = 0.;
            E_Float minCf =  K_CONST::E_MAX_FLOAT;
            E_Float maxCf = -K_CONST::E_MAX_FLOAT;
            for(E_Int noind = 0; noind < sizeOfCloud; noind++)
            {
              sumCf += cfLoc[noind];
              minCf = K_FUNC::E_min(cfLoc[noind],minCf);
              maxCf = K_FUNC::E_max(cfLoc[noind],maxCf);
            }
            if (K_FUNC::E_abs(sumCf-1.)>5.e-2) indicesExtrap.push_back(indR);
            else if ((minCf<-2. || maxCf>2.) && ibm==1) indicesExtrap.push_back(indR); // avoid huge weights for IBMs
	    else if (old == 1) indicesExtrap.push_back(indR); // avoid huge weights for IBMs
            else // interpolate
            {
              for (E_Int posv = 0; posv < nbVars; posv++)
              {
                E_Float* varD = fd->begin(posvarsD[posv]+1);
                E_Float* varR = fieldROut.begin(posvarsR[posv]+1);
                E_Float val = 0.;
                for(E_Int noind = 0; noind < sizeOfCloud; noind++)
                {
                  E_Int indD = listOfCloudPtsPerVertex[noind];
                  val += cfLoc[noind]*varD[indD];                  
                }
                  varR[indR] = val;
              }
            }   
          }          
        }
      }
    }// indR to be interpolated
  }// loop on indR

  //cleanings
  for (E_Int nov =0; nov < nbPtsR; nov++)
    cveR[nov].clear();
  cveR.clear();
  for (E_Int nov = 0; nov < nbEltsR; nov++)
  {
    listOfCloudPtsPerTri[nov].clear();
    cEEN[nov].clear();
  }
  listOfCloudPtsPerTri.clear(); cEEN.clear();

  E_Int nbExtrapPts = indicesExtrap.size();
  E_Float percentageOfExtrap = nbExtrapPts/float(nbPtsR)*100.;
  if (nbExtrapPts > 0)
  {
    printf("Info: projectCloudSolution2Triangle: percentage of extrapolated points = %2.2f%% (" SF_D_ "/" SF_D_ ")\n", percentageOfExtrap, nbExtrapPts, nbPtsR);
    // kdtree of cloud pts
    ArrayAccessor<FldArrayF>* coordAcc = new ArrayAccessor<FldArrayF>(*fd, posxd, posyd, poszd);
    K_SEARCH::KdTree<FldArrayF>* kdt = new K_SEARCH::KdTree<FldArrayF>(*coordAcc);

    for(E_Int noind = 0; noind < nbExtrapPts; noind++)
    {
      E_Int indR = indicesExtrap[noind];
      pt[0] = xtRcv[indR]; pt[1] = ytRcv[indR]; pt[2] = ztRcv[indR];
      E_Int indD = kdt->getClosest(pt);
      for (E_Int posv = 0; posv < nbVars; posv++)
      {
        E_Float* varD = fd->begin(posvarsD[posv]+1);
        E_Float* varR = fieldROut.begin(posvarsR[posv]+1);
        varR[indR] = varD[indD];
      }
    }
    delete kdt; delete coordAcc;
  }
  RELEASESHAREDB(resr, arrayR, fr, cnr);
  RELEASESHAREDB(resd, arrayD, fd, cnd);
  return tpl;
}

// ============================================================================
/*  IN: nuage de points donneur defini par une zone NODE
    IN/OUT: one list of donor indices, one list of donor coefficients */
// ============================================================================
PyObject* K_POST::prepareProjectCloudSolution2Triangle(PyObject* self, PyObject* args)
{
  E_Int dimPb, ibm;
  PyObject *arrayR, *arrayD;
  if (!PYPARSETUPLE_(args, OO_ II_,
                    &arrayD, &arrayR, &dimPb, &ibm))
  {
      return NULL;
  }

  /*-----------------------------------------------*/
  /* Extraction des infos sur le domaine recepteur */
  /*-----------------------------------------------*/
  E_Int imr, jmr, kmr;
  FldArrayF* fr; FldArrayI* cnr;
  char* varStringr; char* eltTyper;
  E_Int resr = K_ARRAY::getFromArray(arrayR, varStringr, fr,
                                     imr, jmr, kmr, cnr, eltTyper, true);
  if (resr != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "projectCloudSolution2Triangle: invalid receptor zone.");
    if (resr == 1) RELEASESHAREDB(resr, arrayR, fr, cnr);
    return NULL;
  }
  if (K_STRING::cmp(eltTyper, "TRI") != 0 && K_STRING::cmp(eltTyper, "BAR") != 0 && 
      K_STRING::cmp(eltTyper, "TRI*") != 0 && K_STRING::cmp(eltTyper, "BAR*") != 0  )
  {
    PyErr_SetString(PyExc_TypeError,
      "projectCloudSolution2Triangle: unstructured receptor zone must be TRI or BAR.");
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    return NULL;
  }

  /*---------------------------------------------*/
  /* Extraction des infos sur le domaine donneur */
  /*---------------------------------------------*/
  E_Int imd, jmd, kmd;
  FldArrayF* fd; FldArrayI* cnd;
  char* varStringd; char* eltTyped;
  E_Int resd = K_ARRAY::getFromArray(arrayD, varStringd, fd,
                                     imd, jmd, kmd, cnd, eltTyped, true);
  if (resd != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "projectCloudSolution2Triangle: 2nd arg is not a valid array.");
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    if (resd == 1) RELEASESHAREDB(resd, arrayD, fd, cnd);
    return NULL;
  }

  if (K_STRING::cmp(eltTyped, "NODE") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
      "projectCloudSolution2Triangle: unstructured donor zone must only be NODE.");
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    RELEASESHAREDB(resd, arrayD, fd, cnd);
    return NULL;
  }
  E_Int posxd = K_ARRAY::isCoordinateXPresent(varStringd);
  E_Int posyd = K_ARRAY::isCoordinateYPresent(varStringd);
  E_Int poszd = K_ARRAY::isCoordinateZPresent(varStringd);
  if (posxd == -1 || posyd == -1 || poszd == -1)
  {
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    RELEASESHAREDB(resd, arrayD, fd, cnd);
    PyErr_SetString(PyExc_TypeError,
                    "projectCloudSolution2Triangle: coordinates not found in donor zone.");
    return NULL;
  }

  E_Int posxr = K_ARRAY::isCoordinateXPresent(varStringr);
  E_Int posyr = K_ARRAY::isCoordinateYPresent(varStringr);
  E_Int poszr = K_ARRAY::isCoordinateZPresent(varStringr);
  if (posxr == -1 || posyr == -1 || poszr == -1)
  {
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    RELEASESHAREDB(resd, arrayD, fd, cnd);
    PyErr_SetString(PyExc_TypeError,
                    "projectCloudSolution2Triangle: coordinates not found in receptor zone.");
    return NULL;
  }
  //Variables a transferer
  E_Int poscd = K_ARRAY::isCellNatureField2Present(varStringd);
  E_Int poscr = K_ARRAY::isCellNatureField2Present(varStringr);
  char* varStringC; // chaine de caractere commune
  E_Int l = strlen(varStringr); varStringC = new char [l+1];
  // les positions demarrent a 1
  vector<E_Int> posvarsD; vector<E_Int> posvarsR;
  K_ARRAY::getPosition(varStringr, varStringd,
                       posvarsR, posvarsD, varStringC);
  delete [] varStringC;
  E_Int sizeVarsD = posvarsD.size();
  E_Int sizeVarsR = posvarsR.size();
  for (E_Int i = 0; i < sizeVarsD; i++) posvarsD[i] -= 1;
  for (E_Int i = 0; i < sizeVarsR; i++) posvarsR[i] -= 1;
  if (poscd != -1) posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), poscd), posvarsD.end());
  if (poscr != -1) posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), poscr), posvarsR.end());
  posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), posxd), posvarsD.end());
  posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), posyd), posvarsD.end());
  posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), poszd), posvarsD.end());
  posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), posxr), posvarsR.end());
  posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), posyr), posvarsR.end());
  posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), poszr), posvarsR.end());
  posxd++; posyd++; poszd++; posxr++; posyr++; poszr++;
  E_Int order = 2; // ordre d'interpolation MLS
  const E_Int sizeBasis = K_FUNC::fact(dimPb+order-1)*1./(K_FUNC::fact(dimPb)*K_FUNC::fact(order-1));

  //1. Determination du(des) triangles contenant le point dnr
  E_Float* xtRcv = fr->begin(posxr);
  E_Float* ytRcv = fr->begin(posyr);
  E_Float* ztRcv = fr->begin(poszr);

  // Creation du kdtree des points de la triangulation
  ArrayAccessor<FldArrayF>* coordAcc = new ArrayAccessor<FldArrayF>(*fr, posxr, posyr, poszr);
  K_SEARCH::KdTree<FldArrayF>* kdt = new K_SEARCH::KdTree<FldArrayF>(*coordAcc);

  // Creation du bboxtree
  E_Int nbEltsR = cnr->getSize();
  E_Int nbPtsR = fr->getSize();

  typedef K_SEARCH::BoundingBox<3>  BBox3DType;
  vector<BBox3DType*> boxes(nbEltsR);// liste des bbox de ts les elements de a2
  K_FLD::FldArrayF bbox(nbEltsR,6);// xmin, ymin, zmin, xmax, ymax, zmax
  K_COMPGEOM::boundingBoxOfUnstrCells(*cnr, xtRcv, ytRcv, ztRcv, bbox);
  E_Float minB[3];  E_Float maxB[3];
  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
  E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
  E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
  for (E_Int et = 0; et < nbEltsR; et++)
  {
    minB[0] = xminp[et]; minB[1] = yminp[et]; minB[2] = zminp[et];
    maxB[0] = xmaxp[et]; maxB[1] = ymaxp[et]; maxB[2] = zmaxp[et];
    boxes[et] = new BBox3DType(minB, maxB);
  }
  // Build the box tree.
  K_SEARCH::BbTree3D bbtree(boxes);
  //projection des points de fd(nuage) sur fr(TRI)
  E_Float xo, yo, zo, rx, ry, rz, rad;
  E_Float pt[3];
  E_Float p0[3]; E_Float p1[3]; E_Float p2[3]; E_Float p[3];
  E_Int indev, indR;
  vector<E_Int> indicesBB; // liste des indices des facettes intersectant la bbox
  E_Int nbPtsD = fd->getSize();

  vector < vector<E_Int> > cveR(nbPtsR);//connectivite vertex/elts
  K_CONNECT::connectEV2VE(*cnr, cveR);

  vector< vector<E_Int> > listOfCloudPtsPerTri(nbEltsR);

  FldArrayI matching(nbPtsR); matching.setAllValuesAt(-1);
  E_Int* matchingP = matching.begin();
  E_Float* xtDnr = fd->begin(posxd);
  E_Float* ytDnr = fd->begin(posyd);
  E_Float* ztDnr = fd->begin(poszd);

  E_Float dist2proj = 0.;

  for (E_Int indD = 0; indD < nbPtsD; indD++)
  {
    // recherche du pt le plus proche P' de P
    pt[0] = xtDnr[indD]; pt[1] = ytDnr[indD]; pt[2] = ztDnr[indD];
    indR = kdt->getClosest(pt);

    // calcul de la bounding box de la sphere de rayon PP'
    rx = pt[0]-xtRcv[indR]; ry = pt[1]-ytRcv[indR]; rz = pt[2]-ztRcv[indR];
    rad = sqrt(rx*rx+ry*ry+rz*rz);
    if (rad < K_CONST::E_GEOM_CUTOFF)//matching
    {matchingP[indR] = indD;}
      
    else //projection 
    {
      minB[0] = pt[0]-rad; minB[1] = pt[1]-rad; minB[2] = pt[2]-rad;
      maxB[0] = pt[0]+rad; maxB[1] = pt[1]+rad; maxB[2] = pt[2]+rad;
      bbtree.getOverlappingBoxes(minB, maxB, indicesBB);
      indev = K_COMPGEOM::projectOrthoPrecond(pt[0], pt[1], pt[2],
                                              xtRcv, ytRcv, ztRcv, indicesBB, *cnr, 
                                              xo, yo, zo, p0, p1, p2, p);
      dist2proj = sqrt((pt[0]-xo)*(pt[0]-xo) + (pt[1]-yo)*(pt[1]-yo) + (pt[2]-zo)*(pt[2]-zo));
      indicesBB.clear();
      if (indev != -1 && !(dist2proj>1.e-6 && ibm==1)) // donor points should be already on the triangle or nearby for IBMs
      {
        listOfCloudPtsPerTri[indev].push_back(indD);
      }
    }
  }
  
  // cleanings
  delete kdt; delete coordAcc;
  E_Int size = boxes.size();
  for (E_Int v = 0; v < size; v++) delete boxes[v];
    
  E_Int sizeMinOfCloud = sizeBasis+1;
  vector< vector<E_Int> > cEEN(nbEltsR);
  K_CONNECT::connectEV2EENbrs(eltTyper, nbPtsR, *cnr, cEEN); 
 
  E_Int axisConst[3]; // =1 si l'axe est une direction constante, 0 sinon
  axisConst[0] = axisConst[1] = axisConst[2] = 0;
  if (dimPb == 2) axisConst[2] = 1;// pas de z
  E_Float radius[3]; // longueurs des axes de l'ellipse
  radius[0] = radius[1] = radius[2] = 1.;
  E_Float axis[9];  // Axes de l'ellipse
  axis[0] = 1.; axis[1] = 0.; axis[2] = 0.;
  axis[3] = 0.; axis[4] = 1.; axis[5] = 0.;
  axis[6] = 0.; axis[7] = 0.; axis[8] = 1.;

  // -----------------------------------------------------------------

  // Donnees d'interpolation
  vector<FldArrayF*> listOfInterpCoefs;
  vector<FldArrayI*> listOfDonorInd; 

  vector<E_Int> indicesExtrap;
  for (E_Int indR = 0; indR < nbPtsR; indR++)
  {
		E_Int indD = matchingP[indR];
    if (indD>-1)
    {
			FldArrayF* donorInterpCoefs = new FldArrayF(1); donorInterpCoefs->setAllValuesAtNull();
			donorInterpCoefs[0] = 1;
			listOfInterpCoefs.push_back(donorInterpCoefs);

			FldArrayI* donorInd = new FldArrayI(1); donorInd->setAllValuesAt(-1);
			donorInd[0] = indD;
			listOfDonorInd.push_back(donorInd);
    }
		else
		{
			vector<E_Int>& eltsVoisins = cveR[indR];
			// create the cloud of donor points
			vector<E_Int> listOfCloudPtsPerVertex;
			for (size_t noe = 0; noe < eltsVoisins.size(); noe++)
			{
				E_Int indEltR = eltsVoisins[noe];// index of one elt with indR as vertex
				vector<E_Int>& projectedPts = listOfCloudPtsPerTri[indEltR];// projected points on this elt
				for (size_t nop = 0; nop < projectedPts.size(); nop++)
					listOfCloudPtsPerVertex.push_back(projectedPts[nop]);
			}
			sort(listOfCloudPtsPerVertex.begin(), listOfCloudPtsPerVertex.end());
			listOfCloudPtsPerVertex.erase(unique(listOfCloudPtsPerVertex.begin(), listOfCloudPtsPerVertex.end()), 
																		listOfCloudPtsPerVertex.end());
			// if the MLS cloud is not big enough : increase it - then closest pt should be used
			E_Int sizeOfCloud = listOfCloudPtsPerVertex.size();
			if (sizeOfCloud < sizeMinOfCloud)
			{
				for (size_t noe = 0; noe < eltsVoisins.size(); noe++)
				{
					E_Int indEltR = eltsVoisins[noe];// index of one elt with indR as vertex
					vector<E_Int>& eltsVoisins2 = cEEN[indEltR];
					for (size_t noe2 = 0; noe2 < eltsVoisins2.size(); noe2++)
					{
						E_Int indEltR2 = eltsVoisins2[noe2];
						vector<E_Int>& projectedPts = listOfCloudPtsPerTri[indEltR2];// projected points on this elt
						for (size_t nop = 0; nop < projectedPts.size(); nop++)
							listOfCloudPtsPerVertex.push_back(projectedPts[nop]);
					}
				}
			}
			sort(listOfCloudPtsPerVertex.begin(), listOfCloudPtsPerVertex.end());
			listOfCloudPtsPerVertex.erase(unique(listOfCloudPtsPerVertex.begin(), listOfCloudPtsPerVertex.end()), 
																		listOfCloudPtsPerVertex.end()); 
			sizeOfCloud = listOfCloudPtsPerVertex.size();

			FldArrayF* donorInterpCoefs = new FldArrayF(sizeOfCloud); donorInterpCoefs->setAllValuesAtNull();
			listOfInterpCoefs.push_back(donorInterpCoefs);
			FldArrayI* donorInd = new FldArrayI(sizeOfCloud); donorInd->setAllValuesAt(-1);
			listOfDonorInd.push_back(donorInd);

			if (sizeOfCloud < sizeMinOfCloud) indicesExtrap.push_back(indR);//extrap
			else
			{
				pt[0] = xtRcv[indR]; pt[1] = ytRcv[indR]; pt[2] = ztRcv[indR];
				vector<E_Float> cfLoc(sizeOfCloud);
				axisConst[0] = 0; axisConst[1] = 0; axisConst[2] = 0;
				E_Float x_max =-K_CONST::E_MAX_FLOAT; 
				E_Float y_max =-K_CONST::E_MAX_FLOAT; 
				E_Float z_max =-K_CONST::E_MAX_FLOAT; 
				E_Float x_min = K_CONST::E_MAX_FLOAT; 
				E_Float y_min = K_CONST::E_MAX_FLOAT; 
				E_Float z_min = K_CONST::E_MAX_FLOAT; 
				for (E_Int noindd = 0; noindd < sizeOfCloud; noindd++)
				{
					E_Int indD = listOfCloudPtsPerVertex[noindd];
					x_max = K_FUNC::E_max(xtDnr[indD],x_max);
					y_max = K_FUNC::E_max(ytDnr[indD],y_max);
					z_max = K_FUNC::E_max(ztDnr[indD],z_max);
					x_min = K_FUNC::E_min(xtDnr[indD],x_min);
					y_min = K_FUNC::E_min(ytDnr[indD],y_min);
					z_min = K_FUNC::E_min(ztDnr[indD],z_min);
				}
				if ( K_FUNC::E_abs(x_max-x_min) < 1.e-6 ) axisConst[0] = 1;
				if ( K_FUNC::E_abs(y_max-y_min) < 1.e-6 ) axisConst[1] = 1;
				if ( K_FUNC::E_abs(z_max-z_min) < 1.e-6 ) axisConst[2] = 1;

				E_Int ok = K_INTERP::getInterpCoefMLS(order, dimPb, sizeBasis, pt, xtDnr, ytDnr, ztDnr, 
																							listOfCloudPtsPerVertex,
																							radius, axis, axisConst, &cfLoc[0]);
				if ( ok != 1)
				{
					indicesExtrap.push_back(indR);
				}
				else 
				{
					E_Int nbNonZeroCfs = 0;
					for (E_Int nocf = 0; nocf < sizeOfCloud; nocf++)
					{
						if ( K_FUNC::E_abs(cfLoc[nocf]) > K_CONST::E_GEOM_CUTOFF) nbNonZeroCfs++;
					}
					if ( nbNonZeroCfs < 2) indicesExtrap.push_back(indR);           
					else
					{
						E_Float sumCf = 0.;
						E_Float minCf =  K_CONST::E_MAX_FLOAT;
						E_Float maxCf = -K_CONST::E_MAX_FLOAT;
						for(E_Int noind = 0; noind < sizeOfCloud; noind++)
						{
							sumCf += cfLoc[noind];
							minCf = K_FUNC::E_min(cfLoc[noind],minCf);
							maxCf = K_FUNC::E_max(cfLoc[noind],maxCf);
						}
						if (K_FUNC::E_abs(sumCf-1.)>5.e-2) indicesExtrap.push_back(indR);
						else if ((minCf<-2. || maxCf>2.) && ibm==1) indicesExtrap.push_back(indR); // avoid huge weights for IBMs
						else // interpolate
						{
							E_Float* donorCf = listOfInterpCoefs[indR]->begin();
							E_Int* donorInd = listOfDonorInd[indR]->begin();
							for(E_Int noind = 0; noind < sizeOfCloud; noind++)
							{
								donorCf[noind] = cfLoc[noind];
								donorInd[noind] = listOfCloudPtsPerVertex[noind];      
							}
						}
					}        
				}
			}
		}
  }// loop on indR

  //cleanings
  for (E_Int nov =0; nov < nbPtsR; nov++)
    cveR[nov].clear();
  cveR.clear();
  for (E_Int nov = 0; nov < nbEltsR; nov++)
  {
    listOfCloudPtsPerTri[nov].clear();
    cEEN[nov].clear();
  }
  listOfCloudPtsPerTri.clear(); cEEN.clear();

  E_Int nbExtrapPts = indicesExtrap.size();
  E_Float percentageOfExtrap = nbExtrapPts/float(nbPtsR)*100.;
  if (nbExtrapPts > 0)
  {
    printf("Info: projectCloudSolution2Triangle: percentage of extrapolated points = %2.2f%% (" SF_D_ "/" SF_D_ ")\n", percentageOfExtrap, nbExtrapPts, nbPtsR);
    // kdtree of cloud pts
    ArrayAccessor<FldArrayF>* coordAcc = new ArrayAccessor<FldArrayF>(*fd, posxd, posyd, poszd);
    K_SEARCH::KdTree<FldArrayF>* kdt = new K_SEARCH::KdTree<FldArrayF>(*coordAcc);

    for(E_Int noind = 0; noind < nbExtrapPts; noind++)
    {
      E_Int indR = indicesExtrap[noind];
      listOfInterpCoefs[indR]->resize(1);
      listOfDonorInd[indR]->resize(1);
      E_Float* donorCf = listOfInterpCoefs[indR]->begin();
      E_Int* donorInd = listOfDonorInd[indR]->begin();
      
      pt[0] = xtRcv[indR]; pt[1] = ytRcv[indR]; pt[2] = ztRcv[indR];
      E_Int indD = kdt->getClosest(pt);
      donorCf[0] = 1;
      donorInd[0] = indD;
    }
    delete kdt; delete coordAcc;
  }

  PyObject * PyListCoefficients = PyList_New(0);
  PyObject * PyListIndices = PyList_New(0);
  for (E_Int indR = 0; indR < nbPtsR; indR++)
  {   
    PyObject* fout = K_NUMPY::buildNumpyArray(*listOfInterpCoefs[indR],1);
    PyList_Append(PyListCoefficients, fout);  Py_DECREF(fout);
    delete listOfInterpCoefs[indR];

    PyObject* donorIndOut = K_NUMPY::buildNumpyArray(*listOfDonorInd[indR],1);
    PyList_Append(PyListIndices, donorIndOut); Py_DECREF(donorIndOut);
    delete listOfDonorInd[indR];
  }

  PyObject* tpl = Py_BuildValue("[OO]", PyListIndices, PyListCoefficients);

  RELEASESHAREDB(resr, arrayR, fr, cnr);
  RELEASESHAREDB(resd, arrayD, fd, cnd);
  return tpl;
}

// ============================================================================
/*  IN: nuage de points donneur defini par une zone NODE
    IN: offset, interpDonor & interpCoef (created with prepareProjectCloudSolution)
    IN/OUT: zone surfacique triangulaire receptrice */
// ============================================================================
PyObject* K_POST::projectCloudSolution2TriangleWithInterpData(PyObject* self, PyObject* args)
{
  E_Int dimPb;
  PyObject *arrayR, *arrayD;
  PyObject *pyOffset, *pyInterpDonor, *pyInterpCoef;
  if (!PYPARSETUPLE_(args, OO_ OOO_ I_,
                    &arrayD, &arrayR, 
                    &pyOffset, &pyInterpDonor, &pyInterpCoef,
                    &dimPb))
  {
      return NULL;
  }

	/*-----------------------------------------------*/
  /* Extraction des donnees d interpolation        */
  /*-----------------------------------------------*/
	FldArrayI* offsetI;
	K_NUMPY::getFromNumpyArray(pyOffset, offsetI, true);
  E_Int* offset  = offsetI->begin();

	FldArrayI* interpDonorI;
  K_NUMPY::getFromNumpyArray(pyInterpDonor, interpDonorI, true);
  E_Int* interpDonor  = interpDonorI->begin();

	FldArrayF* interpCoefF;
  K_NUMPY::getFromNumpyArray(pyInterpCoef, interpCoefF, true);
  E_Float* interpCoef  = interpCoefF->begin();

  /*-----------------------------------------------*/
  /* Extraction des infos sur le domaine recepteur */
  /*-----------------------------------------------*/
  E_Int imr, jmr, kmr;
  FldArrayF* fr; FldArrayI* cnr;
  char* varStringr; char* eltTyper;
  E_Int resr = K_ARRAY::getFromArray(arrayR, varStringr, fr,
                                     imr, jmr, kmr, cnr, eltTyper, true);
  if (resr != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "projectCloudSolution2Triangle: invalid receptor zone.");
    if (resr == 1) RELEASESHAREDB(resr, arrayR, fr, cnr);
    return NULL;
  }
  if (K_STRING::cmp(eltTyper, "TRI") != 0 && K_STRING::cmp(eltTyper, "BAR") != 0 && 
      K_STRING::cmp(eltTyper, "TRI*") != 0 && K_STRING::cmp(eltTyper, "BAR*") != 0  )
  {
    PyErr_SetString(PyExc_TypeError,
      "projectCloudSolution2Triangle: unstructured receptor zone must be TRI or BAR.");
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    return NULL;
  }

  /*---------------------------------------------*/
  /* Extraction des infos sur le domaine donneur */
  /*---------------------------------------------*/
  E_Int imd, jmd, kmd;
  FldArrayF* fd; FldArrayI* cnd;
  char* varStringd; char* eltTyped;
  E_Int resd = K_ARRAY::getFromArray(arrayD, varStringd, fd,
                                     imd, jmd, kmd, cnd, eltTyped, true);
  if (resd != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "projectCloudSolution2Triangle: 2nd arg is not a valid array.");
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    if (resd == 1) RELEASESHAREDB(resd, arrayD, fd, cnd);
    return NULL;
  }

  if (K_STRING::cmp(eltTyped, "NODE") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
      "projectCloudSolution2Triangle: unstructured donor zone must only be NODE.");
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    RELEASESHAREDB(resd, arrayD, fd, cnd);
    return NULL;
  }
  E_Int posxd = K_ARRAY::isCoordinateXPresent(varStringd);
  E_Int posyd = K_ARRAY::isCoordinateYPresent(varStringd);
  E_Int poszd = K_ARRAY::isCoordinateZPresent(varStringd);
  if (posxd == -1 || posyd == -1 || poszd == -1)
  {
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    RELEASESHAREDB(resd, arrayD, fd, cnd);
    PyErr_SetString(PyExc_TypeError,
                    "projectCloudSolution2Triangle: coordinates not found in donor zone.");
    return NULL;
  }

  E_Int posxr = K_ARRAY::isCoordinateXPresent(varStringr);
  E_Int posyr = K_ARRAY::isCoordinateYPresent(varStringr);
  E_Int poszr = K_ARRAY::isCoordinateZPresent(varStringr);
  if (posxr == -1 || posyr == -1 || poszr == -1)
  {
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    RELEASESHAREDB(resd, arrayD, fd, cnd);
    PyErr_SetString(PyExc_TypeError,
                    "projectCloudSolution2Triangle: coordinates not found in receptor zone.");
    return NULL;
  }
  //Variables a transferer
  E_Int poscd = K_ARRAY::isCellNatureField2Present(varStringd);
  E_Int poscr = K_ARRAY::isCellNatureField2Present(varStringr);
  char* varStringC; // chaine de caractere commune
  E_Int l = strlen(varStringr); varStringC = new char [l+1];
  // les positions demarrent a 1
  vector<E_Int> posvarsD; vector<E_Int> posvarsR;
  K_ARRAY::getPosition(varStringr, varStringd,
                       posvarsR, posvarsD, varStringC);
  delete [] varStringC;
  E_Int sizeVarsD = posvarsD.size();
  E_Int sizeVarsR = posvarsR.size();
  for (E_Int i = 0; i < sizeVarsD; i++) posvarsD[i] -= 1;
  for (E_Int i = 0; i < sizeVarsR; i++) posvarsR[i] -= 1;
  if (poscd != -1) posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), poscd), posvarsD.end());
  if (poscr != -1) posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), poscr), posvarsR.end());
  posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), posxd), posvarsD.end());
  posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), posyd), posvarsD.end());
  posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), poszd), posvarsD.end());
  posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), posxr), posvarsR.end());
  posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), posyr), posvarsR.end());
  posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), poszr), posvarsR.end());
  posxd++; posyd++; poszd++; posxr++; posyr++; poszr++;
  E_Int nbVars = posvarsD.size();

  E_Int crsize = cnr->getSize()*cnr->getNfld();
  PyObject* tpl = K_ARRAY::buildArray(fr->getNfld(), varStringr,
    fr->getSize(), cnr->getSize(),-1, eltTyper, false, crsize);
  E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
  K_KCORE::memcpy__(cnnp, cnr->begin(), cnr->getSize()*cnr->getNfld());
  E_Float* ptrFieldOut = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fieldROut(fr->getSize(), fr->getNfld(), ptrFieldOut, true);
  fieldROut = *fr;

  E_Int nbPtsR = fr->getSize();

  for (E_Int indR = 0; indR < nbPtsR; indR++)
  {
		E_Int sizeOfCloud = offset[indR+1]-offset[indR];
		E_Int off = offset[indR];
		for (E_Int posv = 0; posv < nbVars; posv++)
		{
			E_Float* varD = fd->begin(posvarsD[posv]+1);
			E_Float* varR = fieldROut.begin(posvarsR[posv]+1);
			E_Float val = 0.;
			for(E_Int noind = 0; noind < sizeOfCloud; noind++)
			{
				E_Int indD = interpDonor[noind+off];
				val += interpCoef[noind+off]*varD[indD];                  
			}
			varR[indR] = val;
		}
  }// loop on indR

  RELEASESHAREDB(resr, arrayR, fr, cnr);
  RELEASESHAREDB(resd, arrayD, fd, cnd);

	RELEASESHAREDN(pyOffset, offsetI);
	RELEASESHAREDN(pyInterpDonor, interpDonorI);
	RELEASESHAREDN(pyInterpCoef, interpCoefF);

  return tpl;
}
