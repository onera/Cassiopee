/*
    Copyright 2020 Onera.

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

// Computes MLS coefficients for the projection of IBM solution on a triangulated surface

# include "connector.h"
# include "Search/KdTree.h"
# include "Search/BbTree.h"
# include "Fld/ArrayAccessor.h"
using namespace K_FLD;
using namespace std;

#define CLOUDMAX 400

// ============================================================================
/*  IN: nuage de points donneur defini par une zone NODE
    IN : zone surfacique triangulaire receptrice 
    OUT:  donnees d interpolation du nuage sur la surface */
// ============================================================================
PyObject* K_CONNECTOR::setInterpData_IBMWall(PyObject* self, PyObject* args)
{
  E_Int dimPb, order;
  PyObject *arrayR, *arrayD;
  if (!PYPARSETUPLEI(args, "OOll", "OOii", &arrayD, &arrayR, &dimPb, &order))
    return NULL;
  /*-----------------------------------------------*/
  /* Extraction des infos sur le domaine recepteur */
  /*-----------------------------------------------*/
  E_Int imr, jmr, kmr;
  FldArrayF* fr; FldArrayI* cnr;
  char* varStringr; char* eltTyper;
  E_Int resr = K_ARRAY::getFromArray(arrayR, varStringr, fr,
                                      imr, jmr, kmr, cnr, eltTyper);
  if (resr != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setInterpData_IBMWall: invalid receptor zone.");
    return NULL;
  }

  if (K_STRING::cmp(eltTyper, "TRI") != 0 && K_STRING::cmp(eltTyper, "BAR") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setInterpData_IBMWall: unstructured receptor zone must be TRI or BAR.");
    return NULL;
  }

  E_Int posxr = K_ARRAY::isCoordinateXPresent(varStringr);
  E_Int posyr = K_ARRAY::isCoordinateYPresent(varStringr);
  E_Int poszr = K_ARRAY::isCoordinateZPresent(varStringr);
  if (posxr == -1 || posyr == -1 || poszr == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setInterpData_IBMWall: coordinates not found in receptor zone.");
    return NULL;
  }
  posxr++; posyr++; poszr++;

  /*------------------------------------------------*/
  /* Extraction des infos sur les domaines donneurs */
  /*------------------------------------------------*/
  E_Int imd, jmd, kmd;
  FldArrayF* fd; FldArrayI* cnd;
  char* varStringd; char* eltTyped;
  K_ARRAY::getFromArray(arrayD, varStringd, fd,
                        imd, jmd, kmd, cnd, eltTyped);

  E_Int posxd = K_ARRAY::isCoordinateXPresent(varStringd);
  E_Int posyd = K_ARRAY::isCoordinateYPresent(varStringd);
  E_Int poszd = K_ARRAY::isCoordinateZPresent(varStringd);
  if (posxd == -1 || posyd == -1 || poszd == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setInterpData_IBMWall: coordinates not found in donor zone.");
    return NULL;
  }
  posxd++; posyd++; poszd++;

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

  for (E_Int indD = 0; indD < nbPtsD; indD++)
  {
    // printf(" indD %d \n", indD);

    // recherche du pt le plus proche P' de P
    pt[0] = xtDnr[indD]; pt[1] = ytDnr[indD]; pt[2] = ztDnr[indD];
    indR = kdt->getClosest(pt);
    // printf(" indR %d \n", indR);

    // calcul de la bounding box de la sphere de rayon PP'
    rx = pt[0]-xtRcv[indR]; ry = pt[1]-ytRcv[indR]; rz = pt[2]-ztRcv[indR];
    rad = sqrt(rx*rx+ry*ry+rz*rz);
    if (rad < K_CONST::E_GEOM_CUTOFF) //matching
    {matchingP[indR] = indD;}
      
    else //projection 
    {
      minB[0] = pt[0]-rad; minB[1] = pt[1]-rad; minB[2] = pt[2]-rad;
      maxB[0] = pt[0]+rad; maxB[1] = pt[1]+rad; maxB[2] = pt[2]+rad;
      bbtree.getOverlappingBoxes(minB, maxB, indicesBB);
      indev = K_COMPGEOM::projectOrthoPrecond(pt[0], pt[1], pt[2],
                                              xtRcv, ytRcv, ztRcv, indicesBB, *cnr, xo, yo, zo);
      indicesBB.clear();
      if (indev != -1)
      {
        // E_Float dist = (xo-pt[0])*(xo-pt[0])+(yo-pt[1])*(yo-pt[1])+(zo-pt[2])*(zo-pt[2]);
        // // if ( K_FUNC::E_abs(xtRcv[indev]-0.750889)< 1.e-4 && K_FUNC::E_abs(ytRcv[indR]-0.0305362)<1.e-4) 
        // printf("Pt dans le nuage : %g %g %g | projete %g %g %g , size = %d | dist = %g \n", pt[0], pt[1], pt[2], xo,yo,zo,
        //        listOfCloudPtsPerTri[indev].size(), dist);
        listOfCloudPtsPerTri[indev].push_back(indD);
      }
      else printf(" indev = %d\n", indev);
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
  
  FldArrayI* indicesR = new FldArrayI(nbPtsR); indicesR->setAllValuesAtNull();
  FldArrayI* indicesD = new FldArrayI(nbPtsR*(CLOUDMAX+1)); indicesD->setAllValuesAt(-1);
  FldArrayI* itypes   = new FldArrayI(nbPtsR); itypes->setAllValuesAtNull();
  FldArrayF* icoefs   = new FldArrayF(nbPtsR*CLOUDMAX); icoefs->setAllValuesAtNull();
  E_Int*   ptrIndicesD = indicesD->begin();
  E_Float* ptrCoefs = icoefs->begin();
  E_Int* ptrIndicesR = indicesR->begin();  
  E_Int sizeIndDnr=0; E_Int sizeCf=0;

  vector<E_Int> indicesExtrap;
  E_Int noindR=0;
  for (E_Int indR = 0; indR < nbPtsR; indR++)
  {
    E_Int indD = matchingP[indR];

    if (indD>-1) // matching
    {
      ptrIndicesD[sizeIndDnr] = 1;    sizeIndDnr++;
      ptrIndicesD[sizeIndDnr] = indD; sizeIndDnr++;
      ptrCoefs[sizeCf]=1.; sizeCf++;
      ptrIndicesR[noindR]=indR; noindR++;
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
      if ( sizeOfCloud < sizeMinOfCloud)
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
      if ( indR == 899 )
        printf(" xR = %g %g %g | sizeOfCloud = %d \n", xtRcv[indR],ytRcv[indR], ztRcv[indR], sizeOfCloud);
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
        // if ( listOfCloudPtsPerVertex.size() > CLOUDMAX) printf(" size = %d \n", listOfCloudPtsPerVertex.size() );
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
            //check if extrapolated first
            E_Float sumCf = 0.;
            for(E_Int noind = 0; noind < sizeOfCloud; noind++)
              sumCf += cfLoc[noind];
            if ( K_FUNC::E_abs(sumCf-1.)>5.e-2) indicesExtrap.push_back(indR);
            else // not extrapolated -> interpolated
            {
              // Nb de donneurs ds la molecule + indices des donneurs
              ptrIndicesD[sizeIndDnr] = sizeOfCloud; sizeIndDnr++;
              ptrIndicesR[noindR]=indR; noindR++;

              for (E_Int noind = 0; noind < sizeOfCloud; noind++)
              {
                E_Int indD = listOfCloudPtsPerVertex[noind];
                ptrIndicesD[sizeIndDnr] = indD; sizeIndDnr++;                 
                ptrCoefs[sizeCf] = cfLoc[noind]; sizeCf++;
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
  if (nbExtrapPts > 0)
  {
    // kdtree of cloud pts
    ArrayAccessor<FldArrayF>* coordAcc = new ArrayAccessor<FldArrayF>(*fd, posxd, posyd, poszd);
    K_SEARCH::KdTree<FldArrayF>* kdt = new K_SEARCH::KdTree<FldArrayF>(*coordAcc);

    for(E_Int noind = 0; noind < nbExtrapPts; noind++)
    {
      E_Int indR = indicesExtrap[noind];
      pt[0] = xtRcv[indR]; pt[1] = ytRcv[indR]; pt[2] = ztRcv[indR];
      E_Int indD = kdt->getClosest(pt);
      ptrIndicesD[sizeIndDnr] = 1; sizeIndDnr++;
      ptrIndicesD[sizeIndDnr] = indD; sizeIndDnr++;
      ptrCoefs[sizeCf] = 1.; sizeCf++;
      ptrIndicesR[noindR]=indR; noindR++;
    }
    delete kdt; delete coordAcc;
  }
  indicesD->resize(sizeIndDnr); icoefs->resize(sizeCf);

  /*----------------------------------------------------------*/
  /* Output Python objects: creation */
  /*----------------------------------------------------------*/  
  PyObject* PyCellIndicesR = K_NUMPY::buildNumpyArray(*indicesR,1); delete indicesR;
  PyObject* PyCellIndicesD = K_NUMPY::buildNumpyArray(*indicesD,1); delete indicesD;
  PyObject* PyInterpTypes = K_NUMPY::buildNumpyArray(*itypes,1); delete itypes;
  PyObject* PyCoefficients = K_NUMPY::buildNumpyArray(*icoefs,1); delete icoefs;

  PyObject* tpl = Py_BuildValue("[OOOO]", PyCellIndicesR, PyCellIndicesD,  PyInterpTypes, PyCoefficients);
  Py_DECREF(PyCellIndicesD); Py_DECREF(PyInterpTypes); Py_DECREF(PyCoefficients);  Py_DECREF(PyCellIndicesR); 
  return tpl;
}
