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
// Computes MLS coefficients for the projection of IBM solution on a triangulated surface

# include "connector.h"
# include "Nuga/include/KdTree.h"
# include "Nuga/include/BbTree.h"
# include "Nuga/include/ArrayAccessor.h"
using namespace K_FLD;
using namespace std;

#define CLOUDMAX 100

// ============================================================================
/*  IN: nuage de points donneur defini par une zone NODE
    IN : zone surfacique triangulaire receptrice 
    OUT:  donnees d interpolation du nuage sur la surface */
// ============================================================================
PyObject* K_CONNECTOR::setInterpData_IBMWall(PyObject* self, PyObject* args)
{
  E_Int dimPb, order;
  PyObject *arrayR, *arraysD;
  if (!PYPARSETUPLE_(args, OO_ II_,
                     &arraysD, &arrayR, &dimPb, &order))
    return NULL;
  
  if (PyList_Check(arraysD) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpDataIBMWall: 1st argument must be a list of arrays.");
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
    if (resr == 1) RELEASESHAREDS(arrayR, fr);
    PyErr_SetString(PyExc_TypeError,
                    "setInterpData_IBMWall: invalid receptor zone.");
    return NULL;
  }

  if (K_STRING::cmp(eltTyper, "TRI") != 0)//&& K_STRING::cmp(eltTyper, "BAR") != 0)
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
    RELEASESHAREDU(arrayR, fr, cnr);
    PyErr_SetString(PyExc_TypeError,
                    "setInterpData_IBMWall: coordinates not found in receptor zone.");
    return NULL;
  }
  posxr++; posyr++; poszr++;

  /*------------------------------------------------*/
  /* Extraction des infos sur les domaines donneurs */
  /*------------------------------------------------*/
  vector<E_Int> resl;  vector<char*> varStringd;
  vector<FldArrayF*> vectOfDnrZones;
  vector<void*> a2; //ni,nj,nk ou cnt en NS
  vector<void*> a3; //eltType en NS
  vector<void*> a4;
  vector<PyObject*> objs;
  E_Bool skipNoCoord = true;  E_Bool skipStructured = false;
  E_Bool skipUnstructured = false;  E_Bool skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(arraysD, resl, varStringd, vectOfDnrZones,
                                      a2, a3, a4, objs,  
                                      skipDiffVars, skipNoCoord, skipStructured,
                                      skipUnstructured, true);
  E_Int nzonesD = objs.size();
  
  vector<E_Int> posxtd; vector<E_Int> posytd; vector<E_Int> posztd;
  E_Int nptsTotD = 0;
  E_Int nptsMaxD = 0;
  for (E_Int i = 0; i < nzonesD; i++)
  {
    E_Int posxd = K_ARRAY::isCoordinateXPresent(varStringd[i]);
    E_Int posyd = K_ARRAY::isCoordinateYPresent(varStringd[i]);
    E_Int poszd = K_ARRAY::isCoordinateZPresent(varStringd[i]);
    posxd++; posyd++; poszd++;      
    posxtd.push_back(posxd); posytd.push_back(posyd); posztd.push_back(poszd);
    nptsTotD += vectOfDnrZones[i]->getSize();
    nptsMaxD = K_FUNC::E_max(nptsMaxD,vectOfDnrZones[i]->getSize());
  }
  /* End of checks and extractions */

  E_Float* xtRcv = fr->begin(posxr);
  E_Float* ytRcv = fr->begin(posyr);
  E_Float* ztRcv = fr->begin(poszr);

  // Creation of the kdtree of vertices of the tesselation 
  ArrayAccessor<FldArrayF>* coordAcc = new ArrayAccessor<FldArrayF>(*fr, posxr, posyr, poszr);
  K_SEARCH::KdTree<FldArrayF>* kdt = new K_SEARCH::KdTree<FldArrayF>(*coordAcc);
  E_Int nbPtsR = fr->getSize();

  //projection des points de fd(nuage) sur fr(TRI)
  E_Float pt[3];
  E_Int indR;
  vector < vector<E_Int> > cveR(nbPtsR);//connectivite vertex/elts
  K_CONNECT::connectEV2VE(*cnr, cveR);
  vector< vector<E_Int> > listOfDnrIndPerTri(nbPtsR);
  vector< vector<E_Int> > listOfDnrNoZPerTri(nbPtsR);
  nzonesD = vectOfDnrZones.size();
  E_Float rad2;
  E_Int* cn1 = cnr->begin(1);
  E_Int* cn2 = cnr->begin(2);

  vector<E_Int> indicesExtrap;
  //creation du kdtree du nuage de points donneurs
  FldArrayF dnrCoords(nptsTotD,3);
  E_Float* xdnr = dnrCoords.begin(1);
  E_Float* ydnr = dnrCoords.begin(2);
  E_Float* zdnr = dnrCoords.begin(3);
  E_Int noindg = 0;
  FldArrayI indicesGlob(nptsMaxD*nzonesD); indicesGlob.setAllValuesAt(-1);
  for (E_Int nozd =  0; nozd < nzonesD; nozd++)
  {
    E_Float* xd = vectOfDnrZones[nozd]->begin(posxtd[nozd]);
    E_Float* yd = vectOfDnrZones[nozd]->begin(posytd[nozd]);
    E_Float* zd = vectOfDnrZones[nozd]->begin(posztd[nozd]);
    for (E_Int noind = 0; noind < vectOfDnrZones[nozd]->getSize(); noind++)
    {
      xdnr[noindg]=xd[noind]; ydnr[noindg]=yd[noind]; zdnr[noindg]=zd[noind];
      indicesGlob[noindg] = noind + nozd*nptsMaxD;
      noindg++;
    }
  }

  ArrayAccessor<FldArrayF>* coordAccD = new ArrayAccessor<FldArrayF>(dnrCoords, 1,2,3);
  K_SEARCH::KdTree<FldArrayF>* kdtD = new K_SEARCH::KdTree<FldArrayF>(*coordAccD);
  if ( strcmp(eltTyper,"TRI") == 0)
  {
    E_Int* cn3 = cnr->begin(3);

    for (E_Int nozd = 0; nozd < nzonesD; nozd++)
    {
      E_Int nptsD = vectOfDnrZones[nozd]->getSize();
      E_Float* xtDnr = vectOfDnrZones[nozd]->begin(posxtd[nozd]);
      E_Float* ytDnr = vectOfDnrZones[nozd]->begin(posytd[nozd]);
      E_Float* ztDnr = vectOfDnrZones[nozd]->begin(posztd[nozd]);
      for (E_Int indD = 0; indD < nptsD; indD++)
      {
        pt[0] = xtDnr[indD]; pt[1] = ytDnr[indD]; pt[2] = ztDnr[indD];
        indR = kdt->getClosest(pt, rad2);

        // search for the max length of the triangles containing indR 
        vector<E_Int>& eltsVoisins = cveR[indR];
        E_Float projTol2 = 0.;
        for (size_t noetv = 0; noetv < eltsVoisins.size(); noetv++)
        {
          E_Int etv = eltsVoisins[noetv];
          E_Int ind1 = cn1[etv]-1; E_Int ind2 = cn2[etv]-1; E_Int ind3 = cn3[etv]-1;
          E_Float lx = xtRcv[ind1]-xtRcv[ind2];
          E_Float ly = ytRcv[ind1]-ytRcv[ind2];
          E_Float lz = ztRcv[ind1]-ztRcv[ind2];
          E_Float d2 = lx*lx+ly*ly+lz*lz;
          lx = xtRcv[ind3]-xtRcv[ind2];
          ly = ytRcv[ind3]-ytRcv[ind2];
          lz = ztRcv[ind3]-ztRcv[ind2];
          d2 = max(d2, lx*lx+ly*ly+lz*lz);
          lx = xtRcv[ind3]-xtRcv[ind1];
          ly = ytRcv[ind3]-ytRcv[ind1];
          lz = ztRcv[ind3]-ztRcv[ind1];
          d2 = max(d2, lx*lx+ly*ly+lz*lz);        
          projTol2 = max(d2,projTol2);
        }

        // test if the cloud point is on a triangle to which the closest vertex belongs
        // printf(" ind = %d, rad2 = %g %g \n", indD, rad2, projTol2);
        if (rad2 < projTol2)
        {
          listOfDnrIndPerTri[indR].push_back(indD);
          listOfDnrNoZPerTri[indR].push_back(nozd);
        }
      }
    }
  }
  else // receptor is a BAR
  {
    for (E_Int nozd = 0; nozd < nzonesD; nozd++)
    {
      E_Int nptsD = vectOfDnrZones[nozd]->getSize();
      E_Float* xtDnr = vectOfDnrZones[nozd]->begin(posxtd[nozd]);
      E_Float* ytDnr = vectOfDnrZones[nozd]->begin(posytd[nozd]);
      E_Float* ztDnr = vectOfDnrZones[nozd]->begin(posztd[nozd]);
      for (E_Int indD = 0; indD < nptsD; indD++)
      {
        pt[0] = xtDnr[indD]; pt[1] = ytDnr[indD]; pt[2] = ztDnr[indD];
        indR = kdt->getClosest(pt, rad2);

        // search for the max length of the triangles containing indR 
        vector<E_Int>& eltsVoisins = cveR[indR];
        E_Float projTol2 = 0.;
        for (size_t noetv = 0; noetv < eltsVoisins.size(); noetv++)
        {
          E_Int etv = eltsVoisins[noetv];
          E_Int ind1 = cn1[etv]-1; E_Int ind2 = cn2[etv]-1; 
          E_Float lx = xtRcv[ind1]-xtRcv[ind2];
          E_Float ly = ytRcv[ind1]-ytRcv[ind2];
          E_Float lz = ztRcv[ind1]-ztRcv[ind2];
          E_Float d2 = lx*lx+ly*ly+lz*lz;
          projTol2 = max(d2,projTol2);
        }

        // test if the cloud point is on a triangle to which the closest vertex belongs
        // printf(" ind = %d, rad2 = %g %g \n", indD, rad2, projTol2);
        if (rad2 < projTol2)
        {
          listOfDnrIndPerTri[indR].push_back(indD);
          listOfDnrNoZPerTri[indR].push_back(nozd);
        }
      }
    }
  }
  // cleanings
  delete kdt; delete coordAcc;

  /* MLS interpolation data for each vertex indR*/
  const E_Int sizeBasis = K_FUNC::fact(dimPb+order-1)*1./(K_FUNC::fact(dimPb)*K_FUNC::fact(order-1));  
  // printf(" MLS interpolation of order %d: minimum number of cloud points = %d\n", order, sizeBasis+1);

  E_Int sizeMinOfCloud = sizeBasis+1;
  E_Int axisConst[3]; // =1 si l'axe est une direction constante, 0 sinon
  E_Float radius[3]; // longueurs des axes de l'ellipse
  E_Float axis[9];  // Axes de l'ellipse

  vector<FldArrayF*> listOfInterpCoefs;
  vector<FldArrayI*> listOfDnrIndices;
  vector<FldArrayI*> listOfRcvIndices;
  vector<FldArrayI*> listOfDnrTypes;
  vector<E_Int> posPtrDnrIndices(nzonesD);
  vector<E_Int> posPtrCoefs(nzonesD);
  vector<E_Int> posPtrRcvIndices(nzonesD);

  for (E_Int nozd = 0; nozd < nzonesD; nozd++)
  {
    FldArrayI* indicesR = new FldArrayI(nbPtsR); indicesR->setAllValuesAtNull();
    listOfRcvIndices.push_back(indicesR);

    FldArrayI* indicesD = new FldArrayI(nbPtsR*(CLOUDMAX+1)); indicesD->setAllValuesAt(-1);
    listOfDnrIndices.push_back(indicesD);

    FldArrayI* itypes   = new FldArrayI(nbPtsR); itypes->setAllValuesAtNull();
    listOfDnrTypes.push_back(itypes);

    FldArrayF* icoefs   = new FldArrayF(nbPtsR*CLOUDMAX); icoefs->setAllValuesAtNull();
    listOfInterpCoefs.push_back(icoefs);

    posPtrDnrIndices[nozd]=0;
    posPtrCoefs[nozd]=0;
    posPtrRcvIndices[nozd] = 0;
  }

  for (E_Int indR = 0; indR < nbPtsR; indR++)
  {
    pt[0] = xtRcv[indR]; pt[1] = ytRcv[indR]; pt[2] = ztRcv[indR];
    vector<E_Int>& dnrIndices = listOfDnrIndPerTri[indR];
    vector<E_Int>& dnrNoZ = listOfDnrNoZPerTri[indR];
    E_Int sizeOfCloud = dnrIndices.size();
    E_Int isExtrap = 0;
    FldArrayF cfLoc(sizeOfCloud);   
    vector<E_Int> listOfCloudPtsPerVertex(sizeOfCloud);

    if ( sizeOfCloud < sizeMinOfCloud)
    {
      isExtrap = 1;
      indicesExtrap.push_back(indR);
    }
    else
    {      
      FldArrayF cloudCoords(sizeOfCloud,3);
      E_Float* dnrX = cloudCoords.begin(1);
      E_Float* dnrY = cloudCoords.begin(2);
      E_Float* dnrZ = cloudCoords.begin(3);

      for (E_Int noind = 0; noind < sizeOfCloud; noind++)
      {
        E_Int nozDnr = dnrNoZ[noind]; E_Int indDnr = dnrIndices[noind];
        E_Float* xtDnr = vectOfDnrZones[nozDnr]->begin(posxtd[nozDnr]);
        E_Float* ytDnr = vectOfDnrZones[nozDnr]->begin(posytd[nozDnr]);
        E_Float* ztDnr = vectOfDnrZones[nozDnr]->begin(posztd[nozDnr]);

        dnrX[noind] = xtDnr[indDnr]; dnrY[noind] = ytDnr[indDnr]; dnrZ[noind] = ztDnr[indDnr];    
        listOfCloudPtsPerVertex[noind]=noind;
      }
      axisConst[0] = axisConst[1] = axisConst[2] = 0;
      if (dimPb == 2) axisConst[2] = 1;// pas de z
      radius[0] = radius[1] = radius[2] = 1.;
      axis[0] = 1.; axis[1] = 0.; axis[2] = 0.;
      axis[3] = 0.; axis[4] = 1.; axis[5] = 0.;
      axis[6] = 0.; axis[7] = 0.; axis[8] = 1.;
      E_Int ok = K_INTERP::getInterpCoefMLS(order, dimPb, sizeBasis, pt, dnrX, dnrY, dnrZ, 
                                            listOfCloudPtsPerVertex, radius, axis, axisConst, cfLoc.begin());
      listOfCloudPtsPerVertex.clear();

      if ( ok != 1) // extrapolation
      {
        isExtrap = 1;
        indicesExtrap.push_back(indR);
      }
    }//test sizeOfCloud < sizeMin ?
    if (isExtrap==0)
    {
      // Nb de donneurs ds la molecule + indices des donneurs
      vector < vector<E_Int> > vectOfIndDPerDnrZone(nzonesD);

      for (E_Int noind = 0; noind < sizeOfCloud; noind++)
      {
        E_Int nozDnr = dnrNoZ[noind]; 
        vector<E_Int>& vectOfDnrIndices = vectOfIndDPerDnrZone[nozDnr];
        vectOfDnrIndices.push_back(noind);
      }

      for (E_Int nozDnr = 0; nozDnr < nzonesD; nozDnr++)
      {
        vector<E_Int>& vectOfDnrIndices = vectOfIndDPerDnrZone[nozDnr];

        E_Int nbOfDnrIndices=vectOfDnrIndices.size();
        if (nbOfDnrIndices>0)
        {
          E_Int* ptrIndicesR = listOfRcvIndices[nozDnr]->begin();
          E_Int& posIndR = posPtrRcvIndices[nozDnr];
          ptrIndicesR[posIndR] = indR; posIndR++;

          E_Float* ptrCoefs  = listOfInterpCoefs[nozDnr]->begin();
          E_Int& posCf = posPtrCoefs[nozDnr];

          E_Int* ptrIndicesD = listOfDnrIndices[nozDnr]->begin();
          E_Int& posIndD = posPtrDnrIndices[nozDnr];
          ptrIndicesD[posIndD] = vectOfIndDPerDnrZone[nozDnr].size(); posIndD++;

          for (E_Int ii = 0; ii < nbOfDnrIndices; ii++)
          {
            E_Int noind = vectOfDnrIndices[ii];
            ptrIndicesD[posIndD] = dnrIndices[noind]; posIndD++;
            ptrCoefs[posCf] = cfLoc[noind]; posCf++;
          }
        }
      }
    }//fin isExtrap=0
    else
    {
      //recherche du pt le plus proche
      pt[0] = xtRcv[indR]; pt[1] = ytRcv[indR]; pt[2] = ztRcv[indR];
      E_Int noindg = kdtD->getClosest(pt, rad2);

      E_Int nozDorig = E_Int(indicesGlob[noindg]/nptsMaxD);
      E_Int indDorig = (indicesGlob[noindg]-nozDorig*nptsMaxD);
      // if ( indR == 734) 
      // {
      //   printf(" x = %g %g %g : %d \n",pt[0], pt[1], pt[2], isExtrap);
      //   printf(" donneur %g %g %g \n", xdnr[noindg], ydnr[noindg], zdnr[noindg]);
      //   printf(" indD = %d %d \n", indDorig, nozDorig);
      // }
      E_Int* ptrIndicesR = listOfRcvIndices[nozDorig]->begin();
      E_Int& posIndR = posPtrRcvIndices[nozDorig];
      ptrIndicesR[posIndR] = indR; posIndR++;

      E_Float* ptrCoefs  = listOfInterpCoefs[nozDorig]->begin();
      E_Int& posCf = posPtrCoefs[nozDorig];
      
      E_Int* ptrIndicesD = listOfDnrIndices[nozDorig]->begin();
      E_Int& posIndD = posPtrDnrIndices[nozDorig];
      ptrIndicesD[posIndD] = 1; posIndD++;
      ptrIndicesD[posIndD] = indDorig; posIndD++;
      ptrCoefs[posCf] = 1.; posCf++;
    }//fin isExtrap=1
  } // boucle indR
  delete kdtD; delete coordAccD;

  for (E_Int nozd = 0; nozd < nzonesD; nozd++)
  {
    E_Int sizeDnrIndices = posPtrDnrIndices[nozd];
    E_Int sizeCf = posPtrCoefs[nozd];
    E_Int sizeRcvIndices = posPtrRcvIndices[nozd];
    printf(" zone donneuse " SF_D_ " : sizeDnrIndices = " SF_D_ ", sizeCf = " SF_D_ ", sizeRcvIndices = " SF_D_ "\n", nozd, sizeDnrIndices, sizeCf, sizeRcvIndices);
    listOfInterpCoefs[nozd]->resize(sizeCf);
    listOfDnrIndices[nozd]->resize(sizeDnrIndices);
    listOfRcvIndices[nozd]->resize(sizeRcvIndices);
    listOfDnrTypes[nozd]->resize(sizeRcvIndices);
  }   
  RELEASESHAREDU(arrayR, fr, cnr);
  // Nettoyages
  for (E_Int no = 0; no < nzonesD; no++)
  {
    RELEASESHAREDA(resl[no],objs[no], vectOfDnrZones[no],a2[no],a3[no],a4[no]);  
  }
  /*----------------------------------------------------------*/
  /* Ecriture dans des objets Python retournes par la methode */
  /*----------------------------------------------------------*/
  // listes pour stocker les indices de cellules (receveuse et donneuses) par bloc de domaine d'interpolation correspondant
  PyObject* PyListCellIndicesR = PyList_New(0);
  PyObject* PyListCellIndicesD = PyList_New(0);
  // liste pour stocker les coefficients d'interpolation par bloc de domaine d'interpolation correspondant
  PyObject* PyListCoefficients = PyList_New(0);
  // liste pour stocker les types d interpolation par domaine d'interpolation correspondant
  PyObject* PyListInterpTypes = PyList_New(0);

  // ------------------------------
  // Construction des PyArrays
  // ------------------------------
  for (E_Int noz = 0; noz < nzonesD; noz++)
  {   
    //     coefficients d'interpolation
    PyObject* fout = K_NUMPY::buildNumpyArray(*listOfInterpCoefs[noz],1);
    PyList_Append(PyListCoefficients, fout);  Py_DECREF(fout);
    delete listOfInterpCoefs[noz];

    //     donorIndices1D
    PyObject* donorIndOut = K_NUMPY::buildNumpyArray(*listOfDnrIndices[noz],1);
    PyList_Append(PyListCellIndicesD, donorIndOut); Py_DECREF(donorIndOut);
    delete listOfDnrIndices[noz];

     //     receiverIndices1D
    PyObject* rcvIndOut = K_NUMPY::buildNumpyArray(*listOfRcvIndices[noz],1);
    PyList_Append(PyListCellIndicesR, rcvIndOut); Py_DECREF(rcvIndOut);
    delete listOfRcvIndices[noz];

    //     donorType
    PyObject* donorTypeOut = K_NUMPY::buildNumpyArray(*listOfDnrTypes[noz],1);
    PyList_Append(PyListInterpTypes, donorTypeOut); Py_DECREF(donorTypeOut);
    delete listOfDnrTypes[noz];
  } // fin parcours des zones donneuses
  
  PyObject* tpl = Py_BuildValue("[OOOO]", PyListCellIndicesR, PyListCellIndicesD, 
                                PyListInterpTypes, PyListCoefficients);

  Py_DECREF(PyListInterpTypes); 
  Py_DECREF(PyListCoefficients); 
  Py_DECREF(PyListCellIndicesR); 
  Py_DECREF(PyListCellIndicesD); 
  return tpl;
}
