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
// # include <stdio.h>
# include "connector.h"

#define CLOUDMAX 120

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Calcule et stocke les coefficients d'interpolation selon la methode 
   des moindres carres
   CAS SANS DOUBLE WALL 
   IN: receiverArray: points a interpoler definis sous forme de maillage
   IN: donorArrays: maillages donneurs. La localisation du donneur 
        (noeuds/centres/centres etendus) doit �tre effectuee au prealable
   IN: Nature: 0: produit des cellN=0 -> donneur invalide; 
               1: cellN=0 ou 2 -> donneur invalide
   IN: PenalizeBorders: 1: penalite sur le volume des pts ou cellules frontieres
   IN: hook != Py_None: hook sur les kdt associes � donorArrays 
   OUT: [donorBlks,donorInd1D, donorType, coefs, extrapInd1D, orphanInd1D] 
        donorBlks: no du blk donneur, demarre a 0
        donorInd1D: [N1,ind1,...,indN1,N2,indp1,...,indpN2,...] sous forme d un numpy 1D
        donorType: type d'interpolation effectue localement
        coefs: coefficients d'interpolation sous forme d'un numpy 1D
        extrapInd1D: indices des pts extrapoles
        orphanInd1D: indices des pts orphelins */
//=============================================================================
PyObject* K_CONNECTOR::setInterpDataLS(PyObject* self, PyObject* args)
{
  PyObject *receiverArray; 
  PyObject *donorArrays; // domaines d'interpolation
  E_Int Order, Nature, Dim;
  E_Int PenalizeBorders;
  PyObject* hook;

  if (!PYPARSETUPLE_(args, OO_ III_ O_ I_,
                    &receiverArray, &donorArrays, &Order, &Nature, &PenalizeBorders, &hook, &Dim))
  {
      return NULL;
  }
  E_Int order = E_Int(Order);  //Ordre de la formule = ordre maximal dans la base de polynome-1
  if (order < 1 || order > 4)
  {
    printf("Warning: setInterpDataLS: unknown interpolation order.");
    printf(" Set to 3rd order.\n");
    order = 3;
  }
  
  E_Int dimPb = E_Int(Dim); //Dimension du probleme
  E_Int nature = E_Int(Nature);//O: produit des cellN=0 -> donneur invalide; 1 : cellN=0 ou 2 -> donneur invalide
  E_Int penalty = E_Int(PenalizeBorders);//1 : penalite sur le volume des pts ou cellules frontieres

  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  E_Int imr, jmr, kmr;
  FldArrayF* fr; FldArrayI* cnr;
  char* varStringr; char* eltTyper;
  E_Int resr = K_ARRAY::getFromArray(receiverArray, varStringr, fr, 
                                     imr, jmr, kmr, cnr, eltTyper, true); 
  
  // Verif des coordonnees dans la zone a interpoler
  E_Int posxr = K_ARRAY::isCoordinateXPresent(varStringr);
  E_Int posyr = K_ARRAY::isCoordinateYPresent(varStringr);
  E_Int poszr = K_ARRAY::isCoordinateZPresent(varStringr);
  if (posxr == -1 || posyr == -1 || poszr == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setInterpDataLS: 1st arg must contain coordinates.");
    RELEASESHAREDB(resr, receiverArray, fr, cnr); 
    return NULL;
  }
  posxr++; posyr++; poszr++;
 
  /*-------------------------------------------------------*/
  /* Extraction des infos sur les domaines d'interpolation */
  /*-------------------------------------------------------*/
  vector<E_Int> resl;  vector<char*> varString;
  vector<FldArrayF*> fields;
  vector<void*> a2; //ni,nj,nk ou cnt en NS
  vector<void*> a3; //eltType en NS
  vector<void*> a4;
  vector<PyObject*> objs;
  E_Bool skipNoCoord = true;  E_Bool skipStructured = false;
  E_Bool skipUnstructured = false;  E_Bool skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    donorArrays, resl, varString, fields, a2, a3, a4, objs,  
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nDnrZones = objs.size();

  if (isOk == -1)
  {
    RELEASESHAREDB(resr, receiverArray, fr, cnr); 
    for (E_Int no = 0; no < nDnrZones; no++)
      RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
    PyErr_SetString(PyExc_TypeError,
                    "setInterpDataLS: 2nd argument is not valid.");
    return NULL;
  }   
  if (nDnrZones == 0) 
  {
    RELEASESHAREDB(resr, receiverArray, fr, cnr); 
    for (E_Int no = 0; no < nDnrZones; no++)
      RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
    PyErr_SetString(PyExc_TypeError,
                    "setInterpDataLS: no valid donor zone found.");
    return NULL;
  }

  vector<E_Int> posxd; vector<E_Int> posyd; vector<E_Int> poszd; vector<E_Int> poscd;
  for (E_Int no = 0; no < nDnrZones; no++)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString[no]); posx++; posxd.push_back(posx); 
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString[no]); posy++; posyd.push_back(posy); 
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString[no]); posz++; poszd.push_back(posz); 
    E_Int posc = K_ARRAY::isCellNatureField2Present(varString[no]); posc++; poscd.push_back(posc);
  }
  
  // Verifier dans le cas NS que le maillage est NGON
  for (E_Int noz = 0; noz < nDnrZones; noz++)
  {
    if (resl[noz] == 2)
    {
      char* eltType0 = (char*)a3[noz];
      if (K_STRING::cmp(eltType0, "NGON") != 0)
      {
        RELEASESHAREDB(resr, receiverArray, fr, cnr); 
        for (E_Int no = 0; no < nDnrZones; no++)
          RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
        PyErr_SetString(PyExc_TypeError,
                        "setInterpDataLS: unstructured donor zones must be NGON.");
        return NULL;
      }
    }
  }
  
  // Construction de l'arbre de recherche par zone donneuse des points de cellN 1
  vector<K_SEARCH::KdTree<FldArrayF>*> vectOfKdTrees;
  vector<ArrayAccessor<FldArrayF>*> vectOfCoordAcc;
  vector<FldArrayI*> vectOfCorres;
  vector<FldArrayF*> vectOfCellN;
  if (hook == Py_None)
  {
    for (E_Int no = 0; no < nDnrZones; no++)
    {
      E_Float* x = fields[no]->begin(posxd[no]);
      E_Float* y = fields[no]->begin(posyd[no]);
      E_Float* z = fields[no]->begin(poszd[no]);
      E_Float* cellN = fields[no]->begin(poscd[no]);
      E_Int nptsmax = fields[no]->getSize();
      
      FldArrayF* cellN1pts = new FldArrayF(nptsmax, 3);
      E_Float* x2 = cellN1pts->begin(1);
      E_Float* y2 = cellN1pts->begin(2);
      E_Float* z2 = cellN1pts->begin(3);
      FldArrayI* corresF = new FldArrayI(nptsmax);
      E_Int* corres = corresF->begin();

      E_Int pos = 0;
      if (nature == 0)
      {
        for (E_Int ind = 0; ind < nptsmax; ind++)
        {
          if (cellN[ind] != 0.) // On ne garde que les points de cellN=1 ou 2
          {
            x2[pos] = x[ind]; y2[pos] = y[ind]; z2[pos] = z[ind]; 
            corres[pos] = ind;
            pos++;
          }
        }
      }
      else // nature = 1
      {
        for (E_Int ind = 0; ind < nptsmax; ind++)
        {
          if (cellN[ind] == 1.) // On ne garde que les points de cellN=1
          {
            x2[pos] = x[ind]; y2[pos] = y[ind]; z2[pos] = z[ind]; 
            corres[pos] = ind;
            pos++;
          }
        }
      }
      cellN1pts->reAllocMat(pos, 3);
      corresF->resize(pos);
      //corres = corresF-> begin();

      ArrayAccessor<FldArrayF>* coordAcc = new ArrayAccessor<FldArrayF>(*cellN1pts, 1, 2, 3);   
      K_SEARCH::KdTree<FldArrayF>* kdt = new K_SEARCH::KdTree<FldArrayF>(*coordAcc);
      vectOfCoordAcc.push_back(coordAcc);
      vectOfKdTrees.push_back(kdt);
      vectOfCorres.push_back(corresF);
      vectOfCellN.push_back(cellN1pts);
    }
  }
  else //if (hook != Py_None) // hook fourni
  {
    RELEASESHAREDB(resr, receiverArray, fr, cnr); 
    for (E_Int no = 0; no < nDnrZones; no++)
      RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpDataLS: a hook cannot be used for MLS interpolation.");
    return NULL;
  }

  // Cas NS: calcul des connectivites
  vector<vector<vector<E_Int> > > vectOfcVN(nDnrZones); //connectivites Vertex/Vertex voisins  
  vector<vector<vector<E_Int> > > vectOfcVF; //connectivites Vertex/Faces
  vector<vector<vector<E_Int> > > vectOfcEV; //connectivites Element/noeuds
  vector<FldArrayI> vectOfcFE;     //connectivites Faces/Elts
  vector<FldArrayI> vectOfPosFace; //tableaux de position des faces dans la connectivite
  vector<FldArrayI> vectOfPosElt;  //tableaux de position des elements dans la connectivite
  vector<FldArrayI> vectOfDimElt;  //tableaux de la dimension des elements

  for (E_Int noz = 0; noz < nDnrZones; noz++)
  {
    if (resl[noz] == 2)
    {
      FldArrayI& cNG = *(FldArrayI*)a2[noz];
      E_Int* cnp = cNG.begin();
      E_Int npts = fields[noz]->getSize();  // nombre total de noeuds
      E_Int nfaces = cnp[0]; // nombre total de faces
      E_Int ne = cnp[cnp[1]+2];  // nombre total d elements

      // connectivite Vertex/Vertex voisins  
      vector< vector<E_Int> > cVN(npts);
      K_CONNECT::connectNG2VNbrs(cNG, cVN);
      vectOfcVN[noz] = cVN;
      // connectivite Vertex/Faces
      vector< vector<E_Int> > cVF(npts);
      K_CONNECT::connectNG2VF(cNG, cVF);
      vectOfcVF.push_back(cVF);
      // Connectivite Element/Noeuds
      vector< vector<E_Int> > cEV(ne);
      K_CONNECT::connectNG2EV(cNG, cEV);
      vectOfcEV.push_back(cEV);
      // connectivite Faces/Elts
      FldArrayI cFE;
      K_CONNECT::connectNG2FE(cNG, cFE);
      vectOfcFE.push_back(cFE);
      // tableau de position des faces dans la connectivite 
      FldArrayI posFace(nfaces); 
      K_CONNECT::getPosFaces(cNG, posFace);
      vectOfPosFace.push_back(posFace);
      // tableau de position des elements dans la connectivite
      FldArrayI posElt(npts); 
      K_CONNECT::getPosElts(cNG, posElt);
      vectOfPosElt.push_back(posElt);
      // tableau de la dimension des elements
      FldArrayI dimElt(npts); 
      K_CONNECT::getDimElts(cNG, dimElt);
      vectOfDimElt.push_back(dimElt);
    }
    else //si la zone est structuree
    {
      vector< vector<E_Int> > nullv(0);
      FldArrayI nullf(0);
      vectOfcVN[noz] = nullv;
      vectOfcVF.push_back(nullv);
      vectOfcEV.push_back(nullv);
      vectOfcFE.push_back(nullf);
      vectOfPosFace.push_back(nullf);
      vectOfPosElt.push_back(nullf);
      vectOfDimElt.push_back(nullf);
    }
  }

  /*-------------------------------------------------------*/
  /* Calcul des coefficients d'interpolation               */
  /*-------------------------------------------------------*/
  E_Float* xr = fr->begin(posxr);
  E_Float* yr = fr->begin(posyr);
  E_Float* zr = fr->begin(poszr);
  E_Int nbI = fr->getSize();
  
  // Donnees d'interpolation
  vector<FldArrayF*> listOfInterpCoefs;
  vector<FldArrayI*> listOfDonorInd1D;
  vector<FldArrayI*> listOfRcvInd1D;
  vector<FldArrayI*> listOfDonorTypes;
  vector<FldArrayI*> listOfExtrapInd1D;

  // Initialisation des tableaux et dimensionnement a priori
  // On met pour l'instant NCloudPtsMax � CLOUDMAX: la molecule d'interpolation contient a priori CLOUDMAX pts
  E_Int nCloudPtsMax = CLOUDMAX;
  for (E_Int noz = 0; noz < nDnrZones; noz++)
  {
    FldArrayF* donorInterpCoefs = new FldArrayF(nbI*nCloudPtsMax); donorInterpCoefs->setAllValuesAtNull();
    listOfInterpCoefs.push_back(donorInterpCoefs);
    FldArrayI* donorInd1D = new FldArrayI(nbI*(nCloudPtsMax+1)); donorInd1D->setAllValuesAt(-1);
    listOfDonorInd1D.push_back(donorInd1D);
    FldArrayI* rcvInd1D = new FldArrayI(nbI); rcvInd1D->setAllValuesAt(-1);
    listOfRcvInd1D.push_back(rcvInd1D);
    FldArrayI* donorType = new FldArrayI(nbI); donorType->setAllValuesAtNull();
    listOfDonorTypes.push_back(donorType);
    FldArrayI* extrapPoints = new FldArrayI(nbI); extrapPoints->setAllValuesAt(-1);
    listOfExtrapInd1D.push_back(extrapPoints);
  }

  FldArrayI* orphanPts = new FldArrayI(nbI); orphanPts->setAllValuesAt(-1);
  E_Int* orphanPts1D = orphanPts->begin();
  FldArrayI usedDomp(nDnrZones); usedDomp.setAllValuesAtNull();//nb de pts interpoles par domaine donneur
  FldArrayI sizeIndDnrs(nDnrZones); sizeIndDnrs.setAllValuesAtNull();// taille de donorInd1D par domaine donneur 
  FldArrayI sizeCoefs(nDnrZones); sizeCoefs.setAllValuesAtNull();//taille du tableau de coeffs par domaine donneur
  FldArrayI nbExtraps(nDnrZones); nbExtraps.setAllValuesAtNull();

  E_Int nExtrapPts=0; // pour l'instant pas d'extrapolation: comment la qualifier en moindres carres ? 
  E_Int noOrphanPt=0;

  E_Float pt[3]; // pt a interpoler
  E_Int axisConst[3]; // =1 si l'axe est une direction constante, 0 sinon
  // Taille du probleme et direction constante
  // si axisConst=1, direction par prise en compte dans la base
  axisConst[0] = axisConst[1] = axisConst[2] = 0;
  
  if (dimPb == 3)        // pas de direction constante et probleme 3D
  { axisConst[0] = 0; axisConst[1] = 0; axisConst[2] = 0;}
  else if (dimPb == 2)   // pas de direction constante mais probleme 2D
  {
    axisConst[2] = 1; // pas de z
  }
  else if (dimPb == 1)
  {
    axisConst[1] = 1; axisConst[2] = 1; // pas de y/z
  }
        
  // Taille de la base de polynome
  const E_Int sizeBasis = K_FUNC::fact(dimPb+order-1)*1./(K_FUNC::fact(dimPb)*K_FUNC::fact(order-1));

  E_Int isExtrapolated = 0;
  E_Int ok = 0;
  for (E_Int ind = 0; ind < nbI; ind++)
  {
    pt[0] = xr[ind]; pt[1] = yr[ind]; pt[2] = zr[ind];
    
    isExtrapolated = 0;

    // Recherche du meilleur donneur
    E_Int depth = 2; // profondeur du stencil de recherche necessaire pour penalizeBorders
    E_Int no;         //numero de la meilleure zone donneuse
    E_Int indBestDnr; //indice du meilleur point donneur
    E_Float dBestDnr; //distance au meilleur point donneur
    K_INTERP::getBestDonor(dimPb, pt, nDnrZones, fields, vectOfKdTrees, resl,
                           a2, a3, a4, posxd, posyd, poszd, poscd, vectOfCorres, 
                           vectOfcVF, vectOfcEV, vectOfcVN, 
                           vectOfcFE, vectOfPosElt, vectOfPosFace, 
                           vectOfDimElt, depth, penalty,
                           no, indBestDnr, dBestDnr);

    if (no == -1) // point orphelin
    {
      printf("Warning: setInterpDataLS: no valid donor found for point %5.12f %5.12f %5.12f \n",
             pt[0], pt[1], pt[2]);
      printf("Point is marked as orphan\n");
      orphanPts1D[noOrphanPt] = ind; noOrphanPt++;
    }
    else 
    {
      // Recuperation des infos du meilleur donneur
      E_Float* xtDnr = fields[no]->begin(posxd[no]); //x
      E_Float* ytDnr = fields[no]->begin(posyd[no]);
      E_Float* ztDnr = fields[no]->begin(poszd[no]);
      E_Float* cellNtDnr = fields[no]->begin(poscd[no]); // cellN 
   
      // Dimension du maillage donneur
      E_Int ni = 0, nj = 0, nk = 0;
      if (resl[no] == 1)
      {
        ni=*(E_Int*)a2[no];
        nj=*(E_Int*)a3[no];
        nk=*(E_Int*)a4[no];   
      }
    
      // Recherche du nuage de points contenant le pt (x,y,z) 
      // methode de recherche du nuage de points: stencil elliptique
      E_Float radius[3]; // longueurs des axes de l'ellipse
      radius[0] = radius[1] = radius[2] = 1.;
      E_Float axis[9];  // Axes de l'ellipse
      vector<E_Int> indicesIn; //indices des pts situes dans le stencil autour de pt
      if (K_FUNC::fEqual(dBestDnr, 0., 1.e-13) == true) //maillage coincident: pas de MLS 
        indicesIn.push_back(indBestDnr);
      else
      {
        ok = K_INTERP::buildStencil(dimPb, order, nature, sizeBasis,
                                    indBestDnr, resl[no], ni, nj, nk, vectOfcVN[no],
                                    xtDnr, ytDnr, ztDnr, cellNtDnr,
                                    pt, radius, axis, axisConst, indicesIn);
        if ( ok != 1 ) 
        {
          printf("Warning: setInterpDataLS: no valid stencil found for point %5.12f %5.12f %5.12f \n",
                 pt[0], pt[1], pt[2]);
          printf("Point is marked as orphan\n");
          orphanPts1D[noOrphanPt] = ind; noOrphanPt++;
        }
      }
      E_Int nDnrPts = indicesIn.size();

      if (nDnrPts == 1 || (nDnrPts >= sizeBasis && nDnrPts < nCloudPtsMax)) 
      {                                                                
        E_Int noDonorBlk = no;
        E_Int& noInterpPtForBlk = usedDomp[noDonorBlk];// on met une reference car on va modifier cette valeur
        E_Int& sizeIndDnr = sizeIndDnrs[noDonorBlk];//idem
        E_Int& sizeCf = sizeCoefs[noDonorBlk];
        E_Int& noExtrapPtForBlk = nbExtraps[noDonorBlk];
      
        E_Float* donorCf = listOfInterpCoefs[noDonorBlk]->begin();
        E_Int* donorType = listOfDonorTypes[noDonorBlk]->begin();
        E_Int* donorInd1D = listOfDonorInd1D[noDonorBlk]->begin();
        E_Int* rcvInd1D = listOfRcvInd1D[noDonorBlk]->begin();
        E_Int* extrapInd1D = listOfExtrapInd1D[noDonorBlk]->begin();

        // type d'interpolation a appliquer: 0: LS
        donorType[noInterpPtForBlk] = 0;  

        // Nb de donneurs ds la molecule + indices des donneurs
        donorInd1D[sizeIndDnr] = nDnrPts; sizeIndDnr++;

        for (E_Int ii = 0; ii < nDnrPts; ii++)
        {donorInd1D[sizeIndDnr] = indicesIn[ii]; sizeIndDnr++;} 
      
        // Indice du receveur stocke pour le bloc noDonorBlk
        rcvInd1D[noInterpPtForBlk] = ind;

        // Calcul des coefficients d'interpolation par LS
        E_Float cfLoc[nDnrPts];
        for (E_Int ii = 0; ii < nDnrPts; ii++) { cfLoc[ii] = 0.; }
        
      
        if (nDnrPts > 1)//MLS
        {
          ok = K_INTERP::getInterpCoefMLS(order, dimPb, sizeBasis, pt, xtDnr, ytDnr, ztDnr, 
                                          indicesIn, radius, axis, axisConst, cfLoc);
          if ( ok != 1 ) 
          {
            printf("Warning: setInterpDataLS: MLS coefficients cannot be computed for point %5.12f %5.12f %5.12f \n", pt[0], pt[1], pt[2]);
            printf("Point is marked as orphan\n");
            orphanPts1D[noOrphanPt] = ind; noOrphanPt++;
          }
        }
        else if (nDnrPts == 1) cfLoc[0] = 1.;

        // Insertion des coefficients dans le tableau associe au donneur
        for (E_Int ii = 0; ii < nDnrPts; ii++)
        { donorCf[sizeCf] = cfLoc[ii]; sizeCf++;} 

        // Incremente le nb de pts interpoles pour le bloc noDonorBlk
        if (isExtrapolated == 1) 
        { 
          extrapInd1D[noExtrapPtForBlk] = ind; 
          noExtrapPtForBlk++;// = nbExtraps[noDonorBlk]++;
        }
        nExtrapPts++;
        noInterpPtForBlk++;
      }
      else if (nDnrPts < sizeBasis)
      {
        printf("Warning: setInterpDataLS: not enough points in the donor stencil for point %5.12f %5.12f %5.12f \n", pt[0], pt[1], pt[2]);
        printf("Point is marked as orphan\n");
        orphanPts1D[noOrphanPt] = ind; noOrphanPt++;     
      }
      else if (nDnrPts > nCloudPtsMax) 
      {
        printf("Warning: setInterpDataLS: too many points in the donor stencil for point %5.12f %5.12f %5.12f \n", pt[0], pt[1], pt[2]);
        printf("Point is marked as orphan\n");
        orphanPts1D[noOrphanPt] = ind; noOrphanPt++; 
      }
    }
  }
  
  // Redimensionnement des tableaux 
  for (E_Int noz = 0; noz < nDnrZones; noz++)
  {
    E_Int sizecf = sizeCoefs[noz];
    listOfInterpCoefs[noz]->resize(sizecf);
    E_Int sizeInd = sizeIndDnrs[noz];
    listOfDonorInd1D[noz]->resize(sizeInd);
    E_Int nbInterp = usedDomp[noz];
    listOfRcvInd1D[noz]->resize(nbInterp);
    listOfDonorTypes[noz]->resize(nbInterp);
    E_Int nbExtrap = nbExtraps[noz];   
    listOfExtrapInd1D[noz]->resize(nbExtrap);
  }
  // Nettoyages
  for (E_Int no = 0; no < nDnrZones; no++)
  {
    if (hook == Py_None) 
    {
      delete vectOfKdTrees[no];
      delete vectOfCoordAcc[no];
      delete vectOfCorres[no];
      delete vectOfCellN[no];
    }
    // delete vectOfcVN[no];
    
    RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);  
  }
  RELEASESHAREDB(resr, receiverArray, fr, cnr); 
  // Fin nettoyages
  
  // Marquage des cellules orphelines
  orphanPts->resize(noOrphanPt);

  /*----------------------------------------------------------*/
  /* Ecriture dans des objets Python retournes par la methode */
  /*----------------------------------------------------------*/
  // listes pour stocker les indices de cellules (receveuse et donneuses) par bloc de domaine d'interpolation correspondant
  PyObject * PyListCellIndicesR = PyList_New(0);
  PyObject * PyListCellIndicesD = PyList_New(0);
  // liste pour stocker les coefficients d'interpolation par bloc de domaine d'interpolation correspondant
  PyObject * PyListCoefficients = PyList_New(0);
  // liste pour stocker les types d interpolation par domaine d'interpolation correspondant
  PyObject * PyListInterpTypes = PyList_New(0);
 // listes pour marquer les indices des cellules extrapoles
  PyObject * PyListCellIndicesE = PyList_New(0);
  // ------------------------------
  // Construction des PyArrays
  // ------------------------------
  for (E_Int noz = 0; noz < nDnrZones; noz++)
  {   
    //     coefficients d'interpolation
    PyObject* fout = K_NUMPY::buildNumpyArray(*listOfInterpCoefs[noz],1);
    PyList_Append(PyListCoefficients, fout);  Py_DECREF(fout);
    delete listOfInterpCoefs[noz];

    //     donorIndices1D
    PyObject* donorIndOut = K_NUMPY::buildNumpyArray(*listOfDonorInd1D[noz],1);
    PyList_Append(PyListCellIndicesD, donorIndOut); Py_DECREF(donorIndOut);
    delete listOfDonorInd1D[noz];

     //     receiverIndices1D
    PyObject* rcvIndOut = K_NUMPY::buildNumpyArray(*listOfRcvInd1D[noz],1);
    PyList_Append(PyListCellIndicesR, rcvIndOut); Py_DECREF(rcvIndOut);
    delete listOfRcvInd1D[noz];

    //     donorType
    PyObject* donorTypeOut = K_NUMPY::buildNumpyArray(*listOfDonorTypes[noz],1);
    PyList_Append(PyListInterpTypes, donorTypeOut); Py_DECREF(donorTypeOut);
    delete listOfDonorTypes[noz];
    
    //   indices des points extrapoles 
    PyObject* cellIndE = K_NUMPY::buildNumpyArray(*listOfExtrapInd1D[noz],1);
    PyList_Append(PyListCellIndicesE, cellIndE); Py_DECREF(cellIndE);
    delete listOfExtrapInd1D[noz];
  } // fin parcours des zones donneuses
  
  //   indices des points orphelins
  PyObject* cellIndO = K_NUMPY::buildNumpyArray(*orphanPts,1);
  delete orphanPts;
  // listes pour marquer les EXdir (directions pour les pts EX)
  PyObject * PyListEXDir = PyList_New(0); // n est pas utilisee ici
  
  PyObject* tpl = Py_BuildValue("[OOOOOOO]", PyListCellIndicesR, PyListCellIndicesD, 
                                PyListInterpTypes, PyListCoefficients,
                                PyListCellIndicesE, cellIndO, PyListEXDir);
  Py_DECREF(PyListCellIndicesR); 
  Py_DECREF(PyListCellIndicesD);
  Py_DECREF(PyListInterpTypes); 
  Py_DECREF(PyListCoefficients); 
  Py_DECREF(cellIndO); 
  Py_DECREF(PyListCellIndicesE);
  Py_DECREF(PyListEXDir); 
  return tpl;
}
