/*    
    Copyright 2013-2026 ONERA.

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
# include <stdio.h>
# include "connector.h"
# include "KInterp/BlkInterp.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Calcul et stocke les coefficients d'interpolation pour les centres des 
   cellules
   IN: (Nir,Njr): dimensions (i,j) des blocs interpoles (receveurs)
   IN: coordArrays: pts a interpoler
   IN: interpArrays: domaines d interpolation en noeuds
   IN: interpCellN: cellN en centres pour les domaines d interpolation
   IN: isEX: indique si on est en depth=2 (isEX = 0) ou en depth=1 (isEX=1)
   OUT: indices des points interpoles et des cellules donneuses et coefficients d'interpolation 
         ( + liste des directions pour les points EX quand isEX=1)
*/
//=============================================================================
PyObject* K_CONNECTOR::setInterpolations(PyObject* self, PyObject* args)
{
  PyObject *coordArrays; // pts a interpoler: vecteur par domaine d'interpolation pour le double wall
  PyObject *interpArrays, *interpCellN;// domaines d'interpolations
  E_Int Nir, Njr;
  E_Int isEX = 0;
  E_Int Zid;
  E_Float geomCutOff = 1.e-8;
  E_Float cfMax = 30.;
  //E_Float extrapTol;//tolerance sur la  somme des coefs d'extrapolation
  if (!PYPARSETUPLE_(args, II_ OOO_ II_ R_,
                     &Nir, &Njr, &coordArrays, &interpArrays, &interpCellN, &isEX, &Zid, &cfMax))
  {
      return NULL;
  }
  // numero Id de la zone (commence a zero)
  E_Int zid = Zid;

  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  // seulement arrays non structures NODE avec coordonnees ici 
  vector<E_Int> res0;
  vector<char*> structVarString0; vector<char*> unstrVarString0;
  vector<FldArrayF*> structF0; vector<FldArrayF*> unstrF0;
  vector<E_Int> nit0; vector<E_Int> njt0; vector<E_Int> nkt0;
  vector<FldArrayI*> cnt0;
  vector<char*> eltTypet0;
  vector<PyObject*> objs0, obju0;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = true;
  E_Bool skipUnstructured = false;
  E_Bool skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    coordArrays, res0, structVarString0, unstrVarString0,
    structF0, unstrF0, nit0, njt0, nkt0, cnt0, eltTypet0, objs0, obju0, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);

  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setInterpolations: 1st argument is not valid.");
    return NULL;
  }   
  E_Int nzones = unstrF0.size();
  if (nzones == 0) 
  {
    for (E_Int is = 0; is < nzones; is++)
      RELEASESHAREDU(obju0[is], unstrF0[is], cnt0[is]);
    PyErr_SetString(PyExc_TypeError,
                    "setInterpolations: 1st argument is empty.");
    return NULL;
  }
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    if (strcmp(eltTypet0[noz],"NODE") != 0 )
    {
      for (E_Int is = 0; is < nzones; is++)
        RELEASESHAREDU(obju0[is], unstrF0[is], cnt0[is]);
      PyErr_SetString(PyExc_TypeError,
                    "setInterpolations: 1st argument must be unstructured.");
      return NULL;
    }
  }
  
  E_Int posx = K_ARRAY::isCoordinateXPresent(unstrVarString0[0]);
  E_Int posy = K_ARRAY::isCoordinateYPresent(unstrVarString0[0]);
  E_Int posz = K_ARRAY::isCoordinateZPresent(unstrVarString0[0]);
  E_Int posindrcv = K_ARRAY::isNamePresent("indcell",unstrVarString0[0]);
  if (isEX) posindrcv = K_ARRAY::isNamePresent("indcell1",unstrVarString0[0]);
  E_Int posdir = K_ARRAY::isNamePresent("EXdir",unstrVarString0[0]);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    for (E_Int is = 0; is < nzones; is++)
      RELEASESHAREDU(obju0[is], unstrF0[is], cnt0[is]);
    PyErr_SetString(PyExc_TypeError,
                    "setInterpolations: 1st array must contain (x,y,z).");
    return NULL;
  }
  posx++; posy++; posz++; posindrcv++; posdir++;

  if ((isEX)&&(posdir == -1))
  {
    for (E_Int is = 0; is < nzones; is++)
      RELEASESHAREDU(obju0[is], unstrF0[is], cnt0[is]);
    PyErr_SetString(PyExc_TypeError,
                    "setInterpolations: 1st array must contain indirection for EX points.");
    return NULL;
  }
  // dimensions de la grille en centre de cellules du domaine interpole
  //E_Int nir = Nir; E_Int njr = Njr;

  vector<FldArrayF*> vectOfInterpPts;
  vector<FldArrayI*> vectOfIndRcv;
  vector<FldArrayF*> vectOfEXdir;
  for (E_Int v = 0; v < nzones; v++)
  {
    FldArrayF& field = *unstrF0[v];
    E_Int npts = field.getSize();
    FldArrayF* interpPts = new FldArrayF(npts,3);
    FldArrayI* indRcv    = new FldArrayI(npts);
    interpPts->setOneField(field, posx, 1);
    interpPts->setOneField(field, posy, 2);
    interpPts->setOneField(field, posz, 3);
    vectOfInterpPts.push_back(interpPts);
    E_Float* indRcvF = field.begin(posindrcv);
    E_Int* indRcvT = indRcv->begin();
    for (E_Int i = 0; i < npts; i++) indRcvT[i] = E_Int(indRcvF[i]);
    vectOfIndRcv.push_back(indRcv);
    if (isEX)
    {
      FldArrayF* dirEX = new FldArrayF(npts);
      dirEX->setOneField(field, posdir, 1);
      vectOfEXdir.push_back(dirEX);
    }
  }
  for (E_Int is = 0; is < nzones; is++)
    RELEASESHAREDU(obju0[is], unstrF0[is], cnt0[is]); 

  /*-------------------------------------------------------*/
  /* Extraction des infos sur les domaines d'interpolation */
  /*-------------------------------------------------------*/
  // seulement arrays structures avec coordonnees ici 
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt;  vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  skipNoCoord = true;
  skipStructured = false;
  skipUnstructured = true;
  skipDiffVars = true;
  isOk = K_ARRAY::getFromArrays(
    interpArrays, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltTypet, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int structFsize = structF.size();
  if ( isOk == -1 )
  {
    for (E_Int is = 0; is < structFsize; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    PyErr_SetString(PyExc_TypeError,
                    "setInterpolations: 2nd argument is not valid.");
    K_ARRAY::cleanStructFields(vectOfInterpPts); 
    for (unsigned int noz = 0; noz < vectOfIndRcv.size(); noz++)
      delete vectOfIndRcv[noz];
    if (isEX) K_ARRAY::cleanStructFields(vectOfEXdir);
    return NULL;
  }   
  if (structFsize == 0) 
  {
    for (E_Int is = 0; is < structFsize; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    PyErr_SetString(PyExc_TypeError,
                    "setInterpolations: no interpolation domain has been defined.");
    K_ARRAY::cleanStructFields(vectOfInterpPts); 
    for (unsigned int noz = 0; noz < vectOfIndRcv.size(); noz++)
      delete vectOfIndRcv[noz];

    if (isEX) K_ARRAY::cleanStructFields(vectOfEXdir);
    return NULL;
  }
  nzones = structFsize;

  posx = K_ARRAY::isCoordinateXPresent(structVarString[0]); posx++;
  posy = K_ARRAY::isCoordinateYPresent(structVarString[0]); posy++;
  posz = K_ARRAY::isCoordinateZPresent(structVarString[0]); posz++;

  // Extract infos from celln arrays : seulement structures
  vector<E_Int> resc;
  vector<char*> structVarStringc;
  vector<char*> unstrVarStringc;
  vector<FldArrayF*> structFc;
  vector<FldArrayF*> unstrFc;
  vector<E_Int> nitc; 
  vector<E_Int> njtc; 
  vector<E_Int> nktc;
  vector<FldArrayI*> cntc;
  vector<char*> eltTypec;
  vector<PyObject*> objsc, objuc;
  E_Bool skipNoCoordc = false;
  E_Bool skipStructuredc = false;
  E_Bool skipUnstructuredc = true;
  E_Bool skipDiffVarsc = true;
  isOk = K_ARRAY::getFromArrays(
    interpCellN, resc, structVarStringc, unstrVarStringc,
    structFc, unstrFc, nitc, njtc, nktc, cntc, eltTypec, objsc, objuc, 
    skipDiffVarsc, skipNoCoordc, skipStructuredc, skipUnstructuredc, true);
  E_Int nsc = structFc.size();
  if ( isOk == -1 )
  {
    for (E_Int is = 0; is < structFsize; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    for (E_Int is = 0; is < nsc; is++)
      RELEASESHAREDS(objsc[is], structFc[is]);
    PyErr_SetString(PyExc_TypeError,
                    "setInterpolations: 3rd argument is not valid."); 
    K_ARRAY::cleanStructFields(vectOfInterpPts); 
    for (unsigned int noz = 0; noz < vectOfIndRcv.size(); noz++)
      delete vectOfIndRcv[noz];
    if (isEX) K_ARRAY::cleanStructFields(vectOfEXdir);
    return NULL;
  }
  E_Int structFcsize = structFc.size();
  if (structFcsize != nzones)
  {
    for (E_Int is = 0; is < structFsize; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    for (E_Int is = 0; is < nsc; is++)
      RELEASESHAREDS(objsc[is], structFc[is]);
    PyErr_SetString(PyExc_TypeError,
                    "setInterpolations: 2nd and 3rd arguments must be of same size."); 
    K_ARRAY::cleanStructFields(vectOfInterpPts); 
    for (unsigned int noz = 0; noz < vectOfIndRcv.size(); noz++)
      delete vectOfIndRcv[noz];
    if (isEX) K_ARRAY::cleanStructFields(vectOfEXdir);
    return NULL;
  }
  
  // recherche de la variable celln
  E_Int posc = K_ARRAY::isCellNatureField2Present(structVarStringc[0]);
  if ( posc == -1 )
  {
    for (E_Int is = 0; is < structFsize; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    for (E_Int is = 0; is < nsc; is++)
      RELEASESHAREDS(objsc[is], structFc[is]);
    PyErr_SetString(PyExc_TypeError,
                    "setInterpolations: celln variable not found.");
    K_ARRAY::cleanStructFields(vectOfInterpPts);
    for (unsigned int noz = 0; noz < vectOfIndRcv.size(); noz++)
      delete vectOfIndRcv[noz];
    if (isEX) K_ARRAY::cleanStructFields(vectOfEXdir);
    return NULL;
  }
  posc++;
  
  // recuperation des tableaux coords et celln : 
  // celln vaut 0 masque, 1 normal, 2 : interpole
  vector<FldArrayF*> coords;
  vector<FldArrayI*> cellns;
  for (E_Int zone = 0; zone < nzones; zone++)
  {
    FldArrayF& field0 = *structF[zone];
    FldArrayF& fieldc0 = *structFc[zone];
    E_Int npts = field0.getSize();
    E_Int ncells0 = fieldc0.getSize();
    E_Float* cellnp = fieldc0.begin(posc);

    // cellN doit etre un tableau d entier pour la suite
    FldArrayI* celln = new FldArrayI(ncells0);
    E_Int* icellnp = celln->begin();
    for (E_Int i = 0; i < ncells0; i++)
      icellnp[i]=(E_Int)cellnp[i];

    FldArrayF* coord0 = new FldArrayF(npts,3);
    coord0->setOneField(field0, posx, 1);
    coord0->setOneField(field0, posy, 2);
    coord0->setOneField(field0, posz, 3);
    coords.push_back(coord0); cellns.push_back(celln);
  }
  for (E_Int is = 0; is < structFsize; is++)
    RELEASESHAREDS(objst[is], structF[is]);
  for (E_Int is = 0; is < nsc; is++)
    RELEASESHAREDS(objsc[is], structFc[is]);

  /*-------------------------------------------------------*/
  /* Calcul des coefficients d'interpolation               */
  /*-------------------------------------------------------*/
  // InterpDatas structurees 
  // creation des maillages en centres etendus et des interpDatas
  vector<K_KINTERP::KMesh*> listOfKMeshesN;//en centres etendus
  vector<K_KINTERP::KMesh*> listOfExtCentersKMeshes;//en centres etendus
  vector<K_KINTERP::BlkInterpData*> listOfInterpDatas;//en centres etendus
  nzones = coords.size();
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    K_KINTERP::KMesh* kmesh = new K_KINTERP::KMesh(nit[noz], njt[noz], nkt[noz], *coords[noz]);
    listOfKMeshesN.push_back(kmesh);
    K_KINTERP::KMesh* extKMesh = new K_KINTERP::KMesh();
    extKMesh->createExtendedCenterMesh(*kmesh);
    listOfExtCentersKMeshes.push_back(extKMesh);
    K_KINTERP::BlkInterpData* interpData0 = 
      new K_KINTERP::BlkInterpAdt(*extKMesh);  
    listOfInterpDatas.push_back(interpData0);
  }
  InterpData* interpData = new struct InterpData();
  K_KINTERP::BlkInterpData::InterpolationType interpType0 = K_KINTERP::BlkInterpData::O2CF;
  E_Int ncoefs = 8;
  K_KINTERP::BlkInterpData::InterpMeshType interpMeshType = 
    K_KINTERP::BlkInterpData::EXT_CENTERS;

  FldArrayI donorCells(vectOfInterpPts[0]->getSize());

  // working vectors for extrapolated and orphan points 
  vector< vector<E_Int> > vectOfExtrapPts(nzones); 
  vector<E_Int> orphanVect;
  if (isEX == 0) // traitement des centres
    compAndStoreInterpCoefs(geomCutOff, interpType0, interpMeshType,
			    vectOfInterpPts, interpData, donorCells, vectOfExtrapPts, orphanVect,
			    nit, njt, nkt, coords, listOfInterpDatas, cellns, zid, cfMax);
  else //pts EX
    compAndStoreEXInterpCoefs(geomCutOff, interpType0, interpMeshType,
			      vectOfInterpPts, interpData, donorCells, vectOfExtrapPts, orphanVect,
			      nit, njt, nkt, coords, listOfInterpDatas, cellns, zid, cfMax); 

  /*----------------------------------------------------------*/
  /* Ecriture dans des objets Python retournes par la methode */
  /*----------------------------------------------------------*/
  // listes pour stocker les indices de cellules (receveuse et donneuses) par bloc de domaine d'interpolation 
  // correspondant
  PyObject * PyListCellIndicesR = PyList_New(0);
  PyObject * PyListCellIndicesD = PyList_New(0);
  PyObject * PyListCellIndicesDExtC = PyList_New(0); // indices donneurs en centres etendus : pour ecriture ds fichiers elsA

  // listes pour marquer les indices des cellules extrapoles
  PyObject * PyListCellIndicesE = PyList_New(0);
  // listes pour marquer les indices des cellules orphelines
  PyObject * PyListCellIndicesO = PyList_New(0);
  // liste pour stocker le tableau d'indirection pour les points EX
  PyObject * PyListEXDir = PyList_New(0);
  // liste pour stocker les coefficients d'interpolation par bloc de domaine d'interpolation correspondant
  PyObject * PyListCoefficients = PyList_New(0);
  // liste pour stocker les volumes des cellules d'interpolation par bloc de domaine d'interpolation correspondant
  PyObject * PyListVolumes = PyList_New(0);
  // liste pour stocker les types d interpolation par domaine d'interpolation correspondant
  PyObject * PyListInterpTypes = PyList_New(0);

  E_Int nb = 0; E_Int li;
  E_Int c = 0;
  // Parcours des zones d'interpolation
  for (E_Int noz=0; noz<nzones; noz++)
  {
    FldArrayI& indRcv =  *vectOfIndRcv[noz];
    vector<E_Int>& extrapPts = vectOfExtrapPts[noz];
    E_Int nextrapPts = extrapPts.size();
    FldArrayI* extrapPtsP = new FldArrayI(nextrapPts);
    c = 0;
    for (E_Int ie=0; ie < nextrapPts; ie++)
    {
      E_Int inde = extrapPts[ie];
      E_Int indExtrap = indRcv[inde];
      (*extrapPtsP)[ie] = indExtrap; 
      c++;
    }
    
    // Dimensions de la grille en centre de cellules da zone d interpolation
    //E_Int nid = nitc[noz]; E_Int njd = njtc[noz];
    // Tableau des cellules interpolees prenant en compte eventuellement le double_wall pour la zone donneuse noz
    // Creation de tableaux pour stocker les indices, les coefficients d'interpolation, 
    // un tableau d'indirection pour les points EX et les volumes des cellules d'interpolation
    E_Int sizeEachIndirTab = interpData->_sizeOfEachIndirectionTab[noz];
    FldArrayF* donorInterpCoef = new FldArrayF(sizeEachIndirTab,ncoefs);
    FldArrayI* cellIndicesRcv = new FldArrayI(sizeEachIndirTab);
    FldArrayI* cellIndicesDonor = new FldArrayI(sizeEachIndirTab);
    FldArrayI* cellIndicesDonorExtC = new FldArrayI(sizeEachIndirTab);
    FldArrayI* EXDir = new FldArrayI(sizeEachIndirTab);
    FldArrayF* volCell = new FldArrayF(sizeEachIndirTab);
    FldArrayI* interpTypeCell = new FldArrayI(sizeEachIndirTab);

    // pointeurs sur le champ donorInterpCoef de coefficients d'interpolation
    E_Float* donorInterpCoef1=donorInterpCoef->begin(1); 
    E_Float* donorInterpCoef2=donorInterpCoef->begin(2); 
    E_Float* donorInterpCoef3=donorInterpCoef->begin(3); 
    E_Float* donorInterpCoef4=donorInterpCoef->begin(4); 
    E_Float* donorInterpCoef5=donorInterpCoef->begin(5); 
    E_Float* donorInterpCoef6=donorInterpCoef->begin(6); 
    E_Float* donorInterpCoef7=donorInterpCoef->begin(7); 
    E_Float* donorInterpCoef8=donorInterpCoef->begin(8);
    // pointeurs sur le champ cellIndices
    E_Int* cellIndicesRcv1D = cellIndicesRcv->begin(); 
    E_Int* cellIndicesDonorExtC1D = cellIndicesDonorExtC->begin();
    E_Int* cellIndicesDonor1D = cellIndicesDonor->begin();

    // objets Python pour stocker les coefficients d'interpolation, les indices des cellules interpolees,  
    // les indices des cellules donneuses, les indices des cellules orphelines

    // indices 1D receveur et donneur 
    E_Int indRcv1D;
    // Affectation des tableaux
    for (E_Int i=nb; i < nb+interpData->_sizeOfEachIndirectionTab[noz]; i++)
    {
      li = i - nb;
      // cellules interpolees
      indRcv1D = indRcv[interpData->_indirectionPerDom[i]];
      cellIndicesRcv1D[li] = indRcv1D;
      // indirection pour les interfaces
      if (isEX) (*EXDir)[li] = (E_Int)(*vectOfEXdir[noz])[interpData->_indirectionPerDom[i]];
      // cellule d interpolation (indices et volume)
      cellIndicesDonor1D[li] = interpData->_cells[interpData->_indirectionPerDom[i]]; 
      cellIndicesDonorExtC1D[li] = donorCells[interpData->_indirectionPerDom[i]];
      (*volCell)[li]= interpData->_vols[interpData->_indirectionPerDom[i]];
      (*interpTypeCell)[li]= interpData->_interpTypes[interpData->_indirectionPerDom[i]];
      // cellule d interpolation (coefficients)
      donorInterpCoef1[li] = interpData->_coefs(interpData->_indirectionPerDom[i],1);
      donorInterpCoef2[li] = interpData->_coefs(interpData->_indirectionPerDom[i],2);
      donorInterpCoef3[li] = interpData->_coefs(interpData->_indirectionPerDom[i],3);
      donorInterpCoef4[li] = interpData->_coefs(interpData->_indirectionPerDom[i],4);
      donorInterpCoef5[li] = interpData->_coefs(interpData->_indirectionPerDom[i],5);
      donorInterpCoef6[li] = interpData->_coefs(interpData->_indirectionPerDom[i],6);
      donorInterpCoef7[li] = interpData->_coefs(interpData->_indirectionPerDom[i],7);
      donorInterpCoef8[li] = interpData->_coefs(interpData->_indirectionPerDom[i],8);
    }
    // Marquage des cellules orphelines
    FldArrayI orphanPoints(orphanVect.size());
    c = 0;
    for (unsigned int i= 0; i < orphanVect.size(); i++)
    {
      E_Int ind = orphanVect[i];
      E_Int indOrphan = indRcv[ind];
      orphanPoints[c] = indOrphan; c++;
    }

    // construction de PyArray
    // -------------------------------------
    //    coefficients d'interpolation
    PyObject* fout = K_NUMPY::buildNumpyArray(*donorInterpCoef,1);
    PyList_Append(PyListCoefficients, fout); Py_DECREF(fout);
    delete donorInterpCoef;

    //   indices des cellules receveuses
    PyObject* cellIndR = K_NUMPY::buildNumpyArray(*cellIndicesRcv,1);
    PyList_Append(PyListCellIndicesR, cellIndR); Py_DECREF(cellIndR);
    delete cellIndicesRcv;   

    //   indices des cellules donneuses
    PyObject* cellIndD = K_NUMPY::buildNumpyArray(*cellIndicesDonor,1);
    PyList_Append(PyListCellIndicesD, cellIndD); Py_DECREF(cellIndD);
    delete cellIndicesDonor;  

    //   indices des cellules donneuses en centres etendus
    PyObject* cellIndDExtC = K_NUMPY::buildNumpyArray(*cellIndicesDonorExtC,1);
    PyList_Append(PyListCellIndicesDExtC, cellIndDExtC); Py_DECREF(cellIndDExtC);
    delete cellIndicesDonorExtC;  

    //   volumes des cellules donneuses
    PyObject* vout = K_NUMPY::buildNumpyArray(*volCell,1);
    PyList_Append(PyListVolumes, vout); Py_DECREF(vout);
    delete volCell;

    //   interpType des cellules donneuses
    PyObject* tout = K_NUMPY::buildNumpyArray(*interpTypeCell,1);
    PyList_Append(PyListInterpTypes, tout); Py_DECREF(tout);
    delete interpTypeCell;

    //   indices des points extrapoles
    PyObject* cellIndE = K_NUMPY::buildNumpyArray(*extrapPtsP,1);
    PyList_Append(PyListCellIndicesE, cellIndE); Py_DECREF(cellIndE);
    delete extrapPtsP;

    //   indices des points orphelins
    PyObject* cellIndO = K_NUMPY::buildNumpyArray(orphanPoints,1);
    PyList_Append(PyListCellIndicesO, cellIndO); Py_DECREF(cellIndO);
    //delete orphanPoints;

    if (isEX)
    {
      PyObject* pyEXDir = K_NUMPY::buildNumpyArray(*EXDir,1);
      PyList_Append(PyListEXDir, pyEXDir); Py_DECREF(pyEXDir);
    }
    delete EXDir;

    // incrementation du nombre de points interpoles
    nb = nb + interpData->_sizeOfEachIndirectionTab[noz];
  }
  // Construction de l'objet tpl retourne par la methode (liste de 2 listes)
  PyObject* tpl;
  if (isEX)
    tpl = Py_BuildValue("[OOOOOOOOO]", PyListCellIndicesR, PyListCellIndicesD, PyListCellIndicesDExtC,
                        PyListCoefficients, PyListVolumes, PyListInterpTypes, PyListCellIndicesE, 
                        PyListCellIndicesO, PyListEXDir);
  else
    tpl = Py_BuildValue("[OOOOOOOO]", PyListCellIndicesR, PyListCellIndicesD, PyListCellIndicesDExtC,
                        PyListCoefficients, PyListVolumes, PyListInterpTypes, PyListCellIndicesE, 
                        PyListCellIndicesO);
  
  Py_DECREF(PyListCoefficients); Py_DECREF(PyListCellIndicesR); 
  Py_DECREF(PyListCellIndicesD); Py_DECREF(PyListCellIndicesDExtC); 
  Py_DECREF(PyListEXDir);  Py_DECREF(PyListCellIndicesE); Py_DECREF(PyListCellIndicesO); 
  Py_DECREF(PyListVolumes); Py_DECREF(PyListInterpTypes);

  /*----------------------------------------------------------*/
  /* Destruction des objets crees                             */
  /*----------------------------------------------------------*/
  K_ARRAY::cleanStructFields(vectOfInterpPts); 
  for (unsigned int noz = 0; noz < vectOfIndRcv.size(); noz++)
    delete vectOfIndRcv[noz];
  if (isEX) K_ARRAY::cleanStructFields(vectOfEXdir);
    
  E_Int sizev = listOfExtCentersKMeshes.size();
  for (E_Int v = 0; v < sizev; v++)
  {delete listOfKMeshesN[v]; delete listOfExtCentersKMeshes[v]; delete listOfInterpDatas[v];} 
  
  interpData->_cells.malloc(0);
  interpData->_vols.malloc(0);
  interpData->_coefs.malloc(0);
  interpData->_interpTypes.malloc(0);
  interpData->_indirectionPerDom.malloc(0);
  interpData->_sizeOfEachIndirectionTab.malloc(0); 
  delete interpData; 
  nzones = coords.size();
  for (E_Int v = 0; v < nzones; v++)
  {delete coords[v]; delete cellns[v];}
  
  return tpl;
}
//=============================================================================
/* Calcule et stocke les coeff d'interpolations pour les points EX */
//=============================================================================
void K_CONNECTOR::compAndStoreEXInterpCoefs(
  E_Float geomCutOff,
  K_KINTERP::BlkInterpData::InterpolationType interpType0,
  K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
  vector<FldArrayF*>& coordEX, 
  InterpData* interpDataEX, FldArrayI& donorCellsEX, 
  vector< vector<E_Int> >& vectOfExtrapPts, vector<E_Int>& orphanVect,
  vector<E_Int>& nit, vector<E_Int>& njt, vector<E_Int>& nkt,
  vector<FldArrayF*>& coords,
  vector<K_KINTERP::BlkInterpData*>& listOfInterpDatas, 
  vector<FldArrayI*>& cellns, E_Int zid, E_Float cfMax)
{
  E_Int nbI = coordEX[0]->getSize();
  // Array that keep the number of points that are interpolated from a domain
  E_Int numberOfDomains = listOfInterpDatas.size();
  FldArrayI usedDomains(numberOfDomains);
  E_Int* usedDomp = usedDomains.begin();
  usedDomains.setAllValuesAtNull();
  E_Int blocNbSav = 0;
  E_Int blocNbSavExtrap = 0;
  E_Int ncf = 8;
  E_Int nindi = 7;// nfld du tableau indit contenant les indices des cellules d interp + ordre
  E_Float x, y, z;
  short found, foundSav;
  E_Float volSav,volSavExtrap, vol;
  //E_Int interpType;
  //E_Int interpTypeSavExtrap=0;
  //E_Int interpTypeSav=0;
  E_Int ind0, ind0ext;
  E_Float nature;
  E_Int cellSav = 0;
  E_Int cellExtSav = 0;
  E_Int indgSav = 0;
  E_Int indg;
  E_Float test; 
  
  E_Int size = nbI*numberOfDomains;
  E_Int numberOfInterpolatedPoints = 0;
  E_Int numberOfExtrapolatedPoints = 0;
  E_Int numberOfPbPoints = 0;
  E_Int* indirection = new E_Int[size];
  
  // Temp array for vectorization
  FldArrayF coordv(size, 3);
  E_Float* coordx = coordv.begin(1);
  E_Float* coordy = coordv.begin(2);
  E_Float* coordz = coordv.begin(3);
  FldArrayIS foundv(size);
  FldArrayI extrapv(size); extrapv.setAllValuesAt(1);
  FldArrayI indit(size,nindi);
  FldArrayI icv(size);
  FldArrayI jcv(size);
  FldArrayI kcv(size);
  FldArrayF cfv(size, ncf);
  foundv.setAllValuesAtNull();

  // working vectors for periodic interpolation
  vector<E_Int> thetaPEx;
  vector<E_Int> thetaMEx;
  vector<E_Int> thetaPBlk;
  vector<E_Int> thetaMBlk;

  // Allocation
  interpDataEX->_coefs.malloc(nbI, ncf); interpDataEX->_coefs.setAllValuesAtNull();
  interpDataEX->_cells.malloc(nbI);interpDataEX->_cells.setAllValuesAtNull();
  interpDataEX->_vols.malloc(nbI);interpDataEX->_vols.setAllValuesAtNull();
  interpDataEX->_interpTypes.malloc(nbI);interpDataEX->_interpTypes.setAllValuesAtNull();

  E_Int* cellspEX = interpDataEX->_cells.begin();
  E_Float* volspEX = interpDataEX->_vols.begin();
  E_Int* donorCellsEXp = donorCellsEX.begin();
  E_Int* interpTypespEX = interpDataEX->_interpTypes.begin();

  interpDataEX->_indirectionPerDom.malloc(nbI);interpDataEX->_indirectionPerDom.setAllValuesAtNull();
  vector<E_Int> wallPoints;
  vector<E_Float> dist;
  // Calcul de coordv pour vectoriser
  for (E_Int blocNb = 0; blocNb < numberOfDomains; blocNb++)
  {
    E_Float* xc = coordEX[blocNb]->begin(1);
    E_Float* yc = coordEX[blocNb]->begin(2);
    E_Float* zc = coordEX[blocNb]->begin(3);
    for (E_Int i = 0; i < nbI; i++)
    {
      E_Int ind = i;
      E_Int ind2 = i +  nbI*blocNb;
      coordx[ind2] = xc[ind];
      coordy[ind2] = yc[ind];
      coordz[ind2] = zc[ind];
    }
  }
  // Compute interpolation coefficients for all interpolated cells from all interpolation blocks
  for (E_Int blocNb = 0; blocNb < numberOfDomains; blocNb++)
  {
    K_KINTERP::BlkInterpData* blkInterpData = listOfInterpDatas[blocNb];

    // Get the interpolation cell of vector of interpolated points
    // Retourne directement les indices dans indit sous forme de centres
    blkInterpData->getInterpolationCellStructv(coordv, nbI*blocNb, nbI*(blocNb+1),
                                               indit, icv, jcv, kcv, cfv, foundv, extrapv, 
                                               interpMeshType, interpType0);
  }  

  /*-----------------*/
  /* Keep the bests  */
  /*-----------------*/
  E_Int ic, jc, kc;
  E_Int ni1, nj1, nk1, ni1nj1;
  FldArrayF cfSav(ncf);
  FldArrayI indi(nindi);
  FldArrayF cf(ncf);
  E_Int ret = 1;

  for (E_Int i = 0; i < nbI; i++) // for all points that require interpolation
  {
    volSav = K_CONST::E_MAX_FLOAT; vol = volSav;
    foundSav = 0;
    for (E_Int blocNb = 0; blocNb < numberOfDomains; blocNb++)
    {
      FldArrayI& cellNatureField = *cellns[blocNb];
      ni1 = nit[blocNb]-1; nj1 = njt[blocNb]-1; nk1 = nkt[blocNb]-1; ni1nj1 = ni1*nj1;
      indg = i+blocNb*nbI; found = foundv[indg];
      x = coordx[indg]; y = coordy[indg]; z = coordz[indg];
     
      E_Float* xni = coords[blocNb]->begin(1);
      E_Float* yni = coords[blocNb]->begin(2);
      E_Float* zni = coords[blocNb]->begin(3);

      if (found != 0) 
      {
        for (E_Int no = 1; no <= nindi; no++)
          indi[no-1] = indit(indg,no);
        for (E_Int nocf = 1; nocf <= ncf; nocf++)
          cf[nocf-1] = cfv(indg,nocf);

        //interpType = indi[0];
        ind0 = indi[1] + indi[3] * ni1 + indi[5] * ni1nj1;
        ind0ext = icv[indg]-1+(jcv[indg]-1)*(ni1+2)+(kcv[indg]-1)*(ni1+2)*(nj1+2);
               
        K_METRIC::compVolOfStructCell3D(
          nit[blocNb], njt[blocNb], nkt[blocNb], ind0, -1,
          xni, yni, zni, vol);
        ret = extrapv[indg]; if ( ret == 0 ) vol = vol + 500;

        if (vol < volSav)
        {
          // test interpolation pt EX : 
          // - pas de cellule masquee dans la molecule d interpolation, sauf si son cf associe est nul.
          // - somme des cf = 1.
          // -----------------
          // 1/2*cellN*(3-cellN) renvoie 0 si cellN = 0 (pt masque) et 1 si cellN =1 ou 2 (pt calcule ou interpole)
          nature = K_KINTERP::compInterpolatedNatureEX(ni1, nj1, nk1, indi.getSize(), indi.begin(), cf,
                                                       cellNatureField.begin());
          
          if (K_FUNC::fEqual(nature,K_CONST::ONE,geomCutOff) == true) // pas de pt masque dans la cellule d'interp et somme des cf = 1
          {
            indirection[blocNb*nbI+usedDomp[blocNb]] = i;
            volSav = vol; //interpTypeSav = interpType;
            blocNbSav = blocNb;
            indgSav = indg;
            foundSav = found;
            cellSav = ind0;
            cellExtSav = ind0ext;//cellule d interpolation en centres etendus
          }
        }
      }
    }// tous les blocs d interpolation 

    switch (foundSav) // recuperation des points periodiques
    {
      case 1:
        cellspEX[i] = cellSav;
        volspEX[i]=volSav;
        interpTypespEX[i]= 100;//interpTypeSav;
        donorCellsEXp[i] = cellExtSav;
        for (E_Int nf = 1; nf <= ncf; nf++)
          interpDataEX->_coefs(i, nf) = cfv(indgSav, nf);

        usedDomp[blocNbSav]++;
        numberOfInterpolatedPoints++;
        goto endPoint;
      
      case 2: //cellule d'interpolation en +theta
        cellspEX[i] = cellSav;
        volspEX[i]=volSav;
        interpTypespEX[i]=102;//interpTypeSav;
        donorCellsEXp[i] = cellExtSav;
        for (E_Int nf = 1; nf <= ncf; nf++)
          interpDataEX->_coefs(i, nf) = cfv(indgSav, nf);

        thetaPEx.push_back(i);
        thetaPBlk.push_back(blocNbSav);
        
        usedDomp[blocNbSav]++;
        numberOfInterpolatedPoints++;
        goto endPoint;
        
      case 3:
        cellspEX[i] = cellSav;
        volspEX[i]=volSav;
        interpTypespEX[i]=103;//interpTypeSav;
        donorCellsEXp[i] = cellExtSav;
        for (E_Int nf = 1; nf <= ncf; nf++)
          interpDataEX->_coefs(i, nf) = cfv(indgSav, nf);

        thetaMEx.push_back(i);
        thetaMBlk.push_back(blocNbSav);
        
        usedDomp[blocNbSav]++;
        numberOfInterpolatedPoints++;
        goto endPoint;
        
      default: // foundSav = 0
        // Look for extrapolation cell on all blocks
        volSavExtrap = K_CONST::E_MAX_FLOAT;
        for (E_Int blocNb = 0; blocNb <  numberOfDomains; blocNb++)
        {
          indg = i+blocNb*nbI;
          ni1 = nit[blocNb]-1; nj1 = njt[blocNb]-1; nk1 = nkt[blocNb]-1; ni1nj1 = ni1*nj1;
          FldArrayI& cellNatureField = *cellns[blocNb];
          K_KINTERP::BlkInterpData* blkInterpData = listOfInterpDatas[blocNb];

          x = coordx[indg]; y = coordy[indg]; z = coordz[indg];
          found = solveOrphanPoints(interpType0, interpMeshType, blkInterpData, 
                                    nit[blocNb], njt[blocNb], nkt[blocNb], 
                                    cellNatureField, x, y, z, 1, cfMax, cf, indi, ic, jc, kc);
          if (found != 0)
          {
            //interpType = indi[0]; 
            ind0 = indi[1] + indi[3] * ni1 + indi[5] * ni1nj1;
            ind0ext = ic-1+(jc-1)*(ni1+2)+(kc-1)*(ni1+2)*(nj1+2);
         
            E_Float* xni = coords[blocNb]->begin(1);
            E_Float* yni = coords[blocNb]->begin(2);
            E_Float* zni = coords[blocNb]->begin(3);
            K_METRIC::compVolOfStructCell3D(
              nit[blocNb], njt[blocNb], nkt[blocNb], ind0, -1,
              xni, yni, zni, vol);
            ret = extrapv[indg]; 
            if ( ret == 0 ) vol = vol + 1000;

            if (vol < volSavExtrap)
            {
              volSavExtrap = vol;
              blocNbSavExtrap = blocNb;
              //interpTypeSavExtrap = interpType;// extrapolation
              foundSav = found;
              cellSav =  ind0;//cellule d interpolation en centres
              cellExtSav = ind0ext;//cellule d interpolation en centres etendus
              cfSav = cf;
            }
          }
        }
        switch (foundSav)
        {
          case 1:
            cellspEX[i] = cellSav;
            volspEX[i]=volSavExtrap;
            interpTypespEX[i]=100;//interpTypeSavExtrap;
            donorCellsEXp[i] = cellExtSav;
            for (E_Int nf = 1; nf <= ncf; nf++)
              interpDataEX->_coefs(i, nf) = cfSav[nf-1];

            indirection[blocNbSavExtrap*nbI+usedDomp[blocNbSavExtrap]] = i;
            usedDomp[blocNbSavExtrap]++;
            vectOfExtrapPts[blocNbSavExtrap].push_back(i);
            numberOfExtrapolatedPoints++;
            goto endPoint;
              
          case 2:
            cellspEX[i] = cellSav;
            volspEX[i]=volSavExtrap;
            interpTypespEX[i]=102;//interpTypeSavExtrap;
            donorCellsEXp[i] = cellExtSav;
            for (E_Int nf = 1; nf <= ncf; nf++)
              interpDataEX->_coefs(i, nf) = cfSav[nf-1];
                        
            thetaPEx.push_back(i);
            thetaPBlk.push_back(blocNbSavExtrap);
            
            indirection[blocNbSavExtrap*nbI+usedDomp[blocNbSavExtrap]] = i;
            usedDomp[blocNbSavExtrap]++;
            vectOfExtrapPts[blocNbSavExtrap].push_back(i);
            numberOfExtrapolatedPoints++;
            goto endPoint;
              
          case 3:
            cellspEX[i] = cellSav;
            volspEX[i]=volSavExtrap;
            interpTypespEX[i]=103;//interpTypeSavExtrap;
            donorCellsEXp[i] = cellExtSav;
           for (E_Int nf = 1; nf <= ncf; nf++)
              interpDataEX->_coefs(i, nf) = cfSav[nf-1];
   
            thetaMEx.push_back(i);
            thetaMBlk.push_back(blocNbSavExtrap);              
            
            indirection[blocNbSavExtrap*nbI+usedDomp[blocNbSavExtrap]] = i;
            usedDomp[blocNbSavExtrap]++;
            vectOfExtrapPts[blocNbSavExtrap].push_back(i);
            numberOfExtrapolatedPoints++;
            goto endPoint;
          default:
            break;
        }
        break;
    }//switch foundSav

    // Look for extrapolation cell on all blocks
    volSavExtrap = K_CONST::E_MAX_FLOAT;
    for (E_Int blocNb = 0; blocNb < numberOfDomains; blocNb++)
    {
      FldArrayI& cellNatureField = *cellns[blocNb];
      K_KINTERP::BlkInterpData* blkInterpData = listOfInterpDatas[blocNb];
      x = coordx[i+nbI*blocNb];
      y = coordy[i+nbI*blocNb];
      z = coordz[i+nbI*blocNb];
      test = 1.;
      found = blkInterpData->getExtrapolationCell(x, y, z,ic, jc, kc, cf, 
                                                  cellNatureField, 1, test, interpType0, interpMeshType, cfMax);

      if (found != 0)
      {
        E_Float* xni = coords[blocNb]->begin(1);
        E_Float* yni = coords[blocNb]->begin(2);
        E_Float* zni = coords[blocNb]->begin(3);

        ni1 = nit[blocNb]-1; nj1 = njt[blocNb]-1; nk1 = nkt[blocNb]-1; ni1nj1 = ni1*nj1;
        ret = blkInterpData->fromExtendedToStandardCenters(ic, jc, kc, indi, interpType0); 
        //interpType = indi[0];
        ind0 = indi[1] + indi[3]*ni1 + indi[5]*ni1nj1;
        ind0ext = ic-1+(jc-1)*(ni1+2)+(kc-1)*(ni1+2)*(nj1+2);

        K_METRIC::compVolOfStructCell3D(
          nit[blocNb], njt[blocNb], nkt[blocNb], ind0, -1,
          xni, yni, zni, vol);
        if ( ret == 0 ) vol = vol + 1000*test;

        if (vol < volSavExtrap)
        {
          volSavExtrap = vol;
          //interpTypeSavExtrap = interpType + 100 ;//extrapolation
          blocNbSavExtrap = blocNb;
          foundSav = found;
          cellSav = ind0;
          cellExtSav = ind0ext;//cellule d interpolation en centres etendus
          cfSav = cf;
        }
      }
    }
      
    switch (foundSav)
    {
      case 1: // cellule d'extrapolation du bloc initial 
        cellspEX[i] = cellSav;
        volspEX[i]=volSavExtrap;
        interpTypespEX[i]=100;//interpTypeSavExtrap;
        donorCellsEXp[i] = cellExtSav;
        for (E_Int nf = 1; nf <= ncf; nf++)
          interpDataEX->_coefs(i, nf) = cfSav[nf-1];

        indirection[blocNbSavExtrap*nbI+usedDomp[blocNbSavExtrap]] = i;
        usedDomp[blocNbSavExtrap]++;
        vectOfExtrapPts[blocNbSavExtrap].push_back(i);
        numberOfExtrapolatedPoints++;
        goto endPoint;
        
      case  2: // cellule d'extrapolation du bloc en +theta
        cellspEX[i] = cellSav;
        volspEX[i]=volSavExtrap;
        interpTypespEX[i]=102;//interpTypeSavExtrap;
        donorCellsEXp[i] = cellExtSav;
        for (E_Int nf = 1; nf <= ncf; nf++)
          interpDataEX->_coefs(i, nf) = cfSav[nf-1];

        thetaPEx.push_back(i);
        thetaPBlk.push_back(blocNbSavExtrap);
        
        indirection[blocNbSavExtrap*nbI+usedDomp[blocNbSavExtrap]] = i;
        usedDomp[blocNbSavExtrap]++;
        vectOfExtrapPts[blocNbSavExtrap].push_back(i);
        numberOfExtrapolatedPoints++;
        goto endPoint;
        
      case 3: // cellule d'extrapolation du bloc en -theta
        cellspEX[i] = cellSav;
        volspEX[i]=volSavExtrap;
        interpTypespEX[i]=103;//interpTypeSavExtrap;
        donorCellsEXp[i] = cellExtSav;
        for (E_Int nf = 1; nf <= ncf; nf++)
          interpDataEX->_coefs(i, nf) = cfSav[nf-1];

        thetaMEx.push_back(i);
        thetaMBlk.push_back(blocNbSavExtrap);              
        
        indirection[blocNbSavExtrap*nbI+usedDomp[blocNbSavExtrap]] = i;
        usedDomp[blocNbSavExtrap]++;
        vectOfExtrapPts[blocNbSavExtrap].push_back(i);
        numberOfExtrapolatedPoints++;
        goto endPoint;
      default:
        break;
    }
    
    // BIG TROUBLE: cannot interpolate nor extrapolate
    orphanVect.push_back(i);
    numberOfPbPoints++;
    endPoint:;
  } /* all points */
  

  // Recup'
  E_Int compt = 0;
  E_Int* indir1 = interpDataEX->_indirectionPerDom.begin();
  for (E_Int blocNb = 0; blocNb < numberOfDomains; blocNb++)
  {   
    for (E_Int i = 0; i < usedDomp[blocNb]; i++)
    {
      indir1[compt] = indirection[blocNb*nbI+i];
      compt++;
    }
  } 
  delete [] indirection;

  //recuperation des points interpol�s par un bloc periodique
  E_Int sizeP = thetaPEx.size();
  E_Int sizeM = thetaMEx.size();
  FldArrayI periodicEXPtsP(sizeP);
  FldArrayI periodicEXPtsM(sizeM);
  FldArrayI periodicEXBlksP(sizeP);
  FldArrayI periodicEXBlksM(sizeM);

  for (E_Int ii = 0; ii < sizeP; ii++)
  {
    periodicEXPtsP[ii]  = thetaPEx[ii];
    periodicEXBlksP[ii] = thetaPBlk[ii];
  }
  for (E_Int ii = 0; ii < sizeM; ii++)
  {
    periodicEXPtsM[ii]  = thetaMEx[ii];
    periodicEXBlksM[ii] = thetaMBlk[ii];
  }
  
  // // Resume de la situation
  // if (zid == -1)
  //   printf("EX interpolated=%d, extrapolated=%d", numberOfInterpolatedPoints, numberOfExtrapolatedPoints);
  // else
  //   printf("Zone %d: EX interpolated=%d, extrapolated=%d", zid, numberOfInterpolatedPoints, numberOfExtrapolatedPoints);
  // // periodic treatment 
  // if (sizeP != 0)
  //   printf(" Periodic EX points  from +theta block=%d", sizeP);
  // if (sizeM != 0)
  //   printf(" Periodic EX points from -theta block=%d", sizeM);
  // printf(", critical=%d\n", numberOfPbPoints);

  // if (numberOfPbPoints > 0)
  //   printf("# WARNING ! there are EX orphan points !\n");

  interpDataEX->_sizeOfEachIndirectionTab.malloc(numberOfDomains);
  for (E_Int i = 0; i < numberOfDomains; i++)
    interpDataEX->_sizeOfEachIndirectionTab[i] = usedDomp[i];
}

//=============================================================================
/* Calcule et stocke les coeff d interpolations pour les centres
   le vecteur vectOfExtrapPts doit etre dimensionne au prealable par le nombre
   de zones donneuses */
//=============================================================================
void K_CONNECTOR::compAndStoreInterpCoefs(
  E_Float geomCutOff, 
  K_KINTERP::BlkInterpData::InterpolationType interpType0,
  K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
  vector<FldArrayF*>&  interpolatedCenters,
  InterpData* interpData, FldArrayI& donorCells, vector<vector<E_Int> >& vectOfExtrapPts, 
  vector<E_Int>& orphanVect,
  vector<E_Int>& nit, vector<E_Int>& njt,vector<E_Int>& nkt,
  vector<FldArrayF*>& coords,
  vector<K_KINTERP::BlkInterpData*>& listOfInterpDatas, 
  vector<FldArrayI*>& cellns, E_Int zid, E_Float cfMax)
{
  E_Int nbI = interpolatedCenters[0]->getSize();
  // Array that keep the number of points that are interpolated from a domain
  E_Int numberOfDomains = listOfInterpDatas.size();
  FldArrayI usedDomains(numberOfDomains);
  E_Int* usedDomp = usedDomains.begin();
  usedDomains.setAllValuesAtNull();
  
  E_Int blocNbSavExplicit = 0;
  E_Int blocNbSavExtrap = 0;

  E_Int ncf=8;
  E_Int nindi=7;// nfld du tableau indit contenant les indices des cellules d interp + ordre
  E_Float x, y, z;
  short found, foundSav;
  E_Float volSavExplicit, volSavExtrap, vol;
  //E_Int interpType;
  //E_Int interpTypeSav=0, interpTypeSavExtrap=0;
  E_Int ind0, ind0ext;
  E_Float nature;
  E_Int cellSav = 0;
  E_Int cellExtSav = 0;
  E_Int indg;
  E_Float test; 
  
  E_Int size = nbI*numberOfDomains;
  E_Int numberOfInterpolatedPoints = 0;
  E_Int numberOfExtrapolatedPoints = 0;
  E_Int numberOfPbPoints = 0;
  E_Int* indirection = new E_Int[size];
  
  // Temp array for vectorization
  FldArrayF coordv(size, 3);
  E_Float* coordx = coordv.begin(1);
  E_Float* coordy = coordv.begin(2);
  E_Float* coordz = coordv.begin(3);
  FldArrayIS foundv(size);  foundv.setAllValuesAtNull();
  FldArrayI indit(size,nindi);
  FldArrayI icv(size);
  FldArrayI jcv(size);
  FldArrayI kcv(size);
  FldArrayF cfv(size, ncf);
  FldArrayI extrapv(size); extrapv.setAllValuesAt(1);

  // working vectors for periodic interpolation
  vector<E_Int> thetaPPoints;
  vector<E_Int> thetaMPoints;
  vector<E_Int> thetaPBlk;
  vector<E_Int> thetaMBlk;

  // Allocation
  interpData->_coefs.malloc(nbI, ncf); interpData->_coefs.setAllValuesAtNull();
  interpData->_cells.malloc(nbI); interpData->_cells.setAllValuesAtNull();
  interpData->_vols.malloc(nbI); interpData->_vols.setAllValuesAtNull();
  interpData->_interpTypes.malloc(nbI); interpData->_interpTypes.setAllValuesAtNull();

  E_Int* cellsp = interpData->_cells.begin();
  E_Float* volsp = interpData->_vols.begin();
  E_Int* donorCellsp = donorCells.begin();
  E_Int* interpTypesp = interpData->_interpTypes.begin();

  interpData->_indirectionPerDom.malloc(nbI);interpData->_indirectionPerDom.setAllValuesAtNull();
  vector<E_Int> wallPoints;
  vector<E_Float> dist;
  
  // Compute interpolation coefficients for all interpolated cells from all interpolation blocks
  for (E_Int blocNb = 0; blocNb < numberOfDomains; blocNb++)
  {
    E_Float* xc = interpolatedCenters[blocNb]->begin(1);
    E_Float* yc = interpolatedCenters[blocNb]->begin(2);
    E_Float* zc = interpolatedCenters[blocNb]->begin(3);

    for (E_Int ind = 0; ind < nbI; ind++)
    {
      E_Int ind2 = ind +  nbI*blocNb;
      coordx[ind2] = xc[ind];
      coordy[ind2] = yc[ind];
      coordz[ind2] = zc[ind];
    }
    K_KINTERP::BlkInterpData* blkInterpData = listOfInterpDatas[blocNb];
    // Get the interpolation cell of vector of interpolated points
    // Retourne directement les indices dans indit sous forme de centres
    blkInterpData->getInterpolationCellStructv(coordv, nbI*blocNb, nbI*(blocNb+1),
                                               indit, icv, jcv, kcv, cfv, foundv, extrapv, 
                                               interpMeshType, interpType0);
  }  
  /*-----------------*/
  /* Keep the bests  */
  /*-----------------*/
  E_Int ic, jc, kc;
  E_Int ni1, nj1, nk1, ni1nj1;
  FldArrayF cfSav(ncf);
  FldArrayI indi(nindi);
  FldArrayF cf(ncf);
  E_Int ret = 1;

  for (E_Int i = 0; i < nbI; i++) // for all points that require interpolation
  {
    volSavExplicit = K_CONST::E_MAX_FLOAT; 
    volSavExtrap = K_CONST::E_MAX_FLOAT;
    foundSav = 0;
    for (E_Int blocNb = 0; blocNb < numberOfDomains; blocNb++)
    {
      FldArrayI& cellNatureField = *cellns[blocNb];
      ni1 = nit[blocNb]-1; nj1 = njt[blocNb]-1; nk1 = nkt[blocNb]-1; ni1nj1 = ni1*nj1;
      indg = i+blocNb*nbI; 
      x = coordx[indg]; y = coordy[indg]; z = coordz[indg];
      ni1nj1 = ni1*nj1;
      found = foundv[indg];
      E_Float* xni = coords[blocNb]->begin(1);
      E_Float* yni = coords[blocNb]->begin(2);
      E_Float* zni = coords[blocNb]->begin(3);
 
      if (found != 0) 
      {
        for (E_Int no = 1; no <= nindi; no++) indi[no-1] = indit(indg,no);
        for (E_Int nocf = 1; nocf <= ncf; nocf++) cf[nocf-1] = cfv(indg,nocf);
        //interpType = indi[0];
        ind0 = indi[1] + indi[3] * ni1 + indi[5] * ni1nj1;
        ind0ext = icv[indg]-1+(jcv[indg]-1)*(ni1+2)+(kcv[indg]-1)*(ni1+2)*(nj1+2);       
        K_METRIC::compVolOfStructCell3D(
          nit[blocNb], njt[blocNb], nkt[blocNb], ind0, -1,
          xni, yni, zni, vol);
        ret = extrapv[indg]; if (ret == 0) vol += 500.;

        if (vol < volSavExplicit)
        { 
          nature = K_KINTERP::compInterpolatedNature(ni1, nj1, nk1, indi.getSize(), indi.begin(), cf, cellNatureField.begin());
          if (K_FUNC::fEqual(nature,K_CONST::ONE,geomCutOff) == true) // pas de pt masque ou interpole dans la cellule d'interp et somme des cf = 1
          {
            indirection[blocNb*nbI+usedDomp[blocNb]] = i;
            volSavExplicit = vol; //interpTypeSav = interpType;
            blocNbSavExplicit = blocNb;
            //indgSav = indg;
            foundSav = found;
            cellSav = ind0;//cellule d'interpolation en centres
            cellExtSav = ind0ext;//cellule d'interpolation en centres etendus
            cfSav = cf;
          }
        }
      }
    } /* all domains */
    // Interpolable point, explicit
    //printf("%d: found interp=%d\n", indg, foundSav);
    switch (foundSav) // recuperation des points periodiques
    {
      case 1:
        cellsp[i] = cellSav;
        volsp[i] = volSavExplicit;
        interpTypesp[i]=100;//interpTypeSav;
        donorCellsp[i] = cellExtSav;
        for (E_Int nf = 1; nf <= ncf; nf++)
          interpData->_coefs(i, nf) = cfSav[nf-1];
        
        usedDomp[blocNbSavExplicit]++;
        numberOfInterpolatedPoints++;
        goto endPoint;

      case 2: //cellule d'interpolation en +theta
        cellsp[i] = cellSav;
        volsp[i] = volSavExplicit;
        interpTypesp[i]=102;//interpTypeSav;
        donorCellsp[i] = cellExtSav;
        for (E_Int nf = 1; nf <= ncf; nf++) 
      	  interpData->_coefs(i, nf) = cfSav[nf-1];
        
        thetaPPoints.push_back(i);
        thetaPBlk.push_back(blocNbSavExplicit);

        usedDomp[blocNbSavExplicit]++;
        numberOfInterpolatedPoints++;
        goto endPoint;
        
      case 3: //cellule d'interpolation en -theta
        cellsp[i] = cellSav;
        volsp[i] = volSavExplicit;
        interpTypesp[i]=103;//interpTypeSav;
        donorCellsp[i] = cellExtSav;
        for (E_Int nf = 1; nf <= ncf; nf++)
          interpData->_coefs(i, nf) = cfSav[nf-1];
        
        thetaMPoints.push_back(i);
        thetaMBlk.push_back(blocNbSavExplicit);

        usedDomp[blocNbSavExplicit]++;
        numberOfInterpolatedPoints++;
        goto endPoint;

      default: // foundSav = 0 
        // Interpolable, extrapolated
        // Look for extrapolation cell on all blocks
        volSavExtrap = K_CONST::E_MAX_FLOAT;
        for (E_Int blocNb = 0; blocNb < numberOfDomains; blocNb++)
        {
          FldArrayI& cellNatureField = *cellns[blocNb];
          K_KINTERP::BlkInterpData* blkInterpData = listOfInterpDatas[blocNb];

          indg = i+blocNb*nbI;
          ni1 = nit[blocNb]-1; nj1 = njt[blocNb]-1; nk1 = nkt[blocNb]-1; ni1nj1 = ni1*nj1;
          x = coordx[indg]; y = coordy[indg]; z = coordz[indg];              
          found = solveOrphanPoints(interpType0, interpMeshType, blkInterpData, 
                                    nit[blocNb], njt[blocNb], nkt[blocNb], 
                                    cellNatureField, x, y, z, 0, cfMax, cf, indi, ic, jc, kc);
          if (found != 0)
          {
            //interpType = indi[0];
            ind0 = indi[1] + indi[3] * ni1 + indi[5] * ni1nj1;
            ind0ext = ic-1+(jc-1)*(ni1+2)+(kc-1)*(ni1+2)*(nj1+2);

            E_Float* xni = coords[blocNb]->begin(1);
            E_Float* yni = coords[blocNb]->begin(2);
            E_Float* zni = coords[blocNb]->begin(3);
            K_METRIC::compVolOfStructCell3D(
              nit[blocNb], njt[blocNb], nkt[blocNb], ind0, -1,
              xni, yni, zni, vol);
            ret = extrapv[indg]; if (ret == 0) vol += 1000.;
            if (vol < volSavExtrap)
            {
              volSavExtrap = vol;
              //interpTypeSavExtrap = interpType;
              blocNbSavExtrap = blocNb;
              foundSav = found;
              cellSav = ind0; //cellule d'interpolation en centres
              cellExtSav = ind0ext; //cellule d'interpolation en centres etendus
              cfSav = cf;
            }
          }
        }// pour tous les blocs d'interpolation
        //E_Float check = 0.;
        //for (E_Int ii = 0; ii < 8; ii++) check += K_FUNC::E_abs(cfSav[ii]);
        //printf("%d: super extrap %d : %.15g\n", indg, foundSav, check);
        //for (E_Int ii = 0; ii < 8; ii++) printf("%g ", cfSav[ii]);
        //printf("\n");

        switch (foundSav)
        {
          case 1:
            cellsp[i] = cellSav;
            volsp[i] = volSavExtrap;
            interpTypesp[i] = 100; //interpTypeSavExtrap;
            donorCellsp[i] = cellExtSav;
            for (E_Int nf = 1; nf <= ncf; nf++)
              interpData->_coefs(i, nf) = cfSav[nf-1];
              
            indirection[blocNbSavExtrap*nbI+usedDomp[blocNbSavExtrap]] = i;
            usedDomp[blocNbSavExtrap]++;
            vectOfExtrapPts[blocNbSavExtrap].push_back(i);
            numberOfExtrapolatedPoints++;
            goto endPoint;
              
          case 2:
            cellsp[i] = cellSav;
            volsp[i] = volSavExtrap;
            interpTypesp[i] = 102;//interpTypeSavExtrap;
            donorCellsp[i] = cellExtSav;
            for (E_Int nf = 1; nf <= ncf; nf++)
              interpData->_coefs(i,nf) = cfSav[nf-1];
            indirection[blocNbSavExtrap*nbI+usedDomp[blocNbSavExtrap]] = i;
            usedDomp[blocNbSavExtrap]++;
            vectOfExtrapPts[blocNbSavExtrap].push_back(i);
            numberOfExtrapolatedPoints++;
            thetaPPoints.push_back(i);
            thetaPBlk.push_back(blocNbSavExtrap);
            goto endPoint;
              
          case 3:
            cellsp[i] = cellSav;
            volsp[i] = volSavExtrap;
            interpTypesp[i] = 103;//interpTypeSavExtrap;
            donorCellsp[i] = cellExtSav;
            for (E_Int nf = 1; nf <= ncf; nf++)
              interpData->_coefs(i,nf) = cfSav[nf-1];
            indirection[blocNbSavExtrap*nbI+usedDomp[blocNbSavExtrap]] = i;
            usedDomp[blocNbSavExtrap]++;
            vectOfExtrapPts[blocNbSavExtrap].push_back(i);
            numberOfExtrapolatedPoints++;
            thetaMPoints.push_back(i);
            thetaMBlk.push_back(blocNbSavExtrap);
            goto endPoint;
        }

        // Non interpolable, extrapolated
        volSavExtrap = K_CONST::E_MAX_FLOAT;
        for (E_Int blocNb = 0; blocNb < numberOfDomains; blocNb++)
        {
          K_KINTERP::BlkInterpData* blkInterpData = listOfInterpDatas[blocNb];
          FldArrayI& cellNatureField = *cellns[blocNb];

          x = coordx[i+nbI*blocNb];
          y = coordy[i+nbI*blocNb];
          z = coordz[i+nbI*blocNb];
              
          // Perform only explicit extrapolation: ic,jc,kc indice des centres etendus
          test = 1.;
          found = blkInterpData->getExtrapolationCell(x, y, z, ic, jc, kc, 
                                                      cf, cellNatureField, 0, test, interpType0, interpMeshType, cfMax);
          if (found != 0)
          {
            E_Float* xni = coords[blocNb]->begin(1);
            E_Float* yni = coords[blocNb]->begin(2);
            E_Float* zni = coords[blocNb]->begin(3);
            ni1 = nit[blocNb]-1; nj1 = njt[blocNb]-1; nk1 = nkt[blocNb]-1; ni1nj1 = ni1*nj1;
            ret = blkInterpData->fromExtendedToStandardCenters(ic, jc, kc, indi, interpType0); 
            //interpType = indi[0];
            ind0 = indi[1] + indi[3]*ni1 + indi[5]*ni1nj1;
            ind0ext = ic-1+(jc-1)*(ni1+2)+(kc-1)*(ni1+2)*(nj1+2);
          
            K_METRIC::compVolOfStructCell3D(
              nit[blocNb], njt[blocNb], nkt[blocNb], ind0, -1,
              xni, yni, zni, vol);
            if (ret == 0) vol = vol + 1000*test;

            if (vol < volSavExtrap)
            {
              volSavExtrap = vol;
              //interpTypeSavExtrap = interpType;
              blocNbSavExtrap = blocNb;
              foundSav = found;
              cellSav = ind0; // cellule d interpolation en centres
              cellExtSav = ind0ext; //cellule d interpolation en centres etendus
              cfSav = cf;
            }
          } //found
        }//pour tous les bloc d'interpolation

        //check = 0.;
        //for (E_Int ii = 0; ii < 8; ii++) check += K_FUNC::E_abs(cfSav[ii]);
        //printf("%d: last extrap %d : %g\n", indg, foundSav, check);
        switch (foundSav)
        {
          case 1: //cellule d'extrapolation du bloc initial 
            cellsp[i] = cellSav;
            volsp[i] = volSavExtrap;
            interpTypesp[i] = 100;//interpTypeSavExtrap;
            donorCellsp[i] = cellExtSav;
            for (E_Int nf = 1; nf <= ncf; nf++)
              interpData->_coefs(i, nf) = cfSav[nf-1]; 
            indirection[blocNbSavExtrap*nbI+usedDomp[blocNbSavExtrap]] = i;
            usedDomp[blocNbSavExtrap]++;
            vectOfExtrapPts[blocNbSavExtrap].push_back(i);
            numberOfExtrapolatedPoints++;
            goto endPoint;
              
          case 2: //cellule d'extrapolation du bloc en +theta
            cellsp[i] = cellSav;
            volsp[i] = volSavExtrap;
            interpTypesp[i] = 102;//interpTypeSavExtrap;
            donorCellsp[i] = cellExtSav;
            for (E_Int nf = 1; nf <= ncf; nf++)
              interpData->_coefs(i, nf) = cfSav[nf-1]; 
            
            thetaPPoints.push_back(i);
            thetaPBlk.push_back(blocNbSavExtrap);
            
            indirection[blocNbSavExtrap*nbI+usedDomp[blocNbSavExtrap]] = i;
            usedDomp[blocNbSavExtrap]++;
            vectOfExtrapPts[blocNbSavExtrap].push_back(i);
            numberOfExtrapolatedPoints++;
            goto endPoint;
            
          case 3:
            cellsp[i] = cellSav;
            volsp[i] = volSavExtrap;
            interpTypesp[i] = 103;//interpTypeSavExtrap;
            donorCellsp[i] = cellExtSav;
            for (E_Int nf = 1; nf <= ncf; nf++)
              interpData->_coefs(i, nf) = cfSav[nf-1]; 
              
            thetaMPoints.push_back(i);
            thetaMBlk.push_back(blocNbSavExtrap);
            
            indirection[blocNbSavExtrap*nbI+usedDomp[blocNbSavExtrap]] = i;
            usedDomp[blocNbSavExtrap]++;
            vectOfExtrapPts[blocNbSavExtrap].push_back(i);
            numberOfExtrapolatedPoints++;
            goto endPoint;
            
          default:
            break;
        }

        // BIG TROUBLE, cannot interpolate nor extrapolate
        orphanVect.push_back(i);
        numberOfPbPoints++;
    }
    endPoint:;
  } /* all points */
  
  // Recuperation du tableau correct d'indirection
  E_Int compt = 0;
  E_Int* indir1 = interpData->_indirectionPerDom.begin();
  for (E_Int blocNb = 0; blocNb < numberOfDomains; blocNb++)
  {   
    for (E_Int i = 0; i < usedDomp[blocNb]; i++)
    {
      indir1[compt] = indirection[blocNb*nbI+i];
      compt++;
    }
  }
  delete [] indirection;

  //recuperation des points interpol�s par un bloc periodique
  E_Int sizeP = thetaPPoints.size();
  E_Int sizeM = thetaMPoints.size();
  FldArrayI periodicPtsP(sizeP);
  FldArrayI periodicPtsM(sizeM);
  FldArrayI periodicBlksP(sizeP);
  FldArrayI periodicBlksM(sizeM);
   
  for (E_Int ii = 0; ii < sizeP; ii++)
  {
    periodicPtsP[ii]  = thetaPPoints[ii];
    periodicBlksP[ii] = thetaPBlk[ii];
  }
  for (E_Int ii = 0; ii < sizeM; ii++)
  {
    periodicPtsM[ii]  = thetaMPoints[ii];
    periodicBlksM[ii] = thetaMBlk[ii];
  }
  // // Summary of the treatments applied to interpolated points
  // if (zid == -1)
  //   printf("interpolated=%d, extrapolated=%d",numberOfInterpolatedPoints,numberOfExtrapolatedPoints);
  // else
  //   printf("Zone %d: interpolated=%d, extrapolated=%d",zid,numberOfInterpolatedPoints,numberOfExtrapolatedPoints);
  
  // // periodic treatment 
  // if (sizeP != 0)
  //   printf(" Periodic interpolated from +theta block=%d",sizeP);
  // if (sizeM != 0)
  //   printf(" Periodic interpolated from -theta block=%d",sizeM);

  // printf(", critical=%d\n",numberOfPbPoints);
  
  // if (numberOfPbPoints > 0)
  //   printf("# WARNING ! there are orphan points !\n");
    
  interpData->_sizeOfEachIndirectionTab.malloc(numberOfDomains);  
  for (E_Int i = 0; i < numberOfDomains; i++)
    interpData->_sizeOfEachIndirectionTab[i] = usedDomp[i];
}

//=============================================================================
// IN: x,y,z: point coordinates we try to extrapolate
// IN: block: domain tested for extrapolation
// IN: testNature: 0: explicite, 1: pas de cellules masquees
// OUT: indExtrap: cellule d extrapolation en centres
//      cf:        extrapolation coefficients
//=============================================================================
short K_CONNECTOR::solveOrphanPoints(
  K_KINTERP::BlkInterpData::InterpolationType interpType0,
  K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
  K_KINTERP::BlkInterpData* blkInterpData,
  E_Int ni, E_Int nj, E_Int nk, FldArrayI& cellN,
  E_Float x, E_Float y, E_Float z, 
  E_Int testNature, E_Float cfMax,
  FldArrayF& cf, FldArrayI& indExtrap, E_Int& ic, E_Int& jc, E_Int& kc)
{
  E_Int nic = ni+1; E_Int njc = nj+1; E_Int nkc = nk+1;//centres etendus
  E_Int ip, jp, kp;
  E_Int is, js, ks;
  short foundInDomain;
  //E_Int ncf=8, nindi=7;  // nfld du tableau indit contenant les indices des cellules d'interp + ordre

  short found = blkInterpData->getInterpolationCellStruct(
    x, y, z, ic, jc, kc, cf, 
    interpMeshType, interpType0);

  if (found == 0) return 0;

  // Init cf
  cf.setAllValuesAtNull();
  
  // Sauvegarde de la cellule initiale (centre non ext)
  E_Int icsav = ic; E_Int jcsav = jc; E_Int kcsav = kc;

  // ===========
  // Ordre 0
  // ===========
  foundInDomain = blkInterpData->getExtrapolationCoeffForCell(
    x, y, z, ic, jc, kc, cf,
    cellN, testNature, 0, is, js, ks, cfMax,
    interpType0, interpMeshType);
  if (is < 1 || is > nic-1 || js < 1 || js > njc-1 || ks < 1 || ks > nkc-1)
  { is = ic; js = jc; ks = kc; }
  
  if (foundInDomain != 0) goto fin;
  
  ic = is; jc = js; kc = ks;

  foundInDomain = blkInterpData->getExtrapolationCoeffForCell(
    x, y, z, ic, jc, kc, cf,
    cellN, testNature, 0, is, js, ks, cfMax,
    interpType0, interpMeshType);
  if (is < 1 || is > nic-1 || js < 1 || js > njc-1 || ks < 1 || ks > nkc-1)
  { is = ic; js = jc; ks = kc; }

  if (foundInDomain != 0) goto fin;

  // Try neighbouring point of previous one
  if ((jcsav != js) || (kcsav != ks))
  {
    if (is-1 >= 1)
    {
      ic = is-1;
      jc = js;
      kc = ks;
      foundInDomain = 
        blkInterpData->getExtrapolationCoeffForCell(x, y, z, ic, jc, kc, cf,
                                                    cellN, testNature,  0, is, js, ks, cfMax,
                                                    interpType0, interpMeshType);
      if (is < 1 || is > nic-1 || js < 1 || js > njc-1 || ks < 1 || ks > nkc-1)
      { is = ic; js = jc; ks = kc; }

      if (foundInDomain != 0) goto fin;
    }
    else if (is+1 <= nic-1)
    {
      ic = is+1;
      jc = js;
      kc = ks;
      foundInDomain =
        blkInterpData->getExtrapolationCoeffForCell(x, y, z,ic, jc, kc, cf,
                                                    cellN, testNature,  0, is, js, ks, cfMax,
                                                    interpType0, interpMeshType);
      if (is < 1 || is > nic-1 || js < 1 || js > njc-1 || ks < 1 || ks > nkc-1)
      { is = ic; js = jc; ks = kc; }
      if (foundInDomain != 0) goto fin;
    }
  }
  else
  {
    if (js-1 >= 1)
    {
      ic = is;
      jc = js-1;
      kc = ks;
      foundInDomain =
        blkInterpData->getExtrapolationCoeffForCell(x, y, z, ic, jc, kc, cf,
                                                    cellN, testNature, 0, is, js, ks, cfMax,
                                                    interpType0, interpMeshType);
      if (is < 1 || is > nic-1 || js < 1 || js > njc-1 || ks < 1 || ks > nkc-1)
      { is = ic; js = jc; ks = kc; }

      if (foundInDomain != 0) goto fin;
    }
    else if (js+1 <= njc-1)
    {
      ic = is;
      jc = js+1;
      kc = ks;
      foundInDomain =
        blkInterpData->getExtrapolationCoeffForCell(x, y, z, ic, jc, kc, cf,
                                                    cellN, testNature,  0, is, js, ks, cfMax,
                                                    interpType0, interpMeshType);
      if (is < 1 || is > nic-1 || js < 1 || js > njc-1 || ks < 1 || ks > nkc-1)
      { is = ic; js = jc; ks = kc; }
      
      if (foundInDomain != 0) goto fin;
    }
  }
  
  // Final untrapment - we have no choice but test all neighbours
  for (ip = -1; ip <= 1; ip++)
    for (jp = -1; jp <= 1; jp++)
      for (kp = -1; kp <= 1; kp++)
      {
        if (ip != 0 && jp != 0 && kp != 0)
        {
          ic = icsav+ip;
          jc = jcsav+jp;
          kc = kcsav+kp;
          
          if (ic >= 1 && ic <= nic-1 && jc >= 1 && jc <= njc-1 &&
              kc >= 1 && kc <= nkc-1)
          {
            foundInDomain =
              blkInterpData->getExtrapolationCoeffForCell(x, y, z, ic, jc, kc, cf,
                                                          cellN, testNature, 0, is, js, ks, cfMax,
                                                          interpType0, interpMeshType);
      
            if (is < 1 || is > nic-1 || js < 1 || js > njc-1 || ks < 1 || ks > nkc-1)
            { is = ic; js = jc; ks = kc; }
  
            if (foundInDomain != 0) goto fin;
          }
        }
      }

  for (ip = -2; ip <= 2; ip++)
    for (jp = -2; jp <= 2; jp++)
      for (kp = -2; kp <= 2; kp++)
      {
        if ((ip > 1 || ip < -1) && (jp > 1 || jp < -1) && (kp > 1 || kp < -1))
        {
          ic = icsav+ip;
          jc = jcsav+jp;
          kc = kcsav+kp;
          
          if (ic >= 1 && ic <= nic-1 && jc >= 1 && jc <= njc-1 &&
              kc >= 1 && kc <= nkc-1)
          {
            foundInDomain =
              blkInterpData->getExtrapolationCoeffForCell(x, y, z, ic, jc, kc, cf,
                                                          cellN, testNature, 0, is, js, ks, cfMax,
                                                          interpType0, interpMeshType);
  
            if (is < 1 || is > nic-1 || js < 1 || js > njc-1 || ks < 1 || ks > nkc-1)
            { is = ic; js = jc; ks = kc; }
                
            if (foundInDomain != 0) goto fin;
          }
        }
      }

  ic = icsav; jc = jcsav; kc = kcsav;
  
  fin:
  blkInterpData->fromExtendedToStandardCenters(ic, jc, kc, indExtrap, interpType0);
  return foundInDomain;
}
