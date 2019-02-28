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
# include <stdio.h>
# include "connector.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Calcule et stocke les coefficients d'interpolation 
   CAS SANS DOUBLE WALL 
   IN: receiverArray: points a interpoler définis sous forme de maillage
   IN: donorArrays: maillages donneurs. La localisation du donneur 
        (noeuds/centres/centres étendus) doit être effectuee au prealable
   IN: Order: ordre des interpolations (2, 3, 5)
   IN: Nature: 0: produit des cellN=0 -> donneur invalide; 
               1: cellN=0 ou 2 -> donneur invalide
   IN: PenaliseBorders: 1: penalite sur le volume des pts ou cellules frontieres
   IN: allHooks != Py_None: un hook par adt associé à un donor 
   OUT: [donorBlks,donorInd1D, donorType, coefs, extrapInd1D, orphanInd1D] 
        donorBlks: no du blk donneur, démarre à 0
        donorInd1D: indice global (structure), de l elt (NS) du donneur
        donorType: type d interpolation effectué localement
        coefs: coefficients d interpolation, stockés selon le type
        extrapInd1D: indices des pts extrapolés
        orphanInd1D: indices des pts orphelins
*/
//=============================================================================
PyObject* K_CONNECTOR::setInterpData(PyObject* self, PyObject* args)
{
  PyObject *receiverArray; 
  PyObject *donorArrays;// domaines d interpolation
  E_Int Order;
  E_Int Nature;
  E_Int PenalizeBorders;
  E_Int InterpDataType;// 1 : ADT ou 0 : CART
  PyObject* allHooks;
  if (!PYPARSETUPLEI(args, "OOllllO", "OOiiiiO",
                    &receiverArray, &donorArrays, &Order, &Nature, &PenalizeBorders, 
                    &InterpDataType, &allHooks))
  {
      return NULL;
  }
  if ( InterpDataType != 0 && InterpDataType != 1 )
  {
    PyErr_SetString(PyExc_TypeError,
                    "setInterpData: InterpDataType must be 0 for CART or 1 for ADT.");
    return NULL;
  }
  // ordre des interpolations
  E_Int interporder = Order;
  E_Int nature = Nature; // O: produit des cellN=0 -> donneur invalide; 1: cellN=0 ou 2 -> donneur invalide
  E_Int penalty = PenalizeBorders;//1 : penalité sur le volume des pts ou cellules frontieres
  // Interpolation type
  K_INTERP::InterpData::InterpolationType interpType;
  E_Int nindi, ncfmax;
  switch (interporder)
  {
    case 2:
      interpType = K_INTERP::InterpData::O2CF;
      ncfmax = 8; nindi = 1;
      break;

    case 3:
      interpType = K_INTERP::InterpData::O3ABC;
      ncfmax = 9; nindi = 1;
      break;

    case 5:
      interpType = K_INTERP::InterpData::O5ABC;
      ncfmax = 15; nindi = 1;
      break;
        
    default:
      printf("Warning: setInterpData: unknown interpolation order.");
      printf(" Set to 2nd order.\n");
      interpType = K_INTERP::InterpData::O2CF;
      ncfmax = 8; nindi = 1;
      break;
  } 
  FldArrayI indi(nindi); FldArrayF cf(ncfmax);

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
                    "setInterpData: 1st arg must contain coordinates.");
    RELEASESHAREDB(resr, receiverArray, fr, cnr); 
    return NULL;
  }
  posxr++; posyr++; poszr++;
  E_Int posdir = K_ARRAY::isNamePresent("EXdir",varStringr);
  E_Int isEX = 0;
  if (posdir > -1) {isEX=1; posdir++;}
  /*-------------------------------------------------------*/
  /* Extraction des infos sur les domaines d interpolation */
  /*-------------------------------------------------------*/
  vector<E_Int> resl;  vector<char*> varString;
  vector<FldArrayF*> fields;
  vector<void*> a2; //ni,nj,nk ou cnt en NS
  vector<void*> a3; //eltType en NS
  vector<void*> a4;
  vector<PyObject*> objs;
  E_Boolean skipNoCoord = true;  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false;  E_Boolean skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(donorArrays, resl, varString, fields, a2, a3, a4, objs,  
                                      skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nzones = objs.size();

  if (isOk == -1)
  {
    RELEASESHAREDB(resr, receiverArray, fr, cnr); 
    for (E_Int no = 0; no < nzones; no++)
      RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
    PyErr_SetString(PyExc_TypeError,
                    "setInterpData: 2nd argument is not valid.");
    return NULL;
  }   
  if (nzones == 0) 
  {
    RELEASESHAREDB(resr, receiverArray, fr, cnr); 
    for (E_Int no = 0; no < nzones; no++)
      RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
    PyErr_SetString(PyExc_TypeError,
                    "setInterpData: no valid donor zone found.");
    return NULL;
  }
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    if (resl[noz] == 2) 
    {
      char* eltType0 = (char*)a3[noz];
      if (K_STRING::cmp(eltType0,"TETRA")!= 0)
      {
        RELEASESHAREDB(resr, receiverArray, fr, cnr); 
        for (E_Int no = 0; no < nzones; no++)
          RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
        PyErr_SetString(PyExc_TypeError,
                        "setInterpData: unstructured donor zones must be TETRA.");
        return NULL;
      }
    }
    else if (resl[noz] == 1) 
    {
      if (*(E_Int*)a2[noz]<2 || *(E_Int*)a3[noz]<2 ) 
      {
        RELEASESHAREDB(resr, receiverArray, fr, cnr); 
        for (E_Int no = 0; no < nzones; no++)
          RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
        PyErr_SetString(PyExc_TypeError,
                        "setInterpData: structured donor zones must be 3D or nk=1.");
        return NULL;
      }
    }
  }
  E_Int nzonesS = 0; E_Int nzonesU = 0;
  vector<E_Int> posxs; vector<E_Int> posys; vector<E_Int> poszs; vector<E_Int> poscs;
  vector<void*> a5;
  for (E_Int no = 0; no < nzones; no++)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString[no]); posx++;
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString[no]); posy++;
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString[no]); posz++;
    E_Int posc = K_ARRAY::isCellNatureField2Present(varString[no]); posc++;
    if (a4[no] == NULL) nzonesU++;
    else nzonesS++;
    posxs.push_back(posx); 
    posys.push_back(posy); 
    poszs.push_back(posz); 
    poscs.push_back(posc); 
    a5.push_back(NULL);// PAS DE CONNECTIVITE ELTS/ELTS VOISINS
  }
  // Liste des interpDatas
  vector<E_Int> nis; vector<E_Int> njs; vector<E_Int> nks;
  vector<FldArrayI*> cnt;
  vector<K_INTERP::InterpData*> interpDatas;
  E_Int isBuilt;
  // creation des interpDatas
  if (allHooks == Py_None)
  {
    E_Int failed = 0;
    for (E_Int no = 0; no < nzones; no++)
    {
      if ( InterpDataType == 1)
      {
        K_INTERP::InterpAdt* adt = new K_INTERP::InterpAdt(fields[no]->getSize(), 
                                                           fields[no]->begin(posxs[no]),
                                                           fields[no]->begin(posys[no]),
                                                           fields[no]->begin(poszs[no]),
                                                           a2[no], a3[no], a4[no], isBuilt);
        if ( isBuilt == 1 ) interpDatas.push_back(adt);
        else failed = 1;
      }
      else //CART
      {
        if ( resl[no] == 1)
        {
          E_Float* xt = fields[no]->begin(posxs[no]);
          E_Float* yt = fields[no]->begin(posys[no]);
          E_Float* zt = fields[no]->begin(poszs[no]);
          E_Int ni = *(E_Int*)a2[no];
          E_Int nj = *(E_Int*)a3[no]; 
          E_Int nk = *(E_Int*)a4[no];
          E_Float x0 = xt[0]; E_Float y0 = yt[0]; E_Float z0 = zt[0];
          E_Float hi = xt[1]-xt[0];
          E_Float hj = yt[ni]-yt[0];
          E_Float hk = zt[ni*nj]-zt[0];
          K_INTERP::InterpCart* interpCart = new K_INTERP::InterpCart(ni,nj,nk,hi,hj,hk,x0,y0,z0);
          interpDatas.push_back(interpCart);
        }
      }
      if (failed==1)
      {
        RELEASESHAREDB(resr, receiverArray, fr, cnr); 
        for (E_Int no = 0; no < nzones; no++)
          RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
        for (size_t noi = 0; noi < interpDatas.size(); noi++)
          delete interpDatas[noi];
        PyErr_SetString(PyExc_TypeError,
                        "setInterpData: structured donor zones must be 3D or 2D with z=constant.");
        return NULL;
      }
    }
  }
  else //hook fourni pour chq donneur
  {
    if (PyList_Check(allHooks) == false)
    {
      RELEASESHAREDB(resr, receiverArray, fr, cnr); 
      for (E_Int no = 0; no < nzones; no++)
        RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);  
      PyErr_SetString(PyExc_TypeError, 
                      "setInterpData: hook must be a list of hooks.");
      return NULL;
    }
    E_Int nHooks = PyList_Size(allHooks);
    E_Int nDonorZones = fields.size();
    if (nHooks != nDonorZones)
    {
      RELEASESHAREDB(resr, receiverArray, fr, cnr); 
      for (E_Int no = 0; no < nzones; no++)
        RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
      PyErr_SetString(PyExc_TypeError,
                      "setInterpData: size of list of hooks must be equal to the number of donor zones.");
      return NULL;
    }
    for (E_Int noh = 0; noh < nHooks; noh++)
    {
      PyObject* hook = PyList_GetItem(allHooks,noh);
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
      void** packet = (void**) PyCObject_AsVoidPtr(hook);
#else
      void** packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
      E_Int* type = (E_Int*)packet[0];
      if (type[0] != 1) 
      {
        RELEASESHAREDB(resr, receiverArray, fr, cnr); 
        for (E_Int no = 0; no < nzones; no++)
          RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
        PyErr_SetString(PyExc_TypeError,
                        "setInterpData: hook must define an ADT.");
        return NULL;
      }
      if (type[1] != 1) 
      {
        RELEASESHAREDB(resr, receiverArray, fr, cnr); 
        for (E_Int no = 0; no < nzones; no++)
          RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
        PyErr_SetString(PyExc_TypeError,
                        "setInterpData: one ADT must be defined per hook.");
        return NULL;
      }
      interpDatas.push_back((K_INTERP::InterpAdt*)(packet[1])); 
    }
  }
   // cas seulement non structure : ncf a 4 (minimum)
  if (nzonesU != 0)
  {
    if (interporder != 2) printf("Warning: setInterpData: interpolation order is 2 for tetra arrays.\n");
    if (nzonesS == 0) ncfmax = 4;
  }
  /*-------------------------------------------------------*/
  /* Calcul des coefficients d'interpolation               */
  /*-------------------------------------------------------*/
  E_Float* xr = fr->begin(posxr);
  E_Float* yr = fr->begin(posyr);
  E_Float* zr = fr->begin(poszr);
  E_Int nbI = fr->getSize();
  E_Float vol;

  // Donnees d'interpolation
  vector<FldArrayF*> listOfInterpCoefs; // un par thread
  vector<FldArrayI*> listOfDonorInd1D;
  vector<FldArrayI*> listOfRcvInd1D;
  vector<FldArrayI*> listOfDonorTypes;
  vector<FldArrayI*> listOfExtrapInd1D;
  vector<FldArrayF*> listOfEXdirs;
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    FldArrayF* donorInterpCoefs = new FldArrayF(nbI*ncfmax); donorInterpCoefs->setAllValuesAtNull();
    listOfInterpCoefs.push_back(donorInterpCoefs);
    FldArrayI* donorInd1D = new FldArrayI(nbI*2); donorInd1D->setAllValuesAt(-1);
    listOfDonorInd1D.push_back(donorInd1D);
    FldArrayI* rcvInd1D = new FldArrayI(nbI); rcvInd1D->setAllValuesAt(-1);
    listOfRcvInd1D.push_back(rcvInd1D);
    FldArrayI* donorType = new FldArrayI(nbI); donorType->setAllValuesAtNull();
    listOfDonorTypes.push_back(donorType);
    FldArrayI* extrapPoints = new FldArrayI(nbI); extrapPoints->setAllValuesAt(-1);
    listOfExtrapInd1D.push_back(extrapPoints);  
    if (isEX == 1) 
    {
      FldArrayF* EXdirs = new FldArrayF(nbI); EXdirs->setAllValuesAtNull();
      listOfEXdirs.push_back(EXdirs);
    }
  }

  FldArrayI* orphanPts = new FldArrayI(nbI); orphanPts->setAllValuesAt(-1); // un par thread
  E_Int* orphanPts1D = orphanPts->begin();
  E_Int type, isExtrapolated;
  FldArrayI sizeOfIndDonor1D(nzones); sizeOfIndDonor1D.setAllValuesAtNull(); // un par thread
  FldArrayI usedDomp(nzones); usedDomp.setAllValuesAtNull();//nb de pts interpoles par domaine donneur
  FldArrayI sizeCoefs(nzones); sizeCoefs.setAllValuesAtNull();// taille du tableau de coefficients par domaine donneur
  FldArrayI nbExtraps(nzones); nbExtraps.setAllValuesAtNull();
  E_Int nExtrapPts = 0; // par thread
  E_Int noOrphanPt = 0;

  E_Float* EXdir0 = NULL;
  if (isEX == 1) EXdir0 = fr->begin(posdir);
  
  E_Float x, y, z; E_Int noblk = 0;

  for (E_Int ind = 0; ind < nbI; ind++)
  {
    x = xr[ind]; y = yr[ind]; z = zr[ind];
    short ok = K_INTERP::getInterpolationCell(
      x, y, z, interpDatas, fields,
      a2, a3, a4, a5, posxs, posys, poszs, poscs,
      vol, indi, cf, type, noblk, interpType, nature, penalty);   
    isExtrapolated = 0;
    if (ok != 1)
    {
      ok = K_INTERP::getExtrapolationCell(
        x, y, z, interpDatas, fields,
        a2, a3, a4, a5, posxs, posys, poszs, poscs,
        vol, indi, cf, type, noblk, interpType, nature, penalty); 
      if (noblk > 0) isExtrapolated = 1;
    }
    if (noblk > 0)
    {
      E_Int noDonorBlk = noblk-1;
      E_Int& noInterpPtForBlk = usedDomp[noDonorBlk];// on met une reference car on va modifier cette valeur
      E_Int& sizeOfIndDonor1DForBlk = sizeOfIndDonor1D[noDonorBlk];//idem
      E_Int& sizecf = sizeCoefs[noDonorBlk];//idem
      E_Int& noExtrapPtForBlk = nbExtraps[noDonorBlk]; 
      E_Float* donorCf = listOfInterpCoefs[noDonorBlk]->begin();
      E_Int* donorType = listOfDonorTypes[noDonorBlk]->begin();
      E_Int* donorInd1D = listOfDonorInd1D[noDonorBlk]->begin();
      E_Int* rcvInd1D = listOfRcvInd1D[noDonorBlk]->begin();
      E_Int* extrapInd1D = listOfExtrapInd1D[noDonorBlk]->begin();
      E_Float* EXdir = NULL;
      if (isEX == 1)  
      {
        EXdir = listOfEXdirs[noDonorBlk]->begin();
        EXdir[noInterpPtForBlk] = EXdir0[ind];
      }
      donorType[noInterpPtForBlk] = type;         
      rcvInd1D[noInterpPtForBlk] = ind;
      switch (type)
      {
        case 1:
          donorCf[sizecf] = cf[0];
          sizecf += 1;
          donorInd1D[sizeOfIndDonor1DForBlk] = indi[0];
          sizeOfIndDonor1DForBlk++; 
          break;

        case 2:          
          for (E_Int nocf = 0; nocf < 8; nocf++)
            donorCf[sizecf+nocf] = cf[nocf];
          sizecf += 8;
          donorInd1D[sizeOfIndDonor1DForBlk] = indi[0];
          sizeOfIndDonor1DForBlk++; 
          break;

        case 22:
          for (E_Int nocf = 0; nocf < 4; nocf++)
            donorCf[sizecf+nocf] = cf[nocf];
          sizecf += 4;
          donorInd1D[sizeOfIndDonor1DForBlk] = indi[0];
          sizeOfIndDonor1DForBlk++; 
          break;

        case 3: 
          for (E_Int nocf = 0; nocf < 9; nocf++) 
            donorCf[sizecf+nocf] = cf[nocf];
          sizecf += 9;
          donorInd1D[sizeOfIndDonor1DForBlk] = indi[0]; 
          sizeOfIndDonor1DForBlk++;
          break;

        case 4:
          if (a4[noDonorBlk] == NULL ) //non structure
          {
            for (E_Int nocf = 0; nocf < 4; nocf++) 
              donorCf[sizecf+nocf] = cf[nocf];
            sizecf += 4;
            donorInd1D[sizeOfIndDonor1DForBlk] = indi[0];
            sizeOfIndDonor1DForBlk++;
          }
          else // structure
          {
            for (E_Int nocf = 0; nocf < 8; nocf++) 
              donorCf[sizecf+nocf] = cf[nocf];
            sizecf += 8;
            donorInd1D[sizeOfIndDonor1DForBlk] = indi[0];
            sizeOfIndDonor1DForBlk++;
          } 
          break;

        case 5:
          for (E_Int nocf = 0; nocf < 15; nocf++) 
            donorCf[sizecf+nocf] = cf[nocf];
          sizecf += 15;
          donorInd1D[sizeOfIndDonor1DForBlk] = indi[0]; 
          sizeOfIndDonor1DForBlk++;
          break;

        default:
          printf("Error: setInterpData: type not yet implemented.\n");
          exit(0);
      }      
      if (isExtrapolated == 1) 
      { 
        extrapInd1D[noExtrapPtForBlk] = ind; 
        noExtrapPtForBlk++;// = nbExtraps[noDonorBlk]++;
        nExtrapPts++;
      }
      noInterpPtForBlk++;
    }
    else {orphanPts1D[noOrphanPt] = ind; noOrphanPt++;}
  }

  // Redimensionnement des tableaux
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    E_Int sizecf = sizeCoefs[noz];
    E_Int nbInterp = usedDomp[noz];
    E_Int nbExtrap = nbExtraps[noz];   
    E_Int sizeOfIndDonor1DL = sizeOfIndDonor1D[noz];
    listOfInterpCoefs[noz]->resize(sizecf);
    listOfDonorInd1D[noz]->resize(sizeOfIndDonor1DL);
    listOfRcvInd1D[noz]->resize(nbInterp);
    listOfDonorTypes[noz]->resize(nbInterp);
    listOfExtrapInd1D[noz]->resize(nbExtrap);
    if (isEX == 1) listOfEXdirs[noz]->resize(nbInterp);
  }
  // Nettoyages
  for (E_Int no = 0; no < nzones; no++)
  {
    if (allHooks == Py_None) delete interpDatas[no];
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
  PyObject* PyListCellIndicesR = PyList_New(0);
  PyObject* PyListCellIndicesD = PyList_New(0);
  // liste pour stocker les coefficients d'interpolation par bloc de domaine d'interpolation correspondant
  PyObject* PyListCoefficients = PyList_New(0);
  // liste pour stocker les types d interpolation par domaine d'interpolation correspondant
  PyObject* PyListInterpTypes = PyList_New(0);
  // listes pour marquer les indices des cellules extrapoles
  PyObject* PyListCellIndicesE = PyList_New(0);
  // listes pour marquer les EXdir (directions pour les pts EX)
  PyObject* PyListEXDir = PyList_New(0);

  // ------------------------------
  // Construction des PyArrays
  // ------------------------------
  for (E_Int noz = 0; noz < nzones; noz++)
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

    if (isEX == 1)
    {
      PyObject* pyEXDir = K_NUMPY::buildNumpyArray(*listOfEXdirs[noz],1);
      PyList_Append(PyListEXDir, pyEXDir); Py_DECREF(pyEXDir);
      delete listOfEXdirs[noz];
    }

  } // fin parcours des zones donneuses
  
  //   indices des points orphelins
  PyObject* cellIndO = K_NUMPY::buildNumpyArray(*orphanPts,1);
  delete orphanPts;

  PyObject* tpl = Py_BuildValue("[OOOOOOO]", PyListCellIndicesR, PyListCellIndicesD, PyListInterpTypes, PyListCoefficients,
                                PyListCellIndicesE, cellIndO, PyListEXDir);
  Py_DECREF(cellIndO); Py_DECREF(PyListInterpTypes); 
  Py_DECREF(PyListCellIndicesE); Py_DECREF(PyListEXDir); 
  Py_DECREF(PyListCoefficients); Py_DECREF(PyListCellIndicesR); 
  Py_DECREF(PyListCellIndicesD); 
  return tpl;
}

//=============================================================================
/* Calcule et stocke les coefficients d'interpolation 
   CAS AVEC DOUBLE WALL 
   IN: receiverArrays: points à interpoler définis sous forme de maillage     
                       chaque array correspond aux mêmes points mais modifiés par changeWall
                       en fonction des donneurs
   IN: donorArrays: maillages donneurs. La localisation du donneur 
       (noeuds/centres/centres étendus) doit être effectuée au préalable
   !!  receiverArrays et donorArrays doivent être ordonnés de la même manière
   IN: Order: ordre des interpolations (2, 3, 5)
   IN: Nature: 0: produit des cellN=0 -> donneur invalide; 
               1: cellN=0 ou 2 -> donneur invalide
   IN: PenaliseBorders: 1: penalité sur le volume des pts ou cellules frontières
   IN: hook != Py_None: hook sur les adt associés à donorArrays 
   OUT: [donorBlks,donorInd1D, donorType, coefs, extrapInd1D, orphanInd1D] 
        donorBlks: no du blk donneur, démarre à 0
        donorInd1D: indice global (structure), de l elt (NS) du donneur
        donorType: type d interpolation effectué localement
        coefs: coefficients d interpolation, stockés selon le type
        extrapInd1D: indices des pts extrapolés
        orphanInd1D: indices des pts orphelins */
//=============================================================================
PyObject* K_CONNECTOR::setInterpDataDW(PyObject* self, PyObject* args)
{
  PyObject *receiverArrays;
  PyObject *donorArrays;// domaines d'interpolation
  E_Int Order;
  E_Int Nature;
  E_Int PenalizeBorders;
  PyObject* hook;

  if (!PYPARSETUPLEI(args,
                    "OOlllO", "OOiiiO",
                    &receiverArrays, &donorArrays, &Order, &Nature, &PenalizeBorders, &hook))
  {
      return NULL;
  }
  
  // ordre des interpolations
  E_Int interporder = Order;
  E_Int nature = Nature; // O: produit des cellN=0 -> donneur invalide; 1: cellN=0 ou 2 -> donneur invalide
  E_Int penalty = PenalizeBorders;//1 : penalité sur le volume des pts ou cellules frontieres
  // Interpolation type
  K_INTERP::InterpData::InterpolationType interpType;
  E_Int nindi, ncfmax;
  switch (interporder)
  {
    case 2:
      interpType = K_INTERP::InterpData::O2CF;
      ncfmax = 8; nindi = 1;
      break;

    case 3:
      interpType = K_INTERP::InterpData::O3ABC;
      ncfmax = 9; nindi = 1;
      break;

    case 5:
      interpType = K_INTERP::InterpData::O5ABC;
      ncfmax = 15; nindi = 1;
      break;

    default:
      printf("Warning: setInterpDataDW: unknown interpolation order.");
      printf(" Set to 2nd order.\n");
      interpType = K_INTERP::InterpData::O2CF;
      ncfmax = 8; nindi = 1;
      break;
  } 
  FldArrayI indi(nindi); FldArrayF cf(ncfmax);

  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  vector<E_Int> resr;  vector<char*> varStringr;
  vector<FldArrayF*> fieldrs;
  vector<void*> ar2; //ni,nj,nk ou cnt en NS
  vector<void*> ar3; //eltType en NS
  vector<void*> ar4;
  vector<PyObject*> objrs;
  E_Boolean skipNoCoord = true;  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false;  E_Boolean skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    receiverArrays, resr, varStringr, fieldrs, ar2, ar3, ar4, objrs,  
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nzonesr = objrs.size();
  vector<E_Int> posxr; vector<E_Int> posyr; vector<E_Int> poszr; vector<E_Int> posdir;
  E_Int isEX = 0;
  for (E_Int no = 0; no < nzonesr; no++)
  {
    E_Int posxr0 = K_ARRAY::isCoordinateXPresent(varStringr[no]);
    E_Int posyr0 = K_ARRAY::isCoordinateYPresent(varStringr[no]);
    E_Int poszr0 = K_ARRAY::isCoordinateZPresent(varStringr[no]);
    E_Int posdir0 = K_ARRAY::isNamePresent("EXdir",varStringr[no]);
    posxr0++; posyr0++; poszr0++; posdir0++;
    if (posxr0 == -1 || posyr0 == -1 || poszr0 == -1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "setInterpDataDW: 1st arg must contain coordinates.");
      for (E_Int no2 = 0; no2 < nzonesr; no2++)
        RELEASESHAREDA(resr[no2],objrs[no2],fieldrs[no2],ar2[no2],ar3[no2],ar4[no2]);      
      return NULL;
    }
    if (posdir0 != 0) isEX = 1;
    posxr.push_back(posxr0); posyr.push_back(posyr0); poszr.push_back(poszr0); posdir.push_back(posdir0);
  }
 
  /*-------------------------------------------------------*/
  /* Extraction des infos sur les domaines d interpolation */
  /*-------------------------------------------------------*/
  vector<E_Int> resl;  vector<char*> varString;
  vector<FldArrayF*> fields;
  vector<void*> a2; //ni,nj,nk ou cnt en NS
  vector<void*> a3; //eltType en NS
  vector<void*> a4;
  vector<PyObject*> objs;
  isOk = K_ARRAY::getFromArrays(
    donorArrays, resl, varString, fields, a2, a3, a4, objs,  
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nzones = objs.size();
  if (nzonesr != nzones)
  {
    for (E_Int no = 0; no < nzonesr; no++)
      RELEASESHAREDA(resr[no],objrs[no],fieldrs[no],ar2[no],ar3[no],ar4[no]);       
    for (E_Int no = 0; no < nzones; no++)
      RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
    PyErr_SetString(PyExc_TypeError,
                    "setInterpDataDW: 1st and 2nd arguments must be of same length.");
    return NULL;
  }
  if (isOk == -1)
  {
    for (E_Int no = 0; no < nzonesr; no++)
      RELEASESHAREDA(resr[no],objrs[no],fieldrs[no],ar2[no],ar3[no],ar4[no]);
    for (E_Int no = 0; no < nzones; no++)
      RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
    PyErr_SetString(PyExc_TypeError,
                    "setInterpDataDW: 2nd argument is not valid.");
    return NULL;
  }   
  if (nzones == 0) 
  {
    for (E_Int no = 0; no < nzonesr; no++)
      RELEASESHAREDA(resr[no],objrs[no],fieldrs[no],ar2[no],ar3[no],ar4[no]);
    PyErr_SetString(PyExc_TypeError,
                    "setInterpDataDW: no valid donor zone found.");
    return NULL;
  }
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    if (resl[noz] == 2)
    {
      char* eltType0 = (char*)a3[noz];
      if (K_STRING::cmp(eltType0, "TETRA") != 0)
      {
        for (E_Int no = 0; no < nzonesr; no++)
          RELEASESHAREDA(resr[no],objrs[no],fieldrs[no],ar2[no],ar3[no],ar4[no]);
        for (E_Int no = 0; no < nzones; no++)
          RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
        PyErr_SetString(PyExc_TypeError,
                        "setInterpDataDW: unstructured donor arrays must be TETRA.");
        return NULL;
      }
    }
    else if (resl[noz] == 2) 
    {
      if (*(E_Int*)a2[noz]<2 || *(E_Int*)a3[noz]<2 )//|| *(E_Int*)a4[noz]<2) 
      {
        for (E_Int no = 0; no < nzonesr; no++)
          RELEASESHAREDA(resr[no],objrs[no],fieldrs[no],ar2[no],ar3[no],ar4[no]);  
        for (E_Int no = 0; no < nzones; no++)
          RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
        PyErr_SetString(PyExc_TypeError,
                        "setInterpDataDW: structured donor zones must be 3D or 2D with nk=1.");
        return NULL;
      }
    }
  }
  E_Int nzonesS = 0; E_Int nzonesU = 0;
  vector<E_Int> posxs; vector<E_Int> posys; vector<E_Int> poszs; vector<E_Int> poscs;
  vector<void*> a5;
  for (E_Int no = 0; no < nzones; no++)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString[no]); posx++;
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString[no]); posy++;
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString[no]); posz++;
    E_Int posc = K_ARRAY::isCellNatureField2Present(varString[no]); posc++;
    if ( a4[no] == NULL) nzonesU++;
    else nzonesS++;
    posxs.push_back(posx); 
    posys.push_back(posy); 
    poszs.push_back(posz); 
    poscs.push_back(posc); 
    a5.push_back(NULL);// PAS DE CONNECTIVITE ELTS/ELTS VOISINS
  }
  // Liste des interpDatas
  vector<E_Int> nis; vector<E_Int> njs; vector<E_Int> nks;
  vector<FldArrayI*> cnt;
  vector<K_INTERP::InterpData*> interpDatas;
  E_Int isBuilt;
  // creation des interpDatas
  if (hook == Py_None)
  {
    for (E_Int no = 0; no < nzones; no++)
    {
      K_INTERP::InterpAdt* adt = new K_INTERP::InterpAdt(
        fields[no]->getSize(), 
        fields[no]->begin(posxs[no]),
        fields[no]->begin(posys[no]),
        fields[no]->begin(poszs[no]),
        a2[no], a3[no], a4[no], isBuilt);
      if ( isBuilt == 1 ) 
        interpDatas.push_back(adt);
      else 
      {
        for (E_Int no = 0; no < nzonesr; no++)
          RELEASESHAREDA(resr[no],objrs[no],fieldrs[no],ar2[no],ar3[no],ar4[no]);  
        for (E_Int no = 0; no < nzones; no++)
          RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);   
        PyErr_SetString(PyExc_TypeError,
                        "setInterpDataDW: structured donor zones must be 3D or 2D with z=constant.");
        return NULL;
      }
    }
  }
  else //if (hook != Py_None) // hook fourni
  {
#if (PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION < 7) || (PY_MAJOR_VERSION == 3 && PY_MINOR_VERSION < 1)
    void** packet = (void**) PyCObject_AsVoidPtr(hook);
#else
    void** packet = (void**) PyCapsule_GetPointer(hook, NULL);
#endif
    E_Int s1 = fields.size();
    for (E_Int i = 0; i < s1; i++) 
      interpDatas.push_back((K_INTERP::InterpAdt*)(packet[i+1])); 
  }
   // cas seulement non structure : ncfmax a 4 (minimum)
  if (nzonesU != 0)
  {
    if (interporder != 2) printf("Warning: setInterpDataDW: interpolation order is 2 for tetra arrays.\n");
    if (nzonesS == 0) ncfmax = 4;
  }
  /*-------------------------------------------------------*/
  /* Calcul des coefficients d interpolation               */
  /*-------------------------------------------------------*/
  E_Int nbI = fieldrs[0]->getSize();
  E_Float vol;

  // Donnees d interpolation
  vector<FldArrayF*> listOfInterpCoefs;
  vector<FldArrayI*> listOfDonorInd1D;
  vector<FldArrayI*> listOfRcvInd1D;
  vector<FldArrayI*> listOfDonorTypes;
  vector<FldArrayI*> listOfExtrapInd1D;
  vector<FldArrayF*> listOfEXdirs;

  for (E_Int noz = 0; noz < nzones; noz++)
  {
    FldArrayF* donorInterpCoefs = new FldArrayF(nbI*ncfmax); donorInterpCoefs->setAllValuesAtNull();
    listOfInterpCoefs.push_back(donorInterpCoefs);
    FldArrayI* donorInd1D = new FldArrayI(nbI*2); donorInd1D->setAllValuesAt(-1);
    listOfDonorInd1D.push_back(donorInd1D);
    FldArrayI* rcvInd1D = new FldArrayI(nbI); rcvInd1D->setAllValuesAt(-1);
    listOfRcvInd1D.push_back(rcvInd1D);
    FldArrayI* donorType = new FldArrayI(nbI); donorType->setAllValuesAtNull();
    listOfDonorTypes.push_back(donorType);
    FldArrayI* extrapPoints = new FldArrayI(nbI); extrapPoints->setAllValuesAt(-1);
    listOfExtrapInd1D.push_back(extrapPoints);
    if (isEX == 1) 
    {
      FldArrayF* EXdirs = new FldArrayF(nbI); EXdirs->setAllValuesAtNull();
      listOfEXdirs.push_back(EXdirs);
    }
  }
  
  FldArrayI* orphanPts = new FldArrayI(nbI); orphanPts->setAllValuesAt(-1);
  E_Int type, isExtrapolated;
  FldArrayI sizeOfIndDonor1D(nzones); sizeOfIndDonor1D.setAllValuesAtNull();
  FldArrayI usedDomp(nzones); usedDomp.setAllValuesAtNull();//nb de pts interpoles par domaine donneur
  FldArrayI sizeCoefs(nzones); sizeCoefs.setAllValuesAtNull();// taille du tableau de coefficients par domaine donneur
  FldArrayI nbExtraps(nzones); nbExtraps.setAllValuesAtNull();
  // pour l instant un seul indice 1D est stocke par pt interpole, donc pas besoin de specifier la taille de donorInd1D par zoneD  
  E_Int noblk = 0;
  vector<E_Float> xt(nzones); vector<E_Float> yt(nzones); vector<E_Float> zt(nzones);
  E_Float* EXdir0 = NULL;
  if (isEX == 1) EXdir0 = fieldrs[0]->begin(posdir[0]);
  E_Int noOrphanPt=0;
  for (E_Int ind = 0; ind < nbI; ind++)
  {
    for (E_Int noz = 0; noz < nzones; noz++)
    {
      E_Float* xr = fieldrs[noz]->begin(posxr[noz]);
      E_Float* yr = fieldrs[noz]->begin(posyr[noz]);
      E_Float* zr = fieldrs[noz]->begin(poszr[noz]);      
      xt[noz] = xr[ind]; yt[noz] = yr[ind]; zt[noz] = zr[ind]; 
    }
    short ok = K_INTERP::getInterpolationCellDW(
      &xt[0], &yt[0], &zt[0], interpDatas, fields,
      a2, a3, a4, a5, posxs, posys, poszs, poscs,
      vol, indi, cf, type, noblk, interpType, nature, penalty);   
    isExtrapolated = 0;
    if  (ok != 1)
    {
      ok = K_INTERP::getExtrapolationCellDW(
        &xt[0], &yt[0], &zt[0], interpDatas, fields,
        a2, a3, a4, a5, posxs, posys, poszs, poscs,
        vol, indi, cf, type, noblk, interpType, nature, penalty); 
      if (noblk > 0) isExtrapolated = 1;
    }

    if (noblk > 0)
    {
      E_Int noDonorBlk = noblk-1;
      E_Int& noInterpPtForBlk = usedDomp[noDonorBlk];// on met une reference car on va modifier cette valeur
      E_Int& sizeOfIndDonor1DForBlk = sizeOfIndDonor1D[noDonorBlk];
      E_Int& sizecf = sizeCoefs[noDonorBlk];//idem
      E_Int& noExtrapPtForBlk = nbExtraps[noDonorBlk]; 
      E_Float* donorCf = listOfInterpCoefs[noDonorBlk]->begin();
      E_Int* donorType = listOfDonorTypes[noDonorBlk]->begin();
      E_Int* donorInd1D = listOfDonorInd1D[noDonorBlk]->begin();
      E_Int* rcvInd1D = listOfRcvInd1D[noDonorBlk]->begin();
      E_Int* extrapInd1D = listOfExtrapInd1D[noDonorBlk]->begin();
      E_Float* EXdir = NULL;
      if (isEX == 1)  
      {
        EXdir = listOfEXdirs[noDonorBlk]->begin();
        EXdir[noInterpPtForBlk] = EXdir0[ind];
      }
      donorType[noInterpPtForBlk] = type;            
      rcvInd1D[noInterpPtForBlk] = ind;
      switch (type)
      {
        case 1:
          donorCf[sizecf] = cf[0]; sizecf += 1;
          donorInd1D[sizeOfIndDonor1DForBlk+1] = indi[0];
          sizeOfIndDonor1DForBlk++;
          break;

        case 2:
          for (E_Int nocf = 0; nocf < 8; nocf++) 
            donorCf[sizecf+nocf] = cf[nocf];
          sizecf += 8;
          donorInd1D[sizeOfIndDonor1DForBlk] = indi[0]; 
          sizeOfIndDonor1DForBlk++;
          break;

        case 22:
          for (E_Int nocf = 0; nocf < 4; nocf++) 
            donorCf[sizecf+nocf] = cf[nocf];
          sizecf += 4;
          donorInd1D[sizeOfIndDonor1DForBlk] = indi[0]; 
          sizeOfIndDonor1DForBlk++;
          break;

        case 3: 
          for (E_Int nocf = 0; nocf < 9; nocf++) 
            donorCf[sizecf+nocf] = cf[nocf];
          sizecf += 9;
          donorInd1D[sizeOfIndDonor1DForBlk] = indi[0]; 
          sizeOfIndDonor1DForBlk++;
          break;

        case 4:
          if (a4[noblk-1] == NULL)//non structure
          {
            for (E_Int nocf = 0; nocf < 4; nocf++) 
              donorCf[sizecf+nocf] = cf[nocf];
            sizecf += 4;
            donorInd1D[sizeOfIndDonor1DForBlk] = indi[0];
            sizeOfIndDonor1DForBlk++;
          }
          else // structure
          {
            for (E_Int nocf = 0; nocf < 8; nocf++) 
              donorCf[sizecf+nocf] = cf[nocf];
            sizecf += 8;
            donorInd1D[sizeOfIndDonor1DForBlk] = indi[0];
            sizeOfIndDonor1DForBlk++;// cas 2,3,4,5 
          } 
          break;

        case 5:
          for (E_Int nocf = 0; nocf < 15; nocf++) 
            donorCf[sizecf+nocf] = cf[nocf];
          sizecf += 15;
          donorInd1D[sizeOfIndDonor1DForBlk] = indi[0]; 
          sizeOfIndDonor1DForBlk++;
          break;

        default:
          printf("Error: setInterpData: type not yet implemented.\n");
          exit(0);
      }      
      if (isExtrapolated == 1) 
      { 
        extrapInd1D[noExtrapPtForBlk] = ind; 
        noExtrapPtForBlk++;// = nbExtraps[noDonorBlk]++;
      }
      noInterpPtForBlk++;// = usedDomp[noDonorBlk]++;
    }
    else {(*orphanPts)[noOrphanPt] = ind; noOrphanPt++;} 
  }
  // Redimensionnement des tableaux 
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    E_Int sizecf = sizeCoefs[noz];
    E_Int nbInterp = usedDomp[noz];
    E_Int sizeOfIndDonor1DL = sizeOfIndDonor1D[noz];
    E_Int nbExtrap = nbExtraps[noz];    
    listOfInterpCoefs[noz]->resize(sizecf);
    listOfDonorInd1D[noz]->resize(sizeOfIndDonor1DL);
    listOfRcvInd1D[noz]->resize(nbInterp);
    listOfDonorTypes[noz]->resize(nbInterp);
    listOfExtrapInd1D[noz]->resize(nbExtrap);
    if (isEX == 1) listOfEXdirs[noz]->resize(nbInterp);
  }
  // Nettoyages
  for (E_Int no = 0; no < nzones; no++)
  {
    if (hook == Py_None) delete interpDatas[no];
    RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);  
  }
  for (E_Int no = 0; no < nzonesr; no++)
    RELEASESHAREDA(resr[no],objrs[no],fieldrs[no],ar2[no],ar3[no],ar4[no]); 
  // Fin nettoyages
  orphanPts->resize(noOrphanPt);

  /*----------------------------------------------------------*/
  /* Ecriture dans des objets Python retournes par la methode */
  /*----------------------------------------------------------*/
  // listes pour stocker les indices de cellules (receveuse et donneuses) par bloc de domaine d'interpolation correspondant
  PyObject* PyListCellIndicesR = PyList_New(0);
  PyObject* PyListCellIndicesD = PyList_New(0);
  // liste pour stocker les coefficients d'interpolation par bloc de domaine d'interpolation correspondant
  PyObject* PyListCoefficients = PyList_New(0);
  // liste pour stocker les types d'interpolation par domaine d'interpolation correspondant
  PyObject* PyListInterpTypes = PyList_New(0);
  // listes pour marquer les indices des cellules extrapoles
  PyObject * PyListCellIndicesE = PyList_New(0);
  // listes pour marquer les EXdir (directions pour les pts EX)
  PyObject * PyListEXDir = PyList_New(0);

  // ------------------------------
  // Construction des PyArrays
  // ------------------------------
  for (E_Int noz = 0; noz < nzones; noz++)
  {   
    //    coefficients d interpolation
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

    if (isEX == 1)
    {
      PyObject* pyEXDir = K_NUMPY::buildNumpyArray(*listOfEXdirs[noz],1);
      PyList_Append(PyListEXDir, pyEXDir); Py_DECREF(pyEXDir);
      delete  listOfEXdirs[noz];
    }

  } // fin parcours des zones donneuses
  
  //   indices des points orphelins
  PyObject* cellIndO = K_NUMPY::buildNumpyArray(*orphanPts,1);
  delete orphanPts;

  PyObject* tpl = Py_BuildValue("[OOOOOOO]", PyListCellIndicesR, PyListCellIndicesD, PyListInterpTypes, PyListCoefficients,
                                PyListCellIndicesE, cellIndO, PyListEXDir);
  Py_DECREF(cellIndO); Py_DECREF(PyListInterpTypes);
  Py_DECREF(PyListCellIndicesE); Py_DECREF(PyListEXDir); 
  Py_DECREF(PyListCoefficients); Py_DECREF(PyListCellIndicesR); 
  Py_DECREF(PyListCellIndicesD); 
  return tpl;
}
