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
#include "Nuga/include/conservative_chimera.h"

#define CLOUDMAX 50

using namespace std;
using namespace K_FLD;
                                 

#define RELEASEDATA                                                     \
  RELEASESHAREDB(resr, receiverArray, fr, cnr);                         \
  RELEASESHAREDN(numpyIndicesR, IndicesR);                              \
  for (E_Int no = 0; no < nDnrZones; no++)                              \
  { RELEASESHAREDA(resl[no],objs[no],fields[no],a2[no],a3[no],a4[no]);  \
    PyObject* tpld = PyList_GetItem(listOfNumpyIndicesD, no);           \
    RELEASESHAREDN(tpld, listOfIndicesD[no]);  } 

//=============================================================================
/* Calcule et stocke les coefficients d'interpolation par intersection
   OUT: [donorBlks,donorInd1D, donorType, coefs, extrapInd1D, orphanInd1D] 
        donorBlks: no du blk donneur, demarre a 0
        donorInd1D: indice global (structure), de l elt (NS) du donneur
        donorType: type d interpolation effectue localement
        coefs: coefficients d interpolation, stockes selon le type
        extrapInd1D: indices des pts extrapoles
        orphanInd1D: indices des pts orphelins */
//=============================================================================
PyObject* K_CONNECTOR::setInterpDataCons(PyObject* self, PyObject* args)
{
  PyObject *receiverArray; 
  PyObject *donorArrays;// domaines d interpolation
  PyObject *numpyIndicesR;
  PyObject *listOfNumpyIndicesD;   
  if (!PYPARSETUPLE_(args, OOOO_, &receiverArray, &donorArrays, &numpyIndicesR, &listOfNumpyIndicesD)) 
    return NULL;

  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  E_Int imr, jmr, kmr;
  FldArrayF* fr; FldArrayI* cnr;
  char* varStringr; char* eltTyper;
  E_Int resr = K_ARRAY::getFromArray3(receiverArray, varStringr, fr, 
                                      imr, jmr, kmr, cnr, eltTyper); 
  if (resr == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setInterpDataCons: 1st arg is not valid.");
    return NULL;
  }
  else if (resr == 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setInterpDataCons: 1st arg must be unstructured.");
    RELEASESHAREDS(receiverArray, fr); 
    return NULL;
  }    
  else//check NGON
  {
    if (K_STRING::cmp(eltTyper,"NGON") != 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "setInterpDataCons: 1st arg must be NGON.");
      RELEASESHAREDB(resr, receiverArray, fr, cnr); 
      return NULL;
    }
  }
  // Verif des coordonnees dans la zone a interpoler
  E_Int posxr = K_ARRAY::isCoordinateXPresent(varStringr);
  E_Int posyr = K_ARRAY::isCoordinateYPresent(varStringr);
  E_Int poszr = K_ARRAY::isCoordinateZPresent(varStringr);
  if (posxr == -1 || posyr == -1 || poszr == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setInterpDataCons: 1st arg must contain coordinates.");
    RELEASESHAREDB(resr, receiverArray, fr, cnr); 
    return NULL;
  }

  // Extraction of indirection tab for receptor zone 
  FldArrayI* IndicesR;// structured indices corresponding to NGON indices of the receptor zone
  E_Int resi = K_NUMPY::getFromNumpyArray(numpyIndicesR, IndicesR);
  if (resi == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setInterpDataCons: 2nd arg must be a numpy of integers.");
    RELEASESHAREDB(resr, receiverArray, fr, cnr); 
    return NULL;
  }
  // Extraction of indirection tab for donor zones
  if (PyList_Check (listOfNumpyIndicesD) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpDataCons: 4th arg must be a list.");
    RELEASESHAREDB(resr, receiverArray, fr, cnr); 
    RELEASESHAREDN(numpyIndicesR, IndicesR); 
    return NULL;
  }
  vector<FldArrayI*> listOfIndicesD;// list of numpy of structured indices corresponding to NGON indices of all the donor zones
  E_Int sizeOfDnrs = PyList_Size(listOfNumpyIndicesD);
  for (E_Int no = 0; no < sizeOfDnrs; no++)
  {  
    FldArrayI* IndicesD;
    PyObject* tpli = PyList_GetItem(listOfNumpyIndicesD, no);
    resi = K_NUMPY::getFromNumpyArray(tpli, IndicesD);
    if (resi == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "setInterpDataCons: 2nd arg must be a numpy of integers.");
      RELEASESHAREDB(resr, receiverArray, fr, cnr); 
      RELEASESHAREDN(numpyIndicesR, IndicesR);   
      for (E_Int no2 = 0; no2 < no; no2++)
      {
        PyObject* tpli2 = PyList_GetItem(listOfNumpyIndicesD, no2);
        RELEASESHAREDN(tpli2, listOfIndicesD[no2]);   
      }
      return NULL;
    }
    listOfIndicesD.push_back(IndicesD);
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
  E_Bool skipNoCoord = true;  E_Bool skipStructured = true;
  E_Bool skipUnstructured = false;  E_Bool skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(donorArrays, resl, varString, fields, a2, a3, a4, objs,  
                                      skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nDnrZones = objs.size();
  if (isOk == -1)
  {
    RELEASEDATA;
    PyErr_SetString(PyExc_TypeError,
                    "setInterpDataCons: 2nd argument is not valid.");
    return NULL;
  }   
  if (nDnrZones == 0) 
  {
    RELEASEDATA;
    PyErr_SetString(PyExc_TypeError,
                    "setInterpDataCons: no valid donor zone found.");
    return NULL;
  }
  for (E_Int noz = 0; noz < nDnrZones; noz++)
  {
    if (resl[noz] == 2) 
    {
      char* eltType0 = (char*)a3[noz];
      if (K_STRING::cmp(eltType0,"NGON")!= 0)
      {
        RELEASEDATA;
        PyErr_SetString(PyExc_TypeError,
                        "setInterpDataCons: unstructured donor zones must be NGON.");
        return NULL;
      }
    }
  }
  // Check variables
  vector<E_Int> posxd; vector<E_Int> posyd; vector<E_Int> poszd;
  for (E_Int no = 0; no < nDnrZones; no++)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString[no]); posxd.push_back(posx); 
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString[no]); posyd.push_back(posy); 
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString[no]); poszd.push_back(posz); 
  }


  /*-------------------------------------------------------*/
  /* Calcul des coefficients d'interpolation               */
  /*-------------------------------------------------------*/
  E_Int* cngR = cnr->begin();//connectivite ngon recepteur
  E_Int sizeFN = cngR[1];
  E_Int neltsR = cngR[sizeFN+2];


  //E_Int nbI = neltsR; //nb d elts a interpoler (initialisation, a redimensionner en sortie)

  // Initialisation des tableaux et dimensionnement a priori
  // On met pour l'instant NCloudPtsMax a CLOUDMAX: la molecule d'interpolation contient a priori CLOUDMAX pts
  //E_Int nCloudPtsMax = CLOUDMAX;

  vector<E_Int> cellNt(neltsR);
  for (E_Int i = 0; i < neltsR; i++) cellNt[i]=2;

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

  vector<FldArrayI*> listOfDonorInd1D;
  vector<FldArrayI*> listOfDonorTypes;
  vector<FldArrayF*> listOfInterpCoefs;
  vector<FldArrayI*> listOfRcvInd1D;
  FldArrayI isInterp(neltsR); isInterp.setAllValuesAtNull();
  E_Int* isInterpPtr = isInterp.begin();

  E_Int nocf, noi, index;
  E_Int* indicesRcvOrig = IndicesR->begin();
  for (E_Int noz = 0; noz < nDnrZones; noz++)
  {
    E_Int* indicesDnrOrig = listOfIndicesD[noz]->begin();
    FldArrayF* fd = fields[noz];//nodes = coordinates
    FldArrayI* cnd = (FldArrayI*)a2[noz]; 
    vector<E_Int> dindices, xr, roids;
    vector<E_Float> dcoeffs;
    // Interpolation de (xr,yr,zr) par le donneur
    E_Int err = 0;

    // Nuga ne compile pas en DOUBLE_INT
    err = NUGA::P1_CONSERVATIVE::compute_chimera_coeffs(*fr, posxr, posyr, poszr, *cnr, 
                                                        *fd, posxd[noz], posyd[noz], poszd[noz], *cnd, cellNt,
                                                        dindices, dcoeffs, xr, roids);

    if (err)
    {
      FldArrayI* donorInd1D = new FldArrayI(0);
      FldArrayI* donorType = new FldArrayI(0); 
      FldArrayF* coefs = new FldArrayF(0); 
      FldArrayI* rcvInd1D = new FldArrayI(0);
      
      //     donor indices + nb of donors per pt 
      listOfDonorTypes.push_back(donorType);
      listOfDonorInd1D.push_back(donorInd1D);
      listOfRcvInd1D.push_back(rcvInd1D);
      listOfInterpCoefs.push_back(coefs);
      printf("Warning: setInterpDataConservative: failure in computation of interpolation coefficients for donor zone " SF_D_ " (" SF_D_ " zones).\n", noz, nDnrZones);
      // RELEASEDATA; return Py_None;
    }
    else 
    {
      E_Int sizeOfRcvIndices =roids.size();// nb of interpolated points by donor zone noz
      E_Int sizeOfDnrIndices=dindices.size();// size of dnr indices array
      E_Int sizeCf = dcoeffs.size();//size of coefs array
      E_Int sizeOfNbCfsPerPt = xr.size()-1;// size of hash table
      FldArrayI* donorInd1D = new FldArrayI(sizeOfDnrIndices+sizeOfNbCfsPerPt); donorInd1D->setAllValuesAt(-1);
      FldArrayI* donorType = new FldArrayI(sizeOfRcvIndices); donorType->setAllValuesAtNull();//vaudra 0 
      FldArrayF* coefs = new FldArrayF(sizeCf); coefs->setAllValuesAtNull();
      FldArrayI* rcvInd1D = new FldArrayI(sizeOfRcvIndices); rcvInd1D->setAllValuesAtNull();

      E_Int* ptrIndDnr = donorInd1D->begin();
      E_Int* ptrIndRcv = rcvInd1D->begin();
      E_Float* ptrCf = coefs->begin();

      nocf = 0; noi = 0; E_Int nor = 0;
      for (E_Int noindR = 0; noindR < sizeOfRcvIndices; noindR++)
      {
        index = roids[noindR];
        ptrIndRcv[noindR] = indicesRcvOrig[index]-1;
        isInterpPtr[roids[noindR]]=1;
        ptrIndDnr[noi] = xr[noindR+1]-xr[noindR]; 
        nor++; noi++;
        for (E_Int no = xr[noindR]; no < xr[noindR+1]; no++)
        {
          index = dindices[no];          
          ptrIndDnr[noi] = indicesDnrOrig[index]-1;
          ptrCf[no] = dcoeffs[no];
          noi++; nocf++;
        }
      }
      //     donor indices + nb of donors per pt 
      listOfDonorTypes.push_back(donorType);
      listOfDonorInd1D.push_back(donorInd1D);
      listOfRcvInd1D.push_back(rcvInd1D);
      listOfInterpCoefs.push_back(coefs);
    } 
  }
  RELEASEDATA;

  FldArrayI* orphanPts = new FldArrayI(neltsR); orphanPts->setAllValuesAt(-1);
  E_Int* orphanPts1D = orphanPts->begin();
  E_Int nbOrphans = 0;
  isInterpPtr = isInterp.begin();
  for (E_Int nor = 0; nor < neltsR; nor++)
  {
    if (isInterpPtr[nor]==0) 
      {
        orphanPts1D[nbOrphans] = nor;
        nbOrphans++;
      }
  }
  orphanPts->resize(nbOrphans);
  PyObject* cellIndO = K_NUMPY::buildNumpyArray(*orphanPts,1);
  delete orphanPts;

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
  } 

  PyObject* l = Py_BuildValue("[OOOOO]", 
                              PyListCellIndicesR, PyListCellIndicesD, 
                              PyListInterpTypes, PyListCoefficients, cellIndO);
  Py_DECREF(PyListCellIndicesR); 
  Py_DECREF(PyListCellIndicesD);
  Py_DECREF(PyListInterpTypes); 
  Py_DECREF(PyListCoefficients); 
  Py_DECREF(cellIndO); 
  return l;
}
