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
# include <vector>
# include "connector.h"
# include "KInterp/BlkInterp.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* 
   Calcul des transferts Chimere  issus de setInterpolations 
   Le calcul s effectue par paire de bloc interpole/bloc donneur
   On suppose que seules les variables communes aux blocs donneurs et receveurs 
   sont transferees, sauf le cellN  
*/
//=============================================================================
PyObject* K_CONNECTOR::chimeraTransfer(PyObject* self, PyObject* args)
{
  PyObject *pyIndRcv, *pyIndDonor;
  PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs;
  PyObject *arrayR, *arrayD;
  if (!PYPARSETUPLE_(args, OOOO_ OO_,
                    &pyIndRcv, &pyIndDonor, &pyArrayTypes, &pyArrayCoefs, &arrayR, &arrayD))
  {
    return NULL;
  }

  /*--------------------------------------*/
  /* Extraction des indices des receveurs */
  /*--------------------------------------*/
  FldArrayI* rcvPtsI;
  E_Int res = K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI);
  if (res == 0) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "chimeraTransfer: 1st arg must be a numpy of integers.");
    return NULL;
  }
  /*-------------------------------------*/
  /* Extraction des indices des donneurs */
  /*-------------------------------------*/
  FldArrayI* donorPtsI;
  res = K_NUMPY::getFromNumpyArray(pyIndDonor, donorPtsI);
  if (res == 0) 
  {
    RELEASESHAREDN(pyIndRcv, rcvPtsI);
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpTransfers: 2nd arg must be a numpy of integers.");
    return NULL;
  }
  /*----------------------*/
  /* Extraction des types */
  /*----------------------*/
  FldArrayI* typesI;
  res = K_NUMPY::getFromNumpyArray(pyArrayTypes, typesI);
  if (res == 0) 
  {
    RELEASESHAREDN(pyIndRcv, rcvPtsI);
    RELEASESHAREDN(pyIndDonor, donorPtsI);
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpTransfers: 3rd arg must be a numpy of integers.");
    return NULL;
  }

  /*-----------------------*/
  /* Extraction des coefs  */
  /*-----------------------*/
  FldArrayF* donorCoefsF;
  res = K_NUMPY::getFromNumpyArray(pyArrayCoefs, donorCoefsF);

  if (res == 0) 
  {
    RELEASESHAREDN(pyIndRcv, rcvPtsI);
    RELEASESHAREDN(pyIndDonor, donorPtsI);
    RELEASESHAREDN(pyArrayTypes, typesI);
    PyErr_SetString(PyExc_TypeError, 
                    "setInterpTransfers: 4th arg must be a numpy of floats.");
    return NULL;
  }
  // getFromNumpyArray echange la "shape" du FldArray quand il n y a 
  // qu un seul donneur donc on corrige ici.
  if (donorCoefsF->getNfld() == 1)
  {
    donorCoefsF->resize(donorCoefsF->getNfld(),donorCoefsF->getSize());
  }
  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  E_Int imr, jmr, kmr;
  FldArrayF* fr; FldArrayI* cnr;
  char* varStringr; char* eltTyper;
  E_Int resr = K_ARRAY::getFromArray3(arrayR, varStringr, fr, 
                                      imr, jmr, kmr, cnr, eltTyper); 
  if (resr != 1) 
  {
    RELEASESHAREDN(pyIndRcv, rcvPtsI);
    RELEASESHAREDN(pyIndDonor, donorPtsI);
    RELEASESHAREDN(pyArrayTypes, typesI);
    RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
    if (resr == 2) RELEASESHAREDB(resr, arrayR, fr, cnr); 

    PyErr_SetString(PyExc_TypeError,
                    "chimeraTransfer: 5th arg is not a valid array.");
    return NULL; 
  }
  /*---------------------------------------------*/
  /* Extraction des infos sur le domaine donneur */
  /*---------------------------------------------*/
  E_Int imd, jmd, kmd;
  FldArrayF* fd; FldArrayI* cnd;
  char* varStringd; char* eltTyped;
  E_Int resd = K_ARRAY::getFromArray3(arrayD, varStringd, fd, 
                                      imd, jmd, kmd, cnd, eltTyped); 
  if (resd != 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "chimeraTransfers: 6th arg is not a valid array.");
    RELEASESHAREDB(resr, arrayR, fr, cnr); 
    RELEASESHAREDN(pyIndRcv, rcvPtsI);
    RELEASESHAREDN(pyIndDonor, donorPtsI);
    RELEASESHAREDN(pyArrayTypes, typesI);
    RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
    if (resd == 2) RELEASESHAREDB(resd, arrayD, fd, cnd); 
    return NULL; 
  }
  /*---------------------------------------------------*/
  /*  Extrait les positions des variables a transferer */
  /*---------------------------------------------------*/
  vector<E_Int> posvarsD; vector<E_Int> posvarsR;
  E_Int poscd = K_ARRAY::isCellNatureField2Present(varStringd)+1;
  E_Int poscr = K_ARRAY::isCellNatureField2Present(varStringr)+1;
  char* varStringC; // chaine de caractere commune
  E_Int l = strlen(varStringr);
  varStringC = new char [l+1];
  K_ARRAY::getPosition(varStringr, varStringd, posvarsR, posvarsD, varStringC);
  
  delete [] varStringC; 
  if (poscd != 0) 
    posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), poscd), posvarsD.end());
  if (poscr != 0)
    posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), poscr), posvarsR.end());    
  
  /*---------------------------------------------------------------*/
  /* Application de la formule d'interpolation                     */
  /*---------------------------------------------------------------*/
  E_Int api = fr->getApi();
  PyObject* tpl = K_ARRAY::buildArray3(*fr, varStringr, imr, jmr, kmr, api);
  FldArrayF* fieldROut;
  K_ARRAY::getFromArray3(tpl, fieldROut);

  E_Int nbRcvPts = rcvPtsI->getSize();
  E_Int ncoefs = 8;
  E_Int indR, indD, indDLoc;
  E_Int i, j, k, nocf;
  E_Int imdjmd = imd*jmd;
  E_Int* rcvPts = rcvPtsI->begin();
  E_Int* donorPts = donorPtsI->begin();
  E_Int* types = typesI->begin();
  E_Int nfldr0 = posvarsR.size();
  E_Float sumCf;
  vector<E_Float> cfLoc(ncoefs);
  E_Int dim3 = 1; // 3d problem
  if (kmd == 1) dim3 = 0;  // 2d problem
  for (E_Int noind = 0; noind < nbRcvPts; noind++)
  {
    E_Int type = types[noind];
    switch (type)
    {
      case 100:
        for (E_Int no=0; no < ncoefs; no++)
        {
          cfLoc[no] = (*donorCoefsF)(noind,no+1);
        }
 
        indR = rcvPts[noind];
        indD = donorPts[noind];

        k = indD/imdjmd;
        j = (indD-k*imdjmd)/imd;
        i = (indD-j*imd-k*imdjmd);
        for (E_Int eq = 0; eq < nfldr0; eq++)
        {
          E_Float* fieldR = fieldROut->begin(posvarsR[eq]);
          E_Float* fieldD = fd->begin(posvarsD[eq]);
          fieldR[indR] = 0.; nocf = 0; sumCf = 0.;
          for (E_Int kk=0; kk<2; kk++)
            for (E_Int jj=0; jj<2; jj++)
              for (E_Int ii=0; ii<2; ii++)
              {
                indDLoc = (i+ii)+(j+jj)*imd+dim3*(k+kk)*imdjmd;
                sumCf += cfLoc[nocf]; 
                fieldR[indR] += cfLoc[nocf]*fieldD[indDLoc];
                nocf++;
              }
       }
        break;
        
      default:
        printf("Warning: chimeraTransfer: invalid interpolation type; fields are not transfered.\n");
        break;
    }
  }
  // sortie
  RELEASESHAREDS(tpl, fieldROut);
  RELEASESHAREDB(resr, arrayR, fr, cnr); 
  RELEASESHAREDB(resd, arrayD, fd, cnd); 
  RELEASESHAREDN(pyIndRcv, rcvPtsI);
  RELEASESHAREDN(pyIndDonor, donorPtsI);
  RELEASESHAREDN(pyArrayTypes, typesI);
  RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
  return tpl;
}
