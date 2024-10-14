/*
  Copyright 2013-2024 Onera.

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
# include "param_solver.h"

#include <stdlib.h>

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Effectue les transferts par interpolation */
//=============================================================================
PyObject* K_CONNECTOR::setInterpTransfers(PyObject* self, PyObject* args)
{
  PyObject *arrayR, *arrayD;
  PyObject *pyIndRcv, *pyIndDonor;
  PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs;
  PyObject *pyVariables;
  E_Float AngleX, AngleY, AngleZ;
  if (!PYPARSETUPLE_(args, OOOO_ OOO_ TRRR_, 
                     &arrayR, &arrayD,  &pyVariables, &pyIndRcv, &pyIndDonor, &pyArrayTypes, &pyArrayCoefs,
                     &AngleX, &AngleY, &AngleZ))
  {
    return NULL;
  }

  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  E_Int imr, jmr, kmr;
  FldArrayF* fr; FldArrayI* cnr;
  char* varStringR; char* eltTypeR;
  E_Int resr = K_ARRAY::getFromArray(arrayR, varStringR, fr,
                                     imr, jmr, kmr, cnr, eltTypeR, true);
  if (resr != 2 && resr != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setInterpTransfers: 1st arg is not a valid array.");
    return NULL;
  }
  /*---------------------------------------------*/
  /* Extraction des infos sur le domaine donneur */
  /*---------------------------------------------*/
  E_Int imd, jmd, kmd, imdjmd;
  FldArrayF* fd; FldArrayI* cnd;
  char* varStringD; char* eltTypeD;
  E_Int resd = K_ARRAY::getFromArray(arrayD, varStringD, fd,
                                     imd, jmd, kmd, cnd, eltTypeD, true);
  if (resd != 2 && resd != 1)
  {
    PyErr_SetString(PyExc_TypeError,
              "setInterpTransfers: 2nd arg is not a valid array.");
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    return NULL;
  }

  E_Int meshtype = resd; // 1: structure, 2: non structure

  if (resd == 2)
  {
    if (K_STRING::cmp(eltTypeD, "TETRA") != 0 &&
        K_STRING::cmp(eltTypeD, "NGON") != 0)
    {
      PyErr_SetString(PyExc_TypeError,
              "setInterpTransfers: unstructured donor zone must only be TETRA or NGON.");
      RELEASESHAREDB(resr, arrayR, fr, cnr);
      RELEASESHAREDB(resd, arrayD, fd, cnd);
      return NULL;
    }
  }

  /*---------------------------------------------------*/
  /*  Extrait les positions des variables a transferer */
  /*---------------------------------------------------*/
  vector<E_Int> posvarsD; vector<E_Int> posvarsR;
  E_Int posvr, posvd;
  if (PyList_Check(pyVariables) != 0)
  {
    E_Int nvariables = PyList_Size(pyVariables);
    if (nvariables > 0)
    {
      for (int i = 0; i < nvariables; i++)
      {
        PyObject* tpl0 = PyList_GetItem(pyVariables, i);
        if (PyString_Check(tpl0))
        {
          char* varname = PyString_AsString(tpl0);
          posvd = K_ARRAY::isNamePresent(varname, varStringD);
          posvr = K_ARRAY::isNamePresent(varname, varStringR);
          if (posvd != -1 && posvr != -1)
          {
            posvarsD.push_back(posvd+1);
            posvarsR.push_back(posvr+1);
          }
        }
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(tpl0))
        {
          const char* varname = PyUnicode_AsUTF8(tpl0);
          posvd = K_ARRAY::isNamePresent(varname, varStringD);
          posvr = K_ARRAY::isNamePresent(varname, varStringR);
          if (posvd != -1 && posvr != -1)
          {
            posvarsD.push_back(posvd+1);
            posvarsR.push_back(posvr+1);
          }
        }
#endif
        else
          PyErr_Warn(PyExc_Warning, "setInterpTransfers: variable must be a string. Skipped.");
        }
      }
    else // toutes les variables communes sont transferees sauf le cellN
    {
      E_Int poscd = K_ARRAY::isCellNatureField2Present(varStringD)+1;
      E_Int poscr = K_ARRAY::isCellNatureField2Present(varStringR)+1;
      char* varStringC; // chaine de caractere commune
      E_Int l = strlen(varStringR);
      varStringC = new char [l+1];
      K_ARRAY::getPosition(varStringR, varStringD,
                   posvarsR, posvarsD, varStringC);

      delete [] varStringC;
      if (poscd != 0)
        posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), poscd), posvarsD.end());
      if (poscr != 0)
        posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), poscr), posvarsR.end());
    }
    }
  else
  {
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    RELEASESHAREDB(resd, arrayD, fd, cnd);
    PyErr_SetString(PyExc_TypeError,
              "setInterpTransfers: name of transfered variables must be defined by a list.");
    return NULL;
  }

# include "extract_interpD.h"
  /*--------------------------------------*/
  /* Extraction des indices des receveurs */
  /*--------------------------------------*/
  FldArrayI* rcvPtsI;
  E_Int res_rcv = K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI, true);
  nbRcvPts      = rcvPtsI->getSize();
  E_Int* rcvPts = rcvPtsI->begin();

  if (res_donor*res_type*res_coef*res_rcv ==0)
  {
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    RELEASESHAREDB(resd, arrayD, fd, cnd);
    if (res_donor != 0) { RELEASESHAREDN(pyIndDonor  , donorPtsI  );}
    if (res_type  != 0) { RELEASESHAREDN(pyArrayTypes, typesI     );}
    if (res_coef  != 0) { RELEASESHAREDN(pyArrayCoefs, donorCoefsF);}
    if (res_rcv   != 0) { RELEASESHAREDN(pyIndRcv    , rcvPtsI    );}
    PyErr_SetString(PyExc_TypeError,"setInterpTransfers: 4th to 6th arg must be a numpy of integers. 7th arg a numpy floats ");
    return NULL;
  }

  // Extraction de l angle de rotation
  E_Int dirR = 0; E_Float theta = 0.;
  if (K_FUNC::E_abs(AngleX) > 0.) {dirR = 1; theta=AngleX;}
  else if (K_FUNC::E_abs(AngleY) > 0.){dirR=2; theta=AngleY;}
  else if (K_FUNC::E_abs(AngleZ) > 0.) {dirR=3; theta=AngleZ;}
  E_Int posvx=-1, posvy=-1, posvz=-1, posmx=-1, posmy=-1, posmz=-1;
  if (dirR > 0)
  {
    posvx = K_ARRAY::isNamePresent("VelocityX", varStringR);
    posvy = K_ARRAY::isNamePresent("VelocityY", varStringR);
    posvz = K_ARRAY::isNamePresent("VelocityZ", varStringR);
    posmx = K_ARRAY::isNamePresent("MomentumX", varStringR);
    posmy = K_ARRAY::isNamePresent("MomentumY", varStringR);
    posmz = K_ARRAY::isNamePresent("MomentumZ", varStringR);
  }

  E_Int nvars   = fr->getNfld();//nb de champs a interpoler
  E_Int* ptrcnd = cnd->begin();
  E_Int cnNfldD = cnd->getNfld();

  PyObject* tpl;
  if (resr == 1)
  {
    tpl = K_ARRAY::buildArray(nvars, varStringR, imr, jmr, kmr);
  }
  else // unstructured
  {
    E_Int crsize = cnr->getSize()*cnr->getNfld();
    tpl = K_ARRAY::buildArray(nvars, varStringR,
                fr->getSize(), cnr->getSize(),
                -1, eltTypeR, false, crsize);
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cnr->begin(), cnr->getSize()*cnr->getNfld());
  }

  E_Float* frp  = K_ARRAY::getFieldPtr(tpl);
  //FldArrayF fieldROut(nbRcvPts, nvars, frp, true);
  FldArrayF fieldROut( fr->getSize() , nvars, frp, true);
  fieldROut = *fr;

  // Transferts
  // Types valides: 0,1, 2, 3, 4, 5

  nvars = posvarsR.size();

  vector<E_Float*> vectOfRcvFields(nvars);
  vector<E_Float*> vectOfDnrFields(nvars);

  for (E_Int eq = 0; eq < nvars; eq++)
  {
    vectOfRcvFields[eq] = fieldROut.begin(posvarsR[eq]);
    vectOfDnrFields[eq] = fd->begin(posvarsD[eq]);
  }

  ////
  ////
  //  Interpolation parallele
  ////
  ////
# include "commonInterpTransfers_indirect.h"

  // Prise en compte de la periodicite par rotation
  if (dirR != 0)
  {
    # include "includeTransfers.h"
  }

  // sortie
  RELEASESHAREDB(resr, arrayR, fr, cnr);
  RELEASESHAREDB(resd, arrayD, fd, cnd);
  RELEASESHAREDN(pyIndRcv    , rcvPtsI    );
  RELEASESHAREDN(pyIndDonor  , donorPtsI  );
  RELEASESHAREDN(pyArrayTypes, typesI     );
  RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
  return tpl;
}

//=============================================================================
// Idem: in place + from zone
//=============================================================================
PyObject* K_CONNECTOR::_setInterpTransfers(PyObject* self, PyObject* args)
{
  char *GridCoordinates, *FlowSolutionNodes, *FlowSolutionCenters;

  PyObject *zoneR, *zoneD;
  PyObject *pyIndRcv, *pyIndDonor; PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs; PyObject *pyVariables;
  E_Int loc, vartype, compact;
  char* cellNVariable;
  E_Float AngleX, AngleY, AngleZ;

  if (!PYPARSETUPLE_(args, OOOO_ OOO_ III_ SSSS_ RRR_, 
                    &zoneR, &zoneD, &pyVariables, &pyIndRcv,
                    &pyIndDonor, &pyArrayTypes, &pyArrayCoefs, &loc , &vartype, &compact, &cellNVariable,
                    &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters,
                    &AngleX, &AngleY, &AngleZ))
  {
    return NULL;
  }

  // Extraction de l'angle de rotation
  E_Int dirR = 0; E_Float theta = 0.;
  if (K_FUNC::E_abs(AngleX) > 0.) {dirR = 1; theta=AngleX;}
  else if (K_FUNC::E_abs(AngleY) > 0.){dirR=2; theta=AngleY;}
  else if (K_FUNC::E_abs(AngleZ) > 0.) {dirR=3; theta=AngleZ;}

  vector<PyArrayObject*> hook;
  E_Int imdjmd, imd,jmd,kmd, cnNfldD, nvars,ndimdxR, ndimdxD,meshtype;
  E_Float* iptroD; E_Float* iptroR;

# include "extract_interpD.h"

  /*--------------------------------------*/
  /* Extraction des indices des receveurs */
  /*--------------------------------------*/
  FldArrayI* rcvPtsI;
  K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI, true);
  E_Int* rcvPts = rcvPtsI->begin();
  nbRcvPts = rcvPtsI->getSize();

  vector<E_Float*> fieldsR;vector<E_Float*> fieldsD;
  vector<E_Int> posvarsD; vector<E_Int> posvarsR;
  E_Int* ptrcnd;
  char* eltTypeR; char* eltTypeD;

  // codage general (lent ;-) )
  if (compact == 0)
  {
    // recupere les champs du donneur (nodes)
    E_Int cnSizeD;
    char* varStringD;
    vector<E_Int> locsD;
    vector<E_Int*> cnd;
    E_Int resd = K_PYTREE::getFromZone(zoneD, 0, 0, varStringD,
                     fieldsD, locsD, imd, jmd, kmd,
                     cnd, cnSizeD, cnNfldD,
                     eltTypeD, hook,
                     GridCoordinates,
                     FlowSolutionNodes, FlowSolutionCenters);
    
    if (cnd.size() > 0) ptrcnd = cnd[0];
    meshtype = resd; // 1: structure, 2: non structure
    // recupere les champs du receveur (centers)
    E_Int imr, jmr, kmr, cnSizeR, cnNfldR;
    char* varStringR; vector<E_Int> locsR;
    vector<E_Int*> cnr;
    K_PYTREE::getFromZone(zoneR, 0, loc, varStringR,
                fieldsR, locsR, imr, jmr, kmr,
                cnr, cnSizeR, cnNfldR, eltTypeR, hook,
                GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);
    if (varStringD == NULL)
    {
      RELEASESHAREDZ(hook, varStringD, eltTypeD);
      RELEASESHAREDN(pyIndRcv, rcvPtsI);
      RELEASESHAREDN(pyIndDonor, donorPtsI);
      RELEASESHAREDN(pyArrayTypes, typesI);
      RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
      PyErr_SetString(PyExc_TypeError,
              "_setInterpTransfers: no field found in donor zones.");
      return NULL;
    }
    if (varStringR == NULL)
    {
      RELEASESHAREDZ(hook, varStringR, eltTypeR);
      RELEASESHAREDN(pyIndRcv, rcvPtsI);
      RELEASESHAREDN(pyIndDonor, donorPtsI);
      RELEASESHAREDN(pyArrayTypes, typesI);
      RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
      PyErr_SetString(PyExc_TypeError,
              "_setInterpTransfers: no field found in receiver zones.");
      return NULL;
    }

    // Extrait les positions des variables a transferer
    E_Int initAll = false;
    E_Int poscd = K_ARRAY::isNamePresent(cellNVariable, varStringD);
    E_Int poscr = K_ARRAY::isNamePresent(cellNVariable, varStringR);
    E_Int posvr, posvd, posvarcr=-1, posvarcd=-1;

    if (PyList_Check(pyVariables) != 0)
    {
      E_Int nvariables = PyList_Size(pyVariables);
      if (nvariables > 0)
      {
        for (int i = 0; i < nvariables; i++)
        {
          PyObject* tpl0 = PyList_GetItem(pyVariables, i);
          if (PyString_Check(tpl0))
          {
            char* varname = PyString_AsString(tpl0);
            posvd = K_ARRAY::isNamePresent(varname, varStringD);
            posvr = K_ARRAY::isNamePresent(varname, varStringR);
            if (posvd == poscd) posvarcd = posvd;
            if (posvr == poscr) posvarcr = posvr;
            if (posvd != -1 && posvr != -1)
            {
              posvarsD.push_back(posvd);
              posvarsR.push_back(posvr);
            }
          }
#if PY_VERSION_HEX >= 0x03000000
          else if (PyUnicode_Check(tpl0))
          {
            const char* varname = PyUnicode_AsUTF8(tpl0);
            posvd = K_ARRAY::isNamePresent(varname, varStringD);
            posvr = K_ARRAY::isNamePresent(varname, varStringR);
            if (posvd == poscd) posvarcd = posvd;
            if (posvr == poscr) posvarcr = posvr;
            if (posvd != -1 && posvr != -1)
            {
              posvarsD.push_back(posvd);
              posvarsR.push_back(posvr);
            }
          }
#endif
          else
            PyErr_Warn(PyExc_Warning, "_setInterpTransfers: variable must be a string. Skipped.");
        }
      }
      else initAll = true;
    } // ckeck list of variables non empty
    else { initAll = true; }

    E_Int posvx=-1, posvy=-1, posvz=-1, posmx=-1, posmy=-1, posmz=-1;
    if ( dirR > 0 )
    {
      posvx = K_ARRAY::isNamePresent("VelocityX", varStringR);
      posvy = K_ARRAY::isNamePresent("VelocityY", varStringR);
      posvz = K_ARRAY::isNamePresent("VelocityZ", varStringR);
      posmx = K_ARRAY::isNamePresent("MomentumX", varStringR);
      posmy = K_ARRAY::isNamePresent("MomentumY", varStringR);
      posmz = K_ARRAY::isNamePresent("MomentumZ", varStringR);
    }
    if (initAll == true)// all common variables are transfered
    {
      char* varStringC; // chaine de caractere commune
      E_Int l = strlen(varStringR);
      varStringC = new char [l+1];
      // les positions demarrent a 1
      K_ARRAY::getPosition(varStringR, varStringD,
                   posvarsR, posvarsD, varStringC);
      delete [] varStringC;
      E_Int sizeVarsD = posvarsD.size();
      E_Int sizeVarsR = posvarsR.size();

      for (E_Int i = 0; i < sizeVarsD; i++) posvarsD[i] -= 1;
      for (E_Int i = 0; i < sizeVarsR; i++) posvarsR[i] -= 1;
      posvarcr = poscr; posvarcd = poscd;// position de cellNVariable dans la liste des variables a interpoler
    }
    if (posvarcd > -1 && posvarcr > -1) // cellNVariable exists : do not interpolate but specific update
    {
      posvarsD.erase(remove(posvarsD.begin(), posvarsD.end(), posvarcd), posvarsD.end());
      posvarsR.erase(remove(posvarsR.begin(), posvarsR.end(), posvarcr), posvarsR.end());
    }
    delete [] varStringR; delete [] varStringD; delete [] eltTypeR; delete [] eltTypeD;

    // -- no check (perfo) --
    // Transferts
    // Types valides: 0, 1, 2, 3, 4, 5
    nvars = posvarsR.size();

    vector<E_Float*> vectOfRcvFields(nvars);
    vector<E_Float*> vectOfDnrFields(nvars);

    for (E_Int eq = 0; eq < nvars; eq++)
    {
      vectOfRcvFields[eq] = fieldsR[posvarsR[eq]];
      vectOfDnrFields[eq] = fieldsD[posvarsD[eq]];
    }
    // interpolation of all fields
# include "commonInterpTransfers_indirect.h"

    // transfer of cellN variable
    if (posvarcr > -1 && posvarcd > -1) // cellNVariable exists and is transfered specifically
    {
      E_Int indR, type, nocf;
      E_Int indD0, indD, i, j, k, ncfLoc;
      E_Int noi = 0; // compteur sur le tableau d'indices donneur
      E_Int sizecoefs = 0;
      E_Float* cellNR = fieldsR[poscr];
      E_Float* cellND = fieldsD[poscd];
      for (E_Int noind = 0; noind < nbRcvPts; noind++)
      {
        // adressage indirect pour indR
        indR = rcvPts[noind];
# include "commonCellNTransfersStrict.h"
      }
    }
      // Prise en compte de la periodicite par rotation
    if (dirR != 0)
    {
# include "includeTransfers.h"
    }
  }// end of  compact = 0
  else  // compacted fields
  {
    // les variables a transferer sont compactees: on recupere uniquement la premiere et la taille
#include "getfromzonecompact.h"
    if( vartype <= 3 &&  vartype >= 1) nvars =5;
    else                               nvars =6;

    vector<E_Float*> vectOfRcvFields(nvars);
    vector<E_Float*> vectOfDnrFields(nvars);

    for (E_Int eq = 0; eq < nvars; eq++)
    {
      vectOfRcvFields[eq] = iptroR + eq*ndimdxR;
      vectOfDnrFields[eq] = iptroD + eq*ndimdxD;
    }
    //Parallel interpolation of all fields
# include "commonInterpTransfers_indirect.h"

    // Prise en compte de la periodicite par rotation
    if ( dirR != 0 )
    {
      printf("Warning: _setInterpTransfers: no correction for periodicity by rotation can be applied in compacted version.\n");
      //# include "includeTransfers.h"
    }
  }// end of compacted field transfers

  // sortie
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
  RELEASESHAREDN(pyIndRcv    , rcvPtsI    );
  RELEASESHAREDN(pyIndDonor  , donorPtsI  );
  RELEASESHAREDN(pyArrayTypes, typesI     );
  RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
  Py_INCREF(Py_None);
  return Py_None;
}
//=============================================================================
// Idem: in place + from zone + tc compact au niveau zone donneuse
//=============================================================================
PyObject* K_CONNECTOR::__setInterpTransfers(PyObject* self, PyObject* args)
{
  PyObject *zonesR, *zoneD;
  PyObject *pyVariables;
  PyObject *pyParam_int, *pyParam_real;
  E_Int loc, vartype, compact, flagibc, bctype;
  E_Float gamma, cv, muS, Cs, Ts;

  if (!PYPARSETUPLE_(args, OOOO_ O_ IIII_ RRRR_ R_,
                    &zonesR, &zoneD, &pyVariables, &pyParam_int,  &pyParam_real, &vartype, &compact,
                    &flagibc, &bctype, &gamma, &cv, &muS, &Cs, &Ts))
  {
    return NULL;
  }

  /* varType :
     1  : conservatives,
     11 : conservatives + ronutildeSA
     2  : (ro,u,v,w,t)
     21 : (ro,u,v,w,t) + ronutildeSA
     3  : (ro,u,v,w,p)
     31 : (ro,u,v,w,p) + ronutildeSA */
  E_Int varType = E_Int(vartype); 
  E_Int flagIbc = E_Int(flagibc);
  E_Int ibcType =  E_Int(bctype);
  vector<PyArrayObject*> hook;
  E_Int imdjmd, imd,jmd,kmd, cnNfldD, nvars,ndimdxR, ndimdxD,meshtype;
  E_Float* iptroD; E_Float* iptroR;

  E_Int nidom = PyList_Size(zonesR);

  E_Int* ipt_ndimdxR; E_Float** ipt_roR;
  ipt_ndimdxR = new E_Int[nidom];
  ipt_roR     = new E_Float*[nidom];

  E_Float Pr = 0.71;
  PyObject* zone0 = PyList_GetItem(zonesR, 0);
  PyObject* own   = K_PYTREE::getNodeFromName1(zone0 , ".Solver#ownData");
  if (own != NULL)
  {
    PyObject* paramreal0 = K_PYTREE::getNodeFromName1(own, "Parameter_real");
    if (paramreal0 != NULL)
    {
      E_Float* paramreal0val = K_PYTREE::getValueAF(paramreal0, hook);
      Pr = paramreal0val[10];
    }
  }

  /*-------------------------------------*/
  /* Extraction tableau int et real      */
  /*-------------------------------------*/
  FldArrayI* param_int;
  E_Int res_donor = K_NUMPY::getFromNumpyArray(pyParam_int, param_int, true);
  E_Int* ipt_param_int = param_int->begin();
  FldArrayF* param_real;
  res_donor = K_NUMPY::getFromNumpyArray(pyParam_real, param_real, true);
  E_Float* ipt_param_real = param_real->begin();

  E_Int nracID = ipt_param_int[0];
  E_Int nracIBC= ipt_param_int[1];
  E_Int nractot= nracID+nracIBC;

  E_Int nrac, shift_rac;
  if (flagIbc == 0) {nrac = nracID ; shift_rac = 2       ;}
  else              {nrac = nracIBC; shift_rac = 2+nracID;}

  vector<E_Float*> fieldsR;vector<E_Float*> fieldsD;
  vector<E_Int> posvarsD; vector<E_Int> posvarsR;
  E_Int* ptrcnd;
  //char* eltTypeR; char* eltTypeD;


  //Loc restrictif, tous les zones receveuses doivent avoir la meme loc
  loc = ipt_param_int[ 1+2*nractot + 6 ];
# include "getfromzoneDcompact.h"

  if( vartype <= 3 &&  vartype >= 1) nvars =5;
  else                               nvars =6;

  for  (E_Int irac=0; irac< nrac; irac++)
  { 
    E_Int pos_rac   = ipt_param_int[ shift_rac +irac];
    E_Int no_zR     = ipt_param_int[ pos_rac   +  6 ];
    //printf("pos_rac %d %d %d %d \n", pos_rac, no_zR, shift_rac+irac, irac);
    PyObject* zoneR = PyList_GetItem(zonesR, no_zR); // domaine i
#     include "getfromzoneRcompact.h"
    ipt_ndimdxR[irac] = ndimdxR;
    ipt_roR[irac]     = iptroR;
  }

  E_Float** RcvFields = new E_Float*[ nvars];
  E_Float** DnrFields = new E_Float*[ nvars];

  E_Float** vectOfRcvFields = RcvFields;
  E_Float** vectOfDnrFields = DnrFields;

  for  (E_Int irac=0; irac< nrac; irac++)
  {
    for (E_Int eq = 0; eq < nvars; eq++)
    {
      vectOfRcvFields[eq] = ipt_roR[irac] + eq*ipt_ndimdxR[irac];
      vectOfDnrFields[eq] = iptroD        + eq*ndimdxD;
    }

    ////
    //  Interpolation parallele
    ////
    //# include "commonInterpTransfers_indirect.h"
    ////
    imdjmd = imd*jmd;
    E_Int max_thread = min(nvars , E_Int(__NUMTHREADS__));

# pragma omp parallel default(shared) num_threads(max_thread)
    {
#ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num()+1;
    E_Int  Nbre_thread_actif = omp_get_num_threads(); // nombre de thread actif dans cette zone
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
#endif
    // Calcul du nombre de champs a traiter par chaque thread
    E_Int chunk = nvars/Nbre_thread_actif;
    E_Int r = nvars - chunk*Nbre_thread_actif;
    E_Int eq_deb, eq_fin;
    // equations traitees par thread
    if (ithread <= r) { eq_deb = (ithread-1)*(chunk+1)          ; eq_fin = eq_deb + (chunk+1); }
    else              { eq_deb = (chunk+1)*r+(ithread-r-1)*chunk; eq_fin = eq_deb +  chunk;    }

    E_Int pos         = ipt_param_int[ shift_rac + nractot + irac];
    E_Float* ptrCoefs = ipt_param_real + pos;

    pos            = ipt_param_int[ shift_rac + irac];
    E_Int nbRcvPts = ipt_param_int[  pos + 1        ];
    E_Int nbDonPts = ipt_param_int[  pos            ];

    E_Int* types    = ipt_param_int +  pos + 7 + nbRcvPts + nbDonPts;
    E_Int* donorPts = ipt_param_int +  pos + 7             ;
    E_Int* rcvPts   = ipt_param_int +  pos + 7 +   nbDonPts;// donor et receveur inverser car storage donor

    //printf("nbRcvPts %d %d %d %d %d \n", nbRcvPts , types[0], irac, pos, pos + 7 + nbRcvPts + nbDonPts );

    E_Int indR, type;
    E_Int indD0, indD, i, j, k, ncfLoc;
    E_Int noi = 0; // compteur sur le tableau d indices donneur
    E_Int sizecoefs = 0;

    indR = rcvPts[0];
    //printf(" indR00 %d %d %d %d %d \n",nbRcvPts , nbRcvPts, nbDonPts, imd, jmd);
    //printf(" ndimdxR %d  %d  \n", ipt_ndimdxR[irac],ndimdxD );

    for (E_Int noind = 0; noind < nbRcvPts; noind++)
    {
      //
      // adressage indirect pour indR
      //
      indR = rcvPts[noind];
#      include "commonInterpTransfers.h"
      ptrCoefs += sizecoefs;
    }
    }// omp

    if (flagIbc>0)
    {
      E_Int threadmax_sdm  = __NUMTHREADS__;
      E_Int pos      = ipt_param_int[ shift_rac + irac];
      E_Int nbRcvPts = ipt_param_int[  pos + 1        ];
      //E_Int nbDonPts = ipt_param_int[  pos            ];
      E_Int nbInterpD= ipt_param_int[  pos + 2        ];

      E_Int* rcvPts  = ipt_param_int +  pos + 7 +   nbRcvPts;// donor et receveur inverser car storage donor

      pos         = ipt_param_int[ shift_rac + nractot + irac];
      E_Float* ptrCoefs = ipt_param_real + pos;

      E_Int size = (nbRcvPts/threadmax_sdm)+1; // on prend du gras pour gerer le residus
      E_Int    r =  size % 8;
      if (r != 0) size  = size + 8 - r;        // on rajoute du bas pour alignememnt 64bits
      if (bctype <=1 ) size = 0;               // tableau inutile

      FldArrayF  tmp(size*17*threadmax_sdm);
      E_Float* ipt_tmp=  tmp.begin();

      E_Float* xPC     = ptrCoefs + nbInterpD;
      E_Float* xPI     = ptrCoefs + nbInterpD +3*nbRcvPts;
      E_Float* xPW     = ptrCoefs + nbInterpD +6*nbRcvPts;
      E_Float* densPtr = ptrCoefs + nbInterpD +9*nbRcvPts;

      E_Float param_real[30]; 
      param_real[ GAMMA] = gamma;
      param_real[ CVINF] = cv;
      param_real[ XMUL0] = muS;
      param_real[ CS] = Cs;
      param_real[ TEMP0] = Ts;
      param_real[ PRANDT] = Pr;

#     pragma omp parallel default(shared)
      {
        //indice loop pour paralelisation omp
        E_Int ideb, ifin;
#ifdef _OPENMP
        E_Int  ithread           = omp_get_thread_num()+1;
        E_Int  Nbre_thread_actif = omp_get_num_threads(); // nombre de thread actif dans cette zone
#else
        E_Int ithread = 1;
        E_Int Nbre_thread_actif = 1;
#endif

        // Calcul du nombre de champs a traiter par chaque thread
        E_Int chunk = nbRcvPts/Nbre_thread_actif;
        E_Int r = nbRcvPts - chunk*Nbre_thread_actif;
        // pts traitees par thread
        if (ithread <= r)
        { ideb = (ithread-1)*(chunk+1); ifin = ideb + (chunk+1); }
        else { ideb = (chunk+1)*r+(ithread-r-1)*chunk; ifin = ideb + chunk; }
        //  creer 2 zone  para pour threder au Max les loi de paroi
        if (varType == 2 || varType == 21)
          setIBCTransfersCommonVar2(ibcType, rcvPts, nbRcvPts, ideb, ifin, ithread,
                    xPC, xPC+nbRcvPts, xPC+nbRcvPts*2,
                    xPW, xPW+nbRcvPts, xPW+nbRcvPts*2,
                    xPI, xPI+nbRcvPts, xPI+nbRcvPts*2,
                    densPtr,
                    ipt_tmp, size, nvars,
                    param_real,
                    vectOfDnrFields, vectOfRcvFields);
        else {printf("Warning: setInterpTransfers: varType must be 2 or 21 \n");}

      } // Fin zone // omp
    } //ibc

    }//irac


  delete [] ipt_ndimdxR; delete [] ipt_roR;
  delete [] RcvFields; delete [] DnrFields;

  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
  RELEASESHAREDN(pyParam_int    , param_int    );
  RELEASESHAREDN(pyParam_real   , param_real   );

  Py_INCREF(Py_None);
  return Py_None;
}
//=============================================================================
// Idem: in place + from zone + tc compact au niveau base. Valid pour FastS uniquememnt
//=============================================================================
PyObject* K_CONNECTOR::___setInterpTransfers(PyObject* self, PyObject* args)
{
  PyObject *zonesR, *zonesD;

  PyObject *pyVariables, *pydtloc;
  PyObject *pyParam_int, *pyParam_real;
  E_Int vartype, type_transfert, no_transfert, It_target, nstep, nitmax;
  E_Int rk, exploc, num_passage, isWireModel;
  
  if (!PYPARSETUPLE_(args, OOOO_ OO_ IIII_ IIII_ II_,  
                    &zonesR, &zonesD, &pyVariables, &pydtloc, &pyParam_int,  
                    &pyParam_real, &It_target, &vartype,
                    &type_transfert, &no_transfert, &nstep, &nitmax, &rk, 
                    &exploc, &num_passage, &isWireModel))
  {
    return NULL;
  }
  E_Int rank     = -100;
  E_Int it_target=  E_Int(It_target);
  /* varType :
     1  : conservatives,
     11 : conservatives + ronutildeSA
     2  : (ro,u,v,w,t)
     21 : (ro,u,v,w,t) + ronutildeSA
     22 : (ro,u,v,w,t) + ronutildeSA + (gradxRo, gradyRo, gradzRo) + (gradxT, gradyT, gradzT)
     3  : (ro,u,v,w,p)
     31 : (ro,u,v,w,p) + ronutildeSA
     4  : (ro,u,v,w,t) + (Q1,..., QN) LBM
     41 : (ro,u,v,w,t) + (Sxx,...) + (corr_xx,...) + (Q1,...,QN) LBM OVERSET
     ---------------------------------------------
     5  : (ro,u,v,w,t) + (Sxx,...) Couplage NS LBM
     51 : Couplage NS LBM improved ( a coder ) */

  E_Int varType = E_Int(vartype);

  //gestion nombre de pass pour ID et/ou IBC
  E_Int TypeTransfert = E_Int(type_transfert);
  E_Int pass_deb, pass_fin;
  if     (TypeTransfert==0) { pass_deb =1; pass_fin =2; }//ID
  else if(TypeTransfert==1) { pass_deb =0; pass_fin =1; }//IBCD
  else                      { pass_deb =0; pass_fin =2; }//ALL

  E_Int NoTransfert  = E_Int(no_transfert);

  E_Int kmd, cnNfldD, nvars, meshtype, nvars_Pnt2;

  if     ( vartype <= 3 &&  vartype >= 1) nvars =5;
  else if( vartype == 4 ) nvars =32; // LBM transfer, 19 or 27 Qs and 5 macros (32 max in total)
                                     // on majore pour la LBM, car nvar sert uniquememnt a dimensionner taille vector
  else if (vartype == 41) nvars =44; //38;    // LBM Overset : 19 or 27 Q, 5 macro and 6 gradients
  else if( vartype == 5 ) nvars =30;    // Hybrid NSLBM transfer, 5 macros, 19 Qs + 6 gradients
  else                    nvars =6;

  nvars_Pnt2 = 0;
  if (TypeTransfert==0)  isWireModel=0;
  if (isWireModel>0)     nvars_Pnt2=nvars; // *_WM variables (5 for NSLam & 6 for NSTurb) 

  E_Int nidomR   = PyList_Size(zonesR);
  E_Int nidomD   = PyList_Size(zonesD);

  //pointeur pour stocker solution au centre ET au noeud
  E_Int* ipt_ndimdxD; E_Int** ipt_param_intR; E_Int** ipt_cnd;
  E_Float** ipt_roR; E_Float** ipt_roD; E_Float** ipt_roR_vert;  E_Float** ipt_roD_vert; E_Float** ipt_param_realR;
  E_Float** ipt_roR_Pnt2;
  E_Float** ipt_roD_Pnt2;
  
  E_Float** ipt_qR; E_Float** ipt_qD; E_Float** ipt_qR_vert; E_Float** ipt_qD_vert;
  E_Float** ipt_SR; E_Float** ipt_SD; E_Float** ipt_SR_vert; E_Float** ipt_SD_vert;
  E_Float** ipt_psiGR; E_Float** ipt_psiGD; E_Float** ipt_psiGR_vert; E_Float** ipt_psiGD_vert;

  //ipt_ndimdxR      = new E_Int*[nidomR*3];   // on stocke ndimdx  en centre et vertexe
  ipt_param_intR   = new E_Int*[nidomR];

  ipt_roR            = new E_Float*[nidomR*10]; //1
  ipt_roR_vert       = ipt_roR         + nidomR;//2
  ipt_param_realR    = ipt_roR_vert    + nidomR;//3
  ipt_roR_Pnt2       = ipt_param_realR + nidomR;//4
  ipt_qR             = ipt_roR_Pnt2    + nidomR;//5
  ipt_qR_vert        = ipt_qR          + nidomR;//6
  ipt_SR             = ipt_qR_vert     + nidomR;//7
  ipt_SR_vert        = ipt_SR          + nidomR;//8
  ipt_psiGR          = ipt_SR_vert     + nidomR;//9
  ipt_psiGR_vert     = ipt_psiGR       + nidomR;//10

  ipt_ndimdxD        = new E_Int[nidomD*8];  //on stocke ndimdx, imd, jmd, en centre et vertexe, meshtype et cnDfld
  ipt_cnd            = new E_Int*[nidomD];

  ipt_roD            = new E_Float*[nidomD*9];  //1
  ipt_qD             = ipt_roD         + nidomD;//2
  ipt_SD             = ipt_qD          + nidomD;//3
  ipt_psiGD          = ipt_SD          + nidomD;//4
  ipt_roD_vert       = ipt_psiGD       + nidomD;//5
  ipt_qD_vert        = ipt_roD_vert    + nidomD;//6
  ipt_SD_vert        = ipt_qD_vert     + nidomD;//7
  ipt_psiGD_vert     = ipt_SD_vert     + nidomD;//8
  ipt_roD_Pnt2       = ipt_psiGD_vert  + nidomD;//9

  vector<PyArrayObject*> hook;

  FldArrayI* dtloc;
  E_Int res_donor = K_NUMPY::getFromNumpyArray(pydtloc, dtloc, true);
  E_Int* iptdtloc = dtloc->begin();
  /*-------------------------------------*/
  /* Extraction tableau int et real de tc*/
  /*-------------------------------------*/
  FldArrayI* param_int;
  res_donor = K_NUMPY::getFromNumpyArray(pyParam_int, param_int, true);
  E_Int* ipt_param_int = param_int->begin();
  FldArrayF* param_real;
  res_donor = K_NUMPY::getFromNumpyArray(pyParam_real, param_real, true);
  E_Float* ipt_param_real = param_real->begin();

  /*------------------------------------------------*/
  /* RECUPERATION DU NOM DES VARIABLES A TRANSFERER */
  /*------------------------------------------------*/
  char* varname = NULL; char* varname1 = NULL; char* varname2 = NULL; char* varname3 = NULL;
  char* vartmp  = NULL; 
  E_Int nbvar_inlist   = PyList_Size(pyVariables);

  for (E_Int ivar = 0; ivar < nbvar_inlist; ivar++)
  {
    PyObject* tpl0= PyList_GetItem(pyVariables, ivar);
    if (PyString_Check(tpl0)) vartmp = PyString_AsString(tpl0);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(tpl0)) vartmp = (char*)PyUnicode_AsUTF8(tpl0);
#endif
    // printf("varname %s \n", vartmp);
    if     (ivar==0) {varname  = vartmp;}
    else if(ivar==1) {varname1 = vartmp;}  //En LBM, on a besoin d'echanger les macro !!ET!! les Q donc deux varname
    else if(ivar==2) {varname2 = vartmp;}  //Pour le couplage NS-LBM, on a besoin d'echanger les Q !!ET!! les macro !!ET!! les gradients
    else if(ivar==3) {varname3 = vartmp;}
    else {printf("Warning: souci varname setInterpTransfers \n"); }
  }


  //on recupere sol et solcenter ainsi que connectivite et taille zones Donneuses (tc)
  for (E_Int nd = 0; nd < nidomD; nd++)
  {
    PyObject* zoneD = PyList_GetItem(zonesD, nd);
#    include "getfromzoneDcompact_all.h"
  }
  //on recupere sol et solcenter taille zones receuveuses, param_int et param_real (t)
  for (E_Int nd = 0; nd < nidomR; nd++)
  {
    PyObject* zoneR = PyList_GetItem(zonesR, nd);
#     include "getfromzoneRcompact_all.h"
  }


  E_Int nbcomIBC = ipt_param_int[2];
  E_Int nbcomID  = ipt_param_int[3+nbcomIBC];

  E_Int shift_graph = nbcomIBC + nbcomID + 3;

  E_Int threadmax_sdm  = __NUMTHREADS__;
  E_Int ech            = ipt_param_int[ NoTransfert +shift_graph];
  E_Int nrac           = ipt_param_int[ ech +1 ];          //nb total de raccord
  E_Int nrac_inst      = ipt_param_int[ ech +2 ];          //nb total de raccord instationnaire
  E_Int timelevel      = ipt_param_int[ ech +3 ];          //nb de pas de temps stocker pour chaque raccord instationnaire
  E_Int nrac_steady    = nrac - nrac_inst;                 //nb total de raccord stationnaire

  //gestion nombre de pass pour raccord instationnaire
  E_Int pass_inst_deb=0;
  E_Int pass_inst_fin=1;
  if (nrac_inst > 0) pass_inst_fin=2;

  //
  //on optimise les transfert pour implicit local
  //
  E_Int impli_local[nidomR];
  E_Int nssiter  = iptdtloc[0];
  E_Int shift_omp= iptdtloc[11];
  E_Int* ipt_omp = iptdtloc + shift_omp;
  E_Int nbtask   = ipt_omp[nstep-1]; 
  E_Int ptiter   = ipt_omp[nssiter+ nstep-1];

  for (E_Int nd = 0; nd < nidomR; nd++) {impli_local[nd]=0;}//par defaut pas de transfert
  for (E_Int ntask = 0; ntask < nbtask; ntask++)            //transfert sur les zones modifiees � la ssiter nstep
  {
    E_Int pttask = ptiter + ntask*(6+threadmax_sdm*7);
    E_Int nd = ipt_omp[ pttask ];
    impli_local[nd] = 1;
  }
  E_Int maxlevel      =  iptdtloc[ 9];  //transfert sur les zones qui recupere leur valeur interpolees en LBM
  E_Int it_cycl_lbm   =  iptdtloc[10];
  E_Int level_it      =  iptdtloc[12+it_cycl_lbm];
  E_Int max_it        = pow(2, maxlevel-1);

  E_Int level_next_it =  maxlevel;
  if (it_cycl_lbm != max_it -1 ) { level_next_it = iptdtloc[12 +it_cycl_lbm +1];}

  for (E_Int nd = 0; nd < nidomR; nd++)   
    {
       if (  ipt_param_intR[nd][IFLOW] == 4)
       {
          E_Int level = ipt_param_intR[nd][LEVEL];
          if (level <=level_next_it +1 ){impli_local[nd]=1;}
          else                          {impli_local[nd]=0;}
       }
       //printf("implilocal %d %d %d \n", impli_local[nd], nd, it_cycl_lbm);
    }

  E_Int size_autorisation = nrac_steady+1;
  size_autorisation = K_FUNC::E_max(size_autorisation , nrac_inst+1);
  E_Int autorisation_transferts[pass_inst_fin][size_autorisation];

  // printf("nrac = %d, nrac_inst = %d, level= %d, it_target= %d , nitrun= %d \n",  nrac, nrac_inst, timelevel,it_target, NitRun);
  //on dimension tableau travail pour IBC
  E_Int nbRcvPts_mx =0; E_Int ibcTypeMax=0;
  for  (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
  {
    E_Int irac_deb= 0; E_Int irac_fin= nrac_steady;
    if (pass_inst == 1) { irac_deb = ipt_param_int[ ech + 4 + it_target ]; irac_fin = ipt_param_int[ ech + 4 + it_target + timelevel ];}


    for  (E_Int irac=irac_deb; irac< irac_fin; irac++)
    {
      E_Int shift_rac =  ech + 4 + timelevel*2 + irac;
      if (ipt_param_int[ shift_rac+ nrac*10 + 1] > nbRcvPts_mx) nbRcvPts_mx = ipt_param_int[ shift_rac+ nrac*10 + 1];

      if (ipt_param_int[shift_rac+nrac*3] > ibcTypeMax) ibcTypeMax =  ipt_param_int[shift_rac+nrac*3];

      E_Int irac_auto= irac-irac_deb;
      autorisation_transferts[pass_inst][irac_auto]=0;

      if (exploc == 1)  //if(rk==3 && exploc == 2) // Si on est en explicit local, on va autoriser les transferts entre certaines zones seulement en fonction de la ss-ite courante
      {
        E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;
        E_Int levelD = ipt_param_int[debut_rac + 25];
        E_Int levelR = ipt_param_int[debut_rac + 24];
        E_Int cyclD  = nitmax/levelD;

        // Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse
        if (levelD > levelR && num_passage == 1)
        {
          if ( nstep%cyclD==cyclD-1 || (nstep%cyclD==cyclD/2 && (nstep/cyclD)%2==1) )
          {
            autorisation_transferts[pass_inst][irac_auto]=1;
          }
          else {continue;}
        }
        // Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
        else if (levelD < levelR && num_passage == 1)
        {
          if (nstep%cyclD==1 || nstep%cyclD==cyclD/4 || nstep%cyclD== cyclD/2-1 || nstep%cyclD== cyclD/2+1 || nstep%cyclD== cyclD/2+cyclD/4 || nstep%cyclD== cyclD-1)
            { autorisation_transferts[pass_inst][irac_auto]=1; }
          else {continue;}
        }
        // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse
        else if (levelD == levelR && num_passage == 1)
        {
          if (nstep%cyclD==cyclD/2-1 || (nstep%cyclD==cyclD/2 && (nstep/cyclD)%2==0) || nstep%cyclD==cyclD-1)
            { autorisation_transferts[pass_inst][irac_auto]=1; }
          else {continue;}
        }
        // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse (cas du deuxieme passage)
        else if (levelD == ipt_param_int[debut_rac +24] && num_passage == 2)
        {
          if (nstep%cyclD==cyclD/2 && (nstep/cyclD)%2==1)
            { autorisation_transferts[pass_inst][irac_auto]=1; }
          else {continue;}
        }
        else {continue;}
      }
      // Sinon, on autorise les transferts  si la zone donneuse a ete modifiee a l'iteration nstep
      else 
      {
        E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ];
        if (impli_local[NoD]==1) autorisation_transferts[pass_inst][irac_auto]=1;
        //autorisation_transferts[pass_inst][irac_auto]=1;
      }
        
    }
    }

  E_Int size = (nbRcvPts_mx/threadmax_sdm)+1; // on prend du gras pour gerer le residus
  E_Int r =  size % 8;
  if (r != 0) size  = size + 8 - r;           // on rajoute du bas pour alignememnt 64bits
  if (ibcTypeMax <=1 ) size = 0;              // tableau inutile : SP voir avec Ivan

  FldArrayF  tmp(size*17*threadmax_sdm);
  E_Float* ipt_tmp = tmp.begin();
  
  E_Float** RcvFields = new E_Float*[ (nvars+nvars_Pnt2)*threadmax_sdm];
  E_Float** DnrFields = new E_Float*[              nvars*threadmax_sdm];

  //# pragma omp parallel default(shared)  num_threads(1)
# pragma omp parallel default(shared)
  {
      
    E_Float gamma, cv, muS, Cs, Ts, Pr;

#ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num()+1;
    E_Int  Nbre_thread_actif = omp_get_num_threads(); // nombre de thread actif dans cette zone
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
#endif

    E_Int indR, type;
    E_Int indD0, indD, i, j, k, ncfLoc/*, nocf*/, indCoef, noi, sizecoefs, /*Nbchunk,*/ imd, jmd, imdjmd;

    E_Float** vectOfRcvFields = RcvFields + (nvars+nvars_Pnt2)*(ithread-1);
    E_Float** vectOfDnrFields = DnrFields +  nvars*(ithread-1);
    
    //1ere pass_typ: IBC
    //2eme pass_typ: transfert

    for  (E_Int ipass_typ=pass_deb; ipass_typ< pass_fin; ipass_typ++)
    {
      //1ere pass_inst: les raccords fixes
      //2eme pass_inst: les raccords instationnaires
      for  (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
      {
        //printf("pass_inst = %d, level= %d \n",  pass_inst, nrac_inst_level );
        E_Int irac_deb= 0; E_Int irac_fin= nrac_steady;
        if(pass_inst == 1)
        {
          irac_deb = ipt_param_int[ ech + 4 + it_target             ];
          irac_fin = ipt_param_int[ ech + 4 + it_target + timelevel ];
        }

////# pragma omp for schedule(dynamic)
        for  (E_Int irac=irac_deb; irac< irac_fin; irac++)
        {
          E_Int irac_auto= irac-irac_deb;
          if (autorisation_transferts[pass_inst][irac_auto]==1)
          {
            E_Int shift_rac =  ech + 4 + timelevel*2 + irac;
            //printf("ipass_typ = %d, pass_inst= %d, irac=  %d, ithread= %d \n", ipass_typ,pass_inst,irac , ithread );
            // ibc: -1 : raccord ID, 0 : wallslip, 1 slip, 2 : log , 3 : Musker, 4 : outpress
            E_Int ibcType =  ipt_param_int[shift_rac+nrac*3];
            E_Int ibc = 1;
            if (ibcType < 0) ibc = 0;
            if(1-ibc != ipass_typ)  continue;

            E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ];
            E_Int loc      =  ipt_param_int[ shift_rac + nrac*9  +1 ]; //+1 a cause du nrac mpi
            E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ];
            E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 +1 ]; //neq fonction raccord rans/LES
            E_Int rotation =  ipt_param_int[ shift_rac + nrac*14 +1 ]; //flag pour periodicite azimutale

            // COUPLAGE NS LBM - Recupere les solveurs des zones R et D
            E_Int solver_D=2; E_Int solver_R=2;
            if (nvars_loc == 11) {solver_R =4;}
            if (nvars_loc == -5) {solver_D =4; nvars_loc = 5;}
            if (nvars_loc == 19) {solver_D =4; solver_R=4;}
            if (nvars_loc == 27) {solver_D =4; solver_R=4;}

            E_Int overset  =  ipt_param_intR[NoD][LBM_OVERSET];        //flag pour overset en LBM
            if      (nvars_loc==19 && overset==0) nvars_loc = nvars_loc + 5;
            else if (nvars_loc==27 && overset==0) nvars_loc = nvars_loc + 5;
            else if (nvars_loc==19 && overset==1) nvars_loc = nvars_loc + 5 + 6 + 6;
            else if (nvars_loc==27 && overset==1) nvars_loc = nvars_loc + 5 + 6 + 6;
            // cout << nvars_loc << endl;
            // printf("irac=  %d, nvar/nvar_loc= %d %d,  ithread= %d  solverDR %d %d \n",irac , nvars,nvars_loc, ithread,  solver_D, solver_R);

            E_Int levelD   = ipt_param_intR[NoD][LEVEL];
            E_Int levelR   = ipt_param_intR[NoR][LEVEL];

            E_Int meshtype = ipt_ndimdxD[NoD + nidomD*6];
            E_Int cnNfldD  = ipt_ndimdxD[NoD + nidomD*7];
            E_Int* ptrcnd  = ipt_cnd[    NoD           ];

            if (loc == 0)
            {
              printf("Error: transferts optimises non code en vextex " SF_D3_ "\n", shift_rac + nrac*9  +1, NoD, NoR );
              //imd= ipt_ndimdxD[ NoD+ nidomD*4]; jmd= ipt_ndimdxD[ NoD + nidomD*5];
              imd = 0; jmd = 0;
            }
            else
            {
                        /*--------------------------------------------------------------------*/
                        /*                GESTION DES TRANSFERTS EN CENTER                    */
                        /* 2 cas: - si pas de couplage NSLBM, transferts habituels,           */
                        /*        - si couplage NSLBM, raccords a adapter                     */
                        /*--------------------------------------------------------------------*/
                        if (nvars_loc == 5 || nvars_loc == 6) // Transferts NS classiques ou LBM -> NS
                        {
                         for (E_Int eq = 0; eq < nvars_loc; eq++)
                         {
                           vectOfRcvFields[eq] = ipt_roR[ NoR] + eq*ipt_param_intR[ NoR][ NDIMDX ];
                           vectOfDnrFields[eq] = ipt_roD[ NoD] + eq*ipt_param_intR[ NoD][ NDIMDX ];
                         }
                        }
                        else if (nvars_loc == 24 || nvars_loc == 32) // Transferts LBM classiques
                        {
                          // On commence par copier les 5 variables macros
                          for (E_Int eq = 0; eq < 5; eq++)
                          {
                            vectOfRcvFields[eq] = ipt_roR[ NoR] + eq*ipt_param_intR[ NoR][ NDIMDX ];
                            vectOfDnrFields[eq] = ipt_roD[ NoD] + eq*ipt_param_intR[ NoD][ NDIMDX ];
                          }
                          // Puis on copie les fonctions de distribution
                          for (E_Int eq = 5; eq < nvars_loc; eq++)
                          {
                            vectOfRcvFields[eq] = ipt_qR[ NoR] + (eq-5)*ipt_param_intR[ NoR][ NDIMDX ];
                            vectOfDnrFields[eq] = ipt_qD[ NoD] + (eq-5)*ipt_param_intR[ NoD][ NDIMDX ];
                          }
                        }
                        else if (nvars_loc == 11 ) // //Transfert NS -> LBM    
                        {
                          // On commence par copier les 5 variables macros
                          for (E_Int eq = 0; eq < 5; eq++)
                          {
                            vectOfRcvFields[eq] = ipt_roR[ NoR] + eq*ipt_param_intR[ NoR][ NDIMDX ];
                            vectOfDnrFields[eq] = ipt_roD[ NoD] + eq*ipt_param_intR[ NoD][ NDIMDX ];
                          }
                          // Puis on copie les gradients
                          for (E_Int eq = 5; eq < nvars_loc; eq++)
                          {
                            vectOfRcvFields[eq] = ipt_SR[ NoR] + (eq-5)*ipt_param_intR[ NoR][ NDIMDX ];
                            vectOfDnrFields[eq] = ipt_SD[ NoD] + (eq-5)*ipt_param_intR[ NoD][ NDIMDX ];
                          }
                        }
                        else if (nvars_loc == 36 || nvars_loc == 44) // //Transfert LBM  overset   
                        {
                          // On commence par copier les 5 variables macros
                          for (E_Int eq = 0; eq < 5; eq++)
                          {
                            vectOfRcvFields[eq] = ipt_roR[ NoR] + eq*ipt_param_intR[ NoR][ NDIMDX ];
                            vectOfDnrFields[eq] = ipt_roD[ NoD] + eq*ipt_param_intR[ NoD][ NDIMDX ];
                          }
                          // Puis on copie les gradients
                          for (E_Int eq = 5; eq < 11; eq++)
                          {
                            vectOfRcvFields[eq] = ipt_SR[ NoR] + (eq-5)*ipt_param_intR[ NoR][ NDIMDX ];
                            vectOfDnrFields[eq] = ipt_SD[ NoD] + (eq-5)*ipt_param_intR[ NoD][ NDIMDX ];
                          }
                          for (E_Int eq =11; eq < 17; eq++)
                          {
                            vectOfRcvFields[eq] = ipt_psiGR[ NoR] + (eq-11)*ipt_param_intR[ NoR][ NDIMDX ];
                            vectOfDnrFields[eq] = ipt_psiGD[ NoD] + (eq-11)*ipt_param_intR[ NoD][ NDIMDX ];
                          }
                          for (E_Int eq =17; eq < nvars_loc; eq++)
                          {
                            vectOfRcvFields[eq] = ipt_qR[ NoR] + (eq-17)*ipt_param_intR[ NoR][ NDIMDX ];
                            vectOfDnrFields[eq] = ipt_qD[ NoD] + (eq-17)*ipt_param_intR[ NoD][ NDIMDX ];
                          }
                        }
                        if (isWireModel>0)
                        {
                          for (E_Int eq = nvars_loc; eq < nvars_loc+nvars_Pnt2; eq++){
                          vectOfRcvFields[eq] = ipt_roR_Pnt2[ NoR] + (eq-nvars_loc)*ipt_param_intR[ NoR ][ NDIMDX ];
                        }
                     }

                      imd= ipt_param_intR[ NoD ][ NIJK ]; jmd= ipt_param_intR[ NoD ][ NIJK+1];
            }

            imdjmd = imd*jmd;


            ////
            //  Interpolation parallele
            ////
            ////

            E_Int nbRcvPts = ipt_param_int[ shift_rac +  nrac*10 + 1 ];

            E_Int pos;
            pos  = ipt_param_int[ shift_rac + nrac*7 ]     ; E_Int* ntype      = ipt_param_int  + pos;
            pos  = pos +1 + ntype[0]                       ; E_Int* types      = ipt_param_int  + pos;
            pos  = ipt_param_int[ shift_rac + nrac*6      ]; E_Int* donorPts   = ipt_param_int  + pos;
            pos  = ipt_param_int[ shift_rac + nrac*12 + 1 ]; E_Int* rcvPts     = ipt_param_int  + pos;   // donor et receveur inverser car storage donor
            pos  = ipt_param_int[ shift_rac + nrac*8      ]; E_Float* ptrCoefs = ipt_param_real + pos;

            E_Int nbInterpD = ipt_param_int[ shift_rac +  nrac ]; E_Float* xPC=NULL; E_Float* xPI=NULL; E_Float* xPW=NULL; E_Float* densPtr=NULL;
            if (ibc == 1)
            {
              xPC     = ptrCoefs + nbInterpD;
              xPI     = ptrCoefs + nbInterpD +3*nbRcvPts;
              xPW     = ptrCoefs + nbInterpD +6*nbRcvPts;
              densPtr = ptrCoefs + nbInterpD +9*nbRcvPts;
            }

            E_Int ideb      = 0;
            E_Int ifin      = 0;
            E_Int shiftCoef = 0;
            E_Int shiftDonor  = 0;
            
            for (E_Int ndtyp = 0; ndtyp < ntype[0]; ndtyp++)
            {
              type      = types[ifin];

              SIZECF(type, meshtype, sizecoefs);
              ifin =  ifin + ntype[ 1 + ndtyp];

              E_Int pt_deb, pt_fin;

              // Calcul du nombre de champs a traiter par chaque thread
              E_Int size_bc =  ifin-ideb;
              E_Int chunk   =  size_bc/Nbre_thread_actif;
              E_Int r       =  size_bc - chunk*Nbre_thread_actif;
              // pts traitees par thread
              if (ithread <= r)
              { pt_deb = ideb + (ithread-1)*(chunk+1);           pt_fin = pt_deb + (chunk+1); }
              else { pt_deb = ideb + (chunk+1)*r+(ithread-r-1)*chunk; pt_fin = pt_deb + chunk; }

              //Si type 0, calcul sequentiel
              if  ( type == 0 )
              { if (ithread ==1 ){ pt_deb = ideb; pt_fin = ifin;}
                else             { pt_deb = ideb; pt_fin = ideb;}
              }

              // ATTENTION
              //pt_deb = ideb; pt_fin = ifin;

              noi       = shiftDonor;                             // compteur sur le tableau d indices donneur
              indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
    
              E_Int linterp =1;
              E_Int shiftv  =0;
              if      ( isWireModel == 2 ) linterp = 0;
              else if ( isWireModel == 1 )
              {
                if ( ibcType!=141) linterp = 0;
                if ( ibcType==141) shiftv = nvars_loc;
              }
              else if ( isWireModel == 0 and  ibcType==141) linterp = 0;

            if ( (nvars_loc==5 || (ibc==1 && solver_R==4)) && linterp==1 )
            {
#           include "commonInterpTransfers_reorder_5eq.h"
            }
            else if (nvars_loc==6 and linterp== 1)
            {
#           include "commonInterpTransfers_reorder_6eq.h"
            }
            else if(nvars_loc==19)
            {
#           include "commonInterpTransfers_reorder_19eq.h" 
            }
            else
            {
              if (linterp==1) 
              {
#               include "commonInterpTransfers_reorder_neq.h"
              }
            }
            
            // =============================================================
            // RAFFINEMENT DE MAILLAGE LBM/LBM
            // =============================================================
            if (solver_D==4 && solver_R==4 && levelD > levelR)
            {
#             include "includeTransfers_coarse2fine_LBM.h"
            }
            else if (solver_D == 4 && solver_R == 4 && levelD < levelR)
            {
#             include "includeTransfers_fine2coarse_LBM.h"
            }
            // =============================================================

            // =============================================================
            // COUPLAGE NS-LBM: changement d'unite
            // =============================================================
            // Code que pour le cas D3Q19 pour le moment, a adapter au 27
            if (solver_D==4 && solver_R<4)
            {
              // Transfert LBM vers NS: repasse dans unites SI
#             include "includeTransfers_dimLBMtoNS.h"
            }
            else if (solver_D<4 && solver_R==4)
            {
              // Transfert NS vers LBM : adimensionnement
#             include "includeTransfers_dimNStoLBM.h"
            }

            // Prise en compte de la periodicite par rotation
            if (rotation == 1 and isWireModel!=2)
            {
              E_Float* angle = ptrCoefs + nbInterpD;
#          include "includeTransfers_rotation.h"
            }
            

            // ibc
            if (ibc == 1)
            {
                if (isWireModel ==0)
                {
                  setIBCTransfersCommonVar2(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread,
                                  xPC, xPC+nbRcvPts, xPC+nbRcvPts*2,
                                  xPW, xPW+nbRcvPts, xPW+nbRcvPts*2,
                                  xPI, xPI+nbRcvPts, xPI+nbRcvPts*2,
                                  densPtr, 
                                  ipt_tmp, size, nvars,
                                  ipt_param_realR[ NoR ],
                                  vectOfDnrFields, vectOfRcvFields);

                  if (solver_R==4)
                   {
#                   include "includeTransfers_LBM_feq.h"
                   }
                }
                else if(isWireModel==2)
                {
                  if (ibcType==140) // placing values in the tc
                  {
                    for (E_Int noind = 0; noind < pt_fin-pt_deb; noind++)
                    {
                      E_Int indR = rcvPts[noind+pt_deb];
                      (densPtr+nbRcvPts*5 )[noind+pt_deb] = vectOfRcvFields[nvars  ][indR];
                      (densPtr+nbRcvPts*6 )[noind+pt_deb] = vectOfRcvFields[nvars+1][indR];
                      (densPtr+nbRcvPts*7 )[noind+pt_deb] = vectOfRcvFields[nvars+2][indR];
                      (densPtr+nbRcvPts*8)[noind+pt_deb]  = vectOfRcvFields[nvars+3][indR];
                      (densPtr+nbRcvPts*9)[noind+pt_deb]  = vectOfRcvFields[nvars+4][indR];
                      if (nvars==6){ (densPtr+nbRcvPts*10)[noind+pt_deb] = vectOfRcvFields[nvars+5][indR]; }					  
                    }
                  }
                }
              }//ibc
            
            ideb       =  ideb + ntype[ 1 + ndtyp];
            shiftCoef  = shiftCoef  +  ntype[1+ndtyp]*sizecoefs; //shift coef   entre 2 types successif
            shiftDonor = shiftDonor +  ntype[1+ndtyp];           //shift donor entre 2 types successif
              }// type
#pragma omp barrier
          }// autorisation transfert
          }//irac
      }//pass_inst
#pragma omp barrier
      }//ipass
  }// omp


  delete [] ipt_param_intR; delete [] ipt_roR; delete [] ipt_ndimdxD; delete [] ipt_roD; delete [] ipt_cnd;
  delete [] RcvFields; delete [] DnrFields;

  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
  RELEASESHAREDN(pydtloc        , dtloc        );
  RELEASESHAREDN(pyParam_int    , param_int    );
  RELEASESHAREDN(pyParam_real   , param_real   );

  Py_INCREF(Py_None);

  return Py_None;
}

//=============================================================================
// Copy of ___setInterpTransfers() for gradP info
// transfers of ['Density', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature', 'nuSA'] -> 6 vars
// transfers of ['gradx/y/zDensity', 'gradx/y/zTemperature'] -> 6 varsGrad
// varType 22 : tc2/tc -> RCV ZONES
// varType 23 : RCV ZONES -> tc 
//=============================================================================
PyObject* K_CONNECTOR::___setInterpTransfers4GradP(PyObject* self, PyObject* args)
{
  PyObject *zonesR, *zonesD;

  PyObject *pyVariables;
  PyObject *pyParam_int, *pyParam_real;
  E_Int vartype, type_transfert, no_transfert, It_target, nstep, nitmax;
  E_Int rk, exploc, num_passage;
  E_Float gamma, cv, muS, Cs, Ts, Pr;

  if (!PYPARSETUPLE_(args, OOOO_ O_ IIII_ IIII_ I_,  
                    &zonesR, &zonesD, &pyVariables, &pyParam_int,  &pyParam_real, &It_target, &vartype,
                    &type_transfert, &no_transfert, &nstep, &nitmax, &rk, &exploc, &num_passage))
  {
    return NULL;
  }

  E_Int it_target = E_Int(It_target);

  E_Float alpha = 0.;

  /* varType :
     1  : conservatives,
     11 : conservatives + ronutildeSA
     2  : (ro,u,v,w,t)
     21 : (ro,u,v,w,t) + ronutildeSA
     22 : (ro,u,v,w,t) + ronutildeSA + (gradxRo, gradyRo, gradzRo) + (gradxT, gradyT, gradzT)
     3  : (ro,u,v,w,p)
     31 : (ro,u,v,w,p) + ronutildeSA
     4  : (Q1,..., QN)   LBM  */
  E_Int varType = E_Int(vartype);

  //gestion nombre de pass pour ID et/ou IBC
  E_Int TypeTransfert = E_Int(type_transfert);
  E_Int pass_deb, pass_fin;
  if     (TypeTransfert==0) { pass_deb =1; pass_fin =2; }//ID
  else if(TypeTransfert==1) { pass_deb =0; pass_fin =1; }//IBCD
  else                      { pass_deb =0; pass_fin =2; }//ALL

  E_Int NoTransfert  = E_Int(no_transfert);

  E_Int kmd, cnNfldD, nvars, meshtype, nvars_grad;

  if     ( vartype <= 3 &&  vartype >= 1) nvars =5;
  else if( vartype == 4 ) nvars =27;    // on majore pour la LBM, car nvar sert uniquememnt a dimensionner taille vector
  else                    nvars =6;

  if (vartype == 22 || vartype == 23 || vartype == 24)
  {
    nvars = 6;
    nvars_grad = 6;
  }

  E_Int nidomR   = PyList_Size(zonesR);
  E_Int nidomD   = PyList_Size(zonesD);

  //pointeur pour stocker solution au centre ET au noeud
  E_Int* ipt_ndimdxD; E_Int** ipt_param_intR; E_Int** ipt_cnd;
  E_Float** ipt_roR; E_Float** ipt_roD; E_Float** ipt_roR_vert;  E_Float** ipt_roD_vert; E_Float** ipt_param_realR;

  //*************************************
  E_Float** ipt_gradR;
  E_Float** ipt_gradD;
  //*************************************

  //ipt_ndimdxR      = new E_Int*[nidomR*3];   // on stocke ndimdx  en centre et vertexe

  ipt_param_intR   = new E_Int*[nidomR];

  // ipt_roR          = new E_Float*[nidomR*3];
  // ipt_roR_vert     = ipt_roR + nidomR;
  // ipt_param_realR  = ipt_roR_vert + nidomR;

  // ipt_ndimdxD      = new E_Int[nidomD*8];  //on stocke ndimdx, imd, jmd, en centre et vertexe, meshtype et cnDfld
  // ipt_cnd          = new E_Int*[nidomD];

  // ipt_roD          = new E_Float*[nidomD*2];
  // ipt_roD_vert     = ipt_roD + nidomD;


  vector<PyArrayObject*> hook;

  //*************************************
  ipt_roR          = new E_Float*[nidomR*4];
  ipt_roR_vert     = ipt_roR + nidomR;
  ipt_param_realR  = ipt_roR_vert + nidomR;
  ipt_gradR        = ipt_param_realR + nidomR;

  ipt_ndimdxD      = new E_Int[nidomD*8];  //on stocke ndimdx, imd, jmd, en centre et vertexe, meshtype et cnDfld
  ipt_cnd          = new E_Int*[nidomD];

  ipt_roD          = new E_Float*[nidomD*3];
  ipt_roD_vert     = ipt_roD + nidomD;
  ipt_gradD        = ipt_roD_vert + nidomD;
  //*************************************


  /*-------------------------------------*/
  /* Extraction tableau int et real de tc*/
  /*-------------------------------------*/
  FldArrayI* param_int;
  E_Int res_donor = K_NUMPY::getFromNumpyArray(pyParam_int, param_int, true);
  E_Int* ipt_param_int = param_int->begin();
  FldArrayF* param_real;
  res_donor = K_NUMPY::getFromNumpyArray(pyParam_real, param_real, true);
  E_Float* ipt_param_real = param_real->begin();

  //On recupere le nom de la 1ere variable a recuperer
  PyObject* tpl0= PyList_GetItem(pyVariables, 0);
  char* varname = NULL;
  if (PyString_Check(tpl0)) varname = PyString_AsString(tpl0);
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(tpl0)) varname = (char*)PyUnicode_AsUTF8(tpl0);
#endif

  //on recupere sol et solcenter ainsi que connectivite et taille zones Donneuses (tc)
  for (E_Int nd = 0; nd < nidomD; nd++)
  {
    PyObject* zoneD = PyList_GetItem(zonesD, nd);
    // #    include "getfromzoneDcompact_all.h"
    //*************************************
#    include "getfromzoneDcompact_allGradP.h"
    //*************************************
  }

  //on recupere sol et solcenter taille zones receuveuses, param_int et param_real (t)
  for (E_Int nd = 0; nd < nidomR; nd++)
  {
    PyObject* zoneR = PyList_GetItem(zonesR, nd);
    // #     include "getfromzoneRcompact_all.h"
    //*************************************
#    include "getfromzoneRcompact_allGradP.h"
    //*************************************
  }

  E_Int nbcomIBC = ipt_param_int[2];
  E_Int nbcomID  = ipt_param_int[3+nbcomIBC];

  E_Int shift_graph = nbcomIBC + nbcomID + 3;

  E_Int threadmax_sdm  = __NUMTHREADS__;
  E_Int ech            = ipt_param_int[ NoTransfert +shift_graph];
  E_Int nrac           = ipt_param_int[ ech +1 ];          //nb total de raccord
  E_Int nrac_inst      = ipt_param_int[ ech +2 ];          //nb total de raccord instationnaire
  E_Int timelevel      = ipt_param_int[ ech +3 ];          //nb de pas de temps stocker pour chaque raccord instationnaire
  E_Int nrac_steady    = nrac - nrac_inst;                 //nb total de raccord stationnaire


  //gestion nombre de pass pour raccord instationnaire
  E_Int pass_inst_deb=0;
  E_Int pass_inst_fin=1;
  if (nrac_inst > 0) pass_inst_fin=2;

  E_Int size_autorisation = nrac_steady+1;
  size_autorisation = K_FUNC::E_max(size_autorisation , nrac_inst+1);

  E_Int autorisation_transferts[pass_inst_fin][size_autorisation];

  // printf("nrac = %d, nrac_inst = %d, level= %d, it_target= %d , nitrun= %d \n",  nrac, nrac_inst, timelevel,it_target, NitRun);
  //on dimension tableau travail pour IBC
  E_Int nbRcvPts_mx =0; E_Int ibcTypeMax=0;
  for  (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
  {
    E_Int irac_deb= 0; E_Int irac_fin= nrac_steady;
    if (pass_inst == 1) { irac_deb = ipt_param_int[ ech + 4 + it_target ];
    irac_fin = ipt_param_int[ ech + 4 + it_target + timelevel ];}

    for  (E_Int irac=irac_deb; irac< irac_fin; irac++)
    {
      E_Int shift_rac =  ech + 4 + timelevel*2 + irac;
      if (ipt_param_int[ shift_rac+ nrac*10 + 1] > nbRcvPts_mx) nbRcvPts_mx = ipt_param_int[ shift_rac+ nrac*10 + 1];

      if (ipt_param_int[shift_rac+nrac*3] > ibcTypeMax)  ibcTypeMax =  ipt_param_int[shift_rac+nrac*3];


      E_Int irac_auto= irac-irac_deb;
      autorisation_transferts[pass_inst][irac_auto]=0;

      if (exploc == 1)  //if(rk==3 && exploc == 2) // Si on est en explicit local, on va autoriser les transferts entre certaines zones seulement en fonction de la ss-ite courante
      {
        E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;
        E_Int levelD = ipt_param_int[debut_rac + 25];
        E_Int levelR = ipt_param_int[debut_rac + 24];
        E_Int cyclD  = nitmax/levelD;

        // Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse
        if (levelD > levelR && num_passage == 1)
        {
          if ( nstep%cyclD==cyclD-1 || (nstep%cyclD==cyclD/2 && (nstep/cyclD)%2==1) )
          {
            autorisation_transferts[pass_inst][irac_auto]=1;
          }
          else {continue;}
        }
        // Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
        else if (levelD < levelR && num_passage == 1)
        {
          if (nstep%cyclD==1 || nstep%cyclD==cyclD/4 || nstep%cyclD== cyclD/2-1 || nstep%cyclD== cyclD/2+1 || nstep%cyclD== cyclD/2+cyclD/4 || nstep%cyclD== cyclD-1)
            { autorisation_transferts[pass_inst][irac_auto]=1; }
          else {continue;}
        }
        // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse
        else if (levelD == levelR && num_passage == 1)
        {
          if (nstep%cyclD==cyclD/2-1 || (nstep%cyclD==cyclD/2 && (nstep/cyclD)%2==0) || nstep%cyclD==cyclD-1)
            { autorisation_transferts[pass_inst][irac_auto]=1; }
          else {continue;}
        }
        // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse (cas du deuxieme passage)
        else if (levelD == ipt_param_int[debut_rac +24] && num_passage == 2)
        {
          if (nstep%cyclD==cyclD/2 && (nstep/cyclD)%2==1)
            { autorisation_transferts[pass_inst][irac_auto]=1; }
          else {continue;}
        }
        else {continue;}
        }
      // Sinon, on autorise les transferts entre ttes les zones a ttes les ss-ite
      else { autorisation_transferts[pass_inst][irac_auto]=1; }

    }
    }


  E_Int size = (nbRcvPts_mx/threadmax_sdm)+1; // on prend du gras pour gerer le residus
  E_Int r =  size % 8;
  if (r != 0) size  = size + 8 - r;           // on rajoute du bas pour alignememnt 64bits
  if (ibcTypeMax <=1 ) size = 0;              // tableau inutile : SP voir avec Ivan

  FldArrayF  tmp(size*17*threadmax_sdm);
  E_Float* ipt_tmp = tmp.begin();

  //# pragma omp parallel default(shared)  num_threads(1)
# pragma omp parallel default(shared)
  {

#ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num()+1;
    E_Int  Nbre_thread_actif = omp_get_num_threads(); // nombre de thread actif dans cette zone
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
#endif

    E_Int indR, type;
    E_Int indD0, indD, i, j, k, ncfLoc/*, nocf*/, indCoef, noi, sizecoefs, /*Nbchunk,*/ imd, jmd, imdjmd;

    vector<E_Float*> vectOfRcvFields(nvars);
    vector<E_Float*> vectOfDnrFields(nvars);

    //*************************************
    vector<E_Float*> vectOfGradRcvFields(nvars_grad);
    vector<E_Float*> vectOfGradDnrFields(nvars_grad);
    //*************************************

    //1ere pass_typ: IBC
    //2eme pass_typ: transfert

    for (E_Int ipass_typ=pass_deb; ipass_typ< pass_fin; ipass_typ++)
    {
    //1ere pass_inst: les raccords fixes
    //2eme pass_inst: les raccords instationnaires
    for (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
      {
        //printf("pass_inst = %d, level= %d \n",  pass_inst, nrac_inst_level );
        E_Int irac_deb= 0; E_Int irac_fin= nrac_steady;
        if(pass_inst == 1)
        {
          irac_deb = ipt_param_int[ ech + 4 + it_target             ];
          irac_fin = ipt_param_int[ ech + 4 + it_target + timelevel ];
        }

        for  (E_Int irac=irac_deb; irac< irac_fin; irac++)
        {

        E_Int irac_auto= irac-irac_deb;

        if (autorisation_transferts[pass_inst][irac_auto]==1)
        {

            E_Int shift_rac =  ech + 4 + timelevel*2 + irac;

            //printf("ipass_typ = %d, pass_inst= %d, irac=  %d, ithread= %d \n", ipass_typ,pass_inst,irac , ithread );
            // ibc: -1 : raccord ID, 0 : wallslip, 1 slip, 2 : log , 3 : Musker, 4 : outpress
            E_Int ibcType =  ipt_param_int[shift_rac+nrac*3];
            E_Int ibc = 1;
            if (ibcType < 0) ibc = 0;
            if(1-ibc != ipass_typ)  continue;

            E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ];
            E_Int loc      =  ipt_param_int[ shift_rac + nrac*9  +1 ]; //+1 a cause du nrac mpi
            E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ];
            E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 +1 ]; //neq fonction raccord rans/LES
            E_Int rotation =  ipt_param_int[ shift_rac + nrac*14 +1 ]; //flag pour periodicite azimutale

            //printf("irac=  %d, nvar_loc= %d,  ithread= %d \n",irac , nvars_loc, ithread );

            E_Int meshtype = ipt_ndimdxD[NoD + nidomD*6];
            E_Int cnNfldD  = ipt_ndimdxD[NoD + nidomD*7];
            E_Int* ptrcnd  = ipt_cnd[    NoD           ];

            if (loc == 0)
            {
              printf("Error: transferts optimises non code en vextex " SF_D3_ "\n", shift_rac + nrac*9  +1, NoD, NoR );
              //imd= ipt_ndimdxD[ NoD+ nidomD*4]; jmd= ipt_ndimdxD[ NoD + nidomD*5];
              imd = 0; jmd = 0;
            }
            else
            {
              for (E_Int eq = 0; eq < nvars_loc; eq++)
              {
                vectOfRcvFields[eq] = ipt_roR[ NoR] + eq*ipt_param_intR[ NoR ][ NDIMDX ];
                vectOfDnrFields[eq] = ipt_roD[ NoD] + eq*ipt_param_intR[ NoD ][ NDIMDX ];
              }
              imd= ipt_param_intR[ NoD ][ NIJK ]; jmd= ipt_param_intR[ NoD ][ NIJK+1];

              //*************************************
              for (E_Int eq = 0; eq < nvars_grad; eq++)
              {
                vectOfGradRcvFields[eq] = ipt_gradR[ NoR] + eq*ipt_param_intR[ NoR ][ NDIMDX ];
                vectOfGradDnrFields[eq] = ipt_gradD[ NoD] + eq*ipt_param_intR[ NoD ][ NDIMDX ];
              }
            //*************************************
              }

            imdjmd = imd*jmd;

            // std::cout << "loc" << loc << std::endl;
            // std::cout << "Dim trans Py " << " dim  R " << ipt_ndimdxR[NoR] << std::endl;
            // std::cout << "Dim trans Py " << " dim  D " << ipt_ndimdxD[NoD] << std::endl;
            // std::cout << "Dim trans Py " << " dim  imd " << ipt_ndimdxD[NoD + nidomD] << std::endl;
            // std::cout << "Dim trans Py " << " dim  jmd " << ipt_ndimdxD[NoD + nidomD*2] << std::endl;


            ////
            //  Interpolation parallele
            ////
            ////

            E_Int nbRcvPts = ipt_param_int[ shift_rac +  nrac*10 + 1 ];

            E_Int pos;
            pos  = ipt_param_int[ shift_rac + nrac*7 ]     ; E_Int* ntype      = ipt_param_int  + pos;
            pos  = pos +1 + ntype[0]                       ; E_Int* types      = ipt_param_int  + pos;
            pos  = ipt_param_int[ shift_rac + nrac*6      ]; E_Int* donorPts   = ipt_param_int  + pos;
            pos  = ipt_param_int[ shift_rac + nrac*12 + 1 ]; E_Int* rcvPts     = ipt_param_int  + pos;   // donor et receveur inverser car storage donor
            pos  = ipt_param_int[ shift_rac + nrac*8      ]; E_Float* ptrCoefs = ipt_param_real + pos;
            //printf("%d %d\n", ibcType, pos);

            E_Int nbInterpD = ipt_param_int[ shift_rac +  nrac ]; E_Float* xPC=NULL; E_Float* xPI=NULL; E_Float* xPW=NULL; E_Float* densPtr=NULL;
            if (ibc == 1)
            {
              xPC     = ptrCoefs + nbInterpD;
              xPI     = ptrCoefs + nbInterpD +3*nbRcvPts;
              xPW     = ptrCoefs + nbInterpD +6*nbRcvPts;
              densPtr = ptrCoefs + nbInterpD +9*nbRcvPts;
            }

            E_Int ideb      = 0;
            E_Int ifin      = 0;
            E_Int shiftCoef = 0;
            E_Int shiftDonor  = 0;

            for (E_Int ndtyp = 0; ndtyp < ntype[0]; ndtyp++)
            {
            type      = types[ifin];

            SIZECF(type, meshtype, sizecoefs);
            ifin =  ifin + ntype[ 1 + ndtyp];

            /*
            // *      New school: meilleur equilibrage, mais gestion looop dynamique rame...
            // *
            E_Int size_bc =  ifin-ideb;
            E_Int size_min=   16;
            //E_Int chunk = size_bc/Nbre_thread_actif;
            //if     (chunk < size_min && size_bc >= size_min) { chunk = size_min;}
            //else if(chunk < size_min && size_bc <  size_min) { chunk = size_bc ;}
            E_Int chunk = size_min;
            if(size_bc <  size_min) { chunk = size_bc ;}

            if      ( type == 0 ||  chunk <= 0) { Nbchunk = 1;                }
            else if ( chunk > 0)                { Nbchunk = size_bc/chunk;}

            chunk = size_bc/Nbchunk;

            E_Int r = size_bc - chunk*Nbchunk;

            #pragma omp for nowait schedule(dynamic,1)
            for (E_Int nd = 0; nd < Nbchunk; nd++)
            {
            */
            E_Int pt_deb, pt_fin;

            /// oldschool
            // Calcul du nombre de champs a traiter par chaque thread
            E_Int size_bc =  ifin-ideb;
            E_Int chunk   =  size_bc/Nbre_thread_actif;
            E_Int r       =  size_bc - chunk*Nbre_thread_actif;
            // pts traitees par thread
            if (ithread <= r)
              { pt_deb = ideb + (ithread-1)*(chunk+1);           pt_fin = pt_deb + (chunk+1); }
            else { pt_deb = ideb + (chunk+1)*r+(ithread-r-1)*chunk; pt_fin = pt_deb + chunk; }

            //Si type 0, calcul sequentiel
            if      ( type == 0 )
              { if (ithread ==1 ){ pt_deb = ideb; pt_fin = ifin;}
                else             { pt_deb = ideb; pt_fin = ideb;}
              }

            /// newschool suite
            //        if (nd  <  r) { pt_deb = ideb + nd*(chunk+1)               ; pt_fin = pt_deb + (chunk+1); }
            //        else          { pt_deb = ideb +    (chunk+1)*r+(nd-r)*chunk; pt_fin = pt_deb +  chunk;    }

            //printf(" irac= %d, NoR= %d, nvar=  %d, NoD= %d, Rans=  %d, rot= %d, fin= %d, type= %d, ithread= %d \n", irac, NoR, nvars_loc, NoD, pass_inst ,rotation, pt_fin , type,  ithread );
            //if(ithread <=8 && NoD==83 )  printf(" shift %d  %d %d %d  %d %d %d  %d \n", irac, NoR,NoD, ntype[ 1 + ndtyp],pt_deb,pt_fin  , type, ithread );
            //if(ithread <=8 && NoR==114 )  printf(" new   %d  %d %d %d  %d %d %d  %d \n", irac, NoR,NoD, ntype[ 1 + ndtyp],pt_deb,pt_fin  , type, ithread );


            noi       = shiftDonor;                             // compteur sur le tableau d indices donneur
            indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
                        E_Int shiftv =0;
            if (vartype != 24){
              if     (nvars_loc==5)
                {
#           include "commonInterpTransfers_reorder_5eq.h"
                }
              else if(nvars_loc==6)
                {
#           include "commonInterpTransfers_reorder_6eq.h"
                }
              else if(nvars_loc==19)
                {
#           include "commonInterpTransfers_reorder_19eq.h"
                }
              else
                {
#           include "commonInterpTransfers_reorder_neq.h"
                }
              //*************************************
              noi       = shiftDonor;                             // compteur sur le tableau d indices donneur
              indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
#           include "commonInterpTransfers_reorder_6eq_gradP.h"
              //*************************************
            }


            // Prise en compte de la periodicite par rotation
            if (rotation == 1)
              {
                E_Float* angle = ptrCoefs + nbInterpD;
#          include "includeTransfers_rotation.h"
              }
            // ibc
            if (ibc == 1)
              {
                //E_Int nvars = vectOfDnrFields.size();

                alpha    = ipt_param_realR[ NoR ][ ALPHAGRADP ];

                if (vartype == 23){
                  alpha = 1.;
                }

                Pr    = ipt_param_realR[ NoR ][ PRANDT ];
                Ts    = ipt_param_realR[ NoR ][ TEMP0 ];
                Cs    = ipt_param_realR[ NoR ][ CS ];
                muS   = ipt_param_realR[ NoR ][ XMUL0 ];
                cv    = ipt_param_realR[ NoR ][ CVINF ];
                gamma = ipt_param_realR[ NoR ][ GAMMA ];

                if (vartype == 22 || vartype == 23 || vartype == 24)
                {
                  E_Float cvgam = cv*(gamma-1.);

                for (E_Int noind = 0; noind < pt_fin-pt_deb; noind++)
                {
                  E_Int indR = rcvPts[noind+pt_deb];

                  (densPtr+nbRcvPts*7)[noind+pt_deb] = ((vectOfRcvFields[4][indR]*vectOfGradRcvFields[0][indR]+vectOfRcvFields[0][indR]*vectOfGradRcvFields[3][indR])*cvgam)/alpha + (densPtr+nbRcvPts*7)[noind+pt_deb]*(alpha-1.)/alpha;
                  (densPtr+nbRcvPts*8)[noind+pt_deb] = ((vectOfRcvFields[4][indR]*vectOfGradRcvFields[1][indR]+vectOfRcvFields[0][indR]*vectOfGradRcvFields[4][indR])*cvgam)/alpha + (densPtr+nbRcvPts*8)[noind+pt_deb]*(alpha-1.)/alpha;
                  (densPtr+nbRcvPts*9)[noind+pt_deb] = ((vectOfRcvFields[4][indR]*vectOfGradRcvFields[2][indR]+vectOfRcvFields[0][indR]*vectOfGradRcvFields[5][indR])*cvgam)/alpha + (densPtr+nbRcvPts*9)[noind+pt_deb]*(alpha-1.)/alpha;

                }
                }
              }//ibc
            //*
            //        } //chunk
            //*/

            ideb       =  ideb + ntype[ 1 + ndtyp];
            shiftCoef  = shiftCoef  +  ntype[1+ndtyp]*sizecoefs; //shift coef   entre 2 types successif
            shiftDonor = shiftDonor +  ntype[1+ndtyp];           //shift donor entre 2 types successif
              }// type
          }// autorisation transfert
          }//irac
      }//pass_inst
#pragma omp barrier
      }//ipass
  }// omp


  delete [] ipt_param_intR; delete [] ipt_roR; delete [] ipt_ndimdxD; delete [] ipt_roD; delete [] ipt_cnd;

  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
  RELEASESHAREDN(pyParam_int    , param_int    );
  RELEASESHAREDN(pyParam_real   , param_real   );

  Py_INCREF(Py_None);

  return Py_None;
}

