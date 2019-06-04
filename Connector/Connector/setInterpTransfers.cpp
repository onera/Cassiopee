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
#ifdef _MPI
#if defined(_WIN64)
# define __int64 long long
#endif
#include <mpi.h>
#endif

# include "connector.h"

using namespace std;
using namespace K_FLD;

#undef TimeShowsetinterp

#ifdef TimeShowsetinterp
E_Float time_in_0;
E_Float time_out_0;
E_Int rank_intp;
#endif

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
  if (!PYPARSETUPLEF(args, "OOOOOOO(ddd)", "OOOOOOO(fff)",
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
    int nvariables = PyList_Size(pyVariables);
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
          char* varname = PyBytes_AsString(PyUnicode_AsUTF8String(tpl0));
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
  if ( K_FUNC::E_abs(AngleX) > 0.) {dirR = 1; theta=AngleX;}
  else if (K_FUNC::E_abs(AngleY) > 0.){dirR=2; theta=AngleY;}
  else if (K_FUNC::E_abs(AngleZ) > 0.) {dirR=3; theta=AngleZ;} 
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

  if (!PYPARSETUPLE(args,
                    "OOOOOOOlllssssddd", "OOOOOOOiiissssddd",
                    "OOOOOOOlllssssfff", "OOOOOOOiiissssfff",
                    &zoneR, &zoneD, &pyVariables, &pyIndRcv, 
                    &pyIndDonor, &pyArrayTypes, &pyArrayCoefs, &loc , &vartype, &compact, &cellNVariable,
                    &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters, 
                    &AngleX, &AngleY, &AngleZ))
  {
      return NULL;
  }
  // Extraction de l'angle de rotation
  E_Int dirR = 0; E_Float theta = 0.;
  if ( K_FUNC::E_abs(AngleX) > 0.) {dirR = 1; theta=AngleX;}
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
  E_Int* rcvPts  = rcvPtsI->begin();
  nbRcvPts = rcvPtsI->getSize();

  vector<E_Float*> fieldsR;vector<E_Float*> fieldsD;
  vector<E_Int> posvarsD; vector<E_Int> posvarsR;
  E_Int* ptrcnd;
  char* eltTypeR; char* eltTypeD;

  //codage general (lent ;-) )
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
    char* varStringR;  vector<E_Int> locsR;
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
      int nvariables = PyList_Size(pyVariables);
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
            char* varname = PyBytes_AsString(PyUnicode_AsUTF8String(tpl0));
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
      E_Int noi = 0; // compteur sur le tableau d indices donneur
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

  if (!PYPARSETUPLE(args,
                    "OOOOOllllddddd", "OOOOOiiiiddddd",
                    "OOOOOllllfffff", "OOOOOiiiifffff",
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
    { E_Int pos_rac   = ipt_param_int[ shift_rac +irac];
      E_Int no_zR     = ipt_param_int[ pos_rac   +  6 ];
      //printf("pos_rac %d %d %d %d \n", pos_rac, no_zR, shift_rac+irac, irac); 
      PyObject* zoneR = PyList_GetItem(zonesR, no_zR); // domaine i
#     include "getfromzoneRcompact.h"
      ipt_ndimdxR[irac] = ndimdxR;      
      ipt_roR[irac]     = iptroR;
    }
 

  vector<E_Float*> vectOfRcvFields(nvars);
  vector<E_Float*> vectOfDnrFields(nvars);

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
    E_Int max_thread = min(nvars , __NUMTHREADS__);

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

    if(flagIbc>0)
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

      FldArrayF  tmp(size*14*threadmax_sdm);
      E_Float* ipt_tmp=  tmp.begin();

      E_Float* xPC     = ptrCoefs + nbInterpD;
      E_Float* xPI     = ptrCoefs + nbInterpD +3*nbRcvPts;
      E_Float* xPW     = ptrCoefs + nbInterpD +6*nbRcvPts;
      E_Float* densPtr = ptrCoefs + nbInterpD +9*nbRcvPts;
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
        if (varType == 1 || varType == 11)
          setIBCTransfersCommonVar1(ibcType, rcvPts, nbRcvPts, ideb, ifin, ithread, 
                                    xPC, xPC+nbRcvPts, xPC+nbRcvPts*2,
                                    xPW, xPW+nbRcvPts, xPW+nbRcvPts*2,
                                    xPI, xPI+nbRcvPts, xPI+nbRcvPts*2, 
                                    densPtr, densPtr+nbRcvPts, //dens + press
                                    densPtr+nbRcvPts*2, densPtr+nbRcvPts*3, densPtr+nbRcvPts*4, // vx + vy + vz 
                                    densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, // utau + yplus
                                    densPtr+nbRcvPts*7, densPtr+nbRcvPts*8, densPtr+nbRcvPts*9, densPtr+nbRcvPts*10, densPtr+nbRcvPts*11,
                                    ipt_tmp, size,
                                    gamma, cv, muS, Cs, Ts, Pr,
                                    vectOfDnrFields, vectOfRcvFields);
        else if (varType == 2 || varType == 21)
        { 
          setIBCTransfersCommonVar2(ibcType, rcvPts, nbRcvPts, ideb, ifin, ithread,
                                    xPC, xPC+nbRcvPts, xPC+nbRcvPts*2,
                                    xPW, xPW+nbRcvPts, xPW+nbRcvPts*2,
                                    xPI, xPI+nbRcvPts, xPI+nbRcvPts*2, 
                                    densPtr, densPtr+nbRcvPts, //dens + press
                                    densPtr+nbRcvPts*2, densPtr+nbRcvPts*3, densPtr+nbRcvPts*4, // vx + vy + vz 
                                    densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, // utau + yplus
                                    densPtr+nbRcvPts*7, densPtr+nbRcvPts*8, densPtr+nbRcvPts*9, densPtr+nbRcvPts*10, densPtr+nbRcvPts*11,
                                    ipt_tmp, size,
                                    gamma, cv, muS, Cs, Ts, Pr,
                                    vectOfDnrFields, vectOfRcvFields);
        }
        else if (varType == 3 || varType == 31)
          setIBCTransfersCommonVar3(ibcType, rcvPts, nbRcvPts, ideb, ifin, ithread,
                                    xPC, xPC+nbRcvPts, xPC+nbRcvPts*2,
                                    xPW, xPW+nbRcvPts, xPW+nbRcvPts*2,
                                    xPI, xPI+nbRcvPts, xPI+nbRcvPts*2, 
                                    densPtr, densPtr+nbRcvPts, //dens + press
                                    densPtr+nbRcvPts*2, densPtr+nbRcvPts*3, densPtr+nbRcvPts*4, // vx + vy + vz 
                                    densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, // utau + yplus
                                    densPtr+nbRcvPts*7, densPtr+nbRcvPts*8, densPtr+nbRcvPts*9, densPtr+nbRcvPts*10, densPtr+nbRcvPts*11,
                                    ipt_tmp, size,
                                    gamma, cv, muS, Cs, Ts, Pr,
                                    vectOfDnrFields, vectOfRcvFields);

      } // Fin zone // omp
    } //ibc

  }//irac


 delete [] ipt_ndimdxR; delete [] ipt_roR;

 RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL); 
 RELEASESHAREDN(pyParam_int    , param_int    );
 RELEASESHAREDN(pyParam_real   , param_real   );

 Py_INCREF(Py_None);
 return Py_None;
}
//=============================================================================
// Idem: in place + from zone + tc compact au niveau base 
//=============================================================================
PyObject* K_CONNECTOR::___setInterpTransfers(PyObject* self, PyObject* args)
{
  PyObject *zonesR, *zonesD;

  PyObject *pyVariables;
  PyObject *pyParam_int, *pyParam_real;
  E_Int vartype, bctype, type_transfert, no_transfert, It_target, nstep, nitmax;
  E_Int rk, exploc, num_passage;
  E_Float gamma, cv, muS, Cs, Ts;

#ifdef TimeShowsetinterp  
PyObject *timecount;
      if (!PYPARSETUPLE(args,
                        "OOOOOlllllllllldddddO", "OOOOOiiiiiiiiiidddddO",
                        "OOOOOllllllllllfffffO", "OOOOOiiiiiiiiiifffffO",
                        &zonesR, &zonesD, &pyVariables, &pyParam_int,  &pyParam_real, &It_target, &vartype,
                        &bctype, &type_transfert, &no_transfert, &nstep, &nitmax, &rk, &exploc, &num_passage, &gamma, &cv, &muS, &Cs, &Ts, &timecount))
      {
          return NULL;
      }

#else
      if (!PYPARSETUPLE(args,
                        "OOOOOllllllllllddddd", "OOOOOiiiiiiiiiiddddd",
                        "OOOOOllllllllllfffff", "OOOOOiiiiiiiiiifffff",
                        &zonesR, &zonesD, &pyVariables, &pyParam_int,  &pyParam_real, &It_target, &vartype,
                        &bctype, &type_transfert, &no_transfert, &nstep, &nitmax, &rk, &exploc, &num_passage, &gamma, &cv, &muS, &Cs, &Ts))
      {
          return NULL;
      }
#endif

  
  //E_Int init=0;

#ifdef TimeShowsetinterp
  FldArrayF* timecnt;
  E_Float* ipt_timecount = NULL;
  K_NUMPY::getFromNumpyArray(timecount, timecnt, true);
  ipt_timecount = timecnt->begin();

#ifdef _MPI
  MPI_Initialized( &init );
  if(init)
  {
    MPI_Comm_rank (MPI_COMM_WORLD, &rank_intp);  
  }
#endif
#endif


#ifdef TimeShowsetinterp
 if ( rank_intp == 0)
 {
   time_in_0 = omp_get_wtime();
 }
#endif

  //E_Int bcType   = E_Int(bctype);    // 0 : wallslip; 1: noslip; 2: log law of wall; 3: Musker law of wall
  E_Int it_target= E_Int(It_target);

  /* varType : 
     1  : conservatives, 
     11 : conservatives + ronutildeSA 
     2  : (ro,u,v,w,t)
     21 : (ro,u,v,w,t) + ronutildeSA 
     3  : (ro,u,v,w,p)     
     31 : (ro,u,v,w,p) + ronutildeSA */
  E_Int varType = E_Int(vartype); 

  //gestion nombre de pass pour ID et/ou IBC
  E_Int TypeTransfert = E_Int(type_transfert);
  E_Int pass_deb, pass_fin;
  if     (TypeTransfert==0) { pass_deb =1; pass_fin =2; }//ID
  else if(TypeTransfert==1) { pass_deb =0; pass_fin =1; }//IBCD
  else                      { pass_deb =0; pass_fin =2; }//ALL

  E_Int NoTransfert  = E_Int(no_transfert);

  
  E_Int kmd, cnNfldD, nvars, meshtype;

  if( vartype <= 3 &&  vartype >= 1) nvars =5;
  else                               nvars =6;

  E_Int nidomR   = PyList_Size(zonesR);
  E_Int nidomD   = PyList_Size(zonesD);

  //pointeur pour stocker solution au centre ET au noeud 
  E_Int* ipt_ndimdxR; E_Int* ipt_ndimdxD; E_Int** ipt_cnd;
  E_Float** ipt_roR; E_Float** ipt_roD; E_Float** ipt_roR_vert;  E_Float** ipt_roD_vert;

  ipt_ndimdxR      = new E_Int[nidomR*2];   // on stocke ndimdx  en centre et vertexe
  ipt_roR          = new E_Float*[nidomR*2];
  ipt_roR_vert     = ipt_roR + nidomR;

  ipt_ndimdxD      = new E_Int[nidomD*8];  //on stocke ndimdx, imd, jmd, en centre et vertexe, meshtype et cnDfld
  ipt_cnd          = new E_Int*[nidomD];
  ipt_roD          = new E_Float*[nidomD*2];
  ipt_roD_vert     = ipt_roD + nidomD;


  vector<PyArrayObject*> hook ;

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

  /*----------------------------------*/
  /* Get the Shift values for Padding */
  /*----------------------------------*/

  E_Int** ipt_param_int_Shift;

  ipt_param_int_Shift = new E_Int*[nidomR];

  for (E_Int nd = 0; nd < nidomR; nd++)
  {
    
  PyObject* zone = PyList_GetItem(zonesR, nd);  
  PyObject* own  = K_PYTREE::getNodeFromName1(zone , ".Solver#ownData");
  if (own != NULL)
  {
    PyObject* paramint = K_PYTREE::getNodeFromName1(own, "Parameter_int");
    if (paramint != NULL)
    {
      ipt_param_int_Shift[nd] = K_PYTREE::getValueAI(paramint, hook);      
    }
    else
    {
      ipt_param_int_Shift[nd] = NULL;    
    }       
  }  
  else
   {
      ipt_param_int_Shift[nd] = NULL;    
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

  //On recupere le nom de la 1ere variable a recuperer 
  PyObject* tpl0= PyList_GetItem(pyVariables, 0); 
  char* varname = NULL;
  if (PyString_Check(tpl0)) varname = PyString_AsString(tpl0);
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(tpl0)) varname = PyBytes_AsString(PyUnicode_AsUTF8String(tpl0)); 
#endif

  //on recupere sol et solcenter ainsi que connectivite et taille zones Donneuses (tc)
  for (E_Int nd = 0; nd < nidomD; nd++)
  {  
     PyObject* zoneD = PyList_GetItem(zonesD, nd);
#    include "getfromzoneDcompact_all.h"
  }

  //on recupere sol et solcenter et taille zones receuveuses (t)
  for (E_Int nd = 0; nd < nidomR; nd++)
    {  
      PyObject* zoneR = PyList_GetItem(zonesR, nd);
#     include "getfromzoneRcompact_all.h"
    }

  E_Int nbcomIBC = ipt_param_int[1];
  E_Int nbcomID  = ipt_param_int[2+nbcomIBC];
  
  E_Int shift_graph = nbcomIBC + nbcomID + 2;

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

  //E_Int autorisation_transferts[nrac_inst+1][nrac_steady+1];
  E_Int autorisation_transferts[pass_inst_fin][size_autorisation];
  //E_Int** autorisation_transfert = new E_int*[pass_inst_fin]
  //for (E_Int i = 0 ; i < pass_inst_fin; i++)
  //    {
  //	E_Int* autorisation_transfert[i] = new E_Int[size_autorisation];
  //   }


  // printf("nrac = %d, nrac_inst = %d, level= %d, it_target= %d , nitrun= %d \n",  nrac, nrac_inst, timelevel,it_target, NitRun);
  //on dimension tableau travail pour IBC
  E_Int nbRcvPts_mx =0;
  for  (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
  {
      E_Int irac_deb= 0; E_Int irac_fin= nrac_steady;
      if (pass_inst == 1) { irac_deb = ipt_param_int[ ech + 4 + it_target ];
	   irac_fin = ipt_param_int[ ech + 4 + it_target + timelevel ];}


      for  (E_Int irac=irac_deb; irac< irac_fin; irac++)
	{ 
	  E_Int shift_rac =  ech + 4 + timelevel*2 + irac;
	  if( ipt_param_int[ shift_rac+ nrac*10 + 1] > nbRcvPts_mx) nbRcvPts_mx = ipt_param_int[ shift_rac+ nrac*10 + 1];
	  //E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ];
	  //E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ];

    E_Int irac_auto= irac-irac_deb;
	  autorisation_transferts[pass_inst][irac_auto]=0;

	  if(rk==3 && exploc == 2) // Si on est en explicit local, on va autoriser les transferts entre certaines zones seulement en fonction de la ss-ite courante

	    {

	      E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;     
	      E_Int levelD = ipt_param_int[debut_rac + 25];
	      E_Int levelR = ipt_param_int[debut_rac + 24];
	      E_Int cyclD  = nitmax/levelD;

		// Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse   
		if (levelD > levelR && num_passage == 1)		
		  {
		    if (nstep%cyclD==cyclD-1 or nstep%cyclD==cyclD/2 && (nstep/cyclD)%2==1)
		      {
			autorisation_transferts[pass_inst][irac_auto]=1;
			//if( ipt_param_int[ shift_rac+ nrac*10 + 1] > nbRcvPts_mx) nbRcvPts_mx = ipt_param_int[ shift_rac+ nrac*10 + 1];
		      }
		    else {continue;}
		  }
		// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
		else if (levelD < levelR && num_passage == 1) 
		  {
		    if (nstep%cyclD==1 || nstep%cyclD==cyclD/4 || nstep%cyclD== cyclD/2-1 || nstep%cyclD== cyclD/2+1 || nstep%cyclD== cyclD/2+cyclD/4 || nstep%cyclD== cyclD-1)
                         {
			   autorisation_transferts[pass_inst][irac_auto]=1;
			   //if( ipt_param_int[ shift_rac+ nrac*10 + 1] > nbRcvPts_mx) nbRcvPts_mx = ipt_param_int[ shift_rac+ nrac*10 + 1];
			 }
		    else {continue;}
		  }
		// Le pas de temps de la zone donneuse est egal a celui de la zone receveuse
		else if (levelD == levelR && num_passage == 1)
		  {
		    if (nstep%cyclD==cyclD/2-1 || (nstep%cyclD==cyclD/2 && (nstep/cyclD)%2==0) || nstep%cyclD==cyclD-1) 
		      { 
			autorisation_transferts[pass_inst][irac_auto]=1; 
			//if( ipt_param_int[ shift_rac+ nrac*10 + 1] > nbRcvPts_mx) nbRcvPts_mx = ipt_param_int[ shift_rac+ nrac*10 + 1];
		      }
		    else {continue;}
		  }
		// Le pas de temps de la zone donneuse est egal a celui de la zone receveuse (cas du deuxieme passage)   
		else if (levelD == ipt_param_int[debut_rac +24] && num_passage == 2)
		  {
		  if (nstep%cyclD==cyclD/2 && (nstep/cyclD)%2==1)
		    {
		      autorisation_transferts[pass_inst][irac_auto]=1;
		      //if( ipt_param_int[ shift_rac+ nrac*10 + 1] > nbRcvPts_mx) nbRcvPts_mx = ipt_param_int[ shift_rac+ nrac*10 + 1];
		    }
		  else                                            {continue;}
		  }
		else {continue;} 
	    }
           // Sinon, on autorise les transferts entre ttes les zones a ttes les ss-ite
	   else 
	     { 
	       autorisation_transferts[pass_inst][irac_auto]=1; 
	       //if( ipt_param_int[ shift_rac+ nrac*10 + 1] > nbRcvPts_mx) nbRcvPts_mx = ipt_param_int[ shift_rac+ nrac*10 + 1];
	     }
	    
	}
    }

  E_Int size = (nbRcvPts_mx/threadmax_sdm)+1; // on prend du gras pour gerer le residus
  E_Int r =  size % 8;
  if (r != 0) size  = size + 8 - r;           // on rajoute du bas pour alignememnt 64bits
  if (bctype <=1 ) size = 0;                  // tableau inutile : SP voir avec Ivan 

  FldArrayF  tmp(size*14*threadmax_sdm);
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

  //1ere pass_typ: IBC
  //2eme pass_typ: transfert

  for  (E_Int ipass_typ=pass_deb; ipass_typ< pass_fin; ipass_typ++)
  {
   //1ere pass_inst: les raccord fixe
   //2eme pass_inst: les raccord instationnaire
   for  (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
   {
    //printf("pass_inst = %d, level= %d \n",  pass_inst, nrac_inst_level );
    E_Int irac_deb= 0; E_Int irac_fin= nrac_steady;
    if(pass_inst == 1)
    { 
      irac_deb = ipt_param_int[ ech + 4 + it_target             ];
      irac_fin = ipt_param_int[ ech + 4 + it_target + timelevel ];  
    }

    //cout << irac_deb << " "  << endl;
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
      if( 1-ibc != ipass_typ)  continue;

      E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ];
      E_Int loc      =  ipt_param_int[ shift_rac + nrac*9  +1 ]; //+1 a cause du nrac mpi
      E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 +1 ];
      E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 +1 ]; //neq fonction raccord rans/LES
      E_Int rotation =  ipt_param_int[ shift_rac + nrac*14 +1 ]; //flag pour periodicite azimutale

      E_Int meshtype = ipt_ndimdxD[NoD + nidomD*6];
      E_Int cnNfldD  = ipt_ndimdxD[NoD + nidomD*7];
      E_Int* ptrcnd  = ipt_cnd[    NoD           ];

      if (loc == 0)
      {
        for (E_Int eq = 0; eq < nvars_loc; eq++)
        {
         vectOfRcvFields[eq] = ipt_roR_vert[ NoR] + eq*ipt_ndimdxR[ NoR + nidomR  ];
         vectOfDnrFields[eq] = ipt_roD_vert[ NoD] + eq*ipt_ndimdxD[ NoD + nidomD*3];
        }
         imd= ipt_ndimdxD[ NoD+ nidomD*4]; jmd= ipt_ndimdxD[ NoD + nidomD*5];
      }
      else
      {
        for (E_Int eq = 0; eq < nvars_loc; eq++)
        {
         vectOfRcvFields[eq] = ipt_roR[ NoR] + eq*(ipt_ndimdxR[ NoR ] + ipt_param_int_Shift[NoR][66] );
         vectOfDnrFields[eq] = ipt_roD[ NoD] + eq*(ipt_ndimdxD[ NoD ] + ipt_param_int_Shift[NoD][66] );
        }
         imd= ipt_ndimdxD[ NoD+ nidomD  ]; jmd= ipt_ndimdxD[ NoD+ nidomD*2];
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
       //E_Int nbDonPts = ipt_param_int[ shift_rac                ];
  
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

      //printf(" irac= %d, NoR= %d, nvar=  %d, NoD= %d, Rans=  %d, rot= %d, fin= %d, ithread= %d \n", irac, NoR, nvars_loc, NoD, pass_inst ,rotation, pt_fin , ithread );
      //if(ithread <=8 && NoD==83 )  printf(" shift %d  %d %d %d  %d %d %d  %d \n", irac, NoR,NoD, ntype[ 1 + ndtyp],pt_deb,pt_fin  , type, ithread );
      //if(ithread <=8 && NoR==114 )  printf(" new   %d  %d %d %d  %d %d %d  %d \n", irac, NoR,NoD, ntype[ 1 + ndtyp],pt_deb,pt_fin  , type, ithread );

          noi       = shiftDonor;                             // compteur sur le tableau d indices donneur
          indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
          if     (nvars_loc==5)
          {
#           include "commonInterpTransfers_reorder_5eq.h" 
          }
          else if(nvars_loc==6)
          {
#           include "commonInterpTransfers_reorder_6eq.h" 
          }
          else
          {
#           include "commonInterpTransfers_reorder_neq.h"
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
            E_Int nvars = vectOfDnrFields.size();
            if ( (ibcType==2 || (ibcType==3)) && nvars < 6)
            {
              printf("Warning: setIBCTransfersCommonVar: number of variables (<6) inconsistent with ibctype (wall law).\n"); 
            }
            else 
            {
            if (varType == 1 || varType == 11)
              setIBCTransfersCommonVar1(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread,
                                        xPC, xPC+nbRcvPts, xPC+nbRcvPts*2,
                                        xPW, xPW+nbRcvPts, xPW+nbRcvPts*2,
                                        xPI, xPI+nbRcvPts, xPI+nbRcvPts*2, 
                                        densPtr, densPtr+nbRcvPts, //dens + press
                                        densPtr+nbRcvPts*2, densPtr+nbRcvPts*3, densPtr+nbRcvPts*4, // vx + vy + vz 
                                        densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, // utau + yplus
                                        densPtr+nbRcvPts*7, densPtr+nbRcvPts*8, densPtr+nbRcvPts*9, densPtr+nbRcvPts*10, densPtr+nbRcvPts*11,
                                        ipt_tmp, size,
                                        gamma, cv, muS, Cs, Ts, Pr,
                                        vectOfDnrFields, vectOfRcvFields);
            else if (varType == 2 || varType == 21)
            {
              setIBCTransfersCommonVar2(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread,
                                        xPC    , xPC     +nbRcvPts, xPC     +nbRcvPts*2,
                                        xPW    , xPW     +nbRcvPts, xPW     +nbRcvPts*2,
                                        xPI    , xPI     +nbRcvPts, xPI     +nbRcvPts*2, 
                                        densPtr, densPtr+nbRcvPts, //dens + press
                                        densPtr+nbRcvPts*2, densPtr+nbRcvPts*3, densPtr+nbRcvPts*4, // vx + vy + vz 
                                        densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, // utau + yplus
                                        densPtr+nbRcvPts*7, densPtr+nbRcvPts*8, densPtr+nbRcvPts*9, densPtr+nbRcvPts*10, densPtr+nbRcvPts*11,
                                        ipt_tmp, size,
                                        gamma, cv, muS, Cs, Ts, Pr,
                                        vectOfDnrFields, vectOfRcvFields);
            }
            else if (varType == 3 || varType == 31)
              setIBCTransfersCommonVar3(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread,
                                        xPC    , xPC     +nbRcvPts, xPC     +nbRcvPts*2,
                                        xPW    , xPW     +nbRcvPts, xPW     +nbRcvPts*2,
                                        xPI    , xPI     +nbRcvPts, xPI     +nbRcvPts*2, 
                                        densPtr, densPtr+nbRcvPts, //dens + press
                                        densPtr+nbRcvPts*2, densPtr+nbRcvPts*3, densPtr+nbRcvPts*4, // vx + vy + vz 
                                        densPtr+nbRcvPts*5, densPtr+nbRcvPts*6, // utau + yplus
                                        densPtr+nbRcvPts*7, densPtr+nbRcvPts*8, densPtr+nbRcvPts*9, densPtr+nbRcvPts*10, densPtr+nbRcvPts*11,
                                        ipt_tmp, size,
                                        gamma, cv, muS, Cs, Ts, Pr,
                                        vectOfDnrFields, vectOfRcvFields);
          }
          }//ibc          
  //*
  //        } //chunk
  //*/
          ideb       = ideb + ifin;
          shiftCoef  = shiftCoef   +  ntype[1+ndtyp]*sizecoefs; //shift coef   entre 2 types successif
          shiftDonor = shiftDonor +  ntype[1+ndtyp];           //shift donor entre 2 types successif
       }// type 
      }// autorisation transfert
    }//irac
   }//pass_inst
  #pragma omp barrier 
  }//ipass
  }// omp


  delete [] ipt_ndimdxR; delete [] ipt_roR; delete [] ipt_ndimdxD; delete [] ipt_roD; delete [] ipt_cnd;
  delete [] ipt_param_int_Shift;

  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL); 
  RELEASESHAREDN(pyParam_int    , param_int    );
  RELEASESHAREDN(pyParam_real   , param_real   );

  Py_INCREF(Py_None);

#ifdef TimeShowsetinterp  
if ( rank_intp == 0)
    {
      time_out_0 = omp_get_wtime();
      ipt_timecount[2]  =  ipt_timecount[2] + time_out_0 - time_in_0;
    }
#endif

  return Py_None;
}
