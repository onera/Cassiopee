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
# include "connector.h"
using namespace std;
using namespace K_FLD;

//=============================================================================
/* Transfert de champs conservatifs sous forme de numpy */
//=============================================================================
PyObject* K_CONNECTOR::setIBCTransfersD(PyObject* self, PyObject* args)
{
  PyObject *arrayD, *pyIndDonor, *pyArrayTypes, *pyArrayCoefs;
  PyObject *pyArrayXPC, *pyArrayXPI, *pyArrayXPW;
  PyObject *pyArrayYPC, *pyArrayYPI, *pyArrayYPW;
  PyObject *pyArrayZPC, *pyArrayZPI, *pyArrayZPW;
  PyObject *pyArrayDens, *pyArrayPressure;
  PyObject *pyArrayVx, *pyArrayVy, *pyArrayVz; 
  PyObject *pyArrayUtau, *pyArrayYplus;
  E_Int bctype, vartype;
  E_Float gamma, cv, muS, Cs, Ts;

  if (!PYPARSETUPLE(args,
                    "OOOOOOOOOOOOOOOOOOOOllddddd", "OOOOOOOOOOOOOOOOOOOOiiddddd",
                    "OOOOOOOOOOOOOOOOOOOOllfffff", "OOOOOOOOOOOOOOOOOOOOiifffff",
                    &arrayD, &pyIndDonor, &pyArrayTypes, &pyArrayCoefs, 
                    &pyArrayXPC, &pyArrayYPC, &pyArrayZPC,
                    &pyArrayXPW, &pyArrayYPW, &pyArrayZPW,
                    &pyArrayXPI, &pyArrayYPI, &pyArrayZPI,
                    &pyArrayDens, &pyArrayPressure, 
                    &pyArrayVx, &pyArrayVy, &pyArrayVz, 
                    &pyArrayUtau, &pyArrayYplus,
                    &bctype, &vartype, &gamma, &cv, &muS, &Cs, &Ts))
  {
      return NULL;
  }
  E_Int bcType = E_Int(bctype); // 0: wallslip; 1: noslip; 2: log law of wall; 3: Musker law of wall
  /* varType : 
     1  : conservatives, 
     11 : conservatives + ronutildeSA 
     2  : (ro,u,v,w,t)
     21 : (ro,u,v,w,t) + ronutildeSA 
     3  : (ro,u,v,w,p)     
     31 : (ro,u,v,w,p) + ronutildeSA */
  E_Int varType = E_Int(vartype);

  /*---------------------------------------------*/
  /* Extraction des infos sur le domaine donneur */
  /*---------------------------------------------*/
  E_Int imd, jmd, kmd, imdjmd;
  FldArrayF* fd; FldArrayI* cnd;
  char* varStringD; char* eltTypeD;
  E_Int resd = K_ARRAY::getFromArray(arrayD, varStringD, fd, 
                                     imd, jmd, kmd, cnd, eltTypeD, true); 
  E_Int* ptrcnd = NULL;
  E_Int cndSize = 0; E_Int cnNfldD = 0;
  if (resd != 2 && resd != 1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "setIBCTransfersD: 1st arg is not a valid array.");
    return NULL; 
  }
  E_Int meshtype = resd;// 1 : structure, 2 non structure
  if (resd == 2)
  {
    if (strcmp(eltTypeD, "TETRA") != 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "setIBCTransfersD: donor array must be a TETRA if unstructured.");
      RELEASESHAREDB(resd, arrayD, fd, cnd); 
      return NULL; 
    }
    ptrcnd  = cnd->begin();
    cndSize = cnd->getSize();
    cnNfldD = cnd->getNfld();
  }
  E_Int nvars;
  if      ( varType ==  1 || varType ==  2 || varType == 3 )  nvars = 5;
  else if ( varType == 11 || varType == 21 || varType == 31 ) nvars = 6;
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "setIBCTransfersD: varType value is not valid.");
    RELEASESHAREDB(resd, arrayD, fd, cnd); 
    return NULL; 
  }
  if (fd->getNfld() != nvars) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "setIBCTransfersD: number of donor variables is not consistent with varType.");
    RELEASESHAREDB(resd, arrayD, fd, cnd); 
    return NULL; 
  }


# include "extract_interpD.h"

  if (res_donor*res_type*res_coef == 0) 
  {
    RELEASESHAREDB(resd, arrayD, fd, cnd); 
    if (res_donor != 0) { RELEASESHAREDN(pyIndDonor  , donorPtsI  );}
    if (res_type  != 0) { RELEASESHAREDN(pyArrayTypes, typesI     );}
    if (res_coef  != 0) { RELEASESHAREDN(pyArrayCoefs, donorCoefsF);}
    PyErr_SetString(PyExc_TypeError,"setInterpTransfersD: 2nd and 3rd arg must be a numpy of integers. 4th arg a numpy floats ");
    return NULL;
  }

# include "IBC/extract_IBC.h"

  if (okc1*okc2*okc3 == 0)
  {
    RELEASESHAREDB(resd, arrayD, fd, cnd); 
    RELEASESHAREDN(pyIndDonor  , donorPtsI  );
    RELEASESHAREDN(pyArrayTypes, typesI     );
    RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
    if (okc1 != 0) { RELEASESHAREDN(pyArrayXPC  , coordxPC  );}
    if (okc2 != 0) { RELEASESHAREDN(pyArrayYPC  , coordyPC  );}
    if (okc3 != 0) { RELEASESHAREDN(pyArrayZPC  , coordzPC  );}
    PyErr_SetString(PyExc_TypeError, 
                    "setIBCTransfersD: coordinates of corrected points are invalid.");
    return NULL;
  }
  if (okw1*okw2*okw3 == 0 )
  {
    RELEASESHAREDB(resd, arrayD, fd, cnd); 
    RELEASESHAREDN(pyIndDonor  , donorPtsI  );
    RELEASESHAREDN(pyArrayTypes, typesI     );
    RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
    RELEASESHAREDN(pyArrayXPC  , coordxPC  );
    RELEASESHAREDN(pyArrayYPC  , coordyPC  );
    RELEASESHAREDN(pyArrayZPC  , coordzPC  );
    if (okw1 != 0) { RELEASESHAREDN(pyArrayXPW  , coordxPW  );}
    if (okw2 != 0) { RELEASESHAREDN(pyArrayYPW  , coordyPW  );}
    if (okw3 != 0) { RELEASESHAREDN(pyArrayZPW  , coordzPW  );}
    PyErr_SetString(PyExc_TypeError, 
                    "setIBCTransfersD: coordinates of wall points are invalid.");
    return NULL;
  }
  if (oki1*oki2*oki3 == 0 )
  {
    RELEASESHAREDB(resd, arrayD, fd, cnd); 
    RELEASESHAREDN(pyIndDonor  , donorPtsI  );
    RELEASESHAREDN(pyArrayTypes, typesI     );
    RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
    RELEASESHAREDN(pyArrayXPC  , coordxPC  );
    RELEASESHAREDN(pyArrayYPC  , coordyPC  );
    RELEASESHAREDN(pyArrayZPC  , coordzPC  );
    RELEASESHAREDN(pyArrayXPW  , coordxPW  );
    RELEASESHAREDN(pyArrayYPW  , coordyPW  );
    RELEASESHAREDN(pyArrayZPW  , coordzPW  );
    if (oki1 != 0) { RELEASESHAREDN(pyArrayXPI  , coordxPI  );}
    if (oki2 != 0) { RELEASESHAREDN(pyArrayYPI  , coordyPI  );}
    if (oki3 != 0) { RELEASESHAREDN(pyArrayZPI  , coordzPI  );}
    PyErr_SetString(PyExc_TypeError, 
                    "setIBCTransfersD: coordinates of interpolated points are invalid.");
    return NULL;
  }
  if (okD*okP*okVx*okVy*okVz*okU*okY == 0 )
  {
    RELEASESHAREDB(resd, arrayD, fd, cnd); 
    RELEASESHAREDN(pyIndDonor  , donorPtsI  );
    RELEASESHAREDN(pyArrayTypes, typesI     );
    RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
    RELEASESHAREDN(pyArrayXPC  , coordxPC  );
    RELEASESHAREDN(pyArrayYPC  , coordyPC  );
    RELEASESHAREDN(pyArrayZPC  , coordzPC  );
    RELEASESHAREDN(pyArrayXPW  , coordxPW  );
    RELEASESHAREDN(pyArrayYPW  , coordyPW  );
    RELEASESHAREDN(pyArrayZPW  , coordzPW  );
    RELEASESHAREDN(pyArrayXPI  , coordxPI  );
    RELEASESHAREDN(pyArrayYPI  , coordyPI  );
    RELEASESHAREDN(pyArrayZPI  , coordzPI  );
    if (okD != 0) { RELEASESHAREDN(pyArrayDens    , densF    );}
    if (okP != 0) { RELEASESHAREDN(pyArrayPressure, pressF  );}
    if (okVx != 0) { RELEASESHAREDN(pyArrayVx      , vxF     );}
    if (okVy != 0) { RELEASESHAREDN(pyArrayVy      , vyF     );}
    if (okVz != 0) { RELEASESHAREDN(pyArrayVz      , vzF     );}
    if (okU != 0) { RELEASESHAREDN(pyArrayUtau    , utauF   );}
    if (okY != 0) { RELEASESHAREDN(pyArrayYplus   , yplusF  );}
    PyErr_SetString(PyExc_TypeError, 
                    "setIBCTransfersD: Post array are invalid.");
    return NULL;
  }

  PyObject* tpl  = K_ARRAY::buildArray(nvars, varStringD, nbRcvPts, 1, 1);
  E_Float* frp   = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fieldROut(nbRcvPts, nvars, frp, true);
  //
  //utile la mise a zero???
  //
  fieldROut.setAllValuesAtNull();

  // Transferts 
  // Types valides: 2, 3, 4, 5


  vector<E_Float*> vectOfRcvFields(nvars);
  vector<E_Float*> vectOfDnrFields(nvars);
  for (E_Int eq = 0; eq < nvars; eq++)
   {
    vectOfRcvFields[eq] = fieldROut.begin(eq+1);
    vectOfDnrFields[eq] = fd->begin(eq+1); 
   }


  // tableau temporaire pour utiliser la routine commune setIBCTransfersCommon
  FldArrayI rcvPtsI(nbRcvPts); E_Int*  rcvPts = rcvPtsI.begin();

////
////
//  Interpolation parallele
////  
////  
# include "commonInterpTransfers_direct.h" 

    E_Int threadmax_sdm  = __NUMTHREADS__;

    E_Int size = (nbRcvPts/threadmax_sdm)+1; // on prend du gras pour gerer le residus
          size = size + size % 8;                  // on rajoute du bas pour alignememnt 64bits

    FldArrayF  tmp(size*13*threadmax_sdm);
    E_Float* ipt_tmp=  tmp.begin();

#pragma omp parallel default(shared)
 {


  #pragma omp for
  for (E_Int noind = 0; noind < nbRcvPts; noind++) rcvPts[noind] = noind;

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

   if (varType == 1 || varType == 11) 
     setIBCTransfersCommonVar1(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread, 
			      xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, 
            density, pressure, 
            vx, vy, vz, 
            utau, yplus,
            NULL, NULL, NULL, NULL, NULL,
            ipt_tmp, size,
            gamma, cv, muS, Cs, Ts, 0.71,
            vectOfDnrFields, vectOfRcvFields);
  
   else if (varType == 2 || varType == 21) 
     setIBCTransfersCommonVar2(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread,
			      xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, 
            density, pressure, 
            vx, vy, vz, 
            utau, yplus,
            NULL, NULL, NULL, NULL, NULL,
            ipt_tmp, size,
            gamma, cv, muS, Cs, Ts, 0.71,
            vectOfDnrFields, vectOfRcvFields);

   else if (varType == 3 || varType == 31)
     setIBCTransfersCommonVar3(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread,
    			      xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, 
                density, pressure, 
                vx, vy, vz, 
                utau, yplus,
                NULL, NULL, NULL, NULL, NULL,
                ipt_tmp, size,
                gamma, cv, muS, Cs, Ts, 0.71,
                vectOfDnrFields, vectOfRcvFields);

 }//fin omp

  // sortie
  RELEASESHAREDB(resd, arrayD, fd, cnd); 
  BLOCKRELEASEMEMD;
  BLOCKRELEASEMEM2;
  return tpl;
}
//=============================================================================
/* Transfert de champs conservatifs sous forme de numpy 
   From zone- 
   Retourne une liste de numpy directement de champs interpoles */
//=============================================================================
PyObject* K_CONNECTOR::_setIBCTransfersD(PyObject* self, PyObject* args)
{
  PyObject *zoneD, *pyIndDonor, *pyArrayTypes, *pyArrayCoefs;
  PyObject *pyArrayXPC, *pyArrayXPI, *pyArrayXPW;
  PyObject *pyArrayYPC, *pyArrayYPI, *pyArrayYPW;
  PyObject *pyArrayZPC, *pyArrayZPI, *pyArrayZPW;
  PyObject *pyArrayDens, *pyArrayPressure;
  PyObject *pyArrayVx, *pyArrayVy, *pyArrayVz;
  PyObject *pyArrayUtau, *pyArrayYplus;
  PyObject* pyVariables;
  E_Int bctype, vartype, compact;
  E_Float gamma, cv, muS, Cs, Ts;
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;

  if (!PYPARSETUPLE(args,
                    "OOOOOOOOOOOOOOOOOOOOOllldddddsss", "OOOOOOOOOOOOOOOOOOOOOiiidddddsss",
                    "OOOOOOOOOOOOOOOOOOOOOlllfffffsss", "OOOOOOOOOOOOOOOOOOOOOiiifffffsss",
                    &zoneD, &pyVariables, &pyIndDonor, &pyArrayTypes, &pyArrayCoefs, 
                    &pyArrayXPC, &pyArrayYPC, &pyArrayZPC,
                    &pyArrayXPW, &pyArrayYPW, &pyArrayZPW,
                    &pyArrayXPI, &pyArrayYPI, &pyArrayZPI,
                    &pyArrayDens, &pyArrayPressure, 
                    &pyArrayVx, &pyArrayVy, &pyArrayVz, 
                    &pyArrayUtau, &pyArrayYplus,
                    &bctype, &vartype, &compact, &gamma, &cv, &muS, &Cs, &Ts,
                    &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters))
  {
      return NULL;
  }
  E_Int bcType  = E_Int(bctype);  // 0: wallslip; 1: noslip; 2: log law of wall; 3: Musker law of wall
  E_Int varType = E_Int(vartype); // 1:conservatives, 2:(ro,u,v,w,t), 3:(ro,u,v,w,p)

  vector<PyArrayObject*> hook;
  E_Int imdjmd, imd,jmd,kmd, cnNfldD, nvars,ndimdxR, ndimdxD,meshtype;
  E_Float* iptroD; E_Float* iptroR; 

# include "extract_interpD.h"
# include "IBC/extract_IBC.h"

  vector<E_Float*> fieldsD; vector<E_Int> posvarsD; 
  E_Int* ptrcnd;
  char* eltTypeD=NULL; char* varStringD=NULL;
  char* varStringOut = new char[K_ARRAY::VARSTRINGLENGTH];

  //codage general (lent ;-) )
  if (compact==0)
  {   
     /*---------------------------------------------*/
     /* Extraction des infos sur le domaine donneur */
     /*---------------------------------------------*/
     E_Int cnSizeD;
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

     // Extrait les positions des variables a transferer
     E_Int posvd;
     varStringOut[0] = '\0';

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
             if (posvd != -1) 
             { 
               posvarsD.push_back(posvd);
               if (varStringOut[0]=='\0') strcpy(varStringOut,varname);
               else {strcat(varStringOut,","); strcat(varStringOut,varname);}
             }
           }
#if PY_VERSION_HEX >= 0x03000000
           else if (PyUnicode_Check(tpl0)) 
           {
              char* varname = PyBytes_AsString(PyUnicode_AsUTF8String(tpl0));
              posvd = K_ARRAY::isNamePresent(varname, varStringD);      
              if (posvd != -1) 
              { 
                 posvarsD.push_back(posvd);
                if (varStringOut[0]=='\0') strcpy(varStringOut,varname);
                else {strcat(varStringOut,","); strcat(varStringOut,varname);}
              }
           }
#endif
           else 
             PyErr_Warn(PyExc_Warning, "_setIBCTransfersD: variable must be a string. Skipped.");
         }
       }
      }
     nvars = posvarsD.size();
  }
  else // les variables a transferes sont compactes: on recuperes uniquement la premiere et la taille
  {
# include "getfromzonecompactD.h"
  }

  PyObject* tpl = K_ARRAY::buildArray(nvars, varStringOut, nbRcvPts, 1, 1);
  E_Float*  frp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fieldROut(nbRcvPts, nvars, frp, true);
  //
  //utile la mise a zero???
  //
  fieldROut.setAllValuesAtNull();

  vector<E_Float*> vectOfRcvFields(nvars);
  vector<E_Float*> vectOfDnrFields(nvars);

  if (compact==0) 
  {
    for (E_Int eq = 0; eq < nvars; eq++)
      {
       vectOfRcvFields[eq] = fieldROut.begin(eq+1);
       vectOfDnrFields[eq] = fieldsD[ posvarsD[eq] ];
      }
  }
  else
  {
    for (E_Int eq = 0; eq < nvars; eq++)
      {
       vectOfRcvFields[eq] = fieldROut.begin(eq+1);
       vectOfDnrFields[eq] = iptroD + eq*ndimdxD;
      }
  }

////
////
//  Interpolation parallele
////  
////  

     // tableau temporaire pour utiliser la routine commune setIBCTransfersCommon
     FldArrayI rcvPtsI(nbRcvPts); E_Int*  rcvPts = rcvPtsI.begin();
#    include "commonInterpTransfers_direct.h"
  
    E_Int threadmax_sdm  = __NUMTHREADS__;

    E_Int size = (nbRcvPts/threadmax_sdm)+1; // on prend du gras pour gerer le residus
          size = size + size % 8;            // on rajoute du bas pour alignememnt 64bits
    if (bctype <=1 ) size = 0;               // tableau inutile

    FldArrayF  tmp(size*16*threadmax_sdm);
    E_Float* ipt_tmp = tmp.begin();


#pragma omp parallel default(shared)
  {

//     #pragma omp barrier
// barriere inutile car synchro implicit a la prochaine loop parallel
  #pragma omp for
  for (E_Int noind = 0; noind < nbRcvPts; noind++) rcvPts[noind] = noind;

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

   if (varType == 1 || varType == 11) 
     setIBCTransfersCommonVar1(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread, 
			      xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, 
            density, pressure, 
            vx, vy, vz, 
            utau, yplus,
            NULL, NULL, NULL, NULL, NULL,
            ipt_tmp, size,
            gamma, cv, muS, Cs, Ts, 0.71,
            vectOfDnrFields, vectOfRcvFields);
  
   else if (varType == 2 || varType == 21) 
     setIBCTransfersCommonVar2(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread,
			      xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, 
            density, pressure, 
            vx, vy, vz, 
            utau, yplus,
            NULL, NULL, NULL, NULL, NULL,
            ipt_tmp, size,
            gamma, cv, muS, Cs, Ts, 0.71,
            vectOfDnrFields, vectOfRcvFields);

   else if (varType == 3 || varType == 31)
     setIBCTransfersCommonVar3(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread,
    			      xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, 
                density, pressure, 
                vx, vy, vz, 
                utau, yplus,
                NULL, NULL, NULL, NULL, NULL,
                ipt_tmp, size,
                gamma, cv, muS, Cs, Ts, 0.71,
                vectOfDnrFields, vectOfRcvFields);

  } // Fin zone // omp

  // sortie
  delete [] varStringOut;
  RELEASESHAREDZ(hook, varStringD, eltTypeD);
  BLOCKRELEASEMEMD;
  BLOCKRELEASEMEM2;
  return tpl;
}
