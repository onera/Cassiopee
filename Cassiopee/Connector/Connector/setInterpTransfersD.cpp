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

#include "connector.h"
#include "param_solver.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Transfert de champs sous forme de numpy */
//=============================================================================
PyObject* K_CONNECTOR::setInterpTransfersD( PyObject* self, PyObject* args ) 
{
  PyObject *arrayD, *pyIndDonor, *pyArrayTypes, *pyArrayCoefs;
  if ( !PyArg_ParseTuple( args, "OOOO", &arrayD, &pyIndDonor, &pyArrayTypes, &pyArrayCoefs ) ) { return NULL; }

  /*---------------------------------------------*/
  /* Extraction des infos sur le domaine donneur */
  /*---------------------------------------------*/
  E_Int      imd, jmd, kmd, imdjmd;
  FldArrayF* fd;
  FldArrayI* cnd;
  char*      varStringD;
  char*      eltTypeD;
  E_Int      resd = K_ARRAY::getFromArray( arrayD, varStringD, fd, imd, jmd, kmd, cnd, eltTypeD, true );
  if (resd != 2 && resd != 1) 
    {
      PyErr_SetString( PyExc_TypeError, "setInterpTransfersD: 1st arg is not a valid array." );
      return NULL;
    }

  E_Int meshtype = resd;  // 1: structure, 2 non structure

#include "extract_interpD.h"

  if ( res_donor * res_type * res_coef == 0 ) 
    {
      RELEASESHAREDB( resd, arrayD, fd, cnd );
      if ( res_donor != 0 ) { RELEASESHAREDN( pyIndDonor, donorPtsI ); }
      if ( res_type != 0 ) { RELEASESHAREDN( pyArrayTypes, typesI ); }
      if ( res_coef != 0 ) { RELEASESHAREDN( pyArrayCoefs, donorCoefsF ); }
      PyErr_SetString( PyExc_TypeError,
		       "setInterpTransfersD: 2nd and 3rd arg must be a numpy of integers. 4th arg a numpy floats " );
      return NULL;
    }

  E_Int  nvars   = fd->getNfld( );  // nb de champs a interpoler
  E_Int* ptrcnd  = cnd->begin( );
  E_Int  cnNfldD = cnd->getNfld( );

  PyObject* tpl = K_ARRAY::buildArray( nvars, varStringD, nbRcvPts, 1, 1 );
  E_Float*  frp = K_ARRAY::getFieldPtr( tpl );
  FldArrayF fieldROut( nbRcvPts, nvars, frp, true );
  //
  // utile la mise a zero???
  //
  fieldROut.setAllValuesAtNull( );

  // Transferts
  // Types valides: 2, 3, 4, 5

  vector< E_Float* > vectOfRcvFields( nvars );
  vector< E_Float* > vectOfDnrFields( nvars );

  for ( E_Int eq = 0; eq < nvars; eq++ ) 
    {
      vectOfRcvFields[eq] = fieldROut.begin( eq + 1 );
      vectOfDnrFields[eq] = fd->begin( eq + 1 );
    }

  ////
  ////
  //  Interpolation parallele
  ////
  ////
#include "commonInterpTransfers_direct.h"

  // sortie
  RELEASESHAREDB( resd, arrayD, fd, cnd );
  RELEASESHAREDN( pyIndDonor, donorPtsI );
  RELEASESHAREDN( pyArrayTypes, typesI );
  RELEASESHAREDN( pyArrayCoefs, donorCoefsF );
  return tpl;
}

//=============================================================================
/* Transfert de champs sous forme de numpy
   From zone
   Retourne une liste de numpy directement des champs interpoles */
//=============================================================================
PyObject* K_CONNECTOR::_setInterpTransfersD( PyObject* self, PyObject* args )
{
  PyObject *zoneD, *pyIndDonor, *pyArrayTypes, *pyArrayCoefs;
  PyObject *pyVariables;
  char* GridCoordinates;
  char* FlowSolutionNodes;
  char* FlowSolutionCenters;
  char* cellNVariable;

  E_Int vartype, compact;
  if ( !PYPARSETUPLE_(args, OOOO_ O_ II_ SSSS_,
               &zoneD, &pyVariables, &pyIndDonor, &pyArrayTypes,
		       &pyArrayCoefs, &vartype, &compact, &cellNVariable,
		       &GridCoordinates, &FlowSolutionNodes, &FlowSolutionCenters) )
    {
      return NULL;
    }

  vector< PyArrayObject* > hook;
  E_Int imdjmd, imd, jmd, kmd, cnNfldD, nvars, ndimdxD=0, meshtype;
  E_Float * iptroD=NULL;
  E_Int poscd = -1;
  E_Int posvd, posvarcd=-1;

#include "extract_interpD.h"

  vector< E_Float* > fieldsD;
  vector< E_Int >    posvarsD;
  E_Int* ptrcnd;
  char* eltTypeD, *varStringD;
  char* varStringOut = new char[K_ARRAY::VARSTRINGLENGTH];
  // codage general (lent ;-) )
  if (compact == 0) 
    {
      //---------------------------------------------/
      // Extraction des infos sur le domaine donneur /
      //---------------------------------------------/
      E_Int cnSizeD;
      vector<E_Int> locsD;
      vector<E_Int*> cnd;
      meshtype = K_PYTREE::getFromZone( zoneD, 0, 0, varStringD, fieldsD, locsD,
					imd, jmd, kmd, cnd, cnSizeD, cnNfldD,
					eltTypeD, hook, GridCoordinates,
					FlowSolutionNodes, FlowSolutionCenters);
      if (cnd.size()>0) ptrcnd = cnd[0];
      E_Int nfld = fieldsD.size();
      // Extrait les positions des variables a transferer
      E_Int initAll   = false;
      varStringOut[0] = '\0';
      poscd = K_ARRAY::isNamePresent(cellNVariable, varStringD);
      if (PyList_Check(pyVariables) != 0)
        {
	  int nvariables = PyList_Size(pyVariables);
	  if (nvariables > 0) 
            {
	      for (int i = 0; i < nvariables; i++ )
                {
		  PyObject* tpl0 = PyList_GetItem( pyVariables, i );
		  if (PyString_Check(tpl0))
                    {
		      char* varname = PyString_AsString(tpl0);
		      posvd = K_ARRAY::isNamePresent(varname, varStringD);
		      if (posvd == poscd) posvarcd = posvd;
		      if (posvd != -1)
                        {
			  posvarsD.push_back(posvd);
			  if (varStringOut[0] == '\0') strcpy(varStringOut, varname);
			  else
                            {
			      strcat(varStringOut, ",");
			      strcat(varStringOut, varname);
                            }
                        }
                    }
#if PY_VERSION_HEX >= 0x03000000
		  else if (PyUnicode_Check(tpl0))
                    {
		      const char* varname = PyUnicode_AsUTF8(tpl0);
		      posvd = K_ARRAY::isNamePresent(varname, varStringD);
		      if (posvd == poscd) posvarcd = posvd;
		      if (posvd != -1)
                        {
			  posvarsD.push_back(posvd);
			  if (varStringOut[0] == '\0' ) strcpy(varStringOut, varname);
			  else
                            {
			      strcat(varStringOut, ",");
			      strcat(varStringOut, varname);
                            }
                        }
                    }
#endif
		  else
		    PyErr_Warn(PyExc_Warning, "_setInterpTransfersD: variable must be a string. Skipped.");
                }
            }// fin nvariables > 0
	  else  // si [] toutes les variables communes sont transferees sauf le cellN
            {
	      initAll = true;
            }
        }//fin pyVariables==list
      else  // si None toutes les variables communes sont transferees sauf le cellN
        {
	  initAll = true;
        }
      if (initAll == true) //all common variables are transfered
        {
	  vector<char*> tmpvars;
	  K_ARRAY::extractVars(varStringD, tmpvars);
	  for (E_Int i = 0; i < nfld; i++)
            {
	      posvarsD.push_back(i);
	      if (varStringOut[0] == '\0')
		strcpy( varStringOut, tmpvars[i] );
	      else
                {
		  strcat(varStringOut, ",");
		  strcat(varStringOut, tmpvars[i]);
                }
            }
	  for (size_t i = 0; i < tmpvars.size(); i++) delete[] tmpvars[i];
        }// fin initAll==true

      delete[] varStringD;
      delete[] eltTypeD;
      nvars = posvarsD.size();
    }// fin compact = 0
  else  // les variables a transferer sont compactees: on recupere uniquement la premiere et la taille
    {
#include "getfromzonecompactD.h"
    }

  // -- no check (perfo) --
  PyObject* tpl = K_ARRAY::buildArray(nvars, varStringOut, nbRcvPts, 1, 1 );
  E_Float*  frp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fieldROut(nbRcvPts, nvars, frp, true);
  //
  // utile la mise a zero???
  //
  fieldROut.setAllValuesAtNull();

  // Transferts
  // Types valides: 2, 3, 4, 5
  vector< E_Float* > vectOfRcvFields(nvars);
  vector< E_Float* > vectOfDnrFields(nvars);

  /* compact = 0
     si on a specifie la variable cellNVariable, interpolation specifique comme un champ cellN
  */
  if (compact == 0)
    {
      for (E_Int eq = 0; eq < nvars; eq++)
        {
	  vectOfRcvFields[eq] = fieldROut.begin(eq+1);
	  vectOfDnrFields[eq] = fieldsD[posvarsD[eq]];
        }
      //interpolation of all fields
#include "commonInterpTransfers_direct.h"

      //transfer of cellN variable
      if (posvarcd > -1)
        {
	  E_Int indR, type, nocf;
	  E_Int indD0, indD, i, j, k, ncfLoc;
	  E_Int noi = 0; // compteur sur le tableau d indices donneur
	  E_Int sizecoefs = 0;
	  E_Float* cellNR = fieldROut.begin(nvars);
	  E_Float* cellND = fieldsD[poscd];
	  for (E_Int noind = 0; noind < nbRcvPts; noind++)
            {
	      // adressage direct pour indR
	      indR = noind;
# include "commonCellNTransfersStrict.h"   
            }  
        }
    }
  else
    {
      for (E_Int eq = 0; eq < nvars; eq++)
        {
	  vectOfRcvFields[eq] = fieldROut.begin(eq+1);
	  vectOfDnrFields[eq] = iptroD + eq * ndimdxD;
        }
      ///
      ////
      //  Interpolation parallele
      ////
      ////
#include "commonInterpTransfers_direct.h"  
    }

  delete[] varStringOut;
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
  RELEASESHAREDN(pyIndDonor, donorPtsI);
  RELEASESHAREDN(pyArrayTypes, typesI);
  RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
  return tpl;
}

//=============================================================================
// Transfert de champs sous forme de numpy
// From zone
// Retourne une liste de numpy directement des champs interpoles
// in place + from zone + tc compact
//=============================================================================
PyObject* K_CONNECTOR::__setInterpTransfersD(PyObject* self, PyObject* args)
{
  PyObject *zonesR, *zonesD;
  PyObject* pyVariables, *pydtloc;
  PyObject *pyParam_int, *pyParam_real;
  E_Int     vartype, type_transfert, no_transfert, It_target;
  E_Int     nstep, nitmax, rk, exploc, num_passage, rank, isWireModel, isIbmMoving;

  if ( !PYPARSETUPLE_(args, OOOO_ OO_ IIII_ IIII_ IIII_,
		      &zonesR,
		      &zonesD, &pyVariables, &pydtloc, &pyParam_int, &pyParam_real, &It_target, &vartype, &type_transfert,
		      &no_transfert, &nstep, &nitmax, &rk, &exploc, &num_passage, &rank, &isWireModel, &isIbmMoving) ) 
  {
    return NULL;
  }

  //zoneR : zones arbre t
  //zoneD : zones arbre tc
  //param_int/real : arbre tc
  E_Float gamma, cv, muS, Cs, Ts, Pr;
  E_Int it_target= E_Int(It_target);

  E_Int isWireModel_local = 0;
  if (isWireModel == -1)
  {
    isWireModel_local=-1;
    isWireModel=0;
  }

  // gestion nombre de pass pour ID et/ou IBC
  E_Int TypeTransfert = E_Int( type_transfert );
  E_Int pass_deb, pass_fin;
  if     ( TypeTransfert==0) { pass_deb = 1; pass_fin = 2; }  // ID
  else if( TypeTransfert==1) { pass_deb = 0; pass_fin = 1; }  // IBCD
  else                       { pass_deb = 0; pass_fin = 2; }  // ALL

  E_Int NoTransfert = E_Int( no_transfert );

  vector< PyArrayObject* > hook;

  E_Int   kmd, cnNfldD, nvars, meshtype, nvars_Pnt2;

  /* varType determines the variables to transfer
     varType :
       1  : conservatives,
       11 : conservatives + ronutildeSA
       2  : (ro,u,v,w,t)
       21 : (ro,u,v,w,t) + ronutildeSA
       3  : (ro,u,v,w,p)
       31 : (ro,u,v,w,p) + ronutildeSA
       4  : (Q1,..., QN)   LBM
       5  : Couplage NS LBM
  */
  if     ( vartype <= 3 &&  vartype >= 1) nvars =5;  // NSLaminar transfer, only 5 variables
  else if( vartype == 4 ) nvars =32;    // LBM transfer, 19 or 27 Qs and 5 macros (32 max in total)
  else if (vartype == 41) nvars =44; //38;    // LBM Overset : 19 or 27 Q, 5 macro and 6 gradients
  else if( vartype == 5 ) nvars =30;    // Hybrid NSLBM transfer, 5 macros, 19 Qs + 6 gradients
  else                    nvars =6;     // NSTurbulent transfer, 6 variables

  //printf("nvarD %d \n ", nvars);

  E_Int nidomD = PyList_Size( zonesD );

  // pointeur pour stocker solution au centre ET au noeud
  E_Int*    ipt_ndimdxD;
  E_Int**   ipt_cnd;
  E_Float** ipt_roD;
  E_Float** ipt_roD_vert;
  E_Float** ipt_roD_Pnt2;
  E_Float** ipt_qD;
  E_Float** ipt_qD_vert;
  E_Float** ipt_SD;
  E_Float** ipt_SD_vert;
  E_Float** ipt_psiGD;
  E_Float** ipt_psiGD_vert;

  ipt_ndimdxD  = new E_Int[nidomD * 8];  // on stocke ndimdx, imd, jmd, en centre et vertexe, meshtype et cnDfld
  ipt_cnd      = new E_Int*[nidomD];

  ipt_roD      = new E_Float*[nidomD * 9];//1
  ipt_roD_vert   = ipt_roD      + nidomD; //2
  ipt_roD_Pnt2   = ipt_roD_vert + nidomD; //3
  ipt_qD         = ipt_roD_Pnt2 + nidomD; //4
  ipt_qD_vert    = ipt_qD       + nidomD; //5
  ipt_SD         = ipt_qD_vert  + nidomD; //6
  ipt_SD_vert    = ipt_SD       + nidomD; //7
  ipt_psiGD      = ipt_SD_vert  + nidomD; //8
  ipt_psiGD_vert = ipt_psiGD    + nidomD; //9
  /*----------------------------------*/
  /* Get  param_int param_real zone donneuse */
  /*----------------------------------*/

  E_Int nidomR = PyList_Size( zonesR );

  E_Int** ipt_param_intD; E_Float** ipt_param_realD;

  ipt_param_intD   = new E_Int*[nidomR];
  ipt_param_realD  = new E_Float*[nidomR];

  for (E_Int nd = 0; nd < nidomR; nd++)
    {

      //attention en parallel, zone receveuse indisponible, on recupere param_int de la zone donneuse
      PyObject* zoneR = PyList_GetItem(zonesR, nd);

      PyObject* own      = K_PYTREE::getNodeFromName1(zoneR , ".Solver#ownData");
      PyObject* t        = K_PYTREE::getNodeFromName1(own, "Parameter_int");
      ipt_param_intD[nd] = K_PYTREE::getValueAI(t, hook);

      t          = K_PYTREE::getNodeFromName1(own, "Parameter_real");
      ipt_param_realD[nd]= K_PYTREE::getValueAF(t, hook);
    }

  FldArrayI* dtloc;
  E_Int res_donor = K_NUMPY::getFromNumpyArray(pydtloc, dtloc, true);
  E_Int* iptdtloc = dtloc->begin();

  //------------------------------------/
  // Extraction tableau int et real     /
  //------------------------------------/
  FldArrayI* param_int;
  res_donor     = K_NUMPY::getFromNumpyArray( pyParam_int, param_int, true );
  E_Int*     ipt_param_int = param_int->begin( );
  FldArrayF* param_real;
  res_donor               = K_NUMPY::getFromNumpyArray( pyParam_real, param_real, true );
  E_Float* ipt_param_real = param_real->begin( );

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
      if     (ivar==0) {varname  = vartmp;}
      else if(ivar==1) {varname1 = vartmp;}  //En LBM, on a besoin d'echanger les macro !!ET!! les Q donc deux varname
      else if(ivar==2) {varname2 = vartmp;}  //Pour le couplage NS-LBM, on a besoin d'echanger les Q !!ET!! les macro !!ET!! les gradients
      else if(ivar==3) {varname3 = vartmp;}
      else {printf("Warning: souci varname setInterpTransfersD \n"); }
    }

  // on recupere sol et solcenter ainsi que connectivite et taille zones Donneuses (tc)
  for ( E_Int nd = 0; nd < nidomD; nd++ ) {
    PyObject* zoneD = PyList_GetItem( zonesD, nd );
   //printf("avt get fromzone %d %d %d \n ", nidomD, vartype,nd);
#include "getfromzoneDcompact_all.h"
  }


  E_Int sizecomIBC = ipt_param_int[2];
  E_Int sizecomID  = ipt_param_int[3+sizecomIBC];
  E_Int shift_graph = sizecomIBC + sizecomID + 3;

  E_Int threadmax_sdm = __NUMTHREADS__;
  E_Int ech           = ipt_param_int[NoTransfert + shift_graph];
  E_Int nrac          = ipt_param_int[ech + 1];  // nb total de raccord
  E_Int nrac_inst     = ipt_param_int[ech + 2];  // nb total de raccord instationnaire
  E_Int timelevel     = ipt_param_int[ech + 3];  // nb de pas de temps stocker pour chaque raccord instationnaire
  E_Int nrac_steady  = nrac - nrac_inst;                 //nb total de raccord stationnaire


  //gestion nombre de pass pour raccord instationnaire
  E_Int pass_inst_deb=0;
  E_Int pass_inst_fin=1;
  E_Int nrac_inst_level = 0;
  if (nrac_inst > 0) {
    pass_inst_fin=2;
    nrac_inst_level = ipt_param_int[ech + 4 + it_target + timelevel] - ipt_param_int[ech + 4 + it_target] + 1;
  }
  // on dimension tableau travail pour IBC et pour transfert  FldArrayI* dtloc;
  // E_Int nrac_inst_level = ipt_param_int[ech + 4 + it_target + timelevel] - ipt_param_int[ech + 4 + it_target] + 1;
  char* varStringOut    = new char[K_ARRAY::VARSTRINGLENGTH];
  char* varStringOut1   = new char[K_ARRAY::VARSTRINGLENGTH];
  char* varStringOut2   = new char[K_ARRAY::VARSTRINGLENGTH];

  PyObject** list_tpl; PyObject** list_tpl1; PyObject** list_tpl2;
  list_tpl  = new PyObject*[nrac_steady + nrac_inst_level];
  list_tpl1 = new PyObject*[nrac_steady + nrac_inst_level];
  list_tpl2 = new PyObject*[nrac_steady + nrac_inst_level];

  E_Float** frp; E_Float** frp1; E_Float** frp2;
  frp  = new E_Float*[nrac_steady + nrac_inst_level];
  frp1 = new E_Float*[nrac_steady + nrac_inst_level];
  frp2 = new E_Float*[nrac_steady + nrac_inst_level];

  //tableau pour optimisation transfert implicit local. Si zone donneuse non modifiee a la ssiter nstep, on fait rien.
  E_Int impli_local[nidomD];
  E_Int nssiter  = iptdtloc[0];
  E_Int shift_omp= iptdtloc[11];
  E_Int* ipt_omp = iptdtloc + shift_omp;
  E_Int nbtask = ipt_omp[nstep-1]; 
  E_Int ptiter = ipt_omp[nssiter+ nstep-1];

  E_Int impli_local_init=0;
  if (isWireModel_local==-1 || isIbmMoving==1){impli_local_init=1;}
    
  for (E_Int nd = 0; nd < nidomR; nd++) {impli_local[nd]=impli_local_init;}
  for (E_Int ntask = 0; ntask < nbtask; ntask++)
    {
      E_Int pttask = ptiter + ntask*(6+threadmax_sdm*7);
      E_Int nd = ipt_omp[ pttask ];
      impli_local[nd]=1;
    }

  E_Int size_autorisation = nrac_steady+1;
  size_autorisation = K_FUNC::E_max(  size_autorisation , nrac_inst+1);

  E_Int autorisation_transferts[pass_inst_fin][size_autorisation];

  // 1ere pass_inst: les raccord fixe
  // 2eme pass_inst: les raccord instationnaire
  E_Int nbRcvPts_mx = 0;
  E_Int ibcTypeMax  = 0;

  E_Int count_rac   = 0;
  for  (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
    {
      E_Int irac_deb = 0;
      E_Int irac_fin = nrac_steady;

      if ( pass_inst == 1 )
        {
         irac_deb = ipt_param_int[ech + 4 + it_target];
         irac_fin = ipt_param_int[ech + 4 + it_target + timelevel];
        }

      for ( E_Int irac = irac_deb; irac < irac_fin; irac++ ) {
	E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
	E_Int nbRcvPts  = ipt_param_int[shift_rac + nrac * 10 + 1];

	if ( nbRcvPts > nbRcvPts_mx ) nbRcvPts_mx = nbRcvPts;

	if( ipt_param_int[shift_rac+nrac*3] > ibcTypeMax)  ibcTypeMax =  ipt_param_int[shift_rac+nrac*3];

	E_Int ibcType = ipt_param_int[shift_rac + nrac * 3];
	E_Int ibc = 1;
	if ( ibcType < 0) ibc = 0;

	E_Int irac_auto= irac-irac_deb;
	autorisation_transferts[pass_inst][irac_auto]=0;

	if(exploc == 1)// Si on est en explicit local, on va autoriser les transferts entre certaines zones seulement en fonction de la ss-ite courante

	  {
	    E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;

	    E_Int levelD = ipt_param_int[debut_rac + 25];
	    E_Int levelR = ipt_param_int[debut_rac + 24];
	    E_Int cyclD  = nitmax/levelD;

	    // Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse
	    if (levelD > levelR and num_passage == 1)
	      {
		if (nstep%cyclD==cyclD-1 || (nstep%cyclD==cyclD/2 && (nstep/cyclD)%2==1))
		  { autorisation_transferts[pass_inst][irac_auto]=1; }
		else {continue;}
	      }
	    // Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
	    else if (levelD < levelR and num_passage == 1)
	      {
		if (nstep%cyclD==1 or nstep%cyclD==cyclD/4 or nstep%cyclD== cyclD/2-1 or nstep%cyclD== cyclD/2+1 or nstep%cyclD== cyclD/2+cyclD/4 or nstep%cyclD== cyclD-1)
		  { autorisation_transferts[pass_inst][irac_auto]=1; }
		else {continue;}
	      }
	    // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse
	    else if (levelD == levelR and num_passage == 1)
	      {
		if (nstep%cyclD==cyclD/2-1 or (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==0) or nstep%cyclD==cyclD-1)
		  { autorisation_transferts[pass_inst][irac_auto]=1; }
		else {continue;}
	      }
	    // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse (cas du deuxieme passage)
	    else if (levelD ==  ipt_param_int[debut_rac +24] and num_passage == 2)
	      {
		if (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==1) { autorisation_transferts[pass_inst][irac_auto]=1; }
		else {continue;}
	      }
	    else {continue;}

	  }
	// Sinon, on autorise les transferts si zone donneuse modifier a l'iteration courante
	else {
	  E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ];
	  if (impli_local[NoD]==1) {autorisation_transferts[pass_inst][irac_auto]=1;}
          else {continue;}
	}

	// QUOI????? - CBX
	//on skippe les racc IBC  si passe ID  et viceversa: IM
	if      ( TypeTransfert == 0 && ibc == 1 ) { continue;}
	else if ( TypeTransfert == 1 && ibc == 0 ) { continue;}

	E_Int nvars_loc = ipt_param_int[shift_rac + nrac * 13 + 1];  //nbre variable a transferer pour rans/LES

	E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ];
        E_Int overset  =  ipt_param_intD[NoD][LBM_OVERSET];        //flag pour overset en LBM recuper� sur param_int donneuse

        if      (nvars_loc==19 && overset==0) nvars_loc = nvars_loc + 5;
        else if (nvars_loc==19 && overset==1) nvars_loc = nvars_loc + 5 + 6 + 6;

        // COUPLAGE NS LBM - Recupere les solveurs des zones R et D
        E_Int solver_D=2; E_Int solver_R=2;
        if (nvars_loc == 11) {solver_R =4;}
        if (nvars_loc == -5) {solver_D =4; nvars_loc = 5;}


        E_Int nbRcvPts_loc = nbRcvPts;
        if ( nvars_loc == 5 ) 
          {
           if ( strcmp( varname, "Density" ) == 0 )
             {
	      if (isWireModel==0){ strcpy( varStringOut, "Density,VelocityX,VelocityY,VelocityZ,Temperature" ); }
              else               { strcpy( varStringOut, "Density_WM,VelocityX_WM,VelocityY_WM,VelocityZ_WM,Temperature_WM" ); }
             }
           else if ( strcmp( varname, "Density_P1" ) == 0 )
             {
              if (isWireModel==0){ strcpy( varStringOut, "Density_P1,VelocityX_P1,VelocityY_P1,VelocityZ_P1,Temperature_P1" );}
              else{                strcpy( varStringOut, "Density_WM,VelocityX_WM,VelocityY_WM,VelocityZ_WM,Temperature_WM" ); }
             }
           else{ printf("souci setinterpTransferD " SF_D_ " \n", nvars_loc); }
          }
        else if ( nvars_loc == 6 ) 
          {
           if ( strcmp( varname, "Density" ) == 0 )
             {
              if (isWireModel==0){ strcpy( varStringOut, "Density,VelocityX,VelocityY,VelocityZ,Temperature,TurbulentSANuTilde" ); }
              else{                strcpy( varStringOut, "Density_WM,VelocityX_WM,VelocityY_WM,VelocityZ_WM,Temperature_WM,TurbulentSANuTilde_WM" ); }
             }
           else if ( strcmp( varname, "Density_P1" ) == 0 )
             {
              if (isWireModel==0){ strcpy( varStringOut, "Density_P1,VelocityX_P1,VelocityY_P1,VelocityZ_P1,Temperature_P1,TurbulentSANuTilde_P1" );}
              else{                strcpy( varStringOut, "Density_WM,VelocityX_WM,VelocityY_WM,VelocityZ_WM,Temperature_WM,TurbulentSANuTilde_WM" ); }
             }
           else{ printf("souci setinterpTransferD " SF_D_ " \n", nvars_loc); }
          }
        else if ( nvars_loc == 19 ||  nvars_loc == 24 ) 
          { 
            if      ( strcmp( varname, "Density"    ) == 0 ){ strcpy( varStringOut, "Density,VelocityX,VelocityY,VelocityZ,Temperature" ); }
            else if ( strcmp( varname, "Density_P1" ) == 0 ){ strcpy( varStringOut, "Density_P1,VelocityX_P1,VelocityY_P1,VelocityZ_P1,Temperature_P1" );}
            else{ printf("souci setinterpTransferD " SF_D_ " \n", nvars_loc); }

            if ( varname1 != NULL ) // Raccords LBM
              {
                if      ( strcmp( varname1, "Q1")    == 0 ){ strcpy( varStringOut1, "Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16,Q17,Q18,Q19" );}
                else if ( strcmp( varname1, "Q1_M1") == 0 ){ strcpy( varStringOut1, "Q1_M1,Q2_M1,Q3_M1,Q4_M1,Q5_M1,Q6_M1,Q7_M1,Q8_M1,Q9_M1,Q10_M1,Q11_M1,Q12_M1,Q13_M1,Q14_M1,Q15_M1,Q16_M1,Q17_M1,Q18_M1,Q19_M1" );}
                else{ printf("souci setinterpTransferD " SF_D_ " \n", nvars_loc); }
              }
          } 
        else if ( nvars_loc == 27 ||  nvars_loc == 32 ) 
          { 
            if      ( strcmp( varname, "Density"    ) == 0 ){ strcpy( varStringOut, "Density,VelocityX,VelocityY,VelocityZ,Temperature" ); }
            else if ( strcmp( varname, "Density_P1" ) == 0 ){ strcpy( varStringOut, "Density_P1,VelocityX_P1,VelocityY_P1,VelocityZ_P1,Temperature_P1" );}
            else{ printf("souci setinterpTransferD " SF_D_ " \n", nvars_loc); }

            if ( varname1 != NULL ) // Raccords LBM
              {
                if      ( strcmp( varname1, "Q1")    == 0 ){ strcpy( varStringOut1, "Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16,Q17,Q18,Q19,Q20,Q21,Q22,Q23,Q24,Q25,Q26,Q27" );}
                else if ( strcmp( varname1, "Q1_M1") == 0 ){ strcpy( varStringOut1, "Q1_M1,Q2_M1,Q3_M1,Q4_M1,Q5_M1,Q6_M1,Q7_M1,Q8_M1,Q9_M1,Q10_M1,Q11_M1,Q12_M1,Q13_M1,Q14_M1,Q15_M1,Q16_M1,Q17_M1,Q18_M1,Q19_M1,Q20_M1,Q21_M1,Q22_M1,Q23_M1,Q24_M1,Q25_M1,Q26_M1,Q27_M1" );}
                else{ printf("souci setinterpTransferD " SF_D_ " \n", nvars_loc); }
              }
          }
        else if ( nvars_loc == 11 ) 
          { 
            if      ( strcmp( varname, "Density"    ) == 0 ){ strcpy( varStringOut, "Density,VelocityX,VelocityY,VelocityZ,Temperature" ); }
            else if ( strcmp( varname, "Density_P1" ) == 0 ){ strcpy( varStringOut, "Density_P1,VelocityX_P1,VelocityY_P1,VelocityZ_P1,Temperature_P1" );}
            else{ printf("souci setinterpTransferD " SF_D_ " \n", nvars_loc); }

            if ( varname2 != NULL ) // Raccords LBM
              { strcpy( varStringOut2, "Sxx,Sxy,Sxz,Syy,Syz,Szz"); }
          }
        else{ printf("souci aiguillage setinterpTransferD " SF_D_ " \n", nvars_loc); }

        if ( nvars_loc == 5 || nvars_loc == 6 )
          {
            list_tpl[count_rac]  = K_ARRAY::buildArray( nvars_loc, varStringOut, nbRcvPts_loc, 1, 1 );
            frp[count_rac]       = K_ARRAY::getFieldPtr( list_tpl[count_rac] );
	    if(impli_local[NoD]==0) { frp[count_rac][0] =0.;}
          }
        else if ( nvars_loc == 19 || nvars_loc == 24 )
          {
            // TABLEAU QUI VA CONTENIR LES MACROS DONC DE TAILLE 5 (ou 6 si turb)
            list_tpl[count_rac]  = K_ARRAY::buildArray( 5, varStringOut, nbRcvPts_loc, 1, 1 );
            frp[count_rac]       = K_ARRAY::getFieldPtr( list_tpl[count_rac] );
	    if(impli_local[NoD]==0) { frp[count_rac][0] =0.;}
            // TABLEAU QUI VA CONTENIR LES Q DONC DE TAILLE 19
            list_tpl1[count_rac] = K_ARRAY::buildArray( 19, varStringOut1, nbRcvPts_loc, 1, 1 );
            frp1[count_rac]      = K_ARRAY::getFieldPtr( list_tpl1[count_rac] );
	    if(impli_local[NoD]==0) { frp1[count_rac][0] =0.;}
          }
        else if ( nvars_loc == 11 )
          {
            // TABLEAU QUI VA CONTENIR LES MACROS DONC DE TAILLE 5 (ou 6 si turb)
            list_tpl[count_rac]  = K_ARRAY::buildArray( 5, varStringOut, nbRcvPts_loc, 1, 1 );
            frp[count_rac]       = K_ARRAY::getFieldPtr( list_tpl[count_rac] );
	    if(impli_local[NoD]==0) { frp[count_rac][0] =0.;}
            // TABLEAU QUI VA CONTENIR LES GRADIENTS DONC DE TAILLE 6
            list_tpl2[count_rac] = K_ARRAY::buildArray( 6, varStringOut2, nbRcvPts_loc, 1, 1 );
            frp2[count_rac]      = K_ARRAY::getFieldPtr( list_tpl2[count_rac] );
	    if(impli_local[NoD]==0) { frp2[count_rac][0] =0.;}
          }

        count_rac += 1;


      }  // racs
    }      // pass steady et unsteady




  E_Int size              = ( nbRcvPts_mx / threadmax_sdm ) + 1;  // on prend du gras pour gerer le residus
  E_Int r                 = size % 8;
  if ( r != 0 ) size      = size + 8 - r;  // on rajoute du bas pour alignememnt 64bits
  if ( ibcTypeMax <= 1 ) size = 0;             // tableau inutile : SP ; voir avec Ivan

  FldArrayF tmp( size * 17 * threadmax_sdm );
  E_Float*  ipt_tmp = tmp.begin();

  // tableau temporaire pour utiliser la routine commune setIBCTransfersCommon
  FldArrayI rcvPtsI( nbRcvPts_mx );
  E_Int*    rcvPts = rcvPtsI.begin();

  E_Float** RcvFields = new E_Float*[ nvars*threadmax_sdm];
  E_Float** DnrFields = new E_Float*[ nvars*threadmax_sdm];


  //# pragma omp parallel default(shared)  num_threads(1)
#pragma omp parallel default(shared)
  {
#ifdef _OPENMP
    E_Int ithread           = omp_get_thread_num( ) + 1;
    E_Int Nbre_thread_actif = omp_get_num_threads( );  // nombre de thread actif dans cette zone
#else
    E_Int ithread           = 1;
    E_Int Nbre_thread_actif = 1;
#endif

    E_Int type;
    E_Int indD0, indD, i, j, k, ncfLoc, indCoef, noi, sizecoefs, imd, jmd, imdjmd;

    E_Float** vectOfRcvFields = RcvFields + nvars*(ithread-1);
    E_Float** vectOfDnrFields = DnrFields + nvars*(ithread-1);

    // 1ere pass: IBC
    // 2eme pass: transfert
    //
    for ( E_Int ipass_typ = pass_deb; ipass_typ < pass_fin; ipass_typ++ )
    {
      // 1ere pass_inst: les raccord fixe
      // 2eme pass_inst: les raccord instationnaire
	    E_Int count_rac = 0;

      for ( E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
	    {
	    E_Int irac_deb = 0;
	    E_Int irac_fin = nrac_steady;
	    if ( pass_inst == 1 )
	      {
		irac_deb = ipt_param_int[ech + 4 + it_target];
		irac_fin = ipt_param_int[ech + 4 + it_target + timelevel];
	      }
	    for ( E_Int irac = irac_deb; irac < irac_fin; irac++ )
	      {
		E_Int shift_rac = ech + 4 + timelevel * 2 + irac; 

		E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5     ];
		E_Int irac_auto= irac-irac_deb;

                //printf("verif %d %d irac/nrac= %d %d  past_int=%d , pastType= %d \n", autorisation_transferts[pass_inst][irac_auto] , impli_local[NoD], irac, irac_fin, pass_inst, ipass_typ );

		if (autorisation_transferts[pass_inst][irac_auto]==1)
		  {
           
                    //on skippe les racc IBC  si passe ID  et viceversa: IM
                    E_Int ibcType = ipt_param_int[shift_rac+nrac*3];
                    E_Int ibc = 1;
                    if ( ibcType < 0 ) ibc = 0;
                    if ( 1-ibc != ipass_typ ) continue;

                     //printf("ipass= %d, irac= %d, ibc=  %d, envoie vers: %d, size_rac= %d \n", ipass, irac, ibc,
                     //ipt_param_int[ ech ], ipt_param_int[ shift_rac + nrac*10 +1 ]);

                    E_Int NoD       = ipt_param_int[shift_rac + nrac * 5     ];
                    E_Int loc       = ipt_param_int[shift_rac + nrac * 9  + 1];  //+1 a cause du nrac mpi
                    E_Int nbRcvPts  = ipt_param_int[shift_rac + nrac * 10 + 1];
                    E_Int nvars_loc = ipt_param_int[shift_rac + nrac * 13 + 1];  // neq fonction raccord rans/LES
                    E_Int rotation  = ipt_param_int[shift_rac + nrac * 14 + 1];  // flag pour periodicite azymuthal

                    // COUPLAGE NS LBM - Recupere les solveurs des zones R et D
                    E_Int solver_D=2; E_Int solver_R=2;
                    if (nvars_loc == 11) {solver_R =4;}
                    if (nvars_loc == -5) {solver_D =4; nvars_loc = 5;}
                    if (nvars_loc == 19) {solver_D =4; solver_R=4;}

                    E_Int overset  =  ipt_param_intD[NoD][LBM_OVERSET];        //flag pour overset en LBM
                    if      (nvars_loc==19 && overset==0) nvars_loc = nvars_loc + 5;
                    else if (nvars_loc==19 && overset==1) nvars_loc = nvars_loc + 5 + 6 + 6;

                    E_Int  meshtype = ipt_ndimdxD[NoD + nidomD * 6];
                    E_Int  cnNfldD  = ipt_ndimdxD[NoD + nidomD * 7];
                    E_Int* ptrcnd   = ipt_cnd[NoD];

                    // printf("navr_loc %d %d %d \n", nvars_loc, nvars, Rans);

                    if ( loc == 0 ) {
                      printf("transferts optimises pas code en vertex \n");
                      /*for ( E_Int eq = 0; eq < nvars_loc; eq++ ) {
                      vectOfRcvFields[eq] = frp[count_rac] + eq * nbRcvPts;
                      vectOfDnrFields[eq] = ipt_roD_vert[NoD] + eq * ipt_ndimdxD[NoD + nidomD * 3];
                    }*/
                    imd = ipt_ndimdxD[NoD + nidomD * 4];
                    jmd = ipt_ndimdxD[NoD + nidomD * 5];
                    } else {

                        /*--------------------------------------------------------------------*/
                        /*                GESTION DES TRANSFERTS EN CENTER                    */
                        /* 2 cas: - si pas de couplage NSLBM, transferts habituels,           */
                        /*        - si couplage NSLBM, raccords a adapter                     */
                        /*--------------------------------------------------------------------*/
                        // for ( E_Int eq = 0; eq < nvars_loc; eq++ )
                        // {
                        //     vectOfRcvFields[eq] = frp[count_rac] + eq * nbRcvPts;
                        //     vectOfDnrFields[eq] = ipt_roD[NoD]   + eq * ipt_param_intD[NoD][ NDIMDX ];
                        // }
                        if (nvars_loc == 5 || nvars_loc == 6) // Transferts NS classiques ou LBM -> NS
                        {
                          for (E_Int eq = 0; eq < nvars_loc; eq++)
                          {
                             vectOfRcvFields[eq] = frp[count_rac] + eq * nbRcvPts;
                             vectOfDnrFields[eq] = ipt_roD[NoD]   + eq * ipt_param_intD[NoD][ NDIMDX ];
                          }
                        }
                        else if (nvars_loc == 24 || nvars_loc == 32) // Transferts LBM classiques
                        {
                          for (E_Int eq = 0; eq < nvars_loc; eq++)
                          {
                             if (eq < 5) // On commence par copier les 5 variables macros
                             {
                                vectOfRcvFields[eq] = frp[count_rac] + eq * nbRcvPts;
                                vectOfDnrFields[eq] = ipt_roD[NoD]   + eq * ipt_param_intD[NoD][ NDIMDX ];
                             }
                             else // Puis on copie les fonctions de distribution
                             {
                                vectOfRcvFields[eq] = frp1[count_rac] + (eq-5) * nbRcvPts;
                                vectOfDnrFields[eq] = ipt_qD[NoD]     + (eq-5) * ipt_param_intD[NoD][ NDIMDX ];
                             }
                          }
                        }
                        else if (nvars_loc == 11)  //Transfert NS -> LBM
                        {
                          for (E_Int eq = 0; eq < nvars_loc; eq++)
                          {
                             if (eq < 5) // On commence par copier les 5 variables macros
                             {
                                vectOfRcvFields[eq] = frp[count_rac] + eq * nbRcvPts;
                                vectOfDnrFields[eq] = ipt_roD[ NoD]  + eq * ipt_param_intD[NoD][ NDIMDX ];
                             }
                             // else // Puis on copie les gradients
                             // {
                             //    vectOfRcvFields[eq] = frp2[count_rac] + (eq-5) * nbRcvPts;
                             //    vectOfDnrFields[eq] = ipt_SD[NoD]     + (eq-5) * ipt_param_intD[NoD][ NDIMDX ];
                             // }
                          }

                        }

		      imd = ipt_ndimdxD[NoD + nidomD];
		      jmd = ipt_ndimdxD[NoD + nidomD * 2];
                    }

                    imdjmd = imd * jmd;

                    ////
                    //  Interpolation parallele
                    ////
                    ////
                    E_Int pos;
                    pos               = ipt_param_int[shift_rac + nrac * 7];
                    E_Int* ntype      = ipt_param_int + pos;
                    pos               = pos + 1 + ntype[0];
                    E_Int* types      = ipt_param_int + pos;
                    pos               = ipt_param_int[shift_rac + nrac * 6];
                    E_Int* donorPts   = ipt_param_int + pos;
                    pos               = ipt_param_int[shift_rac + nrac * 8];
                    E_Float* ptrCoefs = ipt_param_real + pos;

                    E_Int    nbInterpD = ipt_param_int[shift_rac + nrac];
                    E_Float* xPC       = NULL;
                    E_Float* xPI       = NULL;
                    E_Float* xPW       = NULL;
                    E_Float* densPtr   = NULL;
                    if ( ibc == 1 ) {
		      xPC     = ptrCoefs + nbInterpD;
		      xPI     = ptrCoefs + nbInterpD + 3 * nbRcvPts;
		      xPW     = ptrCoefs + nbInterpD + 6 * nbRcvPts;
		      densPtr = ptrCoefs + nbInterpD + 9 * nbRcvPts;
                    }

                    E_Int ideb        = 0;
                    E_Int ifin        = 0;
                    E_Int shiftCoef   = 0;
                    E_Int shiftDonor = 0;

                    for ( E_Int ndtyp = 0; ndtyp < ntype[0]; ndtyp++ ) {
		      type = types[ifin];

		      SIZECF( type, meshtype, sizecoefs );

		      ifin = ifin + ntype[1 + ndtyp];

		      E_Int pt_deb, pt_fin;

		      /// oldschool
		      // Calcul du nombre de champs a traiter par chaque thread
		      E_Int size_bc = ifin - ideb;
		      E_Int chunk   = size_bc / Nbre_thread_actif;
		      E_Int r       = size_bc - chunk * Nbre_thread_actif;
		      // pts traitees par thread
		      if ( ithread <= r ) {
			pt_deb = ideb + ( ithread - 1 ) * ( chunk + 1 );
			pt_fin = pt_deb + ( chunk + 1 );
		      } else {
			pt_deb = ideb + ( chunk + 1 ) * r + ( ithread - r - 1 ) * chunk;
			pt_fin = pt_deb + chunk;
		      }

		      // Si type 0, calcul sequentiel
		      if ( type == 0 ) {
			if ( ithread == 1 ) {
			  pt_deb = ideb;
			  pt_fin = ifin;
			} else {
			  pt_deb = ideb;
			  pt_fin = ideb;
			}
		      }
		      noi     = shiftDonor;  // compteur sur le tableau d indices donneur
		      indCoef = ( pt_deb - ideb ) * sizecoefs + shiftCoef;

		      E_Int NoR = ipt_param_int[shift_rac + nrac * 11 + 1]; 
			
		      if (isWireModel==1){
			if (ibcType==141){
			  // Interpolation at these points only for _WM variables
			  // Can be included in the commonInterpTransfers_reorder_*eq.h functions
			  // but do not want an if statement inside the loops
#           include "commonInterpTransfersD_reorder_WireModel.h"
			}
			else{
#           include "commonInterpTransfersD_reorder_WireModelZero.h"
			}
			// Only do the Interpolations for ibcType=141 (Wire Model)
			// No other IBC interpolations
			goto IBCType_SkipInterpolation; //looking into this
		      }
		      else{
			if (ibcType==141){
#           include "commonInterpTransfersD_reorder_WireModelZero.h"
			  goto IBCType_SkipInterpolation;
			}
		      }


		      if ( nvars_loc == 5 || (ibc==1 && solver_R==4) )  {
#include "commonInterpTransfersD_reorder_5eq.h"
		      } else if ( nvars_loc == 6 ) {
#include "commonInterpTransfersD_reorder_6eq.h"
		      } else if ( nvars_loc ==19 ) {
#include "LBM/commonInterpTransfersD_reorder_19eq.h"
		      } else {
#include "LBM/commonInterpTransfersD_reorder_neq.h"
		      }

		    IBCType_SkipInterpolation:

                        // COUPLAGE NS-LBM: changement d'unite
                        if (solver_D==4 && solver_R<4)
                        {
                         // Transfert LBM vers NS: repasse dans unites SI
#                        include "includeTransfersD_dimLBMtoNS.h"
                        }
                        else if (solver_D<4 && solver_R==4)
                        {
                         // Transfert NS vers LBM : adimensionnement
#                        include "includeTransfersD_dimNStoLBM.h"
                        }

		      // Prise en compte de la periodicite par rotation
		      if ( rotation == 1 ) {
			E_Float* angle = ptrCoefs + nbInterpD;
#include "includeTransfersD_rotation.h"
		      }

		      // ibc
		      if (ibc == 1 && isWireModel == 0)
                        {
			  //E_Int nvars = vectOfDnrFields.size();

			  Pr    = ipt_param_realD[ NoD ][ PRANDT ];
			  Ts    = ipt_param_realD[ NoD ][ TEMP0 ];
			  Cs    = ipt_param_realD[ NoD ][ CS ];
			  muS   = ipt_param_realD[ NoD ][ XMUL0 ];
			  cv    = ipt_param_realD[ NoD ][ CVINF ];
			  gamma = ipt_param_realD[ NoD ][ GAMMA ];

			  // tableau temporaire pour utiliser la routine commune setIBCTransfersCommon
			  for ( E_Int noind = pt_deb; noind < pt_fin; noind++ ) rcvPts[noind] = noind;

			  setIBCTransfersCommonVar2( ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread, xPC,
						     xPC + nbRcvPts, xPC + nbRcvPts * 2, xPW, xPW + nbRcvPts,
						     xPW + nbRcvPts * 2, xPI, xPI + nbRcvPts, xPI + nbRcvPts * 2,
						     densPtr, 
						     ipt_tmp, size, nvars, ipt_param_realD[NoD], 
						     vectOfDnrFields, vectOfRcvFields ); 

                         if (solver_R==4)
                         {
#                          include "includeTransfersD_LBM_feq.h"
                         }
                      
                        }  // ibc

		      //        } //chunk
		      ideb        =  ideb + ntype[ 1 + ndtyp];
		      shiftCoef   = shiftCoef + ntype[1 + ndtyp] * sizecoefs;  // shift coef   entre 2 types successif
		      shiftDonor = shiftDonor + ntype[1 + ndtyp];            // shift donor entre 2 types successif
                    }                                                            // type

                count_rac += 1;
//#pragma omp barrier
// barriere inutile en MPI: le receveur est local a chaque raccord
//
		  } // autorisation transfert
	      }  // irac
	  }      // pass_inst
#pragma omp barrier
      }  // pass typ (id/ibc)
  }      // omp

  PyObject* infos = PyList_New( 0 );

  count_rac = 0;
  for (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
    {
      E_Int irac_deb = 0;
      E_Int irac_fin = nrac_steady;
      if ( pass_inst == 1 ) 
        { irac_deb = ipt_param_int[ech + 4 + it_target];
     	  irac_fin = ipt_param_int[ech + 4 + it_target + timelevel];
        }

      for ( E_Int irac = irac_deb; irac < irac_fin; irac++ )
        {

	  E_Int irac_auto= irac-irac_deb;
	  if (autorisation_transferts[pass_inst][irac_auto]==1)
            {
	      E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
	      E_Int nbRcvPts  = ipt_param_int[shift_rac + nrac * 10 + 1];

	      E_Int nbRcvPts_loc = nbRcvPts;
	      E_Int NoD          = ipt_param_int[ shift_rac + nrac*5 ];

	      E_Int ibcType = ipt_param_int[shift_rac + nrac * 3];
	      E_Int ibc = 1; 
	      if (ibcType < 0) ibc = 0;
	      if      ( TypeTransfert == 0 && ibc == 1 ) { continue; } 
	      else if ( TypeTransfert == 1 && ibc == 0 ) { continue; }

              E_Int nvars_loc = ipt_param_int[shift_rac + nrac * 13 + 1];
              if (nvars_loc==19) nvars_loc = 24;

              if ( nvars_loc == 5 || nvars_loc == 6 )
              {
	         PyObject* info = PyList_New( 0 );

	         PyObject* Nozone;
	         Nozone = PyInt_FromLong( ipt_param_int[shift_rac + nrac * 11 + 1] );
	         PyList_Append( info, Nozone );  // No Zone receuveuse

                 PyList_Append( info, list_tpl[count_rac] );
                 Py_DECREF( list_tpl[count_rac] );  // tableau data

	         E_Int     PtlistDonor = ipt_param_int[shift_rac + nrac * 12 + 1];
	         E_Int*    ipt_listRcv = ipt_param_int + PtlistDonor;
	         PyObject* listRcv     = K_NUMPY::buildNumpyArray( ipt_listRcv, nbRcvPts_loc, 1, 1 );

	         PyList_Append( info, listRcv );
	         Py_DECREF( listRcv );  // ListReceveur

	         PyList_Append( infos, info );
	         Py_DECREF( info );
              }
              else if ( nvars_loc == 19 || nvars_loc == 24)
              {
                 // ON DOIT GERER DEUX JEUX DE VARIABLES
                   PyObject* Nozone;
                   Nozone = PyInt_FromLong( ipt_param_int[shift_rac + nrac * 11 + 1] );
                   E_Int     PtlistDonor = ipt_param_int[shift_rac + nrac * 12 + 1];
                   E_Int*    ipt_listRcv = ipt_param_int + PtlistDonor;
                   PyObject* listRcv     = K_NUMPY::buildNumpyArray( ipt_listRcv, nbRcvPts_loc, 1, 1 );

                   // 1. Premier jeu de variables : MACROS
                     PyObject* info = PyList_New( 0 );
                     PyList_Append( info, Nozone );  // No Zone receuveuse
                     PyList_Append( info, list_tpl[count_rac] );
                     Py_DECREF( list_tpl[count_rac] );  // tableau data
                     PyList_Append( info, listRcv );
                     PyList_Append( infos, info );
                     Py_DECREF( info );
                   // 2. Premier jeu de variables : Q
                     info = PyList_New( 0 );
                     PyList_Append( info, Nozone );  // No Zone receuveuse
                     PyList_Append( info, list_tpl1[count_rac] );
                     Py_DECREF( list_tpl1[count_rac] );  // tableau data
                     PyList_Append( info, listRcv );
                     PyList_Append( infos, info );
                     Py_DECREF( info );

                   Py_DECREF( listRcv );  // ListReceveur
              }
              else if ( nvars_loc == 11 )
              {
                   // ON DOIT GERER DEUX JEUX DE VARIABLES
                   PyObject* Nozone;
                   Nozone = PyInt_FromLong( ipt_param_int[shift_rac + nrac * 11 + 1] );
                   E_Int     PtlistDonor = ipt_param_int[shift_rac + nrac * 12 + 1];
                   E_Int*    ipt_listRcv = ipt_param_int + PtlistDonor;
                   PyObject* listRcv     = K_NUMPY::buildNumpyArray( ipt_listRcv, nbRcvPts_loc, 1, 1 );

                   // 1. Premier jeu de variables : MACROS
                     PyObject* info = PyList_New( 0 );
                     PyList_Append( info, Nozone );  // No Zone receuveuse
                     PyList_Append( info, list_tpl[count_rac] );
                     Py_DECREF( list_tpl[count_rac] );  // tableau data
                     PyList_Append( info, listRcv );
                     PyList_Append( infos, info );
                     Py_DECREF( info );
                   // 2. Premier jeu de variables : Gradients
                     info = PyList_New( 0 );
                     PyList_Append( info, Nozone );  // No Zone receuveuse
                     PyList_Append( info, list_tpl2[count_rac] );
                     Py_DECREF( list_tpl2[count_rac] );  // tableau data
                     PyList_Append( info, listRcv );
                     PyList_Append( infos, info );
                     Py_DECREF( info );

                   Py_DECREF( listRcv );  // ListReceveur
              }
	      count_rac += 1;
	    }//autorisation transfert
	}// irac
    }// pass_inst

 
  delete[] ipt_param_intD; delete[] ipt_param_realD;
  delete [] RcvFields;  delete [] DnrFields;
  delete[] frp;
  delete[] frp1;
  delete[] frp2;
  delete[] list_tpl;
  delete[] list_tpl1;
  delete[] list_tpl2;
  delete[] ipt_ndimdxD;
  delete[] ipt_roD;
  delete[] ipt_cnd;
  delete[] varStringOut; delete[] varStringOut1; delete[] varStringOut2;
  RELEASESHAREDZ( hook, (char*)NULL, (char*)NULL );
  RELEASESHAREDN(pydtloc        , dtloc        );
  RELEASESHAREDN( pyParam_int, param_int );
  RELEASESHAREDN( pyParam_real, param_real );

  return infos;

}

//=============================================================================
// Transfert de champs sous forme de numpy
// From zone
// Retourne une liste de numpy directement des champs interpoles
// in place + from zone + tc compact + gradP
//=============================================================================
PyObject* K_CONNECTOR::__setInterpTransfersD4GradP(PyObject* self, PyObject* args)
{
  PyObject *zonesR, *zonesD;
  PyObject *pyVariables;
  PyObject *pyParam_int, *pyParam_real;
  E_Int     vartype, type_transfert, no_transfert, It_target;
  E_Int     nstep, nitmax, rk, exploc, num_passage;
  E_Float   gamma, cv, muS, Cs, Ts, Pr;

  if ( !PYPARSETUPLE_( args, OOOO_ O_ IIII_ IIII_ I_,
          &zonesR, &zonesD, &pyVariables, &pyParam_int, &pyParam_real, &It_target, &vartype, 
          &type_transfert, &no_transfert, &nstep, &nitmax, &rk, &exploc, &num_passage)) 
  {
    return NULL;
  }

  E_Int it_target= E_Int(It_target);

  E_Float alpha = 0.;

  // gestion nombre de pass pour ID et/ou IBC
  E_Int TypeTransfert = E_Int( type_transfert );
  E_Int pass_deb, pass_fin;
  if ( TypeTransfert == 0 ) 
  {
    pass_deb = 1;
    pass_fin = 2;
  }  // ID
  else if ( TypeTransfert == 1 ) 
  {
    pass_deb = 0;
    pass_fin = 1;
  }  // IBCD
  else {
    pass_deb = 0;
    pass_fin = 2;
  }  // ALL

  E_Int NoTransfert = E_Int( no_transfert );

  vector< PyArrayObject* > hook;

  // E_Int kmd, cnNfldD, nvars,ndimdxR, ndimdxD,meshtype;
  E_Int   kmd, cnNfldD, nvars, ndimdxD, meshtype, nvars_grad;
  //E_Float* iptroD;

  if     ( vartype <= 3 &&  vartype >= 1) nvars =5;
  else if( vartype == 4 ) nvars =27;    // on majore pour la LBM, car nvar sert uniquememnt a dimensionner taille vector
  else                    nvars =6;

  if (vartype == 22 || vartype == 23 || vartype == 24)
    {
      nvars = 6;
      nvars_grad = 6;
    }

  E_Int nidomD = PyList_Size( zonesD );

  // pointeur pour stocker solution au centre ET au noeud
  E_Int*    ipt_ndimdxD;
  E_Int**   ipt_cnd;
  E_Float** ipt_roD;
  E_Float** ipt_roD_vert;
  E_Float** ipt_gradD;

  ipt_ndimdxD  = new E_Int[nidomD * 8];  // on stocke ndimdx, imd, jmd, en centre et vertexe, meshtype et cnDfld
  ipt_cnd      = new E_Int*[nidomD];
  ipt_roD      = new E_Float*[nidomD * 3];
  ipt_roD_vert = ipt_roD + nidomD;
  ipt_gradD    = ipt_roD_vert + nidomD;

  /*----------------------------------*/
  /* Get  param_int param_real zone donneuse */
  /*----------------------------------*/

  E_Int nidomR = PyList_Size( zonesR );

  E_Int** ipt_param_intD; E_Float** ipt_param_realD;

  ipt_param_intD   = new E_Int*[nidomR];
  ipt_param_realD  = new E_Float*[nidomR];

  for (E_Int nd = 0; nd < nidomR; nd++)
    {

      //attention en parallel, zone receveuse indisponible, on recupere param_int de la zone donneuse
      PyObject* zoneR = PyList_GetItem(zonesR, nd);

      PyObject* own      = K_PYTREE::getNodeFromName1(zoneR , ".Solver#ownData");
      PyObject* t        = K_PYTREE::getNodeFromName1(own, "Parameter_int");
      ipt_param_intD[nd] = K_PYTREE::getValueAI(t, hook);

      t          = K_PYTREE::getNodeFromName1(own, "Parameter_real");
      ipt_param_realD[nd]= K_PYTREE::getValueAF(t, hook);
    }


  //------------------------------------/
  // Extraction tableau int et real     /
  //------------------------------------/
  FldArrayI* param_int;
  E_Int      res_donor     = K_NUMPY::getFromNumpyArray( pyParam_int, param_int, true );
  E_Int*     ipt_param_int = param_int->begin( );
  FldArrayF* param_real;
  res_donor               = K_NUMPY::getFromNumpyArray( pyParam_real, param_real, true );
  E_Float* ipt_param_real = param_real->begin( );

  // On recupere le nom de la 1ere variable a recuperer
  PyObject* tpl0    = PyList_GetItem( pyVariables, 0 );
  char*     varname = NULL;
  if (PyString_Check(tpl0)) varname = PyString_AsString(tpl0);
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(tpl0)) varname = (char*)PyUnicode_AsUTF8(tpl0);
#endif
  // on recupere sol et solcenter ainsi que connectivite et taille zones Donneuses (tc)
  for ( E_Int nd = 0; nd < nidomD; nd++ ) {
    PyObject* zoneD = PyList_GetItem( zonesD, nd );
#include "getfromzoneDcompact_allGradP.h"
  }

  E_Int nbcomIBC = ipt_param_int[2];
  E_Int nbcomID  = ipt_param_int[3+nbcomIBC];

  E_Int shift_graph = nbcomIBC + nbcomID + 3;

  E_Int threadmax_sdm = __NUMTHREADS__;
  E_Int ech           = ipt_param_int[NoTransfert + shift_graph];
  E_Int nrac          = ipt_param_int[ech + 1];  // nb total de raccord
  E_Int nrac_inst     = ipt_param_int[ech + 2];  // nb total de raccord instationnaire
  E_Int timelevel     = ipt_param_int[ech + 3];  // nb de pas de temps stocker pour chaque raccord instationnaire
  E_Int nrac_steady  = nrac - nrac_inst;                 //nb total de raccord stationnaire


  //gestion nombre de pass pour raccord instationnaire
  E_Int pass_inst_deb=0;
  E_Int pass_inst_fin=1;
  E_Int nrac_inst_level = 0;
  if (nrac_inst > 0) {
    pass_inst_fin=2;
    nrac_inst_level = ipt_param_int[ech + 4 + it_target + timelevel] - ipt_param_int[ech + 4 + it_target] + 1;
  }
  // on dimension tableau travail pour IBC et pour transfert
  // E_Int nrac_inst_level = ipt_param_int[ech + 4 + it_target + timelevel] - ipt_param_int[ech + 4 + it_target] + 1;
  char* varStringOut         = new char[K_ARRAY::VARSTRINGLENGTH];
  char* varStringOutGrad     = new char[K_ARRAY::VARSTRINGLENGTH];
  PyObject** list_tpl;
  list_tpl = new PyObject*[nrac_steady + nrac_inst_level];

  PyObject** list_tplGrad;
  list_tplGrad = new PyObject*[nrac_steady + nrac_inst_level];

  E_Float** frp;
  frp = new E_Float*[nrac_steady + nrac_inst_level];

  E_Float** frpGrad;
  frpGrad = new E_Float*[nrac_steady + nrac_inst_level];

  E_Int size_autorisation = nrac_steady+1;
  size_autorisation = K_FUNC::E_max(  size_autorisation , nrac_inst+1);

  E_Int autorisation_transferts[pass_inst_fin][size_autorisation];

  // 1ere pass_inst: les raccord fixe
  // 2eme pass_inst: les raccord instationnaire
  E_Int nbRcvPts_mx = 0;
  E_Int ibcTypeMax  = 0;
  E_Int count_rac   = 0;




  for  (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
    {
      E_Int irac_deb = 0;
      E_Int irac_fin = nrac_steady;

      if ( pass_inst == 1 )
        {
	  irac_deb = ipt_param_int[ech + 4 + it_target];
	  irac_fin = ipt_param_int[ech + 4 + it_target + timelevel];
        }

      for ( E_Int irac = irac_deb; irac < irac_fin; irac++ ) {
	E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
	E_Int nbRcvPts  = ipt_param_int[shift_rac + nrac * 10 + 1];

	if ( nbRcvPts > nbRcvPts_mx ) nbRcvPts_mx = nbRcvPts;

	if( ipt_param_int[shift_rac+nrac*3] > ibcTypeMax)  ibcTypeMax =  ipt_param_int[shift_rac+nrac*3];

	E_Int ibcType = ipt_param_int[shift_rac + nrac * 3];
	E_Int ibc = 1;
	if ( ibcType < 0) ibc = 0;

	E_Int irac_auto= irac-irac_deb;
	autorisation_transferts[pass_inst][irac_auto]=0;

        if(exploc == 1) // if(rk==3 and exploc == 2) // Si on est en explicit local, on va autoriser les transferts entre certaines zones seulement en fonction de la ss-ite courante

          {

	    E_Int debut_rac = ech + 4 + timelevel*2 + 1 + nrac*16 + 27*irac;

	    E_Int levelD = ipt_param_int[debut_rac + 25];
	    E_Int levelR = ipt_param_int[debut_rac + 24];
	    E_Int cyclD  = nitmax/levelD;

	    // Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse
	    if (levelD > levelR and num_passage == 1)
	      {
		if (nstep%cyclD==cyclD-1 or nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==1)
		  { autorisation_transferts[pass_inst][irac_auto]=1; }
		else {continue;}
	      }
	    // Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
	    else if (levelD < levelR and num_passage == 1)
	      {
		if (nstep%cyclD==1 or nstep%cyclD==cyclD/4 or nstep%cyclD== cyclD/2-1 or nstep%cyclD== cyclD/2+1 or nstep%cyclD== cyclD/2+cyclD/4 or nstep%cyclD== cyclD-1)
		  { autorisation_transferts[pass_inst][irac_auto]=1; }
		else {continue;}
	      }
	    // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse
	    else if (levelD == levelR and num_passage == 1)
	      {
		if (nstep%cyclD==cyclD/2-1 or (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==0) or nstep%cyclD==cyclD-1)
		  { autorisation_transferts[pass_inst][irac_auto]=1; }
		else {continue;}
	      }
	    // Le pas de temps de la zone donneuse est egal a celui de la zone receveuse (cas du deuxieme passage)
	    else if (levelD ==  ipt_param_int[debut_rac +24] and num_passage == 2)
	      {
		if (nstep%cyclD==cyclD/2 and (nstep/cyclD)%2==1) { autorisation_transferts[pass_inst][irac_auto]=1; }
		else {continue;}
	      }
	    else {continue;}

          }
	// Sinon, on autorise les transferts entre ttes les zones a ttes les ss-ite
        else
          { autorisation_transferts[pass_inst][irac_auto]=1; }

	// QUOI????? - CBX
	if ( TypeTransfert == 0 && ibc == 1 ) {
	  continue;
	} else if ( TypeTransfert == 1 && ibc == 0 ) {
	  continue;
	}

	E_Int nvars_loc            = ipt_param_int[shift_rac + nrac * 13 + 1];  //nbre variable a transferer pour rans/LES

	if ( strcmp( varname, "Density" ) == 0 ) {
	  if ( nvars_loc == 5 )
	    strcpy( varStringOut, "Density,VelocityX,VelocityY,VelocityZ,Temperature" );
	  else
	    strcpy( varStringOut, "Density,VelocityX,VelocityY,VelocityZ,Temperature,TurbulentSANuTilde" );
	}
	else if ( strcmp( varname, "Density_P1" ) == 0 ) {
	  if ( nvars_loc == 5 )
	    strcpy( varStringOut, "Density_P1,VelocityX_P1,VelocityY_P1,VelocityZ_P1,Temperature_P1" );
	  else
	    strcpy( varStringOut,
		    "Density_P1,VelocityX_P1,VelocityY_P1,VelocityZ_P1,Temperature_P1,TurbulentSANuTilde_P1" );
	}
	else { //LBM
	  if ( nvars_loc == 9 )
	    strcpy( varStringOut, "Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9" );
	  else if ( nvars_loc == 19 )
	    strcpy( varStringOut, "Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16,Q17,Q18,Q19" );
	  else if ( nvars_loc == 27 )
	    strcpy( varStringOut, "Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,Q12,Q13,Q14,Q15,Q16,Q17,Q18,Q19,Q20,Q21,Q22,Q23,Q24,Q25,Q26,Q27" );
	}

	strcpy( varStringOutGrad, "gradxDensity,gradyDensity,gradzDensity,gradxTemperature,gradyTemperature,gradzTemperature" );

	list_tpl[count_rac] = K_ARRAY::buildArray( nvars_loc, varStringOut, nbRcvPts, 1, 1 );
	frp[count_rac]      = K_ARRAY::getFieldPtr( list_tpl[count_rac] );

	list_tplGrad[count_rac] = K_ARRAY::buildArray( nvars_grad, varStringOutGrad, nbRcvPts, 1, 1 );
	frpGrad[count_rac]      = K_ARRAY::getFieldPtr( list_tplGrad[count_rac] );

	count_rac += 1;


      }  // racs
    }      // pass steady et unsteady



  E_Int size              = ( nbRcvPts_mx / threadmax_sdm ) + 1;  // on prend du gras pour gerer le residus
  E_Int r                 = size % 8;
  if ( r != 0 ) size      = size + 8 - r;  // on rajoute du bas pour alignememnt 64bits
  if ( ibcTypeMax <= 1 ) size = 0;             // tableau inutile : SP ; voir avec Ivan

  FldArrayF tmp( size * 17 * threadmax_sdm );
  E_Float*  ipt_tmp = tmp.begin();

  // tableau temporaire pour utiliser la routine commune setIBCTransfersCommon
  FldArrayI rcvPtsI( nbRcvPts_mx );
  E_Int*    rcvPts = rcvPtsI.begin();



  //# pragma omp parallel default(shared)  num_threads(1)
#pragma omp parallel default(shared)
  {
#ifdef _OPENMP
    E_Int ithread           = omp_get_thread_num( ) + 1;
    E_Int Nbre_thread_actif = omp_get_num_threads( );  // nombre de thread actif dans cette zone
#else
    E_Int ithread           = 1;
    E_Int Nbre_thread_actif = 1;
#endif

    E_Int type;
    E_Int indD0, indD, i, j, k, ncfLoc, indCoef, noi, sizecoefs, imd, jmd, imdjmd;

    vector< E_Float* > vectOfRcvFields( nvars );
    vector< E_Float* > vectOfDnrFields( nvars );

    vector< E_Float* > vectOfGradRcvFields( nvars_grad );
    vector< E_Float* > vectOfGradDnrFields( nvars_grad );

    // 1ere pass: IBC
    // 2eme pass: transfert
    //
    for ( E_Int ipass_typ = pass_deb; ipass_typ < pass_fin; ipass_typ++ )
      {
	// 1ere pass_inst: les raccord fixe
	// 2eme pass_inst: les raccord instationnaire
	E_Int count_rac = 0;

	for ( E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
	  {
	    E_Int irac_deb = 0;
	    E_Int irac_fin = nrac_steady;
	    if ( pass_inst == 1 )
	      {
		irac_deb = ipt_param_int[ech + 4 + it_target];
		irac_fin = ipt_param_int[ech + 4 + it_target + timelevel];
	      }
	    for ( E_Int irac = irac_deb; irac < irac_fin; irac++ )
	      {
		E_Int irac_auto= irac-irac_deb;
		if (autorisation_transferts[pass_inst][irac_auto]==1)
		  {
                    E_Int shift_rac = ech + 4 + timelevel * 2 + irac;

                    E_Int ibcType = ipt_param_int[shift_rac+nrac*3];


                    E_Int ibc = 1;
                    if ( ibcType < 0 ) ibc = 0;


                    // printf("ipass= %d, irac= %d, ibc=  %d, envoie vers: %d, size_rac= %d \n", ipass, irac, ibc,
                    // ipt_param_int[ ech ], ipt_param_int[ shift_rac + nrac*10 +1 ]);
                    if ( 1-ibc!= ipass_typ ) continue;

                    E_Int NoD       = ipt_param_int[shift_rac + nrac * 5     ];
                    E_Int loc       = ipt_param_int[shift_rac + nrac * 9  + 1];  //+1 a cause du nrac mpi
                    E_Int nbRcvPts  = ipt_param_int[shift_rac + nrac * 10 + 1];
                    E_Int nvars_loc = ipt_param_int[shift_rac + nrac * 13 + 1];  // neq fonction raccord rans/LES
                    E_Int rotation  = ipt_param_int[shift_rac + nrac * 14 + 1];  // flag pour periodicite azymuthal


                    E_Int  meshtype = ipt_ndimdxD[NoD + nidomD * 6];
                    E_Int  cnNfldD  = ipt_ndimdxD[NoD + nidomD * 7];
                    E_Int* ptrcnd   = ipt_cnd[NoD];

                    // printf("navr_loc %d %d %d \n", nvars_loc, nvars, Rans);

                    if ( loc == 0 ) {
		      printf("transferts optimises pas code en vertex \n");
		      /*for ( E_Int eq = 0; eq < nvars_loc; eq++ ) {
			vectOfRcvFields[eq] = frp[count_rac] + eq * nbRcvPts;
			vectOfDnrFields[eq] = ipt_roD_vert[NoD] + eq * ipt_ndimdxD[NoD + nidomD * 3];
                        }*/
		      imd = ipt_ndimdxD[NoD + nidomD * 4];
		      jmd = ipt_ndimdxD[NoD + nidomD * 5];
                    } else {

		      for ( E_Int eq = 0; eq < nvars_loc; eq++ ) {
			vectOfRcvFields[eq] = frp[count_rac] + eq * nbRcvPts;
			vectOfDnrFields[eq] = ipt_roD[NoD] + eq *ipt_param_intD[NoD][ NDIMDX ];
		      }

		      imd = ipt_ndimdxD[NoD + nidomD];
		      jmd = ipt_ndimdxD[NoD + nidomD * 2];

		      for (E_Int eq = 0; eq < nvars_grad; eq++)
			{
			  vectOfGradRcvFields[eq] = frpGrad[count_rac] + eq * nbRcvPts;
			  vectOfGradDnrFields[eq] = ipt_gradD[NoD] + eq*ipt_param_intD[NoD][ NDIMDX ];
			}
                    }

                    imdjmd = imd * jmd;

                    ////
                    //  Interpolation parallele
                    ////
                    ////
                    E_Int pos;
                    pos               = ipt_param_int[shift_rac + nrac * 7];
                    E_Int* ntype      = ipt_param_int + pos;
                    pos               = pos + 1 + ntype[0];
                    E_Int* types      = ipt_param_int + pos;
                    pos               = ipt_param_int[shift_rac + nrac * 6];
                    E_Int* donorPts   = ipt_param_int + pos;
                    pos               = ipt_param_int[shift_rac + nrac * 8];
                    E_Float* ptrCoefs = ipt_param_real + pos;

                    E_Int    nbInterpD = ipt_param_int[shift_rac + nrac];
                    E_Float* xPC       = NULL;
                    E_Float* xPI       = NULL;
                    E_Float* xPW       = NULL;
                    E_Float* densPtr   = NULL;
                    if ( ibc == 1 ) 
                    {
                      xPC     = ptrCoefs + nbInterpD;
                      xPI     = ptrCoefs + nbInterpD + 3 * nbRcvPts;
                      xPW     = ptrCoefs + nbInterpD + 6 * nbRcvPts;
                      densPtr = ptrCoefs + nbInterpD + 9 * nbRcvPts;
                    }

                    E_Int ideb        = 0;
                    E_Int ifin        = 0;
                    E_Int shiftCoef   = 0;
                    E_Int shiftDonor = 0;

                    for ( E_Int ndtyp = 0; ndtyp < ntype[0]; ndtyp++ ) {
		      type = types[ifin];

		      SIZECF( type, meshtype, sizecoefs );

		      ifin = ifin + ntype[1 + ndtyp];

		      E_Int pt_deb, pt_fin;

		      /// oldschool
		      // Calcul du nombre de champs a traiter par chaque thread
		      E_Int size_bc = ifin - ideb;
		      E_Int chunk   = size_bc / Nbre_thread_actif;
		      E_Int r       = size_bc - chunk * Nbre_thread_actif;
		      // pts traitees par thread
		      if ( ithread <= r ) {
			pt_deb = ideb + ( ithread - 1 ) * ( chunk + 1 );
			pt_fin = pt_deb + ( chunk + 1 );
		      } else {
			pt_deb = ideb + ( chunk + 1 ) * r + ( ithread - r - 1 ) * chunk;
			pt_fin = pt_deb + chunk;
		      }

		      // Si type 0, calcul sequentiel
		      if ( type == 0 ) {
			if ( ithread == 1 ) {
			  pt_deb = ideb;
			  pt_fin = ifin;
			} else {
			  pt_deb = ideb;
			  pt_fin = ideb;
			}
		      }

		      noi     = shiftDonor;  // compteur sur le tableau d indices donneur
		      indCoef = ( pt_deb - ideb ) * sizecoefs + shiftCoef;

		      //E_Int NoR = ipt_param_int[shift_rac + nrac * 11 + 1];



		      //if (ipt_param_int[ech]==0) printf("No rac= %d , NoR= %d, NoD= %d, Ntype= %d, ptdeb= %d, ptfin= %d, NptD= %d, neq= %d, skip= %d, rank= %d, dest= %d,  thread= %d\n",
		      //irac, NoR,NoD, ntype[ 1 + ndtyp],pt_deb,pt_fin ,
		      //ipt_param_int[ shift_rac + nrac*10+1  ], ipt_param_int[ shift_rac + nrac*13+1  ], ipt_param_int[ shift_rac + nrac*15+1  ],
		      //rank, ipt_param_int[ ech  ], ithread );


		      if ( nvars_loc == 5 ) {
#include "commonInterpTransfersD_reorder_5eq.h"
		      } else if ( nvars_loc == 6 ) {
#include "commonInterpTransfersD_reorder_6eq.h"
		      } else if ( nvars_loc ==19 ) {
#include "commonInterpTransfersD_reorder_19eq.h"
		      } else {
#include "commonInterpTransfersD_reorder_neq.h"
		      }

		      noi       = shiftDonor;                             
		      indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
#include "commonInterpTransfersD_reorder_6eq_gradP.h"

		      // Prise en compte de la periodicite par rotation
		      if ( rotation == 1 ) {
			E_Float* angle = ptrCoefs + nbInterpD;
#include "includeTransfersD_rotation.h"
		      }

		      // ibc
		      if (ibc == 1)
                        {
			  //E_Int nvars = vectOfDnrFields.size();

			  alpha    = ipt_param_realD[ NoD ][ ALPHAGRADP ];

			  if (vartype == 23){
			    alpha = 1.;
			  }

			  Pr    = ipt_param_realD[ NoD ][ PRANDT ];
			  Ts    = ipt_param_realD[ NoD ][ TEMP0 ];
			  Cs    = ipt_param_realD[ NoD ][ CS ];
			  muS   = ipt_param_realD[ NoD ][ XMUL0 ];
			  cv    = ipt_param_realD[ NoD ][ CVINF ];
			  gamma = ipt_param_realD[ NoD ][ GAMMA ];

			  // tableau temporaire pour utiliser la routine commune setIBCTransfersCommon
			  for ( E_Int noind = pt_deb; noind < pt_fin; noind++ ) rcvPts[noind] = noind;

			  if (vartype == 22 || vartype == 23 || vartype == 24){
			    E_Float cvgam = cv*(gamma-1.);

			    for (E_Int noind = 0; noind < pt_fin-pt_deb; noind++)
			      {
				E_Int indR = rcvPts[noind+pt_deb];

				(densPtr+nbRcvPts*7)[noind+pt_deb] = ((vectOfRcvFields[4][indR]*vectOfGradRcvFields[0][indR]+vectOfRcvFields[0][indR]*vectOfGradRcvFields[3][indR])*cvgam)/alpha + (densPtr+nbRcvPts*7)[noind+pt_deb]*(alpha-1.)/alpha;
				(densPtr+nbRcvPts*8)[noind+pt_deb] = ((vectOfRcvFields[4][indR]*vectOfGradRcvFields[1][indR]+vectOfRcvFields[0][indR]*vectOfGradRcvFields[4][indR])*cvgam)/alpha + (densPtr+nbRcvPts*8)[noind+pt_deb]*(alpha-1.)/alpha;
				(densPtr+nbRcvPts*9)[noind+pt_deb] = ((vectOfRcvFields[4][indR]*vectOfGradRcvFields[2][indR]+vectOfRcvFields[0][indR]*vectOfGradRcvFields[5][indR])*cvgam)/alpha + (densPtr+nbRcvPts*9)[noind+pt_deb]*(alpha-1.)/alpha;
                                      
			      }
			  }

                        }  // ibc

		      //        } //chunk
		      ideb        =  ideb + ntype[ 1 + ndtyp];
		      shiftCoef   = shiftCoef + ntype[1 + ndtyp] * sizecoefs;  // shift coef   entre 2 types successif
		      shiftDonor = shiftDonor + ntype[1 + ndtyp];            // shift donor entre 2 types successif
                    }                                                            // type

                    count_rac += 1;

		  } // autorisation transfert

	      }  // irac
	  }      // pass_inst
#pragma omp barrier
      }  // pass
  }      // omp

  PyObject* infos     = PyList_New( 0 );
  PyObject* infosGrad = PyList_New( 0 );
  PyObject* allInfos  = PyList_New( 0 );

  count_rac = 0;
  for (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++)
    {
      E_Int irac_deb = 0;
      E_Int irac_fin = nrac_steady;
      if ( pass_inst == 1 ) {
	irac_deb = ipt_param_int[ech + 4 + it_target];
	irac_fin = ipt_param_int[ech + 4 + it_target + timelevel]; }

      for ( E_Int irac = irac_deb; irac < irac_fin; irac++ )
	{

	  E_Int irac_auto= irac-irac_deb;
	  if (autorisation_transferts[pass_inst][irac_auto]==1)

	    {

	      E_Int shift_rac = ech + 4 + timelevel * 2 + irac;
	      E_Int nbRcvPts  = ipt_param_int[shift_rac + nrac * 10 + 1];
	      //E_Int NoD       = ipt_param_int[shift_rac + nrac * 5     ];
	      if ( nbRcvPts > nbRcvPts_mx ) nbRcvPts_mx = nbRcvPts;


	      E_Int ibcType = ipt_param_int[shift_rac + nrac * 3];
	      E_Int ibc = 1;
	      if (ibcType < 0) ibc = 0;
	      if      ( TypeTransfert == 0 && ibc == 1 ) { continue; }
	      else if ( TypeTransfert == 1 && ibc == 0 ) { continue; }
	      // E_Int Rans=  ipt_param_int[ shift_rac + nrac*13 +1 ]; //flag raccord rans/LES
	      // E_Int nvars_loc = nvars;
	      // if(Rans==1) nvars_loc = 5;
	      PyObject* info     = PyList_New( 0 );
	      PyObject* infoGrad = PyList_New( 0 );

	      PyObject* Nozone;
	      Nozone = PyInt_FromLong( ipt_param_int[shift_rac + nrac * 11 + 1] );
	      PyList_Append( info, Nozone );  // No Zone receuveuse
	      PyList_Append( infoGrad, Nozone );  // No Zone receuveuse

	      PyList_Append( info, list_tpl[count_rac] );
	      Py_DECREF( list_tpl[count_rac] );  // tableau data

	      PyList_Append( infoGrad, list_tplGrad[count_rac] );
	      Py_DECREF( list_tplGrad[count_rac] );  // tableau data

	      E_Int     PtlistDonor = ipt_param_int[shift_rac + nrac * 12 + 1];
	      E_Int*    ipt_listRcv = ipt_param_int + PtlistDonor;
	      PyObject* listRcv     = K_NUMPY::buildNumpyArray( ipt_listRcv, nbRcvPts, 1, 1 );

	      PyList_Append( info, listRcv );
	      PyList_Append( infoGrad, listRcv );
	      Py_DECREF( listRcv );  // ListReceveur

	      PyList_Append( infos, info );
	      PyList_Append( infosGrad, infoGrad );
	      Py_DECREF( info );
	      Py_DECREF( infoGrad );
	      count_rac += 1;
	    }//autorisation transfert
	}// irac
    }// pass_inst

  PyList_Append( allInfos, infos );  
  PyList_Append( allInfos, infosGrad );    
  Py_DECREF( infos );
  Py_DECREF( infosGrad );

  delete[] ipt_param_intD; delete[] ipt_param_realD;
  delete[] list_tpl;
  delete[] frp;
  delete[] list_tplGrad;
  delete[] frpGrad;
  delete[] ipt_ndimdxD;
  delete[] ipt_roD;
  delete[] ipt_cnd;
  delete[] varStringOut;
  delete[] varStringOutGrad;
  RELEASESHAREDZ( hook, (char*)NULL, (char*)NULL );
  RELEASESHAREDN( pyParam_int, param_int );
  RELEASESHAREDN( pyParam_real, param_real );

  // return infos, infosGrad;
  return allInfos;

}
