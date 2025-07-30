/*    
      Copyright 2013-2020 Onera.

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

using namespace std;
using namespace K_FLD;

//
// LBM
//

// Getting Intersection of Q with immersed boundary
PyObject* K_CONNECTOR::___setQintersectionLBM(PyObject* self, PyObject* args){
  PyObject *zonesR;
  PyObject *zonesD;
  PyObject *pyVariables;
  PyObject *pyParam_int;
  PyObject *pyParam_real;
  E_Int nvars;
  E_Int vartype, type_transfert, no_transfert,It_target, nstep, nitmax;
  E_Int rk, exploc, num_passage;
  if (!PYPARSETUPLE_(args, OOOO_ O_ IIII_ IIII_ II_,
                    &zonesR, &zonesD, &pyVariables, &pyParam_int, &pyParam_real,&nvars, &vartype,
                    &It_target,&type_transfert, &no_transfert, &nstep, &nitmax, &rk, &exploc, &num_passage))
    {
      return NULL;
    }
  E_Int it_target= E_Int(It_target);
  E_Int varType = E_Int(vartype);
  
  //gestion nombre de pass pour ID et/ou IBC
  E_Int TypeTransfert = E_Int(type_transfert);
  E_Int pass_deb, pass_fin;
  if     (TypeTransfert==0) { pass_deb =1; pass_fin =2; }//ID
  else if(TypeTransfert==1) { pass_deb =0; pass_fin =1; }//IBCD
  else                      { pass_deb =0; pass_fin =2; }//ALL
  E_Int NoTransfert  = E_Int(no_transfert);
  
  E_Int kmd, cnNfldD, meshtype;

  E_Int nidomR   = PyList_Size(zonesR);
  E_Int nidomD   = PyList_Size(zonesD);

  //pointeur pour stocker solution au centre ET au noeud 
  E_Int** ipt_param_intR;
  
  E_Float** ipt_param_realR;
  E_Float** ipt_roR;
  E_Float** ipt_distR;

  E_Int*  ipt_ndimdxD;
  E_Int** ipt_cnd;
  
  E_Float** ipt_roD;
  E_Float** ipt_roD_vert;
  ipt_param_intR   = new E_Int*[nidomR];
  
  ipt_param_realR  = new E_Float*[nidomR*3];   //1
  ipt_roR          = ipt_param_realR + nidomR; //2
  ipt_distR        = ipt_roR         + nidomR; //3
      
  ipt_ndimdxD      = new E_Int[nidomD*8];  
  ipt_cnd          = new E_Int*[nidomD];

  ipt_roD          = new E_Float*[nidomD*8];
  ipt_roD_vert     = ipt_roD + nidomD;
  vector<PyArrayObject*> hook;
  
  /*-------------------------------------*/
  /* Extraction tableau int et real de tc*/
  /*-------------------------------------*/
  FldArrayI* param_int;E_Int res_donor = K_NUMPY::getFromNumpyArray(pyParam_int, param_int);
  E_Int* ipt_param_int = param_int->begin();
  
  FldArrayF* param_real;res_donor = K_NUMPY::getFromNumpyArray(pyParam_real, param_real);
  E_Float* ipt_param_real = param_real->begin();

  //On recupere le nom de la 1ere variable a recuperer 
  PyObject* tpl0= PyList_GetItem(pyVariables, 0); 
  char* varname = NULL;
  if (PyString_Check(tpl0)) varname = PyString_AsString(tpl0);
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(tpl0)) varname = (char*)PyUnicode_AsUTF8(tpl0); 
#endif
 
  // CB DBX
  E_Float** ipt_qD_vert; E_Float** ipt_SD_vert; E_Float** ipt_psiGD_vert;
  E_Float** ipt_qD; E_Float** ipt_SD; E_Float** ipt_psiGD;
  ipt_qD_vert    = ipt_roD_vert + nidomD;
  ipt_SD_vert    = ipt_qD_vert + nidomD;
  ipt_psiGD_vert = ipt_SD_vert + nidomD;
  ipt_qD         = ipt_psiGD_vert + nidomD;
  ipt_SD         = ipt_qD + nidomD;
  ipt_psiGD      = ipt_SD + nidomD;

  char* vartmp; char* varname1; char* varname2; char* varname3;
  E_Int nbvar_inlist = PyList_Size(pyVariables);
  for (E_Int ivar = 0; ivar < nbvar_inlist; ivar++)
    {
      PyObject* tpl0= PyList_GetItem(pyVariables, ivar);
      if (PyString_Check(tpl0)) vartmp = PyString_AsString(tpl0);
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(tpl0)) vartmp = (char*)PyUnicode_AsUTF8(tpl0);
#endif
      //printf("varname %s \n", vartmp);
      if     (ivar==0) {varname  = vartmp;}
      else if(ivar==1) {varname1 = vartmp;}  //En LBM, on a besoin d'echanger les macro !!ET!! les Q donc deux varname
      else if(ivar==2) {varname2 = vartmp;}  //Pour le couplage NS-LBM, on a besoin d'echanger les Q !!ET!! les macro !!ET!! les gradients
      else if(ivar==3) {varname3 = vartmp;}
      else {printf("Warning: souci varname setInterpTransfers \n"); }
    }
  // ENDCB

  for (E_Int nd = 0; nd < nidomD; nd++){  
    PyObject* zoneD = PyList_GetItem(zonesD, nd);
#    include "getfromzoneDcompact_all.h"
  }
  for (E_Int nd = 0; nd < nidomR; nd++){  
    PyObject* zoneR = PyList_GetItem(zonesR, nd);
#     include "getfromzoneRcompact_allLBM_Qonly.h"
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

  E_Int autorisation_transferts[pass_inst_fin][size_autorisation];

  //on dimension tableau travail pour IBC
  E_Int nbRcvPts_mx =0;
  E_Int ibcTypeMax=0;
  for  (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++){
    E_Int irac_deb= 0;
    E_Int irac_fin= nrac_steady;
    
    if (pass_inst == 1) {
      irac_deb = ipt_param_int[ ech + 4 + it_target ];
      irac_fin = ipt_param_int[ ech + 4 + it_target + timelevel ];
    }
    
    for  (E_Int irac=irac_deb; irac< irac_fin; irac++){ 
      E_Int shift_rac =  ech + 4 + timelevel*2 + irac;

      if( ipt_param_int[ shift_rac+ nrac*10 ] > nbRcvPts_mx)
	nbRcvPts_mx = ipt_param_int[ shift_rac+ nrac*10 ];
      if( ipt_param_int[shift_rac+nrac*3] > ibcTypeMax)
	ibcTypeMax =  ipt_param_int[shift_rac+nrac*3];
      
      E_Int irac_auto= irac-irac_deb;
      autorisation_transferts[pass_inst][irac_auto]=0;
      // No if statements if(rk==3 && exploc == 2){} so going directly with the else which is right below 
      autorisation_transferts[pass_inst][irac_auto]=1; 
    }
  }

  E_Int size = (nbRcvPts_mx/threadmax_sdm)+1; // on prend du gras pour gerer le residus
  E_Int r =  size % 8;
  if (r != 0) size  = size + 8 - r;           // on rajoute du bas pour alignememnt 64bits
  if (ibcTypeMax <=1 ) size = 0;              // tableau inutile : SP voir avec Ivan 

  FldArrayF  tmp(size*14*threadmax_sdm);
  E_Float* ipt_tmp = tmp.begin();

  //# pragma omp parallel default(shared)  num_threads(1)
# pragma omp parallel default(shared) num_threads(1)
  {

#ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num()+1;
    E_Int  Nbre_thread_actif = omp_get_num_threads(); // nombre de thread actif dans cette zone
#else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
#endif

    E_Int indR, type;
    E_Int indD0, indD, i, j, k, ncfLoc, indCoef, noi, sizecoefs, imd, jmd, imdjmd;

    E_Float* RcvFields;
    
    for  (E_Int ipass_typ=pass_deb; ipass_typ< pass_fin; ipass_typ++){
      for  (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++){
	E_Int irac_deb= 0;
	E_Int irac_fin= nrac_steady;
	if(pass_inst == 1){ 
	  irac_deb = ipt_param_int[ ech + 4 + it_target             ];
	  irac_fin = ipt_param_int[ ech + 4 + it_target + timelevel ];  
	}
	
	for  (E_Int irac=irac_deb; irac< irac_fin; irac++){
	  E_Int irac_auto= irac-irac_deb;
	  if (autorisation_transferts[pass_inst][irac_auto]==1){
	    E_Int shift_rac =  ech + 4 + timelevel*2 + irac;
	    E_Int NoD       =  ipt_param_int[ shift_rac + nrac*5  ];
	    E_Int loc       =  ipt_param_int[ shift_rac + nrac*9  ]; //+1 a cause du nrac mpi
	    E_Int NoR       =  ipt_param_int[ shift_rac + nrac*11 ];
	    E_Int nvars_loc =  ipt_param_int[ shift_rac + nrac*13 ]; //neq fonction raccord rans/LES
	   
	    E_Int meshtype = ipt_ndimdxD[NoD + nidomD*6];
	    E_Int cnNfldD  = ipt_ndimdxD[NoD + nidomD*7];
	    E_Int* ptrcnd  = ipt_cnd[    NoD           ];
	    E_Float* ipt_distQ ;
	    E_Int*   ipt_intrQ ;
	    E_Int*   ipt_cvel  ;
	    if (loc == 0){
	      printf("Error: transferts optimises non code en vextex " SF_D3_ "\n", shift_rac + nrac*9  +1, NoD, NoR ); 
	      imd = 0; jmd = 0;
	    }
	    else {
	      imd= ipt_param_intR[ NoD ][ NIJK ];
	      jmd= ipt_param_intR[ NoD ][ NIJK+1];
	      RcvFields = ipt_distR[NoR] ;

	      //[LBM]	      
	      ipt_distQ = &ipt_param_realR[NoR][ipt_param_intR[NoR][PT_LBM_IBC_DIST]];
	      ipt_intrQ = &ipt_param_intR [NoR][ipt_param_intR[NoR][PT_LBM_IBC_DIR]];
	      ipt_cvel  = &ipt_param_intR [NoR][ipt_param_intR[NoR][PT_LBM_Cs]];
	    }
	    
	    imdjmd = imd*jmd;

	    ////
	    //  Determine values needed for intersection of Q's
	    ////
	    
	    E_Int nbRcvPts = ipt_param_int[ shift_rac +  nrac*10 ];
  	    E_Int pos;
	    pos  = ipt_param_int[ shift_rac + nrac*7 ] ; E_Int* ntype      = ipt_param_int  + pos;
	    pos  = pos +1 + ntype[0]                   ; E_Int* types      = ipt_param_int  + pos;
	    pos  = ipt_param_int[ shift_rac + nrac*6  ]; E_Int* donorPts   = ipt_param_int  + pos;
	    pos  = ipt_param_int[ shift_rac + nrac*12 ]; E_Int* rcvPts     = ipt_param_int  + pos;   
	    pos  = ipt_param_int[ shift_rac + nrac*8  ]; E_Float* ptrCoefs = ipt_param_real + pos;	    
	    
	    E_Int nbInterpD = ipt_param_int[ shift_rac +  nrac ];
	    E_Float* xPC    =NULL;
	    xPC     = ptrCoefs + nbInterpD;


	    E_Int ideb      = 0;
	    E_Int ifin      = 0;
       
	    for (E_Int ndtyp = 0; ndtyp < ntype[0]; ndtyp++){ 
	      type      = types[ifin];

	      SIZECF(type, meshtype, sizecoefs);
	      ifin =  ifin + ntype[ 1 + ndtyp];

	      E_Int pt_deb, pt_fin;

	      // Calcul du nombre de champs a traiter par chaque thread
	      E_Int size_bc =  ifin-ideb;
	      E_Int chunk   =  size_bc/Nbre_thread_actif;
	      E_Int r       =  size_bc - chunk*Nbre_thread_actif;
	      // pts traitees par thread
	      if (ithread <= r){
		pt_deb = ideb + (ithread-1)*(chunk+1);
		pt_fin = pt_deb + (chunk+1);
	      }
	      else {
		pt_deb = ideb + (chunk+1)*r+(ithread-r-1)*chunk;
		pt_fin = pt_deb + chunk;
	      } 

	      //Si type 0, calcul sequentiel
	      if( type == 0 ){
		if (ithread ==1 ){ pt_deb = ideb; pt_fin = ifin;}
		else             { pt_deb = ideb; pt_fin = ideb;}
	      }
	      if(varType == 44 || varType == 42){
		E_Int isnside = 0;
		if (varType == 42){isnside = 1;}
		setIBCTransfersCommonVarQtagonly(rcvPts  , nbRcvPts    , pt_deb        ,
						 pt_fin  , ithread     , 
						 xPC     , xPC+nbRcvPts, xPC+nbRcvPts*2,
						 ipt_cvel, ipt_intrQ   ,
						 nvars   , imd         , jmd           ,
						 ipt_tmp , size        ,
						 ipt_param_realR[ NoR ][ LBM_zlim ],
						 isnside,
						 RcvFields);
	      }
	      else{
		printf("ipt_param_realR[ NoR ][ LBM_zlim ]=%f\n",ipt_param_realR[ NoR ][ LBM_zlim ]);
		setIBCTransfersCommonVarQ(rcvPts  , nbRcvPts    , pt_deb        ,
					  pt_fin  , ithread     , 
					  xPC     , xPC+nbRcvPts, xPC+nbRcvPts*2,
					  ipt_cvel, ipt_intrQ   , ipt_distQ     ,
					  nvars   , imd         , jmd           ,
					  ipt_tmp , size        ,
					  ipt_param_realR[ NoR ][ LBM_zlim ],
					  RcvFields);
	      }
	      ideb       = ideb + ifin;
	    }// type 
	  }// autorisation transfert
	}//irac
      }//pass_inst
#pragma omp barrier 
    }//ipass
  }// omp

  delete [] ipt_param_intR;
  delete [] ipt_param_realR;
  delete [] ipt_ndimdxD;
  delete [] ipt_cnd;
  delete [] ipt_roD;
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
  RELEASESHAREDN(pyParam_int    , param_int    );
  RELEASESHAREDN(pyParam_real   , param_real   );
  Py_INCREF(Py_None);
  return Py_None;
}

// General set interp transfer for LBM
PyObject* K_CONNECTOR::___setInterpTransfersLBM(PyObject* self, PyObject* args){
  PyObject *zonesR, *zonesD;

  PyObject *pyVariables;
  PyObject *pyParam_int, *pyParam_real;
  E_Int vartype, type_transfert, no_transfert, It_target, nstep, nitmax;
  E_Int rk, exploc, num_passage;
  E_Int global_IBM,nitrun,lbm_octree;
  E_Float gamma, cv, muS, Cs, Ts, Pr;
  
  if (!PYPARSETUPLE(args,
		    "OOOOOllllllllllll", "OOOOOiiiiiiiiiiii",
		    "OOOOOllllllllllll", "OOOOOiiiiiiiiiiii",
		    &zonesR, &zonesD, &pyVariables, &pyParam_int,  &pyParam_real, &It_target, &vartype,
		    &type_transfert, &no_transfert, &nstep, &nitmax, &rk, &exploc, &num_passage,&global_IBM,&nitrun,&lbm_octree)){
    return NULL;
  }

  E_Int it_target= E_Int(It_target);
  /* varType : 
     1  : conservatives, 
     11 : conservatives + ronutildeSA 
     2  : (ro,u,v,w,t)
     21 : (ro,u,v,w,t) + ronutildeSA 
     3  : (ro,u,v,w,p)     
     31 : (ro,u,v,w,p) + ronutildeSA 
     4  : (ro,u,v,w,t) & (Q1,..., QN)    LBM  Note: There are more options for LBM, refer to warmup in /FastLBM/Pytree.py
  */     

  E_Int varType = E_Int(vartype); 

  //gestion nombre de pass pour ID et/ou IBC
  E_Int TypeTransfert = E_Int(type_transfert);
  E_Int pass_deb, pass_fin;
  if     (TypeTransfert==0) { pass_deb =1; pass_fin =2; }//ID
  else if(TypeTransfert==1) { pass_deb =0; pass_fin =1; }//IBCD
  else                      { pass_deb =0; pass_fin =2; }//ALL

  E_Int NoTransfert  = E_Int(no_transfert);

  //printf("vartype=" SF_D_ " \n",vartype);
  E_Int kmd, cnNfldD, nvars, meshtype;
  E_Int nvars_macro_local;

  E_Int isLBM                 = 0;
  E_Int isIBM                 = 0;
  E_Int iscalcMacro_LBM       = 0;
  E_Int isneedQstar           = 0;
  
  
  //Legend for vartype for LBM
  // 4   :: default behavior for D3Q19
  // 40  :: default behavior for D3Q19  but calculate macro at ghost cell
  // 41  :: default behavior for D3Q27

  // IBM
  // 42  :: default behavior for D3Q19 Tiwari & Vanka
  // 421 :: default behavior for D3Q19 ProLB approach(-ish)  
  // 44  :: default behavior for D3Q19 penalization - Bounce Back
  // 45  :: default behavior for D3Q19 Bouzidis et al linear
  // 46  :: default behavior for D3Q19 Mei

  if ( vartype <= 3 &&  vartype >= 1){
    nvars =5;
  }
  else if( vartype == 4  ){ // [LBM]
    nvars =19;
    isLBM =1;
    isIBM =0;
    if (lbm_octree==1){
      iscalcMacro_LBM = 1;
    }
  }
  else if( vartype == 40 ){ // [LBM]
    nvars =19;
    isLBM =1;
    isIBM =0;
    iscalcMacro_LBM = 1;
  }
  else if( vartype == 41 ){ // [LBM]
    nvars =27;
    isLBM =1;
    isIBM =0;
  }
  else if( vartype == 42 || vartype == 421){ // [LBM]
    // setting iscalcMacro_LBM as interpolation of macros is needed for IBC
    // the macro calculated at ghost cells is not and should not used !!!
    nvars =19;
    isLBM =1;
    isIBM =1;    
    iscalcMacro_LBM = 1;
  }
  else if( vartype == 44 || vartype == 45 || vartype == 46 || vartype == 47){ // [LBM]
    nvars       = 19;
    isLBM       = 1;
    isIBM       = 1;
    isneedQstar = 1;
    if(vartype == 46){
      iscalcMacro_LBM = 1;
    }
  }
  else{
    nvars =6;
  }
    
  nvars_macro_local=0;
  if(iscalcMacro_LBM==1){ // [LBM]
    nvars_macro_local=5;
  }
  E_Int nidomR   = PyList_Size(zonesR);
  E_Int nidomD   = PyList_Size(zonesD);

  //pointeur pour stocker solution au centre ET au noeud 
  E_Int* ipt_ndimdxD;
  E_Int** ipt_param_intR;
  E_Int** ipt_cnd;
  
  E_Float** ipt_roR;
  E_Float** ipt_roD;
  
  E_Float** ipt_roR_vert;
  E_Float** ipt_roD_vert;
  
  E_Float** ipt_param_realR;

  // [LBM]
  E_Float** ipt_macros_localR;
  E_Float** ipt_macros_localD;

  E_Float** ipt_Qneq_localR;
  E_Float** ipt_Qneq_localD;

  E_Float** ipt_Qstar_localR;
  E_Float** ipt_Qstar_localD;

  E_Float** ipt_Qm1_localD;
  E_Float** ipt_macrosm1_localD;
  
  E_Float** ipt_cellNIBC_R;

  //ipt_ndimdxR      = new E_Int*[nidomR*3];   // on stocke ndimdx  en centre et vertexe
  ipt_param_intR   = new E_Int*[nidomR];


  ipt_roR             = new E_Float*[nidomR*7]     ;   //1
  ipt_roR_vert        = ipt_roR           + nidomR ;   //2
  ipt_param_realR     = ipt_roR_vert      + nidomR ;   //3
  ipt_cellNIBC_R      = ipt_param_realR   + nidomR ;   //4
  ipt_macros_localR   = ipt_cellNIBC_R    + nidomR ;   //5
  ipt_Qneq_localR     = ipt_macros_localR + nidomR ;   //6
  ipt_Qstar_localR    = ipt_Qneq_localR   + nidomR ;   //7
    
  ipt_ndimdxD      = new E_Int[nidomD*8];  //on stocke ndimdx, imd, jmd, en centre et vertexe, meshtype et cnDfld
  ipt_cnd          = new E_Int*[nidomD];

  ipt_roD            = new E_Float*[nidomD*7];     //1
  ipt_roD_vert       = ipt_roD           + nidomD; //2
  ipt_macros_localD  = ipt_roD_vert      + nidomD; //3
  ipt_Qneq_localD    = ipt_macros_localD + nidomD; //4
  ipt_Qstar_localD   = ipt_Qneq_localD   + nidomD; //5
  ipt_Qm1_localD     = ipt_Qstar_localD  + nidomD; //6
  ipt_macrosm1_localD= ipt_Qm1_localD    + nidomD; //7
  
  vector<PyArrayObject*> hook;
  /*-------------------------------------*/ 
  /* Extraction tableau int et real de tc*/
  /*-------------------------------------*/
  FldArrayI* param_int;
  E_Int res_donor = K_NUMPY::getFromNumpyArray(pyParam_int, param_int);
  E_Int* ipt_param_int = param_int->begin();
  FldArrayF* param_real;
  res_donor = K_NUMPY::getFromNumpyArray(pyParam_real, param_real);
  E_Float* ipt_param_real = param_real->begin();
  //On recupere le nom de la 1ere variable a recuperer 
  PyObject* tpl0= PyList_GetItem(pyVariables, 0); 
  char* varname = NULL;
  if (PyString_Check(tpl0)) varname = PyString_AsString(tpl0); 
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(tpl0)) varname = (char*)PyUnicode_AsUTF8(tpl0); 
#endif
  //on recupere sol et solcenter ainsi que connectivite et taille zones Donneuses (tc)
  for (E_Int nd = 0; nd < nidomD; nd++){  
    PyObject* zoneD = PyList_GetItem(zonesD, nd); 
#    include "getfromzoneDcompact_allLBM.h"
  }

  //on recupere sol et solcenter taille zones receuveuses, param_int et param_real (t)
  for (E_Int nd = 0; nd < nidomR; nd++){  
    PyObject* zoneR = PyList_GetItem(zonesR, nd);
#     include "getfromzoneRcompact_allLBM.h"
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

  E_Int autorisation_transferts[pass_inst_fin][size_autorisation];
  E_Int ntab_int    =18;

  //on dimension tableau travail pour IBC
  E_Int nbRcvPts_mx =0; E_Int ibcTypeMax=0;
  for  (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++){
    E_Int irac_deb= 0; E_Int irac_fin= nrac_steady;
    if (pass_inst == 1) { irac_deb = ipt_param_int[ ech + 4 + it_target ];
      irac_fin = ipt_param_int[ ech + 4 + it_target + timelevel ];}
    
    for  (E_Int irac=irac_deb; irac< irac_fin; irac++){
      
      E_Int shift_rac =  ech + 4 + timelevel*2 + irac;
      if( ipt_param_int[ shift_rac+ nrac*10] > nbRcvPts_mx) nbRcvPts_mx = ipt_param_int[ shift_rac+ nrac*10];
      if( ipt_param_int[shift_rac+nrac*3] > ibcTypeMax)  ibcTypeMax =  ipt_param_int[shift_rac+nrac*3];
      E_Int irac_auto= irac-irac_deb;
      autorisation_transferts[pass_inst][irac_auto]=0;

      // Si on est en explicit local, on va autoriser les transferts entre certaines zones seulement en fonction de la ss-ite courante
      //printf("rk,exploc=%d %d\n",rk,exploc);
      if(rk==3 && exploc == 2){
 	E_Int debut_rac = ech + 4 + timelevel*2 + nrac*ntab_int + 27*irac;     
 	E_Int levelD = ipt_param_int[debut_rac + 25];
 	E_Int levelR = ipt_param_int[debut_rac + 24];
 	E_Int cyclD  = nitmax/levelD;
	
 	// Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse   
 	if (levelD > levelR && num_passage == 1){
 	  if ( nstep%cyclD==cyclD-1 || (nstep%cyclD==cyclD/2 && (nstep/cyclD)%2==1) ){
 	    autorisation_transferts[pass_inst][irac_auto]=1;
 	  }
 	  else {continue;}
 	}
 	// Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
 	else if (levelD < levelR && num_passage == 1){
 	  if (nstep%cyclD==1 || nstep%cyclD==cyclD/4 || nstep%cyclD== cyclD/2-1 || nstep%cyclD== cyclD/2+1 || nstep%cyclD== cyclD/2+cyclD/4 || nstep%cyclD== cyclD-1)
 	    { autorisation_transferts[pass_inst][irac_auto]=1; }
 	  else {continue;}
 	}
 	// Le pas de temps de la zone donneuse est egal a celui de la zone receveuse
 	else if (levelD == levelR && num_passage == 1){
 	  if (nstep%cyclD==cyclD/2-1 || (nstep%cyclD==cyclD/2 && (nstep/cyclD)%2==0) || nstep%cyclD==cyclD-1) 
 	    { autorisation_transferts[pass_inst][irac_auto]=1; }
 	  else {continue;}
 	}
 	// Le pas de temps de la zone donneuse est egal a celui de la zone receveuse (cas du deuxieme passage)   
 	else if (levelD == ipt_param_int[debut_rac +24] && num_passage == 2){
 	  if (nstep%cyclD==cyclD/2 && (nstep/cyclD)%2==1)
 	    { autorisation_transferts[pass_inst][irac_auto]=1; }
 	  else {continue;}
 	}
 	else {continue;} 
      }
      else if(exploc == 19){
	E_Int debut_rac = ech + 4 + timelevel*2 + nrac*ntab_int + 27*irac;     
	E_Int levelD = ipt_param_int[debut_rac + 25];
	E_Int levelR = ipt_param_int[debut_rac + 24];
	E_Int cyclD  = nitmax/levelD;
	// Le pas de temps de la zone donneuse est plus petit que celui de la zone receveuse   
	if (lbm_octree == 1){
	    if (levelD > levelR && num_passage == 1){
	      if ( nstep==2 ){
		autorisation_transferts[pass_inst][irac_auto]=1;
	      }
	      else {continue;}
	    }
	    // Le pas de temps de la zone donneuse est plus grand que celui de la zone receveuse
	    else if (levelD < levelR && num_passage == 1){
	      autorisation_transferts[pass_inst][irac_auto]=1;
	    }
	}
 	// Le pas de temps de la zone donneuse est egal a celui de la zone receveuse
 	else if (levelD == levelR && num_passage == 1){
 	  autorisation_transferts[pass_inst][irac_auto]=1; 
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

  FldArrayF  tmp(size*14*threadmax_sdm);
  E_Float* ipt_tmp = tmp.begin();

  E_Float** RcvFields = new E_Float*[ nvars*threadmax_sdm];
  E_Float** DnrFields = new E_Float*[ nvars*threadmax_sdm];

  // [LBM]
  E_Float meax       = 0;
  E_Float meay       = 0;
  E_Float meaz       = 0;
  
  E_Float param_realLoc[30]; 
  param_realLoc[ GAMMA] = gamma;
  param_realLoc[ CVINF] = cv;
  param_realLoc[ XMUL0] = muS;
  param_realLoc[ CS] = Cs;
  param_realLoc[ TEMP0] = Ts;
  param_realLoc[ PRANDT] = Pr;

# pragma omp parallel default(shared) num_threads(1) reduction(+:meax,meay,meaz)
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

    // [LBM ==> Q or Qstar]
    vector<E_Float*> vectOfRcvFieldsLoc(nvars);
    E_Float** vectOfRcvFields = RcvFields + nvars*(ithread-1);
    E_Float** vectOfDnrFields = DnrFields + nvars*(ithread-1);

    // [LBM ==> Macro properties]
    vector<E_Float*> vectOfmacroRcvFields(nvars_macro_local);
    vector<E_Float*> vectOfmacroDnrFields(nvars_macro_local);

    // [LBM ==> Qneq]
    vector<E_Float*> vectOfQneqRcvFields(nvars);
    vector<E_Float*> vectOfQneqDnrFields(nvars);

    // [LBM ==> Qstar]
    vector<E_Float*> vectOfQstarRcvFields(nvars);
    vector<E_Float*> vectOfQstarDnrFields(nvars);

    // [LBM ==> QM1 & Macro M1]
    vector<E_Float*> vectOfQm1DnrFields(nvars);
    vector<E_Float*> vectOfmacrom1DnrFields(nvars);

    // [LBM ==> cellNIBC]
    vector<E_Float*> vectOfcellNIBMRcvFields(3);

    
    //1ere pass_typ: IBC
    //2eme pass_typ: transfert

    for  (E_Int ipass_typ=pass_deb; ipass_typ< pass_fin; ipass_typ++){
      //1ere pass_inst: les raccords fixes
      //2eme pass_inst: les raccords instationnaires
      for  (E_Int pass_inst=pass_inst_deb; pass_inst< pass_inst_fin; pass_inst++){
 	E_Int irac_deb= 0; E_Int irac_fin= nrac_steady;
 	if(pass_inst == 1){ 
 	  irac_deb = ipt_param_int[ ech + 4 + it_target             ];
 	  irac_fin = ipt_param_int[ ech + 4 + it_target + timelevel ];  
 	}

 	for  (E_Int irac=irac_deb; irac< irac_fin; irac++){
 	  E_Int irac_auto= irac-irac_deb;
	  E_Int debut_rac = ech + 4 + timelevel*2 + nrac*ntab_int + 27*irac; 
	  E_Int levelD = ipt_param_int[debut_rac + 25];
	  E_Int levelR = ipt_param_int[debut_rac + 24];
	  
 	  if (autorisation_transferts[pass_inst][irac_auto]==1){
	    //printf("AUTHORIZED :: levelR,levelD=%d %d\n", levelR,levelD);
 	    E_Int shift_rac    =  ech + 4 + timelevel*2 + irac;

 	    // ibcType values
 	    //-1 : raccord ID
 	    // 0 : wallslip
 	    // 1 : noslip
 	    // 2 : log
 	    // 3 : Musker
 	    // 4 : outpress
 	    E_Int ibcType =  ipt_param_int[shift_rac+nrac*3];
	    
 	    E_Int ibc = 1;
 	    if (ibcType < 0) ibc = 0;
 	    if(1-ibc != ipass_typ)  continue;
 	    E_Int NoD      =  ipt_param_int[ shift_rac + nrac*5  ];
 	    E_Int loc      =  ipt_param_int[ shift_rac + nrac*9  ];
 	    E_Int NoR      =  ipt_param_int[ shift_rac + nrac*11 ];
 	    E_Int nvars_loc=  ipt_param_int[ shift_rac + nrac*13 ]; //neq fonction raccord rans/LES
 	    E_Int rotation =  ipt_param_int[ shift_rac + nrac*14 ]; //flag pour periodicite azimutale

	    

	   
 	    E_Int meshtype = ipt_ndimdxD[NoD + nidomD*6];
 	    E_Int cnNfldD  = ipt_ndimdxD[NoD + nidomD*7];
 	    E_Int* ptrcnd  = ipt_cnd[    NoD           ];

	    // [LBM]
	    E_Int   num_bcs;
	    E_Int*  ipt_cvel ;
 	    E_Int*  ipt_intrQ ;
 	    E_Int*  ipt_cminus ;
 	    E_Int*  ipt_bcs ;
	    
 	    E_Float c0_local          = ipt_param_realR[NoR][LBM_c0];
	    E_Float c02               = 1./(c0_local*c0_local);
	    E_Float c04               = c02*c02;
	    E_Float gamma_precon_inv  = 1./ipt_param_realR[NoR][LBM_gamma_precon];
	    E_Float tau_relax         = ipt_param_realR[NoR][LBM_difcoef]*gamma_precon_inv*c02+0.5;
	    E_Float tau_relaxD        = ipt_param_realR[NoD][LBM_difcoef]*gamma_precon_inv*c02+0.5;
	    E_Float tau_relaxR        = ipt_param_realR[NoR][LBM_difcoef]*gamma_precon_inv*c02+0.5;
	    E_Float* ipt_wdist ;
	    E_Float* ipt_distQ ;	    
	    E_Float* ipt_H2H3;

	    // IMPORTANT NOTE:: calc of tau_relax assumes constant nu & not in sponge layers
	    
	    
 	    if (loc == 0){
 	      printf("Error: transferts optimises non code en vextex " SF_D3_ "\n", shift_rac + nrac*9  +1, NoD, NoR ); 
 	      //imd= ipt_ndimdxD[ NoD+ nidomD*4]; jmd= ipt_ndimdxD[ NoD + nidomD*5];
 	      imd = 0; jmd = 0;
 	    }
 	    else {
 	      imd= ipt_param_intR[ NoD ][ NIJK ]; jmd= ipt_param_intR[ NoD ][ NIJK+1];
 	      for (E_Int eq = 0; eq < nvars_loc; eq++){
		vectOfRcvFieldsLoc[eq] = ipt_roR[ NoR] + eq*ipt_param_intR[ NoR ][ NDIMDX ];
		vectOfRcvFields[eq] = ipt_roR[ NoR] + eq*ipt_param_intR[ NoR ][ NDIMDX ];
        vectOfDnrFields[eq] = ipt_roD[ NoD] + eq*ipt_param_intR[ NoD ][ NDIMDX ];
 	      }
	      
	      // [LBM]
 	      if(iscalcMacro_LBM==1){
 		for (E_Int eq = 0; eq < nvars_macro_local; eq++){
 		  vectOfmacroRcvFields[eq] = ipt_macros_localR[NoR] + eq*ipt_param_intR[ NoR ][ NDIMDX ];
 		  vectOfmacroDnrFields[eq] = ipt_macros_localD[NoD] + eq*ipt_param_intR[ NoD ][ NDIMDX ];
 		}
 		ipt_cvel  = &ipt_param_intR [NoR][ipt_param_intR[NoR][PT_LBM_Cs]];
 	      }
	      if(lbm_octree==1){
		for (E_Int eq = 0; eq < nvars_loc; eq++){
 		  vectOfQneqRcvFields[eq] = ipt_Qneq_localR[ NoR] + eq*ipt_param_intR[ NoR ][ NDIMDX ];
 		  vectOfQneqDnrFields[eq] = ipt_Qneq_localD[ NoD] + eq*ipt_param_intR[ NoD ][ NDIMDX ];
		  vectOfQm1DnrFields[eq]  = ipt_Qm1_localD[ NoD]  + eq*ipt_param_intR[ NoD ][ NDIMDX ];
		}
		for (E_Int eq = 0; eq < nvars_macro_local; eq++){
		  vectOfmacrom1DnrFields[eq]= ipt_macrosm1_localD[NoD] + eq*ipt_param_intR[ NoD ][ NDIMDX ];
 		}
		ipt_wdist   = &ipt_param_realR[NoR][ipt_param_intR[NoR][PT_LBM_Ws]];	
	      }		
 	      if(isIBM==1){
 		for (E_Int eq = 0; eq < nvars_loc; eq++){
 		  vectOfQneqRcvFields[eq] = ipt_Qneq_localR[ NoR] + eq*ipt_param_intR[ NoR ][ NDIMDX ];
 		  vectOfQneqDnrFields[eq] = ipt_Qneq_localD[ NoD] + eq*ipt_param_intR[ NoD ][ NDIMDX ];
 		}
		// [for HRR]
		for (E_Int eq = 0; eq < 3; eq++){
		  vectOfcellNIBMRcvFields[eq] = ipt_cellNIBC_R[ NoR] + eq*ipt_param_intR[ NoR ][ NDIMDX ];
		}
	        num_bcs     =  ipt_param_intR [NoR][LBM_NQ_BC];
		ipt_intrQ   = &ipt_param_intR [NoR][ipt_param_intR[NoR][PT_LBM_IBC_DIR]];
 		ipt_cvel    = &ipt_param_intR [NoR][ipt_param_intR[NoR][PT_LBM_Cs]];
 		ipt_cminus  = &ipt_param_intR [NoR][ipt_param_intR[NoR][PT_LBM_Cminus]]; 		
 		ipt_bcs     = &ipt_param_intR [NoR][ipt_param_intR[NoR][PT_LBM_BC]];
		
 		ipt_wdist   = &ipt_param_realR[NoR][ipt_param_intR[NoR][PT_LBM_Ws]];		
 		ipt_distQ   = &ipt_param_realR[NoR][ipt_param_intR[NoR][PT_LBM_IBC_DIST]]; 				
		ipt_H2H3    = &ipt_param_realR [NoR][ipt_param_intR[NoR][PT_LBM_H2H3]];
		
		if(isneedQstar==1){
		  for (E_Int eq = 0; eq < nvars_loc; eq++){
		    vectOfQstarRcvFields[eq] = ipt_Qstar_localR[ NoR] + eq*ipt_param_intR[ NoR ][ NDIMDX ];
		    vectOfQstarDnrFields[eq] = ipt_Qstar_localD[ NoD] + eq*ipt_param_intR[ NoD ][ NDIMDX ];
		  }
		}		
	      }
 	    }
	    
 	    imdjmd = imd*jmd;

 	    ////
 	    //  Interpolation parallele
 	    ////
	    
 	    E_Int nbRcvPts = ipt_param_int[ shift_rac +  nrac*10 ];

 	    E_Int pos;
 	    pos  = ipt_param_int[ shift_rac + nrac*7 ]; E_Int* ntype      = ipt_param_int  + pos;
 	    pos  = pos +1 + ntype[0]                  ; E_Int* types      = ipt_param_int  + pos;
 	    pos  = ipt_param_int[ shift_rac + nrac*6 ]; E_Int* donorPts   = ipt_param_int  + pos;
 	    pos  = ipt_param_int[ shift_rac + nrac*12]; E_Int* rcvPts     = ipt_param_int  + pos;   // donor et receveur inverser car storage donor
 	    pos  = ipt_param_int[ shift_rac + nrac*8 ]; E_Float* ptrCoefs = ipt_param_real + pos;

    
	    E_Int nbInterpD = ipt_param_int[ shift_rac +  nrac ];
 	    E_Float* xPC=NULL; E_Float* xPI=NULL; E_Float* xPW=NULL; E_Float* densPtr=NULL;
 	    if (ibc == 1){ 
 	      xPC     = ptrCoefs + nbInterpD;
 	      xPI     = ptrCoefs + nbInterpD +3*nbRcvPts;
 	      xPW     = ptrCoefs + nbInterpD +6*nbRcvPts;
 	      // [LBM] pointer below is also used for Q's for LBM
 	      densPtr = ptrCoefs + nbInterpD +9*nbRcvPts;
	    }

 	    E_Int ideb        = 0;
 	    E_Int ifin        = 0;
 	    E_Int shiftCoef   = 0;
 	    E_Int shiftDonor  = 0;
        E_Int shiftv = 0;
     
 	    for (E_Int ndtyp = 0; ndtyp < ntype[0]; ndtyp++){ 
 	      type      = types[ifin];

 	      SIZECF(type, meshtype, sizecoefs);
 	      ifin =  ifin + ntype[ 1 + ndtyp];

 	      E_Int pt_deb, pt_fin;

 	      /// oldschool
 	      // Calcul du nombre de champs a traiter par chaque thread
 	      E_Int size_bc =  ifin-ideb;
 	      E_Int chunk   =  size_bc/Nbre_thread_actif;
 	      E_Int r       =  size_bc - chunk*Nbre_thread_actif;
 	      // pts traitees par thread
 	      if (ithread <= r){
 		pt_deb = ideb + (ithread-1)*(chunk+1);
 		pt_fin = pt_deb + (chunk+1);
 	      }
 	      else {
 		pt_deb = ideb + (chunk+1)*r+(ithread-r-1)*chunk;
 		pt_fin = pt_deb + chunk;
 	      } 

 	      //Si type 0, calcul sequentiel
 	      if( type == 0 ){
 		if (ithread ==1 ){ pt_deb = ideb; pt_fin = ifin;}
 		else             { pt_deb = ideb; pt_fin = ideb;}
 	      }
	      
 	      noi       = shiftDonor;                             // compteur sur le tableau d indices donneur
 	      indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
	      //[LBM]
	      // IMPORTANT NOTE:: for BB & IBB I cannot calculate the image point value as it supresses the real values
	      // at the location of the ghost ndoes
	      if (ibc == 1 && (global_IBM<4 || global_IBM==21 || global_IBM==31 || global_IBM==41)){
		goto skipping_interpolation;
	      }
 	      if (nvars_loc==5){
		printf("nvars_loc==5\n");
#           include "commonInterpTransfers_reorder_5eq.h" 
 	      }
 	      else if(nvars_loc==6){
		printf("nvars_loc==6\n");
#           include "commonInterpTransfers_reorder_6eq.h" 
 	      }
 	      else if(nvars_loc==19)  {
		//printf("nvars_loc=19\n");
		// loop unrolling for efficiency
#           include "commonInterpTransfers_reorder_19eqLBM.h"
		// [AJ] OCTREE
//		if (lbm_octree==1){
//		  printf("interpolation LBM setInterTransfer \n");
////		  noi       = shiftDonor;     
////		  indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
////#           include "commonInterpTransfers_reorder_19eqLBM_neq_octree.h"
// 		  noi       = shiftDonor;        
//		  indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
//#           include "commonInterpTransfers_reorder_5eqLBM_macro.h" // seems correct
//		  noi       = shiftDonor;        
//		  indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
//#           include "commonInterpTransfers_reorder_19eqLBM_octree.h"
//		}
		if(vartype==40){
		  noi       = shiftDonor;        
		  indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
#           include "commonInterpTransfers_reorder_5eqLBM.h"
		}
		if (vartype==42 or vartype==421){
		  noi       = shiftDonor;     
		  indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
#           include "commonInterpTransfers_reorder_19eqLBM_neq.h"		  
		  noi       = shiftDonor;        
		  indCoef   = (pt_deb-ideb)*sizecoefs +  shiftCoef;
#           include "commonInterpTransfers_reorder_5eqLBM.h"
		}
	      }
 	      else {
#           include "commonInterpTransfers_reorder_neqLBM.h"		
 	      }

    
 	      // Prise en compte de la periodicite par rotation
 	      if (rotation == 1){
 		E_Float* angle = ptrCoefs + nbInterpD;
#          include "includeTransfers_rotation.h"
 	      }

	      skipping_interpolation:
 	      // IBC    
 	      if (ibc == 1){
		
 	      	Pr    = ipt_param_realR[ NoR ][ PRANDT ];
 	      	Ts    = ipt_param_realR[ NoR ][ TEMP0 ];
 	      	Cs    = ipt_param_realR[ NoR ][ CS ];
 	      	muS   = ipt_param_realR[ NoR ][ XMUL0 ];
 	      	cv    = ipt_param_realR[ NoR ][ CVINF ];
 	      	gamma = ipt_param_realR[ NoR ][ GAMMA ];
 		//1) rho
 		//2) pressure
 		//3) vel u
 		//4) vel v
 		//5) vel w
 		//6) utau
 		//7) yplus
 		//8) stagnation enthalpy
 		//9) stagnation pressure
 		//10) dirx
 		//11) diry
 		//12) dirz
 		//13) q
 		//14) qstar
 		//15) qneq
		  
 	      	if (varType == 2 || varType == 21){
 	      	  setIBCTransfersCommonVar2(ibcType, rcvPts, nbRcvPts, pt_deb, pt_fin, ithread,
 	      				    xPC    , xPC     +nbRcvPts, xPC     +nbRcvPts*2,
 	      				    xPW    , xPW     +nbRcvPts, xPW     +nbRcvPts*2,
 	      				    xPI    , xPI     +nbRcvPts, xPI     +nbRcvPts*2, 
 	      				    densPtr, 
 	      				    ipt_tmp, size, nvars,
                                            param_realLoc,
 	      				    vectOfDnrFields, vectOfRcvFields);
		   }
 	      	else if (varType == 42){
		  //printf("==== IBM Tiwari & Vanka =%d\n",varType);
 	      	  ////order of added pts
 	      	  ////1)Q      densPtr+nbRcvPts*12
 	      	  ////2)Qstar  densPtr+nbRcvPts*12+  nvars*nbRcvPts  
 	      	  ////3)Qneq   densPtr+nbRcvPts*12+2*nvars*nbRcvPts 
	      	  //
 	      	  // Interpolation for Macro
 		  setIBCTransfersCommonVar2LBM(ibcType           , rcvPts             , nbRcvPts          ,
					       pt_deb            , pt_fin             , ithread           ,
					       xPC               , xPC     +nbRcvPts  , xPC     +nbRcvPts*2,
					       xPW               , xPW     +nbRcvPts  , xPW     +nbRcvPts*2,
					       xPI               , xPI     +nbRcvPts  , xPI     +nbRcvPts*2, 
					       densPtr           , densPtr+nbRcvPts   , 
					       densPtr+nbRcvPts*2, densPtr+nbRcvPts*3 , densPtr+nbRcvPts*4 , 
					       densPtr+nbRcvPts*5, densPtr+nbRcvPts*6 , 
					       densPtr+nbRcvPts*7, densPtr+nbRcvPts*8 ,
					       densPtr+nbRcvPts*9, densPtr+nbRcvPts*10, densPtr+nbRcvPts*11,
					       ipt_tmp           , size               , gamma              ,
					       cv                , muS                , Cs                 ,
					       Ts                , Pr                 , ipt_intrQ          ,
					       nvars_loc         ,
					       vectOfmacroDnrFields, vectOfmacroRcvFields);

 		  E_Int shift_qstart = nbRcvPts*12+nvars*nbRcvPts;
 		  setIBCTransfersCommonVar5(ibcType             , rcvPts           , nbRcvPts          ,
					    pt_deb              , pt_fin           , ithread           ,
					    ipt_cvel            , ipt_wdist        , ipt_param_realR[ NoR ][ LBM_c0 ],
					    ipt_param_realR[ NoR ][ LBM_difcoef ],
					    densPtr+shift_qstart, varType          , ipt_intrQ         , gamma_precon_inv,
					    ipt_H2H3            , imd              , imdjmd            ,
					    vectOfmacroRcvFields ,vectOfQneqRcvFields  ,vectOfRcvFieldsLoc,
					    vectOfcellNIBMRcvFields);

 		}
		else if (varType == 421){
		  E_Int ibcType_local = 11;
		  //printf("==== IBM Feng et al. 2019 \n");
		  // Interpolation for Macro
		  setIBCTransfersCommonVar2LBM(ibcType_local     , rcvPts             , nbRcvPts          ,
					       pt_deb            , pt_fin             , ithread           ,
					       xPC               , xPC     +nbRcvPts  , xPC     +nbRcvPts*2,
					       xPW               , xPW     +nbRcvPts  , xPW     +nbRcvPts*2,
					       xPI               , xPI     +nbRcvPts  , xPI     +nbRcvPts*2, 
					       densPtr           , densPtr+nbRcvPts   , 
					       densPtr+nbRcvPts*2, densPtr+nbRcvPts*3 , densPtr+nbRcvPts*4 , 
					       densPtr+nbRcvPts*5, densPtr+nbRcvPts*6 , 
					       densPtr+nbRcvPts*7, densPtr+nbRcvPts*8 ,
					       densPtr+nbRcvPts*9, densPtr+nbRcvPts*10, densPtr+nbRcvPts*11,
					       ipt_tmp           , size               , gamma              ,
					       cv                , muS                , Cs                 ,
					       Ts                , Pr                 , ipt_intrQ          ,
					       nvars_loc         ,
					       vectOfmacroDnrFields, vectOfmacroRcvFields);

		  ibcType_local = 1; // [AJ] 1 = extrapolation of Qneq & 11 = HRR (reconstruction)
		  E_Int shift_qstart = nbRcvPts*12+nvars*nbRcvPts;
 		  setIBCTransfersCommonVar5(ibcType_local       , rcvPts           , nbRcvPts          ,
					    pt_deb              , pt_fin           , ithread           ,
					    ipt_cvel            , ipt_wdist        , ipt_param_realR[ NoR ][ LBM_c0 ],
					    ipt_param_realR[ NoR ][ LBM_difcoef ],
					    densPtr+shift_qstart, varType          , ipt_intrQ         , gamma_precon_inv,
					    ipt_H2H3            , imd              , imdjmd            ,
					    vectOfmacroRcvFields ,vectOfQneqRcvFields  ,vectOfRcvFieldsLoc,
					    vectOfcellNIBMRcvFields );
		  
		}
 		else if (varType == 44){
		  //printf("==== IBM Penalization");
		  setIBCTransfersCommonVar44(rcvPts  , nbRcvPts     , pt_deb,
					     pt_fin  , ithread      ,
 					     xPC     , xPC +nbRcvPts, xPC +nbRcvPts*2    ,
 					     ipt_cvel, ipt_intrQ    , ipt_cminus         ,
					     ipt_bcs , num_bcs      , densPtr+nbRcvPts*12,
 					     ipt_tmp , size         , meax,
					     meay, meaz, 
 					     vectOfRcvFieldsLoc,vectOfQstarRcvFields);
 		}
		else if (varType == 45){
		  //printf("==== IBM Bouzidis et al. \n");
		  setIBCTransfersCommonVar45(rcvPts     , nbRcvPts     , pt_deb            ,
					     pt_fin     , ithread      ,
 					     xPC        , xPC +nbRcvPts, xPC +nbRcvPts*2   ,
 					     ipt_cvel   , ipt_intrQ    , ipt_distQ         ,
					     ipt_cminus , imd          , jmd               ,
					     ipt_bcs    , num_bcs      , densPtr+nbRcvPts*12,
 					     ipt_tmp    , size         , meax,
					     meay, meaz, 
 					     vectOfRcvFieldsLoc,vectOfQstarRcvFields);
 		}
		else if (varType == 46){
		  //printf("==== IBM Mei \n");
		  setIBCTransfersCommonVar46(rcvPts           , nbRcvPts     , pt_deb             ,
					     pt_fin           , ithread      ,			  
 					     xPC              , xPC +nbRcvPts, xPC +nbRcvPts*2    ,
 					     ipt_cvel         , ipt_intrQ    , ipt_distQ          ,
					     ipt_cminus       , imd          , jmd                ,
					     tau_relax        , ipt_wdist    , c0_local           ,
					     gamma_precon_inv , imd          , imdjmd             ,
					     ipt_bcs          , num_bcs      , densPtr+nbRcvPts*12,
 					     ipt_tmp          , size         , meax,
					     meay             , meaz         , 
 					     vectOfRcvFieldsLoc,vectOfQstarRcvFields,vectOfmacroRcvFields);
 		}
		else if (varType == 47){
		  //printf("==== IBM Zhao \n");
		  setIBCTransfersCommonVar47(rcvPts     , nbRcvPts     , pt_deb            ,
					     pt_fin     , ithread      ,
 					     xPC        , xPC +nbRcvPts, xPC +nbRcvPts*2   ,
 					     ipt_cvel   , ipt_intrQ    , ipt_distQ         ,
					     ipt_cminus , imd          , jmd               ,
					     tau_relax  , ipt_wdist    , c0_local          ,
					     ipt_bcs    , num_bcs      , densPtr+nbRcvPts*12,
 					     ipt_tmp    , size         , meax,
					     meay       , meaz         , 
 					     vectOfRcvFieldsLoc,vectOfQstarRcvFields);
 		}
		
 	      }
 	      ideb       = ideb + ifin;
 	      shiftCoef  = shiftCoef  +  ntype[1+ndtyp]*sizecoefs; //shift coef   entre 2 types successif
 	      shiftDonor = shiftDonor +  ntype[1+ndtyp];           //shift donor entre 2 types successif
 	    }// type 
 	  }// autorisation transfert
 	}//irac
      }//pass_inst
#pragma omp barrier 
    }//ipass
  }// omp

  /*
  E_Float maxval_mea= max({abs(meax),abs(meay),abs(meaz)});
  if(maxval_mea>1e-9  && TypeTransfert==2) 
  {
    FILE* outfile;
    outfile=fopen("mea_data_setinterp.txt", "a");
    fprintf(outfile,"%d %e %e %e\n",nitrun,meax,meay,meaz);
    fclose(outfile);
  }*/
  
  delete [] RcvFields; delete [] DnrFields;
  delete [] ipt_param_intR;
  delete [] ipt_roR;
  delete [] ipt_ndimdxD;
  delete [] ipt_roD;
  delete [] ipt_cnd;
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
  RELEASESHAREDN(pyParam_int    , param_int    );
  RELEASESHAREDN(pyParam_real   , param_real   );
  Py_INCREF(Py_None);

  return Py_None;
}
