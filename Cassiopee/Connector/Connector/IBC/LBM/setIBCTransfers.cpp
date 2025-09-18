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
# include <math.h>
using namespace std;
using namespace K_FLD;

//# include "IBC/commonLaws.h"
//# define NUTILDE_FERRARI 2

// on suppose les variables compactees densPtr, pressPtr, 
E_Int K_CONNECTOR::setIBCTransfersCommonVar2LBM(E_Int bctype    , E_Int* rcvPts    , E_Int& nbRcvPts,
						E_Int& ideb     , E_Int& ifin      , E_Int& ithread ,
						E_Float* xPC    , E_Float* yPC     , E_Float* zPC   ,
						E_Float* xPW    , E_Float* yPW     , E_Float* zPW   ,
						E_Float* xPI    , E_Float* yPI     , E_Float* zPI   , 
						E_Float* densPtr, E_Float* pressPtr, 		    
						E_Float* vxPtr  , E_Float* vyPtr   , E_Float* vzPtr , 
						E_Float* utauPtr, E_Float* yplusPtr,
						E_Float* d1     , E_Float* d2      ,
						E_Float* d3     , E_Float* d4      , E_Float* d5    ,
						E_Float* tmp    , E_Int& size      , E_Float gamma  , 
						E_Float cv      , E_Float muS      , E_Float Cs     ,
						E_Float Ts      , E_Float Pr       , E_Int* Qdir    ,
						E_Int& nvarsQ   ,
						vector<E_Float*>& vectOfDnrFields, vector<E_Float*>& vectOfRcvFields,
						E_Int nbptslinelets, E_Float* linelets, E_Int* indexlinelets ){
  
  /* lois de paroi */
  E_Float roext, uext, pext, text, muext, yext, yplus, yibc;
  E_Float uscaln, un, vn, wn, ut, vt, wt, utau, utauv, utau0, umod;
  E_Float aa, bb, dd, fp, tp, f1v;
  E_Float expy, denoml10,ax,l1,l2, l3;
  E_Float ucible0, ucible, vcible, wcible, nutilde, signibc, twall, rowall, muwall;
  E_Float rhocible;
  E_Int npass;
  
  // Lois de paroi: criteres d'arret pour estimer le frottement par Newton
  E_Float newtoneps = 1.e-7; // critere d'arret pour u+
  E_Float newtonepsnutilde = 1.e-10; // critere d arret pour nutilde
  E_Float newtonepsprime = 1.e-12;// critere d'arret pour la derivee  
  E_Float cvgam = cv*(gamma-1.);
  E_Float cvgaminv = 1./(cvgam);
  E_Float coefSuth = muS * (1.+Cs/Ts);
  E_Float Tsinv = 1./Ts;
  E_Float kappa = 0.4; // Constante de Von Karman
  E_Float kappainv = 1./kappa;
  E_Float cc = 5.2; //pour la loi log
  E_Float one_third = 1./3.; 
  /* fin parametres loi de parois*/

  E_Int nvars  = vectOfDnrFields.size();
  
  E_Float a0,a1,a2,b0,b1,b2,n0,n1,n2;
  E_Float gamma_local, eta_local;
  E_Float normb, rho, u, v, w, rho_inv;
  E_Float vnc, alpha, beta, alphasbeta;
  E_Float* roOut = vectOfRcvFields[0];// ro
  E_Float* uOut  = vectOfRcvFields[1];// u
  E_Float* vOut  = vectOfRcvFields[2];// v
  E_Float* wOut  = vectOfRcvFields[3];// w
  E_Float* tOut  = vectOfRcvFields[4];// temperature
  E_Float* varSAOut = NULL;

  // ODE-based wall model

  E_Float Cv1cube = pow(7.1,3);
  E_Int nmax = 20;

  E_Float L2norm ; 
  E_Float L2norm0;  
  
  E_Float ynm,ynp,dy,dym,dyp,nutm,nutp;
  E_Float xim,xi,xip,m;
  E_Float nutrm,nutrp,nutr;

  E_Float* ipt_u1dold = NULL;
  E_Float* yline        = NULL;
  E_Float* u_line       = NULL;
  E_Float* nutilde_line = NULL;
  E_Float* matm_line    = NULL;
  E_Float* mat_line     = NULL;
  E_Float* matp_line    = NULL;
  E_Float* alphasbeta_line = NULL;
  
  FldArrayF u1dold(nbptslinelets);
  ipt_u1dold = u1dold.begin();

  if (nbptslinelets != 0){
    yline           = linelets;
    u_line          = linelets + nbRcvPts*nbptslinelets;
    nutilde_line    = linelets + nbRcvPts*nbptslinelets*2;
    matm_line       = linelets + nbRcvPts*nbptslinelets*3;
    mat_line        = linelets + nbRcvPts*nbptslinelets*4;
    matp_line       = linelets + nbRcvPts*nbptslinelets*5;
    alphasbeta_line = linelets + nbRcvPts*nbptslinelets*6;    
  }

  if (nvars == 6) varSAOut = vectOfRcvFields[5]; // nutildeSA
  
  if (bctype == 1){
#ifdef _OPENMP4
#pragma omp simd
#endif 
    for (E_Int noind = 0; noind < ifin-ideb; noind++){
      E_Int indR = rcvPts[noind+ideb];	
      if(Qdir[0+noind*nvarsQ] == 1){
	rho_inv = 1./roOut[indR]; 
	// vitesse
	u = uOut[indR]*rho_inv; 
	v = vOut[indR]*rho_inv;
	w = wOut[indR]*rho_inv;
# include "IBC/commonBCType1.h"
	uOut[indR] = ucible*roOut[indR]; 
	vOut[indR] = vcible*roOut[indR]; 
	wOut[indR] = wcible*roOut[indR];      
	if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

	pressPtr[noind + ideb] = roOut[indR]* tOut[indR]*cvgam; 
	densPtr[noind + ideb]  = roOut[indR];	
	vxPtr  [noind + ideb]  = uOut[indR];
	vyPtr  [noind + ideb]  = vOut[indR];
	vzPtr  [noind + ideb]  = wOut[indR];
      }
    }
  }
  else  if (bctype == 11){
    // Quadratic interpolation for density
    // Better than linear (bctype==1)
#ifdef _OPENMP4
#pragma omp simd
#endif 
    for (E_Int noind = 0; noind < ifin-ideb; noind++){
      E_Int indR = rcvPts[noind+ideb];	
      if(Qdir[0+noind*nvarsQ] == 1){
	rho_inv = 1./roOut[indR]; 
	// vitesse
	u   = uOut[indR]*rho_inv; 
	v   = vOut[indR]*rho_inv;
	w   = wOut[indR]*rho_inv;
	rho = roOut[indR];
# include "IBC/LBM/commonBCType11Dens.h"
	uOut[indR]  = ucible*rhocible; 
	vOut[indR]  = vcible*rhocible; 
	wOut[indR]  = wcible*rhocible;
	roOut[indR] = rhocible;      
	if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

	pressPtr[noind + ideb] = roOut[indR]* tOut[indR]*cvgam;
	densPtr[noind + ideb]  = roOut[indR];
	vxPtr  [noind + ideb]  = uOut[indR];
	vyPtr  [noind + ideb]  = vOut[indR];
	vzPtr  [noind + ideb]  = wOut[indR];
      }
    }
  }
  else{
    printf("Warning !!! setIBCTransfersCommonVar2LBM: bcType " SF_D_ " not implemented.\n", bctype);
    return 0;
  }
  return 1;
}

//==================================================================================
// Regular LBM - determine Q at ghost cell
//==================================================================================
E_Int K_CONNECTOR::setIBCTransfersCommonVar5(E_Int bctype    , E_Int* rcvPts, E_Int& nbRcvPts,
					     E_Int& ideb     , E_Int& ifin  , E_Int& ithread ,
					     E_Int* ca       , E_Float* w   , E_Float cs     , E_Float nu     ,
					     E_Float* QlocPtr, E_Int varType, E_Int* Qdir    , E_Float gamma_precon_inv,
					     E_Float* H2H3   , E_Int ni     , E_Int ninj     , 
					     vector<E_Float*>& macrosfields, vector<E_Float*>& Qneqfields  , vector<E_Float*>& Qinfields,
					     vector<E_Float*>& cellN_IBC_LBM){

  E_Int nvars = Qneqfields.size();
  
  // variables for commonGeom
  E_Float a0,a1,a2;
  E_Float b0,b1,b2;
  E_Float n0,n1,n2;
  E_Float normb, vn;
  E_Float alpha, beta, alphasbeta;

  // variables for flow varialbes
  E_Float* qOut;
  E_Float* rho_;
  E_Float* rho_u;
  E_Float* rho_v;
  E_Float* rho_w;
  E_Float* qneq;
  E_Float* qeq;
  
  E_Float q_local;
  E_Float rho;
  E_Float Ux;
  E_Float Uy;
  E_Float Uz;
  E_Float qneq_local;
  E_Float qeq_local;
  E_Float udotc,udotu;
  E_Float c02,c04;

  // HRR
  E_Int l,l1,l12,l2,l22,l3,l4,l32,l42,l5,l6,l52,l62;
  E_Int p;
  E_Float sigma,um_sigma,coef_coll,c0,c1,c2,c_fd;
  E_Float coef  ,coef1 ,coef2 ,coef3 ,coef4 ,coef5 ,coef6;
  E_Float coef7 ,coef8 ,coef9 ,coef10,coef11,coef12,coef13;
  E_Float coef14,coef15,coef16,coef17  ;
  E_Float a1fd_1,a1fd_2,a1fd_3,a1fd_4,a1fd_5;
  E_Float a1fd_6,a1fd_7,a1fd_8,a1fd_9;
  E_Float a1hrr_1,a1hrr_2,a1hrr_3,a1hrr_4,a1hrr_5;
  E_Float a1hrr_6,a1hrr_7,a1hrr_8,a1hrr_9;
  E_Float a1pr_1,a1pr_2,a1pr_3,a1pr_4,a1pr_5;
  E_Float a1pr_6,a1pr_7,a1pr_8,a1pr_9;
  E_Float a3neq_1,a3neq_2,a3neq_3;
  E_Float a3neq_4,a3neq_5,a3neq_6,a3neq_7;
  E_Float a4neq_1,a4neq_2,a4neq_3;
  E_Float a4neq_4,a4neq_5,a4neq_6;
  E_Float a5neq_1,a5neq_2,a5neq_3;
  E_Float a6neq_1;
  E_Float fd_prt1 , fd_prt2 , fd_prt3;
  E_Float fd_prt12, fd_prt22, fd_prt32;
  E_Float fd_a,fd_b,fd_c;
  E_Float* cellNIBCx;
  E_Float* cellNIBCy;
  E_Float* cellNIBCz;
  
  rho_  = macrosfields[0];
  rho_u = macrosfields[1];
  rho_v = macrosfields[2];
  rho_w = macrosfields[3];

  cellNIBCx=cellN_IBC_LBM[0];
  cellNIBCy=cellN_IBC_LBM[1];
  cellNIBCz=cellN_IBC_LBM[2];
  
  c02 =1./(cs*cs);
  c04 =c02*c02; 
  // no slip
  if (bctype == 1){
#ifdef _OPENMP4
#pragma omp simd
#endif    
    for (E_Int noind = 0; noind < ifin-ideb; noind++){
      E_Int indR = rcvPts[noind+ideb];
      if(Qdir[0+noind*nvars] == 1){
	rho  = 1./rho_[indR];
	Ux   = rho_u[indR]*rho;
	Uy   = rho_v[indR]*rho;
	Uz   = rho_w[indR]*rho;
	udotu =Ux*Ux+Uy*Uy+Uz*Uz;      
	for (E_Int f_index = 0; f_index <nvars; f_index++){
	  qneq        = Qneqfields[f_index];
	  qOut        = Qinfields[f_index];

	  qneq_local  = qneq[indR];

	  udotc =(Ux*ca[f_index         ]+
		  Uy*ca[f_index+nvars   ]+
		  Uz*ca[f_index+2*+nvars]);

	  qeq_local= w[f_index]*rho_[indR]*(1.+
					    c02*udotc+
					    0.5*c04*udotc*udotc*gamma_precon_inv-
					    0.5*c02*udotu*gamma_precon_inv);
	
	  
	  qOut[indR] = qeq_local+qneq_local;
	  //BELOW IS NOT CORRECT
	  //QlocPtr[noind+ideb] = qOut[indR];
	}
      }      
    }
  }
  else if (bctype == 11){
    E_Int shift_h3 =nvars*9;
    E_Int shift_h3b=nvars*(9+7);
    E_Float taug   = nu*c02+0.5;
    E_Float offeq;
    
    sigma     = 0.99999995;
    um_sigma  = 1. - sigma; 
    coef      = -2.*taug*cs*cs;  
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++){
      E_Int indR = rcvPts[noind+ideb];
      l   = indR;
      l1  = indR + 1;
      l2  = indR - 1;
      l12 = indR + 2;
      l22 = indR - 2;
      l3  = indR + ni;
      l4  = indR - ni;
      l32 = indR + 2*ni;
      l42 = indR - 2*ni;
      l5  = indR + ninj;
      l6  = indR - ninj;
      l52 = indR + 2*ninj;
      l62 = indR - 2*ninj;
      if(Qdir[0+noind*nvars] == 1){
	rho  = 1./rho_[indR];
	Ux   = rho_u[indR]*rho;
	Uy   = rho_v[indR]*rho;
	Uz   = rho_w[indR]*rho;
	udotu =Ux*Ux+Uy*Uy+Uz*Uz;
# include "IBC/LBM/HRR_stresstensor.h"
	for (E_Int f_index = 0; f_index <nvars; f_index++){
	  qeq        = Qneqfields[f_index];

	  udotc =(Ux*ca[f_index         ]+
		  Uy*ca[f_index+nvars   ]+
		  Uz*ca[f_index+2*+nvars]);

	  coef1 = Ux*Ux*Uy + Uy*Uz*Uz;
	  coef2 = Ux*Uz*Uz + Ux*Uy*Uy;
	  coef3 = Uy*Uy*Uz + Ux*Ux*Uz;
	  coef4 = Ux*Ux*Uy - Uy*Uz*Uz;
	  coef5 = Ux*Uz*Uz - Ux*Uy*Uy;
	  coef6 = Uy*Uy*Uz - Ux*Ux*Uz;

	  qeq[indR]= w[f_index]*rho_[indR]*(1.+
					     c02*udotc+
					     0.5*c04*udotc*udotc-
					     0.5*c02*udotu+
					     13.5*(H2H3[        shift_h3]*coef1 + 
						   H2H3[1*nvars+shift_h3]*coef2 +
						   H2H3[2*nvars+shift_h3]*coef3)+
					     4.5 *(H2H3[3*nvars+shift_h3]*coef4 +
						   H2H3[4*nvars+shift_h3]*coef5 +
						   H2H3[5*nvars+shift_h3]*coef6));
	
	    }
# include "IBC/LBM/HRR_H2.h"
	for (E_Int f_index = 0; f_index <nvars; f_index++){
	  qeq  = Qneqfields[f_index];
	  qOut = Qinfields[f_index];
	  qneq_local = w[f_index]*(4.5*(H2H3[(1-1)*nvars+f_index]*a1hrr_1    +
				  H2H3[(2-1)*nvars+f_index]*a1hrr_2+H2H3[(3-1)*nvars+f_index]*a1hrr_3 +
				  H2H3[(4-1)*nvars+f_index]*a1hrr_4+H2H3[(5-1)*nvars+f_index]*a1hrr_5 +
				  H2H3[(6-1)*nvars+f_index]*a1hrr_6+H2H3[(7-1)*nvars+f_index]*a1hrr_7 +
				  H2H3[(8-1)*nvars+f_index]*a1hrr_8+H2H3[(9-1)*nvars+f_index]*a1hrr_9)+
			     13.5*(H2H3[(1-1)*nvars+f_index+shift_h3b]*coef1   +
				   H2H3[(2-1)*nvars+f_index+shift_h3b]*coef2   +
				   H2H3[(3-1)*nvars+f_index+shift_h3b]*coef3)  +
			     4.5*(H2H3[(4-1)*nvars+f_index+shift_h3b]*coef4    +
				  H2H3[(5-1)*nvars+f_index+shift_h3b]*coef5    +
				  H2H3[(6-1)*nvars+f_index+shift_h3b]*coef6)   );
	  qOut[indR] = qeq[indR]+qneq_local;
	}
      }
    }      
  }
  else{
    printf("Warning !!! setIBCTransfersCommonVar5: bcType " SF_D_ " not implemented.\n", bctype);
    return 0;
  }
  return 1;
}



// IBC for penalization
E_Int K_CONNECTOR::setIBCTransfersCommonVar44(E_Int* rcvPts, E_Int& nbRcvPts, E_Int& ideb     , 
					      E_Int& ifin  , E_Int& ithread                   ,
					      E_Float* xPC , E_Float* yPC   , E_Float* zPC    ,
					      E_Int* ca    , E_Int* Qdir    , E_Int* cminus   ,
					      E_Int* ptbcs , E_Int num_bcs  , E_Float* QlocPtr, 
					      E_Float* tmp , E_Int& size    , E_Float& meax   ,
					      E_Float& meay, E_Float& meaz  ,
					      vector<E_Float*>& vectOfRcvFields,
					      vector<E_Float*>& vectOfQstarRcvFields){
  
  E_Int nvars = vectOfRcvFields.size();
  E_Int ii,jj;
  
  // variables for flow varialbes
  E_Float* qin;
  E_Float* qout;
  // no slip 
  for (E_Int noind = 0; noind < ifin-ideb; noind++){
    E_Int indR = rcvPts[noind+ideb];
    for (E_Int f_index = 1; f_index <nvars; f_index++){
      if(Qdir[f_index+noind*nvars] == 1){
	jj  = cminus[f_index];
	qin = vectOfQstarRcvFields[f_index];
	qout= vectOfRcvFields[jj];
	qout[indR] = qin[indR];
	// Momentum Exchange Algorithm
# include "IBC/LBM/mea_calc.h"	
      }
    }
  }
  return 1;
}

E_Int K_CONNECTOR::setIBCTransfersCommonVar45(E_Int* rcvPts, E_Int& nbRcvPts, E_Int& ideb     , 
					      E_Int& ifin  , E_Int& ithread ,		      
					      E_Float* xPC , E_Float* yPC   , E_Float* zPC    ,
					      E_Int* ca    , E_Int* Qdir    , E_Float* Qdist  ,
					      E_Int* cminus, E_Int& imd     , E_Int& jmd      ,
					      E_Int* ptbcs , E_Int num_bcs  , E_Float* QlocPt , 
					      E_Float* tmp , E_Int& size    , E_Float& meax   ,
					      E_Float& meay, E_Float& meaz  ,
					      vector<E_Float*>& vectOfRcvFields,
					      vector<E_Float*>& vectOfQstarRcvFields){
  E_Int nvars = vectOfRcvFields.size();
  E_Int ii,jj;
  E_Int shift;
  // variables for commonGeom
  E_Float a0,a1,a2;
  E_Float b0,b1,b2;
  E_Float n0,n1,n2;
  E_Float normb, vn;
  E_Float alpha, beta, alphasbeta;

  // variables for flow varialbes
  E_Int shift_xf;
  E_Int indRf;
  E_Float* qin;
  E_Float* qin2;
  E_Float* qout;
  E_Float qdist_local;
  // no slip
  for (E_Int noind = 0; noind < ifin-ideb; noind++){
    E_Int indR  = rcvPts[noind+ideb];    
    for (E_Int f_index = 1; f_index <nvars; f_index++){
      if(Qdir[f_index+noind*nvars] == 1){
	jj = cminus[f_index];
	shift_xf =ca[f_index        ]    +
	          ca[f_index+  nvars]*imd +
	          ca[f_index+2*nvars]*imd*jmd;	
	indRf = rcvPts[noind+ideb]-shift_xf;
	qin   = vectOfQstarRcvFields[f_index];
	qin2  = vectOfQstarRcvFields[jj];
	qout  = vectOfRcvFields[jj];	
	qdist_local = Qdist[f_index+noind*nvars] ;
	if(qdist_local<0.5){
	  qout[indR] = 2.*qdist_local*qin[indR]+(1.-2*qdist_local)*qin[indRf];
	}
	else{
	  qout[indR] = 1./(2.*qdist_local)*(qin[indR]+(2*qdist_local-1.)*qin2[indR]);
	}
	// Momentum Exchange Algorithm
# include "IBC/LBM/mea_calc.h"
      }      
    }
  }

  return 1;
}

E_Int K_CONNECTOR::setIBCTransfersCommonVar46(E_Int* rcvPts             , E_Int& nbRcvPts, E_Int& ideb     , 
					      E_Int& ifin               , E_Int& ithread ,
					      E_Float* xPC              , E_Float* yPC   , E_Float* zPC    ,
					      E_Int* ca                 , E_Int* Qdir    , E_Float* Qdist  ,
					      E_Int* cminus             , E_Int& imd     , E_Int& jmd      ,
					      E_Float& tau_relax        , E_Float* w     , E_Float& cs     ,
					      E_Float& gamma_precond_inv, E_Int& ni      , E_Int& ninj    ,
					      E_Int* ptbcs              , E_Int num_bcs  , E_Float* QlocPtr, 
					      E_Float* tmp              , E_Int& size    , E_Float& meax   ,
					      E_Float& meay             , E_Float& meaz  ,
					      vector<E_Float*>& vectOfRcvFields,
					      vector<E_Float*>& vectOfQstarRcvFields,
					      vector<E_Float*>& vectOfmacroRcvFields){

  E_Int nvars    = vectOfRcvFields.size();
  E_Float* rho_  = vectOfmacroRcvFields[0];// ro
  E_Float* rho_u = vectOfmacroRcvFields[1];// u
  E_Float* rho_v = vectOfmacroRcvFields[2];// v
  E_Float* rho_w = vectOfmacroRcvFields[3];// w
  
  E_Int ii,jj;
  E_Int shift, shift_xf;
  // variables for commonGeom
  E_Float a0,a1,a2;
  E_Float b0,b1,b2;
  E_Float n0,n1,n2;
  E_Float normb, vn;
  E_Float alpha, beta, alphasbeta;

  // variables for flow varialbes
  E_Float* qin;
  E_Float* qin2;
  E_Float* qout;
  E_Float qdist_local;
  E_Float chi_local;
  E_Float c02 = 1./(cs*cs);
  E_Float c04 = c02*c02;
  E_Float feq_local;
  E_Float udotu;
  E_Float udotc_s;
  E_Float udotc_b;
  E_Float rho  , Ux  , Uy  , Uz  ;
  E_Float rho_f, Ux_f, Uy_f, Uz_f;
  E_Float u_s,v_s,w_s;

  // no slip
  for (E_Int noind = 0; noind < ifin-ideb; noind++){
    E_Int indR  = rcvPts[noind+ideb];
    E_Float tau_freq = 1./tau_relax;

    rho=rho_ [indR];
    Ux =rho_u[indR]/rho;
    Uy =rho_v[indR]/rho;
    Uz =rho_w[indR]/rho;

    for (E_Int f_index = 1; f_index <nvars; f_index++){
      if(Qdir[f_index+noind*nvars] == 1){
	jj = cminus[f_index];
	
	shift_xf=(ca[jj        ]   +
		  ca[jj+  nvars]*ni+
		  ca[jj+2*nvars]*ninj);

	rho_f=rho_ [indR+shift_xf];
	Ux_f =rho_u[indR+shift_xf]/rho_f;
	Uy_f =rho_v[indR+shift_xf]/rho_f;
	Uz_f =rho_w[indR+shift_xf]/rho_f;
	
	qin   = vectOfQstarRcvFields[f_index];
	qin2  = vectOfQstarRcvFields[jj];
	qout  = vectOfRcvFields[jj];
	
	qdist_local = Qdist[f_index+noind*nvars];	
	
	if(qdist_local<0.5){
	  // Original FH
	  chi_local  = tau_freq*(2.*qdist_local-1.)/(1.-tau_freq);
	  u_s        = Ux;
	  v_s        = Uy;
	  w_s        = Uz;
	  
	  // Mei,Luo & Shyy Improvement
	  chi_local  = tau_freq*(2.*qdist_local-1.)/(1.-2.*tau_freq);
	  u_s        = Ux_f;
	  v_s        = Uy_f;
	  w_s        = Uz_f;
	}
	else{
	  chi_local  = tau_freq*(2.*qdist_local-1.);
	  u_s        = 1./qdist_local*((qdist_local-1.)*Ux);
	  v_s        = 1./qdist_local*((qdist_local-1.)*Uy);
	  w_s        = 1./qdist_local*((qdist_local-1.)*Uz);
	}
	
	udotc_s=( ca[f_index        ]*u_s +
		  ca[f_index+  nvars]*v_s +
		  ca[f_index+2*nvars]*w_s );
	udotc_b=( ca[f_index        ]*Ux +
		  ca[f_index+  nvars]*Uy +
		  ca[f_index+2*nvars]*Uz );
	udotu=( Ux*Ux + Uy*Uy + Uz*Uz );
	
	feq_local = w[f_index]*rho*(1.0                     +
				    udotc_s*c02             +
				    0.5*udotc_b*udotc_b*c04*gamma_precond_inv -
				    0.5*udotu*c02*gamma_precond_inv );

	  
	qout[indR] = (1.-chi_local)*qin[indR]+chi_local*feq_local;
	// Momentum Exchange Algorithm
#include "IBC/LBM/mea_calc.h"
      }
    }
  } 
  return 1;
}


E_Int K_CONNECTOR::setIBCTransfersCommonVar47(E_Int* rcvPts     , E_Int& nbRcvPts, E_Int& ideb     , 
					      E_Int& ifin       , E_Int& ithread ,
					      E_Float* xPC      , E_Float* yPC   , E_Float* zPC    ,
					      E_Int* ca         , E_Int* Qdir    , E_Float* Qdist  ,
					      E_Int* cminus     , E_Int& imd     , E_Int& jmd      ,
					      E_Float& tau_relax, E_Float* w     , E_Float& cs     ,
					      E_Int* ptbcs      , E_Int num_bcs  , E_Float* QlocPtr, 
					      E_Float* tmp      , E_Int& size    , E_Float& meax   ,
					      E_Float& meay     , E_Float& meaz  ,
					      vector<E_Float*>& vectOfRcvFields,
					      vector<E_Float*>& vectOfQstarRcvFields){

  E_Int nvars      = vectOfRcvFields.size();
  E_Int ii,jj;
  // variables for flow variables
  E_Float* qin;
  E_Float* qin2;
  E_Float* qout;
  E_Float* qout2;
  E_Float qdist;
  // no slip
  for (E_Int noind = 0; noind < ifin-ideb; noind++){
    E_Int indR  = rcvPts[noind+ideb];    
    for (E_Int f_index = 1; f_index <nvars; f_index++){
      if(Qdir[f_index+noind*nvars] == 1){
	jj    = cminus[f_index];
	
	qin   = vectOfQstarRcvFields[f_index];
	qin2  = vectOfQstarRcvFields[jj];	
	qout  = vectOfRcvFields[jj];
	qout2 = vectOfRcvFields[f_index];
	
	qdist = Qdist[f_index+noind*nvars] ;
	qout[indR]  = 1./(1.+2.*qdist)*((2.*qdist)*qin2[indR]+qout2[indR]);
	// Momentum Exchange Algorithm
#include "IBC/LBM/mea_calc.h"
      }
    }
  }
  return 1;
}




// Determine Q's that intersect immersed boundary for penalization
E_Int K_CONNECTOR::setIBCTransfersCommonVarQtagonly(E_Int* rcvPts, E_Int& nbRcvPts, E_Int& ideb,
						    E_Int& ifin  , E_Int& ithread ,
						    E_Float* xPC , E_Float* yPC   , E_Float* zPC,
						    E_Int* cvel  , E_Int* Qdir    , E_Int& nvars,
						    E_Int& imd   , E_Int& jmd     , E_Float* tmp,
						    E_Int& size  ,
						    E_Float& zlimit, E_Int& isinside,
						    E_Float*& RcvFields){

  E_Int imdjmd = imd*jmd;
  E_Int is,js,ks;
  E_Int indR2;
  E_Float eps=1e-10;
  if(isinside==0){
    for (E_Int noind = 0; noind < ifin-ideb; noind++){
      E_Int indR  = rcvPts[noind+ideb];      
      for (E_Int f_index = 1; f_index <nvars; f_index++){	
	is = cvel[f_index        ];
	js = cvel[f_index+  nvars];
	ks = cvel[f_index+2*nvars];
	indR2 = indR + is + js*imd + ks*imdjmd;
	Qdir[f_index+noind*nvars] = 0;
	if(abs(zPC[noind+ideb])<zlimit){	  
	  if(RcvFields[indR]*RcvFields[indR2]<0+eps && RcvFields[indR]>0-eps){
	    Qdir[f_index+noind*nvars] = 1;
	  }
	}
      }      
    }
  }
  else{
    E_Int f_index=0;
    for (E_Int noind = 0; noind < ifin-ideb; noind++){
      E_Int indR  = rcvPts[noind+ideb];
      Qdir[f_index+noind*nvars] = 0;
      if(RcvFields[indR]<0 && abs(zPC[noind+ideb])<zlimit){
	Qdir[f_index+noind*nvars] = 1;      
      }      
    }
  }
  return 1;
}


// Determine Q's that intersect immersed boundary for ghost method approach
E_Int K_CONNECTOR::setIBCTransfersCommonVarQ(E_Int* rcvPts  , E_Int& nbRcvPts, E_Int& ideb   ,
					     E_Int& ifin    , E_Int& ithread ,
					     E_Float* xPC   , E_Float* yPC   , E_Float* zPC  ,
					     E_Int* cvel    , E_Int* Qdir    , E_Float* Qdist,
					     E_Int& nvars   , E_Int& imd     , E_Int& jmd    ,
					     E_Float* tmp   , E_Int& size    ,
					     E_Float& zlimit,
					     E_Float*& RcvFields){
  E_Float distmp;
  E_Int imdjmd = imd*jmd;
  E_Int is,js,ks;
  E_Int indR2;
  E_Float eps=1e-17;
  for (E_Int noind = 0; noind < ifin-ideb; noind++){
    E_Int indR  = rcvPts[noind+ideb];
    if(RcvFields[indR]>0+eps && abs(zPC[noind+ideb])<zlimit){
      Qdir[0+noind*nvars] = 1;
    }
    for (E_Int f_index = 1; f_index <nvars; f_index++){
      is = cvel[f_index        ];
      js = cvel[f_index+  nvars];
      ks = cvel[f_index+2*nvars];
      indR2 = indR + is + js*imd + ks*imdjmd;
      Qdir [f_index+noind*nvars] = 0;
      Qdist[f_index+noind*nvars] = 0.;
      
      if(abs(zPC[noind+ideb])<zlimit){
	if(RcvFields[indR]*RcvFields[indR2]<0.+eps &&
	 RcvFields[indR]>0+eps){
	  Qdir[f_index+noind*nvars] = 1;
	  distmp     = RcvFields[indR]/(RcvFields[indR]-RcvFields[indR2]);
	  Qdist[f_index+noind*nvars] = distmp;
	}
      }
    }
  }
  return 1;
}



//// [AJ] NEED TO CHECK THIS
//// Regular LBM (Interpolate Qstar,Qneq)
//E_Int K_CONNECTOR::save_target_values_LBM(E_Int* rcvPts, E_Int& nbRcvPts, E_Int& ideb, E_Int& ifin,  E_Int& ithread,
//					  E_Float* QlocPtr, E_Int& save_interp,
//					  E_Float* tmp, E_Int& size,
//					  vector<E_Float*>& vectOfRcvFields){
//
//  E_Int nvars = vectOfRcvFields.size();
//  // variables for flow varialbes
//
//  E_Float* qOut;
//  E_Int count_local;
//  // no slip
//  count_local = 0;
//  for (E_Int noind = 0; noind < ifin-ideb; noind++){
//    E_Int indR = rcvPts[noind+ideb];
//    for (E_Int f_index = 0; f_index <nvars; f_index++){
//      E_Int ptr_indx = noind+ideb+f_index*nvars;
//      qOut                = vectOfRcvFields[f_index];
//      count_local = count_local + 1;
//      if(save_interp==1){
//	QlocPtr[ptr_indx] = count_local;//qOut[indR];
//	printf("%f \n",count_local);
//      }
//      else{
//	printf("%f \n",QlocPtr[ptr_indx]);
//	//qOut[indR]=QlocPtr[ptr_indx];
//      }
//    }      
//  }
//  
//  return 1;
//}



