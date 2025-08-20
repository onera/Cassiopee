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
# include "param_solver.h"
# include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
using namespace std;
using namespace K_FLD;

# include "IBC/commonLaws.h"
# define NUTILDE_FERRARI 2

//=============================================================================
//Retourne -2: incoherence entre meshtype et le type d interpolation
//         -1: type invalide
//          1: ok
// Entree/Sortie: des variables conservatives ( + ronutildeSA )
//=============================================================================
E_Int K_CONNECTOR::setIBCTransfersCommonVar1(
               E_Int bctype,
               E_Int* rcvPts, E_Int& nbRcvPts, E_Int& ideb, E_Int& ifin, E_Int& ithread,
               E_Float* xPC, E_Float* yPC, E_Float* zPC,
               E_Float* xPW, E_Float* yPW, E_Float* zPW,
               E_Float* xPI, E_Float* yPI, E_Float* zPI,
               E_Float* densPtr, E_Float* pressPtr,
               E_Float* vxPtr, E_Float* vyPtr, E_Float* vzPtr,
               E_Float* utauPtr, E_Float* yplusPtr, E_Float* kcurvPtr,
               E_Float* d1, E_Float* d2, E_Float* d3, E_Float* d4, E_Float* d5,
               E_Float* tmp, E_Int& size,
               E_Float gamma, E_Float cv, E_Float muS, E_Float Cs, E_Float Ts, E_Float Pr,
               vector<E_Float*>& vectOfDnrFields, vector<E_Float*>& vectOfRcvFields)
{
  E_Float gam1 = gamma-1.;
  /* lois de paroi */
  E_Float roext, uext, pext, text, muext, yext, yplus, yibc;
  E_Float uscaln, un, vn, wn, ut, vt, wt, utau0, umod;
  E_Float expy, denoml10, ax, l1, l2, l3;
  E_Float fp, tp;
  E_Float ucible0, ucible, vcible, wcible, signibc, twall, rowall, muwall;
  //Lois de paroi : criteres d arret pour estimer le frottement par Newton
  E_Float newtoneps = 1.e-7; // critere d arret pour u+
  //E_Float newtonepsnutilde = 1.e-10; // critere d'arret pour nutilde
  E_Float newtonepsprime = 1.e-12;// critere d'arret pour la derivee
  E_Float cvgaminv = 1./(cv*gam1);
  E_Float coefSuth = muS * (1.+Cs/Ts);
  E_Float Tsinv = 1./Ts;
  E_Float kappa = 0.4; // Constante de Von Karman
  E_Float kappainv = 1./kappa;
  E_Float cc = 5.2;//pour la loi log
  E_Float one_third = 1./3.;
  /* fin parametres loi de parois */

  E_Int nvars = vectOfDnrFields.size();
  //E_Int nbRcvPts = rcvPtsI.getSize();

  E_Float a0,a1,a2,b0,b1,b2,n0,n1,n2;
  E_Float normb, u,v,w, roE;
  E_Float vnc, alpha, beta, alphasbeta;
  E_Float* roOut  = vectOfRcvFields[0];// ro
  E_Float* rouOut = vectOfRcvFields[1];// rou
  E_Float* rovOut = vectOfRcvFields[2];// rov
  E_Float* rowOut = vectOfRcvFields[3];// row
  E_Float* roEOut = vectOfRcvFields[4];// roE
  E_Float* varSAOut = NULL;
  if (nvars == 6) varSAOut = vectOfRcvFields[5]; // ronutildeSA

  //E_Int* rcvPts = rcvPtsI.begin();
  // if ( (bctype==2 || (bctype==3)) && nvars < 6)
  // {
  //   printf("Warning: setIBCTransfersCommonVar1: number of variables (<6) inconsistent with bctype.\n");
  //   return 0;
  // }
  if (bctype == 0) // wallslip
  {
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      E_Int indR = rcvPts[noind+ideb];

      roext  = roOut[indR];     //ro
      roE    = roEOut[indR];    //roE
      E_Float roinv = 1./roext;
      // vitesse
      u = rouOut[indR]*roinv; v = rovOut[indR]*roinv; w = rowOut[indR]*roinv;

# include "IBC/commonBCType0.h"

      rouOut[indR] = ucible*roext; rovOut[indR] = vcible*roext; rowOut[indR] = wcible*roext;
      // pression
      pext            = gam1*(roE-0.5*roext*(u*u+v*v+w*w));//PI
      roEOut[indR]    = pext/gam1+0.5*roext*(ucible*ucible+vcible*vcible+wcible*wcible);
      pressPtr[noind +ideb] = pext;
      densPtr[ noind +ideb] = roext;

      vxPtr[noind+ideb] = ucible;
      vyPtr[noind+ideb] = vcible;
      vzPtr[noind+ideb] = wcible;

      if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;
    }
  }
  else if (bctype == 1) // adherence
  {
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      E_Int indR = rcvPts[noind+ideb];

      roext  = roOut[indR];      //ro
      roE    = roEOut[indR];     //roE
      E_Float roinv = 1./roext;
      // vitesse
      u = rouOut[indR]*roinv; v = rovOut[indR]*roinv; w = rowOut[indR]*roinv;

#     include "IBC/commonBCType1.h"

      rouOut[indR] = ucible*roext; rovOut[indR] = vcible*roext; rowOut[indR] = wcible*roext;
      // pression
      pext         = gam1*(roE-0.5*roext*(u*u+v*v+w*w));//PI
      roEOut[indR] = pext/gam1+0.5*roext*(ucible*ucible+vcible*vcible+wcible*wcible);

      pressPtr[noind +ideb] = pext;
      densPtr[ noind +ideb] = roext;

      vxPtr[noind+ideb] = ucible;
      vyPtr[noind+ideb] = vcible;
      vzPtr[noind+ideb] = wcible;

      //utau pas calcule ni y+
      //
      if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;
    }
  }
  else if (bctype == 2) // loi de paroi en log
  {
#   include "IBC/pointer.h"

    E_Int err  = 0;
    E_Int skip = 0;
    //initilisation parametre geometrique et utau
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        //E_Int indR = rcvPts[noind];
        E_Int indR = rcvPts[noind+ideb];

        roext = roOut[indR]; // densite du point interpole
        roE   = roEOut[indR];    //roE

        E_Float roinv = 1./roext;
        // vitesse du pt ext
        u = rouOut[indR]*roinv;
        v = rovOut[indR]*roinv;
        w = rowOut[indR]*roinv;

        pext = gam1*(roE-0.5*roext*(u*u+v*v+w*w));//PI
        text = pext / roext * cvgaminv;

#       include "IBC/commonLogLaw_init.h"
        // out= utau  et err
      }

      // Newton pour utau
#     include "IBC/commonLogLaw_Newton.h"

      //initialisation Newton SA  + vitesse cible
#if NUTILDE_FERRARI == 0
#     include "IBC/commonLogLaw_cible.h"
#else
#     include "nutilde_Ferrari.h"
#endif
      if (nvars == 6)
      {
        // Newton pour nut
#if NUTILDE_FERRARI == 0
#       include "IBC/nutildeSA_Newton.h"
#endif

        // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
          E_Int indR = rcvPts[noind+ideb];

          // For Post (temperature du point image text, densite et pression du point de paroi: densPtr,pressPtr)
          text = press_vec[noind]/ ro_vec[noind]* cvgaminv;

          twall = text  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
          densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
          pressPtr[noind+ideb]= press_vec[noind ];

          roE   = roEOut[indR];//roE

          E_Float roinv = 1./ro_vec[noind];
          // vitesse du pt ext
          u = rouOut[indR]*roinv;
          v = rovOut[indR]*roinv;
          w = rowOut[indR]*roinv;

          // Densite corrigee (CC-Busemann)
          roOut[indR]  = press_vec[noind ]/tcible_vec[noind]*cvgaminv;

          rouOut[indR] = ucible_vec[noind]*roOut[indR];
          rovOut[indR] = vcible_vec[noind]*roOut[indR];
          rowOut[indR] = wcible_vec[noind]*roOut[indR];

          // roEOut[indR] = roE + 0.5*ro_vec[noind]*( ucible_vec[noind]*ucible_vec[noind] - u*u
          //                                         +vcible_vec[noind]*vcible_vec[noind] - v*v
          //                                         +wcible_vec[noind]*wcible_vec[noind] - w*w );

          roEOut[indR] = roE +   0.5*roOut[indR]*( ucible_vec[noind]*ucible_vec[noind]
                   +vcible_vec[noind]*vcible_vec[noind]
                   +wcible_vec[noind]*wcible_vec[noind]  )
                  - 0.5*ro_vec[noind]*( u*u+v*v+w*w );

          varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind]*ro_vec[noind];                                    //nutilde*signibc

          vxPtr[noind+ideb] = ucible_vec[noind];
          vyPtr[noind+ideb] = vcible_vec[noind];
          vzPtr[noind+ideb] = wcible_vec[noind];
        }
      }
      else //5eq
      {
        // mise a jour des variables
#ifdef _OPENMP4
#pragma omp simd
#endif
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
          E_Int indR = rcvPts[noind+ideb];

          // For Post (temperature du point image text, densite et pression du point de paroi: densPtr,pressPtr)
          text = press_vec[noind]/ ro_vec[noind]* cvgaminv;

          twall = text  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
          densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
          pressPtr[noind+ideb]= press_vec[noind ];

          E_Float roinv = 1./ro_vec[noind];
          // vitesse du pt ext
          u = rouOut[indR]*roinv;
          v = rovOut[indR]*roinv;
          w = rowOut[indR]*roinv;

          // Densite corrigee (CC-Busemann)
          roOut[indR]  = press_vec[noind ]/tcible_vec[noind]*cvgaminv;

          rouOut[indR] = ucible_vec[noind]*roOut[indR];
          rovOut[indR] = vcible_vec[noind]*roOut[indR];
          rowOut[indR] = wcible_vec[noind]*roOut[indR];

          roE          = roEOut[indR];//roE
          roEOut[indR] = roE +   0.5*roOut[indR]*( ucible_vec[noind]*ucible_vec[noind]
                   +vcible_vec[noind]*vcible_vec[noind]
                   +wcible_vec[noind]*wcible_vec[noind]  )
                  - 0.5*ro_vec[noind]*( u*u+v*v+w*w );

          // roEOut[indR] = roE + 0.5*ro_vec[noind]*( ucible_vec[noind]*ucible_vec[noind] - u*u
          //                                         +vcible_vec[noind]*vcible_vec[noind] - v*v
          //                                         +wcible_vec[noind]*wcible_vec[noind] - w*w );

          vxPtr[noind+ideb] = ucible_vec[noind];
          vyPtr[noind+ideb] = vcible_vec[noind];
          vzPtr[noind+ideb] = wcible_vec[noind];
        }
      }
    }
  else if (bctype == 3)// loi de paroi Musker
  {
#   include "IBC/pointer.h"

    E_Int err  = 0;
    E_Int skip = 0;
    //initialisation parametre geometrique et utau
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      //E_Int indR = rcvPts[noind];
      E_Int indR = rcvPts[noind+ideb];

      roext = roOut[indR]; // densite du point interpole
      roE   = roEOut[indR];    //roE

      E_Float roinv = 1./roext;
      // vitesse du pt ext
      u = rouOut[indR]*roinv;
      v = rovOut[indR]*roinv;
      w = rowOut[indR]*roinv;

      pext = gam1*(roE-0.5*roext*(u*u+v*v+w*w));//PI
      text = pext / roext * cvgaminv;

#     include "IBC/commonMuskerLaw_init.h"
      // out= utau  et err
    }

    // Newton pour utau
#   include "IBC/commonMuskerLaw_Newton.h"

    //initialisation Newton SA  + vitesse cible
#if NUTILDE_FERRARI == 0
#   include "IBC/commonMuskerLaw_cible.h"
#else
#   include "IBC/nutilde_Ferrari.h"
#endif
    if (nvars == 6)
    {
      err  = 0; skip = 0;
      // Newton pour nut
#if NUTILDE_FERRARI == 0
#     include "IBC/nutildeSA_Newton.h"
#endif

      // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];

        // For Post (temperature du point image text, densite et pression du point de paroi: densPtr,pressPtr)
        text = press_vec[noind]/ ro_vec[noind]* cvgaminv;

        twall = text  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
        densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
        pressPtr[noind+ideb]= press_vec[noind ];

        roE   = roEOut[indR];//roE

        E_Float roinv = 1./ro_vec[noind];
        // vitesse du pt ext
        u = rouOut[indR]*roinv;
        v = rovOut[indR]*roinv;
        w = rowOut[indR]*roinv;

        // Densite corrigee (CC-Busemann)
        roOut[indR]  = press_vec[noind ]/tcible_vec[noind]*cvgaminv;

        rouOut[indR] = ucible_vec[noind]*roOut[indR];
        rovOut[indR] = vcible_vec[noind]*roOut[indR];
        rowOut[indR] = wcible_vec[noind]*roOut[indR];

        // roEOut[indR] = roE + 0.5*ro_vec[noind]*( ucible_vec[noind]*ucible_vec[noind] - u*u
        //                                         +vcible_vec[noind]*vcible_vec[noind] - v*v
        //                                         +wcible_vec[noind]*wcible_vec[noind] - w*w );
        roEOut[indR] = roE +   0.5*roOut[indR]*( ucible_vec[noind]*ucible_vec[noind]
                   +vcible_vec[noind]*vcible_vec[noind]
                   +wcible_vec[noind]*wcible_vec[noind]  )
                   - 0.5*ro_vec[noind]*( u*u+v*v+w*w );

        varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind]*ro_vec[noind];                                    //nutilde*signibc

        vxPtr[noind+ideb] = ucible_vec[noind];
        vyPtr[noind+ideb] = vcible_vec[noind];
        vzPtr[noind+ideb] = wcible_vec[noind];
      }
    }
    else //5eq
    { // mise a jour des variables
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];

        // For Post (temperature du point image text, densite et pression du point de paroi: densPtr,pressPtr)
        text = press_vec[noind]/ ro_vec[noind]* cvgaminv;

        twall = text  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
        densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
        pressPtr[noind+ideb]= press_vec[noind ];

        roE   = roEOut[indR];//roE

        E_Float roinv = 1./ro_vec[noind];
        // vitesse du pt ext
        u = rouOut[indR]*roinv;
        v = rovOut[indR]*roinv;
        w = rowOut[indR]*roinv;

        // Densite corrigee (CC-Busemann)
        roOut[indR]  = press_vec[noind ]/tcible_vec[noind]*cvgaminv;

        rouOut[indR] = ucible_vec[noind]*roOut[indR];
        rovOut[indR] = vcible_vec[noind]*roOut[indR];
        rowOut[indR] = wcible_vec[noind]*roOut[indR];

        roEOut[indR] = roE +   0.5*roOut[indR]*( ucible_vec[noind]*ucible_vec[noind]
                   +vcible_vec[noind]*vcible_vec[noind]
                   +wcible_vec[noind]*wcible_vec[noind]  )
                  - 0.5*ro_vec[noind]*( u*u+v*v+w*w );

        // roEOut[indR] = roE + 0.5*ro_vec[noind]*( ucible_vec[noind]*ucible_vec[noind] - u*u
        //                                         +vcible_vec[noind]*vcible_vec[noind] - v*v
        //                                         +wcible_vec[noind]*wcible_vec[noind] - w*w );

        vxPtr[noind+ideb] = ucible_vec[noind];
        vyPtr[noind+ideb] = vcible_vec[noind];
        vzPtr[noind+ideb] = wcible_vec[noind];
      }
    }
  }
  else if (bctype == 4) // outpres
  {
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      E_Int indR = rcvPts[noind+ideb];
      densPtr[noind+ideb] = roOut[indR];
      pext = pressPtr[noind+ideb];//pression sur la frontiere
      E_Float roinv = 1./roOut[indR];
      u = rouOut[indR]*roinv;
      v = rovOut[indR]*roinv;
      w = rowOut[indR]*roinv;
      roEOut[indR] = pext/gam1+0.5*roOut[indR]*(u*u+v*v+w*w);
      vxPtr[noind+ideb] = u;
      vyPtr[noind+ideb] = v;
      vzPtr[noind+ideb] = w;
    }
  }
  else
  {
    printf("Warning !!! setIBCTransfersCommonVar1: bcType " SF_D_ " not implemented.\n", bctype);
    return 0;
  }
  return 1;
}

//=============================================================================
//Retourne -2: incoherence entre meshtype et le type d'interpolation
//         -1: type invalide
//          1: ok
// Entree/Sortie:  (ro,u,v,w,t) ( + nutildeSA )
//=============================================================================
E_Int K_CONNECTOR::setIBCTransfersCommonVar2(
               E_Int bctype,
               E_Int* rcvPts, E_Int& nbRcvPts, E_Int& ideb, E_Int& ifin, E_Int& ithread,
               E_Float* xPC, E_Float* yPC, E_Float* zPC,
               E_Float* xPW, E_Float* yPW, E_Float* zPW,
               E_Float* xPI, E_Float* yPI, E_Float* zPI, 
               E_Float* densPtr, 
               E_Float* tmp, E_Int& size, E_Int& nvars, 
               E_Float*  param_real,
               //vector<E_Float*>& vectOfDnrFields, vector<E_Float*>& vectOfRcvFields,
               E_Float** vectOfDnrFields, E_Float** vectOfRcvFields,
               E_Int nbptslinelets, E_Float* linelets, E_Int* indexlinelets)
{
  E_Float Pr           = param_real[ PRANDT ];
  E_Float Ts           = param_real[ TEMP0 ];
  E_Float Tinf         = param_real[ TINF ];
  E_Float Pinf         = param_real[ PINF ];
  E_Float Roinf        = param_real[ ROINF ];
  E_Float Cs           = param_real[ CS ];
  E_Float muS          = param_real[ XMUL0 ];
  E_Float cv           = param_real[ CVINF ];
  E_Float gamma        = param_real[ GAMMA ];
  E_Float K_wire       = param_real[ KWire ];
  E_Float Delta_V_wire = param_real[ DeltaVWire ];
  E_Float Diam_wire    = param_real[ DiameterWire ];
  E_Float Ct_WM        = param_real[ CtWire ];
  E_Float R_gas        = Pinf/(Roinf*Tinf);


  E_Int bctypeLocal; 
  E_Int motionType = param_real[MotionType];
  //[AJ] Keep for now
  //E_Float transpeed[3]    = {param_real[TransSpeed],param_real[TransSpeed+1],param_real[TransSpeed+2]};
  //E_Float axispnt[3]      = {param_real[AxisPnt],param_real[AxisPnt+1],param_real[AxisPnt+2]};
  //E_Float axisvec[3]      = {param_real[AxisVec],param_real[AxisVec+1],param_real[AxisVec+2]};
  //E_Float omg             = param_real[OMG];  

  E_Float cmx,cmy,cmz;
  E_Float kvcmx,kvcmy,kvcmz;
  E_Float tmp_x,tmp_y,tmp_z;
  E_Float uGrid_local,vGrid_local,wGrid_local;
  E_Float normalVelGrid_local;
  E_Int c_ale = K_FUNC::E_max(E_Int(0), K_FUNC::E_min(E_Int(1), motionType));
    
  E_Float* pressPtr = densPtr + 1*nbRcvPts;
  E_Float* vxPtr    = densPtr + 2*nbRcvPts;
  E_Float* vyPtr    = densPtr + 3*nbRcvPts; 
  E_Float* vzPtr    = densPtr + 4*nbRcvPts;

  E_Float* utauPtr = NULL;
  E_Float* yplusPtr = NULL;
  E_Float* kcurvPtr = NULL;

  E_Float* d1 = NULL;
  E_Float* d2 = NULL;
  E_Float* d3 = NULL;
  E_Float* d4 = NULL;
  E_Float* d5 = NULL;

  E_Float* tempPtr      = NULL;
  E_Float* tempExtraPtr = NULL;

  E_Float* gradxPressPtr = NULL;
  E_Float* gradyPressPtr = NULL;
  E_Float* gradzPressPtr = NULL;

  E_Float* gradxUPtr = NULL;
  E_Float* gradyUPtr = NULL;
  E_Float* gradzUPtr = NULL;

  E_Float* gradxVPtr = NULL;
  E_Float* gradyVPtr = NULL;
  E_Float* gradzVPtr = NULL;

  E_Float* gradxWPtr = NULL;
  E_Float* gradyWPtr = NULL;
  E_Float* gradzWPtr = NULL;

  E_Float* motionPtr   = NULL;
  E_Float* transpeedPtrX = NULL;
  E_Float* transpeedPtrY = NULL;
  E_Float* transpeedPtrZ = NULL;
  E_Float* axispntPtrX = NULL;
  E_Float* axispntPtrY = NULL;
  E_Float* axispntPtrZ = NULL;
  E_Float* axisvecPtrX = NULL;
  E_Float* axisvecPtrY = NULL;
  E_Float* axisvecPtrZ = NULL;
  E_Float* omgPtr   = NULL;

  E_Float* y_linePtr = NULL;
  E_Float* u_linePtr = NULL;
  E_Float* nutilde_linePtr = NULL;
  E_Float* psi_linePtr = NULL;
  E_Float* matm_linePtr = NULL;
  E_Float* mat_linePtr = NULL;
  E_Float* matp_linePtr = NULL;
  E_Float* alphasbeta_linePtr = NULL;
  E_Float* index_linePtr = NULL;

  // bctype = 3 for all Musker, SA, MuskerLin, & SALin to avoid adding
  // more bctype conditions in if statements.
  // bctypeLocal will be kept for a flag switch for SA (32), MuskerLin (331), & SALin (332).
  // These are in development and will be added in the near future.
  bctypeLocal = bctype;
  if (bctypeLocal == 32 || bctypeLocal == 331 || bctypeLocal == 332)
  {
    bctype=3;
  }
  
  if (motionType==3)
  {
    E_Int shift_var=0;
    // log, Musker, TBLE, MuskerMob, Pohlhausen, Thwaites - also have utau & yplus - need the shift
    if (bctype == 2 || bctype == 3 || bctype == 6 || bctype == 7 || bctype == 8 || bctype == 9) shift_var=2;
      
    motionPtr = densPtr + (14+shift_var)*nbRcvPts;

    transpeedPtrX=densPtr + (15+shift_var)*nbRcvPts;
    transpeedPtrY=densPtr + (16+shift_var)*nbRcvPts;
    transpeedPtrZ=densPtr + (17+shift_var)*nbRcvPts;

    axispntPtrX=densPtr + (18+shift_var)*nbRcvPts;
    axispntPtrY=densPtr + (19+shift_var)*nbRcvPts;
    axispntPtrZ=densPtr + (20+shift_var)*nbRcvPts;

    axisvecPtrX=densPtr + (21+shift_var)*nbRcvPts;
    axisvecPtrY=densPtr + (22+shift_var)*nbRcvPts;
    axisvecPtrZ=densPtr + (23+shift_var)*nbRcvPts;

    omgPtr = densPtr + (24+shift_var)*nbRcvPts;
  }

  if (bctype == 11)
  {
    nbptslinelets = param_real[ NBPTS_LINELETS ];
  }

  if (nbptslinelets == 0 && bctype == 11) //TBLE_FULL -> Musker
  {
    bctype = 3; //Musker
  }

  if ( bctype == 0 || bctype == 1 || bctype == 4 || bctype == 140 || bctype == 141)
  {;}
  else if (bctype==100)//slip + curvature radius
  {
    kcurvPtr = densPtr+5*nbRcvPts;
  }
  else if (bctype == 2 || bctype == 3 || bctype == 6 || bctype == 7 || bctype == 8 || bctype == 9)// log, Musker, TBLE, MuskerMob, Pohlhausen, Thwaites
  {
    utauPtr  = densPtr+5*nbRcvPts;
    yplusPtr = densPtr+6*nbRcvPts;
  }
  else if (bctype == 5)// injection
  {
    d1 = densPtr+5*nbRcvPts;
    d2 = densPtr+6*nbRcvPts;
    d3 = densPtr+7*nbRcvPts;
    d4 = densPtr+8*nbRcvPts;
    d5 = densPtr+9*nbRcvPts;
  }
  else if (bctype == 10) // Mafzal
  {
    utauPtr       = densPtr+5*nbRcvPts;
    yplusPtr      = densPtr+6*nbRcvPts;

    gradxPressPtr = densPtr+7*nbRcvPts;
    gradyPressPtr = densPtr+8*nbRcvPts;
    gradzPressPtr = densPtr+9*nbRcvPts;

    //E_Int   mafzalMode    = param_real[ MAFZAL_MODE ];
    //E_Float alphaGradP    = param_real[ ALPHAGRADP ];
    //nbptslinelets         = param_real[ NBPTS_LINELETS ];
    // std::cout << "mafzalMode = " << mafzalMode << " alpha = " << alphaGradP << " nbpts linelets = " << nbptslinelets << std::endl;
  }
  else if (bctype == 11) // TBLE-FULL
  {
    utauPtr       = densPtr+5*nbRcvPts;
    yplusPtr      = densPtr+6*nbRcvPts;

    gradxPressPtr = densPtr+7*nbRcvPts;
    gradyPressPtr = densPtr+8*nbRcvPts;
    gradzPressPtr = densPtr+9*nbRcvPts;

    gradxUPtr = densPtr+10*nbRcvPts;
    gradyUPtr = densPtr+11*nbRcvPts;
    gradzUPtr = densPtr+12*nbRcvPts;

    gradxVPtr = densPtr+13*nbRcvPts;
    gradyVPtr = densPtr+14*nbRcvPts;
    gradzVPtr = densPtr+15*nbRcvPts;

    gradxWPtr = densPtr+16*nbRcvPts;
    gradyWPtr = densPtr+17*nbRcvPts;
    gradzWPtr = densPtr+18*nbRcvPts;

    y_linePtr          = densPtr+nbRcvPts*(19+nbptslinelets*0);
    u_linePtr          = densPtr+nbRcvPts*(19+nbptslinelets*1);
    nutilde_linePtr    = densPtr+nbRcvPts*(19+nbptslinelets*2);
    psi_linePtr        = densPtr+nbRcvPts*(19+nbptslinelets*3);
    matm_linePtr       = densPtr+nbRcvPts*(19+nbptslinelets*4);
    mat_linePtr        = densPtr+nbRcvPts*(19+nbptslinelets*5);
    matp_linePtr       = densPtr+nbRcvPts*(19+nbptslinelets*6);
    alphasbeta_linePtr = densPtr+nbRcvPts*(19+nbptslinelets*7);
    index_linePtr      = densPtr+nbRcvPts*(19+nbptslinelets*7+1);
  }
  else if (bctype == 12 || bctype == 13 )// isothermal or heat flux
  {
    tempPtr      = densPtr+5*nbRcvPts;
    tempExtraPtr = densPtr+6*nbRcvPts;
  }
  else 
  {
    printf("Warning !!! setIBCTransfersCommonVar2: bcType " SF_D_ " not implemented.\n", bctype);
    return 0;
  }

  /* lois de paroi */
  E_Float roext, uext, pext, text, muext, yext, yplus, yibc, eta, delta;
  E_Float uext_wall,uext_image;
  E_Float uscaln_wall,uscaln_image;
  E_Float gradxPext, gradyPext, gradzPext;
  E_Float gradxUext, gradyUext, gradzUext;
  E_Float gradxVext, gradyVext, gradzVext;
  E_Float gradxWext, gradyWext, gradzWext;
  E_Float uscaln, un, vn, wn, ut, vt, wt, utau, utauv, utau0, umod;  
  E_Float utinf, vtinf, wtinf;
  E_Float t_thwaites, a_thwaites, b_thwaites, c_thwaites, lambda_thwaites, m_thwaites;
  E_Float ut_thwaites, vt_thwaites, wt_thwaites, utau_thwaites, x_T, alpha_T;
  E_Float ut_transition, vt_transition, wt_transition;
  E_Float aa, bb, dd, fp, tp, f1v;
  E_Float expy, denoml10, ax, px, kx, y1, y2, l1, l2, l3, l4, l5, l11, l12, l13;
  E_Float ag11, ag12, ag13, bg10, bg11, bg12, bg13, bg14, cg10, cg11, cg12, cg13, cg14;
  E_Float ag22, ag23, bg20, bg21, bg22, bg23, bg24, cg20, cg21, cg22, cg23, cg24;
  E_Float ucible0, ucible, vcible, wcible, tcible, nutilde, signibc, twall, rowall, muwall;
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

  //E_Int nvars    = vectOfDnrFields.size();
  //E_Int nvarsRcv = vectOfRcvFields.size();

  E_Float a0,a1,a2,b0,b1,b2,n0,n1,n2;
  E_Float t0,t1,t2;
  E_Float normb, ro, u, v, w, t;
  E_Float vnc, alpha, beta, alphasbeta;
  E_Float* roOut = vectOfRcvFields[0];// ro
  E_Float* uOut  = vectOfRcvFields[1];// u
  E_Float* vOut  = vectOfRcvFields[2];// v
  E_Float* wOut  = vectOfRcvFields[3];// w
  E_Float* tOut  = vectOfRcvFields[4];// temperature


  E_Float* varSAOut = NULL;
  
  //---------------------------------
  // ODE-based wall model
  //---------------------------------
  E_Float Cv1cube = pow(7.1,3);
  E_Int nmax      = 20;
  
  E_Float L2norm ;
  E_Float L2norm0;

  E_Float ynm,ynp,dy,dym,dyp,nutm,nutp;
  E_Float xim,xi,xip,m;
  E_Float nutrm,nutrp,nutr;

  E_Float* ipt_u1dold      = NULL;
  E_Float* yline           = NULL;
  E_Float* u_line          = NULL;
  E_Float* nutilde_line    = NULL;
  E_Float* matm_line       = NULL;
  E_Float* mat_line        = NULL;
  E_Float* matp_line       = NULL;
  E_Float* alphasbeta_line = NULL;

  FldArrayF u1dold(nbptslinelets);
  ipt_u1dold = u1dold.begin();

  if (nbptslinelets != 0 && (bctype == 6))
  {
    yline           = linelets;
    u_line          = linelets + nbRcvPts*nbptslinelets;
    nutilde_line    = linelets + nbRcvPts*nbptslinelets*2;
    matm_line       = linelets + nbRcvPts*nbptslinelets*3;
    mat_line        = linelets + nbRcvPts*nbptslinelets*4;
    matp_line       = linelets + nbRcvPts*nbptslinelets*5;
    alphasbeta_line = linelets + nbRcvPts*nbptslinelets*6;
  }
  //---------------------------------

  if (nvars == 6) varSAOut = vectOfRcvFields[5]; // nutildeSA

  // if ( (bctype==2 || (bctype==3)) && nvars < 6)
  // {
  //   printf("Warning: setIBCTransfersCommonVar2: number of variables (<6) inconsistent with bctype.\n");
  //   return 0;
  // }

  if (bctype == 100)//wallslip + curvature radius
  {
#ifdef _OPENMP4
#pragma omp simd
#endif 
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      E_Int indR = rcvPts[noind+ideb];

      // vitesse
      u = uOut[indR];
      v = vOut[indR];
      w = wOut[indR];
#     include "IBC/commonBCType0.h"
      uOut[indR] = ucible; 
      vOut[indR] = vcible; 
      wOut[indR] = wcible;
      if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

      pressPtr[noind + ideb] = roOut[indR]* tOut[indR]*cvgam; 
      densPtr[ noind + ideb] = roOut[indR];

      vxPtr[noind+ideb] = uOut[indR];
      vyPtr[noind+ideb] = vOut[indR];
      vzPtr[noind+ideb] = wOut[indR];
    }
  }
  else if (bctype == 0) // wallslip
  {
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      E_Int indR = rcvPts[noind+ideb];

      // vitesse
      u = uOut[indR];
      v = vOut[indR];
      w = wOut[indR];

      //[AJ]
#     include "IBC/commonIBCmotionAbs2Rel.h"  
    
#     include "IBC/commonBCType0.h"

      uOut[indR] = ucible;
      vOut[indR] = vcible;
      wOut[indR] = wcible;

      //[AJ]
#     include "IBC/commonIBCmotionRel2Abs.h"
    
      if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

      pressPtr[noind + ideb] = roOut[indR]* tOut[indR]*cvgam;
      densPtr[ noind + ideb] = roOut[indR];

      vxPtr[noind+ideb] = uOut[indR];
      vyPtr[noind+ideb] = vOut[indR];
      vzPtr[noind+ideb] = wOut[indR];
    }
  }
  else if (bctype == 1) // adherence (lineaire)
  {
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      E_Int indR = rcvPts[noind+ideb];

      // vitesse
      u = uOut[indR];
      v = vOut[indR];
      w = wOut[indR];

      //[AJ]
#     include "IBC/commonIBCmotionAbs2Rel.h"  	    

#     include "IBC/commonBCType1.h"
    
      uOut[indR] = ucible;
      vOut[indR] = vcible;
      wOut[indR] = wcible;

      //[AJ]
#     include "IBC/commonIBCmotionRel2Abs.h"

      if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

      pressPtr[noind + ideb] = roOut[indR]* tOut[indR]*cvgam;
      densPtr[ noind + ideb] = roOut[indR];

      vxPtr[noind+ideb] = uOut[indR];
      vyPtr[noind+ideb] = vOut[indR];
      vzPtr[noind+ideb] = wOut[indR];

    }
  }
  else if (bctype == 2) // loi de paroi log
  {
#   include "IBC/pointer.h"

    E_Int err  = 0;
    E_Int skip = 0;
    //initialisation parametre geometrique et utau
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      //E_Int indR = rcvPts[noind];
      E_Int indR = rcvPts[noind+ideb];

      roext = roOut[indR]; // densite du point interpole
      text  = tOut[indR];  // pression du point interpole
      pext  = text*roext*cvgam;

      // vitesse du pt ext
      u = uOut[indR];
      v = vOut[indR];
      w = wOut[indR];

#     include "IBC/commonLogLaw_init.h"
      // out= utau  et err
    }

    // Newton pour utau
#   include "IBC/commonLogLaw_Newton.h"

    // initialisation Newton SA  + vitesse cible
#if NUTILDE_FERRARI == 0
#    include "IBC/commonLogLaw_cible.h"
#elif NUTILDE_FERRARI == 1
#    include "IBC/nutilde_Ferrari.h"
#else
#    include "IBC/nutilde_Ferrari_adim.h"
#endif

    if (nvars == 6)
    {
      // Newton pour mut
#if NUTILDE_FERRARI == 0
#     include "IBC/nutildeSA_Newton.h"
#endif
      // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];

        // For Post (tOut temperature du point image en entree, pt corrige en sortie)

        twall = tOut[indR]  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
        densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
        pressPtr[noind+ideb]= press_vec[noind ];

        roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
        uOut[indR]     = ucible_vec[noind]; vOut[indR] = vcible_vec[noind]; wOut[indR]  = wcible_vec[noind];
        tOut[indR]     = tcible_vec[noind];
        varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind];                                         //nutilde*signibc

        vxPtr[noind+ideb] = uOut[indR];
        vyPtr[noind+ideb] = vOut[indR];
        vzPtr[noind+ideb] = wOut[indR];
      }
    }
    else //5eq
    {
      // mise a jour des variables
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];

        roOut[indR]= press_vec[noind ]/tcible_vec[noind]*cvgaminv;

        uOut[indR] = ucible_vec[noind]; vOut[indR] = vcible_vec[noind]; wOut[indR] = wcible_vec[noind];
        tOut[indR] = tcible_vec[noind];

        vxPtr[noind+ideb] = uOut[indR];
        vyPtr[noind+ideb] = vOut[indR];
        vzPtr[noind+ideb] = wOut[indR];
      }
    }
  }
  else if (bctype == 3) // loi de paroi Musker
  {
#   include "IBC/pointer.h" 

    E_Int err  = 0;
    E_Int skip = 0; 
    //initialisation parametre geometrique et utau
#ifdef _OPENMP4
#pragma omp simd
#endif 
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      //E_Int indR = rcvPts[noind];
      E_Int indR = rcvPts[noind+ideb];
 
      roext = roOut[indR]; // densite @ Image Pnt
      text  = tOut[indR];  // pression @ Image Pnt
      pext  = text*roext*cvgam;

      // vitesse @ Image Pnt
      u = uOut[indR];
      v = vOut[indR]; 
      w = wOut[indR];
    
#     include "IBC/commonIBCmotionAbs2Rel.h"
      //Tangential and Normal velocities: need relative velocity 
#     include "IBC/commonMuskerLaw_init.h"
      // out= utau  et err
    }  

    // Newton pour utau
#   include "IBC/commonMuskerLaw_Newton.h" 

    //initialisation Newton SA  + vitesse cible
#if NUTILDE_FERRARI == 0
#   include "IBC/commonMuskerLaw_cible.h"
#elif NUTILDE_FERRARI == 1
#   include "IBC/nutilde_Ferrari.h"
#else
#   include "IBC/nutilde_Ferrari_adim.h"
#endif
    if (nvars == 6)
    {
      // Newton pour mut
#if NUTILDE_FERRARI == 0
#     include "IBC/nutildeSA_Newton.h" 
#endif
      // mise a jour des variables
#ifdef _OPENMP4
#pragma omp simd
#endif 
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];

        // For Post (tOut temperature du point image en entree, pt corrige en sortie)
        twall = tOut[indR] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
        densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
        pressPtr[noind+ideb]= press_vec[noind ];

        // Mise a jour pt corrige
        roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;       
        uOut[indR]     = ucible_vec[noind];
        vOut[indR]     = vcible_vec[noind];
        wOut[indR]     = wcible_vec[noind];
        tOut[indR]     = tcible_vec[noind];
        varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind];  //nutilde*signibc

#       include "IBC/commonIBCmotionRel2Abs.h"  

        vxPtr[noind+ideb] = uOut[indR];
        vyPtr[noind+ideb] = vOut[indR];
        vzPtr[noind+ideb] = wOut[indR];

        // printf("OUT WALL LAW: %f %f %f %f\n",uOut[indR],vOut[indR],wOut[indR],varSAOut[indR]);
      }
    }
    else //5eq 
    {
      // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif 
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];

        // For Post (tOut temperature du point image)
        twall = tOut[indR]  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
        densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
        pressPtr[noind+ideb]= press_vec[noind ];

        roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;   
        uOut[indR]     = ucible_vec[noind];
        vOut[indR]     = vcible_vec[noind];
        wOut[indR]     = wcible_vec[noind];
        tOut[indR]     = tcible_vec[noind];

#       include "IBC/commonIBCmotionRel2Abs.h"  

        vxPtr[noind+ideb] = uOut[indR];
        vyPtr[noind+ideb] = vOut[indR];
        vzPtr[noind+ideb] = wOut[indR];
      }
    }
  }
  else if (bctype == 4) // outpres
  {
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      E_Int indR = rcvPts[noind+ideb];

      if (densPtr[noind+ideb] > 0.) 
      {
        roOut[indR] = densPtr[noind+ideb];
      }
      else
      {
        densPtr[noind+ideb] = -roOut[indR];
      }
      if (pressPtr[noind+ideb] > 0.) 
      {
        tOut[indR] = pressPtr[noind+ideb]/(roOut[indR]*cvgam);//pext/(roext*cvgam)
      }
    
      //tOut[indR] = pressPtr[noind+ideb]/(roOut[indR]*cvgam);//pext/(roext*cvgam)
      //densPtr[noind+ideb] = roOut[indR];
    
      vxPtr[noind+ideb] = uOut[indR];
      vyPtr[noind+ideb] = vOut[indR];
      vzPtr[noind+ideb] = wOut[indR];
    }
  }
  else if (bctype == 5) // inj
  {
    //printf("injection\n");
    if (d1 == NULL)
    {
      printf("Warning: IBC inj with no data.\n"); return 0; // no data
    }

    E_Int indR, newtonmax, nitnwt;
    E_Float d0x, d0y, d0z, gam, gam1, gam2, rgp, tolnewton, c4, c5, c6, c0, c1, c2, c3;
    E_Float roi, ui, vi, wi, Ti, pi, ro0, u0, v0, w0, T0, p0, wni, ri, residug;
    E_Float tnx, tny, tnz, norm, vni, roc0, usd0n, usd0n2, qen, wn0, wng;
    E_Float rog, vg, Tg, pg, ha, pa, b, vng, dwng, f, df;
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      indR = rcvPts[noind+ideb];
      // std::cout<< "********************************************************************" << std::endl;
      // printf("ro,p stored %f %f\n", densPtr[noind+ideb], pressPtr[noind+ideb]);
      // printf("u, v, w: %f %f %f\n", vxPtr[noind+ideb], vyPtr[noind+ideb], vzPtr[noind+ideb]);
      // printf("u,y stored %f %f\n", utauPtr[noind+ideb], yplusPtr[noind+ideb]);
      // printf("state: %f %f %f %f %f\n", d1[noind+ideb], d2[noind+ideb], d3[noind+ideb], d4[noind+ideb], d5[noind+ideb]);
      // std::cout<< "********************************************************************" << std::endl;
      d0x = d3[noind+ideb];
      d0y = d4[noind+ideb];
      d0z = d5[noind+ideb];
    
      ha = d1[noind+ideb];
      pa = d2[noind+ideb];

      gam = gamma;
      gam1 = gam - 1.;
      //gam1_1 = 1. / gam1;
      gam2 = gam / gam1;
      rgp = cv*gam1;
      newtonmax = 40;
      tolnewton = 1.e-6;
      //c4 = 5. / 6.;
      //c5 = 2. / 6.;
      //c6 = -1. / 6.;
      //c0 = 1./c6;
      //c1 =-(c4 + c5)*c0;
      //c2 =- c6*c0;
      //c3 = (2.- c5- c4)*c0;

      // Init newton
      roi = roOut[indR];
      ui  = uOut[indR];
      vi  = vOut[indR];
      wi  = wOut[indR];
      Ti  = tOut[indR];
      pi  = roi*rgp*Ti;

      ro0 = roOut[indR];
      u0  = uOut[indR];
      v0  = vOut[indR];
      w0  = wOut[indR];
      T0  = tOut[indR];
      p0  = ro0*rgp*T0;

      // normale ext normalisee, a inverser en Euler
      tnx = xPI[noind+ideb]-xPW[noind+ideb];
      tny = yPI[noind+ideb]-yPW[noind+ideb];
      tnz = zPI[noind+ideb]-zPW[noind+ideb];
      tnx = -tnx; tny = -tny; tnz = -tnz; // sure?
      norm = sqrt(tnx*tnx+tny*tny+tnz*tnz);
      norm = 1./K_FUNC::E_max(norm, 1.e-12);
      tnx = tnx*norm;
      tny = tny*norm;
      tnz = tnz*norm;

      vni  = ui*tnx + vi*tny + wi*tnz;
      roc0 = sqrt(ro0*gam*p0); // rho*c

      // sans dimension: produit scalaire direction vitesse . normale
      usd0n  = 1./(d0x*tnx + d0y*tny + d0z*tnz);
      usd0n2 = usd0n*usd0n;

      qen = 0.; // a ajouter en ALE

      // ...   Inner caracteristic variable
      // ...   Relative normal velocity
      wni = vni - qen;
      ri  = pi + roc0*wni;

      // ...   Newton Initialization for the relative normal velocity
      wn0  = u0*tnx + v0*tny + w0*tnz  - qen;
      wng  = wn0;

      // resolution Newton
      residug  = 1.e+20;
      nitnwt   = 0;

      while (residug > tolnewton && nitnwt < newtonmax)
      {
        nitnwt  += 1;
        residug = 0.;

        // LINEARISATION a ajouter??
        b = 1. - ((wng+qen)*(wng+qen))*usd0n2/(2.*ha);

        b = K_FUNC::E_max(0.2,b); //cutoff robustesse
        // p = Pi (1 -U^2/(2CpTi))^(gamma/(gamma-1))
        pg  = pa*pow(b,gam2);

        //      nan = isnan(pg)
        //      if(nan)   write(*,*)'fuck Nan inflow_newton',pa(li),b,gam2

        rog = gam2*pg/(ha*b);

        f    = pg + roc0*wng - ri;
        df   = roc0 - rog*(wng+qen)*usd0n2;

        dwng = -f/df;
        wng = wng + dwng;

        residug = K_FUNC::E_max(residug, K_FUNC::E_abs(dwng/wng));
      }

      // LINEARISATION A ajouter??
      vng = (wng+qen);
      //...      Absolute velocity module
      vg  = vng*usd0n;

      b = 1. - (vng*vng)*usd0n2/(2.*ha);
      b = K_FUNC::E_max(0.2,b); //cutoff robustesse

      pg   = pa* pow(b,gam2);
      rog  = gam2*pg/(ha*b);

      Tg   = pg/(rog*rgp);

      roOut[indR] = rog;
      uOut[indR] = vg*d0x;
      vOut[indR] = vg*d0y;
      wOut[indR] = vg*d0z;
      tOut[indR] = Tg;

      // update des grandeurs pour post
      densPtr[noind+ideb] = rog;
      pressPtr[noind+ideb] = pg;
      vxPtr[noind+ideb] = uOut[indR];
      vyPtr[noind+ideb] = vOut[indR];
      vzPtr[noind+ideb] = wOut[indR];
    }
  }
  else if (bctype == 6) // TBLE
  {
    if (nbptslinelets != 0)
    {

#     include "IBC/pointer.h"

      E_Int err  = 0;
      E_Int skip = 0;
      E_Int init = 1;

      if ( (alphasbeta_line[ ideb ] != 0.0) || ideb == ifin)   {init=0;}

      if (init)
      {

        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {

          E_Float* yline1d   = yline + (noind + ideb)*nbptslinelets;

          a0 = xPC[noind+ideb]-xPW[noind+ideb];
          a1 = yPC[noind+ideb]-yPW[noind+ideb];
          a2 = zPC[noind+ideb]-zPW[noind+ideb];

          E_Float dist = sqrt(a0*a0 + a1*a1 + a2*a2);

          E_Int imin  = 0;
          E_Int imax  = nbptslinelets;
          E_Int isearch    = nbptslinelets/2;

          while ( (imax - imin) > 1  )
          {
            if ( dist <= yline1d[isearch])
            {
              imax = isearch;
              isearch = imin + (imax-imin)/2;
            }
            else
            {
              imin = isearch;
              isearch = imin + (imax-imin)/2;
            }
          }

        alphasbeta_line[noind + ideb] = (dist-yline1d[imin] )/(yline1d[imax] -yline1d[imin]);
        indexlinelets[ noind + ideb ] = imax;
      }


      // }
      // Fill nutilde in linelets

      // Test Mixing-length
      for ( E_Int iline = 0 ; iline < nbptslinelets; iline++)
      {
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
          E_Int   indl   = (noind + ideb)*nbptslinelets   + iline;
          E_Int   indR   = rcvPts[noind+ideb];

          roext = roOut[indR]; // densite du point interpole
          text  = tOut[indR];  // pression du point interpole
          pext  = text*roext*cvgam;

          // vitesse du pt ext
          u = uOut[indR];
          v = vOut[indR];
          w = wOut[indR];

#         include "IBC/commonMuskerLaw_init_tble.h"
        }

        // Newton pour utau
#       include "IBC/commonMuskerLaw_Newton.h"

        //initialisation Newton SA  + vitesse cible
#       include "IBC/commonMuskerLaw_cible_tble.h"

        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {

          E_Int   indl   = (noind + ideb)*nbptslinelets   + iline;//*nbptslinelets;
          E_Int   indR   = rcvPts[noind+ideb];

          // Initialisation linelets
          // u_line[indl]     =  uOut[indR]*abs(yline[indl]/yline[nbptslinelets-1]);

          u_line[indl]     =  sqrt(ucible_vec[noind]*ucible_vec[noind] +
                 vcible_vec[noind]*vcible_vec[noind] +
                 wcible_vec[noind]*wcible_vec[noind] );//(uOut[indR] )*abs(alphasbeta) ;//ucible_vec[noind];

          // nutilde_line[indl] = 0.0;//aa_vec[noind]*sign_vec[noind]*uext_vec[noind];  //nu
          nutilde_line[indl] = (varSAOut[indR] )*K_FUNC::E_abs(yline[indl]/yline[nbptslinelets-1]) ;
        }

    }
  } // END INIT


  // 1D TBLE for u (no grad p nor convective terms... for now)  NEED TO CORRECT V and W

  for (E_Int noind = 0; noind < ifin-ideb; noind++)
  {
    // for (E_Int loop = 0; loop < 400; ++loop)
    // {

    E_Int indR = rcvPts[noind+ideb];

    E_Float* yline1d   = yline        + (noind + ideb)*nbptslinelets;
    E_Float* nutilde1d = nutilde_line + (noind + ideb)*nbptslinelets;
    E_Float* u1d       = u_line       + (noind + ideb)*nbptslinelets;
    E_Float* mat       = mat_line     + (noind + ideb)*nbptslinelets;
    E_Float* matm      = matm_line    + (noind + ideb)*nbptslinelets;
    E_Float* matp      = matp_line    + (noind + ideb)*nbptslinelets;

#   include "IBC/commonGeom.h"

    ro_vec[noind] = roOut[indR]; // densite du point interpole
    text  = tOut[indR];  // pression du point interpole
    pext  = text*ro_vec[noind]*cvgam;

    // vitesse du pt ext
    u = uOut[indR];
    v = vOut[indR];
    w = wOut[indR];

    // Loi de Sutherland -> viscosite au point interpole
    mu_vec[noind] = coefSuth * sqrt(K_FUNC::E_abs(text)*Tsinv) / (1.+Cs/text);

    // uscaln: u au point interpole scalaire la normale au pt interpole
    uscaln = u*n0 + v*n1 + w*n2;

    // composante normale de la vitesse
    un = uscaln*n0;
    vn = uscaln*n1;
    wn = uscaln*n2;

    //composante tangentielle de la vitesse au pt interpole
    ut = u-un;
    vt = v-vn;
    wt = w-wn;
    // uext: norme de la composante tangentielle de la vitesse externe
    uext = sqrt(ut*ut+vt*vt+wt*wt);
    uext_vec[noind] = std::max(uext, 1.e-12);

    E_Float nu = mu_vec[noind]/roOut[indR];

    u1d[0] = 0.0;
    u1d[nbptslinelets-1] = uext_vec[noind];

    npass = 0;
    L2norm = 1.0;
    L2norm0= 1.0;

    nmax = 30;


    for (E_Int j = 0; j < nbptslinelets; ++j)
    {
      ipt_u1dold[j] = 0.0;//u1d[j];
    }


    // while ( (L2norm/L2norm0 >= 0.01) && (npass < nmax))
    while ( (L2norm >= 1.0e-1) && (npass < nmax))
    {
#     include "IBC/tble_1d.h"

      // if ((npass == nmax) && (ithread==1))
      // {
      //   std::cout << "TBLE doesn't converge !  " << L2norm << std::endl;
      //   std::cout << "Residual Delta(utble)/utble :" << L2norm << " utble :" << u1d[j_cvg] << " j_cvg :" << j_cvg <<  std::endl;
      // }

      // if (ithread==1)
      // {
      //   std::cout << "Residual utble = " << L2norm << " utble :" << u1d[j_cvg] << " j_cvg :" << j_cvg <<  std::endl;
      // }
      utau_vec[noind ] = sqrt(K_FUNC::E_abs(nu*u1d[1])/yline1d[1]);

      // ODE 1D SA for nutilde

      spalart_1d_(ithread,yline1d, matm,mat,matp, nutilde1d ,u1d,kappa,
            nu,varSAOut[indR],nbptslinelets,kappa);

    }

    if ((npass == nmax) && (ithread==1))
    {
      std::cout << "TBLE doesn't converge !  " << L2norm << std::endl;
      // std::cout << "Residual Delta(utble)/utble :" << L2norm << " utble :" << u1d[j_cvg] << " j_cvg :" << j_cvg <<  std::endl;
    }

  }


  // mise a jour des variables
#ifdef _OPENMP4
#pragma omp simd
#endif
  for (E_Int noind = 0; noind < ifin-ideb; noind++)
  {

    E_Int indR = rcvPts[noind+ideb];

    E_Float* nutilde1d = nutilde_line + (noind + ideb)*nbptslinelets;
    E_Float* u1d       = u_line + (noind + ideb)*nbptslinelets;


    // vitesse du pt ext
    u = uOut[indR];
    v = vOut[indR];
    w = wOut[indR];

    pext  = tOut[indR]*roOut[indR]*cvgam;

#   include "IBC/commonLaws1.h"

    ucible_vec[noind]   = alphasbeta*un; // init : normal component of velocity is linearly reconstructed
    vcible_vec[noind]   = alphasbeta*vn;
    wcible_vec[noind]   = alphasbeta*wn;

    umod = (u1d[indexlinelets[noind + ideb]] - u1d[indexlinelets[noind + ideb]-1])*alphasbeta_line[noind + ideb] + u1d[indexlinelets[noind + ideb]-1];

    ucible0 = signibc/uext * umod;
    ucible_vec[noind] += ucible0 * ut; // vitesse tangentielle pour le pt IBC
    vcible_vec[noind] += ucible0 * vt;
    wcible_vec[noind] += ucible0 * wt;

    // For Post (tOut temperature du point image en entree, pt corrige en sortie)

    twall = tOut[indR]  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext*uext); // Crocco-Busemann()
    densPtr[noind+ideb] = pext/twall*cvgaminv;
    pressPtr[noind+ideb]= pext;

    tcible_vec[noind] = tOut[indR] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext*uext - umod*umod); // Crocco-Busemann

    // Mise a jour pt corrige
    roOut[indR]    = pext/tcible_vec[noind]*cvgaminv;
    uOut[indR]     = ucible_vec[noind];
    vOut[indR]     = vcible_vec[noind];
    wOut[indR]     = wcible_vec[noind];
    tOut[indR]     = tcible_vec[noind];
    varSAOut[indR] = (nutilde1d[indexlinelets[noind + ideb]] - nutilde1d[indexlinelets[noind + ideb]-1])*alphasbeta_line[noind + ideb] + nutilde1d[indexlinelets[noind + ideb]-1]; //aa_vec[noind]*sign_vec[noind]*uext_vec[noind];  //nutilde*signibc

    vxPtr[noind+ideb] = uOut[indR];
    vyPtr[noind+ideb] = vOut[indR];
    vzPtr[noind+ideb] = wOut[indR];

  }


  } // nbptslinelets

  else  // premier call effectue par fillGhostCell (pas de TBLE +SA ) --> on renseigne les PC par du Musker + Mix Length

  {

#   include "IBC/pointer.h"

    E_Int err  = 0;
    E_Int skip = 0;
    //initialisation parametre geometrique et utau
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        //E_Int indR = rcvPts[noind];
        E_Int indR = rcvPts[noind+ideb];

        roext = roOut[indR]; // densite du point interpole
        text  = tOut[indR];  // pression du point interpole
        pext  = text*roext*cvgam;

        // vitesse du pt ext
        u = uOut[indR];
        v = vOut[indR];
        w = wOut[indR];
        //printf("IN WALL LAW: %f %f %f %f %f \n",roext, text, u,v,w);
#       include "IBC/commonMuskerLaw_init.h"
        // out= utau  et err
      }

    // Newton pour utau
#    include "IBC/commonMuskerLaw_Newton.h"

    //initialisation Newton SA  + vitesse cible
#if NUTILDE_FERRARI == 0
#    include "IBC/commonMuskerLaw_cible.h"
#elif NUTILDE_FERRARI == 1
#    include "IBC/nutilde_Ferrari.h"
#else
#    include "IBC/nutilde_Ferrari_adim.h"
#endif
    if (nvars == 6)
      {
        // Newton pour mut
#if NUTILDE_FERRARI == 0
#       include "IBC/nutildeSA_Newton.h"
#endif
        // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      E_Int indR = rcvPts[noind+ideb];

      // For Post (tOut temperature du point image en entree, pt corrige en sortie)

      twall = tOut[indR] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
      densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
      pressPtr[noind+ideb]= press_vec[noind ];

      // Mise a jour pt corrige
      roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
      uOut[indR]     = ucible_vec[noind];
      vOut[indR]     = vcible_vec[noind];
      wOut[indR]     = wcible_vec[noind];
      tOut[indR]     = tcible_vec[noind];
      varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind];  //nutilde*signibc

      vxPtr[noind+ideb] = uOut[indR];
      vyPtr[noind+ideb] = vOut[indR];
      vzPtr[noind+ideb] = wOut[indR];

      // printf("OUT WALL LAW: %f %f %f %f\n",uOut[indR],vOut[indR],wOut[indR],varSAOut[indR]);
    }
      }
    else //5eq
      {
        // mise a jour des variables
#ifdef _OPENMP4
#pragma omp simd
#endif
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      E_Int indR = rcvPts[noind+ideb];

      // For Post (tOut temperature du point image)
      twall = tOut[indR]  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
      densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
      pressPtr[noind+ideb]= press_vec[noind ];

      roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
      uOut[indR]     = ucible_vec[noind];
      vOut[indR]     = vcible_vec[noind];
      wOut[indR]     = wcible_vec[noind];
      tOut[indR]     = tcible_vec[noind];

      vxPtr[noind+ideb] = uOut[indR];
      vyPtr[noind+ideb] = vOut[indR];
      vzPtr[noind+ideb] = wOut[indR];

    }
      }

  } 
    }      
  else if (bctype == 7) // loi de paroi adh paroi rotation
    {
#   include "IBC/pointer.h" 

      E_Int err  = 0;
      E_Int skip = 0; 

      E_Float teta_out = param_real[ROT_TETA];
      E_Float tetap    = param_real[ROT_TETAP];
      //E_Float teta     = teta_out;
      E_Float teta     = 0;


    E_Float cay,caz,ctheta, stheta,vx,vy,vz,vn_paroi;
    stheta = sin(teta);
    ctheta = cos(teta);
#ifdef _OPENMP4
       #pragma omp simd
#endif
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
         E_Int indR = rcvPts[noind+ideb];

         cay = -stheta*(zPW[noind+ideb] - param_real[ROT_CENTER+2]) + ctheta*(yPW[noind+ideb] - param_real[ROT_CENTER+1]);
         caz =  stheta*(yPW[noind+ideb] - param_real[ROT_CENTER+1]) + ctheta*(zPW[noind+ideb] - param_real[ROT_CENTER+2]);
# include "IBC/commonGeom.h"
         vx  =  0;
         vy  = -tetap*caz;
         vz  =  tetap*cay;

         //composante normale de la vitesse paroi
         vn_paroi = vx*n0 + vy*n1 + vz*n2;
  
         //composante tangentielle de la vitesse paroi au pt interpole
         vx = vx-vn_paroi*n0;
         vy = vy-vn_paroi*n1;
         vz = vz-vn_paroi*n2;


         // vitesse relative paroi
         u = uOut[indR]-vx;
         v = vOut[indR]-vy; 
         w = wOut[indR]-vz;
         //if (noind == 0){printf("avt %f %f %f %f \n",vOut[indR],vy,wOut[indR],vz );}
         //
         vn_paroi = u*n0 + v*n1 + w*n2;         

         ucible = (u-vn_paroi*n0)*alphasbeta;// u du pt corrige
         vcible = (v-vn_paroi*n1)*alphasbeta;// v du pt corrige
         wcible = (w-vn_paroi*n2)*alphasbeta;// w du pt corrige


         //E_Float uext2    = vx*vx+vy*vy+vz*vz;
         //E_Float pressure = tOut[indR]*roOut[indR];
         //tOut[indR]    = tOut[indR] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext2); // Crocco-Busemann
         //roOut[indR]   = pressure/tOut[indR];

         uOut[indR] = ucible+vx;
         vOut[indR] = vcible+vy;
         wOut[indR] = wcible+vz;
         //printf("apr %f %f \n",vOut[indR],wOut[indR] );
         if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

         pressPtr[noind + ideb] = roOut[indR]* tOut[indR]*cvgam;
         densPtr[ noind + ideb] = roOut[indR];

         vxPtr[noind+ideb] = uOut[indR];
         vyPtr[noind+ideb] = vOut[indR];
         vzPtr[noind+ideb] = wOut[indR];

        }

    }//bctype 
  else if (bctype == 8) // loi de paroi Pohlhausen
    {
#   include "IBC/pointer.h"

      E_Int err  = 0;
      E_Int skip = 0;
      //initialisation parametre geometrique et utau
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
  {
    //E_Int indR = rcvPts[noind];
    E_Int indR = rcvPts[noind+ideb];

    roext = roOut[indR]; // densite du point interpole
    text  = tOut[indR];  // pression du point interpole
    pext  = text*roext*cvgam;

    // vitesse du pt ext
    u = uOut[indR];
    v = vOut[indR];
    w = wOut[indR];
#       include "IBC/commonPohlhausenLaw_init.h"
  }

#    include "IBC/commonPohlhausenLaw_cible.h"

      if (nvars == 6)
  {
    // mise a jour des variables
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];

        // For Post (tOut temperature du point image en entree, pt corrige en sortie)
        twall = tOut[indR] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
        densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
        pressPtr[noind+ideb]= press_vec[noind ];

        // Mise a jour pt corrige
        roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
        uOut[indR]     = ucible_vec[noind];
        vOut[indR]     = vcible_vec[noind];
        wOut[indR]     = wcible_vec[noind];
        tOut[indR]     = tcible_vec[noind];
        varSAOut[indR] = 0.; // laminar

        vxPtr[noind+ideb] = uOut[indR];
        vyPtr[noind+ideb] = vOut[indR];
        vzPtr[noind+ideb] = wOut[indR];

      }
  }
      else //5eq
  {
    // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];

        // For Post (tOut temperature du point image)
        twall = tOut[indR]  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
        densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
        pressPtr[noind+ideb]= press_vec[noind ];

        roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
        uOut[indR]     = ucible_vec[noind];
        vOut[indR]     = vcible_vec[noind];
        wOut[indR]     = wcible_vec[noind];
        tOut[indR]     = tcible_vec[noind];

        vxPtr[noind+ideb] = uOut[indR];
        vyPtr[noind+ideb] = vOut[indR];
        vzPtr[noind+ideb] = wOut[indR];
      }
  }

    }//bctype
  else if (bctype == 9) // loi de paroi Thwaites
    {
#   include "IBC/pointer.h"

      E_Int err  = 0;
      E_Int skip = 0;
      E_Float* matp;

      //initialisation parametre geometrique et utau
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
  {
    //E_Int indR = rcvPts[noind];
    E_Int indR = rcvPts[noind+ideb];

    roext = roOut[indR]; // densite du point interpole
    text  = tOut[indR];  // pression du point interpole
    pext  = text*roext*cvgam;

    // vitesse du pt ext
    u = uOut[indR];
    v = vOut[indR];
    w = wOut[indR];
#       include "IBC/commonThwaitesLaw_init.h"
    // out= utau  et err
  }

#    include "IBC/commonThwaitesLaw_cible.h"

      if (nvars == 6)
  {
    // mise a jour des variables
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];

        // For Post (tOut temperature du point image en entree, pt corrige en sortie)
        twall = tOut[indR] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
        densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
        pressPtr[noind+ideb]= press_vec[noind ];

        // Mise a jour pt corrige
        roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
        uOut[indR]     = ucible_vec[noind];
        vOut[indR]     = vcible_vec[noind];
        wOut[indR]     = wcible_vec[noind];
        tOut[indR]     = tcible_vec[noind];
        varSAOut[indR] = 0.;

        vxPtr[noind+ideb] = uOut[indR];
        vyPtr[noind+ideb] = vOut[indR];
        vzPtr[noind+ideb] = wOut[indR];

      }
  }
      else //5eq
  {
    // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];

        // For Post (tOut temperature du point image)
        twall = tOut[indR]  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
        densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
        pressPtr[noind+ideb]= press_vec[noind ];

        roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
        uOut[indR]     = ucible_vec[noind];
        vOut[indR]     = vcible_vec[noind];
        wOut[indR]     = wcible_vec[noind];
        tOut[indR]     = tcible_vec[noind];

        vxPtr[noind+ideb] = uOut[indR];
        vyPtr[noind+ideb] = vOut[indR];
        vzPtr[noind+ideb] = wOut[indR];

      }
  }
    }//bctype
  else if (bctype == 10) // loi de paroi Mafzal
    {
#   include "IBC/pointer.h"

      E_Float MafzalMode = 3; // param_real[ MAFZAL_MODE ];

      E_Int err  = 0;
      E_Int skip = 0;
      E_Float tgradU = 0.;
      E_Float ngradU = 0.;
      E_Float unext = 0.;
      //initialisation parametre geometrique et utau
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
  {
    //E_Int indR = rcvPts[noind];
    E_Int indR = rcvPts[noind+ideb];

    roext = roOut[indR]; // densite du point interpole
    text  = tOut[indR];  // pression du point interpole
    pext  = text*roext*cvgam;

    // vitesse du pt ext
    u = uOut[indR];
    v = vOut[indR];
    w = wOut[indR];

    // gradient de pression au point ext
    gradxPext = gradxPressPtr[noind+ideb];
    gradyPext = gradyPressPtr[noind+ideb];
    gradzPext = gradzPressPtr[noind+ideb];

#       include "IBC/commonMafzalLaw_init.h"
  }

      // PREMIERE PASSE MUSKER POUR UTAU_ORI #####################################
#    include "IBC/commonMuskerLaw_Newton.h"
      err  = 0;
      skip = 0;
      E_Float alphaMafMus;

#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
  {
    utau0 = utauOri_vec[noind];
    utauOri_vec[noind] = utau_vec[noind];
    utau_vec[noind] = utau0;

#  include "IBC/mafzal_vec.h"
  }
      // #########################################################################

      // Newton pour utau -> Mafzal
#    include "IBC/commonMafzalLaw_Newton.h"

      //initialisation Newton SA  + vitesse cible
#if NUTILDE_FERRARI == 0
#    include "IBC/commonMafzalLaw_cible.h"
#elif NUTILDE_FERRARI == 1
#    include "IBC/nutilde_Ferrari_Mafzal.h"
#else
#    include "IBC/nutilde_Ferrari_adim_Mafzal.h"
#endif
      if (nvars == 6)
  {
    // Newton pour mut
#if NUTILDE_FERRARI == 0
#       include "IBC/nutildeSA_Newton.h"
#endif
    // mise a jour des variables
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];

        // For Post (tOut temperature du point image en entree, pt corrige en sortie)
        twall = tOut[indR] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
        densPtr[noind+ideb]  = press_vec[noind]/twall*cvgaminv;
        pressPtr[noind+ideb] = press_vec[noind];

        // Mise a jour pt corrige
        roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
        uOut[indR]     = ucible_vec[noind];
        vOut[indR]     = vcible_vec[noind];
        wOut[indR]     = wcible_vec[noind];
        tOut[indR]     = tcible_vec[noind];
        varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind];  //nutilde*signibc

        vxPtr[noind+ideb] = uOut[indR];
        vyPtr[noind+ideb] = vOut[indR];
        vzPtr[noind+ideb] = wOut[indR];
      }
  }
      else //5eq
  {
    // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];

        // For Post (tOut temperature du point image)
        twall = tOut[indR]  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
        densPtr[noind+ideb]  = press_vec[noind]/twall*cvgaminv;
        pressPtr[noind+ideb] = press_vec[noind];

        roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
        uOut[indR]     = ucible_vec[noind];
        vOut[indR]     = vcible_vec[noind];
        wOut[indR]     = wcible_vec[noind];
        tOut[indR]     = tcible_vec[noind];

        vxPtr[noind+ideb] = uOut[indR];
        vyPtr[noind+ideb] = vOut[indR];
        vzPtr[noind+ideb] = wOut[indR];
      }
  }

    }//bctype
  else if (bctype == 11) // TBLE_FULL
    {

#   include "IBC/pointer.h"

      E_Int err  = 0;
      E_Int skip = 0;
      E_Int init = 1;
      E_Int index_int = 0;
      E_Float index_float = 0.;

      E_Float tgradU = 0.;
      E_Float ngradU = 0.;
      E_Float unext = 0.;

#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
  {
    E_Int indR = rcvPts[noind+ideb];

    roext = roOut[indR]; 
    text  = tOut[indR];  
    pext  = text*roext*cvgam;

    u = uOut[indR];
    v = vOut[indR];
    w = wOut[indR];

    // gradient de pression au point ext
    gradxPext = gradxPressPtr[noind+ideb];
    gradyPext = gradyPressPtr[noind+ideb];
    gradzPext = gradzPressPtr[noind+ideb];

    // gradients de vitesse au point ext
    gradxUext = gradxUPtr[noind+ideb];
    gradyUext = gradyUPtr[noind+ideb];
    gradzUext = gradzUPtr[noind+ideb];

    gradxVext = gradxVPtr[noind+ideb];
    gradyVext = gradyVPtr[noind+ideb];
    gradzVext = gradzVPtr[noind+ideb];

    gradxWext = gradxWPtr[noind+ideb];
    gradyWext = gradyWPtr[noind+ideb];
    gradzWext = gradzWPtr[noind+ideb];

#     include "IBC/commonTbleFull_init.h"
  }

      // Newton pour utau
#   include "IBC/commonMuskerLaw_Newton.h"

      if( (alphasbeta_linePtr[ ideb ] != 0.0) || ideb == ifin)   {init=0;}

      if(init)
  {

    for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {

        E_Float* yline1d   = y_linePtr + (noind + ideb)*nbptslinelets;

        a0 = xPC[noind+ideb]-xPW[noind+ideb];
        a1 = yPC[noind+ideb]-yPW[noind+ideb];
        a2 = zPC[noind+ideb]-zPW[noind+ideb];

        E_Float dist = sqrt(a0*a0 + a1*a1 + a2*a2);

        E_Int imin  = 0;
        E_Int imax  = nbptslinelets;
        E_Int isearch    = nbptslinelets/2;

        while ( (imax - imin) > 1  )
    {
      if ( dist <= yline1d[isearch])
        {
          imax = isearch;
          isearch = imin + (imax-imin)/2;
        }
      else
        {
          imin = isearch;
          isearch = imin + (imax-imin)/2;
        }
    }

        alphasbeta_linePtr[noind + ideb] = (dist-yline1d[imin])/(yline1d[imax]-yline1d[imin]);
        index_float = imax;
        index_linePtr[noind + ideb] = index_float;
      }

    for ( E_Int iline = 0 ; iline < nbptslinelets; iline++)
      {
        for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {

      E_Int   indl   = (noind + ideb)*nbptslinelets + iline;
      E_Int   indR   = rcvPts[noind+ideb];

#     include "IBC/commonLaws1.h" 
      yplus = utau_vec[noind]*yplus_vec[noind]/yibc*y_linePtr[indl];

      denoml10 = yplus*yplus-8.15*yplus+86.;
      denoml10 = denoml10*denoml10;
              
      umod = utau_vec[noind]*(5.424*atan((2.*yplus-8.15)/16.7) + log10(pow(yplus+10.6,9.6)/denoml10) - 3.52);
      umod = K_FUNC::E_abs(umod);

      u_linePtr[indl] = umod;
      psi_linePtr[indl] = umod/uext_vec[noind]; // Berger's law
      // psi_linePtr[indl] = pow(umod/uext_vec[noind], 2); // Marc's law
      nutilde_linePtr[indl] = (varSAOut[indR] )*K_FUNC::E_abs(y_linePtr[indl]/y_linePtr[nbptslinelets-1]);

    }

      }
  }

      else{
  for ( E_Int iline = 0 ; iline < nbptslinelets; iline++)
    {
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {

    E_Int   indl   = (noind + ideb)*nbptslinelets + iline;
    E_Int   indR   = rcvPts[noind+ideb];

#     include "IBC/commonLaws1.h" 
    yplus = utau_vec[noind]*yplus_vec[noind]/yibc*y_linePtr[indl];

    denoml10 = yplus*yplus-8.15*yplus+86.;
    denoml10 = denoml10*denoml10;
            
    umod = utau_vec[noind]*(5.424*atan((2.*yplus-8.15)/16.7) + log10(pow(yplus+10.6,9.6)/denoml10) - 3.52);
    umod = K_FUNC::E_abs(umod);

    psi_linePtr[indl] = umod/uext_vec[noind]; // Berger's law
    // psi_linePtr[indl] = pow(umod/uext_vec[noind], 2); // Marc's law
        }

    }
      } // END INIT


      // TBLE FULL + ODE 1D SA for nutilde

      for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {

    E_Int indR = rcvPts[noind+ideb];

    E_Float* yline1d   = y_linePtr        + (noind + ideb)*nbptslinelets;
    E_Float* u1d       = u_linePtr        + (noind + ideb)*nbptslinelets;
    E_Float* nutilde1d = nutilde_linePtr  + (noind + ideb)*nbptslinelets;
    E_Float* psild     = psi_linePtr      + (noind + ideb)*nbptslinelets;
    E_Float* mat       = mat_linePtr      + (noind + ideb)*nbptslinelets;
    E_Float* matm      = matm_linePtr     + (noind + ideb)*nbptslinelets;
    E_Float* matp      = matp_linePtr     + (noind + ideb)*nbptslinelets;

    E_Float nu = mu_vec[noind]/ro_vec[noind];

    utauOri_vec[noind] = utau_vec[noind];

    u1d[0] = 0.0;
    u1d[nbptslinelets-1] = uext_vec[noind];

    psild[0] = 0.0;
    psild[nbptslinelets-1] = 1.;

    npass = 0;
    L2norm = 1.0;
    L2norm0= 1.0;

    nmax = 30;


    for (E_Int j = 0; j < nbptslinelets; ++j)
      {
        ipt_u1dold[j] = 0.0;
      }

    while ( (L2norm >= 1.e-1) && (npass < nmax))
      {

#           include "IBC/tble_1d_full.h"

        utau_vec[noind ] = sqrt(K_FUNC::E_abs(nu*u1d[1])/yline1d[1]);

        spalart_1d_(ithread,yline1d,matm,mat,matp,nutilde1d,u1d,kappa,
        nu,varSAOut[indR],nbptslinelets,kappa);
      }

    if ((npass == nmax) && (ithread==1))
      {
        std::cout << "FULL TBLE doesn't converge ! " << L2norm << std::endl;
      }

        }

      // mise a jour des variables
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {

          E_Int indR = rcvPts[noind+ideb];

          E_Float* u1d       = u_linePtr        + (noind + ideb)*nbptslinelets;
          E_Float* nutilde1d = nutilde_linePtr  + (noind + ideb)*nbptslinelets;

          index_int = index_linePtr[noind + ideb];
          umod = (u1d[index_int] - u1d[index_int-1])*alphasbeta_linePtr[noind + ideb] + u1d[index_int-1];

          yplus            = utau_vec[noind]*yplus_vec[noind];
          yplus_vec[noind] = yplus;

          ucible0 = sign_vec[noind] * umod;
          ucible_vec[noind] += ucible0 * ut_vec[noind]; 
          vcible_vec[noind] += ucible0 * vt_vec[noind];
          wcible_vec[noind] += ucible0 * wt_vec[noind];

          // uext: norme de la composante tangentielle de la vitesse externe
          uext = sqrt(ut_vec[noind]*ut_vec[noind]+vt_vec[noind]*vt_vec[noind]+wt_vec[noind]*wt_vec[noind]);
          uext = std::max(uext, 1.e-12);

          // For Post (tOut temperature du point image en entree, pt corrige en sortie)
          twall = tOut[indR] + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
          densPtr[noind+ideb] = press_vec[noind]/twall*cvgaminv;
          pressPtr[noind+ideb]= press_vec[noind];

          // Mise a jour pt corrige
          roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
          uOut[indR]     = ucible_vec[noind];
          vOut[indR]     = vcible_vec[noind];
          wOut[indR]     = wcible_vec[noind];
          tOut[indR]     = tcible_vec[noind];

          index_int = index_linePtr[noind + ideb];
          varSAOut[indR] = (nutilde1d[index_int] - nutilde1d[index_int-1])*alphasbeta_linePtr[noind + ideb] + nutilde1d[index_int-1]; 

          vxPtr[noind+ideb] = uOut[indR];
          vyPtr[noind+ideb] = vOut[indR];
          vzPtr[noind+ideb] = wOut[indR];

        }
    }
  else if (bctype == 12) // isothermal - prescribed wall temp
    {
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
    E_Int indR = rcvPts[noind+ideb];

    // get values at image point (stored in FlowSolution#Centers)
    u = uOut[indR];     
    v = vOut[indR];     
    w = wOut[indR];     
    t = tOut[indR]; 

    //Calc values at target points
# include "IBC/commonBCType1.h"
    tcible = t*alphasbeta + (1-alphasbeta)*tempExtraPtr[ noind + ideb];

    // set values at target point(stored in FlowSolution#Centers)
    uOut[indR] = ucible; 
    vOut[indR] = vcible; 
    wOut[indR] = wcible; 
    tOut[indR] = tcible; 
    if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

    // set values in tc
    pressPtr[noind + ideb] = roOut[indR]* tOut[indR]*cvgam; 
    densPtr[ noind + ideb] = roOut[indR];                   
    tempPtr[ noind + ideb] = tOut[indR];                    
    vxPtr[noind+ideb] = uOut[indR];                         
    vyPtr[noind+ideb] = vOut[indR];                         
    vzPtr[noind+ideb] = wOut[indR];                         
        }
    }
  else if (bctype == 13) // heat flux - prescribed heat flux
    {
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
    E_Int indR = rcvPts[noind+ideb];

    // get values at image point (stored in FlowSolution#Centers)
    u = uOut[indR];     
    v = vOut[indR];     
    w = wOut[indR];     
    t = tOut[indR]; 

    //Calc values at target points
# include "IBC/commonBCType1.h"
    if (int(copysign(1.0,alpha)) == int(copysign(1.0,beta)))
      {      
        //Target points are inside the fluid
        // second order one sided difference w/ first order one sided different
        // ∂T/∂n=q=(T_T - T_W)/(x_T-x_W) -> T_W=T_T-q*(x_T-x_W)
        // use above in
        //  ∂T/∂n=q=(-T_IP*(x_T-x_W)^2+T_T*(x_IP-x_W)^2-T_W*[(x_IP-x_W)^2-(x_T-x_W)^2])/[(x_T-x_W)(x_I-x_W)(x_I-x_T)]
        E_Float absquared=pow((alpha+beta),2);
        E_Float asquared =pow(alpha,2);
        E_Float tmp_var  =absquared-asquared;
        tcible = (tempExtraPtr[ noind + ideb]*alpha*beta*(alpha+beta)+asquared*t-alpha*tempExtraPtr[ noind + ideb]*tmp_var)/(absquared-tmp_var);
      }
    else
      {
        //Target points are inside the solid
        // second order central difference
        // ∂T/∂n=q=(T_IP - T_T)/(x_IP-x_T) -> T_T=T_IP-q*(x_IP-x_T) 
        tcible = t-(abs(alpha)+abs(beta))*tempExtraPtr[ noind + ideb];
      }
    

    // set values at target point(stored in FlowSolution#Centers)
    uOut[indR] = ucible; 
    vOut[indR] = vcible; 
    wOut[indR] = wcible; 
    tOut[indR] = tcible; 
    if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

    // set values in tc
    pressPtr[noind + ideb] = roOut[indR]* tOut[indR]*cvgam; 
    densPtr[ noind + ideb] = roOut[indR];                   
    tempPtr[ noind + ideb] = tOut[indR];                    
    vxPtr[noind+ideb] = uOut[indR];                         
    vyPtr[noind+ideb] = vOut[indR];                         
    vzPtr[noind+ideb] = wOut[indR];                         
        }
    }
  else if (bctype == 140) // Wire Model - M. Terracol & E. Monaha 2021 - Numerical Wire Mesh Model for the Simulation of Noise Reduction Devices
    {
#   include "IBC/pointer.h"
      E_Int err  = 0;
      E_Int skip = 0;
      
      E_Float ro_Pnt2, u_Pnt2, v_Pnt2, w_Pnt2, t_Pnt2;
      E_Float p_w,p_t,s_w,t_t;
      E_Float ro_up,u_up,v_up,w_up;
      E_Float u_n_w,u_t_w;
      E_Float v_n_w,v_t_w;
      E_Float w_n_w,w_t_w;

      E_Float* roOut_Pnt2    = densPtr +  5*nbRcvPts;// ro          
      E_Float* uOut_Pnt2     = densPtr +  6*nbRcvPts;// u           
      E_Float* vOut_Pnt2     = densPtr +  7*nbRcvPts;// v            
      E_Float* wOut_Pnt2     = densPtr +  8*nbRcvPts;// w           
      E_Float* tOut_Pnt2     = densPtr +  9*nbRcvPts;// temperature 
      E_Float* varSAOut_Pnt2 = densPtr + 10*nbRcvPts;// pseudoviscosity nu_tilde 
        
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
  {
    E_Int indR    = rcvPts[noind+ideb];
        
    // get values at image point (stored in FlowSolution#Centers)
    // These image point values are those for tc
    // 1pnt
    ro            = roOut[indR]; 
    u             = uOut[indR];     
    v             = vOut[indR];     
    w             = wOut[indR];     
    t             = tOut[indR];

    // These image point values are those for tc2
    // Value at A2 for the wire model (second image points)
    // 2pnt
    ro_Pnt2 = roOut_Pnt2[noind+ ideb]; 
    u_Pnt2  = uOut_Pnt2[noind+ ideb];     
    v_Pnt2  = vOut_Pnt2[noind+ ideb];     
    w_Pnt2  = wOut_Pnt2[noind+ ideb];     
    t_Pnt2  = tOut_Pnt2[noind+ ideb];
    
# include "IBC/commonGeom.h"
    s_w = u*n0+v*n1+w*n2;
    s_w = copysign(1,s_w);
    // Flow -->
    //           |
    // s_w=-1    |  s_w=1
    //           |
    //   x   n<--|
    //   ↑
    //  target pnt
    if (s_w<0){
      ro_up=ro;
      u_up=u;
      v_up=v;
      w_up=w;
    }
    else{
      ro_up=ro_Pnt2;
      u_up=u_Pnt2;
      v_up=v_Pnt2;
      w_up=w_Pnt2;
    }

    p_w    = 0.5*R_gas*(ro_Pnt2*t_Pnt2+ro*t);
    p_t    = p_w-s_w*0.25*K_wire*ro_up*(u_up*u_up+v_up*v_up+w_up*w_up);
    tcible = p_t/(ro_up*R_gas);
    twall  = p_w/(ro_up*R_gas);

    roOut[indR]= ro_up;
    tOut[indR] = tcible;

    //save in tc 
    pressPtr[noind + ideb] = roOut[indR]* tOut[indR]*cvgam;	  
    densPtr[ noind + ideb] = roOut[indR];

    // vec is to save the value for the following loop of the target points only
    // save wall values as interpolation is of interpolated & wall are used to determine
    // target values
    ro_vec[noind]=ro_up;	  
    mu_vec[noind]=coefSuth * sqrt(K_FUNC::E_abs(twall)*Tsinv) / (1.+Cs/twall);

    // UP TO HERE CORRECT	  
    
    //n0=-n0;n1=-n1;n2=-n2; // needed to get normals in correct direction
# include "IBC/normalTangentVelocity.h"
    t0=ut/uext;t1=vt/uext;t2=wt/uext;
    
    uext_image  =uext;
    uscaln_image=uscaln;
    
    u=u_up;v=v_up;w=w_up;	
# include "IBC/normalTangentVelocity.h"
    uext_wall  = uext;
    uscaln_wall= uscaln;
    uext_wall  = Delta_V_wire*uext_wall;
    

    //uext      = alphasbeta*uext_image  +(1.-alphasbeta)*uext_wall;
    //uscaln    = alphasbeta*uscaln_image+(1.-alphasbeta)*uscaln_wall;
    
    //ucible = uscaln*n0+uext*t0;
    //vcible = uscaln*n1+uext*t1;
    //wcible = uscaln*n2+uext*t2;

    ucible      = alphasbeta*uOut[indR]  +(1.-alphasbeta)*(uscaln_wall*n0+uext_wall*t0);
    vcible      = alphasbeta*vOut[indR]  +(1.-alphasbeta)*(uscaln_wall*n1+uext_wall*t1);
    wcible      = alphasbeta*wOut[indR]  +(1.-alphasbeta)*(uscaln_wall*n2+uext_wall*t2);

    //if (yPC[noind+ideb]>0 && yPC[noind+ideb]<0.0007){
    //  printf("xPC[noind+ideb] uscaln_wall uscaln_target:: %g %g %g %g %g\n",xPC[noind+ideb],uscaln_wall,uscaln_target,ucible,vcible);
    //}	  
    
    uOut[indR] = ucible; 
    vOut[indR] = vcible; 
    wOut[indR] = wcible; 	  

    //save in tc 
    vxPtr[noind+ideb] = uOut[indR];                         
    vyPtr[noind+ideb] = vOut[indR];                         
    vzPtr[noind+ideb] = wOut[indR];

    //vec is to save the value for the following loop of the target points only	  
    ut_vec[noind]=uext_wall*t0;
    vt_vec[noind]=uext_wall*t1;	
    wt_vec[noind]=uext_wall*t2;

    ucible_vec[noind]=uscaln_wall*n0;
    vcible_vec[noind]=uscaln_wall*n1;	
    wcible_vec[noind]=uscaln_wall*n2;
  }
    
      if (nvars==6){
  E_Int count = 0;
  for (E_Int noind = 0; noind < ifin-ideb; noind++)
  {
    E_Int indR    = rcvPts[noind+ideb];

    //|u_w|=sqrt(un,w²+ut,w²)
    E_Float norm_wall_vel = sqrt(ucible_vec[noind]*ucible_vec[noind]+vcible_vec[noind]*vcible_vec[noind]+wcible_vec[noind]*wcible_vec[noind]+
               ut_vec[noind]    *ut_vec[noind]    +vt_vec[noind]    *vt_vec[noind]    +wt_vec[noind]    *wt_vec[noind]);
    
    E_Float nu_t_local    = Ct_WM*Diam_wire*norm_wall_vel; // mu_t = rho*C_t*d*|u_w| & nu_t=mu_t/rho

    nutcible_vec[noind] = nu_t_local;
    nutilde             = K_FUNC::E_abs( nutcible_vec[noind] );
    aa_vec[noind]       = nutilde; //to start first iteration
    ut_vec[noind]       = nutilde; // Sauvegarde du nutilde, au cas ou newton non convergent.
                                   // save  ṽ(O)=v_t (v=nu)  	  
# include "IBC/fnutilde_vec.h"
    
  }
  // Newton pour ṽ
#       include "IBC/nutildeSA_Newton.h"
  
  for (E_Int noind = 0; noind < ifin-ideb; noind++)
  {
    E_Int indR    = rcvPts[noind+ideb];
# include "IBC/commonGeom.h"	  
    varSAOut_Pnt2[noind+ ideb]=nutcible_vec[noind];
    varSAOut[indR]=alphasbeta*varSAOut[indR]  +(1.-alphasbeta)*(aa_vec[noind]);
  }
      }
    }
  else if (bctype == 141) // Wire Model prt2 - just interpolation of image points & placed at control points
    {;}      
  else
    {
      printf("Warning !!! setIBCTransfersCommonVar2: bcType " SF_D_ " not implemented.\n", bctype);
      return 0;
    }

 WireMeshSkip:
  return 1;
}

//=============================================================================
//Retourne -2 : incoherence entre meshtype et le type d interpolation
//         -1 : type invalide
//          1 : ok
// Entree/Sortie :  (ro,u,v,w,p) ( + ronutildeSA )
//=============================================================================
E_Int K_CONNECTOR::setIBCTransfersCommonVar3(
               E_Int bctype,
               E_Int* rcvPts, E_Int& nbRcvPts, E_Int& ideb, E_Int& ifin,  E_Int& ithread,
               E_Float* xPC, E_Float* yPC, E_Float* zPC,
               E_Float* xPW, E_Float* yPW, E_Float* zPW,
               E_Float* xPI, E_Float* yPI, E_Float* zPI,
               E_Float* densPtr, E_Float* pressPtr,
               E_Float* vxPtr, E_Float* vyPtr, E_Float* vzPtr,
               E_Float* utauPtr, E_Float* yplusPtr, E_Float* kcurvPtr,
               E_Float* d1, E_Float* d2, E_Float* d3, E_Float* d4, E_Float* d5,
               E_Float* tmp, E_Int& size,
               E_Float gamma, E_Float cv, E_Float muS, E_Float Cs, E_Float Ts, E_Float Pr,
               vector<E_Float*>& vectOfDnrFields, vector<E_Float*>& vectOfRcvFields)
{
  /* lois de paroi */
  E_Float roext, uext, pext, text, muext, yext, yplus, yibc;
  E_Float uscaln, un, vn, wn, ut, vt, wt, utauv, utau0, umod;
  E_Float bb, fp, tp;
  E_Float expy, denoml10,ax,l1,l2, l3;
  E_Float ucible0, ucible, vcible, wcible, signibc, twall, rowall, muwall;
  //Lois de paroi : criteres d arret pour estimer le frottement par Newton
  E_Float newtoneps = 1.e-7; // critere d arret pour u+
  //E_Float newtonepsnutilde = 1.e-10; // critere d arret pour nutilde
  E_Float newtonepsprime = 1.e-12;// critere d arret pour la derivee
  E_Float cvgaminv = 1./(cv*(gamma-1.));
  E_Float coefSuth = muS * (1.+Cs/Ts);
  E_Float Tsinv = 1./Ts;
  E_Float kappa = 0.4; // Constante de Von Karman
  E_Float kappainv = 1./kappa;
  E_Float cc = 5.2;//pour la loi log
  E_Float one_third = 1./3.;

  /* fin parametres loi de parois */

  E_Int nvars = vectOfDnrFields.size();

  //E_Float a[3], b[3], n[3];
  E_Float a0,a1,a2,b0,b1,b2,n0,n1,n2;
  E_Float normb, u,v,w;
  E_Float vnc, alpha, beta, alphasbeta;
  E_Float* roOut = vectOfRcvFields[0];// density
  E_Float* uOut = vectOfRcvFields[1];// ux
  E_Float* vOut = vectOfRcvFields[2];// uy
  E_Float* wOut = vectOfRcvFields[3];// uz
  E_Float* pOut = vectOfRcvFields[4];// pressure
  E_Float* varSAOut = NULL;
  if ( nvars == 6 ) varSAOut = vectOfRcvFields[5];//ronutildeSA
  // if ( (bctype==2 || (bctype==3)) && nvars < 6)
  // {
  //   printf("Warning: setIBCTransfersCommonVar3: number of variables (<6) inconsistent with bctype.\n");
  //   return 0;
  // }
  if (bctype == 0) // wallslip
    {
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
    E_Int indR = rcvPts[noind+ideb];

    // vitesse
    u = uOut[indR];
    v = vOut[indR];
    w = wOut[indR];
# include "IBC/commonBCType0.h"
    uOut[indR] = ucible;
    vOut[indR] = vcible;
    wOut[indR] = wcible;
    if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

    pressPtr[noind +ideb] = pOut[indR];
    densPtr[ noind +ideb] = roOut[indR];

    vxPtr[noind+ideb] = uOut[indR];
    vyPtr[noind+ideb] = vOut[indR];
    vzPtr[noind+ideb] = wOut[indR];

        }
    }
  else if (bctype == 1) // adherence
    {
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
    E_Int indR = rcvPts[noind+ideb];

    // vitesse
    u = uOut[indR];
    v = vOut[indR];
    w = wOut[indR];
# include "IBC/commonBCType1.h"
    uOut[indR] = ucible;
    vOut[indR] = vcible;
    wOut[indR] = wcible;
    if (nvars == 6) varSAOut[indR] = varSAOut[indR]*alphasbeta;

    pressPtr[noind +ideb] = pOut[indR];
    densPtr[ noind +ideb] = roOut[indR];

    vxPtr[noind+ideb] = uOut[indR];
    vyPtr[noind+ideb] = vOut[indR];
    vzPtr[noind+ideb] = wOut[indR];
        }
    }
  else if (bctype == 2) // loi de paroi en log
    {
#   include "IBC/pointer.h"

      E_Int err  = 0;
      E_Int skip = 0;
      //initilisation parametre geometrique et utau
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
  {
    E_Int indR = rcvPts[noind+ideb];

    roext= roOut[indR]; // densite du point interpole
    pext = pOut[indR];  // pression du point interpole
    text = pext / roext * cvgaminv;

    // vitesse du pt ext
    u = uOut[indR];
    v = vOut[indR];
    w = wOut[indR];

#       include "IBC/commonLogLaw_init.h"
    // out= utau  et err
  }

      // Newton pour utau
#   include "IBC/commonLogLaw_Newton.h"

      //initialisation Newton SA  + vitesse cible
#if NUTILDE_FERRARI == 0
#   include "IBC/commonLogLaw_cible.h"
#else
#   include "IBC/nutilde_Ferrari.h"
#endif
      if (nvars == 6)
  {
    // Newton pour mut
#if NUTILDE_FERRARI == 0
#       include "IBC/nutildeSA_Newton.h"
#endif
    // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];

        // For Post (temperature du point image text, densite et pression du point de paroi: densPtr,pressPtr)
        text = press_vec[noind]/ ro_vec[noind]* cvgaminv;
        twall = text  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
        densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
        pressPtr[noind+ideb]= press_vec[noind ];

        roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
        uOut[indR]     = ucible_vec[noind]; vOut[indR] = vcible_vec[noind]; wOut[indR] = wcible_vec[noind];
        varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind];                                         //nutilde*signibc

        vxPtr[noind+ideb] = uOut[indR];
        vyPtr[noind+ideb] = vOut[indR];
        vzPtr[noind+ideb] = wOut[indR];
      }
  }
      else //5eq
  {
    // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];

        // For Post (temperature du point image text, densite et pression du point de paroi: densPtr,pressPtr)
        text = press_vec[noind]/ ro_vec[noind]* cvgaminv;
        twall = text  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
        densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
        pressPtr[noind+ideb]= press_vec[noind ];

        roOut[indR] = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
        uOut[indR] = ucible_vec[noind]; vOut[indR] = vcible_vec[noind]; wOut[indR] = wcible_vec[noind];

        vxPtr[noind+ideb] = uOut[indR];
        vyPtr[noind+ideb] = vOut[indR];
        vzPtr[noind+ideb] = wOut[indR];
      }
  }
    }
  else if (bctype == 3)// loi de paroi de Musker
    {
#   include "IBC/pointer.h"

      E_Int err  = 0;
      E_Int skip = 0;
      //initilisation parametre geometrique et utau
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
  {
    E_Int indR = rcvPts[noind+ideb];

    roext= roOut[indR]; // densite du point interpole
    pext = pOut[indR];  // pression du point interpole
    text = pext / roext * cvgaminv;

    // vitesse du pt ext
    u = uOut[indR];
    v = vOut[indR];
    w = wOut[indR];

#       include "IBC/commonMuskerLaw_init.h"
    // out= utau  et err
  }

      // Newton pour utau
#   include "IBC/commonMuskerLaw_Newton.h"

      //initialisation Newton SA  + vitesse cible
#if NUTILDE_FERRARI == 0
#   include "IBC/commonMuskerLaw_cible.h"
#else
#   include "IBC/nutilde_Ferrari.h"
#endif
      if (nvars == 6)
  {
    // Newton pour mut
#if NUTILDE_FERRARI == 0
#       include "IBC/nutildeSA_Newton.h"
#endif
    // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];

        uOut[indR]     = ucible_vec[noind]; vOut[indR] = vcible_vec[noind]; wOut[indR] = wcible_vec[noind];
        varSAOut[indR] = aa_vec[noind]*sign_vec[noind]*uext_vec[noind];                                         //nutilde*signibc
      }
  }
      else //5eq
  {
    // mise a jour des variable
#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];

        // For Post (temperature du point image text, densite et pression du point de paroi: densPtr,pressPtr)
        text = press_vec[noind]/ ro_vec[noind]* cvgaminv;
        twall = text  + 0.5*pow(Pr,one_third)/(cv*gamma)*(uext_vec[noind]*uext_vec[noind]); // Crocco-Busemann
        densPtr[noind+ideb] = press_vec[noind ]/twall*cvgaminv;
        pressPtr[noind+ideb]= press_vec[noind ];

        roOut[indR]    = press_vec[noind ]/tcible_vec[noind]*cvgaminv;
        uOut[indR] = ucible_vec[noind]; vOut[indR] = vcible_vec[noind]; wOut[indR] = wcible_vec[noind];

        vxPtr[noind+ideb] = uOut[indR];
        vyPtr[noind+ideb] = vOut[indR];
        vzPtr[noind+ideb] = wOut[indR];
      }
  }
    }
  else if (bctype == 4) // outpres
    {
#ifdef _OPENMP4
#pragma omp simd
#endif
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
        {
          E_Int indR = rcvPts[noind+ideb];
          densPtr[noind+ideb] = roOut[indR];
          vxPtr[noind+ideb] = uOut[indR];
          vyPtr[noind+ideb] = vOut[indR];
          vzPtr[noind+ideb] = wOut[indR];
        }
    }
  else
    {
      printf("Warning !!! setIBCTransfersCommonVar3: bcType " SF_D_ " not implemented.\n", bctype);
      return 0;
    }
  return 1;
}
//=============================================================================
/* Effectue les transfers IBC */
//=============================================================================
//Stephanie: fonction absente du pytree?
PyObject* K_CONNECTOR::setIBCTransfers(PyObject* self, PyObject* args)
{
  PyObject *arrayR, *arrayD;
  PyObject *pyVariables;
  PyObject *pyIndRcv, *pyIndDonor;
  PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs;
  PyObject *pyArrayXPC, *pyArrayXPI, *pyArrayXPW;
  PyObject *pyArrayYPC, *pyArrayYPI, *pyArrayYPW;
  PyObject *pyArrayZPC, *pyArrayZPI, *pyArrayZPW;
  PyObject *pyArrayDens;
  E_Int bctype;
  E_Int vartype;
  E_Float gamma, cv, muS, Cs, Ts;

  if (!PYPARSETUPLE_(args,
                    OOOO_ OOOO_ OOOO_ OOOO_ O_ II_ RRRR_ R_, 
                    &arrayR, &arrayD,  &pyVariables,
                    &pyIndRcv, &pyIndDonor, &pyArrayTypes, &pyArrayCoefs,
                    &pyArrayXPC, &pyArrayYPC, &pyArrayZPC,
                    &pyArrayXPW, &pyArrayYPW, &pyArrayZPW,
                    &pyArrayXPI, &pyArrayYPI, &pyArrayZPI,
                    &pyArrayDens, 
                    &bctype, &vartype, &gamma, &cv, &muS, &Cs, &Ts))
    {
      return NULL;
    }


  E_Int bcType = E_Int(bctype); // 0 : wallslip; 1: noslip; 2: log law of wall; 3: Musker law of wall
  /* varType :
     1  : conservatives,
     11 : conservatives + ronutildeSA
     2  : (ro,u,v,w,t)
     21 : (ro,u,v,w,t) + ronutildeSA
     3  : (ro,u,v,w,p)
     31 : (ro,u,v,w,p) + ronutildeSA */
  E_Int varType = E_Int(vartype);

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
          "setIBCTransfers: 1st arg is not a valid array.");
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
  E_Int* ptrcnd = NULL;
  E_Int cndSize = 0; E_Int cnNfldD = 0;
  if (resd != 2 && resd != 1)
    {
      PyErr_SetString(PyExc_TypeError,
          "setIBCTransfers: 2nd arg is not a valid array.");
      RELEASESHAREDB(resr, arrayR, fr, cnr);
      return NULL;
    }
  E_Int meshtype = resd;// 1 : structure, 2 non structure
  if (resd == 2)
    {
      if (strcmp(eltTypeD,"TETRA") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
        "setIBCTransfers: donor array must be a TETRA if unstructured.");
    RELEASESHAREDB(resr, arrayR, fr, cnr);
    RELEASESHAREDB(resd, arrayD, fd, cnd);
    return NULL;
  }
      ptrcnd  = cnd->begin();
      cndSize = cnd->getSize();
      cnNfldD = cnd->getNfld();
    }

  E_Int nvars;
  if ( varType == 1 || varType == 2 || varType == 3 )
    nvars = 5;
  else if ( varType == 11 || varType == 21 || varType == 31 )
    nvars = 6;
  else
    {
      RELEASESHAREDB(resr, arrayR, fr, cnr);
      RELEASESHAREDB(resd, arrayD, fd, cnd);
      PyErr_SetString(PyExc_TypeError,
          "setIBCTransfers: varType value is not valid.");
      return NULL;
    }

# include "extract_interpD.h"
  /*--------------------------------------*/
  /* Extraction des indices des receveurs */
  /*--------------------------------------*/
  FldArrayI* rcvPtsI;
  E_Int res_rcv = K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI);
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
      PyErr_SetString(PyExc_TypeError,"setInterpTransfersD: 4th to 6th arg must be a numpy of integers. 7th arg a numpy floats ");
      return NULL;
    }

# include "IBC/extract_IBC.h"

  if (okc1*okc2*okc3 == 0 )
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
  if (okD == 0 )
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
      PyErr_SetString(PyExc_TypeError, 
          "setIBCTransfersD: Post array are invalid.");
      return NULL;
    }


  // Transferts
  // Types valides: 2, 3, 4, 5
  PyObject* tpl;
  if (resr == 1)
    {
      tpl = K_ARRAY::buildArray(fr->getNfld(), varStringR, imr, jmr, kmr);
    }
  else // unstructured
    {
      E_Int crsize = cnr->getSize()*cnr->getNfld();
      tpl = K_ARRAY::buildArray(fr->getNfld(), varStringR,
        fr->getSize(), cnr->getSize(),
        -1, eltTypeR, false, crsize);
      E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
      K_KCORE::memcpy__(cnnp, cnr->begin(), cnr->getSize()*cnr->getNfld());
    }
  E_Float* ptrFieldOut = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fieldROut(fr->getSize(), fr->getNfld(), ptrFieldOut, true);
  fieldROut = *fr;

  //vector<E_Float*> vectOfDnrFields;
  //vector<E_Float*> vectOfRcvFields;
  E_Float** RcvFields = new E_Float*[ nvars*10];
  E_Float** DnrFields = new E_Float*[ nvars*10];
  E_Float** vectOfRcvFields = RcvFields;
  E_Float** vectOfDnrFields = DnrFields;

  E_Int posvr, posvd;
  // Extrait les positions des variables a transferer
  E_Int nfoundvar = 0;
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
          vectOfDnrFields[nfoundvar]= fd->begin(posvd+1);
          vectOfRcvFields[nfoundvar]= fieldROut.begin(posvr+1);
          //vectOfDnrFields.push_back(fd->begin(posvd+1));
          //vectOfRcvFields.push_back(fieldROut.begin(posvr+1));
          nfoundvar += 1;
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
          vectOfDnrFields[nfoundvar]= fd->begin(posvd+1);
          vectOfRcvFields[nfoundvar]= fieldROut.begin(posvr+1);
          //vectOfDnrFields.push_back(fd->begin(posvd+1));
          //vectOfRcvFields.push_back(fieldROut.begin(posvr+1));
          nfoundvar += 1;
        }
    }
#endif
        else
    PyErr_Warn(PyExc_Warning, "setIBCTransfers: variable must be a string. Skipped.");
      }
  }// nvariables > 0
    }
  if (nfoundvar != nvars)
    {
      RELEASESHAREDB(resr, arrayR, fr, cnr);
      RELEASESHAREDB(resd, arrayD, fd, cnd);
      BLOCKRELEASEMEM;
      BLOCKRELEASEMEM2;
      PyErr_SetString(PyExc_TypeError,
          "setIBCTransfers: number of variables is inconsistent with varType.");
      return NULL;
    }

  E_Float param_real[30]; 
  param_real[ GAMMA] = gamma;
  param_real[ CVINF] = cv;
  param_real[ XMUL0] = muS;
  param_real[ CS] = Cs;
  param_real[ TEMP0] = Ts;
  param_real[ PRANDT] = 0.71;

  ////
  ////
  //  Interpolation parallele
  ////
  ////
# include "commonInterpTransfers_indirect.h"


  E_Int threadmax_sdm  = __NUMTHREADS__;

  E_Int size = (nbRcvPts/threadmax_sdm)+1; // on prend du gras pour gerer le residus
  E_Int    r =  size % 8;
  size = size + r;                  // on rajoute du bas pour alignememnt 64bits

  FldArrayF  tmp(size*17*threadmax_sdm);
  E_Float* ipt_tmp=  tmp.begin();

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

    if (varType == 2 || varType == 21) 
      setIBCTransfersCommonVar2(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread,
        xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, 
        density, 
        ipt_tmp, size, nvars,
        param_real,
        vectOfDnrFields, vectOfRcvFields);
    else {printf(" setIBCTransfers: only valid for vartype=2 or 21 \n");}
 
  } // Fin zone // omp

  delete [] RcvFields;  delete [] DnrFields;
  // sortie
  RELEASESHAREDB(resr, arrayR, fr, cnr);
  RELEASESHAREDB(resd, arrayD, fd, cnd);
  BLOCKRELEASEMEM;
  BLOCKRELEASEMEM2;
  return tpl;
}
//=============================================================================
/* Effectue les transfers IBC en in-place - 
   !!!Warning !!! densPtr doit etre compacte et contient les autres variables IBC
   selon le bcType */
//=============================================================================
PyObject* K_CONNECTOR::_setIBCTransfers(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  PyObject *zoneR, *zoneD;
  PyObject *pyVariables; 
  PyObject *pyIndRcv, *pyIndDonor;
  PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs;
  PyObject *pyArrayXPC, *pyArrayXPI, *pyArrayXPW;
  PyObject *pyArrayYPC, *pyArrayYPI, *pyArrayYPW;
  PyObject *pyArrayZPC, *pyArrayZPI, *pyArrayZPW;
  PyObject *pyArrayDens;
  E_Int bctype, loc, vartype, compact;
  E_Float gamma, cv, muS, Cs, Ts;

  if (!PYPARSETUPLE_(args,
                    OOOO_ OOOO_ OOOO_ OOOO_ O_ IIII_ RRRR_ R_ SSS_,
                    &zoneR, &zoneD, &pyVariables,
                    &pyIndRcv  , &pyIndDonor, &pyArrayTypes, &pyArrayCoefs, 
                    &pyArrayXPC, &pyArrayYPC, &pyArrayZPC,
                    &pyArrayXPW, &pyArrayYPW, &pyArrayZPW,
                    &pyArrayXPI, &pyArrayYPI, &pyArrayZPI, 
                    &pyArrayDens, 
                    &bctype, &loc, &vartype, &compact, &gamma, &cv, &muS, &Cs, &Ts,
                    &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters))
    {
      return NULL;
    }

  vector<PyArrayObject*> hook;
  
  // Get Prandtl (FastS)
  E_Float Pr = 0.71;
  PyObject* own   = K_PYTREE::getNodeFromName1(zoneR, ".Solver#ownData");
  if (own != NULL)
    {
      PyObject* paramreal0 = K_PYTREE::getNodeFromName1(own, "Parameter_real");
      if (paramreal0 != NULL)
  {
    E_Float* paramreal0val = K_PYTREE::getValueAF(paramreal0, hook);
    Pr = paramreal0val[10];
  }
    }
  E_Int bcType = E_Int(bctype); // 0 : glissement; 1: adherence; 2: loi log; 3: loi de Musker

  /* varType : 
     1  : conservatives, 
     11 : conservatives + ronutildeSA 
     2  : (ro,u,v,w,t)
     21 : (ro,u,v,w,t) + ronutildeSA 
     3  : (ro,u,v,w,p)     
     31 : (ro,u,v,w,p) + ronutildeSA   */
  E_Int nvars;
  E_Int varType = E_Int(vartype);
  if (varType ==  1 || varType ==  2 || varType ==  3) 
    nvars = 5;
  else if (varType == 11 || varType == 21 || varType == 31)
    nvars = 6;
  else 
    {
      PyErr_SetString(PyExc_TypeError, 
          "_setIBCTransfers: varType value is not valid.");
      return NULL;
    }

  // recupere les champs du donneur (nodes)
  E_Int imdjmd, imd, jmd, kmd, cnNfldD, ndimdxR, ndimdxD, meshtype;;
  E_Float* iptroD; E_Float* iptroR;

# include "extract_interpD.h"

  /*--------------------------------------*/
  /* Extraction des indices des receveurs */
  /*--------------------------------------*/
  FldArrayI* rcvPtsI;
  K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI);
  E_Int* rcvPts  = rcvPtsI->begin();
  nbRcvPts       = rcvPtsI->getSize();  

# include "IBC/extract_IBC.h"

  vector<E_Float*> fieldsR;vector<E_Float*> fieldsD;

  E_Float** RcvFields = new E_Float*[ nvars];
  E_Float** DnrFields = new E_Float*[ nvars];
  E_Float** vectOfRcvFields = RcvFields;
  E_Float** vectOfDnrFields = DnrFields;
  E_Int* ptrcnd=NULL;
  char* eltTypeR; char* eltTypeD;
  //codage general (lent ;-) )
  if (compact == 0)
    {// recupere les champs du donneur (nodes)
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
      // recupere les champs du receveur 
      E_Int imr, jmr, kmr, cnSizeR, cnNfldR;
      char* varStringR; vector<E_Int> locsR;
      vector<E_Int*> cnr;
      K_PYTREE::getFromZone(zoneR, 0, loc, varStringR,
                            fieldsR, locsR, imr, jmr, kmr,
                            cnr, cnSizeR, cnNfldR, eltTypeR, hook,
                            GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);

      // -- no check (perfo) --
      // Transferts 
      // Types valides: 2, 3, 4, 5
      E_Int posvr, posvd;

      // Extrait les positions des variables a transferer
      E_Int nfoundvar = 0;
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
        vectOfRcvFields[nfoundvar]= fieldsR[posvr];
        vectOfDnrFields[nfoundvar]= fieldsD[posvd];
        //vectOfDnrFields.push_back(fieldsD[posvd]);
        //vectOfRcvFields.push_back(fieldsR[posvr]);
        nfoundvar += 1;
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
        vectOfRcvFields[nfoundvar]= fieldsR[posvr];
        vectOfDnrFields[nfoundvar]= fieldsD[posvd];
        //vectOfDnrFields.push_back(fieldsD[posvd]);
        //vectOfRcvFields.push_back(fieldsR[posvr]);
        nfoundvar += 1;
      }
        }
#endif
      else
        {
          PyErr_Warn(PyExc_Warning, "_setIBCTransfers: variable must be a string. Skipped.");
        }
            
    }
      }
  }
      delete [] varStringR; delete [] varStringD; delete [] eltTypeR; delete [] eltTypeD;

      if (nfoundvar != nvars)
  {
    RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
    BLOCKRELEASEMEM;
    BLOCKRELEASEMEM2;
    PyErr_SetString(PyExc_TypeError, 
        "_setIBCTransfers: number of variables is inconsistent with varType.");
    return NULL;
  }
    }  

  // les variables a transferer sont compactees: on recupere uniquement la premiere et la taille 
  else
    {
# include "getfromzonecompact.h"
      for (E_Int eq = 0; eq < nvars; eq++)
  {
    vectOfRcvFields[eq]= iptroR + eq*ndimdxR;
    vectOfDnrFields[eq]= iptroD + eq*ndimdxD;
  }
    }

  ////
  ////
  ////
  //  Interpolation parallele
  ////  
  ////  
  ////  
# include "commonInterpTransfers_indirect.h"

  E_Int threadmax_sdm  = __NUMTHREADS__;

  E_Int size = (nbRcvPts/threadmax_sdm)+1; // on prend du gras pour gerer le residus
  E_Int    r =  size % 8;
  if (r != 0) size  = size + 8 - r;        // on rajoute du bas pour alignememnt 64bits
  //printf("size %d %d \n", size, r);
  if (bctype <=1 ) size = 0;               // tableau inutile

  FldArrayF  tmp(size*17*threadmax_sdm);
  E_Float* ipt_tmp = tmp.begin();

  E_Float param_real[30]; 
  param_real[ GAMMA] = gamma;
  param_real[ CVINF] = cv;
  param_real[ XMUL0] = muS;
  param_real[ CS] = Cs;
  param_real[ TEMP0] = Ts;
  param_real[ PRANDT] = Pr;

#    pragma omp parallel default(shared)
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

    if (varType == 2 || varType == 21){
      setIBCTransfersCommonVar2(bcType, rcvPts, nbRcvPts, ideb, ifin, ithread,
        xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI, 
        density, 
        ipt_tmp, size, nvars,
        param_real,
        vectOfDnrFields, vectOfRcvFields); 

    }
    else 
      {printf("Warning: _setIBCTransfers only valid for vartype=2 or 21 \n");}
  } // Fin zone // omp

  delete [] RcvFields;  delete [] DnrFields;
  // sortie
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
  BLOCKRELEASEMEM;
  BLOCKRELEASEMEM2;
  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
// 
//=============================================================================
PyObject* K_CONNECTOR::_setIBCTransfersForPressureGradientsOrder1(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  PyObject *zoneR, *zoneD;
  PyObject *pyVariables;
  PyObject *pyIndRcv, *pyIndDonor;
  PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs;
  PyObject *pyArrayPressure;
  PyObject *pyArrayGradxP , *pyArrayGradyP , *pyArrayGradzP;
  E_Int loc;

  if (!PYPARSETUPLE_(args,
                    OOOO_ OOOO_ OOO_ I_ SSS_,
                    &zoneR, &zoneD, &pyVariables,
                    &pyIndRcv  , &pyIndDonor, &pyArrayTypes, &pyArrayCoefs,
                    &pyArrayPressure,
                    &pyArrayGradxP , &pyArrayGradyP , &pyArrayGradzP,
                    &loc,
                    &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters))
    {
      return NULL;
    }

  vector<PyArrayObject*> hook;

  E_Int nvars = 4;

  E_Int imdjmd, imd, jmd, kmd, cnNfldD, ndimdxR, ndimdxD, meshtype;
  E_Float* iptroD; E_Float* iptroR;

  # include "extract_interpD.h"
  FldArrayI* rcvPtsI;
  K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI);
  E_Int* rcvPts  = rcvPtsI->begin();
  nbRcvPts       = rcvPtsI->getSize();

  FldArrayF* pressF;
  E_Int okP = K_NUMPY::getFromNumpyArray(pyArrayPressure, pressF);
  E_Float* pressure = pressF->begin();

  FldArrayF* gradxPressF; FldArrayF* gradyPressF; FldArrayF* gradzPressF;
  E_Int okGxP = K_NUMPY::getFromNumpyArray(pyArrayGradxP, gradxPressF);
  E_Int okGyP = K_NUMPY::getFromNumpyArray(pyArrayGradyP, gradyPressF);
  E_Int okGzP = K_NUMPY::getFromNumpyArray(pyArrayGradzP, gradzPressF);
  E_Float* gradxP = gradxPressF->begin();
  E_Float* gradyP = gradyPressF->begin();
  E_Float* gradzP = gradzPressF->begin();

  vector<E_Float*> fieldsR;vector<E_Float*> fieldsD;
  vector<E_Float*> vectOfDnrFields(nvars); vector<E_Float*> vectOfRcvFields(nvars);
  E_Int* ptrcnd;
  char* eltTypeR; char* eltTypeD;

  // Data are not compacted beforehand
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

  meshtype = resd; // 1: structured, 2: unstructured
  E_Int imr, jmr, kmr, cnSizeR, cnNfldR;
  char* varStringR; vector<E_Int> locsR;
  vector<E_Int*> cnr;
  K_PYTREE::getFromZone(zoneR, 0, loc, varStringR,
                        fieldsR, locsR, imr, jmr, kmr,
                        cnr, cnSizeR, cnNfldR, eltTypeR, hook,
                        GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);

  E_Int posvr, posvd;

  E_Int nfoundvar = 0;
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
          if (posvr != -1)
          {
            vectOfRcvFields[nfoundvar]= fieldsR[posvr];
            vectOfDnrFields[nfoundvar]= fieldsD[posvd];
            nfoundvar += 1;
          }
        }
        #if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(tpl0))
        {
          const char* varname = PyUnicode_AsUTF8(tpl0);
          posvd = K_ARRAY::isNamePresent(varname, varStringD);
          posvr = K_ARRAY::isNamePresent(varname, varStringR);
          if (posvr != -1)
          {
            vectOfRcvFields[nfoundvar]= fieldsR[posvr];
            vectOfDnrFields[nfoundvar]= fieldsD[posvd];
            nfoundvar += 1;
          }
        }
        #endif
        else
        {
          PyErr_Warn(PyExc_Warning, "_setIBCTransfers: variable must be a string. Skipped.");
        }
      }
    }
  }

  delete [] varStringR; delete [] varStringD; delete [] eltTypeR; delete [] eltTypeD;

  # include "commonInterpTransfers_indirect.h"

  # pragma omp parallel default(shared)
  {
    E_Int ideb, ifin;
    #ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num()+1;
    E_Int  Nbre_thread_actif = omp_get_num_threads(); 
    #else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
    #endif
    E_Int chunk = nbRcvPts/Nbre_thread_actif;
    E_Int r = nbRcvPts - chunk*Nbre_thread_actif;
    if (ithread <= r)
    { 
      ideb = (ithread-1)*(chunk+1); ifin = ideb + (chunk+1); 
    }
    else 
    { 
      ideb = (chunk+1)*r+(ithread-r-1)*chunk; ifin = ideb + chunk; 
    }

    // 0 : Pressure
    // 1 / 2 / 3 : gradxPressure  / gradyPressure  / gradzPressure
    // 4 / 5 / 6 : gradxxPressure / gradyxPressure / gradzxPressure
    // 7 / 8 / 9 : gradxyPressure / gradyyPressure / gradyyPressure
    // 10/11 /12 : gradxzPressure / gradyzPressure / gradzzPressure

    #ifdef _OPENMP4
    #pragma omp simd
    #endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      E_Int indR = rcvPts[noind+ideb];
      
      pressure[noind+ideb] = vectOfRcvFields[0][indR];

      gradxP[noind+ideb] = vectOfRcvFields[1][indR];
      gradyP[noind+ideb] = vectOfRcvFields[2][indR];
      gradzP[noind+ideb] = vectOfRcvFields[3][indR];
    }
  }

  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);

  RELEASESHAREDN(pyIndRcv, rcvPtsI);
  RELEASESHAREDN(pyIndDonor, donorPtsI);
  RELEASESHAREDN(pyArrayTypes, typesI);
  RELEASESHAREDN(pyArrayCoefs, donorCoefsF);

  RELEASESHAREDN(pyArrayPressure, pressF);
  RELEASESHAREDN(pyArrayGradxP, gradxPressF);
  RELEASESHAREDN(pyArrayGradyP, gradyPressF);
  RELEASESHAREDN(pyArrayGradzP, gradzPressF);

  Py_INCREF(Py_None);
  return Py_None;
}
//=============================================================================
// 
//=============================================================================
PyObject* K_CONNECTOR::_setIBCTransfersForPressureGradientsOrder2(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  PyObject *zoneR, *zoneD;
  PyObject *pyVariables;
  PyObject *pyIndRcv, *pyIndDonor;
  PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs;
  PyObject *pyArrayPressure;
  PyObject *pyArrayGradxP , *pyArrayGradyP , *pyArrayGradzP;
  PyObject *pyArrayGradxxP, *pyArrayGradxyP, *pyArrayGradxzP;
  PyObject *pyArrayGradyxP, *pyArrayGradyyP, *pyArrayGradyzP;
  PyObject *pyArrayGradzxP, *pyArrayGradzyP, *pyArrayGradzzP;
  E_Int loc;

  if (!PYPARSETUPLE_(args,
                    OOOO_ OOOO_ OOOO_ OOOO_ OOOO_ I_ SSS_,
                    &zoneR, &zoneD, &pyVariables,
                    &pyIndRcv  , &pyIndDonor, &pyArrayTypes, &pyArrayCoefs,
                    &pyArrayPressure,
                    &pyArrayGradxP , &pyArrayGradyP , &pyArrayGradzP,
                    &pyArrayGradxxP, &pyArrayGradxyP, &pyArrayGradxzP,
                    &pyArrayGradyxP, &pyArrayGradyyP, &pyArrayGradyzP,
                    &pyArrayGradzxP, &pyArrayGradzyP, &pyArrayGradzzP,
                    &loc,
                    &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters))
    {
      return NULL;
    }

  vector<PyArrayObject*> hook;

  E_Int nvars = 13;

  E_Int imdjmd, imd, jmd, kmd, cnNfldD, ndimdxR, ndimdxD, meshtype;
  E_Float* iptroD; E_Float* iptroR;

  # include "extract_interpD.h"
  FldArrayI* rcvPtsI;
  K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI);
  E_Int* rcvPts  = rcvPtsI->begin();
  nbRcvPts       = rcvPtsI->getSize();

  FldArrayF* pressF;
  E_Int okP = K_NUMPY::getFromNumpyArray(pyArrayPressure, pressF);
  E_Float* pressure = pressF->begin();

  FldArrayF* gradxPressF; FldArrayF* gradyPressF; FldArrayF* gradzPressF;
  E_Int okGxP = K_NUMPY::getFromNumpyArray(pyArrayGradxP, gradxPressF);
  E_Int okGyP = K_NUMPY::getFromNumpyArray(pyArrayGradyP, gradyPressF);
  E_Int okGzP = K_NUMPY::getFromNumpyArray(pyArrayGradzP, gradzPressF);
  E_Float* gradxP = gradxPressF->begin();
  E_Float* gradyP = gradyPressF->begin();
  E_Float* gradzP = gradzPressF->begin();

  FldArrayF* gradxxPressF; FldArrayF* gradxyPressF; FldArrayF* gradxzPressF;
  E_Int okGxxP = K_NUMPY::getFromNumpyArray(pyArrayGradxxP, gradxxPressF);
  E_Int okGxyP = K_NUMPY::getFromNumpyArray(pyArrayGradxyP, gradxyPressF);
  E_Int okGxzP = K_NUMPY::getFromNumpyArray(pyArrayGradxzP, gradxzPressF);
  E_Float* gradxxP = gradxxPressF->begin();
  E_Float* gradxyP = gradxyPressF->begin();
  E_Float* gradxzP = gradxzPressF->begin();

  FldArrayF* gradyxPressF; FldArrayF* gradyyPressF; FldArrayF* gradyzPressF;
  E_Int okGyxP = K_NUMPY::getFromNumpyArray(pyArrayGradyxP, gradyxPressF);
  E_Int okGyyP = K_NUMPY::getFromNumpyArray(pyArrayGradyyP, gradyyPressF);
  E_Int okGyzP = K_NUMPY::getFromNumpyArray(pyArrayGradyzP, gradyzPressF);
  E_Float* gradyxP = gradyxPressF->begin();
  E_Float* gradyyP = gradyyPressF->begin();
  E_Float* gradyzP = gradyzPressF->begin();

  FldArrayF* gradzxPressF; FldArrayF* gradzyPressF; FldArrayF* gradzzPressF;
  E_Int okGzxP = K_NUMPY::getFromNumpyArray(pyArrayGradzxP, gradzxPressF);
  E_Int okGzyP = K_NUMPY::getFromNumpyArray(pyArrayGradzyP, gradzyPressF);
  E_Int okGzzP = K_NUMPY::getFromNumpyArray(pyArrayGradzzP, gradzzPressF);
  E_Float* gradzxP = gradzxPressF->begin();
  E_Float* gradzyP = gradzyPressF->begin();
  E_Float* gradzzP = gradzzPressF->begin();

  vector<E_Float*> fieldsR; vector<E_Float*> fieldsD;
  vector<E_Float*> vectOfDnrFields(nvars); vector<E_Float*> vectOfRcvFields(nvars);
  E_Int* ptrcnd;
  char* eltTypeR; char* eltTypeD;

  // Data are not compacted beforehand
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

  meshtype = resd; // 1: structured, 2: unstructured
  E_Int imr, jmr, kmr, cnSizeR, cnNfldR;
  char* varStringR; vector<E_Int> locsR;
  vector<E_Int*> cnr;
  K_PYTREE::getFromZone(zoneR, 0, loc, varStringR,
                        fieldsR, locsR, imr, jmr, kmr,
                        cnr, cnSizeR, cnNfldR, eltTypeR, hook,
                        GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);

  E_Int posvr, posvd;

  E_Int nfoundvar = 0;
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
          if (posvr != -1)
          {
            vectOfRcvFields[nfoundvar]= fieldsR[posvr];
            vectOfDnrFields[nfoundvar]= fieldsD[posvd];
            nfoundvar += 1;
          }
        }
        #if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(tpl0))
        {
          const char* varname = PyUnicode_AsUTF8(tpl0);
          posvd = K_ARRAY::isNamePresent(varname, varStringD);
          posvr = K_ARRAY::isNamePresent(varname, varStringR);
          if (posvr != -1)
          {
            vectOfRcvFields[nfoundvar]= fieldsR[posvr];
            vectOfDnrFields[nfoundvar]= fieldsD[posvd];
            nfoundvar += 1;
          }
        }
        #endif
        else
        {
          PyErr_Warn(PyExc_Warning, "_setIBCTransfers: variable must be a string. Skipped.");
        }
      }
    }
  }

  delete [] varStringR; delete [] varStringD; delete [] eltTypeR; delete [] eltTypeD;

  # include "commonInterpTransfers_indirect.h"

  # pragma omp parallel default(shared)
  {
    E_Int ideb, ifin;
    #ifdef _OPENMP
    E_Int  ithread           = omp_get_thread_num()+1;
    E_Int  Nbre_thread_actif = omp_get_num_threads(); 
    #else
    E_Int ithread = 1;
    E_Int Nbre_thread_actif = 1;
    #endif
    E_Int chunk = nbRcvPts/Nbre_thread_actif;
    E_Int r = nbRcvPts - chunk*Nbre_thread_actif;
    if (ithread <= r)
    { 
      ideb = (ithread-1)*(chunk+1); ifin = ideb + (chunk+1); 
    }
    else 
    { 
      ideb = (chunk+1)*r+(ithread-r-1)*chunk; ifin = ideb + chunk; 
    }

    // 0 : Pressure
    // 1 / 2 / 3 : gradxPressure  / gradyPressure  / gradzPressure
    // 4 / 5 / 6 : gradxxPressure / gradyxPressure / gradzxPressure
    // 7 / 8 / 9 : gradxyPressure / gradyyPressure / gradyyPressure
    // 10/11 /12 : gradxzPressure / gradyzPressure / gradzzPressure

    #ifdef _OPENMP4
    #pragma omp simd
    #endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      E_Int indR = rcvPts[noind+ideb];
        
      pressure[noind+ideb] = vectOfRcvFields[0][indR];

      gradxP[noind+ideb] = vectOfRcvFields[1][indR];
      gradyP[noind+ideb] = vectOfRcvFields[2][indR];
      gradzP[noind+ideb] = vectOfRcvFields[3][indR];

      gradxxP[noind+ideb] = vectOfRcvFields[4][indR];
      gradxyP[noind+ideb] = vectOfRcvFields[7][indR];
      gradxzP[noind+ideb] = vectOfRcvFields[10][indR];

      gradyxP[noind+ideb] = vectOfRcvFields[5][indR];
      gradyyP[noind+ideb] = vectOfRcvFields[8][indR];
      gradyzP[noind+ideb] = vectOfRcvFields[11][indR];

      gradzxP[noind+ideb] = vectOfRcvFields[6][indR];
      gradzyP[noind+ideb] = vectOfRcvFields[9][indR];
      gradzzP[noind+ideb] = vectOfRcvFields[12][indR];
      }
  }

  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);

  RELEASESHAREDN(pyIndRcv, rcvPtsI);
  RELEASESHAREDN(pyIndDonor, donorPtsI);
  RELEASESHAREDN(pyArrayTypes, typesI);
  RELEASESHAREDN(pyArrayCoefs, donorCoefsF);

  RELEASESHAREDN(pyArrayPressure, pressF);
  RELEASESHAREDN(pyArrayGradxP, gradxPressF);
  RELEASESHAREDN(pyArrayGradyP, gradyPressF);
  RELEASESHAREDN(pyArrayGradzP, gradzPressF);

  RELEASESHAREDN(pyArrayGradxxP, gradxxPressF);
  RELEASESHAREDN(pyArrayGradxyP, gradxyPressF);
  RELEASESHAREDN(pyArrayGradxzP, gradxzPressF);

  RELEASESHAREDN(pyArrayGradyxP, gradyxPressF);
  RELEASESHAREDN(pyArrayGradyyP, gradyyPressF);
  RELEASESHAREDN(pyArrayGradyzP, gradyzPressF);

  RELEASESHAREDN(pyArrayGradzxP, gradzxPressF);
  RELEASESHAREDN(pyArrayGradzyP, gradzyPressF);
  RELEASESHAREDN(pyArrayGradzzP, gradzzPressF);

  Py_INCREF(Py_None);
  return Py_None;
}
//=============================================================================
// Copy of _setIBCTransfers for gradP info
// tc/tc2 -> RCV ZONES
// Called by XOD._setIBCTransfers4GradP()
//=============================================================================
PyObject* K_CONNECTOR::_setIBCTransfers4GradP(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  PyObject *zoneR, *zoneD;
  PyObject *pyVariables;
  PyObject *pyIndRcv, *pyIndDonor;
  PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs;
  PyObject *pyArrayPressure;
  PyObject *pyArrayGradxP, *pyArrayGradyP, *pyArrayGradzP;
  E_Int bctype, loc, vartype, compact;
  E_Float gamma, cv, muS, Cs, Ts, alpha;

  if (!PYPARSETUPLE_(args,
                    OOOO_ OOOO_ OOO_ IIII_ RRRR_ RR_ SSS_,
                    &zoneR, &zoneD, &pyVariables,
                    &pyIndRcv  , &pyIndDonor, &pyArrayTypes, &pyArrayCoefs,
                    &pyArrayPressure,
                    &pyArrayGradxP, &pyArrayGradyP, &pyArrayGradzP,
                    &bctype    , &loc       , &vartype   , &compact   ,&gamma, &cv, &muS, &Cs, &Ts, &alpha,
                    &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters))
  {
    return NULL;
  }

  vector<PyArrayObject*> hook;

  E_Int bcType = E_Int(bctype);

  E_Int nvars;
  E_Int varType = E_Int(vartype);

  nvars = 4;

  // recupere les champs du donneur (nodes)
  E_Int imdjmd, imd, jmd, kmd, cnNfldD, ndimdxR, ndimdxD, meshtype;;
  E_Float* iptroD; E_Float* iptroR;

  # include "extract_interpD.h"

  /*--------------------------------------*/
  /* Extraction des indices des receveurs */
  /*--------------------------------------*/
  FldArrayI* rcvPtsI;
  K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI);
  E_Int* rcvPts  = rcvPtsI->begin();
  nbRcvPts       = rcvPtsI->getSize();

  FldArrayF* pressF;
  E_Int okP = K_NUMPY::getFromNumpyArray( pyArrayPressure, pressF);
  E_Float* pressure = pressF->begin();

  FldArrayF* gradxPressF; FldArrayF* gradyPressF; FldArrayF* gradzPressF;
  E_Int okGxP = K_NUMPY::getFromNumpyArray(pyArrayGradxP, gradxPressF);
  E_Int okGyP = K_NUMPY::getFromNumpyArray(pyArrayGradyP, gradyPressF);
  E_Int okGzP = K_NUMPY::getFromNumpyArray(pyArrayGradzP, gradzPressF);
  E_Float* gradxP = gradxPressF->begin();
  E_Float* gradyP = gradyPressF->begin();
  E_Float* gradzP = gradzPressF->begin();

  vector<E_Float*> fieldsR;vector<E_Float*> fieldsD;
  vector<E_Float*> vectOfDnrFields(nvars); vector<E_Float*> vectOfRcvFields(nvars);
  E_Int* ptrcnd;
  char* eltTypeR; char* eltTypeD;
  //codage general (lent ;-) )
  if (compact == 0)
    {// recupere les champs du donneur (nodes)
      // std::cout << "coucou PAS COMPACT" << std::endl;
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

      // printf("%s\n", varStringD);

      if (cnd.size() > 0) ptrcnd = cnd[0];

      meshtype = resd; // 1: structure, 2: non structure
      // recupere les champs du receveur
      E_Int imr, jmr, kmr, cnSizeR, cnNfldR;
      char* varStringR; vector<E_Int> locsR;
      vector<E_Int*> cnr;
      K_PYTREE::getFromZone(zoneR, 0, loc, varStringR,
                            fieldsR, locsR, imr, jmr, kmr,
                            cnr, cnSizeR, cnNfldR, eltTypeR, hook,
                            GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);

      // printf("%s\n", varStringR);


      // -- no check (perfo) --
      // Transferts
      // Types valides: 2, 3, 4, 5
      E_Int posvr, posvd;

      // Extrait les positions des variables a transferer
      E_Int nfoundvar = 0;
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
          if (posvr != -1)
      {
        vectOfRcvFields[nfoundvar]= fieldsR[posvr];
        vectOfDnrFields[nfoundvar]= fieldsD[posvd];
        nfoundvar += 1;
      }
        }
          #if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(tpl0))
        {
          const char* varname = PyUnicode_AsUTF8(tpl0);
          posvd = K_ARRAY::isNamePresent(varname, varStringD);
          posvr = K_ARRAY::isNamePresent(varname, varStringR);
          if (posvr != -1)
      {
        vectOfRcvFields[nfoundvar]= fieldsR[posvr];
        vectOfDnrFields[nfoundvar]= fieldsD[posvd];
        nfoundvar += 1;
      }
        }
          #endif
      else
        {
          PyErr_Warn(PyExc_Warning, "_setIBCTransfers: variable must be a string. Skipped.");
        }

    }
      }
  }

      delete [] varStringR; delete [] varStringD; delete [] eltTypeR; delete [] eltTypeD;

    }

  // les variables a transferes sont compactes: on recuperes uniquement la premiere et la taille
  else
    {
      // std::cout << "coucou COMPACT" << std::endl;
      # include "getfromzonecompact.h"
      for (E_Int eq = 0; eq < nvars; eq++)
  {
    vectOfRcvFields[eq]= iptroR + eq*ndimdxR;
    vectOfDnrFields[eq]= iptroD + eq*ndimdxD;
  }
    }

    # include "commonInterpTransfers_indirect.h"

  # pragma omp parallel default(shared)
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

    E_Float cvgam = cv*(gamma-1.);

    // 0 : Pressure
    // 1 / 2 / 3 : gradxPressure  / gradyPressure  / gradzPressure

    #ifdef _OPENMP4
    #pragma omp simd
    #endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];
        
        pressure[noind+ideb] = vectOfRcvFields[0][indR];

        gradxP[noind+ideb] = vectOfRcvFields[1][indR];
        gradyP[noind+ideb] = vectOfRcvFields[2][indR];
        gradzP[noind+ideb] = vectOfRcvFields[3][indR];
      }

  } // Fin zone // omp
  // sortie
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);

  RELEASESHAREDN(pyIndRcv, rcvPtsI);
  RELEASESHAREDN(pyIndDonor, donorPtsI);
  RELEASESHAREDN(pyArrayTypes, typesI);
  RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
  // BLOCKRELEASEMEM;

  RELEASESHAREDN(pyArrayPressure, pressF);
  RELEASESHAREDN(pyArrayGradxP, gradxPressF);
  RELEASESHAREDN(pyArrayGradyP, gradyPressF);
  RELEASESHAREDN(pyArrayGradzP, gradzPressF);
  // BLOCKRELEASEMEM2;

  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
// Copy of _setIBCTransfers for gradP HO info
// tc/tc2 -> RCV ZONES
// Called by XOD._setIBCTransfers4GradP2()
//=============================================================================
PyObject* K_CONNECTOR::_setIBCTransfers4GradP2(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  PyObject *zoneR, *zoneD;
  PyObject *pyVariables;
  PyObject *pyIndRcv, *pyIndDonor;
  PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs;
  PyObject *pyArrayPressure;
  PyObject *pyArrayGradxP,  *pyArrayGradyP,  *pyArrayGradzP;
  PyObject *pyArrayGradxxP, *pyArrayGradxyP, *pyArrayGradxzP;
  PyObject *pyArrayGradyxP, *pyArrayGradyyP, *pyArrayGradyzP;
  PyObject *pyArrayGradzxP, *pyArrayGradzyP, *pyArrayGradzzP;
  E_Int bctype, loc, vartype, compact;
  E_Float gamma, cv, muS, Cs, Ts, alpha;

  if (!PYPARSETUPLE_(args,
                    OOOO_ OOOO_ OOOO_ OOOO_ OOOO_ IIII_ RRRR_ RR_ SSS_,
                    &zoneR, &zoneD, &pyVariables,
                    &pyIndRcv  , &pyIndDonor, &pyArrayTypes, &pyArrayCoefs,
                    &pyArrayPressure,
                    &pyArrayGradxP,  &pyArrayGradyP,  &pyArrayGradzP,
                    &pyArrayGradxxP, &pyArrayGradxyP, &pyArrayGradxzP,
                    &pyArrayGradyxP, &pyArrayGradyyP, &pyArrayGradyzP,
                    &pyArrayGradzxP, &pyArrayGradzyP, &pyArrayGradzzP,
                    &bctype    , &loc       , &vartype   , &compact   ,&gamma, &cv, &muS, &Cs, &Ts, &alpha,
                    &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters))
    {
      return NULL;
    }

  vector<PyArrayObject*> hook;

  E_Int bcType = E_Int(bctype);

  E_Int nvars;
  E_Int varType = E_Int(vartype);

  nvars = 13;

  // recupere les champs du donneur (nodes)
  E_Int imdjmd, imd, jmd, kmd, cnNfldD, ndimdxR, ndimdxD, meshtype;;
  E_Float* iptroD; E_Float* iptroR;

# include "extract_interpD.h"

  /*--------------------------------------*/
  /* Extraction des indices des receveurs */
  /*--------------------------------------*/
  FldArrayI* rcvPtsI;
  K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI);
  E_Int* rcvPts  = rcvPtsI->begin();
  nbRcvPts       = rcvPtsI->getSize();

  FldArrayF* pressF;
  E_Int okP = K_NUMPY::getFromNumpyArray( pyArrayPressure, pressF);
  E_Float* pressure = pressF->begin();

  FldArrayF* gradxPressF; FldArrayF* gradyPressF; FldArrayF* gradzPressF;
  E_Int okGxP = K_NUMPY::getFromNumpyArray(pyArrayGradxP, gradxPressF);
  E_Int okGyP = K_NUMPY::getFromNumpyArray(pyArrayGradyP, gradyPressF);
  E_Int okGzP = K_NUMPY::getFromNumpyArray(pyArrayGradzP, gradzPressF);
  E_Float* gradxP = gradxPressF->begin();
  E_Float* gradyP = gradyPressF->begin();
  E_Float* gradzP = gradzPressF->begin();

  FldArrayF* gradxxPressF; FldArrayF* gradxyPressF; FldArrayF* gradxzPressF;
  E_Int okGxxP = K_NUMPY::getFromNumpyArray(pyArrayGradxxP, gradxxPressF);
  E_Int okGxyP = K_NUMPY::getFromNumpyArray(pyArrayGradxyP, gradxyPressF);
  E_Int okGxzP = K_NUMPY::getFromNumpyArray(pyArrayGradxzP, gradxzPressF);
  E_Float* gradxxP = gradxxPressF->begin();
  E_Float* gradxyP = gradxyPressF->begin();
  E_Float* gradxzP = gradxzPressF->begin();

  FldArrayF* gradyxPressF; FldArrayF* gradyyPressF; FldArrayF* gradyzPressF;
  E_Int okGyxP = K_NUMPY::getFromNumpyArray(pyArrayGradyxP, gradyxPressF);
  E_Int okGyyP = K_NUMPY::getFromNumpyArray(pyArrayGradyyP, gradyyPressF);
  E_Int okGyzP = K_NUMPY::getFromNumpyArray(pyArrayGradyzP, gradyzPressF);
  E_Float* gradyxP = gradyxPressF->begin();
  E_Float* gradyyP = gradyyPressF->begin();
  E_Float* gradyzP = gradyzPressF->begin();

  FldArrayF* gradzxPressF; FldArrayF* gradzyPressF; FldArrayF* gradzzPressF;
  E_Int okGzxP = K_NUMPY::getFromNumpyArray(pyArrayGradzxP, gradzxPressF );
  E_Int okGzyP = K_NUMPY::getFromNumpyArray(pyArrayGradzyP, gradzyPressF);
  E_Int okGzzP = K_NUMPY::getFromNumpyArray(pyArrayGradzzP, gradzzPressF);
  E_Float* gradzxP = gradzxPressF->begin();
  E_Float* gradzyP = gradzyPressF->begin();
  E_Float* gradzzP = gradzzPressF->begin();

  vector<E_Float*> fieldsR;vector<E_Float*> fieldsD;
  vector<E_Float*> vectOfDnrFields(nvars); vector<E_Float*> vectOfRcvFields(nvars);
  E_Int* ptrcnd;
  char* eltTypeR; char* eltTypeD;
  //codage general (lent ;-) )
  if (compact == 0)
    {// recupere les champs du donneur (nodes)
      // std::cout << "coucou PAS COMPACT" << std::endl;
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
      // recupere les champs du receveur
      E_Int imr, jmr, kmr, cnSizeR, cnNfldR;
      char* varStringR; vector<E_Int> locsR;
      vector<E_Int*> cnr;
      K_PYTREE::getFromZone(zoneR, 0, loc, varStringR,
                            fieldsR, locsR, imr, jmr, kmr,
                            cnr, cnSizeR, cnNfldR, eltTypeR, hook,
                            GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);

      // -- no check (perfo) --
      // Transferts
      // Types valides: 2, 3, 4, 5
      E_Int posvr, posvd;

      // Extrait les positions des variables a transferer
      E_Int nfoundvar = 0;
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
          if (posvr != -1)
      {
        vectOfRcvFields[nfoundvar]= fieldsR[posvr];
        vectOfDnrFields[nfoundvar]= fieldsD[posvd];
        nfoundvar += 1;
      }
        }
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(tpl0))
        {
          const char* varname = PyUnicode_AsUTF8(tpl0);
          posvd = K_ARRAY::isNamePresent(varname, varStringD);
          posvr = K_ARRAY::isNamePresent(varname, varStringR);
          if (posvr != -1)
      {
        vectOfRcvFields[nfoundvar]= fieldsR[posvr];
        vectOfDnrFields[nfoundvar]= fieldsD[posvd];
        nfoundvar += 1;
      }
        }
#endif
      else
        {
          PyErr_Warn(PyExc_Warning, "_setIBCTransfers: variable must be a string. Skipped.");
        }

    }
      }
  }

      delete [] varStringR; delete [] varStringD; delete [] eltTypeR; delete [] eltTypeD;

    }

  // les variables a transferes sont compactes: on recuperes uniquement la premiere et la taille
  else
    {
      // std::cout << "coucou COMPACT" << std::endl;
# include "getfromzonecompact.h"
      for (E_Int eq = 0; eq < nvars; eq++)
  {
    vectOfRcvFields[eq]= iptroR + eq*ndimdxR;
    vectOfDnrFields[eq]= iptroD + eq*ndimdxD;
  }
    }

# include "commonInterpTransfers_indirect.h"

#    pragma omp parallel default(shared)
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

    E_Float cvgam = cv*(gamma-1.);

    // 0 : Pressure
    // 1 / 2 / 3 : gradxPressure  / gradyPressure  / gradzPressure
    // 4 / 5 / 6 : gradxxPressure / gradyxPressure / gradzxPressure
    // 7 / 8 / 9 : gradxyPressure / gradyyPressure / gradyyPressure
    // 10/11 /12 : gradxzPressure / gradyzPressure / gradzzPressure

#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
  E_Int indR = rcvPts[noind+ideb];
     
  pressure[noind+ideb] = vectOfRcvFields[0][indR];

  gradxP[noind+ideb] = vectOfRcvFields[1][indR];
  gradyP[noind+ideb] = vectOfRcvFields[2][indR];
  gradzP[noind+ideb] = vectOfRcvFields[3][indR];

  gradxxP[noind+ideb] = vectOfRcvFields[4][indR];
  gradxyP[noind+ideb] = vectOfRcvFields[7][indR];
  gradxzP[noind+ideb] = vectOfRcvFields[10][indR];

  gradyxP[noind+ideb] = vectOfRcvFields[5][indR];
  gradyyP[noind+ideb] = vectOfRcvFields[8][indR];
  gradyzP[noind+ideb] = vectOfRcvFields[11][indR];

  gradzxP[noind+ideb] = vectOfRcvFields[6][indR];
  gradzyP[noind+ideb] = vectOfRcvFields[9][indR];
  gradzzP[noind+ideb] = vectOfRcvFields[12][indR];
      }

  } // Fin zone // omp
  // sortie
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);

  RELEASESHAREDN(pyIndRcv, rcvPtsI);
  RELEASESHAREDN(pyIndDonor, donorPtsI);
  RELEASESHAREDN(pyArrayTypes, typesI);
  RELEASESHAREDN(pyArrayCoefs, donorCoefsF);

  RELEASESHAREDN(pyArrayPressure, pressF);
  RELEASESHAREDN(pyArrayGradxP, gradxPressF);
  RELEASESHAREDN(pyArrayGradyP, gradyPressF);
  RELEASESHAREDN(pyArrayGradzP, gradzPressF);

  RELEASESHAREDN(pyArrayGradxxP, gradxxPressF);
  RELEASESHAREDN(pyArrayGradxyP, gradxyPressF);
  RELEASESHAREDN(pyArrayGradxzP, gradxzPressF);

  RELEASESHAREDN(pyArrayGradyxP, gradyxPressF);
  RELEASESHAREDN(pyArrayGradyyP, gradyyPressF);
  RELEASESHAREDN(pyArrayGradyzP, gradyzPressF);

  RELEASESHAREDN(pyArrayGradzxP, gradzxPressF);
  RELEASESHAREDN(pyArrayGradzyP, gradzyPressF);
  RELEASESHAREDN(pyArrayGradzzP, gradzzPressF);

  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
// Copy of _setIBCTransfers for gradP info - OPTI + COMPACT 
// RCV ZONES -> tc
// transfers of ['Density', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature', 'nuSA'] -> 6 vars
// transfers of ['gradx/y/zDensity', 'gradx/y/zTemperature'] -> 6 varsGrad
// Called by XOD._setIBCTransfers4GradP3()
//=============================================================================
PyObject* K_CONNECTOR::_setIBCTransfers4GradP3(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  PyObject *zoneR;
  PyObject *pyVariables;
  PyObject *pyIndRcv, *pyIndDonor;
  PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs;
  PyObject *pyArrayGradxP, *pyArrayGradyP, *pyArrayGradzP;
  E_Int bctype, loc, vartype, compact;
  E_Float gamma, cv, muS, Cs, Ts, alpha;

  if (!PYPARSETUPLE_(args,
                    OOOO_ OOOO_ O_ IIII_ RRRR_ RR_ SSS_,
                    &zoneR, &pyVariables,
                    &pyIndRcv  , &pyIndDonor, &pyArrayTypes, &pyArrayCoefs,
                    &pyArrayGradxP, &pyArrayGradyP, &pyArrayGradzP,
                    &bctype    , &loc       , &vartype   , &compact   ,&gamma, &cv, &muS, &Cs, &Ts, &alpha,
                    &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters))
    {
      return NULL;
    }

  vector<PyArrayObject*> hook;

  E_Int bcType = E_Int(bctype);

  E_Int nvars, nvars_grad;
  E_Int varType = E_Int(vartype);

  nvars = 6;
  nvars_grad = 6;

  // recupere les champs du donneur (nodes)
  E_Int imdjmd, imd, jmd, kmd, ndimdxR, meshtype;
  E_Float* iptroR;
  E_Float* iptgradR;

# include "extract_interpD.h"

  /*--------------------------------------*/
  /* Extraction des indices des receveurs */
  /*--------------------------------------*/
  FldArrayI* rcvPtsI;
  K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI);
  E_Int* rcvPts  = rcvPtsI->begin();
  nbRcvPts       = rcvPtsI->getSize();

  FldArrayF* gradxPressF; FldArrayF* gradyPressF; FldArrayF* gradzPressF;
  E_Int okGxP = K_NUMPY::getFromNumpyArray(pyArrayGradxP, gradxPressF);
  E_Int okGyP = K_NUMPY::getFromNumpyArray(pyArrayGradyP, gradyPressF);
  E_Int okGzP = K_NUMPY::getFromNumpyArray(pyArrayGradzP, gradzPressF);
  E_Float* gradxP = gradxPressF->begin();
  E_Float* gradyP = gradyPressF->begin();
  E_Float* gradzP = gradzPressF->begin();

  vector<E_Float*> vectOfRcvFields(nvars);
  vector<E_Float*> vectOfGradRcvFields(nvars_grad);

  // les variables a transferes sont compactes: on recuperes uniquement la premiere et la taille
  //##############################
  PyObject* solR;
  PyObject* t;
  char* type; E_Int s, s0, s1;  E_Int* d;

  PyObject* tpl0= PyList_GetItem(pyVariables, 0);
  char* varname = NULL;
  if (PyString_Check(tpl0)) varname = PyString_AsString(tpl0);
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(tpl0)) varname = (char*)PyUnicode_AsUTF8(tpl0);
#endif

  if  (loc==0) { solR = K_PYTREE::getNodeFromName1(zoneR , "FlowSolution"        ); }
  else  { solR = K_PYTREE::getNodeFromName1(zoneR , "FlowSolution#Centers"); }
  t = K_PYTREE::getNodeFromName1(solR, varname);
  iptroR = K_PYTREE::getValueAF(t, hook);

  if  (loc==0) { solR = K_PYTREE::getNodeFromName1(zoneR , "FlowSolution"        ); }
  else  { solR = K_PYTREE::getNodeFromName1(zoneR , "FlowSolution#Centers"); }
  t = K_PYTREE::getNodeFromName1(solR, "gradxDensity");
  iptgradR = K_PYTREE::getValueAF(t, hook);

  // get type
  t =  K_PYTREE::getNodeFromName1(zoneR, "ZoneType");
  type =  K_PYTREE::getValueS(t, s, hook);
  // get dims zone receveuse
  d  =  K_PYTREE::getValueAI(zoneR, s0, s1, hook);

  if  (K_STRING::cmp(type, s, "Structured") == 0)
    {
      E_Int shift = 0; if(loc == 1) shift = 3;
      if (s0 == 1) { ndimdxR= d[0+shift]; }
      else if (s0 == 2) { ndimdxR= d[0+shift]*d[1+shift]; } 
      else if (s0 == 3) { ndimdxR= d[0+shift]*d[1+shift]*d[2+shift]; } 
    }
  else // non structure
    {
      ndimdxR= d[0]* d[1]; // npoint, nelements
    }
  //##############################

  for (E_Int eq = 0; eq < nvars; eq++)
    {
      vectOfRcvFields[eq]= iptroR + eq*ndimdxR;
    }

  for (E_Int eq = 0; eq < nvars_grad; eq++)
    {
      vectOfGradRcvFields[eq]= iptgradR + eq*ndimdxR;
    }

  // # include "commonInterpTransfers_indirect.h"

#    pragma omp parallel default(shared)
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

    E_Float cvgam = cv*(gamma-1.);

#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
  E_Int indR = rcvPts[noind+ideb];

  gradxP[noind+ideb] = ((vectOfRcvFields[4][indR]*vectOfGradRcvFields[0][indR]+vectOfRcvFields[0][indR]*vectOfGradRcvFields[3][indR])*cvgam)/alpha + gradxP[noind+ideb]*(alpha-1.)/alpha;
  gradyP[noind+ideb] = ((vectOfRcvFields[4][indR]*vectOfGradRcvFields[1][indR]+vectOfRcvFields[0][indR]*vectOfGradRcvFields[4][indR])*cvgam)/alpha + gradyP[noind+ideb]*(alpha-1.)/alpha;
  gradzP[noind+ideb] = ((vectOfRcvFields[4][indR]*vectOfGradRcvFields[2][indR]+vectOfRcvFields[0][indR]*vectOfGradRcvFields[5][indR])*cvgam)/alpha + gradzP[noind+ideb]*(alpha-1.)/alpha;
      }

  } // Fin zone // omp
  // sortie
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);

  RELEASESHAREDN(pyIndRcv, rcvPtsI);
  RELEASESHAREDN(pyIndDonor, donorPtsI);
  RELEASESHAREDN(pyArrayTypes, typesI);
  RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
  // BLOCKRELEASEMEM;

  RELEASESHAREDN(pyArrayGradxP, gradxPressF);
  RELEASESHAREDN(pyArrayGradyP, gradyPressF);
  RELEASESHAREDN(pyArrayGradzP, gradzPressF);
  // BLOCKRELEASEMEM2;

  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
// Update gradP info in IBCD zones given three numpy gradx/y/zP_new
// Called by XOD._setIBCTransfers4GradP3()
//=============================================================================
PyObject* K_CONNECTOR::_setIBCTransfers4GradP4(PyObject* self, PyObject* args)
{
  PyObject *pyArrayGradxP_new, *pyArrayGradyP_new, *pyArrayGradzP_new;
  PyObject *pyArrayGradxP, *pyArrayGradyP, *pyArrayGradzP;

  if (!PYPARSETUPLE_(args, OOOO_ OO_,
                    &pyArrayGradxP_new, &pyArrayGradyP_new, &pyArrayGradzP_new,
                    &pyArrayGradxP, &pyArrayGradyP, &pyArrayGradzP))
    {
      return NULL;
    }

  FldArrayF* gradxPressF; FldArrayF* gradyPressF; FldArrayF* gradzPressF;
  E_Int okGxP = K_NUMPY::getFromNumpyArray(pyArrayGradxP, gradxPressF);
  E_Int okGyP = K_NUMPY::getFromNumpyArray(pyArrayGradyP, gradyPressF);
  E_Int okGzP = K_NUMPY::getFromNumpyArray(pyArrayGradzP, gradzPressF);
  E_Float* gradxP = gradxPressF->begin();
  E_Float* gradyP = gradyPressF->begin();
  E_Float* gradzP = gradzPressF->begin();

  FldArrayF* gradxPress_newF; FldArrayF* gradyPress_newF; FldArrayF* gradzPress_newF;
  E_Int okGxP_new = K_NUMPY::getFromNumpyArray(pyArrayGradxP_new, gradxPress_newF);
  E_Int okGyP_new = K_NUMPY::getFromNumpyArray(pyArrayGradyP_new, gradyPress_newF);
  E_Int okGzP_new = K_NUMPY::getFromNumpyArray(pyArrayGradzP_new, gradzPress_newF);
  E_Float* gradxP_new = gradxPress_newF->begin();
  E_Float* gradyP_new = gradyPress_newF->begin();
  E_Float* gradzP_new = gradzPress_newF->begin();

  E_Int nbRcvPts = gradzPress_newF->getSize();

#    pragma omp parallel default(shared)
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

#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      gradxP[noind+ideb] = gradxP_new[noind+ideb];
      gradyP[noind+ideb] = gradyP_new[noind+ideb];
      gradzP[noind+ideb] = gradzP_new[noind+ideb];
    }

  } // Fin zone // omp
  // sortie

  RELEASESHAREDN(pyArrayGradxP, gradxPressF);
  RELEASESHAREDN(pyArrayGradyP, gradyPressF);
  RELEASESHAREDN(pyArrayGradzP, gradzPressF);
  RELEASESHAREDN(pyArrayGradxP_new, gradxPress_newF);
  RELEASESHAREDN(pyArrayGradyP_new, gradyPress_newF);
  RELEASESHAREDN(pyArrayGradzP_new, gradzPress_newF);

  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
// Copy of _setIBCTransfers for gradP + gradVelocity info
// tc/tc2 -> RCV ZONES
// Called by XOD._setIBCTransfers4FULLTBLE()
//=============================================================================
PyObject* K_CONNECTOR::_setIBCTransfers4FULLTBLE(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  PyObject *zoneR, *zoneD;
  PyObject *pyVariables;
  PyObject *pyIndRcv, *pyIndDonor;
  PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs;
  PyObject *pyArrayPressure;
  PyObject *pyArrayU, *pyArrayV, *pyArrayW;
  PyObject *pyArrayGradxP, *pyArrayGradyP, *pyArrayGradzP;
  PyObject *pyArrayGradxU, *pyArrayGradyU, *pyArrayGradzU;
  PyObject *pyArrayGradxV, *pyArrayGradyV, *pyArrayGradzV;
  PyObject *pyArrayGradxW, *pyArrayGradyW, *pyArrayGradzW;
  E_Int bctype, loc, vartype, compact;
  E_Float gamma, cv, muS, Cs, Ts, alpha;

  if (!PYPARSETUPLE_(args,
                    OOOO_ OOOO_ OOOO_ OOOO_ OOOO_ OOO_ IIII_ RRRR_ RR_ SSS_,
                    &zoneR, &zoneD, &pyVariables,
                    &pyIndRcv  , &pyIndDonor, &pyArrayTypes, &pyArrayCoefs,
                    &pyArrayPressure,
                    &pyArrayU, &pyArrayV, &pyArrayW,
                    &pyArrayGradxP, &pyArrayGradyP, &pyArrayGradzP,
                    &pyArrayGradxU, &pyArrayGradyU, &pyArrayGradzU,
                    &pyArrayGradxV, &pyArrayGradyV, &pyArrayGradzV,
                    &pyArrayGradxW, &pyArrayGradyW, &pyArrayGradzW,
                    &bctype    , &loc       , &vartype   , &compact   ,&gamma, &cv, &muS, &Cs, &Ts, &alpha,
                    &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters))
  {
    return NULL;
  }

  vector<PyArrayObject*> hook;

  E_Int bcType = E_Int(bctype);

  E_Int nvars;
  E_Int varType = E_Int(vartype);

  nvars = 20;

  // recupere les champs du donneur (nodes)
  E_Int imdjmd, imd, jmd, kmd, cnNfldD, ndimdxR, ndimdxD, meshtype;;
  E_Float* iptroD; E_Float* iptroR;

# include "extract_interpD.h"

  /*--------------------------------------*/
  /* Extraction des indices des receveurs */
  /*--------------------------------------*/
  FldArrayI* rcvPtsI;
  K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI);
  E_Int* rcvPts  = rcvPtsI->begin();
  nbRcvPts       = rcvPtsI->getSize();

  // # include "IBC/extract_IBC.h"

  FldArrayF* pressF;
  E_Int okP = K_NUMPY::getFromNumpyArray( pyArrayPressure, pressF);
  E_Float* pressure = pressF->begin();

  FldArrayF* UF; FldArrayF* VF; FldArrayF* WF;
  E_Int okU = K_NUMPY::getFromNumpyArray(pyArrayU , UF);
  E_Int okV = K_NUMPY::getFromNumpyArray(pyArrayV , VF);
  E_Int okW = K_NUMPY::getFromNumpyArray(pyArrayW , WF);
  E_Float* U = UF->begin();
  E_Float* V = VF->begin();
  E_Float* W = WF->begin();

  FldArrayF* gradxPressF; FldArrayF* gradyPressF; FldArrayF* gradzPressF;
  E_Int okGxP = K_NUMPY::getFromNumpyArray(pyArrayGradxP, gradxPressF);
  E_Int okGyP = K_NUMPY::getFromNumpyArray(pyArrayGradyP, gradyPressF);
  E_Int okGzP = K_NUMPY::getFromNumpyArray(pyArrayGradzP, gradzPressF);
  E_Float* gradxP = gradxPressF->begin();
  E_Float* gradyP = gradyPressF->begin();
  E_Float* gradzP = gradzPressF->begin();

  FldArrayF* gradxVelocityXF; FldArrayF* gradyVelocityXF; FldArrayF* gradzVelocityXF;
  E_Int okGxU = K_NUMPY::getFromNumpyArray(pyArrayGradxU, gradxVelocityXF);
  E_Int okGyU = K_NUMPY::getFromNumpyArray(pyArrayGradyU, gradyVelocityXF);
  E_Int okGzU = K_NUMPY::getFromNumpyArray(pyArrayGradzU, gradzVelocityXF);
  E_Float* gradxU = gradxVelocityXF->begin();
  E_Float* gradyU = gradyVelocityXF->begin();
  E_Float* gradzU = gradzVelocityXF->begin();

  FldArrayF* gradxVelocityYF; FldArrayF* gradyVelocityYF; FldArrayF* gradzVelocityYF;
  E_Int okGxV = K_NUMPY::getFromNumpyArray(pyArrayGradxV, gradxVelocityYF);
  E_Int okGyV = K_NUMPY::getFromNumpyArray(pyArrayGradyV, gradyVelocityYF);
  E_Int okGzV = K_NUMPY::getFromNumpyArray(pyArrayGradzV, gradzVelocityYF);
  E_Float* gradxV = gradxVelocityYF->begin();
  E_Float* gradyV = gradyVelocityYF->begin();
  E_Float* gradzV = gradzVelocityYF->begin();

  FldArrayF* gradxVelocityZF; FldArrayF* gradyVelocityZF; FldArrayF* gradzVelocityZF;
  E_Int okGxW = K_NUMPY::getFromNumpyArray(pyArrayGradxW, gradxVelocityZF);
  E_Int okGyW = K_NUMPY::getFromNumpyArray(pyArrayGradyW, gradyVelocityZF);
  E_Int okGzW = K_NUMPY::getFromNumpyArray(pyArrayGradzW, gradzVelocityZF);
  E_Float* gradxW = gradxVelocityZF->begin();
  E_Float* gradyW = gradyVelocityZF->begin();
  E_Float* gradzW = gradzVelocityZF->begin();

  // std::cout << "COUCOU 3188" << std::endl;

  vector<E_Float*> fieldsR;vector<E_Float*> fieldsD;
  vector<E_Float*> vectOfDnrFields(nvars); vector<E_Float*> vectOfRcvFields(nvars);
  E_Int* ptrcnd=NULL;
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
    // recupere les champs du receveur
    E_Int imr, jmr, kmr, cnSizeR, cnNfldR;
    char* varStringR; vector<E_Int> locsR;
    vector<E_Int*> cnr;
    K_PYTREE::getFromZone(zoneR, 0, loc, varStringR,
                            fieldsR, locsR, imr, jmr, kmr,
                            cnr, cnSizeR, cnNfldR, eltTypeR, hook,
                            GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);


    // -- no check (perfo) --
    // Transferts
    // Types valides: 2, 3, 4, 5
    E_Int posvr, posvd;

    // Extrait les positions des variables a transferer
    E_Int nfoundvar = 0;
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
            if (posvr != -1)
            {
              vectOfRcvFields[nfoundvar]= fieldsR[posvr];
              vectOfDnrFields[nfoundvar]= fieldsD[posvd];
              nfoundvar += 1;
            }
          }
#if PY_VERSION_HEX >= 0x03000000
          else if (PyUnicode_Check(tpl0))
          {
            const char* varname = PyUnicode_AsUTF8(tpl0);
            posvd = K_ARRAY::isNamePresent(varname, varStringD);
            posvr = K_ARRAY::isNamePresent(varname, varStringR);
            if (posvr != -1)
            {
              vectOfRcvFields[nfoundvar]= fieldsR[posvr];
              vectOfDnrFields[nfoundvar]= fieldsD[posvd];
              nfoundvar += 1;
            }
          }
#endif
        else
        {
          PyErr_Warn(PyExc_Warning, "_setIBCTransfers: variable must be a string. Skipped.");
        }

      }
    }
  }

  delete [] varStringR; delete [] varStringD; delete [] eltTypeR; delete [] eltTypeD;

  }

  // les variables a transferes sont compactes: on recuperes uniquement la premiere et la taille
  else
  {
#   include "getfromzonecompact.h"
    for (E_Int eq = 0; eq < nvars; eq++)
    {
      vectOfRcvFields[eq]= iptroR + eq*ndimdxR;
      vectOfDnrFields[eq]= iptroD + eq*ndimdxD;
    }
  }

# include "commonInterpTransfers_indirect.h"

# pragma omp parallel default(shared)
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

    E_Float cvgam = cv*(gamma-1.);

    //0 : Density
    //1 : Temperature
    //2 / 3 / 4 : gradxDensity / gradyDensity / gradzDensity
    //5 / 6 / 7 : gradxTemperature / gradyTemperature / gradzTemperature
    //8 / 9 / 10 : gradxVelocityX / ...
    //11 / 12 / 13 : gradxVelocityY / ...
    //14 / 15 / 16 : gradxVelocityZ / ...
    //17 / 18 / 19 : VelocityX / ...

    // #ifdef _OPENMP4
    //    #pragma omp simd
    // #endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      E_Int indR = rcvPts[noind+ideb];

      pressure[noind+ideb] = vectOfRcvFields[0][indR]*vectOfRcvFields[1][indR]*cvgam;

      U[noind+ideb] = vectOfRcvFields[17][indR];
      V[noind+ideb] = vectOfRcvFields[18][indR];
      W[noind+ideb] = vectOfRcvFields[19][indR];

      if (alpha <= 1.)
      {
        gradxP[noind+ideb] = ((vectOfRcvFields[1][indR]*vectOfRcvFields[2][indR]+vectOfRcvFields[0][indR]*vectOfRcvFields[5][indR])*cvgam)*alpha;
        gradyP[noind+ideb] = ((vectOfRcvFields[1][indR]*vectOfRcvFields[3][indR]+vectOfRcvFields[0][indR]*vectOfRcvFields[6][indR])*cvgam)*alpha;
        gradzP[noind+ideb] = ((vectOfRcvFields[1][indR]*vectOfRcvFields[4][indR]+vectOfRcvFields[0][indR]*vectOfRcvFields[7][indR])*cvgam)*alpha;

        gradxU[noind+ideb] = (vectOfRcvFields[8][indR])*alpha;  gradyU[noind+ideb] = (vectOfRcvFields[9][indR])*alpha;  gradzU[noind+ideb] = (vectOfRcvFields[10][indR])*alpha;
        gradxV[noind+ideb] = (vectOfRcvFields[11][indR])*alpha; gradyV[noind+ideb] = (vectOfRcvFields[12][indR])*alpha; gradzV[noind+ideb] = (vectOfRcvFields[13][indR])*alpha;
        gradxW[noind+ideb] = (vectOfRcvFields[14][indR])*alpha; gradyW[noind+ideb] = (vectOfRcvFields[15][indR])*alpha; gradzW[noind+ideb] = (vectOfRcvFields[16][indR])*alpha;
      }
      else
      {
        gradxP[noind+ideb] = ((vectOfRcvFields[1][indR]*vectOfRcvFields[2][indR]+vectOfRcvFields[0][indR]*vectOfRcvFields[5][indR])*cvgam)/alpha + gradxP[noind+ideb]*(alpha-1.)/alpha;
        gradyP[noind+ideb] = ((vectOfRcvFields[1][indR]*vectOfRcvFields[3][indR]+vectOfRcvFields[0][indR]*vectOfRcvFields[6][indR])*cvgam)/alpha + gradyP[noind+ideb]*(alpha-1.)/alpha;
        gradzP[noind+ideb] = ((vectOfRcvFields[1][indR]*vectOfRcvFields[4][indR]+vectOfRcvFields[0][indR]*vectOfRcvFields[7][indR])*cvgam)/alpha + gradzP[noind+ideb]*(alpha-1.)/alpha;

        gradxU[noind+ideb] = (vectOfRcvFields[8][indR])/alpha  + gradxU[noind+ideb]*(alpha-1.)/alpha;  
        gradyU[noind+ideb] = (vectOfRcvFields[9][indR])/alpha  + gradyU[noind+ideb]*(alpha-1.)/alpha;    
        gradzU[noind+ideb] = (vectOfRcvFields[10][indR])/alpha + gradzU[noind+ideb]*(alpha-1.)/alpha;  

        gradxV[noind+ideb] = (vectOfRcvFields[11][indR])/alpha + gradxV[noind+ideb]*(alpha-1.)/alpha;  
        gradyV[noind+ideb] = (vectOfRcvFields[12][indR])/alpha + gradyV[noind+ideb]*(alpha-1.)/alpha;    
        gradzV[noind+ideb] = (vectOfRcvFields[13][indR])/alpha + gradzV[noind+ideb]*(alpha-1.)/alpha;  

        gradxW[noind+ideb] = (vectOfRcvFields[14][indR])/alpha + gradxW[noind+ideb]*(alpha-1.)/alpha;  
        gradyW[noind+ideb] = (vectOfRcvFields[15][indR])/alpha + gradyW[noind+ideb]*(alpha-1.)/alpha;    
        gradzW[noind+ideb] = (vectOfRcvFields[16][indR])/alpha + gradzW[noind+ideb]*(alpha-1.)/alpha;  
      }

      // for (E_Int toto = 2; toto < 17; toto++){
      //   vectOfRcvFields[toto][indR] = toto;
      // }

      // vectOfRcvFields[16][indR] = 2.;

    }
  } // Fin zone // omp
  // sortie
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);

  RELEASESHAREDN(pyIndRcv, rcvPtsI);
  RELEASESHAREDN(pyIndDonor, donorPtsI);
  RELEASESHAREDN(pyArrayTypes, typesI);
  RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
  // BLOCKRELEASEMEM;

  RELEASESHAREDN(pyArrayPressure, pressF);
  RELEASESHAREDN(pyArrayU, UF);
  RELEASESHAREDN(pyArrayV, VF);
  RELEASESHAREDN(pyArrayW, WF);
  RELEASESHAREDN(pyArrayGradxP, gradxPressF);
  RELEASESHAREDN(pyArrayGradyP, gradyPressF);
  RELEASESHAREDN(pyArrayGradzP, gradzPressF);
  RELEASESHAREDN(pyArrayGradxU, gradxVelocityXF);
  RELEASESHAREDN(pyArrayGradyU, gradyVelocityXF);
  RELEASESHAREDN(pyArrayGradzU, gradzVelocityXF);
  RELEASESHAREDN(pyArrayGradxV, gradxVelocityYF);
  RELEASESHAREDN(pyArrayGradyV, gradyVelocityYF);
  RELEASESHAREDN(pyArrayGradzV, gradzVelocityYF);
  RELEASESHAREDN(pyArrayGradxW, gradxVelocityZF);
  RELEASESHAREDN(pyArrayGradyW, gradyVelocityZF);
  RELEASESHAREDN(pyArrayGradzW, gradzVelocityZF);
  // BLOCKRELEASEMEM2;

  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
// Copy of _setIBCTransfers for gradP + gradVelocity info
// RCV ZONES -> tc
// Called by XOD._setIBCTransfers4FULLTBLE2()
//=============================================================================
PyObject* K_CONNECTOR::_setIBCTransfers4FULLTBLE2(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  PyObject *zoneR, *zoneD;
  PyObject *pyVariables;
  PyObject *pyIndRcv, *pyIndDonor;
  PyObject *pyArrayTypes;
  PyObject *pyArrayCoefs;
  PyObject *pyArrayPressure;
  PyObject *pyArrayU, *pyArrayV, *pyArrayW;
  PyObject *pyArrayGradxP, *pyArrayGradyP, *pyArrayGradzP;
  PyObject *pyArrayGradxU, *pyArrayGradyU, *pyArrayGradzU;
  PyObject *pyArrayGradxV, *pyArrayGradyV, *pyArrayGradzV;
  PyObject *pyArrayGradxW, *pyArrayGradyW, *pyArrayGradzW;
  E_Int bctype, loc, vartype, compact;
  E_Float gamma, cv, muS, Cs, Ts, alpha;

  if (!PYPARSETUPLE_(args,
                    OOOO_ OOOO_ OOOO_ OOOO_ OOOO_ OOO_ IIII_ RRRR_ RR_ SSS_,
                    &zoneR, &zoneD, &pyVariables,
                    &pyIndRcv  , &pyIndDonor, &pyArrayTypes, &pyArrayCoefs,
                    &pyArrayPressure,
                    &pyArrayU, &pyArrayV, &pyArrayW,
                    &pyArrayGradxP, &pyArrayGradyP, &pyArrayGradzP,
                    &pyArrayGradxU, &pyArrayGradyU, &pyArrayGradzU,
                    &pyArrayGradxV, &pyArrayGradyV, &pyArrayGradzV,
                    &pyArrayGradxW, &pyArrayGradyW, &pyArrayGradzW,
                    &bctype    , &loc       , &vartype   , &compact   ,&gamma, &cv, &muS, &Cs, &Ts, &alpha,
                    &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters))
    {
      return NULL;
    }

  vector<PyArrayObject*> hook;

  //E_Int bcType = E_Int(bctype);

  E_Int nvars;
  E_Int varType = E_Int(vartype);

  nvars = 20;

  // recupere les champs du donneur (nodes)
  E_Int imdjmd, imd, jmd, kmd, cnNfldD, ndimdxR, ndimdxD, meshtype;;
  E_Float* iptroD; E_Float* iptroR;

# include "extract_interpD.h"

  /*--------------------------------------*/
  /* Extraction des indices des receveurs */
  /*--------------------------------------*/
  FldArrayI* rcvPtsI;
  K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI);
  E_Int* rcvPts  = rcvPtsI->begin();
  nbRcvPts       = rcvPtsI->getSize();

  // # include "IBC/extract_IBC.h"

  FldArrayF* pressF;
  E_Int okP = K_NUMPY::getFromNumpyArray(pyArrayPressure, pressF);
  E_Float* pressure = pressF->begin();

  FldArrayF* UF; FldArrayF* VF; FldArrayF* WF;
  E_Int okU = K_NUMPY::getFromNumpyArray(pyArrayU, UF);
  E_Int okV = K_NUMPY::getFromNumpyArray(pyArrayV, VF);
  E_Int okW = K_NUMPY::getFromNumpyArray(pyArrayW, WF);
  //E_Float* U = UF->begin();
  //E_Float* V = VF->begin();
  //E_Float* W = WF->begin();

  FldArrayF* gradxPressF; FldArrayF* gradyPressF; FldArrayF* gradzPressF;
  E_Int okGxP = K_NUMPY::getFromNumpyArray(pyArrayGradxP, gradxPressF);
  E_Int okGyP = K_NUMPY::getFromNumpyArray(pyArrayGradyP, gradyPressF);
  E_Int okGzP = K_NUMPY::getFromNumpyArray(pyArrayGradzP, gradzPressF);
  E_Float* gradxP = gradxPressF->begin();
  E_Float* gradyP = gradyPressF->begin();
  E_Float* gradzP = gradzPressF->begin();

  FldArrayF* gradxVelocityXF; FldArrayF* gradyVelocityXF; FldArrayF* gradzVelocityXF;
  E_Int okGxU = K_NUMPY::getFromNumpyArray(pyArrayGradxU, gradxVelocityXF);
  E_Int okGyU = K_NUMPY::getFromNumpyArray(pyArrayGradyU, gradyVelocityXF);
  E_Int okGzU = K_NUMPY::getFromNumpyArray(pyArrayGradzU, gradzVelocityXF);
  E_Float* gradxU = gradxVelocityXF->begin();
  E_Float* gradyU = gradyVelocityXF->begin();
  E_Float* gradzU = gradzVelocityXF->begin();

  FldArrayF* gradxVelocityYF; FldArrayF* gradyVelocityYF; FldArrayF* gradzVelocityYF;
  E_Int okGxV = K_NUMPY::getFromNumpyArray(pyArrayGradxV, gradxVelocityYF);
  E_Int okGyV = K_NUMPY::getFromNumpyArray(pyArrayGradyV, gradyVelocityYF);
  E_Int okGzV = K_NUMPY::getFromNumpyArray(pyArrayGradzV, gradzVelocityYF);
  E_Float* gradxV = gradxVelocityYF->begin();
  E_Float* gradyV = gradyVelocityYF->begin();
  E_Float* gradzV = gradzVelocityYF->begin();

  FldArrayF* gradxVelocityZF; FldArrayF* gradyVelocityZF; FldArrayF* gradzVelocityZF;
  E_Int okGxW = K_NUMPY::getFromNumpyArray(pyArrayGradxW, gradxVelocityZF);
  E_Int okGyW = K_NUMPY::getFromNumpyArray(pyArrayGradyW, gradyVelocityZF);
  E_Int okGzW = K_NUMPY::getFromNumpyArray(pyArrayGradzW, gradzVelocityZF);
  E_Float* gradxW = gradxVelocityZF->begin();
  E_Float* gradyW = gradyVelocityZF->begin();
  E_Float* gradzW = gradzVelocityZF->begin();

  vector<E_Float*> fieldsR;vector<E_Float*> fieldsD;
  vector<E_Float*> vectOfDnrFields(nvars); vector<E_Float*> vectOfRcvFields(nvars);
  E_Int* ptrcnd;
  char* eltTypeR; char* eltTypeD;
  //codage general (lent ;-) )
  if (compact == 0)
  {
    // recupere les champs du donneur (nodes)
    // std::cout << "COUCOU 5702" << std::endl;
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
    // recupere les champs du receveur
    E_Int imr, jmr, kmr, cnSizeR, cnNfldR;
    char* varStringR; vector<E_Int> locsR;
    vector<E_Int*> cnr;
    K_PYTREE::getFromZone(zoneR, 0, loc, varStringR,
                          fieldsR, locsR, imr, jmr, kmr,
                          cnr, cnSizeR, cnNfldR, eltTypeR, hook,
                          GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);

    // -- no check (perfo) --
    // Transferts
    // Types valides: 2, 3, 4, 5
    E_Int posvr, posvd;

    // Extrait les positions des variables a transferer
    E_Int nfoundvar = 0;
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
            if (posvr != -1)
            {
              vectOfRcvFields[nfoundvar]= fieldsR[posvr];
              vectOfDnrFields[nfoundvar]= fieldsD[posvd];
              nfoundvar += 1;
            }
          }
#if PY_VERSION_HEX >= 0x03000000
          else if (PyUnicode_Check(tpl0))
          {
            const char* varname = PyUnicode_AsUTF8(tpl0);
            posvd = K_ARRAY::isNamePresent(varname, varStringD);
            posvr = K_ARRAY::isNamePresent(varname, varStringR);
            if (posvr != -1)
            {
              vectOfRcvFields[nfoundvar]= fieldsR[posvr];
              vectOfDnrFields[nfoundvar]= fieldsD[posvd];
              nfoundvar += 1;
            }
          }
#endif
          else
          {
            PyErr_Warn(PyExc_Warning, "_setIBCTransfers: variable must be a string. Skipped.");
          }
        }
      }
    }

    delete [] varStringR; delete [] varStringD; delete [] eltTypeR; delete [] eltTypeD;
  }

  // les variables a transferes sont compactes: on recuperes uniquement la premiere et la taille
  else
  {
    // std::cout << "COUCOU 5788" << std::endl;
# include "getfromzonecompact.h"
    for (E_Int eq = 0; eq < nvars; eq++)
    {
      vectOfRcvFields[eq]= iptroR + eq*ndimdxR;
      vectOfDnrFields[eq]= iptroD + eq*ndimdxD;
    }
  }

  // # include "commonInterpTransfers_indirect.h"
  // std::cout << "COUCOU 5796" << std::endl;

# pragma omp parallel default(shared)
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

    E_Float cvgam = cv*(gamma-1.);

    //0 : Density
    //1 : Temperature
    //2 / 3 / 4 : gradxDensity / gradyDensity / gradzDensity
    //5 / 6 / 7 : gradxTemperature / gradyTemperature / gradzTemperature
    //8 / 9 / 10 : gradxVelocityX / ...
    //11 / 12 / 13 : gradxVelocityY / ...
    //14 / 15 / 16 : gradxVelocityZ / ...
    //17 / 18 / 19 : VelocityX / ...

    // #ifdef _OPENMP4
    //    #pragma omp simd
    // #endif

    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      E_Int indR = rcvPts[noind+ideb];

      // pressure[noind+ideb] = vectOfRcvFields[0][indR]*vectOfRcvFields[1][indR]*cvgam;

      //gradxP[noind+ideb] = ((vectOfRcvFields[1][indR]*vectOfRcvFields[2][indR]+vectOfRcvFields[0][indR]*vectOfRcvFields[5][indR])*cvgam);
      //gradyP[noind+ideb] = ((vectOfRcvFields[1][indR]*vectOfRcvFields[3][indR]+vectOfRcvFields[0][indR]*vectOfRcvFields[6][indR])*cvgam);
      //gradzP[noind+ideb] = ((vectOfRcvFields[1][indR]*vectOfRcvFields[4][indR]+vectOfRcvFields[0][indR]*vectOfRcvFields[7][indR])*cvgam);
     
      gradxP[noind+ideb] = ((vectOfRcvFields[1][indR]*vectOfRcvFields[2][indR]+vectOfRcvFields[0][indR]*vectOfRcvFields[5][indR])*cvgam)/alpha + gradxP[noind+ideb]*(alpha-1.)/alpha;
      gradyP[noind+ideb] = ((vectOfRcvFields[1][indR]*vectOfRcvFields[3][indR]+vectOfRcvFields[0][indR]*vectOfRcvFields[6][indR])*cvgam)/alpha + gradyP[noind+ideb]*(alpha-1.)/alpha;
      gradzP[noind+ideb] = ((vectOfRcvFields[1][indR]*vectOfRcvFields[4][indR]+vectOfRcvFields[0][indR]*vectOfRcvFields[7][indR])*cvgam)/alpha + gradzP[noind+ideb]*(alpha-1.)/alpha;

      // gradxU[noind+ideb] = (vectOfRcvFields[8][indR]);  gradyU[noind+ideb] = (vectOfRcvFields[9][indR]);  gradzU[noind+ideb] = (vectOfRcvFields[10][indR]);
      // gradxV[noind+ideb] = (vectOfRcvFields[11][indR]); gradyV[noind+ideb] = (vectOfRcvFields[12][indR]); gradzV[noind+ideb] = (vectOfRcvFields[13][indR]);
      // gradxW[noind+ideb] = (vectOfRcvFields[14][indR]); gradyW[noind+ideb] = (vectOfRcvFields[15][indR]); gradzW[noind+ideb] = (vectOfRcvFields[16][indR]);
     
    }

  } // Fin zone // omp

  // # include "commonInterpTransfers_indirect.h"
  // sortie
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);

  RELEASESHAREDN(pyIndRcv, rcvPtsI);
  RELEASESHAREDN(pyIndDonor, donorPtsI);
  RELEASESHAREDN(pyArrayTypes, typesI);
  RELEASESHAREDN(pyArrayCoefs, donorCoefsF);
  // BLOCKRELEASEMEM;

  RELEASESHAREDN(pyArrayPressure, pressF);
  RELEASESHAREDN(pyArrayU, UF);
  RELEASESHAREDN(pyArrayV, VF);
  RELEASESHAREDN(pyArrayW, WF);
  RELEASESHAREDN(pyArrayGradxP, gradxPressF);
  RELEASESHAREDN(pyArrayGradyP, gradyPressF);
  RELEASESHAREDN(pyArrayGradzP, gradzPressF);
  RELEASESHAREDN(pyArrayGradxU, gradxVelocityXF);
  RELEASESHAREDN(pyArrayGradyU, gradyVelocityXF);
  RELEASESHAREDN(pyArrayGradzU, gradzVelocityXF);
  RELEASESHAREDN(pyArrayGradxV, gradxVelocityYF);
  RELEASESHAREDN(pyArrayGradyV, gradyVelocityYF);
  RELEASESHAREDN(pyArrayGradzV, gradzVelocityYF);
  RELEASESHAREDN(pyArrayGradxW, gradxVelocityZF);
  RELEASESHAREDN(pyArrayGradyW, gradyVelocityZF);
  RELEASESHAREDN(pyArrayGradzW, gradzVelocityZF);
  // BLOCKRELEASEMEM2;

  Py_INCREF(Py_None);
  return Py_None;
}


//=============================================================================
// Get the values at target points in RCV and put them in the tc (local - as in the copy of the tc locations)
//=============================================================================
PyObject* K_CONNECTOR::_WM_getVal2tc(PyObject* self, PyObject* args)
{
  PyObject *zoneR;
  PyObject *pyVariables;
  PyObject *pyIndRcv;
  PyObject *pyArrayDensWM, *pyArrayVelXWM, *pyArrayVelYWM, *pyArrayVelZWM, *pyArrayTempWM, *pyArraySaNuWM;
  E_Int     loc,nvars;
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;

  if (!PYPARSETUPLE_(args,
                    OOOO_ OOOO_ O_ II_ SSS_,
                    &zoneR, &pyVariables, &pyIndRcv, 
                    &pyArrayDensWM, &pyArrayVelXWM, &pyArrayVelYWM, &pyArrayVelZWM, &pyArrayTempWM, &pyArraySaNuWM,
                    &loc,&nvars,
                    &GridCoordinates,  &FlowSolutionNodes, &FlowSolutionCenters)){
    return NULL;
  }

  vector<PyArrayObject*> hook;

  // recupere les champs du donneur (nodes)
  E_Int imdjmd, imd, jmd, kmd, ndimdxR, meshtype;
  E_Float* iptroR;

  /*--------------------------------------*/
  /* Extraction des indices des receveurs */
  /*--------------------------------------*/
  FldArrayI* rcvPtsI;
  K_NUMPY::getFromNumpyArray(pyIndRcv, rcvPtsI);
  E_Int* rcvPts  = rcvPtsI->begin();
  E_Int nbRcvPts = rcvPtsI->getSize();

  FldArrayF* dens_F; FldArrayF* velx_F; FldArrayF* vely_F;
  FldArrayF* velz_F; FldArrayF* temp_F; FldArrayF* sanu_F;
  E_Int okdens = K_NUMPY::getFromNumpyArray(pyArrayDensWM, dens_F);
  E_Int okvelx = K_NUMPY::getFromNumpyArray(pyArrayVelXWM, velx_F);
  E_Int okvely = K_NUMPY::getFromNumpyArray(pyArrayVelYWM, vely_F);
  E_Int okvelz = K_NUMPY::getFromNumpyArray(pyArrayVelZWM, velz_F);
  E_Int oktemp = K_NUMPY::getFromNumpyArray(pyArrayTempWM, temp_F);
  E_Int oksanu = K_NUMPY::getFromNumpyArray(pyArraySaNuWM, sanu_F);
  E_Float* dens = dens_F->begin();
  E_Float* velx = velx_F->begin();
  E_Float* vely = vely_F->begin();
  E_Float* velz = velz_F->begin();
  E_Float* temp = temp_F->begin();
  E_Float* sanu = sanu_F->begin();

  vector<E_Float*> vectOfRcvFields(nvars);

  // les variables a transferes sont compactes: on recuperes uniquement la premiere et la taille
  //##############################
  PyObject* solR;
  PyObject* t;
  char* type; E_Int s, s0, s1;  E_Int* d;

  PyObject* tpl0= PyList_GetItem(pyVariables, 0);
  char* varname = NULL;
  if (PyString_Check(tpl0)) varname = PyString_AsString(tpl0);
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(tpl0)) varname = (char*)PyUnicode_AsUTF8(tpl0);
#endif

  if  (loc==0) { solR = K_PYTREE::getNodeFromName1(zoneR , "FlowSolution"); }
  else  { solR = K_PYTREE::getNodeFromName1(zoneR , "FlowSolution#Centers"); }
  t = K_PYTREE::getNodeFromName1(solR, "Density_WM");
  iptroR = K_PYTREE::getValueAF(t, hook);

  // get type
  t =  K_PYTREE::getNodeFromName1(zoneR, "ZoneType");
  type =  K_PYTREE::getValueS(t, s, hook);
  // get dims zone receveuse
  d  =  K_PYTREE::getValueAI(zoneR, s0, s1, hook);

  if  (K_STRING::cmp(type, s, "Structured") == 0)
  {
    E_Int shift = 0; if(loc == 1) shift = 3;
    if (s0 == 1) { ndimdxR= d[0+shift]; }
    else if (s0 == 2) { ndimdxR= d[0+shift]*d[1+shift]; } 
    else if (s0 == 3) { ndimdxR= d[0+shift]*d[1+shift]*d[2+shift]; } 
  }
  else
  {
    // non structure
    ndimdxR= d[0]* d[1]; // npoint, nelements
  }
  //##############################

  for (E_Int eq = 0; eq < nvars; eq++)
  {
    vectOfRcvFields[eq]= iptroR + eq*ndimdxR;
  }

# pragma omp parallel default(shared)
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
    if (ithread <= r){ ideb = (ithread-1)*(chunk+1); ifin = ideb + (chunk+1); }
    else { ideb = (chunk+1)*r+(ithread-r-1)*chunk; ifin = ideb + chunk; }

#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      E_Int indR = rcvPts[noind+ideb];
      dens[noind+ideb] = vectOfRcvFields[0][indR];
      velx[noind+ideb] = vectOfRcvFields[1][indR];
      vely[noind+ideb] = vectOfRcvFields[2][indR];
      velz[noind+ideb] = vectOfRcvFields[3][indR];
      temp[noind+ideb] = vectOfRcvFields[4][indR];
    }
    if (nvars==6)
    {
      for (E_Int noind = 0; noind < ifin-ideb; noind++)
      {
        E_Int indR = rcvPts[noind+ideb];
        sanu[noind+ideb] = vectOfRcvFields[5][indR];
      }
    }
  } // Fin zone // omp
  // sortie
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);

  RELEASESHAREDN(pyIndRcv     , rcvPtsI  );
  RELEASESHAREDN(pyArrayDensWM, dens_F );
  RELEASESHAREDN(pyArrayVelXWM, velx_F );
  RELEASESHAREDN(pyArrayVelYWM, vely_F );
  RELEASESHAREDN(pyArrayVelZWM, velz_F );
  RELEASESHAREDN(pyArrayTempWM, temp_F );
  RELEASESHAREDN(pyArraySaNuWM, sanu_F );

  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
// Set the value at the target points in the tc (local) into the real tc (as in the tree itself)
//=============================================================================
PyObject* K_CONNECTOR::_WM_setVal2tc(PyObject* self, PyObject* args)
{
  PyObject *pyArraydens_new, *pyArrayvelx_new, *pyArrayvely_new;
  PyObject *pyArrayvelz_new, *pyArraytemp_new, *pyArraysanu_new;
  PyObject *pyArraydens    , *pyArrayvelx    , *pyArrayvely;
  PyObject *pyArrayvelz    , *pyArraytemp    , *pyArraysanu;

  if (!PYPARSETUPLE_(args, OOOO_ OOOO_ OOOO_,
                    &pyArraydens_new, &pyArrayvelx_new, &pyArrayvely_new,
                    &pyArrayvelz_new, &pyArraytemp_new, &pyArraysanu_new,
                    &pyArraydens    , &pyArrayvelx    , &pyArrayvely    ,
                    &pyArrayvelz    , &pyArraytemp    , &pyArraysanu    ))
    {
      return NULL;
    }

  FldArrayF* densF; FldArrayF* velxF; FldArrayF* velyF;
  FldArrayF* velzF; FldArrayF* tempF; FldArrayF* sanuF;
  E_Int okdens = K_NUMPY::getFromNumpyArray(pyArraydens, densF);
  E_Int okvelx = K_NUMPY::getFromNumpyArray(pyArrayvelx, velxF);
  E_Int okvely = K_NUMPY::getFromNumpyArray(pyArrayvely, velyF);
  E_Int okvelz = K_NUMPY::getFromNumpyArray(pyArrayvelz, velzF);
  E_Int oktemp = K_NUMPY::getFromNumpyArray(pyArraytemp, tempF);
  E_Int oksanu = K_NUMPY::getFromNumpyArray(pyArraysanu, sanuF);
  E_Float* dens = densF->begin();
  E_Float* velx = velxF->begin();
  E_Float* vely = velyF->begin();
  E_Float* velz = velzF->begin();
  E_Float* temp = tempF->begin();
  E_Float* sanu = sanuF->begin();

  FldArrayF* densF_new; FldArrayF* velxF_new; FldArrayF* velyF_new;
  FldArrayF* velzF_new; FldArrayF* tempF_new; FldArrayF* sanuF_new;
  okdens = K_NUMPY::getFromNumpyArray(pyArraydens_new, densF_new);
  okvelx = K_NUMPY::getFromNumpyArray(pyArrayvelx_new, velxF_new);
  okvely = K_NUMPY::getFromNumpyArray(pyArrayvely_new, velyF_new);
  okvelz = K_NUMPY::getFromNumpyArray(pyArrayvelz_new, velzF_new);
  oktemp = K_NUMPY::getFromNumpyArray(pyArraytemp_new, tempF_new);
  oksanu = K_NUMPY::getFromNumpyArray(pyArraysanu_new, sanuF_new);
  E_Float* dens_new = densF_new->begin();
  E_Float* velx_new = velxF_new->begin();
  E_Float* vely_new = velyF_new->begin();
  E_Float* velz_new = velzF_new->begin();
  E_Float* temp_new = tempF_new->begin();
  E_Float* sanu_new = sanuF_new->begin();

  E_Int nbRcvPts = densF_new->getSize();

# pragma omp parallel default(shared)
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

#ifdef _OPENMP4
#pragma omp simd
#endif
    for (E_Int noind = 0; noind < ifin-ideb; noind++)
    {
      dens[noind+ideb] = dens_new[noind+ideb];
      velx[noind+ideb] = velx_new[noind+ideb];
      vely[noind+ideb] = vely_new[noind+ideb];
      velz[noind+ideb] = velz_new[noind+ideb];
      temp[noind+ideb] = temp_new[noind+ideb];
      sanu[noind+ideb] = sanu_new[noind+ideb];
    }
  } // Fin zone // omp
  // sortie

  RELEASESHAREDN(pyArraydens, densF);
  RELEASESHAREDN(pyArrayvelx, velxF);
  RELEASESHAREDN(pyArrayvely, velyF);
  RELEASESHAREDN(pyArrayvelz, velzF);
  RELEASESHAREDN(pyArraytemp, tempF);
  RELEASESHAREDN(pyArraysanu, sanuF);

  RELEASESHAREDN(pyArraydens_new, densF_new);
  RELEASESHAREDN(pyArrayvelx_new, velxF_new);
  RELEASESHAREDN(pyArrayvely_new, velyF_new);
  RELEASESHAREDN(pyArrayvelz_new, velzF_new);
  RELEASESHAREDN(pyArraytemp_new, tempF_new);
  RELEASESHAREDN(pyArraysanu_new, sanuF_new);

  Py_INCREF(Py_None);
  return Py_None;
}

