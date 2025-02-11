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
#ifndef _POST_STREAM_H_
#define _POST_STREAM_H_

# include "../post.h"
#define FldArrayF K_FLD::FldArrayF
#define FldArrayI K_FLD::FldArrayI

/* Define a tracer */
struct tracer
{
    E_Float x,y,z;
    tracer* left;
    tracer* right;
};

using namespace std;

namespace K_POST
{
/* Dans toutes les routines ci-dessous :
   IN: listOfStructInterpData: liste des interpData des grilles structurees
   IN: listOfStructFields associes
   IN: listOfStructVelocities associes
   IN: nis,njs,nks: dimensions des grilles structurees
   IN: posxs, posys, poszs, poscs: positions de x,y,z,cellN dans listOfStructFields
   IN: listOfUnstrInterpData: liste des interpData des grilles TETRA
   IN: listOfUnstrFields associes
   IN: listOfUnstrVelocities associes
   IN: connectu: connectivite des grilles TETRA
   IN: posxu, posyu, poszu, poscu: positions de x,y,z,cellN dans listOfUnstrFields
   IN: interpType : type d'interpolation: necessaire uniquement en structure */

// 1 - methodes communes a tous 
/* Calcul du pas de temps initial 
   IN: noblk: no du bloc d'interpolation
   IN: indi: indices des pts d'interpolation (+indi[0] : info structure/non structure)
   OUT: dt0: pas de temps initial */
  void compInitialStep(E_Int noblk, E_Int type, FldArrayI& indi, 
                       vector<E_Int>& nis, vector<E_Int>& njs, 
                       vector<E_Int>& nks, vector<E_Int>& posxs, 
                       vector<E_Int>& posys, vector<E_Int>& poszs,
                       vector<FldArrayF*>& structFields, 
                       vector<FldArrayF*>& structVelocities,
                       vector<E_Int>& posxu, vector<E_Int>& posyu, 
                       vector<E_Int>& poszu,  
                       vector<FldArrayF*>& unstrFields,
                       vector<FldArrayI*>& connectu,
                       vector<FldArrayF*>& unstrVelocities,
                       E_Float& dt0);

/* 
   Calcul d un point de la ligne de courant
   IN : nopt : numero du pt  dans le tableau streamPt
   IN : (x0,y0,z0) : premier pt a interpoler
   IN : noblk : numero du bloc d interpolation
   IN : indi : indices des sommets de la cellule d interpolation
   IN : cf : coefs d interpolation 
   IN/OUT : streamPt : tableau des pts de la ligne de courant */
  void 
  compStreamPtFields(E_Int nopt, E_Float x0, E_Float y0, E_Float z0, 
                     E_Int noblk, E_Int type, FldArrayI& indi, FldArrayF& cf,
                     vector<E_Int>& nis, vector<E_Int>& njs, 
                     vector<E_Int>& nks, vector<E_Int>& posxs, 
                     vector<E_Int>& posys, vector<E_Int>& poszs,
                     vector<FldArrayF*>& structFields,
                     vector<E_Int>& posxu, vector<E_Int>& posyu, 
                     vector<E_Int>& poszu,  
                     vector<FldArrayF*>& unstrFields,
                     vector<FldArrayI*>& connectu,
                     FldArrayF& streamPt,
                     K_INTERP::InterpData::InterpolationType interpType);

/* Initialisation de la ligne de courant : vitesse et dt calcules
   IN : (xp,yp,xp) : coordonnees du pt X0
   IN : connectSurf, (xSurf, ySurf, zSurf), sizeSurf : connection, coordonnees et taille de l'array
        contenant les surfaces 2D
   OUT : up, vp, wp : vitesse du pt X0
   OUT : dt : pas de temps calcule
   OUT : indi : indices de la cellule d interpolation de X0
   OUT : cf : coefs d interpolation correspondants
   OUT : streamPts : mise a jour pour le point X0
   Retourne le numero du bloc d interpolation (demarre a 1) pour le nouveau pt, 
   et 0 si pas interpolable */
  short 
  initStreamLine(E_Float xp, E_Float yp, E_Float zp,
                 vector<K_INTERP::InterpData*>& listOfStructInterpData, 
                 vector<FldArrayF*>& listOfStructFields,
                 vector<FldArrayF*>& listOfStructVelocities,
                 vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
                 vector<E_Int>& posxs, vector<E_Int>& posys, 
                 vector<E_Int>& poszs, vector<E_Int>& poscs,
                 vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
                 vector<FldArrayF*>& listOfUnstrFields,
                 vector<FldArrayF*>& listOfUnstrVelocities,
                 vector<FldArrayI*>& connectu,
                 vector<E_Int>& posxu, vector<E_Int>& posyu, 
                 vector<E_Int>& poszu, vector<E_Int>& poscu,
                 FldArrayI& connectSurf, 
                 E_Float* xSurf, E_Float* ySurf, E_Float* zSurf, 
                 E_Int sizeSurf, 
                 E_Float& up, E_Float& vp, E_Float& wp, E_Float& dt,
                 FldArrayI& indi, FldArrayF& cf, FldArrayF& streamPts, 
                 K_INTERP::InterpData::InterpolationType interpType);

//------------------------------------------------------------------
//2 - methodes specifiques a streamLine

/* calcul des points constituant la ligne de courant 
   IN: (x0,y0,z0): point initial de la ligne
   IN: connectSurf, (xSurf, ySurf, zSurf), sizeSurf: connection, coordonnees et taille de l'array
        contenant les surfaces 2D
   OUT: streamPoints: tableau des points de la ligne de courant: 
   coordonnees + champs initiaux */
  E_Int computeStreamLineElts(
    E_Float x0, E_Float y0, E_Float z0, 
    vector<K_INTERP::InterpData*>& listOfStructInterpData, 
    vector<FldArrayF*>& listOfStructFields,
    vector<FldArrayF*>& listOfStructVelocities,
    vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks, 
    vector<E_Int>& posxs, vector<E_Int>& posys, vector<E_Int>& poszs, 
    vector<E_Int>& poscs,
    vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
    vector<FldArrayF*>& listOfUnstrFields,
    vector<FldArrayF*>& listOfUnstrVelocities,
    vector<FldArrayI*>& connectu,
    vector<E_Int>& posxu, vector<E_Int>& posyu, vector<E_Int>& poszu, 
    vector<E_Int>& poscu,
    FldArrayI& connectSurf, E_Float* xSurf, E_Float* ySurf, E_Float* zSurf, E_Int sizeSurf, 
    FldArrayF& streamPts);

/* Calcule les coordonnees du pt X(n+1), ses cf et indi et la vitesse U(n+1)
   Si necessaire, le pas dt est modifie 
   IN : connectSurf, (xSurf, ySurf, zSurf), sizeSurf : connection, coordonnees et taille de l'array
        contenant les surfaces 2D
   IN/OUT : xp, yp, zp : coordonnees du pt X(n+1) 
   IN/OUT : up, vp, wp : son vecteur vitesse aussi calcule ici
   IN/OUT : cfp : ses nvx coefs d interpolation 
   IN/OUT : indip : nvelle cellule d interpolation 
   IN/OUT : dt : pas de temps eventuellement modifie
   Retourne le numero du bloc d interpolation (demarre a 1) pour le nouveau pt, 
   et 0 si pas interpolable */
  short updateStreamLinePoint(
    E_Float& xp, E_Float& yp, E_Float& zp,
    E_Float& up, E_Float& vp, E_Float& wp,
    E_Int& type, FldArrayI& indip, FldArrayF& cfp,
    E_Float& dt, 
    vector<K_INTERP::InterpData*>& listOfStructInterpData,
    vector<FldArrayF*>& listOfStructFields,
    vector<FldArrayF*>& listOfStructVelocities,
    vector<E_Int>& nis, vector<E_Int>& njs,
    vector<E_Int>& nks, vector<E_Int>& posxs, 
    vector<E_Int>& posys, vector<E_Int>& poszs, 
    vector<E_Int>& poscs, 
    vector<K_INTERP::InterpData*>& listOfUnstrInterpData,
    vector<FldArrayF*>& listOfUnstrFields,
    vector<FldArrayF*>& listOfUnstrVelocities,
    vector<FldArrayI*>& connectu,
    vector<E_Int>& posxu, vector<E_Int>& posyu, 
    vector<E_Int>& poszu, vector<E_Int>& poscu, 
    FldArrayI& connectSurf, E_Float* xSurf, E_Float* ySurf, E_Float* zSurf, E_Int sizeSurf, 
    K_INTERP::InterpData::InterpolationType interpType);

/* Calcul des 4 coefs de Runge-Kutta4 
   IN : xp,yp,zp : coordonnees du pt X(n)
   IN : up,vp,wp : vitesse du pt X(n)
   IN : connectSurf, (xSurf, ySurf, zSurf), sizeSurf : connection, coordonnees et taille de l'array
        contenant les surfaces 2D
   OUT : dt : pas de temps modifie eventuellement
   OUT : xn, yn, zn : coordonnees du point X(n+1)
   Retourne 0 si echec  */
  short compRungeKutta4(E_Float xp, E_Float yp, E_Float zp,
                        E_Float up, E_Float vp, E_Float wp, 
                        E_Float& dt, E_Float& xn, E_Float& yn, E_Float& zn,
                        vector<K_INTERP::InterpData*>& listOfStructInterpData,
                        vector<FldArrayF*>& listOfStructFields,
                        vector<FldArrayF*>& listOfStructVelocities,
                        vector<E_Int>& nis,vector<E_Int>& njs,
                        vector<E_Int>& nks, vector<E_Int>& posxs, 
                        vector<E_Int>& posys, vector<E_Int>& poszs, 
                        vector<E_Int>& poscs, 
                        vector<K_INTERP::InterpData*>& listOfUnstrInterpData,
                        vector<FldArrayF*>& listOfUnstrFields,
                        vector<FldArrayF*>& listOfUnstrVelocities,
                        vector<FldArrayI*>& connectu,
                        vector<E_Int>& posxu, vector<E_Int>& posyu, 
                        vector<E_Int>& poszu, vector<E_Int>& poscu, 
                        FldArrayI& connectSurf, E_Float* xSurf, E_Float* ySurf, E_Float* zSurf, E_Int sizeSurf, 
                        K_INTERP::InterpData::InterpolationType interpType);

//------------------------------------------------------------------
// 3 - methodes specifiques a streamRibbon

/* calcul des points constituant la ligne de courant 
   IN : (xp,yp,zp) : point initial de la ligne
   IN : (nxp,nyp,nzp) : normale donnant la taille/orientation du ruban 
   OUT : streamPoints1 : tableau des points de la ligne de courant : 
   coordonnees + champs initiaux 
   OUT : streamPts2 : pts en face de streamPts1 ds le ruban*/
  E_Int computeStreamRibbonElts(
    E_Float xp, E_Float yp, E_Float zp, 
    E_Float nxp, E_Float nyp, E_Float nzp,
    vector<K_INTERP::InterpData*>& listOfStructInterpData, 
    vector<FldArrayF*>& listOfStructFields,
    vector<FldArrayF*>& listOfStructVelocities,
    vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks, 
    vector<E_Int>& posxs, vector<E_Int>& posys, vector<E_Int>& poszs, 
    vector<E_Int>& poscs,
    vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
    vector<FldArrayF*>& listOfUnstrFields,
    vector<FldArrayF*>& listOfUnstrVelocities,
    vector<FldArrayI*>& connectu,
    vector<E_Int>& posxu, vector<E_Int>& posyu, vector<E_Int>& poszu, 
    vector<E_Int>& poscu,
    FldArrayF& streamPts1, FldArrayF& streamPts2);

/* Initialisation du ruban de courant : vitesse et dt calcules
   IN : (xp,yp,xp) : coordonnees du pt X0
   OUT : up, vp, wp : vitesse du pt X0
   OUT : nxp, nyp, nzp : normale au pt X0 de (up,vp,wp)
   OUT : indi : indices de la cellule d interpolation de X0
   OUT : cf : coefs d interpolation correspondants
   OUT : dt : pas de temps calcule
   OUT : streamPts1 : mise a jour pour le point X0
   OUT : streamPts2 : mise a jour pour le point X0'
   Retourne le numero du bloc d interpolation (demarre a 1) pour le nouveau
   pt, et 0 si pas interpolable */
  short initStreamRibbon(
    E_Float xp, E_Float yp, E_Float zp,
    vector<K_INTERP::InterpData*>& listOfStructInterpData, 
    vector<FldArrayF*>& listOfStructFields,
    vector<FldArrayF*>& listOfStructVelocities,
    vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
    vector<E_Int>& posxs, vector<E_Int>& posys, 
    vector<E_Int>& poszs, vector<E_Int>& poscs,
    vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
    vector<FldArrayF*>& listOfUnstrFields,
    vector<FldArrayF*>& listOfUnstrVelocities,
    vector<FldArrayI*>& connectu,
    vector<E_Int>& posxu, vector<E_Int>& posyu, 
    vector<E_Int>& poszu, vector<E_Int>& poscu,
    E_Float& up, E_Float& vp, E_Float& wp, E_Float& dt,
    E_Float& nxp, E_Float& nyp, E_Float& nzp,
    FldArrayI& indi, FldArrayF& cf, 
    FldArrayF& streamPts1, FldArrayF& streamPts2, 
    K_INTERP::InterpData::InterpolationType interpType);

/* Initialisation de la normale a la ligne de courant
   IN : (x0,y0,z0) : point de depart P0 
   IN : (u0,v0,w0) : vecteur vitesse associe
   OUT : (nx0,ny0,nz0) : normale en X0 a U0. Si Np n est pas orthogonale
   a U0, alors elle est redressee. Si echec : retourne -1 */
  short initNormalToStreamLine(
    const E_Float x0, const E_Float y0, const E_Float z0,
    const E_Float u0, const E_Float v0, const E_Float w0,
    E_Float& nx0, E_Float& ny0, E_Float& nz0);

/* Calcule le second point XP' a partir de XP 
   Retourne 0 si pas interpolable, le no du bloc d interpolation sinon
   IN : nopt : numero du pt XP' a calculer dans streamPt2
   IN : xp, yp, zp : coordonnees du pt XP
   IN : nxp, nyp, nzp : normale a UP au pt XP
   IN/OUT : streamPts2 mis à jour au pt nopt */
  short compSecondPoint(
    const E_Int nopt,
    const E_Float xp, const E_Float yp, const E_Float zp,
    const E_Float nxp, const E_Float nyp, const E_Float nzp,
    vector<K_INTERP::InterpData*>& listOfStructInterpData, 
    vector<FldArrayF*>& listOfStructFields,
    vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
    vector<E_Int>& posxs, vector<E_Int>& posys, 
    vector<E_Int>& poszs, vector<E_Int>& poscs,
    vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
    vector<FldArrayF*>& listOfUnstrFields,
    vector<FldArrayI*>& connectu,
    vector<E_Int>& posxu, vector<E_Int>& posyu, 
    vector<E_Int>& poszu, vector<E_Int>& poscu,
    FldArrayF& streamPts2,  
    K_INTERP::InterpData::InterpolationType interpType);

/* Calcule les coordonnees du pt X(n+1), ses cf et indi et la vitesse U(n+1)
   Si necessaire, le pas dt est modifie 
   Retourne 1 si le nouveau pt est interpolable, 0 sinon.
   IN : niter=n+1 : numero des pts X(n+1) et X(n+1)' a calculer dans streamPt1 et streamPt2
   IN : noblkp : numero du bloc d interpolation
   IN : xp, yp, zp : coordonnees du pt X(n+1)
   IN : up, vp, wp : vitesse du pt U(n+1)
   IN : nxp, nyp, nzp : normale a U(n+1) au pt X(n+1)
   IN : thetap : angle entre U(n-1) et U(n) 
   IN : indip, cfp : donnees d interpolation pour le pt X(p)
   IN/OUT : dt modifie si necessaire
   IN/OUT : streamPts1 mis à jour au pt X(n+1)  
   IN/OUT : streamPts2 mis à jour au pt X'(n+1) 
   OUT: noblkn: le numero du bloc d interpolation (demarre a 1) 
   pour le nouveau pt, et 0 si pas interpolable 
   OUT: typen: nouveau type d'interpolation
*/
  short updateStreamRibbonPoints(
    E_Int niter, E_Int noblkp, E_Int typep,
    E_Float& xp, E_Float& yp, E_Float& zp,
    E_Float& up, E_Float& vp, E_Float& wp,
    E_Float& nxp, E_Float& nyp, E_Float& nzp,
    E_Float& thetap,
    FldArrayI& indip, FldArrayF& cfp, E_Float& dt, 
    FldArrayF& streamPts1, FldArrayF& streamPts2,
    vector<K_INTERP::InterpData*>& listOfStructInterpData, 
    vector<FldArrayF*>& listOfStructFields,
    vector<FldArrayF*>& listOfStructVelocities,
    vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
    vector<E_Int>& posxs, vector<E_Int>& posys, 
    vector<E_Int>& poszs, vector<E_Int>& poscs,
    vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
    vector<FldArrayF*>& listOfUnstrFields,
    vector<FldArrayF*>& listOfUnstrVelocities,
    vector<FldArrayI*>& connectu,
    vector<E_Int>& posxu, vector<E_Int>& posyu, 
    vector<E_Int>& poszu, vector<E_Int>& poscu,
    K_INTERP::InterpData::InterpolationType interpType,
    E_Int& noblkn, E_Int& typen);

/* Calcul des 4 coefs de Runge-Kutta4 
   IN : noblkp : no du bloc d interpolation de Xp dans listOfInterpData
   IN : xp,yp,zp : coordonnees du pt X(p)
   IN : thetap : angle 
   IN : up,vp,wp : vitesse du pt X(p)
   IN : indip, cfp : infos d interpolation pour le pt X(p)
   OUT : dt : pas de temps modifie eventuellement
   OUT : xn, yn, zn : coordonnees du point X(p+1)
   OUT : thetan : theta du point X(p+1)
   Retourne 0 si echec  */
  short compRungeKutta4ForRibbon(const E_Int noblkp, E_Int typep,
                                 E_Float xp, E_Float yp, E_Float zp, E_Float thetap,
                                 E_Float up, E_Float vp, E_Float wp, 
                                 FldArrayI& indip, FldArrayF& cfp,
                                 E_Float& dt, E_Float& xn, E_Float& yn, E_Float& zn, E_Float& thetan,
                                 vector<K_INTERP::InterpData*>& listOfStructInterpData, 
                                 vector<FldArrayF*>& listOfStructFields,
                                 vector<FldArrayF*>& listOfStructVelocities,
                                 vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
                                 vector<E_Int>& posxs, vector<E_Int>& posys, 
                                 vector<E_Int>& poszs, vector<E_Int>& poscs,
                                 vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
                                 vector<FldArrayF*>& listOfUnstrFields,
                                 vector<FldArrayF*>& listOfUnstrVelocities,
                                 vector<FldArrayI*>& connectu,
                                 vector<E_Int>& posxu, vector<E_Int>& posyu, 
                                 vector<E_Int>& poszu, vector<E_Int>& poscu,
                                 K_INTERP::InterpData::InterpolationType interpType);

/* Calcul du coefficient ki du sous pas i de Runge Kutta 4 pour l equation 
   en theta : dtheta/dt = 0.5*(rotv . v/normv) 
   IN : noblk : numero du bloc d interpolation de XP
   IN : up,vp,wp : vitesse au pt XP
   IN : indi : indices des pts de la cellule d interpolation
   OUT : kcoef : coef ki du sous pas i Runge-Kutta pour l edo en dtheta */
  short getThetaRKCoef(const E_Int noblk, E_Int type,
                       const E_Float up, const E_Float vp, const E_Float wp,
                       const FldArrayI& indi, 
                       vector<FldArrayF*>& listOfStructVelocities,
                       vector<FldArrayF*>& listOfStructFields,  
                       vector<E_Int>& nis, vector<E_Int>& njs,
                       vector<E_Int>& nks, vector<E_Int>& posxs, 
                       vector<E_Int>& posys, vector<E_Int>& poszs, 
                       vector<FldArrayF*>& listOfUnstrVelocities,
                       vector<FldArrayF*>& listOfUnstrFields,
                       vector<FldArrayI*>& connectu, vector<E_Int>& posxu, 
                       vector<E_Int>& posyu, vector<E_Int>& poszu, 
                       E_Float& kcoef);

/* Calcul de la connectivite */
  void buildConnectivity(FldArrayF& field1, FldArrayF& field2,
                         FldArrayF& field, FldArrayI& cn);

/* Calcul de la normale au pt Xn, connaissant ses coordonnees, son theta,
   et la normale a l instant precedent, en Xp */
  void compNormalToStreamLine(const E_Float xn, const E_Float yn, 
                              const E_Float zn, const E_Float thetan, 
                              const E_Float nxp, const E_Float nyp, 
                              const E_Float nzp,
                              E_Float& nxn, E_Float& nyn, E_Float& nzn);

/* Determine la liste des arrays structures, dont un des champs est la 
   vitesse. Retourne la vitesse sur ces arrays
   IN: signe: parcours de la streamLine ds 1 sens ou l autre +/-1
   IN: listes 'In' des donnees sur les arrays structures
   OUT: listes 'Out': seuls sont conserves les arrays ayant les  variables du vecteur
   OUT: vectorf: liste correspondante des champs des variables du vecteur*/
  E_Int extractVectorFromStructArrays(
    E_Float signe,
    vector<E_Int>& niIn, vector<E_Int>& njIn, vector<E_Int>& nkIn,
    vector<E_Int>& posxIn, vector<E_Int>& posyIn, vector<E_Int>& poszIn,
    vector<E_Int>& poscIn, vector<char*>& varStringIn, 
    vector<FldArrayF*>& fieldsIn, 
    vector<K_INTERP::InterpData*>& structInterpDataIn,
    vector<E_Int>& niOut, vector<E_Int>& njOut, vector<E_Int>& nkOut,
    vector<E_Int>& posxOut, vector<E_Int>& posyOut, vector<E_Int>& poszOut,
    vector<E_Int>& poscOut, vector<char*>& varStringOut, 
    vector<FldArrayF*>& fieldsOut, 
    vector<K_INTERP::InterpData*>& structInterpDataOut,
    vector<FldArrayF*>& vectorf, vector<char*>& vnames);

/* Determine la liste des arrays non structures TETRA, dont un des champs
   est la vitesse. Retourne la vitesse sur ces arrays
   IN: signe: parcours de la streamLine ds 1 sens ou l autre +/-1 
   IN: listes 'In' des donnees sur les arrays tetra
   OUT: listes 'Out': seuls sont conserves les arrays ayant les  variables du vecteur
   OUT: vectorf: liste correspondante des champs des variables du vecteur */
  E_Int extractVectorFromUnstrArrays(
    E_Float signe, vector<E_Int>& posxIn, vector<E_Int>& posyIn, 
    vector<E_Int>& poszIn, vector<E_Int>& poscIn, 
    vector<char*>& varStringIn, vector<FldArrayF*>& fieldsIn, 
    vector<FldArrayI*>& cntIn, vector<char*>& eltTypeIn,
    vector<K_INTERP::InterpData*>& interpDataIn,
    vector<E_Int>& posxOut, vector<E_Int>& posyOut, 
    vector<E_Int>& poszOut, vector<E_Int>& poscOut, 
    vector<char*>& varStringOut, vector<FldArrayF*>& fieldsOut, 
    vector<FldArrayI*>& cntOut, vector<char*>& eltTypeOut,
    vector<K_INTERP::InterpData*>& interpDataOut,
    vector<FldArrayF*>& vectorf, vector<char*>& vnames);

/* Initialisation de la ligne de courant : vitesse et dt calcules
   IN: (xp,yp,xp): coordonnees du pt X0
   OUT: up, vp, wp: vitesse du pt X0
   OUT: dt: pas de temps calcule
   OUT: indi: indices de la cellule d interpolation de X0
   OUT: cf: coefs d interpolation correspondants
   OUT: streamPts: mise a jour pour le point X0 */
  void 
  initStreamSurf(E_Float xp, E_Float yp, E_Float zp,
                 vector<K_INTERP::InterpData*>& listOfStructInterpData, 
                 vector<FldArrayF*>& listOfStructFields,
                 vector<FldArrayF*>& listOfStructVelocities,
                 vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks,
                 vector<E_Int>& posxs, vector<E_Int>& posys, 
                 vector<E_Int>& poszs, vector<E_Int>& poscs,
                 vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
                 vector<FldArrayF*>& listOfUnstrFields,
                 vector<FldArrayF*>& listOfUnstrVelocities,
                 vector<FldArrayI*>& connectu,
                 vector<E_Int>& posxu, vector<E_Int>& posyu, 
                 vector<E_Int>& poszu, vector<E_Int>& poscu,
                 E_Float& up, E_Float& vp, E_Float& wp, E_Float& dt,
                 FldArrayI& indi, FldArrayF& cf, FldArrayF& streamPts, 
                 K_INTERP::InterpData::InterpolationType interpType);


/* creation du front d'avancement pour le calcul des nappes de courant 
   IN  : npts             : nombre de points constituant le BAR-array
   IN  : xBAR, yBAR, zBAR : coordonnees issus du BAR-array entre par l'utilisateur 
   IN  : cn               : connectivites issues du BAR-array entre par l'utilisateur 
   OUT : front            : front d'avancement correspondant au BAR-array */
void createFront(E_Float* xBAR, E_Float* yBAR, E_Float* zBAR, FldArrayI& cn, 
                 vector<tracer*>& front, E_Int npts);

/* Progression du front 
   IN     : npts             : nombre de points constituant le BAR-array
   IN/OUT : front            : front d'avancement
   IN/OUT : tleft            : extremite gauche du front d'avancement 
   IN/OUT : tright           : extremite droite du front d'avancement  
   OUT    : nt               : nombre de triangles constituant la nappe de courant */
void advanceFront(vector<tracer*> front, tracer* tleft, tracer* tright, E_Int npts, E_Int& nt,
                  FldArrayF* field, FldArrayI* cn,
                  vector<K_INTERP::InterpData*>& structInterpDatas, 
                  vector<FldArrayF*>& structFields, 
                  vector<FldArrayF*>& structVelocities,
                  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks, 
                  vector<E_Int>& posxs, vector<E_Int>& posys, vector<E_Int>& poszs, 
                  vector<E_Int>& poscs,
                  vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
                  vector<FldArrayF*>& listOfUnstrFields,
                  vector<FldArrayF*>& listOfUnstrVelocities,
                  vector<FldArrayI*>& connectu,
                  vector<E_Int>& posxu, vector<E_Int>& posyu, vector<E_Int>& poszu, 
                  vector<E_Int>& poscu);

/* Avance du "ribbon" en allant de gauche a droite
   IN     : npts             : nombre de points constituant le BAR-array
   IN/OUT : t                : traceur du "ribbon" concerne
   OUT    : nt               : nombre de triangles constituant la nappe de courant */
void advanceRibbonLeft(tracer* t, E_Int npts, E_Int& nt,
                       FldArrayF* field, FldArrayI* cn,
                       vector<K_INTERP::InterpData*>& structInterpDatas, 
                       vector<FldArrayF*>& structFields, 
                       vector<FldArrayF*>&structVelocities,
                       vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks, 
                       vector<E_Int>& posxs, vector<E_Int>& posys, vector<E_Int>& poszs, 
                       vector<E_Int>& poscs,
                       vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
                       vector<FldArrayF*>& listOfUnstrFields,
                       vector<FldArrayF*>& listOfUnstrVelocities,
                       vector<FldArrayI*>& connectu,
                       vector<E_Int>& posxu, vector<E_Int>& posyu, vector<E_Int>& poszu, 
                       vector<E_Int>& poscu);

/* Avance du "ribbon" en allant de droite a gauche
   IN     : npts             : nombre de points constituant le BAR-array
   IN/OUT : t                : traceur du "ribbon" concerne
   OUT    : nt               : nombre de triangles constituant la nappe de courant */
void advanceRibbonRight(tracer* t, E_Int npts, E_Int& nt,
                        FldArrayF* field, FldArrayI* cn,
                        vector<K_INTERP::InterpData*>& structInterpDatas, 
                        vector<FldArrayF*>& structFields, 
                        vector<FldArrayF*>&structVelocities,
                        vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks, 
                        vector<E_Int>& posxs, vector<E_Int>& posys, vector<E_Int>& poszs, 
                        vector<E_Int>& poscs,
                        vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
                        vector<FldArrayF*>& listOfUnstrFields,
                        vector<FldArrayF*>& listOfUnstrVelocities,
                        vector<FldArrayI*>& connectu,
                        vector<E_Int>& posxu, vector<E_Int>& posyu, vector<E_Int>& poszu, 
                        vector<E_Int>& poscu);

/* Ecriture d'un triangle (coordonnees et connectivites) de la nappe de courant 
IN     : (f1,f2,f3)       : champ au 3 points du triangle
IN/OUT : nt               : nombre de triangle constituant la nappe de courant
OUT    : field            : tableau de champs pour le triangle
OUT    : cn               : tableau des connectivites pour le triangle
*/
void writeTriangle(FldArrayF& field, FldArrayI& cn, 
                   FldArrayF& f1, FldArrayF& f2, FldArrayF& f3,
                   E_Int& nt);
}

#undef FldArrayF
#undef FldArrayI
#endif
