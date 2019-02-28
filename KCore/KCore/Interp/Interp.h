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
#ifndef _KCORE_INTERP_H_
#define _KCORE_INTERP_H_

# include "Def/DefTypes.h"
# include "Def/DefFunction.h"
# include "Fld/FldArray.h"
# include "Array/Array.h"
# include "Search/KdTree.h"
# include "Fld/ArrayAccessor.h"
# include <vector>
# include <list>
# define FldArrayF K_FLD::FldArrayF
# define FldArrayI K_FLD::FldArrayI
# define FldArrayIS K_FLD::FldArrayIS
# include "Interp/InterpData.h"
# include "Interp/InterpAdt.h"
# include "Interp/InterpCart.h"

namespace K_INTERP
{

/* IN: x,y,z: point a interpoler.
   IN: interpDatas: les interpData des domaines donneurs
   IN: fields: des champs des domaines donneurs
   IN: STRUCTURE: a1=nit, a2=njt, a3=nkt
   IN: NON STRUCTURE: a1=cEV, a4=cEEN ou NULL
   IN: posxt, posyt, poszt, posct: position de x,y,z,c dans les domaines 
   donneurs
   IN: interpType: type d'interpolation
   IN: nature: 0 (pas de cellN=0 dans la cellule d'interpolation)
               1 (pas de cellN=0 ni cellN=2 dans la cellule d'interpolation)
   IN: penalty: 0: les bords ne sont pas penalise dans la recherche
                1: les bords sont penalise (seront pris en dernier)
   OUT: voli: volume de la cellule d'interpolation
   OUT: donorIndices: indices de la cellule d'interpolation
   OUT: donorCoefs: coeff d'interpolation
   OUT: type: type d'interpolation reellement effectuee
        0: echec, 1: coincident, 2: O2CF, 3: O3ABC, 
        4: O2 avec indice de l'element pour la cellule d'interpolation, 
        5: O5ABC 
   OUT: noDonorBlk: no du bloc donneur utilise dans l'interpolation  */
  short getInterpolationCell(
    E_Float x, E_Float y, E_Float z,
    InterpData* interpData,
    FldArrayF* fields,
    void* a1, void* a2, void* a3, void* a4, 
    E_Int posxt, E_Int posyt, E_Int poszt, E_Int posct,
    E_Float& voli, FldArrayI& donorIndices, FldArrayF& donorCoefs,
    E_Int& type, E_Int& noDonorBlk,
    InterpData::InterpolationType interpType=InterpData::O2CF,
    E_Int nature=0, E_Int penalty=0);

  short getExtrapolationCell(
    E_Float x, E_Float y, E_Float z,
    InterpData* interpData,
    FldArrayF* fields,
    void* a1, void* a2, void* a3, void* a4, 
    E_Int posxt, E_Int posyt, E_Int poszt, E_Int posct,
    E_Float& voli, FldArrayI& donorIndices, FldArrayF& donorCoefs,
    E_Int& type, E_Int& noDonorBlk,
    InterpData::InterpolationType interpType=InterpData::O2CF,
    E_Int nature=0, E_Int penalty=0, E_Float constraint=40.,
    E_Int extrapOrder=1);
  
  short getInterpolationCell(
    E_Float x, E_Float y, E_Float z,
    std::vector<InterpData*>& interpData,
    std::vector<FldArrayF*>& fields,
    std::vector<void*>& a1, std::vector<void*>& a2, std::vector<void*>& a3, 
    std::vector<void*>& a4, 
    std::vector<E_Int>& posxt, std::vector<E_Int>& posyt, 
    std::vector<E_Int>& poszt, std::vector<E_Int>& posct,
    E_Float& voli, FldArrayI& donorIndices, FldArrayF& donorCoefs,
    E_Int& type, E_Int& noDonorBlk,
    InterpData::InterpolationType interpType=InterpData::O2CF,
    E_Int nature=0, E_Int penalty=0);
  
  short getExtrapolationCell(
    E_Float x, E_Float y, E_Float z,
    std::vector<InterpData*>& interpDatas,
    std::vector<FldArrayF*>& fields,
    std::vector<void*>& a1, std::vector<void*>& a2, std::vector<void*>& a3, 
    std::vector<void*>& a4, 
    std::vector<E_Int>& posxt, std::vector<E_Int>& posyt, 
    std::vector<E_Int>& poszt, std::vector<E_Int>& posct,
    E_Float& voli, FldArrayI& donorIndices, FldArrayF& donorCoefs,
    E_Int& type, E_Int& noDonorBlk,
    InterpData::InterpolationType interpType=InterpData::O2CF,
    E_Int nature=0,E_Int penalty=0, E_Float constraint=40.,
    E_Int extrapOrder=1);

/* Meme routine que getInterpolationCell, mais on envoie
   les images de (x,y,z) projete sur les parois des differents
   domaines donneurs dans xt,yt,zt.
   Les pts de coordonnees (xt[noz],yt[noz],zt[noz]) sont interpoles depuis 
   l'interpData interpDatas[noz] uniquement. Ce cas est appliqué pour calculer
   les coefficients d'interpolation dans le cas double wall.
   Attention: l'ordre de (xt,yt,zt) et de interpDatas doit etre le meme 
   donorIndices et donorCoefs sont redimensionnes */
  short getInterpolationCellDW(
    E_Float* xt, E_Float* yt, E_Float* zt,
    std::vector<InterpData*>& interpData,
    std::vector<FldArrayF*>& fields,
    std::vector<void*>& a1, std::vector<void*>& a2, std::vector<void*>& a3, 
    std::vector<void*>& a4, 
    std::vector<E_Int>& posxt, std::vector<E_Int>& posyt, 
    std::vector<E_Int>& poszt, std::vector<E_Int>& posct,
    E_Float& voli, FldArrayI& donorIndices, FldArrayF& donorCoefs,
    E_Int& type, E_Int& noDonorBlk,
    InterpData::InterpolationType interpType=InterpData::O2CF,
    E_Int nature=0, E_Int penalty=0);
  
  short getExtrapolationCellDW(
    E_Float* xt, E_Float* yt, E_Float* zt,
    std::vector<InterpData*>& interpDatas,
    std::vector<FldArrayF*>& fields,
    std::vector<void*>& a1, std::vector<void*>& a2, std::vector<void*>& a3, 
    std::vector<void*>& a4, 
    std::vector<E_Int>& posxt, std::vector<E_Int>& posyt, 
    std::vector<E_Int>& poszt, std::vector<E_Int>& posct,
    E_Float& voli, FldArrayI& donorIndices, FldArrayF& donorCoefs,
    E_Int& type, E_Int& noDonorBlk,
    InterpData::InterpolationType interpType=InterpData::O2CF,
    E_Int nature=0,E_Int penalty=0, E_Float constraint=40.,
    E_Int extrapOrder=1);
  
  short compOneInterpolatedValue(E_Int* indi, FldArrayF& cf,
                                 E_Float* f0,  void* a2, void* a3, void* a4, 
                                 E_Int type, E_Float& val);
  
  short compInterpolatedValues(E_Int* indi, FldArrayF& cf,
                               FldArrayF& f0,  void* a2, void* a3, void* a4, 
                               E_Int ind, E_Int type, FldArrayF& f);
  /* meme fonction que precedente mais interpolation sur les variables de numero posvars0
  posvars0 demarre a 1.  */
  short compInterpolatedValues(E_Int* indi, FldArrayF& cf,
                               FldArrayF& f0,  void* a2, void* a3, void* a4, 
                               E_Int ind, E_Int type, FldArrayF& f, std::vector<E_Int>& posvars0);

  short compInterpolatedField(E_Int* indi, FldArrayF& cf,
                              FldArrayF& f0,  void* a2, void* a3, void* a4, 
                              E_Int type, FldArrayF& f);

  /* Fonctions internes */
  E_Int getInterpolationData(
    E_Float x, E_Float y, E_Float z,
    InterpData* interpData, void* c1, void* c2, 
    E_Int meshtype, E_Int ni, E_Int nj, E_Int nk,
    E_Float* xl, E_Float* yl, E_Float* zl, E_Float* cellN,
    E_Int& isBorder, E_Int& type, FldArrayI& indi, FldArrayF& cf, 
    InterpData::InterpolationType interpType, E_Int nature);  
  
  E_Int getExtrapolationData(
    E_Float x, E_Float y, E_Float z,
    InterpData* interpData, void* c1,  void* c2, 
    E_Int meshtype, E_Int ni, E_Int nj, E_Int nk,
    E_Float* xl, E_Float* yl, E_Float* zl, E_Float* cellN,
    E_Int& isBorder, E_Int& type, FldArrayI& indi, FldArrayF& cf, 
    InterpData::InterpolationType interpType,
    E_Int nature, E_Float constraint, E_Int extrapOrder);

  short compLagrangeCoefs(E_Float x, E_Float y, E_Float z,
                          E_Int ic, E_Int jc, E_Int kc,
                          E_Int ni, E_Int nj, E_Int nk,
                          E_Float* xl, E_Float* yl, E_Float* zl,
                          FldArrayF& cf, 
                          InterpData::InterpolationType interpType);


  /* Interpolation MLS a partir d un nuage de points donneur, en un point donne
     IN: order: ordre de la formule
     IN: dimPb: dimension du probleme
     IN: pt: point ou l'on calcule les coefficients d'interpolation
     IN: xtDnr, ytDnr, ztDnr: coordonnees x, y, z du maillage donneur
     IN: dnrIndices: indices des points du stencil donneur
     IN: radius: longueurs des demi-axes de l'ellipse
     IN: axis: direction des demi-axes de l'ellipse
     IN: axisConst: direction(s) constante(s)
     OUT: cfloc: coefficients d'interpolation 
     Retourne 1 : succes, 
     -1 : pas de base de polynomes creee
     -2 : matrice singuliere
     -3 : matrice non SDP
     -4 : matrice B non creee */ 
  E_Int getInterpCoefMLS(E_Int order, E_Int dimPb, E_Int sizeBasis,
                         E_Float* pt, 
                         E_Float* xtDnr, E_Float* ytDnr, E_Float* ztDnr, 
                         std::vector<E_Int>& dnrIndices, 
                         E_Float*radius, E_Float* axis, E_Int* axisConst, 
                         E_Float* cfLoc);
  
  E_Int matrixB(E_Int nDnrPts, E_Int sizeBasis, std::vector<E_Int>& dnrIndices, 
                E_Float *xtDnr, E_Float *ytDnr, E_Float *ztDnr, 
                E_Int order, E_Int dimPb, 
                E_Float *pt, E_Float *radius, E_Float *axis, E_Int *axisConst, 
                E_Float *B);
  E_Int matrixA(E_Int nDnrPts, E_Int sizeBasis, std::vector<E_Int>& indicesIn, 
                E_Float *xtDnr, E_Float *ytDnr, E_Float *ztDnr, 
                E_Int order, E_Int dimPb,
                E_Float *pt, E_Float *radius, E_Float *axis, E_Int *axisConst, 
                E_Float *A);
  
  E_Int polyBasis(E_Int order, E_Int dimPb, E_Int sizeBasis, 
                  E_Int *axisConst, E_Float* X,
                  std::vector<E_Float>& basis);
  E_Float weightFunction(E_Float *pt, E_Float *xi, 
                         E_Float *radius, E_Float *axis);

  E_Float distEllipse(E_Float x, E_Float y, E_Float z,
                      E_Float *center, E_Float *radius, E_Float *axis);

  /* Functions used for setInterpDataLS */
  E_Int buildStencil(E_Int dimPb, E_Int order, E_Int nature, E_Int sizeBasis,
                     E_Int center, E_Int resl, E_Int ni, E_Int nj, E_Int nk,
                     std::vector<std::vector<E_Int> >& vectOfcVN,
                     E_Float *xtDnr, E_Float *ytDnr, E_Float *ztDnr, E_Float *cellNtDnr,
                     E_Float *pt, E_Float *radius, E_Float *axis, E_Int *axisConst, 
                     std::vector<E_Int>& indicesIn);
  
  void structStencil(E_Int dimPb, E_Int indi, E_Int depth, E_Int ni, E_Int nj, E_Int nk,
                     std::vector<E_Int>& indicesIn);
  void NGONStencil(E_Int dimPb, E_Int indi, E_Int depth, std::vector< std::vector<E_Int> >& cVN, 
                   E_Float* xtDnr, E_Float* ytDnr, E_Float* ztDnr,  std::vector<E_Int>& indicesIn);

  void findRadius(E_Int dimPb, std::vector<E_Int> &indicesIn,  E_Int depth, E_Int order,
                  E_Float *pt, E_Float *xtDnr, E_Float *ytDnr, E_Float *ztDnr, 
                  E_Float *axis, E_Float *radius, E_Int *axisConst);

  void indicesInEllipse(E_Float *pt, E_Float *radius, E_Float *axis, 
                        E_Float* xtDnr, E_Float* ytDnr, E_Float* ztDnr, 
                        std::vector<E_Int> &indicesIn);
  void eraseCellN0or2(E_Float *cellNtDnr, std::vector<E_Int> &indicesIn);
  void eraseCellN0(E_Float *cellNtDnr, std::vector<E_Int> &indicesIn);
  void eraseCellN(E_Int nature,
                  E_Float *cellNtDnr, std::vector<E_Int> &indicesIn);

  E_Float distEllipse2(E_Float x, E_Float y, E_Float z,
                       E_Float *center, E_Float *radius, E_Float *axis);
  void OBbox(E_Int dimPb, E_Int n, E_Float *x, E_Float *y, E_Float *z, 
             E_Float *axis, E_Float *bbox);  
  void PCA(E_Int dimPb, E_Int n, E_Float *x, E_Float *y, E_Float *z, 
           E_Float *axis);

  void getBestDonor(
    E_Int dimPb, E_Float* pt, E_Int nDnrZones, 
    std::vector<FldArrayF*>& fields, 
    std::vector<K_SEARCH::KdTree<FldArrayF>*>& vectOfKdTrees,
    std::vector<E_Int>& resl,
    std::vector<void*>& a2, std::vector<void*>& a3, std::vector<void*>& a4,
    std::vector<E_Int>& posxs, std::vector<E_Int>& posys, 
    std::vector<E_Int>& poszs,
    std::vector<E_Int>& poscs, std::vector<FldArrayI*>& vectOfCorres, 
    std::vector<std::vector<std::vector<E_Int> > >& vectOfcVF,
    std::vector<std::vector<std::vector<E_Int> > >& vectOfcEV,
    std::vector<std::vector<std::vector<E_Int> > >& vectOfcVN,
    std::vector<FldArrayI>& vectOfcFE, 
    std::vector<FldArrayI>& vectOfPosElt, 
    std::vector<FldArrayI>& vectOfPosFace,
    std::vector<FldArrayI>& vectOfDimElt,
    E_Int depth, E_Int penalty,
    E_Int& noBest, E_Int& indBestDnr, E_Float& dBestDnr);
  E_Int isOnEdge(E_Int indice, E_Int ni, E_Int nj, E_Int nk, E_Int depth, E_Float *cellN);
  
  E_Int isOnEgdeStruct(E_Int ind, E_Int ni, E_Int nj, E_Int nk);

  E_Int isOnEgdeNGON(E_Int indV, std::vector<std::vector<E_Int> >& cVF, FldArrayI& cFE);

  E_Int penalizeBorderNGON(E_Int dimPb, E_Int indV, E_Int depth, 
                           E_Float *cellN, 
                           E_Float* xtDnr, E_Float* ytDnr, E_Float* ztDnr,
                           std::vector<std::vector<E_Int> >& cVN, 
                           std::vector<std::vector<E_Int> >& cVF, 
                           FldArrayI& cFE);
  E_Int penalizeBorderStruct(E_Int dimPb, E_Int indV, E_Int depth, 
                             E_Float *cellN,
                             E_Int ni, E_Int nj, E_Int nk);

  E_Float compMaxLengthOf1stNgonNbr(
    E_Int ind, E_Float* x, E_Float* y, E_Float* z,
    std::vector<std::vector<E_Int> >& cVN);  

  E_Float compMaxLengthOf1stStructNbr(
    E_Int dimPb, E_Int ind, E_Float* x, E_Float* y, E_Float* z,
    E_Int ni, E_Int nj, E_Int nk);

  /* calcule les elements contenant le sommet d indice indV pour un NGON conforme
     IN: indv: indice du vertex
     IN: cVF : connectivite Vertex/Faces
     IN: cFE : connectivite Faces/Elements
     OUT : ENbrs : indices des elements contenant indV  */
  void connectNG2ENbrs(E_Int indV, std::vector< std::vector<E_Int> >& cVF, 
                       FldArrayI& cFE,
                       std::vector<E_Int>& ENbrs);
}
#undef FldArrayF
#undef FldArrayI
#undef FldArrayIS 
#endif
