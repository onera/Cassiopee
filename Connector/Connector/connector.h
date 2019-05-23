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

#ifndef _CONNECTOR_CONNECTOR_H_
#define _CONNECTOR_CONNECTOR_H_
# include <locale>
# include <cctype>
# include "kcore.h"
# include "KInterp/BlkInterp.h"

# define SIZECF(type,meshtype,sizecf){                                  \
    if      (type == 1                 ) sizecf=1;                      \
    else if (type == 2 && meshtype == 1) sizecf=8;                      \
    else if (type == 3 && meshtype == 1) sizecf=9;                      \
    else if (type == 4 && meshtype == 1) sizecf=8;                      \
    else if (type == 4 && meshtype == 2) sizecf=4;                      \
    else if (type == 5 && meshtype == 1) sizecf=15;                     \
    else if (type == 22 && meshtype == 1) sizecf=4;                     \
    else sizecf=-1;                                                     \
  }

/*# define COMMONINTERPTRANSFERS(adr){                                         \
  E_Int chunk = nvars/Nbre_thread_actif;                                     \
  E_Int r = nvars - chunk*Nbre_thread_actif;                                 \
  E_Int eq_deb, eq_fin;                                                      \
  if (ithread <= r)                                                          \
  { eq_deb = (ithread-1)*(chunk+1);  eq_fin = eq_deb + (chunk+1);  }         \
  else                                                                       \
  { eq_deb = (chunk+1)*r+(ithread-r-1)*chunk; eq_fin = eq_deb + chunk; }     \
  E_Float* ptrCoefs = donorCoefsF->begin();                                  \
  E_Int indR, type;                                                          \
  E_Int indD0, indD, i, j, k, ncfLoc, nocf;                                  \
  E_Int noi = 0;                                                             \
  E_Int sizecoefs = 0;                                                       \
                                                                             \
  for (E_Int noind = 0; noind < nbRcvPts; noind++)                           \
  {                                                                          \
    indR = adr;                                                              \
# include "commonInterpTransfers.h"                                      \
    ptrCoefs += sizecoefs;                                                   \
  }                                                                          \
 }       
*/
                                                                    

# define BLOCKRELEASEMEM \
  RELEASESHAREDN(pyIndRcv, rcvPtsI);                    \
  RELEASESHAREDN(pyIndDonor, donorPtsI);                \
  RELEASESHAREDN(pyArrayTypes, typesI);                 \
  RELEASESHAREDN(pyArrayXPC, coordxPC);                 \
  RELEASESHAREDN(pyArrayYPC, coordyPC);                 \
  RELEASESHAREDN(pyArrayZPC, coordzPC);                 \
  RELEASESHAREDN(pyArrayXPW, coordxPW);                 \
  RELEASESHAREDN(pyArrayYPW, coordyPW);                 \
  RELEASESHAREDN(pyArrayZPW, coordzPW);                 \
  RELEASESHAREDN(pyArrayXPI, coordxPI);                 \
  RELEASESHAREDN(pyArrayYPI, coordyPI);                 \
  RELEASESHAREDN(pyArrayZPI, coordzPI);                 \
  RELEASESHAREDN(pyArrayCoefs, donorCoefsF);            

# define BLOCKRELEASEMEMD \
  RELEASESHAREDN(pyIndDonor, donorPtsI);                \
  RELEASESHAREDN(pyArrayTypes, typesI);                 \
  RELEASESHAREDN(pyArrayXPC, coordxPC);                 \
  RELEASESHAREDN(pyArrayYPC, coordyPC);                 \
  RELEASESHAREDN(pyArrayZPC, coordzPC);                 \
  RELEASESHAREDN(pyArrayXPW, coordxPW);                 \
  RELEASESHAREDN(pyArrayYPW, coordyPW);                 \
  RELEASESHAREDN(pyArrayZPW, coordzPW);                 \
  RELEASESHAREDN(pyArrayXPI, coordxPI);                 \
  RELEASESHAREDN(pyArrayYPI, coordyPI);                 \
  RELEASESHAREDN(pyArrayZPI, coordzPI);                 \
  RELEASESHAREDN(pyArrayCoefs, donorCoefsF);            

#define BLOCKRELEASEMEM2 \
  if ( pyArrayPressure != Py_None)                      \
    RELEASESHAREDN(pyArrayPressure, pressF);            \
  if ( pyArrayDens != Py_None)                          \
    RELEASESHAREDN(pyArrayDens, densF);                 \
  if ( pyArrayVx != Py_None)                            \
    RELEASESHAREDN(pyArrayVx, vxF);                     \
  if ( pyArrayVy != Py_None)                            \
    RELEASESHAREDN(pyArrayVy, vyF);                     \
  if ( pyArrayVz != Py_None)                            \
    RELEASESHAREDN(pyArrayVz, vzF);                     \
  if ( pyArrayUtau != Py_None)                          \
    RELEASESHAREDN(pyArrayUtau, utauF);                 \
  if ( pyArrayYplus != Py_None)                         \
    RELEASESHAREDN(pyArrayYplus, yplusF);                              


extern "C"
{
  void spalart_1d_(E_Int& ithread, E_Float* y, E_Float* matm,E_Float* mat,E_Float* matp,E_Float* nutilde, E_Float* utble, E_Float& pdtc, E_Float& nu, E_Float& nutildeext, E_Int& jmax, E_Float& kappa);
}  

namespace K_CONNECTOR 
{
/* Interpolation datas.*/
  struct InterpData
  {
      /*
        Interpolation coefficients on fine grid. The coarse cells
        use another method for interpolation.
      */
      K_FLD::FldArrayF _coefs;
      /*
        Interpolation cells per interpolation domains.
      */
      K_FLD::FldArrayI _cells;    
      /*
        Volume of interpolation cells per interpolation domains.
      */
      K_FLD::FldArrayF _vols;    
      /*
        Indirection per interpolations domains
      */
      K_FLD::FldArrayI _indirectionPerDom;
      /*
        Number of interpolated cells per interpolations domains.
      */
      K_FLD::FldArrayI _sizeOfEachIndirectionTab;
      /*
        Interpolation type stored for each interpolated point
      */
      K_FLD::FldArrayI _interpTypes;
  };

/* X-Ray mask */
  struct XRayPlane
  {
      E_Float xmin;
      E_Float ymin;
      E_Float zmin;
      E_Float xmax;
      E_Float ymax;
      E_Float zmax;
      E_Int ni, nj;
      E_Float hi, hj;
      K_FLD::FldArrayI indir; // Compact storage indirection
      K_FLD::FldArrayF Z; // Compact Z storage (increasing order)
      std::vector<E_Float>* tempZ; // temporary storage of Z
  };

/* C est a l appelant de detruire les planes*/
  E_Int compCharacteristics(E_Int isNot, E_Int elevationDir, 
                            E_Int dim1, E_Int dim2, 
                            E_Float tol, E_Float delta,
                            std::vector<E_Int>& posxt, 
                            std::vector<E_Int>& posyt, 
                            std::vector<E_Int>& poszt,
                            std::vector<K_FLD::FldArrayF*>& fieldt, 
                            std::vector<K_FLD::FldArrayI*>& cns,
                            std::list<XRayPlane*>& planes,
                            E_Float& xmin, E_Float& ymin, E_Float& zmin,
                            E_Float& xmax, E_Float& ymax, E_Float& zmax);

/*xmin, ymin, zmin : min de la bbox globale
  Retourne comp : le nb de pts ambigus */
  E_Int computeZ(E_Int elevationDir,  
                 E_Float xmin, E_Float ymin, E_Float xmax, E_Float ymax, 
                 K_FLD::FldArrayF& epsilon, 
                 E_Float* xtb, E_Float* ytb, E_Float* ztb, 
                 K_FLD::FldArrayI& cnb, struct XRayPlane* p);

/*xmina,ymina, zmina : min de la bbox globale
  Retourne comp : le nb de pts ambigus */
  E_Int triangleOnPlane(E_Int elevationDir, 
                        E_Float xmina, E_Float ymina, 
                        E_Float xmaxa, E_Float ymaxa,
                        K_FLD::FldArrayF& epsilon, 
                        E_Float* xt, E_Float* yt, E_Float* zt, 
                        E_Int ind1, E_Int ind2, E_Int ind3,
                        struct XRayPlane* p);
/* Compute intersection between triangle defined by (x0, y0, z0), 
   (x1, y1, z1), (x2, y2, z2) and ray (x,y). Result in z if found.
   Return value is : 0 non-intersection
   1 intersection
   -1 ambiguous intersection, face is along z
   -2 intersection on a node of triangle*/
  E_Int compIntersect(E_Int elevationDir,
                      E_Float x0, E_Float y0, E_Float z0,
                      E_Float x1, E_Float y1, E_Float z1,
                      E_Float x2, E_Float y2, E_Float z2,
                      E_Float x, E_Float y, 
                      E_Float zmin, E_Float zmax,
                      E_Float& z);
  void
  compactZ(E_Int isNot, E_Int elevationDir, E_Float delta, E_Float tol,
           E_Float xmin, E_Float ymin, E_Float zmin, 
           E_Float xmax, E_Float ymax, E_Float zmax,             
           XRayPlane* p);

/* add delta to intersection points ordinates */
  void addDeltaToZ(E_Int isNot, E_Int elevationDir, E_Float delta, 
                   E_Float xmin, E_Float ymin, E_Float zmin, 
                   E_Float xmax, E_Float ymax, E_Float zmax,             
                   XRayPlane* p);

/* mask X-Ray delta*/
  E_Int holeExpansionStruct(E_Int elevationDir, E_Int blankingType, 
                            E_Int isNot, E_Float delta,
                            std::list<XRayPlane*>& planes,
                            E_Int ni, E_Int nj, E_Int nk, 
                            E_Int posx, E_Int posy, E_Int posz,
                            K_FLD::FldArrayF& field, 
                            K_FLD::FldArrayI& cellNatFld);
  E_Int holeExpansionUnstr(E_Int elevationDir, E_Int blankingType, 
                           E_Int isNot, E_Float delta,
                           std::list<XRayPlane*>& planes,
                           E_Int posx, E_Int posy, E_Int posz,
                           K_FLD::FldArrayF& field, K_FLD::FldArrayI& cn,
                           K_FLD::FldArrayI& cellNatFld);

/* Calcule la liste des points interpoles :
   IN : elevationDir : 2 = 2D
   IN: type : blankingType 
   IN: depth=1 ou 2 
   IN: nic, njc, nkc: dimensions de cellNatFld
   IN: cellNatFld: cellNatureField
   OUT: tableau des indices de points interpoles 
   OUT: if type = 0 diri and dirj are built and define the sign of delta expansion */

  void compListOfInterpolatedPoints(
    E_Int elevationDir, E_Int type, E_Int depth,
    E_Int nic, E_Int njc, E_Int nkc,
    K_FLD::FldArrayI& cellNatFld, 
    K_FLD::FldArrayI& listOfInterpolatedPoints,
    K_FLD::FldArrayI& diri, K_FLD::FldArrayI& dirj);

  void compXRayMaskInfo(E_Int elevationDir, std::list<XRayPlane*>& planes, 
                        K_FLD::FldArrayF& coord);

/* Construit le tableau des priorites pour toutes les zones. Vaut 0
   par defaut pour toutes les zones */
  E_Int getPriorities(PyObject* priorities, K_FLD::FldArrayIS& prios);

/* comparaison de la taille des cellules interpolables/d interpolation  
   des blk1 et blk2. Interpolation en centres etendus effectuee. 
   Retourne les celln modifies
   IN: ni1, nj1, nk1: liste des dimensions des grilles en centres 
   IN: msh1: kmesh en noeuds
   IN: liste des interpdatas basees sur les centres etendus
   retourne les celln modifies */
  void compareInterpCells(
    E_Int ni1, E_Int nj1, E_Int nk1, 
    E_Int nie1, E_Int nje1, E_Int nke1, K_FLD::FldArrayF* extCenters1,
    K_INTERP::InterpAdt* interpData1,
    E_Float* xc1, E_Float* yc1, E_Float* zc1, E_Float* celln1,
    E_Float* vol1, K_FLD::FldArrayI& interpCells1, K_FLD::FldArrayI& tag1,
    E_Int ni2, E_Int nj2, E_Int nk2, 
    E_Int nie2, E_Int nje2, E_Int nke2, K_FLD::FldArrayF* extCenters2,
    K_INTERP::InterpAdt* interpData2,
    E_Float* xc2, E_Float* yc2, E_Float* zc2, E_Float* celln2,
    E_Float* vol2, K_FLD::FldArrayI& interpCells2, K_FLD::FldArrayI& tag2, E_Int isDW);


/* Modification du cellN des deux zones base sur le critere de 
   masquage de la cellule de plus grand volume */
  void modifyCellNWithVolCriterion(
    E_Int ni1, E_Int nj1, E_Int nk1, E_Float* vol1,
    E_Int nie1, E_Int nje1, E_Int nke1, K_FLD::FldArrayF* extCenters1,
    K_INTERP::InterpAdt* interpData1,
    E_Float* xc1, E_Float* yc1, E_Float* zc1, E_Float* celln1,
    E_Int ni2, E_Int nj2, E_Int nk2, E_Float* vol2,
    E_Int nie2, E_Int nje2, E_Int nke2, K_FLD::FldArrayF* extCenters2,
    K_INTERP::InterpAdt* interpData2,
    E_Float* xc2, E_Float* yc2, E_Float* zc2, E_Float* celln2, E_Int isDW);

  void modifyCellNWithPriority(
    E_Int prio1, E_Int ni1, E_Int nj1, E_Int nk1, E_Float* vol1,
     E_Int nie1, E_Int nje1, E_Int nke1, K_FLD::FldArrayF* extCenters1,
    K_INTERP::InterpAdt* interpData1,
    E_Float* xc1, E_Float* yc1, E_Float* zc1, E_Float* celln1,
    E_Int prio2, E_Int ni2, E_Int nj2, E_Int nk2,E_Float* vol2,
    E_Int nie2, E_Int nje2, E_Int nke2, K_FLD::FldArrayF* extCenters2,
    K_INTERP::InterpAdt* interpData2,
    E_Float* xc2, E_Float* yc2, E_Float* zc2, E_Float* celln2, E_Int isDW);

  void modifyCellnForDonorCells(
    K_FLD::FldArrayI& tag1, E_Float* celln1, K_FLD::FldArrayI& tag2, E_Float* celln2, 
    K_FLD::FldArrayI& interpCells2);

/* blank cells. Return the modified cell nature fields
   IN: delta: distance to the body for blanked cells
   IN: dim1,dim2: dimensions for XRay plane
   IN: nit, njt, nkt: dimensions of meshes (nodes)
   IN: nibt, njbt, nkbt: dimensions of body arrays
   IN: coord: mesh  coordinates in nodes 
   IN/OUT: celln: list of cellnaturefields : centers
   IN: fieldsb: list of coordinates of arrays defining the mask
   IN: cnb: connectivite associee
*/
  void blankCellsStruct(E_Int elevationDir, E_Int isNot, E_Int blankingType,
                        E_Float delta, E_Float tol, E_Int dim1, E_Int dim2,
                        std::vector<E_Int>& posxt, std::vector<E_Int>& posyt, 
                        std::vector<E_Int>& poszt,
                        std::vector<E_Int>& nit, std::vector<E_Int>& njt, 
                        std::vector<E_Int>& nkt, 
                        std::vector<K_FLD::FldArrayF*>& blankedCoords, 
                        std::vector<K_FLD::FldArrayI*>& cellns,
                        std::vector<E_Int>& posxb, std::vector<E_Int>& posyb, 
                        std::vector<E_Int>& poszb,
                        std::vector<K_FLD::FldArrayF*>& fieldsb,
                        std::vector<K_FLD::FldArrayI*>& cnb);

/* Meme fonction mais avec des maillages non structures. Attention, dans 
   ce cas, seul le blankingType 0 s'applique.
   IN: blankedCoords: coordonnees en noeuds, cnt connectivite associee
   IN/OUT: cellns: cellN defini aux sommets des elements
   IN: fieldsb: corps servant pour le masque 
*/
  void 
  blankCellsUnstr(E_Int elevationDir, E_Int isNot, E_Int blankingType,
                  E_Float delta, E_Float tol, E_Int dim1, E_Int dim2,
                  std::vector<E_Int>& posxt,  std::vector<E_Int>& posyt, 
                  std::vector<E_Int>& poszt,
                  std::vector<K_FLD::FldArrayF*>& blankedCoords, std::vector<K_FLD::FldArrayI*>& cnt,
                  std::vector<K_FLD::FldArrayI*>& cellns, std::vector<K_FLD::FldArrayI*>& cntc,
                  std::vector<E_Int>& posxb, std::vector<E_Int>& posyb, 
                  std::vector<E_Int>& poszb, std::vector<K_FLD::FldArrayF*>& fieldsb,
                  std::vector<K_FLD::FldArrayI*>& cnb);

  E_Int searchForBlankedCellsStruct(E_Int elevationDir, E_Int blankingType,
                                    E_Int isNot, E_Float delta,
                                    std::list<XRayPlane*>& planes,
                                    E_Int ni, E_Int nj, E_Int nk,
                                    E_Int posx, E_Int posy, E_Int posz,
                                    K_FLD::FldArrayF& field,
                                    K_FLD::FldArrayI& blankedCell);
  E_Int searchForBlankedCellsUnstr(E_Int elevationDir, E_Int blankingType,
                                   E_Int isNot, E_Float delta,
                                   std::list<XRayPlane*>& planes,
                                   E_Int posx, E_Int posy, E_Int posz,
                                   K_FLD::FldArrayF& field, K_FLD::FldArrayI& cn,
                                   K_FLD::FldArrayI& blankedCell);
  /* Blank cells for strand grids*/
  E_Int blankAboveCellsStruct(E_Int indcell1, E_Int indcell2, 
                              E_Int nic, E_Int njc, E_Int nkc, E_Int nicnjc, 
                              E_Float* cellN);
  E_Int blankAboveCellsPenta(E_Int etg, E_Int etd, K_FLD::FldArrayI& cn,
                             E_Float* cellN);
  E_Int blankAboveCellsHexa(E_Int etg, E_Int etd, K_FLD::FldArrayI& cn,
                            E_Float* cellN);
  E_Int blankIntersectingCellsStruct(
    E_Float eps, E_Int posx, E_Int posy, E_Int posz, E_Int posc,
    std::vector<E_Int>& nit, std::vector<E_Int>& njt, std::vector<E_Int>& nkt, 
    std::vector<K_FLD::FldArrayF*>& structF, 
    std::vector<E_Int>& nict, std::vector<E_Int>& njct, 
    std::vector<E_Int>& nkct, 
    std::vector<K_FLD::FldArrayF*>& structFc);
  
  E_Int blankInvalidCellsHexa( 
    E_Float eps, E_Int posx, E_Int posy, E_Int posz, E_Int posc,
    std::vector<K_FLD::FldArrayI*>& cnt, std::vector<K_FLD::FldArrayF*>& unstrF, 
    std::vector<K_FLD::FldArrayI*>& cnct, std::vector<K_FLD::FldArrayF*>& unstrFc);
  E_Int blankInvalidCellsPenta( 
    E_Float eps, E_Int posx, E_Int posy, E_Int posz, E_Int posc,
    std::vector<K_FLD::FldArrayI*>& cnt, std::vector<K_FLD::FldArrayF*>& unstrF, 
    std::vector<K_FLD::FldArrayI*>& cnct, std::vector<K_FLD::FldArrayF*>& unstrFc);
  E_Int blankInvalidCellsStruct(
    E_Float eps, E_Int posx, E_Int posy, E_Int posz, E_Int posc,
    std::vector<E_Int>& nit, std::vector<E_Int>& njt, std::vector<E_Int>& nkt, 
    std::vector<K_FLD::FldArrayF*>& structF, 
    std::vector<E_Int>& nict, std::vector<E_Int>& njct,
    std::vector<E_Int>& nkct, 
    std::vector<K_FLD::FldArrayF*>& structFc);

  E_Int blankIntersectingCellsPenta( 
    E_Float eps, E_Int posx, E_Int posy, E_Int posz, E_Int posc,
    std::vector<K_FLD::FldArrayI*>& cnt, 
    std::vector<K_FLD::FldArrayF*>& unstrF, 
    std::vector<K_FLD::FldArrayI*>& cnct, 
    std::vector<K_FLD::FldArrayF*>& unstrFc);
  E_Int blankIntersectingCellsHexa( 
    E_Float eps, E_Int posx, E_Int posy, E_Int posz, E_Int posc,
    std::vector<K_FLD::FldArrayI*>& cnt, 
    std::vector<K_FLD::FldArrayF*>& unstrF, 
    std::vector<K_FLD::FldArrayI*>& cnct, 
    std::vector<K_FLD::FldArrayF*>& unstrFc);
  /* For HEXA/PENTA mesh generated by extrusion, determines the number of 
     nk layers in the normal direction to the initial surface 
     IN: npts: nb de noeuds dans le maillage
     IN: cEV: connectivite elts->noeuds
     retourne 0 si pb, or nk the number of layers */
  E_Int getNbOfHexaLayers(E_Int npts, K_FLD::FldArrayI& cEV);
  E_Int getNbOfPentaLayers(E_Int npts, K_FLD::FldArrayI& cEV);

/* Verification de la coherence des dimensions entre coords et celln*/
  E_Int checkDimensions(E_Int im, E_Int jm, E_Int km, 
                        E_Int imc, E_Int jmc, E_Int kmc);

/* Modification du celln sur les parois doubly defined lorsque le 
   point n est pas interpolable
   IN : xc, yc, zc : maillage en centres de taille imc,jmc,kmc */
  E_Int modifyCellNForDoublyDefined(
    E_Int imc, E_Int jmc, E_Int kmc, E_Int depth, E_Int* range, 
    E_Float* xc, E_Float* yc, E_Float* zc, E_Float* celln,
    E_Int posxt, E_Int posyt, E_Int poszt, E_Int posct, 
    std::vector<E_Int>& nit, std::vector<E_Int>& njt, 
    std::vector<E_Int>& nkt, std::vector<K_FLD::FldArrayF*>& structF, 
    std::vector<E_Int>& nitc, std::vector<E_Int>& njtc, 
    std::vector<E_Int>& nktc, std::vector<K_FLD::FldArrayF*>& structFc);

/*---------------------------------------------------*/
/* Routines internes a testInterpolationCellCenters  */ 
/*---------------------------------------------------*/
/* recherche des cellules d'interpolations valides et calcule les coeff d interpolation 
   IN: interpPts: pts a interpoler 
   IN: nit, njt, nkt: dimensions des domaines d interpolation en noeuds
   IN: coords: coordonnees en noeuds des domaines d interpolation
   IN: cellns: celln des domaines d interpolation - en centres
   OUT: donorCells: indices des cellules donneuses dans le maillage en centres etendus
   vectOfExtrapPts doit etre dimensionne au nombre de domaines d interpolation */
  void compAndStoreInterpCoefs(
    E_Float geomCutOff,
    K_KINTERP::BlkInterpData::InterpolationType interpType,
    K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
    std::vector<K_FLD::FldArrayF*>& interpPts, InterpData* interpData, K_FLD::FldArrayI& donorCells, 
    std::vector< std::vector<E_Int> >& vectOfExtrapPts, std::vector<E_Int>& orphanVect,
    std::vector<E_Int>& nit, std::vector<E_Int>& njt, std::vector<E_Int>& nkt,
    std::vector<K_FLD::FldArrayF*>& coords,
    std::vector<K_KINTERP::BlkInterpData*>& listOfInterpDatas, 
    std::vector<K_FLD::FldArrayI*>& cellns, E_Int zid, E_Float cfMax);
/* recherche des pts EX valides et calcule les coeff d interpolation 
   IN: interpPts: pts EX a interpoler 
   IN: nit, njt, nkt: dimensions des domaines d interpolation en noeuds
   IN: coords: coordonnees en noeuds des domaines d interpolation
   IN: cellns: celln des domaines d interpolation - en centres
   OUT: donorCells: indices des cellules donneuses dans le maillage en centres etendus
*/
  void compAndStoreEXInterpCoefs(
    E_Float geomCutOff, 
    K_KINTERP::BlkInterpData::InterpolationType interpType,
    K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
    std::vector<K_FLD::FldArrayF*>& interpPts, InterpData* interpDataEX, K_FLD::FldArrayI& donorCellsEX, 
    std::vector< std::vector<E_Int> >& vectOfExtrapPts, std::vector<E_Int>& orphanVectEX,
    std::vector<E_Int>& nit, std::vector<E_Int>& njt, std::vector<E_Int>& nkt,
    std::vector<K_FLD::FldArrayF*>& coords,
    std::vector<K_KINTERP::BlkInterpData*>& listOfInterpDatas, 
    std::vector<K_FLD::FldArrayI*>& cellns, E_Int zid, E_Float cfMax);

  /* calcul des transferts Chimere
   IN: lArraysDonor: liste des champs des domaines donneurs 
   IN: lArraysCoef: liste des coefficients d interpolation 
   IN/OUT: lArraysRcv: liste des champs interpoles, sur lesquels le transfert est applique 
*/
  void chimeraTransfer(
    std::vector<K_FLD::FldArrayF*>& lArraysDonor,
    std::vector<K_FLD::FldArrayF*>& lArraysCoef,
    std::vector<K_FLD::FldArrayF*>& lArraysRcv);

  /* OUT: cf: coefficients non directionnels (O2CF)
          indExtrap: tableau des indices de cellule d interpolation 
          indExtrap est de taille 7 pour interpType = O2CF*/
  short solveOrphanPoints(K_KINTERP::BlkInterpData::InterpolationType interpType,
                          K_KINTERP::BlkInterpData::InterpMeshType interpMeshType,
                          K_KINTERP::BlkInterpData* blkInterpData,
                          E_Int ni, E_Int nj, E_Int nk, K_FLD::FldArrayI& cellN,
                          E_Float x, E_Float y, E_Float z, 
                          E_Int testNature, E_Float cfMax,
                          K_FLD::FldArrayF& cf, K_FLD::FldArrayI& indExtrap, 
                          E_Int& ic, E_Int& jc, E_Int& kc);
  E_Float compInterpolatedNature(E_Int ni, E_Int nj, E_Int nk,
                                 E_Int indiSize, E_Int* indi,
                                 K_FLD::FldArrayF& cf, E_Int* cellNp,
                                 E_Int interpType);
  /* changeWall basé sur la courbure : modifie le maillage z à partir 
     d'une liste de parois sur lesquelles projeter. Ne projette que des 
     points de cellN = 2
     IN: imc, jmc, kmc: dimensions du maillage z a projeter 
     IN: cellN: champ associe, vaut 0, 1 ou 2
     IN: indicesw: indices des premiers points pres de la paroi
     IN: dirw1: vaut 1 si le pt est sur la paroi i=1, -1 si sur la paroi im, 0 sinon
     IN: dirw2: vaut 1 si le pt est sur la paroi j=1, -1 si sur la paroi jm, 0 sinon
     IN: dirw3: vaut 1 si le pt est sur la paroi k=1, -1 si sur la paroi km, 0 sinon
     IN: hmaxw: hmax pour les pts de indicesw 
     IN: unstrF, cnt: maillages TRI des surfaces sur lesquelles projeter 
     IN: posxt, posyt, poszt: position de x, y et z dans ces maillages TRI
     IN: xc, yc, zc: maillage à projeter en centres, les coordonnées sont modifiées si cellN = 2 
     OUT: xc2, yc2, zc2 : maillage projete resultant */
  void changeWall(E_Int imc, E_Int jmc, E_Int kmc, E_Float* cellN, 
                  E_Int nbCentersW, E_Float* indicesw, E_Float* dirw1, E_Float* dirw2, E_Float* dirw3, E_Float* hmaxw, 
                  std::vector<E_Int> posxt, std::vector<E_Int> posyt, std::vector<E_Int> poszt, 
                  std::vector<E_Int> posht, std::vector<E_Int> posct,
                  std::vector<K_FLD::FldArrayI*>& cnt, std::vector<K_FLD::FldArrayF*>& unstrF,
                  E_Float* xc, E_Float* yc, E_Float* zc,
                  E_Float* xc2, E_Float* yc2, E_Float* zc2,
                  E_Float planartol=0.);
  
  void shiftAbovePoints(E_Int imc, E_Int jmc, E_Int kmc,
                        E_Int dir, E_Int indA, E_Int iA, E_Int jA, E_Int kA,
                        E_Float xa, E_Float ya, E_Float za, 
                        E_Float deltax, E_Float deltay, E_Float deltaz,
                        E_Float* xc, E_Float* yc, E_Float* zc, E_Float* cellN,
                        E_Float* xc2, E_Float* yc2, E_Float* zc2);
  /**/
  void changeWallEX(E_Int nEXPts, E_Float* xEX, E_Float* yEX, E_Float* zEX,
                    E_Float* neighbourCellsp1, E_Float* neighbourCellsp2, E_Float* dirEXp, E_Float* nodeMinp,
                    E_Int im, E_Int jm, E_Int km, E_Float* xn, E_Float* yn, E_Float* zn,
                    E_Int imc, E_Int jmc, E_Int kmc, E_Float* xc, E_Float* yc, E_Float* zc, E_Float* cellnc,
                    E_Int nbCentersW, E_Float* indicesw, E_Float* dirw1, E_Float* dirw2, E_Float* dirw3, E_Float* hmaxw,
                    std::vector<E_Int> posxt, std::vector<E_Int> posyt, std::vector<E_Int> poszt, 
                    std::vector<E_Int> posht, std::vector<E_Int> posct,
                    std::vector<K_FLD::FldArrayI*>& cnt, std::vector<K_FLD::FldArrayF*>& unstrF,
                    E_Float planartol=0.);

  /* Calcul du centre d interface au dessus du centre d interface dont le noeud min est indn
     IN: dirEX : direction du pt EX (interface i constante => 1)
     IN: indn: indice du noeud min associe au pt EX
     IN: indc: indice du centre dont provient le pt EX
     IN: dir: direction de recherche du above pt associe a indcell
     IN: im, jm: dimensions du maillage en noeuds
     IN: imc, jmc : dimensions du maillage en centres
     IN: xn,yn,zn: coordonnees des noeuds
     OUT: xb,yb,zb: coordonnees du pt EX au dessus de indEX */  
  void compAbovePointEX(E_Int dirEX, E_Int indn, E_Int indc, E_Int dir,
                        E_Int im, E_Int jm, E_Int imc, E_Int jmc,
                        E_Float* xn, E_Float* yn, E_Float* zn,
                        E_Float& xb, E_Float& yb, E_Float& zb);

  /* Elimination des points interpoles non necessaires (passage en cellN = 0) 
     selon le nombre de rangees de cellules d'interpolation autorise depth 
     si dir=0: molecule directionelle, si dir=1: molecule complete */
  E_Int getRidOfInterpPoints(E_Float* celln, 
                             E_Int im, E_Int jm, E_Int km, 
                             E_Int depth, E_Int dir);

  /* Determines the depth layers of interpolated points (nodes or centers) */ 
  void searchMaskInterpolatedNodesUnstr(E_Int depth, K_FLD::FldArrayI& cnEV,
                                        K_FLD::FldArrayI& blankedCells,
                                        K_FLD::FldArrayI& cellN);
  void searchMaskInterpolatedCellsNGON(E_Int depth, K_FLD::FldArrayI& cNG,
                                       K_FLD::FldArrayI& blankedCells,
                                       K_FLD::FldArrayI& cellN);
  void searchMaskInterpolatedCellsUnstr(char* eltType, 
                                        E_Int depth, K_FLD::FldArrayI& cnEV,
                                        K_FLD::FldArrayI& blankedCells,
                                        K_FLD::FldArrayI& cellN);
  void searchMaskInterpolatedCellsStruct(E_Int imc, E_Int jmc, E_Int kmc, E_Int depth,
                                         E_Int dir,
                                         K_FLD::FldArrayI& blankedCells,
                                         K_FLD::FldArrayI& cellN);

  /* Functions used for gatherMatching functions */
  void compIncrement(E_Int indwA1, E_Int imw1, E_Float* oppositeWins, E_Float* oppositePts,
                      E_Int& inci, E_Int& incj, E_Int& inciopp, E_Int& incjopp);

  void compTrirac(E_Int im1, E_Int jm1, E_Int im2, E_Int jm2, 
                  E_Int typewin1, E_Int inc1, E_Int inc2, 
                  E_Int typewin2, E_Int incm1, E_Int incm2,
                  std::vector<E_Int>& rac1, std::vector<E_Int>& rac2, std::vector<E_Int>& rac3);
  
  void compTrirac2D(E_Int im1, E_Int jm1, E_Int im2, E_Int jm2, 
                    E_Int typewin1, E_Int typewin2, E_Int inc1, E_Int incm1,
                    std::vector<E_Int>& rac1, std::vector<E_Int>& rac2, std::vector<E_Int>& rac3);

  E_Int signature(E_Int r1, E_Int r2, E_Int r3, E_Int var);

  void getIndicesInBlk(E_Int isw1, E_Int iew1, E_Int jsw1, E_Int jew1, 
                       E_Int imw1, E_Int jmw1, E_Int typewin1, 
                       E_Int im1, E_Int jm1, E_Int km1,
                       E_Int& imin, E_Int& imax, E_Int& jmin, E_Int& jmax, E_Int& kmin, E_Int& kmax);

  /* Transferts IBC avec variables conservatives en entree/sortie */
  E_Int setIBCTransfersCommonVar1(E_Int bctype,
                                  E_Int* rcvPtsI, E_Int& nbRcvPts, E_Int& ideb, E_Int& ifin, E_Int& ithread,
                                  E_Float* xPC, E_Float* yPC, E_Float* zPC,
                                  E_Float* xPW, E_Float* yPW, E_Float* zPW,
                                  E_Float* xPI, E_Float* yPI, E_Float* zPI, 
                                  E_Float* densPtr, E_Float* pressPtr, 
                                  E_Float* vxPtr, E_Float* vyPtr, E_Float* vzPtr, 
                                  E_Float* utauPtr, E_Float* yplusPtr,
                                  E_Float* d1, E_Float* d2, E_Float* d3, E_Float* d4, E_Float* d5,
                                  E_Float* tmp, E_Int&  size,
                                  E_Float gamma, E_Float cv, E_Float muS, E_Float Cs, E_Float Ts, E_Float Pr,
                                  std::vector<E_Float*>& WIn,
                                  std::vector<E_Float*>& WOut);

  /* Transferts IBC avec variables (ro,u,v,w,t) en entree/sortie */
  E_Int setIBCTransfersCommonVar2(E_Int bctype, 
                                  E_Int* rcvPtsI, E_Int& nbRcvPts, E_Int& ideb, E_Int& ifin, E_Int& ithread,
                                  E_Float* xPC, E_Float* yPC, E_Float* zPC,
                                  E_Float* xPW, E_Float* yPW, E_Float* zPW,
                                  E_Float* xPI, E_Float* yPI, E_Float* zPI, 
                                  E_Float* densPtr, E_Float* pressPtr, 
                                  E_Float* vxPtr, E_Float* vyPtr, E_Float* vzPtr, 
                                  E_Float* utauPtr, E_Float* yplusPtr,
                                  E_Float* d1, E_Float* d2, E_Float* d3, E_Float* d4, E_Float* d5,
                                  E_Float* tmp, E_Int&  size,
                                  E_Float gamma, E_Float cv, E_Float muS, E_Float Cs, E_Float Ts, E_Float Pr,
                                  std::vector<E_Float*>& WIn, std::vector<E_Float*>& WOut, 
                                  E_Int nbptslinelets=0, E_Float* linelets=NULL, E_Int* indexlinelets=NULL);

  /* Transferts IBC avec variables (ro,u,v,w,p) en entree/sortie */
  E_Int setIBCTransfersCommonVar3(E_Int bctype,
                                  E_Int* rcvPtsI, E_Int& nbRcvPts, E_Int& ideb, E_Int& ifin, E_Int& ithread,
                                  E_Float* xPC, E_Float* yPC, E_Float* zPC,
                                  E_Float* xPW, E_Float* yPW, E_Float* zPW,
                                  E_Float* xPI, E_Float* yPI, E_Float* zPI, 
                                  E_Float* densPtr, E_Float* pressPtr, 
                                  E_Float* vxPtr, E_Float* vyPtr, E_Float* vzPtr,
                                  E_Float* utauPtr, E_Float* yplusPtr, 
                                  E_Float* d1, E_Float* d2, E_Float* d3, E_Float* d4, E_Float* d5,
                                  E_Float* tmp, E_Int&  size,
                                  E_Float gamma, E_Float cv, E_Float muS, E_Float Cs, E_Float Ts, E_Float Pr,
                                  std::vector<E_Float*>& WIn, std::vector<E_Float*>& WOut);

  /* For setInterpDataGC : a mettre dans KCore ? */
  /*Determination des indices des points a modifier 
    rcvIndices doit deja etre alloue 
    IN : loc = 0 : noeuds, 1 : centres
    IN : dir: -1,1,-2,2,-3,3
    IN : dim : dimensions du pb 
    IN : win=[wmin,...,wkmax] indices min,max de la fenetre a traiter, dimensionnee
    dans le maillage ghost cells
    IN : rindwin : decalage des pts dans la direction de la fenetre
    IN : imr, jmr, kmr : dimensions de la zone (en noeuds meme si loc = 1)
    OUT : indices globaux des rindwin pts a partir de la fenetre win*/
  void getRcvIndices(E_Int imr, E_Int jmr, E_Int kmr,
                     E_Int wimin, E_Int wimax, 
                     E_Int wjmin, E_Int wjmax, 
                     E_Int wkmin, E_Int wkmax,
                     E_Int rindwin,
                     E_Int dir, E_Int dim, E_Int loc,
                     K_FLD::FldArrayI& rcvIndices);  

  /* Fonctions pour les lois de paroi */
  E_Float logf(E_Float x, E_Float a, E_Float b, E_Float c, E_Float d);
  E_Float logfprime(E_Float x, E_Float a, E_Float d) ;
  E_Float musker(E_Float x, E_Float a, E_Float b);
  E_Float muskerprime(E_Float x, E_Float a, E_Float b);
  E_Float fnutilde(E_Float nutilde, E_Float nut, E_Float rho, E_Float xmu);
  E_Float fnutildeprime(E_Float nutilde, E_Float nut, E_Float rho, E_Float xmu);

  /*---------------------------------------------------*/
  PyObject* getIBMPtsBasic(PyObject* self, PyObject* args);
  PyObject* getIBMPtsWithFront(PyObject* self, PyObject* args);
  PyObject* getIBMPtsWithoutFront(PyObject* self, PyObject* args);
  PyObject* optimizeOverlap(PyObject* self, PyObject* args);
  PyObject* maximizeBlankedCells( PyObject* self, PyObject* args );
  PyObject* blankCells( PyObject* self, PyObject* args);
  PyObject* _blankCells( PyObject* self, PyObject* args);
  PyObject* blankCellsTetra( PyObject* self, PyObject* args);
  PyObject* createTetraMask( PyObject* self, PyObject* args);
  PyObject* deleteTetraMask( PyObject* self, PyObject* args); 
  PyObject* createTriMask( PyObject* self, PyObject* args);
  PyObject* deleteTriMask( PyObject* self, PyObject* args); 
  PyObject* maskXRay(PyObject* self, PyObject* args);
  PyObject* getIntersectingDomainsAABB(PyObject* self, PyObject* args);
  PyObject* applyBCOverlapStruct(PyObject* self, PyObject* args);
  PyObject* applyBCOverlapsNG(PyObject* self, PyObject* args);
  PyObject* setDoublyDefinedBC(PyObject* self, PyObject* args);
  PyObject* getOversetHolesInterpCellCenters(PyObject* self, PyObject* args);
  PyObject* getOversetHolesInterpNodes(PyObject* self, PyObject* args);
  PyObject* _getOversetHolesInterpCellCenters(PyObject* self, PyObject* args);
  PyObject* _getOversetHolesInterpNodes(PyObject* self, PyObject* args);
  PyObject* getEXPoints(PyObject* self, PyObject* args);
  PyObject* setInterpolations(PyObject* self, PyObject* args);
  PyObject* setInterpData(PyObject* self, PyObject* args);
  PyObject* setInterpDataDW(PyObject* self, PyObject* args);
  PyObject* setInterpDataForGC(PyObject* self, PyObject* args);
  PyObject* setInterpDataLS(PyObject* self, PyObject* args);
  PyObject* setInterpDataCons(PyObject* self, PyObject* args);
  PyObject* writeCoefs(PyObject* self, PyObject* args);
  PyObject* chimeraTransfer(PyObject* self, PyObject* args);
  PyObject* transferFields(PyObject* self, PyObject* args);
  PyObject* initNuma(PyObject* self, PyObject* args);
  PyObject* setInterpTransfers(PyObject* self, PyObject* args);// en pratique non appelee de PyTree
  PyObject* _setInterpTransfers(PyObject* self, PyObject* args);
  PyObject* __setInterpTransfers(PyObject* self, PyObject* args);
  PyObject* ___setInterpTransfers(PyObject* self, PyObject* args);
  PyObject* setInterpTransfersD(PyObject* self, PyObject* args);// en pratique non appelee de PyTree
  PyObject* _setInterpTransfersD(PyObject* self, PyObject* args);
  PyObject* __setInterpTransfersD(PyObject* self, PyObject* args);
  PyObject* getInterpolatedPoints(PyObject* self, PyObject* args);
  PyObject* getInterpolatedPointsZ(PyObject* self, PyObject* args);
  PyObject* changeWall(PyObject* self, PyObject* args);
  PyObject* changeWallEX(PyObject* self, PyObject* args);
  PyObject* modifyBorders(PyObject* self, PyObject* args);
  PyObject* blankIntersectingCells(PyObject* self, PyObject* args);
  PyObject* cellN2OversetHolesStruct(PyObject* self, PyObject* args);
  PyObject* cellN2OversetHolesUnStruct(PyObject* self, PyObject* args);
  PyObject* identifyMatching(PyObject* self, PyObject* args);
  PyObject* identifyMatchingP(PyObject* self, PyObject* args);
  PyObject* identifyMatchingNM(PyObject* self, PyObject* args);
  PyObject* identifyDegenerated(PyObject* self, PyObject* args);
  PyObject* gatherMatching(PyObject* self, PyObject* args);
  PyObject* gatherMatchingNM(PyObject* self, PyObject* args);
  PyObject* gatherMatchingNGon(PyObject* self, PyObject* args);
  PyObject* _getEmptyBCInfoNGON(PyObject* self, PyObject* args);
  PyObject* gatherDegenerated(PyObject* self, PyObject* args);
  PyObject* setIBCTransfers(PyObject* self, PyObject* args);
  PyObject* _setIBCTransfers(PyObject* self, PyObject* args);
  PyObject* setIBCTransfersD(PyObject* self, PyObject* args);
  PyObject* _setIBCTransfersD(PyObject* self, PyObject* args);
  PyObject* getExtrapAbsCoefs(PyObject* self, PyObject* args);
  PyObject* _updateNatureForIBM(PyObject* self, PyObject* args);//on a zone, in place
  PyObject* indiceToCoord2(PyObject* self, PyObject* args);//on a zone, in place
  PyObject* correctCoeffList(PyObject* self, PyObject* args);//on a zone, in place
  PyObject* _blankClosestTargetCells(PyObject* self, PyObject* args);
}
#endif
