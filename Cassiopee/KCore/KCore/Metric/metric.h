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

#ifndef _KCORE_METRIC_H
#define _KCORE_METRIC_H
# include "Def/DefTypes.h"
# include "Fld/FldArray.h"
# include "Def/DefFunction.h"

namespace K_METRIC
{
  /* Calcul des centres des interfaces pour une grille structuree definie 
     en noeuds.
     IN: im, im, km: Number of mesh vertices along %i, %j, %k
     IN: inti: Total number of interfaces along the %i direction
     IN: intij: Total number of interfaces along the %i direction and %j direction
     IN: x, y, z: Vertex coordinates
     OUT: cix, ciy, ciz: interface centers
  */
  void compCenterInterface(
    const E_Int im, const E_Int jm, const E_Int km,
    const E_Int inti, const E_Int intij,
    const E_Float* x, const E_Float* y, const E_Float* z,
    E_Float* cix, E_Float* ciy, E_Float* ciz);
  
  /* Surface vector and normal to the interfaces for structured grids.
     IN: im, im, km: Number of mesh vertices along %i, %j, %k
     IN: inti: Total number of interfaces along the %i direction
     IN: intij: Total number of interfaces along the %i direction and %j direction
     IN: x, y, z: Vertex coordinates
     OUT: surfx, surfy, surfz: Surface vector
     OUT: snorm: Norm of the surface vector
  */
  void compIntSurf(
    const E_Int im, const E_Int jm, const E_Int km,
    const E_Int inti, const E_Int intij,
    const E_Float* x, const E_Float* y, const E_Float* z,
    E_Float* surfx, E_Float* surfy, E_Float* surfz,
    E_Float* snorm);

  /* Calcul des vecteurs surfaces des 6 interfaces d'une cellule 3D structuree.
     IN: ind: indice du premier sommet de la cellule
     IN: ni, nj, nk: dimensions du maillage en noeuds
     IN: xt, yt, zt: Vertex coordinates
     OUT: surf: vecteur surface de la cellule
  */
  void compIntSurfOfCell(
    const E_Int ind, const E_Int ni, const E_Int nj, const E_Int nk,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* surf);
  
  /* Calcul d'une longueur caracteristique d'une cellule structuree: 
     moyenne des longueurs des cotes de la cellule
     IN: ni, nj, nk: Number of mesh vertices along %i, %j, %k
     IN: indA: First vertex of the cell
     IN: xt, yt, zt: Vertex coordinates
     OUT: meanl: Cell mean length
  */
  void compMeanLengthOfStructCell(
    const E_Int ni, const E_Int nj, const E_Int nk, E_Int indA,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float& meanl);

  /* Calcul d'une longueur caracteristique d'une cellule tetra: 
     moyenne des longueurs des cotes de la cellule
     IN: indA, indB, indC, indD: Sommets du tetraedre
     IN: xt, yt, zt: Vertex coordinates
     OUT: meanl: Cell mean length
  */
  void compMeanLengthOfTetraCell(
    const E_Int indA, const E_Int indB, const E_Int indC, const E_Int indD,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float& meanl);

  /* Calcul d'une longueur caracteristique de la cellule : minimum  des
     longueurs des cotes de la cellule
     IN: ni, nj, nk: Number of mesh vertices along %i, %j, %k
     IN: indA: Index of the first cell vertex
     IN: xt, yt, zt: Vertex coordinates
     OUT: minl: Minimum length of the cell
  */
  void compMinLengthOfCell(
    const E_Int ni, const E_Int nj, const E_Int nk, E_Int indA,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float& minl);
  
  /* Calcul d'une longueur caracteristique de la cellule : minimum  des
     longueurs des cotes de la cellule
     IN: indA, indB, indC, indD: Sommets du tetraedre
     IN: xt, yt, zt: Vertex coordinates
     OUT: minl: Minimum length of the cell
  */
  void compMinLengthOfTetraCell(
    const E_Int indA, const E_Int indB, const E_Int indC, const E_Int indD,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float& minl);
  
  /* Calcul les volumes pour des elements NGon
     IN: (xt, yt, zt): pointeurs sur les coordonnees du maillage
     IN: cn: connectivite NGon
     OUT: volp: pointeur sur le tableau des volumes calcules aux centres des elements
  */
  E_Int compVolNGon(
    const E_Float* xt, const E_Float* yt, const E_Float* zt, 
    K_FLD::FldArrayI& cn, E_Float* volp);
  
  E_Int compNGonVolOfElement(
    const E_Float* xt, const E_Float* yt, const E_Float* zt, 
    K_FLD::FldArrayI& cn, E_Int indE, std::vector<std::vector<E_Int> > cnEV, 
    K_FLD::FldArrayI& posElt, K_FLD::FldArrayI& posFace, 
    K_FLD::FldArrayI& dimElt, E_Float& vol);

  E_Int compSurfNGon(
    const E_Float* xt, const E_Float* yt, const E_Float* zt, 
    K_FLD::FldArrayI& cn, 
    E_Float* sxp, E_Float* syp,  E_Float* szp);
  
  /* Calcul des surfaces orientees des faces et la norme associee
     On suppose que le NGON est deja correctement oriente
     IN: (xt, yt, zt): pointeurs sur les coordonnees du maillage
     IN:  cnp: pointeur sur la connectivite NGon
     OUT: sxp, syp, szp, snp: surface orientee calculee pour les faces et norme associee
     Return 0 (OK), 1 (Failed)
  */
  E_Int compNGonFacesSurf(
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    K_FLD::FldArrayI& cn,
    E_Float* sxp, E_Float* syp,  E_Float* szp, E_Float* snp, 
    K_FLD::FldArrayI* cFE=NULL);  

  /* Attention cette routine est 2D: pas de calcul des normales sur un hexaedre!
    IN: ni, nj: Number of mesh vertices along %i, %j
    IN: xt, yt, zt: Vertex coordinates
    OUT: nxt, nyt, nzt: Norm of the surface vector
  */
  void compNormStructSurf(
    const E_Int ni, const E_Int nj,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* nxt, E_Float* nyt, E_Float* nzt);

  /* Calcul les vecteurs normaux aux triangles
     IN: cn: connectivite elts-noeuds
     IN: xt, yt, zt: coordonnees x, y, z des pts de la grille
     OUT: nsurf: normales aux facettes
  */
  void compNormUnstructSurf(
    K_FLD::FldArrayI& cn, const char* eltType,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* nsurfx, E_Float* nsurfy, E_Float* nsurfz);

  /* Compute barycenter of cells.
     This version does not use minOfIntvoid
     IN: im, im, km: Number of mesh vertices along %i, %j, %k
     IN: xt, yt, zt: Vertex coordinates
     OUT: bary: Cell centers
  */
  void compStructCellCenter(
    const E_Int im, const E_Int jm, const E_Int km,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* bary);

  /* Calcul du volume de toutes les cellules et des surface des interfaces.
     CAS STRUCTURE
     IN: ni, nj, nk: Number of mesh vertices along %i, %j, %k
     IN: nbInti, nbIntj, nbIntk: Number of interfaces along the %i, %j, %k directions
     IN: x, y, z: Vertex coordinates
     OUT: vol: Volume of the elements
     OUT: surfx, surfy, surfz: Surface vectors
     OUT: snorm: Norm of the surface vectors
     OUT: cix, ciy, ciz: Interface centers
  */
  void compMetricStruct(
    const E_Int ni, const E_Int nj, const E_Int nk,
    const E_Int nbInti, const E_Int nbIntj, const E_Int nbIntk,
    const E_Float* x, const E_Float* y, const E_Float* z,
    E_Float* vol, E_Float* surfx, E_Float* surfy, E_Float* surfz,
    E_Float* snorm, E_Float* cix, E_Float* ciy, E_Float* ciz);

  /* Calcul de la surface pour une grille surfacique structuree.
     IN: ni, nj, nk: Number of mesh vertices along %i, %j, %k
     IN: xt, yt, zt: Vertex coordinates
     OUT: surface: Mesh area
  */
  void compSurfStruct2D(
    const E_Int ni, const E_Int nj, const E_Int nk,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* surface);

  /* Calcul de la longueur entre chaque sommet pour une ligne structuree.
     IN: ni, nj, nk: Number of mesh vertices along %i, %j, %k
     IN: xt, yt, zt: Vertex coordinates
     OUT: length: Longueur entre chaque sommet
  */
  void compSurfStruct1D(
    const E_Int ni, const E_Int nj, const E_Int nk,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* length);

  /* Calcul de la surface pour une grille surfacique structuree.
     IN: ni, nj, nk: Number of mesh vertices along %i, %j, %k
     IN: xt, yt, zt: Vertex coordinates
     OUT: surface: Aire
  */
  void compSurfOfStructCell(
    const E_Int ni, const E_Int nj, const E_Int nk, const E_Int indcell,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float& surface);

  //============================================================================
  // Compute cell barycenter for a ME mesh.
  // IN: cn: Element-Node connectivity
  // IN: xt, yt, zt: Vertex coordinates
  // OUT: xb, yb, zb: Barycenter coordinates
  //============================================================================
  void compUnstructCellCenter(
    K_FLD::FldArrayI& cn,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* xb, E_Float* yb, E_Float* zb);

  /* Calcul des centres des interfaces pour des mailles non structurees.
     IN: nedges: nb de facettes par elemt
     IN: cn: Element-Node connectivity
     IN: xt, yt, zt: Vertex coordinates
     OUT: xint, yint, zint: Coordonnees du centre des facettes
  */
  void compUnstructCenterInt(
    K_FLD::FldArrayI& cn, const char* eltType,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* xint, E_Float* yint, E_Float* zint);

  /* Calcul du volume de toutes les cellules et des surfaces des interfaces
     CAS NON STRUCTURE
     IN: npts: nb de pts du maillage
     IN: nelts: nb d elements
     IN: nedges: nb de facettes par elemt
     IN: nnodes: nb de noeuds par elemt
     IN: cn: nb de noeuds par elemt
     IN: coordx, coordy, coordz: coordonnees x, y, z des pts de la grille
     OUT: snx, sny, snz: normales aux facettes %x, %y, %z
     OUT: surf: aires des facettes
     OUT: vol: volume des cellules
  */
  void compMetricUnstruct(
    K_FLD::FldArrayI& cn, const char* eltType,
    const E_Float* coordx, const E_Float* coordy, const E_Float* coordz,
    E_Float* snx, E_Float* sny, E_Float* snz, E_Float* surf, E_Float* vol);

  /* Calcul des normales pour un maillage BE.
     Les normales aux surfaces sont orientees vers l'exterieur de l'element.
     IN: xt, yt, zt: pointeurs sur les coordonnees du maillage
  */
  void compSurfUnstruct(
    K_FLD::FldArrayI& cn, const char* eltType,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* surfnx, E_Float* surfny, E_Float* surfnz, E_Float* surface);

  // depend du type d'element considere 
  void compTriSurf(
    K_FLD::FldArrayI& cm, const E_Int fcOffset,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* surfnx, E_Float* surfny, E_Float* surfnz, E_Float* surface);
  void compQuadSurf(
    K_FLD::FldArrayI& cm, const E_Int fcOffset,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* surfnx, E_Float* surfny, E_Float* surfnz, E_Float* surface);
  void compTetraSurf(
    K_FLD::FldArrayI& cm, const E_Int fcOffset,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* surfnx, E_Float* surfny, E_Float* surfnz, E_Float* surface);
  void compPyraSurf(
    K_FLD::FldArrayI& cm, const E_Int fcOffset,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* surfnx, E_Float* surfny, E_Float* surfnz, E_Float* surface);
  void compPentaSurf(
    K_FLD::FldArrayI& cm, const E_Int fcOffset,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* surfnx, E_Float* surfny, E_Float* surfnz, E_Float* surface);
  void compHexaSurf(
    K_FLD::FldArrayI& cm, const E_Int fcOffset,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* surfnx, E_Float* surfny, E_Float* surfnz, E_Float* surface);

  /* Calcul de la longueur entre chaque sommet pour une ligne non structuree

  */
  void compUnstructSurf1d(
    K_FLD::FldArrayI& cn, const char* eltType,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float* length);

  /* Calcule l aire d'une cellule d un maillage surfacique nk=1.
     N'est pas necessairement dans le plan 
     On rentre soit l indice de la cellule indcell, soit indnode l indice du premier point
     d indices i et j min de la cellule. Si indnode est different de -1, c'est lui qui prime
  */
  void compVolOfStructCell2D(
    const E_Int ni, const E_Int nj,
    const E_Int indcell, E_Int indnode,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float& area);

  void compVolOfStructCell3D(
    const E_Int ni, const E_Int nj, const E_Int nk,
    const E_Int indcell, E_Int indnode,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float& vol);

  //============================================================================
  // Calcul du volume des elements pour un maillage Multi-Elements.
  // IN: cn: Element-Node connectivity
  // IN: eltType: Element names
  // IN: xint, yint, zint: coordonnees du centre des facettes
  // IN: snx, sny, snz: normales aux facettes %x, %y, %z
  // OUT: vol: volume des cellules
  //============================================================================
  void compUnstructVol(
    K_FLD::FldArrayI& cn, const char* eltType,
    const E_Float* xint, const E_Float* yint, const E_Float* zint,
    const E_Float* snx, const E_Float* sny, const E_Float* snz, E_Float* vol);

  // Compute cell volumes for NGons
  void compute_face_center_and_area(E_Int id, E_Int stride, E_Int *pn,
    E_Float *x, E_Float *y, E_Float *z, E_Float *fc, E_Float *fa);
  
  E_Int compVolNGonImad(E_Float *x, E_Float *y, E_Float *z,
    K_FLD::FldArrayI &cn, E_Float *cellVols);

  void compute_cell_volume(
    E_Int, K_FLD::FldArrayI &,
    E_Int *, E_Int *, E_Int *, E_Int *,
    E_Float *, E_Float *, E_Float *, E_Float &, E_Int refIdx=0
  );
  
  void compute_face_centers_and_areas(K_FLD::FldArrayI &cn, E_Float *x,
    E_Float *y, E_Float *z, E_Float *fcenters, E_Float *fareas);
  
  void compute_cell_centers_and_vols(K_FLD::FldArrayI &cn, E_Float *x,
    E_Float *y, E_Float *z, E_Int *owner, E_Int *neigh, E_Float *fcenters,
    E_Float *fareas, E_Float *cx, E_Float *cy, E_Float *cz, E_Float *volumes);
  
  E_Int compute_gradients_ngon(K_FLD::FldArrayI &cn, E_Float *x, E_Float *y,
    E_Float *z, E_Int *owner, E_Int *neigh, E_Float *centers,
    const std::vector<E_Float *> &flds, std::vector<E_Float *> &Gs);
  
  void compute_cell_center_and_volume(E_Int id, E_Int stride,
    E_Int *pf, E_Float *x, E_Float *y, E_Float *z, E_Float *fc, E_Float *fa,
    E_Int *owner, E_Float &cx, E_Float &cy, E_Float &cz, E_Float &vol);
}
#endif
