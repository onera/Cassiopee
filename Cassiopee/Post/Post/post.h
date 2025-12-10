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
#ifndef _POST_POST_H_
#define _POST_POST_H_

# include <locale>
# include <cctype>
# include "kcore.h"
# include <vector>

namespace K_POST
{
  PyObject* extractPoint(PyObject* self, PyObject* args);
  PyObject* extractPlane(PyObject* self, PyObject* args);
  PyObject* extractMesh(PyObject* self, PyObject* args);
  PyObject* projectCloudSolution2Triangle(PyObject* self, PyObject* args);
  PyObject* prepareProjectCloudSolution2Triangle(PyObject* self, PyObject* args);
  PyObject* projectCloudSolution2TriangleWithInterpData(PyObject* self, PyObject* args);
  PyObject* coarsen(PyObject* self, PyObject* args);
  PyObject* refine(PyObject* self, PyObject* args);
  PyObject* refineButterfly(PyObject* self, PyObject* args);
  PyObject* selectCells(PyObject* self, PyObject* args);
  PyObject* selectCellsBoth(PyObject* self, PyObject* args);
  PyObject* selectCells3(PyObject* self, PyObject* args);
  PyObject* selectCellCenters(PyObject* self, PyObject* args);
  PyObject* selectCellCentersBoth(PyObject* self, PyObject* args);
  PyObject* selectInteriorFaces(PyObject* self, PyObject* args);
  PyObject* selectExteriorFaces(PyObject* self, PyObject* args);
  PyObject* selectExteriorFacesStructured(PyObject* self, PyObject* args);
  PyObject* selectExteriorElts(PyObject* self, PyObject* args);
  PyObject* exteriorEltsStructured(PyObject* self, PyObject* args);
  PyObject* frontFaces(PyObject* self, PyObject* args);
  PyObject* integ(PyObject* self, PyObject* args);
  PyObject* integ2(PyObject* self, PyObject* args);
  PyObject* integVect(PyObject* self, PyObject* args);
  PyObject* integNorm(PyObject* self, PyObject* args);
  PyObject* integNormProduct(PyObject* self, PyObject* args);
  PyObject* integMoment(PyObject* self, PyObject* args);
  PyObject* integMomentNorm(PyObject* self, PyObject* args);
  PyObject* zipperF(PyObject* self, PyObject* args);
  PyObject* usurpF(PyObject* self, PyObject* args);
  PyObject* computeVariables(PyObject* self,PyObject* args);
  PyObject* computeVariables2(PyObject* self,PyObject* args);
  PyObject* computeGrad(PyObject* self,PyObject* args);
  PyObject* computeGrad2NGon(PyObject* self,PyObject* args);
  PyObject* computeGrad2Struct(PyObject* self,PyObject* args);
  PyObject* computeGradLSQ(PyObject *self, PyObject *args);
  PyObject* computeNormGrad(PyObject* self, PyObject* args);
  PyObject* computeDiv(PyObject* self, PyObject* args);
  PyObject* computeDiv2NGon(PyObject* self, PyObject* args);
  PyObject* computeDiv2Struct(PyObject* self, PyObject* args);
  PyObject* computeCurl(PyObject* self, PyObject* args);
  PyObject* computeNormCurl(PyObject* self, PyObject* args);
  PyObject* computeDiff(PyObject* self, PyObject* args);
  PyObject* perlinNoise(PyObject* self, PyObject* args);
  PyObject* compStreamLine(PyObject* self, PyObject* args);
  PyObject* comp_stream_line(PyObject* self, PyObject* args);
  PyObject* compStreamRibbon(PyObject* self, PyObject* args);
  PyObject* comp_streamribbon(PyObject* self, PyObject* args);
  PyObject* compStreamSurf(PyObject* self, PyObject* args);
  PyObject* isoLine(PyObject* self, PyObject* args);
  PyObject* isoSurf(PyObject* self, PyObject* args);
  PyObject* isoSurfMC(PyObject* self, PyObject* args);
  PyObject* isoSurfNGon(PyObject* self, PyObject* args);
  PyObject* computeIndicatorValue(PyObject* self, PyObject* args);
  PyObject* enforceIndicatorNearBodies(PyObject* self, PyObject* args);
  PyObject* enforceIndicatorForFinestLevel(PyObject* self, PyObject* args);
  PyObject* enforceIndicatorForCoarsestLevel(PyObject* self, PyObject* args);
  PyObject* sharpEdges(PyObject* self, PyObject* args);
  PyObject* silhouette(PyObject* self, PyObject* args);

// ******************************* INTEGRATION ****************************** //
// ============================================================================
// Compute surface integral of the field F, coordinates 
// and field have the same size
//   I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
//   Aire(ABCD) = ||AB^AC||/2+||DB^DC||/2
// ============================================================================
  void integStructNodeCenter2D(
    const E_Int ni, const E_Int nj,
    const E_Float* ratio, const E_Float* surf, const E_Float* field,
    E_Float& result);
// ============================================================================
// Compute linear integral of field F, structured 1D case
//   I(AB) = Length(AB)*(F(A)+F(B))/2
// ============================================================================
  void integStructNodeCenter1D(
    const E_Int ni,
    const E_Float* ratio, const E_Float* length, const E_Float* field,
    E_Float& result);
// ============================================================================
// Compute surface integral of field F, coordinates in nodes,
// field defined in centers, structured case
// IN: ni1, nj1 : dim en centres
// ============================================================================
  void integStructCellCenter2D(
    const E_Int ni1, const E_Int nj1,
    const E_Float* ratio, const E_Float* surf, const E_Float* field,
    E_Float& result);

// ============================================================================
// Compute linear integral of field F, node/center case, 1D
// ============================================================================
  void integStructCellCenter1D(
    const E_Int ni,
    const E_Float* ratio, const E_Float* length, const E_Float* field,
    E_Float& result);

// ============================================================================
// Compute surface integral of the field F, unstructured case
// Coordinates and field have the same size
//   I(ABC) = Aire(ABC)*(F(A)+F(B)+F(C))/3        TRI
//   I(ABC) = Aire(ABC)*(F(A)+F(B)+F(C)+F(D))/4   QUAD
//   Aire(ABC) = ||AB^AC||/2
// ============================================================================
  void integUnstructNodeCenter(
    FldArrayI& cn,
    const E_Float* ratio, const E_Float* surf, const E_Float* field,
    E_Float& result);

// ============================================================================
// Compute surface integral of field F, coordinates in nodes,
// field defined in centers, unstructured case
// ============================================================================
  void integUnstructCellCenter(
    const E_Int nelts, const E_Float* ratio, 
    const E_Float* surf, const E_Float* field,
    E_Float& result);

/*
  Compute the surface integral of field F
  posx, posy, posz : positions de x,y,z dans coord
  ni*nj*nk: dimension of coordinate array (coord)
  ratio: to have the better value in overlap mesh case
  (equal to 1 if not defined)
  center2node: 0 if coord and F have the same size
  1 if coord is in nodes and F in centers
  resultat: integration result, same size as F variable number
*/
  E_Int integStruct2D(E_Int ni, E_Int nj, E_Int nk,
                      E_Int center2node, E_Int posx, E_Int posy, E_Int posz,
                      FldArrayF& coord, FldArrayF& F,
                      FldArrayF& ratio, FldArrayF& resultat);
/*
  Compute the linear integral of field F
  ni*nj*nk: dimension of coordinate array (coord)
  posx, posy, posz: positions de x,y,z dans coord
  ratio: to have the better value in overlap mesh case
  (equal to 1 if not defined)
  center2node: 0 if coord and F have the same size
  1 if coord is in nodes and F in centers
  resultat: integration result, same size as F variable number
*/
  E_Int integStruct1D(E_Int ni, E_Int nj, E_Int nk,
                      E_Int center2node, E_Int posx, E_Int posy, E_Int posz,
                      FldArrayF& coord, FldArrayF& F,
                      FldArrayF& ratio, FldArrayF& resultat);
/*
  Compute the surface integral of field F
  posx, posy, posz: positions de x,y,z dans coord
  cn: connectivite: element ->noeud
  coord: coordonnees  des noeuds
  ratio: ponderation ds integ. Taille : celle de xt
  F: champ a integrer
  center2node : 0 if coord and F have the same size
  1 if coord is in nodes and F in centers
  resultat: integration result, same size as F variable number
*/
  E_Int integUnstruct2D(E_Int center2node,
                        E_Int posx, E_Int posy, E_Int posz,
                        FldArrayI& cn, const char* eltType, FldArrayF& coord,
                        FldArrayF& F, FldArrayF& ratio,
                        FldArrayF& resultat);
/*
  Calcul de l'integrale de F sur une surface "BAR"
  posx, posy, posz: positions de x,y,z dans coord
  cn: connectivite: element ->noeud
  coord: coordonnees  des noeuds
  ratio: ponderation ds integ. Taille: celle de xt
  F: champ a integrer
  center2node: 0 if coord and F have the same size
  1 if coord is in nodes and F in centers
  resultat: integration result, same size as F variable number
*/
  E_Int integUnstruct1D(E_Int center2node,
                        E_Int posx, E_Int posy, E_Int posz,
                        FldArrayI& cn, const char* eltType, FldArrayF& coord,
                        FldArrayF& F, FldArrayF& ratio,
                        FldArrayF& resultat);
/*
   Compute the surface integral of field F * normal vect(n)
   ni*nj*nk: dimension of coordinate array (coord)
   ratio: to have the better value in overlap mesh case
   (equal to 1 if not defined)
   center2node: 0 if coord and F have the same size
   1 if coord is in nodes and F in centers
   resultat: integration result, same size as F variable number*3
*/
  E_Int integNormStruct2D(E_Int ni, E_Int nj, E_Int nk,
                          E_Int center2node,
                          E_Int posx, E_Int posy, E_Int posz,
                          FldArrayF& coord, FldArrayF& F,
                          FldArrayF& ratio, FldArrayF& resultat);

/*
  Compute the surface integral of field F * normal vect(n)
  cn: connectivity
  ratio: to have the better value in overlap mesh case
  (equal to 1 if not defined)
  center2node: 0 if coord and F have the same size
  1 if coord is in nodes and F in centers
  resultat: integration result, same size as F variable number*3
*/
  E_Int integNormUnstruct2D(E_Int center2node,
                            E_Int posx, E_Int posy, E_Int posz,
                            FldArrayI& cn, const char* eltType, FldArrayF& coord,
                            FldArrayF& F, FldArrayF& ratio,
                            FldArrayF& resultat);

/*
   Compute the surface integral of scalar product field
   vect(F) * normal vect(n)
   ni*nj*nk: dimension of coordinate array (coord)
   ratio: to have the better value in overlap mesh case
   (equal to 1 if not defined)
   center2node: 0 if coord and F have the same size
   1 if coord is in nodes and F in centers
   resultat: integration result, same size as F variable number
*/
  E_Int integNormProdStruct2D(E_Int ni, E_Int nj, E_Int nk,
                              E_Int center2node,
                              E_Int posx, E_Int posy, E_Int posz,
                              FldArrayF& coord, FldArrayF& F,
                              FldArrayF& ratio, E_Float& resultat);

/*
   Compute the surface integral of scalar product field
   vect(F) * normal vect(n)
   cn: connectivity
   ratio: to have the better value in overlap mesh case
   (equal to 1 if not defined)
   center2node: 0 if coord and F have the same size
   1 if coord is in nodes and F in centers
   resultat: integration result, same size as F variable number
*/
  E_Int integNormProdUnstruct2D(E_Int center2node,
                               E_Int posx, E_Int posy, E_Int posz,
                               FldArrayI& cn, const char* eltType, FldArrayF& coord,
                               FldArrayF& F, FldArrayF& ratio,
                               E_Float& resultat);
/*
   Compute the surface integral of moment (OM^F)
   ni*nj*nk: dimension of coordinate array (coord)
   ratio: to have the better value in overlap mesh case
   (equal to 1 if not defined)
   center2node: 0 if coord and F have the same size
   1 if coord is in nodes and F in centers
   resultat: integration result, same size as F variable number
*/
  E_Int integMomentStruct2D(E_Int ni, E_Int nj, E_Int nk,
                            E_Int center2node,
                            E_Int posx, E_Int posy, E_Int posz,
                            E_Float cx, E_Float cy, E_Float cz,
                            FldArrayF& coord, FldArrayF& F,
                            FldArrayF& ratio, FldArrayF& resultat);
/*
   Compute the linear integral of moment (OM^F)
   ni*nj*nk: dimension of coordinate array (coord)
   ratio: to have the better value in overlap mesh case
   (equal to 1 if not defined)
   center2node: 0 if coord and F have the same size
   1 if coord is in nodes and F in centers
   resultat: integration result, same size as F variable number
*/
  E_Int integMomentStruct1D(E_Int ni, E_Int nj, E_Int nk,
                            E_Int center2node,
                            E_Int posx, E_Int posy, E_Int posz,
                            E_Float cx, E_Float cy, E_Float cz,
                            FldArrayF& coord, FldArrayF& F,
                            FldArrayF& ratio, FldArrayF& resultat);
/*
   Compute the surface integral of moment (OM^F)
   cn: connectivity
   ratio: to have the better value in overlap mesh case
   (equal to 1 if not defined)
   center2node : 0 if coord and F have the same size
   1 if coord is in nodes and F in centers
   resultat: integration result, same size as F variable number
*/
  E_Int integMomentUnstruct2D(E_Int center2node,
                              E_Int posx, E_Int posy, E_Int posz,
                              E_Float cx, E_Float cy, E_Float cz,
                              FldArrayI& cn, const char* eltType, FldArrayF& coord,
                              FldArrayF& F, FldArrayF& ratio,
                              FldArrayF& resultat);
/*
   Compute the linear integral of moment (OM^F)
   cn: connectivity
   ratio: to have the better value in overlap mesh case
   (equal to 1 if not defined)
   center2node: 0 if coord and F have the same size
   1 if coord is in nodes and F in centers
   resultat: integration result, same size as F variable number
*/
  E_Int integMomentUnstruct1D(E_Int center2node,
                              E_Int posx, E_Int posy, E_Int posz,
                              E_Float cx, E_Float cy, E_Float cz,
                              FldArrayI& cn, const char* eltType, FldArrayF& coord,
                              FldArrayF& F, FldArrayF& ratio,
                              FldArrayF& resultat);
/*
   Compute the surface integral of moment (OM^F.vect(n))
   ni*nj*nk: dimension of coordinate array (coord)
   ratio: to have the better value in overlap mesh case
   (equal to 1 if not defined)
   center2node: 0 if coord and F have the same size
   1 if coord is in nodes and F in centers
   resultat: integration result, same size as F variable number
*/
  E_Int integMomentNormStruct2D(E_Int ni, E_Int nj, E_Int nk,
                                E_Int center2node,
                                E_Int posx, E_Int posy, E_Int posz,
                                E_Float cx, E_Float cy, E_Float cz,
                                FldArrayF& coord,
                                FldArrayF& F, FldArrayF& ratio,
                                FldArrayF& resultat);
/*
   Compute the surface integral of moment (OM^F.vect(n))
   cn: connectivity
   ratio: to have the better value in overlap mesh case
   (equal to 1 if not defined)
   center2node: 0 if coord and F have the same size
   1 if coord is in nodes and F in centers
   resultat: integration result, same size as F variable number
*/
  E_Int integMomentNormUnstruct2D(E_Int center2node,
                                  E_Int posx, E_Int posy, E_Int posz,
                                  E_Float cx, E_Float cy, E_Float cz,
                                  FldArrayI& cn, const char* eltType, FldArrayF& coord,
                                  FldArrayF& F, FldArrayF& ratio,
                                  FldArrayF& resultat);

// ============================================================================
// Compute surface integral of the moment M (OM^F), coordinates 
//     and field have the same size
//     I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
//     Aire(ABCD) = ||AB^AC||/2+||DB^DC||/2
// ============================================================================
  void integMomentStructNodeCenter2D(
    const E_Int ni, const E_Int nj,
    const E_Float cx, const E_Float cy, const E_Float cz,
    const E_Float* ratio, const E_Float* xt, const E_Float* yt,
    const E_Float* zt, const E_Float* surf,
    const E_Float* vx, const E_Float* vy, const E_Float* vz,
    E_Float* result);

//=============================================================================
// Compute linear integral of the moment M (OM^F), coordinates 
//     and field have the same size
//     I(AB) = LENGTH(ABCD)*(F(A)+F(B))/2
// ============================================================================
  void integMomentStructNodeCenter1D(
    const E_Int ni, const E_Float cx, const E_Float cy,
    const E_Float cz, const E_Float* ratio,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* length,
    const E_Float* vx, const E_Float* vy, const E_Float* vz,
    E_Float* result);

//=============================================================================
// Compute surface integral of the moment M (OM^F), coordinates 
//     are defined in nodes and F is defined in center
//=============================================================================
  void integMomentStructCellCenter2D(
    const E_Int ni, const E_Int nj,
    const E_Float cx, const E_Float cy, const E_Float cz,
    const E_Float* ratio, const E_Float* xt, const E_Float* yt,
    const E_Float* zt, const E_Float* surf, const E_Float* vx,
    const E_Float* vy, const E_Float* vz, E_Float* result);

//=============================================================================
// Compute linear integral of the moment M (OM^F), coordinates 
//     are defined in nodes and F is defined in center
//=============================================================================
  void integMomentStructCellCenter1D(
    const E_Int ni, const E_Float cx, const E_Float cy, const E_Float cz,
    const E_Float* ratio, const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* length, const E_Float* vx, const E_Float* vy, const E_Float* vz,
    E_Float* result);

// ============================================================================
// Compute surface integral of the moment M (OM^F.vect(n)), coordinates 
//       and F have the same size
// ============================================================================
  void integMomentNormStructNodeCenter2D(
    const E_Int ni, const E_Int nj, const E_Float cx,
    const E_Float cy, const E_Float cz, const E_Float* ratio, const E_Float* xt,
    const E_Float* yt, const E_Float* zt, const E_Float* sx, const E_Float* sy,
    const E_Float* sz, const E_Float* field, E_Float* result);
  
//=============================================================================
// Compute surface integral of the moment M (OM^F.vect(n)), coordinates 
//     are defined in nodes and F is defined in center
//=============================================================================
  void integMomentNormStructCellCenter2D(
    const E_Int ni, const E_Int nj,
    const E_Float cx, const E_Float cy, const E_Float cz, const E_Float* ratio,
    const E_Float* xt, const E_Float* yt, const E_Float* zt, const E_Float* sx,
    const E_Float* sy, const E_Float* sz, const E_Float* field, E_Float* result);

// ============================================================================
// Compute surface integral of the field F.vect(n), coordinates 
//     and field have the same size
//     I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
//     Aire(ABCD) = ||AB^AC||/2+||DB^DC||/2
// ============================================================================
  void integNormProdStructNodeCenter2D(
    const E_Int ni, const E_Int nj, const E_Float* ratio,
    const E_Float* sx, const E_Float* sy, const E_Float* sz,
    const E_Float* vx, const E_Float* vy, const E_Float* vz,
    E_Float& result);
  
//==============================================================================
// Compute surface integral of the product vect(F).vect(n), coordinates 
// are defined in nodes and F is defined in nodes (center-based formulation)
//==============================================================================
  void integNormProdStructCellCenter2D(
    const E_Int ni, const E_Int nj, const E_Float* ratio,
    const E_Float* sx, const E_Float* sy, const E_Float* sz,
    const E_Float* vx, const E_Float* vy, const E_Float* vz,
    E_Float& result);

// ============================================================================
// Compute surface integral of the field F.vect(n), coordinates 
//      and field have the same size
//      I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
//      Aire(ABCD) = ||AB^AC||/2+||DB^DC||/2
//=============================================================================
  void integNormStructNodeCenter2D(
    const E_Int ni, const E_Int nj, const E_Float* ratio,
    const E_Float* sx, const E_Float* sy, const E_Float* sz,
    const E_Float* field, E_Float* result);

//=============================================================================
// Compute surface integral of the field F.vect(n), coordinates 
//     are defined in nodes and F is defined in center
//=============================================================================
  void integNormStructCellCenter2D(
    const E_Int ni, const E_Int nj,
    const E_Float* ratio, const E_Float* sx, const E_Float* sy,
    const E_Float* sz, const E_Float* field, E_Float* result);

// ============================================================================
// Compute surface integral of the moment M (OM^F.vect(n)), coordinates 
// and F have the same size
// ============================================================================
  void integMomentNormUnstructNodeCenter(
      FldArrayI& cn,
      const E_Float cx, const E_Float cy, const E_Float cz, 
      const E_Float* ratio,
      const E_Float* xt, const E_Float* yt, const E_Float* zt,
      const E_Float* sx, const E_Float* sy, const E_Float* sz, 
      const E_Float* field, E_Float* result);

// ============================================================================
// Compute linear integral of the moment.norm (OM^F.n), coordinates 
// are defined in nodes and F is defined in center, unstructured case
// ============================================================================
  void integMomentNormUnstructCellCenter(
    FldArrayI& cn,
    const E_Float cx, const E_Float cy, const E_Float cz,
    const E_Float* ratio,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* sx, const E_Float* sy, const E_Float* sz, 
    const E_Float* field, E_Float* result);

// ============================================================================
// Compute surface integral of the field F.vect(n), coordinates 
// and field have the same size
// I(ABC) = Aire(ABC)*(F(A)+F(B)+F(C))/3        - TRI
// I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)F(D))/4  - QUAD
// Aire(ABCD) = ||AB^AC||/2
// ============================================================================
  void integNormUnstructNodeCenter(
    FldArrayI& cn,
    const E_Float *ratio,
    const E_Float *sx, const E_Float *sy, const E_Float *sz, 
    const E_Float *field, E_Float *result);

// ============================================================================
// Compute surface integral of the field F, coordinates are defined
// in nodes and F is defined in center, unstructured case
// ============================================================================
  void integNormUnstructCellCenter(
    const E_Int nelts, const E_Float *ratio,
    const E_Float *nsurfx, const E_Float *nsurfy, const E_Float *nsurfz,
    const E_Float *field, E_Float *result);

// ============================================================================
// Compute surface integral of the moment M (OM^F), coordinates 
// and field have the same size
// I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
// Aire(ABCD) = ||AB^AC||/2 + ||DB^DC||/2
// ============================================================================
  void integMomentUnstructNodeCenter(
    FldArrayI& cn,
    const E_Float cx, const E_Float cy, const E_Float cz, const E_Float* ratio,
    const E_Float* xt, const E_Float* yt, const E_Float* zt, const E_Float* surf,
    const E_Float* vx, const E_Float* vy, const E_Float* vz, E_Float* result);

// ============================================================================
// Compute surface integral of the moment M (OM^F)
// coordinates are defined in nodes and F is defined in center (unstructured)
// ============================================================================
  void integMomentUnstructCellCenter(
    FldArrayI& cn,
    const E_Float cx, const E_Float cy, const E_Float cz,
    const E_Float* ratio, const E_Float* xt, const E_Float* yt,
    const E_Float* zt, const E_Float* surf,
    const E_Float* vx, const E_Float* vy, const E_Float* vz,
    E_Float* result);

// ============================================================================
// Compute surface integral of the field F.vect(n), coordinates 
//     and field have the same size
//     I(ABC) = Aire(ABC) * (F(A) + F(B) + F(C)) / 3
//     Aire(ABC) = ||AB ^ AC|| / 2
// ============================================================================
  void integNormProdUnstructNodeCenter(
    FldArrayI& cn,
    const E_Float* ratio,
    const E_Float* sx, const E_Float* sy, const E_Float* sz,
    const E_Float* vx, const E_Float* vy, const E_Float* vz,
    E_Float& result);

// ============================================================================
// Compute surface integral of the field F, coordinates are defined
// in nodes and F is defined in center, unstructured case
// ============================================================================
  void integNormProdUnstructCellCenter(
    const E_Int nbt, const E_Float* ratio,
    const E_Float* sx, const E_Float* sy, const E_Float* sz,
    const E_Float* vx, const E_Float* vy, const E_Float* vz,
    E_Float& result);

// ============================================================================
// Etant donnes n champs definis aux noeuds d une grille 3D, 
// calcul des champs aux centres des interfaces de la grille  
// fint = 0.25*(fa+fb+fc+fd)
// ============================================================================
  void compIntField(
    const E_Int ni, const E_Int nj, const E_Int nk,
    const E_Int nfld, const E_Float* f,
    E_Float* fint);

// ============================================================================
// Compute the values of the vector (f1,f2,f3) at interfaces
// ============================================================================
  void compIntFieldV(
    const E_Int ni, const E_Int nj, const E_Int nk,
    const E_Float* f1, const E_Float* f2, const E_Float* f3,
    E_Float* fint1, E_Float* fint2, E_Float* fint3);


// ***************************** COMPUTE FIELDS ***************************** //
/* Extrait de la chaine vars0 les variables a calculer. Une verification
   est effectuee sur les noms de variables. La chaine varStringOut est
   aussi construite pour l array de sortie
   IN: vars0: chaine contenant les variables a extraire
   OUT: vars: vecteur contenant les variables a calculer
   OUT: varStringOut: chaine de variables calculees pour l array de sortie
   retourne 0 si aucune variable n'a ete trouvee
*/
  short checkAndExtractVariables(char* vars0, std::vector<char*>& vars,
                                 char* varStringOut);

/* Calcule les variables composees a partir des grandeurs conservatives.
 * Retourne 0 si erreur.*/
  E_Int computeCompVariables(const FldArrayF& f, const E_Int posro,
                             const E_Int posrou, const E_Int posrov,
                             const E_Int posrow, const E_Int posroe,
                             const E_Float gamma, const E_Float rgp,
                             const E_Float s0,
                             const E_Float betas, const E_Float Cs,
                             std::vector<char*>& vars,
                             FldArrayF& fnew);

/* Calcule les variables composees a partir des grandeurs ro,u,T.
 * Retourne 0 si erreur.*/
  E_Int computeCompVariables2(const FldArrayF& f, const E_Int posro,
                              const E_Int posu, const E_Int posv,
                              const E_Int posw, const E_Int post,
                              const E_Float gamma, const E_Float rgp,
                              const E_Float s0,
                              const E_Float betas, const E_Float Cs,
                              std::vector<char*>& vars,
                              FldArrayF& fnew);

/* Calcule les variables composees a partir des grandeurs conservatives.
 * en association de la fonction computeVariables2
 * Retourne 0 si erreur.*/
  E_Int computeCompVars(const FldArrayF& f,  const E_Int posnew,
                        char* varname,       const E_Int posro,
                        const E_Int posrou,  const E_Int posrov,
                        const E_Int posrow,  const E_Int posroe,
                        const E_Float gamma, const E_Float rgp,
                        const E_Float s0,    const E_Float betas,
                        const E_Float Cs);

/* Calcule les variables composees a partir des grandeurs primitives.
 * en association de la fonction computeVariables2
 * Retourne 0 si erreur.*/
  E_Int computeCompVars2(const FldArrayF& f,  const E_Int posnew,
                         char* varnew,  const E_Int posro,
                         const E_Int posu,    const E_Int posv,
                         const E_Int posw,    const E_Int post,
                         const E_Float gamma, const E_Float rgp,
                         const E_Float s0,    const E_Float betas,
                         const E_Float Cs);

/* Calcule la vitesse absolue */
  void computeVelocity(const FldArrayF& f,
                       const E_Int posro, const E_Int posrou,
                       const E_Int posrov, const E_Int posrow,
                       FldArrayF& velo);

/* Calcule la pression statique */
  void computePressure(const FldArrayF& f,
                       const E_Int posro, const E_Int posrou,
                       const E_Int posrov, const E_Int posrow,
                       const E_Int posroe, const E_Float gamma,
                       FldArrayF& velo, FldArrayF& press);

/* Calcule la temperature statique */
  void computeTemperature(const FldArrayF& f,
                          const E_Int posro, const E_Int posrou,
                          const E_Int posrov, const E_Int posrow,
                          const E_Int posroe, const E_Float gamma,
                          const E_Float rgp,
                          FldArrayF& velo, FldArrayF& press,
                          FldArrayF& temp);

/* calcul de la viscosite du fluide */
  void computeMu(const FldArrayF& f,
                 const E_Int posro, const E_Int posrou,
                 const E_Int posrov, const E_Int posrow,
                 const E_Int posroe,
                 const E_Float gamma, const E_Float rgp,
                 const E_Float betas, const E_Float Cs,
                 FldArrayF& velo, FldArrayF& press,
                 FldArrayF& temp, FldArrayF& mu);

/* Calcule l'entropie s = s0 + rgp*gam/(gam-1)*log(T) - rgp*log(P)
   ou s0 = sref -  rgp*gam/(gam-1)*log(Tref) + rgp*log(Pref)*/
  void computeEntropy(const FldArrayF& f,
                      const E_Int posro, const E_Int posrou,
                      const E_Int posrov, const E_Int posrow,
                      const E_Int posroe, const E_Float gamma,
                      const E_Float rgp, const E_Float s0,
                      FldArrayF& velo,
                      FldArrayF& temp, FldArrayF& press,
                      FldArrayF& s);

/* Calcule l'enthalpie */
  void computeEnthalpy(const FldArrayF& f,
                       const E_Int posro, const E_Int posrou,
                       const E_Int posrov, const E_Int posrow,
                       const E_Int posroe, const E_Float gamma,
                       FldArrayF& velo, FldArrayF& press,
                       FldArrayF& h);

/* Calcul du Mach */
  void computeMach(const FldArrayF& f,
                   const E_Int posro, const E_Int posrou,
                   const E_Int posrov, const E_Int posrow,
                   const E_Int posroe, const E_Float gamma,
                   FldArrayF& velo, FldArrayF& press,
                   FldArrayF& mach);

/* Zip: lecture des options et retourne les valeurs associ�es:
   overlapTol: tolerance geometrique de recouvrement
   matchTol: tolerance geometrique pour les raccords coincidents */
  void readZipperOptions(PyObject* optionList,
                         E_Float& overlapTol,
                         E_Float& matchTol);

// ***************************** GRAD, DIV, CURL **************************** //
/* Creation de la chaine de caracteres pour la  fonction computeGrad
   IN: varString: "x,y,z, var1..." avec var1... variables calculees
   OUT: varStringOut "gradxvar1, gradyvar1, gradzvar1...." */
  void computeGradVarsString(char* varString, char*& varStringOut);
  
/* Creation de la chaine de caracteres pour les fonctions computeDiv/computeDiv2
   IN: varString: "x,y,z, var1..." avec var1... variables calculees
   OUT: varStringOut "divvar1, ...." */
  void computeDivVarsString(char* varString, char*& varStringOut);

  /* gradx, grady, gradz must be allocated previously */
  E_Int computeGradStruct(
    const E_Int ni, const E_Int nj, const E_Int nk,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* field,
    E_Float* gradx, E_Float* grady, E_Float* gradz);

  E_Int computeGradUnstruct(
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    FldArrayI& cn, const char* eltType,
    const E_Float* field,
    E_Float* gradx, E_Float* grady, E_Float* gradz);

  E_Int computeGradNGon(
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* fp, FldArrayI& cn,
    E_Float* gradx, E_Float* grady, E_Float* gradz);

  E_Int compute_gradients_ngon(
    FldArrayI &cn, E_Float *x, E_Float *y, E_Float *z,
    E_Int *owner, E_Int *neigh, E_Float *centers,
    const std::vector<E_Float *> &flds, std::vector<E_Float *> &Gs);
  
// ============================================================================
// Calcul de la divergence d'un champ defini aux noeuds d une grille surfacique
// structuree
// retourne la divergence defini aux centres des cellules
// IN
// ni, nj, nk : dimensions de la grille aux noeuds
// xt, yt, zt : coordonnees des noeuds de la grille
// fieldX, fieldY, fieldZ : champ defini aux noeuds auquel on applique div
// OUT
// div : divergence du champ vectoriel
// ============================================================================
  void compDivStruct2D(
    const E_Int ni, const E_Int nj, const E_Int nk,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* fieldX, const E_Float* fieldY, const E_Float* fieldZ,
    E_Float* div);

  void compDivStruct3D(
    const E_Int ni, const E_Int nj, const E_Int nk,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* fieldX, const E_Float* fieldY, const E_Float* fieldZ,
    E_Float* div);

// ============================================================================
// Calcul de la divergence d'un champ défini aux noeuds d'une grille non structuree
// Retourne la divergence définie aux centres des éléments
// ============================================================================
  void compDivUnstruct2D(
    FldArrayI& cn, const char* eltType,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* fieldX, const E_Float* fieldY, const E_Float* fieldZ,
    E_Float* div);

  void compDivUnstruct3D(
    FldArrayI& cn, const char* eltType,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* fieldX, const E_Float* fieldY, const E_Float* fieldZ,
    E_Float* div);

// ============================================================================
// Calcul du gradient d'un champ defini aux noeuds d une grille structuree
// retourne le gradient defini aux centres des cellules
// IN
// ni: dimensions de la grille aux noeuds
// nbcell: nb de cellules
// xt, yt, zt: coordonnees des noeuds de la grille
// field: champ defini aux noeuds auquel on applique grad
// OUT
// gradx, grady, gradz: gradient de field aux centres des cellules
// ============================================================================
  void compGradStruct1D(
    const E_Int ni,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* field,
    E_Float* gradx, E_Float* grady, E_Float* gradz);
    
  void compGradStruct2D(
    const E_Int ni, const E_Int nj, const E_Int nk,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* field,
    E_Float* gradx, E_Float* grady, E_Float* gradz);
    
  void compGradStruct3D(
    const E_Int ni, const E_Int nj, const E_Int nk,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* field,
    E_Float* gradx, E_Float* grady, E_Float* gradz);

// ============================================================================
// Calcul du rotationnel d'un vecteur defini aux noeuds d une grille surfacique
// Retourne le rotationnel defini aux centres des cellules
//  IN: ni,nj,nk: dimensions du maillage en noeuds
//  IN: xt,yt,zt: coordonnees de la grille
//  IN: u: vecteur dont le rotationnel est a calculer
//  OUT: rotu: rotationnel de u aux centres des cellules
// ============================================================================
  void compCurlStruct2D(
    const E_Int ni, const E_Int nj, const E_Int nk,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* ux, const E_Float* uy, const E_Float* uz,
    E_Float* rotux, E_Float* rotuy, E_Float* rotuz);
  
  void compCurlStruct3D(
    const E_Int ni, const E_Int nj, const E_Int nk,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* ux, const E_Float* uy, const E_Float* uz,
    E_Float* rotx, E_Float* roty, E_Float* rotz);

// ============================================================================
// Calcul du rotationnel moyen d un champ (u,v,w) sur une cellule
// attention cette routine est uniquement 3d
// IN: ind: indice du premier sommet de la cellule
// IN: ni,nj,nk: dimensions du maillage en noeuds
// IN: velo: vecteur dont le rotationnel est a calculer. Defini sur le maillage
// OUT: rotu,rotv,rotw: rotationnel moyen de (u,v,w) sur la cellule
// ============================================================================
  void compMeanCurlOfStructCell(
    const E_Int ind, const E_Int ni, const E_Int nj, const E_Int nk,
    const E_Float* velox, const E_Float* veloy, const E_Float* veloz,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float& rotu, E_Float& rotv, E_Float& rotw);

// ============================================================================
// Calcul du rotationnel moyen d un champ (u,v,w) sur une cellule non structuree
// ============================================================================
  void compMeanCurlOfUnstructCell(
    E_Int noet, FldArrayI& cn, const char* eltType,
    const E_Float* ux, const E_Float* uy, const E_Float* uz,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    E_Float& rotx, E_Float& roty, E_Float& rotz);

  /* Calcul du gradient d'un champ defini aux noeuds d une grille non structuree
   retourne le gradient defini aux centres des elts.
   IN: xt, yt, zt: coordonnees x, y, z des pts de la grille
   IN: cn: connectivite elts-noeuds
   IN: eltType: list of BE element types forming the ME mesh
   IN: field: champ defini aux noeuds auquel on applique grad
   OUT: gradx, grady, gradz: gradient de field %x, %y, %z
  */
  void compGradUnstruct1D(
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    FldArrayI& cn, const char* eltType, const E_Float* field,
    E_Float* gradx, E_Float* grady, E_Float* gradz);

  void compGradUnstruct2D(
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    FldArrayI& cn, const char* eltType, const E_Float* field,
    E_Float* gradx, E_Float* grady, E_Float* gradz);

  void compGradUnstruct3D(
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    FldArrayI& cn, const char* eltType, const E_Float* field,
    E_Float* gradx, E_Float* grady, E_Float* gradz);

  /* Div */
  E_Int computeDivStruct(
    E_Int ni, E_Int nj, E_Int nk,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* fieldX, const E_Float* fieldY, const E_Float* fieldZ,
    E_Float* div);

  E_Int computeDivUnstruct(
    FldArrayI& cn, const char* eltType,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* fieldX, const E_Float* fieldY, const E_Float* fieldZ,
    E_Float* div);

  E_Int computeDivNGon(
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* fpx, const E_Float* fpy, const E_Float* fpz, FldArrayI& cn,
    E_Float* div);

  PyObject* computeGrad2Struct2D(E_Int ni, E_Int nj, E_Int nic, E_Int njc,
                                 const char* varStringOut, E_Float* cellNp,
                                 E_Float* xt, E_Float* yt, E_Float* zt,
                                 FldArrayF& fc, FldArrayF& faceField,
                                 E_Int* cellG, E_Int* cellD,
                                 PyObject* indices, PyObject* field);
  PyObject* computeGrad2Struct3D(E_Int ni, E_Int nj, E_Int nk,
                                 E_Int nic, E_Int njc, E_Int nkc,
                                 const char* varStringOut, E_Float* cellNp,
                                 E_Float* xt, E_Float* yt, E_Float* zt,
                                 FldArrayF& fc, FldArrayF& faceField,
                                 E_Int* cellG, E_Int* cellD,
                                 PyObject* indices, PyObject* field);
  PyObject* computeDiv2Struct2D(E_Int ni, E_Int nj, E_Int nic, E_Int njc,
                                E_Int ixyz, const char* varStringOut, E_Float* cellNp,
                                E_Float* xt, E_Float* yt, E_Float* zt,
                                FldArrayF& fc, FldArrayF& faceField,
                                E_Int* cellG, E_Int* cellD,
                                PyObject* indices, PyObject* fieldX,
                                PyObject* fieldY, PyObject* fieldZ);
  PyObject* computeDiv2Struct3D(E_Int ni, E_Int nj, E_Int nk,
                                E_Int nic, E_Int njc, E_Int nkc,
                                const char* varStringOut, E_Float* cellNp,
                                E_Float* xt, E_Float* yt, E_Float* zt,
                                FldArrayF& fc, FldArrayF& faceField,
                                E_Int* cellG, E_Int* cellD,
                                PyObject* indices, PyObject* fieldX,
                                PyObject* fieldY, PyObject* fieldZ);

  /* Curl */
  E_Int computeCurlStruct(
    const E_Int ni, const E_Int nj, const E_Int nk,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* ux, const E_Float* uy, const E_Float* uz,
    E_Float* rotx, E_Float* roty, E_Float* rotz);

  E_Int computeCurlUnstruct(
    FldArrayI& cn, const char* eltType,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* ux, const E_Float* uy, const E_Float* uz,
    E_Float* rotx, E_Float* roty, E_Float* rotz);

  E_Int computeCurlNGon(
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* fxp, const E_Float* fyp, const E_Float* fzp, FldArrayI& cn,
    E_Float* curlx, E_Float* curly, E_Float* curlz);

// ============================================================================
// Calcul du rotationnel d'un vecteur défini aux noeuds d'une grille 
// non structurée
// retourne le rotationnel défini aux centres des éléments
// ============================================================================
  void compCurlUnstruct2D(
    FldArrayI& cn, const char* eltType,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* ux, const E_Float* uy, const E_Float* uz,
    E_Float* rotx, E_Float* roty, E_Float* rotz);
    
  void compCurlUnstruct3D(
    FldArrayI& cn, const char* eltType,
    const E_Float* xt, const E_Float* yt, const E_Float* zt,
    const E_Float* ux, const E_Float* uy, const E_Float* uz,
    E_Float* rotx, E_Float* roty, E_Float* rotz);


// ******************************* NODE2CENTER ****************************** //
  /* Convertit les noeuds en centres pour les array structures */
  E_Int node2centerStruct(FldArrayF& FNode,
                          E_Int ni, E_Int nj, E_Int nk,
                          E_Int cellN, E_Int mod,
                          FldArrayF& FCenter);
  /* Convertit les noeuds en centres pour les array non-structures */
  E_Int node2centerUnstruct(FldArrayF& FNode,
                            FldArrayI& c,
                            E_Int cellN, E_Int mod,
                            FldArrayF& FCenter);
  /* Calcul des valeurs d un champ aux faces des elements a partir des 
    des valeurs aux noeuds.
    IN: cn: connectivite elts-noeuds
    IN: fieldn: champs aux noeuds
    OUT: fieldf: champs aux centres des faces
  */
  void compUnstrNodes2Faces(FldArrayI& cn, const char* eltType,
                            const E_Float* fieldn,
                            E_Float* fieldf);
  /* Convertit les centres en noeuds pour les array structures */
  E_Int center2nodeStruct(FldArrayF& FCenter,
                          E_Int ni, E_Int nj, E_Int nk,
                          E_Int cellN, E_Int mod,
                          E_Int posx, E_Int posy, E_Int posz,
                          FldArrayF& FNode,
                          E_Int& nin, E_Int& njn, E_Int& nkn);

// ************************************************************************** //
  PyObject* exteriorFacesStructured(char* varString, FldArrayF& f,
                                    E_Int ni, E_Int nj, E_Int nk,
                                    PyObject* Indices);

  PyObject* exteriorFacesBasic(char* varString, FldArrayF& f,
                               FldArrayI& cn, char* eltType,
                               PyObject* Indices);

  PyObject* selectExteriorFacesME(char* varString, FldArrayF& f,    
                                  FldArrayI& cn, char* eltType,
                                  PyObject* indices);                         

  PyObject* selectExteriorFacesNGon3D(char* varString, FldArrayF& f,
                                      FldArrayI& cn,
                                      PyObject* Indices);
  PyObject* selectExteriorFacesNGon2D(char* varString, FldArrayF& f,
                                      FldArrayI& cn,
                                      PyObject* Indices);

  void selectExteriorFacesStruct2D(E_Int ni, E_Int nj, E_Int nk,
                                   FldArrayF& f, char* varString,
                                   PyObject* list);
  void selectExteriorFacesStruct3D(E_Int ni, E_Int nj, E_Int nk,
                                   FldArrayF& f, char* varString,
                                   PyObject* list);

  PyObject* selectExteriorEltsBasic(FldArrayF& f, FldArrayI& cn,
                                    char* eltType, char* varString);
  PyObject* selectExteriorEltsBasic2(FldArrayF& f, FldArrayI& cn,
                                     char* eltType, char* varString,
                                     E_Int posx, E_Int posy, E_Int posz);
  PyObject* selectExteriorEltsBasic3(FldArrayF& f, FldArrayI& cn,
                                     char* eltType, char* varString,
                                     E_Int posx, E_Int posy, E_Int posz);
  PyObject* selectExteriorEltsME(FldArrayF& f, FldArrayI& cn, 
                                 char* eltType, char* varString);
  PyObject* selectExteriorEltsNGon(FldArrayF& f, FldArrayI& cn,
                                   char* varString);

  /* exteriorFacesBasic par topologie.
     IN: nfaces: nb de faces par elt ns, nvertex: nb de sommets
     IN: f: tableau des coordonnees sur la grille ns
     IN: cn: connectivite elt-noeud
     OUT: fext: tableau des noeuds ext,
     OUT: cnext: connectivite arete - noeud pour la frontiere ext
     retourne 0 si erreur interne
  */
  short exteriorFacesBasic(E_Int nfaces, E_Int nvertex,
                           FldArrayF& f, FldArrayI& cn,
                           FldArrayF& fext, FldArrayI& cnext);
  /* idem: autre version topologique */
  short exteriorFacesBasic3(FldArrayF& f, FldArrayI& cn, char*eltType,
                            FldArrayI& cnext,
                            bool boolIndir, PyObject* indicesFaces);
  /* idem: version geometrique */
  short exteriorFacesBasic2(E_Int nfaces, E_Int nvertex,
                            FldArrayF& f, FldArrayI& cn,
                            E_Int posx, E_Int posy, E_Int posz,
                            FldArrayF& fext, FldArrayI& cnext,
                            PyObject* Indices);
  short buildFaceInfo(E_Int et, FldArrayI& cn, FldArrayI& face);
  short buildFaceInfo2(FldArrayI& cn, FldArrayI& face);
  short testCommonFaces(FldArrayI& face1, FldArrayI& face2, FldArrayIS& tag1);
/*
   Corps de la fonction de fusion des elements en fonction de la taille
   des cellules recouvrantes
   IN/OUT: connect: connectivite triangle de la surface a traiter
   IN/OUT: field: coordonnees de la surface a traiter
   IN: posx, posy, posz: position de x,y,z dans field
   IN: indic: indique si l elt est a deraffiner
   IN: eps: tolerance geometrique
*/
  void mergeElements(FldArrayI& connect, FldArrayF& field,
                     E_Int posx, E_Int posy, E_Int posz,
                     E_Float argqual,
                     FldArrayIS& indic, E_Float& eps);

/* Peut etre a mettre ds KCompGeom ?? */
/* Retourne 1 si le triangle noet est degenere*/
  E_Int isDegenerated(E_Int noet, FldArrayI& cn,
                      E_Float* xt, E_Float* yt, E_Float* zt);

/* Determination des sommets externes.
   extNodes vaut 1 si sommet externe, 0 sinon */
  void getExternNodes(FldArrayI& cn, FldArrayF& coord,
                      FldArrayI& extNodes);

/* Determine si les pts ind1 et ind2 sont a fusionner
   IN: et: element a fusionner
   IN: ind1: pt candidat a la fusion (source)
   IN: ind2: pt candidat a la fusion (destination)
   IN: xt, yt, zt: coordonnees des pts
   IN cVE: connectivite Sommet-Elements
   IN: eps: tolerance geometrique de fusion
   IN: indic
   OUT: indo1, indo2: pts trouves pr la fusion
*/
  E_Int testFusion(E_Int et, E_Int ind1, E_Int ind2, FldArrayI& connect,
                   FldArrayIS& indic,
                   E_Float* xt, E_Float* yt, E_Float* zt,
                   std::vector< std::vector<E_Int> >& cVE, E_Float eps,
                   E_Int& indo1, E_Int& indo2);

/* Calcul de teta1+teta2: teta1 = (BC,BA), teta2 = (BA,BD) */
  E_Float computeAngle(E_Int indA, E_Int indB, E_Int indC, E_Int indD,
                       E_Float* xt, E_Float* yt, E_Float* zt);

  /* Retourne 0 si la galette autour de A n est pas convexe */
  E_Int isConvex(E_Int indA, FldArrayI& connect, std::vector<E_Int> & cVE,
                 E_Float* xt, E_Float* yt, E_Float* zt);

  /* Determination du meilleur candidat a la fusion */
  void selectBestCandidatesForMerge(E_Int et0, FldArrayI& selectedVertices,
                                    FldArrayI& connect, E_Float argqual,
                                    E_Float* xt, E_Float* yt, E_Float* zt,
                                    std::vector< std::vector<E_Int> >& cVE,
                                    E_Int& indo1, E_Int& indo2);

  /* Raffine en fonction de indic */
  void refineElements(FldArrayF& f, FldArrayI& cn, FldArrayIS& indic);
  /* Refine partour avec l'interpolation butterfly */
  void refineButterfly(FldArrayF& f, FldArrayI& cn, E_Float w,
    FldArrayF*& fo, FldArrayI*& cno);

  /* Interpole 2 sommets (pour isoLine, isoSurf, isoSurfMC) */
  void vertexInterp(E_Int nfld, E_Float value,
                    FldArrayF& f, E_Int poscellN,
                    E_Float f0, E_Float f1,
                    E_Int ind0, E_Int ind1,
                    FldArrayF& fiso, E_Int& npts);

  /* IsoNode sur un maillage BAR */
  void doIsoNode(FldArrayF& f, FldArrayI& cn, E_Int posf, E_Float value,
                 E_Int poscellN, FldArrayF& fiso, FldArrayI& ciso);
  /* Isoline sur un maillage triangulaire surfacique */
  void doIsoLine(FldArrayF& f, FldArrayI& cn, E_Int posf, E_Float value,
               E_Int poscellN, FldArrayF& fiso, FldArrayI& ciso);
  /* Isosurface dans un maillage tetra (sortie TRI) */
  void doIsoSurf(FldArrayF& f, FldArrayI& cn, E_Int posf, E_Float value,
                 E_Int poscellN, FldArrayF& fiso, FldArrayI& ciso);
  /* Isosurface dans un maillage NGON (sortie TRI) */
  void doIsoSurfNGon(FldArrayF& f, FldArrayI& cn, E_Int posf, E_Float value,
                     E_Int poscellN, FldArrayF& fiso, FldArrayI& ciso);
  /* Isosurface dans un maillage hexa (sortie TRI) */
  void doIsoSurfMC(FldArrayF& f, FldArrayI& cn, E_Int posf, E_Float value,
                   E_Int poscellN, FldArrayF& fiso, FldArrayI& ciso);
  /* Isosurface dans un maillage hexa (sortie QUAD) */
  void doIsoSurfMCQuads(FldArrayF& f, FldArrayI& cn, E_Int posf, E_Float value,
                        E_Int poscellN, FldArrayF& fiso, FldArrayI& ciso);

  /* Construit partiellement une isoSurf pour les cellules taggees par tagC = 0*/
  void createMCQuadsForStructZone(E_Int ni, E_Int nj, E_Int nk,
                                  FldArrayF& f, E_Int poscellN, E_Float* fp,
                                  E_Int* tagC,
                                  E_Int npts0, E_Int nquads0,
                                  E_Int& npts, E_Int& nquad,
                                  FldArrayF& fiso, FldArrayI& ciso);
  void createMCTrisForUnstrZone(FldArrayF& f, FldArrayI& cn, E_Int poscellN,
                                E_Float* fp, E_Int* tagC,
                                E_Int npts0, E_Int ntri0,
                                E_Int& npts, E_Int& ntri,
                                FldArrayF& fiso, FldArrayI& ciso);

  /* Fonctions necessaires a computeDiff */
  void computeDiffForVertex(E_Int npts, FldArrayI& cn, E_Float* field,
                            E_Float* difft);
  void computeDiffForVertexWithCellN(E_Int npts, FldArrayI& cn,
                                     E_Float* field,
                                     E_Float* cellN, E_Float* difft);

  /* Fonctions necessaires a frontFaces */
  PyObject* frontFacesUnstructured(char* varString, FldArrayF& f,
                                   FldArrayI& cn, char* eltType,
                                   FldArrayF& tag);

  short buildFaceInfo(E_Int et, FldArrayI& cn, FldArrayI& face);
}
#undef FldArrayF
#undef FldArrayI
#undef FldArrayIS

#endif
