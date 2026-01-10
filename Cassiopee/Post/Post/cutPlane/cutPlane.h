/*    
    Copyright 2013-2026 ONERA.

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
#ifndef _POST_CUTPLANE_CUTPLANE_H_
#define _POST_CUTPLANE_CUTPLANE_H_

# include "../post.h"
# include <vector>
# include "PlaneIntersection.h"

#define FldArrayF K_FLD::FldArrayF
#define FldArrayI K_FLD::FldArrayI

namespace K_POST
{
/* Dans toutes les routines ci dessous: 
   IN: structInterpDatas: liste des interpDatas structurees
   IN: nis, njs, nks: liste des dimensions de chq array structure
   IN: posxs, posys, poszs, poscs: position de x,y,z,cellN dans structFields
   IN: structFields: champs des grilles structurees (x,y,z, ev cellN, ...)
   IN: unstrInterpDatas: liste des interpDatas non structurees
   IN: connectu: liste des connectivites TETRA
   IN: posxu, posyu, poszu: position de x,y,z,cellN ds unstrFields
   IN: unstrFields: champs des grilles TETRA (x,y,z,(celln),...)
*/

/* Calcul des intersections des grilles structurees avec le plan:
   IN: coefa, coefb, coefc, coefd: coefficients de l eq. de plan
   IN: tagS: tag des cellules structurees a intersecter. tag=1: a intersecter 
   IN: tagU: tag des cellules non structurees a intersecter. tag=1: a intersecter 
   OUT: vectOfIntersectPts: liste des pts d intersection pour chq array
   OUT: volOfIntersectPts: volume de la cellule d'interpolation correspondante 
   pour chq pt d intersection. Par array.*/
void compIntersectionWithPlane(
  E_Float coefa, E_Float coefb, E_Float coefc, E_Float coefd,
  std::vector<K_INTERP::InterpData*>& structInterpDatas,
  std::vector<E_Int>& nis, std::vector<E_Int>& njs, std::vector<E_Int>& nks,
  std::vector<E_Int>& posxs, std::vector<E_Int>& posys, 
  std::vector<E_Int>& poszs, std::vector<E_Int>& poscs,
  std::vector<FldArrayF*>& structFields, std::vector<FldArrayI*>& tagS, 
  std::vector<K_INTERP::InterpData*>& unstrInterpDatas,
  std::vector<FldArrayI*>& connectu,
  std::vector<E_Int>& posxu, std::vector<E_Int>& posyu, 
  std::vector<E_Int>& poszu, std::vector<E_Int>& poscu,
  std::vector<FldArrayF*>& unstrFields, std::vector<FldArrayI*>& tagU, 
  std::vector<FldArrayF*>& vectOfIntersectPts,
  std::vector<FldArrayF*>& vectOfInterpCellVol,
  K_INTERP::InterpData::InterpolationType interpType);

/* Select points from overlapping zones: 
   selected one is the one whose cell volume  is the smallest one 
   OUT : selectedPts : tableau des pts selectionnes.
*/
void selectPointsInOverlappingZones( 
  std::vector<E_Int>& nis, std::vector<E_Int>& njs,std::vector<E_Int>& nks,
  std::vector<E_Int>& posxs, std::vector<E_Int>& posys, 
  std::vector<E_Int>& poszs, std::vector<E_Int>& poscs,
  std::vector<K_INTERP::InterpData*>& structInterpDatas,
  std::vector<FldArrayF*>& structFields,
  std::vector<FldArrayI*>& connectu,
  std::vector<E_Int>& posxu, std::vector<E_Int>& posyu, 
  std::vector<E_Int>& poszu, std::vector<E_Int>& poscu,
  std::vector<K_INTERP::InterpData*>& unstrInterpDatas,
  std::vector<FldArrayF*>& unstrFields,
  std::vector<FldArrayF*>& vectOfIntersectPts,
  std::vector<FldArrayF*>& volOfIntersectPts,
  FldArrayF& selectedPts,
  K_INTERP::InterpData::InterpolationType interpType);

/* Compute the triangulation of a given set of points.
   IN : coefa, coefb, coefc, coefd : coefficients defining the plane the
        points belong to. Plane equation : 
        coefa * x + coefb * y + coefc * z + coefd = 0.
   IN : structInterpDatas : liste des interpDatas de tous les blocs
   IN/OUT : field : selectedpoints for triangulation
   OUT : connect : connectivity of the set of points by triangles.
*/
void makeTriangulation(
  E_Float coefa, E_Float coefb,E_Float coefc, E_Float coefd,
  std::vector<E_Int>& nis, std::vector<E_Int>& njs, 
  std::vector<E_Int>& nks, std::vector<E_Int>& posxs,
  std::vector<E_Int>& posys, std::vector<E_Int>& poszs,
  std::vector<E_Int>& poscs, std::vector<FldArrayF*>& structF,
  std::vector<K_INTERP::InterpData*>& structInterpDatas,
  std::vector<FldArrayI*>& connectu, 
  std::vector<E_Int>& posxu, std::vector<E_Int>& posyu, 
  std::vector<E_Int>& poszu, std::vector<E_Int>& poscu, 
  std::vector<FldArrayF*>& unstrF,
  std::vector<K_INTERP::InterpData*>& unstrInterpDatas,
  FldArrayF& field, FldArrayI& connect,
  K_INTERP::InterpData::InterpolationType interpType);


/* Test if the barycenter of the triangle is interpolated or not 
   if no interpolation point with cellN > 0 is found the triangle is
   blanked*/
void removeTrianglesWithBlankedCenters(
  std::vector<E_Int>& nis, std::vector<E_Int>& njs, 
  std::vector<E_Int>& nks, std::vector<E_Int>& posxs,
  std::vector<E_Int>& posys, std::vector<E_Int>& poszs,
  std::vector<E_Int>& poscs, std::vector<FldArrayF*>& structF,
  std::vector<K_INTERP::InterpData*>& structInterpDatas,
  std::vector<FldArrayI*>& connectu, 
  std::vector<E_Int>& posxu, std::vector<E_Int>& posyu, 
  std::vector<E_Int>& poszu, std::vector<E_Int>& poscu, 
  std::vector<FldArrayF*>& unstrF,
  std::vector<K_INTERP::InterpData*>& unstrInterpDatas,
  FldArrayF& field, FldArrayI& connect,
  K_INTERP::InterpData::InterpolationType interpType);

// /* Remove triangles with one vertex that is blanked (masked)  */
// void removeTrianglesWithBlankedVertices(FldArrayF& field, 
//                                         std::list<Triangle*>&triangles);
}//fin namespace

#undef FldArrayF
#undef FldArrayI
#endif
