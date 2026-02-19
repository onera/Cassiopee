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

#ifndef _POST_CUTPLANE_PLANE_INTERSECTION_H_
#define _POST_CUTPLANE_PLANE_INTERSECTION_H_

# include "cutPlane.h"
# include <vector>

#define FldArrayF K_FLD::FldArrayF
#define FldArrayI K_FLD::FldArrayI

namespace K_POST
{
/* Given the points defined by coord, for each segment [i,i+1] or [j,j+1] ...
   search for the intersection with the plane. If the point is found, store it
   in plane Field.
   In: coef of the plane equation coefaX+ coefbY + coefcZ + coefd=0 
       ni, nj,nk size of the mesh in each direction
       posx : position of x coordinate in coord
       posy : position of y coordinate in coord
       posz : position of z coordinate in coord
       posc : position of celln in coord. posc = 0 if no celln found
       field: coordinates + solution  
       tagC : tag = 1 si la cellule est a intersecter
   Out : storage of the fields of intersection points of mesh with plane
   Retourne 0 si erreur interne
*/ 
short computeStructIntersectionWithPlane( 
  K_INTERP::InterpData* interpData,
  K_INTERP::InterpData::InterpolationType interpType,
  E_Float coefa, E_Float coefb, E_Float coefc,
  E_Float coefd, E_Int ni, E_Int nj, E_Int nk,
  E_Int posx, E_Int posy, E_Int posz, E_Int posc,
  FldArrayF& coord,   FldArrayI& tagC,
  FldArrayF& intersectPts,
  FldArrayF& volOfIntersectPts);

/* Calcul de l intersection des grilles non structurees TETRA avec le plan
   Retourne la liste des pts d intersection et le volume de la cellule 
   d interpolation correspondante 
   tagC = 1 : on recherche l intersection pour la cellule */
void computeUnstrIntersectionWithPlane(
  E_Float coefa, E_Float coefb, E_Float coefc, E_Float coefd, 
  K_INTERP::InterpData* interpData,
  K_INTERP::InterpData::InterpolationType interpType,
  FldArrayI& connect,
  E_Int posx, E_Int posy, E_Int posz, E_Int posc, FldArrayF& field, 
  FldArrayI& tagC,
  FldArrayF& intersectPts, FldArrayF& volOfIntersectPts);

/* Given 1- 2 vertices A(coord(indA,.)) and Bi(coord(indp[i],.)) 
         2- the plane equation coefficients coefa, coefb...
         3- the cellNatureField (if it exists : cellN != -1)
   Compute the intersection pt H coordinates of (AB) and the plane 
   If (AB) is in the plane insert A and B, 
   else insert H if and only if k in [0,1], where k is such that  AH = k.AB 
   cnt is the number for the intersectPts array of the last point inserted*/
void searchStructIntersectForSegment(
  K_INTERP::InterpData* interpData,
  K_INTERP::InterpData::InterpolationType interpType,
  E_Float coefa, E_Float coefb, E_Float coefc, E_Float coefd,
  E_Int ni, E_Int nj, E_Int nk, E_Int indA, 
  E_Int posx, E_Int posy, E_Int posz, E_Int poscelln, 
  E_Int& cnt, E_Int indp, FldArrayF& field,
  FldArrayF& intersectPts, FldArrayF& volOfIntersectPts);

/* Calcul de l intersection entre une arete [indA, indB] du tetra no 'et'
   et le plan/
   IN : indA, indB : indices ds field des numeros des sommets de l arete
   IN : field : champs defini aux sommets de la grille tetra
   IN : connect : connectivite elt->sommet de la grille
   IN : cellVol : volume de la cellule 'et'
   IN : posx, posy, posz, posc : positions de x,y,z,celln ds field
   IN : interpData : interpData tetra de la grille
   IN/OUT :  cnt : compteur ds intersectPts
   IN:OUT : intersectPts : tableau des pts d intersections trouves
   IN/OUT : volOfIntersectPts : volume de la cellule d interpolation */
void searchUnstrIntersectForSegment(
  E_Float coefa, E_Float coefb, E_Float coefc, E_Float coefd,
  E_Int indA, E_Int indB, E_Int posx, E_Int posy, E_Int posz, E_Int posc,
  E_Float cellVol, FldArrayI& connect, FldArrayF& field, 
  K_INTERP::InterpData* interpData, 
  K_INTERP::InterpData::InterpolationType interpType,
  E_Int& cnt, FldArrayF& intersectPts, FldArrayF& volOfIntersectPts);
} //fin namespace
#undef FldArrayF
#undef FldArrayI
#endif
