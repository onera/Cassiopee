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

# include "CompGeom/compGeom.h"

//=============================================================================
/* Verifie si un quad est convexe ou non
   IN: x1,y1,z1,...: coord du quad dans un ordre consecutif
   OUT: retourne 0: quad convexe
   retourne i: quad non convexe au sommet i
   retourne -1: quad degenere
*/
//=============================================================================
E_Int K_COMPGEOM::checkQuadConvexity(E_Float* pt1, E_Float* pt2, 
                                     E_Float* pt3, E_Float* pt4)
{
  E_Float vx12, vy12, vz12, vx14, vy14, vz14, vx32, vy32, vz32;
  E_Float vx34, vy34, vz34, vx21, vy21, vz21, vx23, vy23, vz23;
  E_Float vx41, vy41, vz41, vx43, vy43, vz43, px1, py1, pz1;
  E_Float px2, py2, pz2, scal;

  E_Float x1 = pt1[0]; E_Float y1 = pt1[1]; E_Float z1 = pt1[2];
  E_Float x2 = pt2[0]; E_Float y2 = pt2[1]; E_Float z2 = pt2[2];
  E_Float x3 = pt3[0]; E_Float y3 = pt3[1]; E_Float z3 = pt3[2];
  E_Float x4 = pt4[0]; E_Float y4 = pt4[1]; E_Float z4 = pt4[2];

  // Premier decoupage
  // Vecteur 1->2
  vx12 = x2-x1; vy12 = y2-y1; vz12 = z2-z1;
  // Vecteur 1->4
  vx14 = x4-x1; vy14 = y4-y1; vz14 = z4-z1;
  // Vecteur 3->2
  vx32 = x2-x3; vy32 = y2-y3; vz32 = z2-z3;
  // Vecteur 3->4
  vx34 = x4-x3; vy34 = y4-y3; vz34 = z4-z3;
  // Produit vectoriel 1
  px1 = vy12*vz14 - vz12*vy14;
  py1 = vz12*vx14 - vx12*vz14;
  pz1 = vx12*vy14 - vy12*vx14;
  // Produit vectoriel 2
  px2 = vy34*vz32 - vz34*vy32;
  py2 = vz34*vx32 - vx34*vz32;
  pz2 = vx34*vy32 - vy34*vx32;
  // Produit scalaire
  scal = px1*px2 + py1*py2 + pz1*pz2;
  if (scal < 0) return 1;
 
  // Deuxieme decoupage
  // Vecteur 2->1
  vx21 = x1-x2; vy21 = y1-y2; vz21 = z1-z2;
  // Vecteur 2->3
  vx23 = x3-x2; vy23 = y3-y2; vz23 = z3-z2;
  // Vecteur 4->1
  vx41 = x1-x4; vy41 = y1-y4; vz41 = z1-z4;
  // Vecteur 4->3
  vx43 = x3-x4; vy43 = y3-y4; vz43 = z3-z4;
  // Produit vectoriel 1
  px1 = vy21*vz23 - vz21*vy23;
  py1 = vz21*vx23 - vx21*vz23;
  pz1 = vx21*vy23 - vy21*vx23;
  // Produit vectoriel 2
  px2 = vy43*vz41 - vz43*vy41;
  py2 = vz43*vx41 - vx43*vz41;
  pz2 = vx43*vy41 - vy43*vx41;
  // Produit scalaire
  scal = px1*px2 + py1*py2 + pz1*pz2;
  if (scal < 0) return 2;
  return 0;
}
