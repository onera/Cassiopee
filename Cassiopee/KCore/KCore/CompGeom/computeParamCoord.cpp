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

# include "CompGeom/compGeom.h"
# include <stdio.h>

//=============================================================================
// Calcul les coordonnees barycentriques d'un point P dans un triangle
// IN: p1, p2, p3: les 3 pts du triangle
// IN: P: le point en question
// OUT: s, t: coordonnees barycentriques
// Retourne 1 si le pt est dans le triangle, 0 sinon.
//=============================================================================
E_Int K_COMPGEOM::computeParamCoord(E_Float* p1, E_Float* p2, E_Float* p3,
                                    E_Float* P,
                                    E_Float& s, E_Float& t)
{
  E_Float dist2;
  E_Bool in;
  E_Float xp, yp, zp;
  distanceToTriangle(p1, p2, p3, P, 1, 
                     dist2, in, 
                     xp, yp, zp,
                     s, t);
  if (dist2 > 1.e-12)
    printf("Warning: computeParamCoord: point P is not in triangle plane.\n");
  if (in == true) return 1;
  else return 0;
}
