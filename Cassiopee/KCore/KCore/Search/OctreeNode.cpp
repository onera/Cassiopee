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

# include "Search/OctreeNode.h"

//=============================================================================
/* Cree un noeud d'octree correspondant a une cellule cartesienne de frontiere
   xmin,ymin,zmin,xmin+dh,ymin+dh,zmin+dh
   Les pointeurs sur les 8 eventuels voisins suivants sont stockes */
//=============================================================================
K_SEARCH::OctreeNode::OctreeNode(
  E_Float xmin, E_Float ymin, E_Float zmin, 
  E_Float dh, E_Int level,
  E_Int ind1, E_Int ind2, E_Int ind3, E_Int ind4,
  E_Int ind5, E_Int ind6, E_Int ind7, E_Int ind8):
  _xmin(xmin), _ymin(ymin), _zmin(zmin), _dh(dh), _level(level), 
  _ind1(ind1), _ind2(ind2), _ind3(ind3), _ind4(ind4),
  _ind5(ind5), _ind6(ind6), _ind7(ind7), _ind8(ind8)
{
  _next1 = NULL; _next2 = NULL;
  _next3 = NULL; _next4 = NULL;
  _next5 = NULL; _next6 = NULL;
  _next7 = NULL; _next8 = NULL;
  _voisin1 = NULL; _voisin2 = NULL;
  _voisin3 = NULL; _voisin4 = NULL;
  _voisin5 = NULL; _voisin6 = NULL;
}

//=============================================================================
K_SEARCH::OctreeNode::~OctreeNode()
{
}

//=============================================================================
