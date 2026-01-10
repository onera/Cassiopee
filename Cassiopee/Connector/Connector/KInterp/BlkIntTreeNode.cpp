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

#include "BlkInterp.h"

//=============================================================================
K_KINTERP::BlkIntTreeNode::BlkIntTreeNode()
{
   _ind = 0;
   for (E_Int i = 0; i < 6; i++) _BB[i] = 0.;
   _left = NULL;
   _right = NULL;
}

//=============================================================================
K_KINTERP::BlkIntTreeNode::BlkIntTreeNode(
  E_Int ind,
  E_Float xmax, E_Float ymax, E_Float zmax,
  E_Float xmin, E_Float ymin, E_Float zmin)
{
  _ind = ind;
  _BB[0] = xmax;
  _BB[1] = ymax;
  _BB[2] = zmax;
  _BB[3] = xmin;
  _BB[4] = ymin;
  _BB[5] = zmin;
  _left = NULL;
  _right = NULL;
}

//=============================================================================
K_KINTERP::BlkIntTreeNode::~BlkIntTreeNode()
{
}

//=============================================================================



