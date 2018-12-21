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

#include "Slot1D.h"

//=============================================================================
Slot1D::Slot1D(int no, int gridPosI, int gridPosJ, double bgBlend)
{
  _no = no;
  _gridPosI = gridPosI;
  _gridPosJ = gridPosJ;
  _bgBlend = bgBlend;
}

//=============================================================================
Slot1D::~Slot1D()
{
  unsigned int nz = _zones.size();
  for (unsigned int i = 0; i < nz; i++) delete _zones[i];
}
