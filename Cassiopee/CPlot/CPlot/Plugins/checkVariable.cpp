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
#include "../Data.h"

//=============================================================================
/*
  check if the variable varName is a variable of the zone.
  Return the no of this variable in the list.
  If variable is not found return -1.
*/
//=============================================================================
E_Int Data::checkVariable(E_Int zone, const char* varName)
{
  E_Int i;
  Zone* z = _zones[zone];
  for (i = 0; i < z->nfield; i++)
  {
    if (strcmp(z->varnames[i], varName) == 0) return i;
  }
  return -1;
}
