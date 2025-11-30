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
#include "Connect/connect.h"
#include "String/kstring.h"

//=============================================================================
// Get the dimensionality of a BE/ME connectivity based on its element types
// IN: eltType: Type of element
// OUT: dim: dimensionality
//=============================================================================
E_Int K_CONNECT::getDimME(const char* eltType)
{
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);
  E_Int dim = getDimME(eltTypes);
  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
  return dim;
}

E_Int K_CONNECT::getDimME(std::vector<char*> eltTypes)
{
  E_Int dim = 3;
  if (K_STRING::cmp(eltTypes[0], 4, "NODE") == 0) dim = 0;
  else if (K_STRING::cmp(eltTypes[0], 3, "BAR") == 0) dim = 1;
  else if (K_STRING::cmp(eltTypes[0], 3, "TRI") == 0 ||
           K_STRING::cmp(eltTypes[0], 4, "QUAD") == 0) dim = 2;
  return dim;
}
