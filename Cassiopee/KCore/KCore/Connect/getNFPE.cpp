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
// Get the number of facets per element type of a Multiple Element connectivity
// Note that 'faces' of 1D and 2D elements are vertices and edges, respectively,
// if expandToLowerDim is set to True
// IN: eltType: Type of basic element
// OUT: nfpe: number of facets per element type
// Return: error index, 0 (ok), 1 (error)
//=============================================================================
E_Int K_CONNECT::getNFPE(
  std::vector<E_Int>& nfpe, const char* eltType,
  E_Bool expandToLowerDim
)
{
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);
  E_Int nc = eltTypes.size();
  nfpe.resize(nc);

  for (E_Int ic = 0; ic < nc; ic++)
  {
    char* eltTypeConn = eltTypes[ic];

    if (K_STRING::cmp(eltTypeConn, "BAR") == 0 || 
        K_STRING::cmp(eltTypeConn, "BAR*") == 0)
    {
      if (expandToLowerDim) nfpe[ic] = 2;
      else nfpe[ic] = 0;
    }
    else if (K_STRING::cmp(eltTypeConn, "TRI") == 0 ||
             K_STRING::cmp(eltTypeConn, "TRI*") == 0)
    {
      if (expandToLowerDim) nfpe[ic] = 3;
      else nfpe[ic] = 1;
    }
    else if (K_STRING::cmp(eltTypeConn, "QUAD") == 0 || 
             K_STRING::cmp(eltTypeConn, "QUAD*") == 0)
    {
      if (expandToLowerDim) nfpe[ic] = 4;
      else nfpe[ic] = 1;
    }
    else if (K_STRING::cmp(eltTypeConn, "TETRA") == 0 || 
             K_STRING::cmp(eltTypeConn, "TETRA*") == 0)
    {
      nfpe[ic] = 4;
    }
    else if (K_STRING::cmp(eltTypeConn, "PYRA") == 0 || 
             K_STRING::cmp(eltTypeConn, "PYRA*") == 0)
    {
      nfpe[ic] = 5;
    }
    else if (K_STRING::cmp(eltTypeConn, "PENTA") == 0 || 
             K_STRING::cmp(eltTypeConn, "PENTA*") == 0)
    {
      nfpe[ic] = 5;
    }
    else if (K_STRING::cmp(eltTypeConn, "HEXA") == 0 || 
             K_STRING::cmp(eltTypeConn, "HEXA*") == 0) 
    {
      nfpe[ic] = 6;
    }
    else
    {
      printf("Error in getNFPE: element type unknown or not covered.");
      for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
      return 1;
    };
  }
  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
  return 0;
}
