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

#include "Connect/connect.h"
#include "String/kstring.h"
#include <vector>

//=============================================================================
/* Get all facets of a basic element
  IN: eltType: Type of basic element
  IN: allow_degenerated: whether to allow degenerated elements (in which case
                         the last vertex of a degenerated facet is a repetition
                         of the first)
  OUT: facets: vertex indices for each of the facets - indexing starts at 1
  Return: error index 0 (ok), 1 (error) */
//=============================================================================
E_Int K_CONNECT::getEVFacets(std::vector<std::vector<E_Int> >& facets,
                             const char* eltType, E_Bool allow_degenerated)
{
  E_Int ierr = 0;
  facets.clear();

  if (K_STRING::cmp(eltType, "BAR") == 0 || K_STRING::cmp(eltType, "BAR*") == 0)
  {
    facets.push_back({1}); facets.push_back({2});
  }
  else if (K_STRING::cmp(eltType, "TRI") == 0 || K_STRING::cmp(eltType, "TRI*") == 0)
  {
    facets.push_back({1, 2}); facets.push_back({2, 3});
    facets.push_back({3, 1});
  }
  else if (K_STRING::cmp(eltType, "QUAD") == 0 || K_STRING::cmp(eltType, "QUAD*") == 0)
  {
    facets.push_back({1, 2}); facets.push_back({2, 3});
    facets.push_back({3, 4}); facets.push_back({4, 1});
  }
  else if (K_STRING::cmp(eltType, "TETRA") == 0 || K_STRING::cmp(eltType, "TETRA*") == 0)
  {
    facets.push_back({1, 3, 2}); facets.push_back({1, 2, 4});
    facets.push_back({2, 3, 4}); facets.push_back({3, 1, 4});
  }
  else if (K_STRING::cmp(eltType, "PYRA") == 0 || K_STRING::cmp(eltType, "PYRA*") == 0)
  {
    facets.push_back({1, 4, 3, 2});
    if (allow_degenerated)
    {
      facets.push_back({1, 2, 5, 1}); facets.push_back({2, 3, 5, 2});
      facets.push_back({3, 4, 5, 3}); facets.push_back({4, 1, 5, 4});
    }
    else
    {
      facets.push_back({1, 2, 5}); facets.push_back({2, 3, 5});
      facets.push_back({3, 4, 5}); facets.push_back({4, 1, 5});
    }
  }
  else if (K_STRING::cmp(eltType, "PENTA") == 0 || K_STRING::cmp(eltType, "PENTA*") == 0)
  {
    facets.push_back({1, 2, 5, 4}); facets.push_back({2, 3, 6, 5});
    facets.push_back({3, 1, 4, 6});
    if (allow_degenerated)
    {
      facets.push_back({1, 3, 2, 1}); facets.push_back({4, 5, 6, 4});
    }
    else
    {
      facets.push_back({1, 3, 2}); facets.push_back({4, 5, 6});
    }
  }
  else if (K_STRING::cmp(eltType, "HEXA") == 0 || K_STRING::cmp(eltType, "HEXA*") == 0)
  {
    facets.push_back({1, 4, 3, 2}); facets.push_back({1, 2, 6, 5});
    facets.push_back({2, 3, 7, 6}); facets.push_back({3, 4, 8, 7});
    facets.push_back({1, 5, 8, 4}); facets.push_back({5, 6, 7, 8});
  }
  else ierr = 1;
  return ierr;
}
