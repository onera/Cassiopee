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

// ============================================================================
// Reorder the vertex indices of a ME connectivity such that each facet normal
// is pointing outward.
// The input connectivity must be such that the image points of the 'base' facet
// are correct.
// Return 
//     0 if no vertex reordering took place (incl. when the coordinates are
//                                           missing in varString).
//     1 if vertex reordering took place.
// ============================================================================
E_Int K_CONNECT::reorderUnstruct(const char* varString, FldArrayF& f, 
                                 FldArrayI& cn, const char* eltType)
{
  // Get dimensionality
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);
  E_Int dim = K_CONNECT::getDimME(eltTypes);
  if(dim < 2) return 0;

  // Check that coordinates are present
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_WarnEx(PyExc_Warning,
                 "reorderUnstruct: coords must be present in array.", 1);
    return 0;
  }
  posx++; posy++; posz++;

  E_Float* xp = f.begin(posx);
  E_Float* yp = f.begin(posy);
  E_Float* zp = f.begin(posz);

  // List of vertex indices forming each of the element's facets, for each
  // connectivity
  std::vector<std::vector<std::vector<E_Int> > > facetspc(nc);

  // List the pairs of vertex indices to swap for each connectivity if a face
  // normal is found to be pointing inward (i.e., in the element's volume)
  std::vector<std::vector<std::pair<E_Int, E_Int> > > ind2Swap(nc);
  for (E_Int ic = 0; ic < nc; ic++)
  {
    if (K_STRING::cmp(eltTypes[ic], 3, "TRI") == 0)
      ind2Swap[ic].push_back({2, 3});
    else if (K_STRING::cmp(eltTypes[ic], 4, "QUAD") == 0)
      ind2Swap[ic].push_back({2, 4});
    else if (K_STRING::cmp(eltTypes[ic], 5, "TETRA") == 0)
      ind2Swap[ic].push_back({2, 3});
    else if (K_STRING::cmp(eltTypes[ic], 4, "PYRA") == 0)
      ind2Swap[ic].push_back({2, 4});
    else if (K_STRING::cmp(eltTypes[ic], 5, "PENTA") == 0)
    {
      ind2Swap[ic].push_back({2, 3});
      ind2Swap[ic].push_back({5, 6});
    }
    else if (K_STRING::cmp(eltTypes[ic], 4, "HEXA") == 0)
    {
      ind2Swap[ic].push_back({2, 4});
      ind2Swap[ic].push_back({6, 8});
    }

    // Fill the list of facets for this BE
    E_Int ierr = K_CONNECT::getEVFacets(facetspc[ic], eltTypes[ic]);
    if (ierr == 1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "reorderUnstruct: element type not supported.");
      for (E_Int ic = 0; ic < nc; ic++) delete [] eltTypes[ic];
      return -1;
    }
  }
  
  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];

  E_Int reoriented = 0;

  #pragma omp parallel reduction(+: reoriented)
  {
    E_Int nelts, nvpe;
    E_Int ind, ind1, ind2, ind3;
    E_Float xb, yb, zb;
    E_Float x1, y1, z1, x2, y2, z2, x3, y3, z3;
    E_Float x12, y12, z12, x13, y13, z13;
    E_Float xf, yf, zf, nx, ny, nz, xbf, ybf, zbf;
    E_Float res;
    
    // Loop over each element of all connectivities
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn.getConnect(ic));
      nelts = cm.getSize();
      nvpe = cm.getNfld();
      const std::vector<E_Int>& facet0 = facetspc[ic][0];

      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        // Compute the position of the cell barycenter (xb, yb, zb)
        xb = 0.; yb = 0.; zb = 0.;
        for (E_Int j = 1; j <= nvpe; j++)
        {
          ind = cm(i, j)-1;
          xb += xp[ind]; yb += yp[ind]; zb += zp[ind];
        }
        xb /= nvpe; yb /= nvpe; zb /= nvpe;

        // Compute the facet normal (nx, ny, nz) using the first 3 points of the
        // first facet, facet0
        ind1 = cm(i, facet0[0]) - 1;
        ind2 = cm(i, facet0[1]) - 1;
        ind3 = cm(i, facet0[2]) - 1;
        x1 = xp[ind1]; y1 = yp[ind1]; z1 = zp[ind1];
        x2 = xp[ind2]; y2 = yp[ind2]; z2 = zp[ind2];
        x3 = xp[ind3]; y3 = yp[ind3]; z3 = zp[ind3];
        x12 = x2 - x1; x13 = x3 - x1;
        y12 = y2 - y1; y13 = y3 - y1;
        z12 = z2 - z1; z13 = z3 - z1;
        K_MATH::cross(x12, y12, z12, x13, y13, z13, nx, ny, nz);

        // Compute the face center and the vector from element barycenter to the
        // face center (xbf, ybf, zbf)
        xf = K_MATH::ONE_THIRD * (x1 + x2 + x3);
        yf = K_MATH::ONE_THIRD * (y1 + y2 + y3);
        zf = K_MATH::ONE_THIRD * (z1 + z2 + z3);
        xbf = xf - xb; ybf = yf - yb; zbf = zf - zb;

        // If dot(n, bf) is positive, the facet normal is pointing outward
        // which is the correct orientation
        K_MATH::dot(nx, ny, nz, xbf, ybf, zbf, res);

        if (res < 0.)
        {
          reoriented = true;
          for (size_t j = 0; j < ind2Swap[ic].size(); j++)
          {
            std::swap(
              cm(i, ind2Swap[ic][j].first),
              cm(i, ind2Swap[ic][j].second)
            );
          }
        }
      }
    }
  }

  if (reoriented > 0) reoriented = 1;
  return reoriented;
}
