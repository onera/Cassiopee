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

# include "post.h"

// ============================================================================
// Compute surface integral of the moment M (OM^F.vect(n)), coordinates 
// and F have the same size
// ============================================================================
void K_POST::integMomentNormUnstruct(
  FldArrayI& cn, const char* eltType,
  const E_Float cx, const E_Float cy, const E_Float cz, 
  const E_Float* ratio, const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* sx, const E_Float* sy, const E_Float* sz, 
  const E_Float* field, E_Float* result
)
{
  E_Float resx, resy, resz;
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  resx = 0.0;
  resy = 0.0;
  resz = 0.0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    E_Int ind1, ind2, ind3;
    E_Float f, f1, f2, f3;
    E_Float mx, my, mz, sx0, sy0, sz0;
    E_Float centerx, centery, centerz;
  
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    
    if (strcmp(eltTypes[ic], "TRI") == 0)
    {
      for (E_Int i = 0; i < nelts; i++)
      {
        ind1 = cm(i, 1) - 1;
        ind2 = cm(i, 2) - 1;
        ind3 = cm(i, 3) - 1;

        f1 = ratio[ind1] * field[ind1];
        f2 = ratio[ind2] * field[ind2];
        f3 = ratio[ind3] * field[ind3];

        sx0 = sx[i];
        sy0 = sy[i];
        sz0 = sz[i];

        f = K_CONST::ONE_THIRD * (f1 + f2 + f3);
        centerx = K_CONST::ONE_THIRD * (xt[ind1] + xt[ind2] + xt[ind3]) - cx;
        centery = K_CONST::ONE_THIRD * (yt[ind1] + yt[ind2] + yt[ind3]) - cy;
        centerz = K_CONST::ONE_THIRD * (zt[ind1] + zt[ind2] + zt[ind3]) - cz;

        mx = centery * sz0 - centerz * sy0;
        my = centerz * sx0 - centerx * sz0;
        mz = centerx * sy0 - centery * sx0;

        resx += f * mx;
        resy += f * my;
        resz += f * mz;
      }
    }
    else if (strcmp(eltTypes[ic], "QUAD") == 0)
    {
      E_Int ind4;
      E_Float f4;

      for (E_Int i = 0; i < nelts; i++)
      {
        ind1 = cm(i, 1) - 1;
        ind2 = cm(i, 2) - 1;
        ind3 = cm(i, 3) - 1;
        ind4 = cm(i, 4) - 1;

        f1 = ratio[ind1] * field[ind1];
        f2 = ratio[ind2] * field[ind2];
        f3 = ratio[ind3] * field[ind3];
        f4 = ratio[ind4] * field[ind4];

        sx0 = sx[i];
        sy0 = sy[i];
        sz0 = sz[i];

        f = K_CONST::ONE_FOURTH * (f1 + f2 + f3 + f4);
        centerx = K_CONST::ONE_FOURTH * (xt[ind1] + xt[ind2] + xt[ind3] + xt[ind4]) - cx;
        centery = K_CONST::ONE_FOURTH * (yt[ind1] + yt[ind2] + yt[ind3] + yt[ind4]) - cy;
        centerz = K_CONST::ONE_FOURTH * (zt[ind1] + zt[ind2] + zt[ind3] + zt[ind4]) - cz;

        mx = centery * sz0 - centerz * sy0;
        my = centerz * sx0 - centerx * sz0;
        mz = centerx * sy0 - centery * sx0;

        resx += f * mx;
        resy += f * my;
        resz += f * mz;
      }
    }
    else
    {
      fprintf(stderr, "Error: in K_POST::integMomentNormUnstruct.\n");
      fprintf(stderr, "Unsupported type of element, %s.\n", eltTypes[ic]);
    }
  }

  result[0] = resx;
  result[1] = resy;
  result[2] = resz;

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

// ============================================================================
// Compute linear integral of the moment.norm (OM^F.n), coordinates 
// are defined in nodes and F is defined in center, unstructured case
// ============================================================================
void K_POST::integMomentNormUnstructNodeCenter(
  FldArrayI& cn, const char* eltType,
  const E_Float cx, const E_Float cy, const E_Float cz,
  const E_Float* ratio, const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* sx, const E_Float* sy, const E_Float* sz, 
  const E_Float* field, E_Float* result
)
{
  E_Float resx, resy, resz;
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  resx = 0.0;
  resy = 0.0;
  resz = 0.0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    E_Int ind1, ind2, ind3;
    E_Float f;
    E_Float mx, my, mz, sx0, sy0, sz0;
    E_Float centerx, centery, centerz;
  
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    
    if (strcmp(eltTypes[ic], "TRI") == 0)
    {
      for (E_Int i = 0; i < nelts; i++)
      {
        f = ratio[i] * field[i];

        ind1 = cm(i, 1) - 1;
        ind2 = cm(i, 2) - 1;
        ind3 = cm(i, 3) - 1;

        sx0 = sx[i];
        sy0 = sy[i];
        sz0 = sz[i];

        centerx = xt[ind1] + xt[ind2] + xt[ind3];
        centerx = K_CONST::ONE_THIRD * centerx - cx;

        centery = yt[ind1] + yt[ind2] + yt[ind3];
        centery = K_CONST::ONE_THIRD * centery - cy;

        centerz = zt[ind1] + zt[ind2] + zt[ind3];
        centerz = K_CONST::ONE_THIRD * centerz - cz;

        mx = centery * sz0 - centerz * sy0;
        my = centerz * sx0 - centerx * sz0;
        mz = centerx * sy0 - centery * sx0;

        resx += f * mx;
        resy += f * my;
        resz += f * mz;
      }
    }
    else if (strcmp(eltTypes[ic], "QUAD") == 0)
    {
      E_Int ind4;

      for (E_Int i = 0; i < nelts; i++)
      {
        f = ratio[i] * field[i];

        ind1 = cm(i, 1) - 1;
        ind2 = cm(i, 2) - 1;
        ind3 = cm(i, 3) - 1;
        ind4 = cm(i, 4) - 1;

        sx0 = sx[i];
        sy0 = sy[i];
        sz0 = sz[i];

        centerx = xt[ind1] + xt[ind2] + xt[ind3] + xt[ind4];
        centerx = K_CONST::ONE_FOURTH * centerx - cx;

        centery = yt[ind1] + yt[ind2] + yt[ind3] + yt[ind4];
        centery = K_CONST::ONE_FOURTH * centery - cy;

        centerz = zt[ind1] + zt[ind2] + zt[ind3] + zt[ind4];
        centerz = K_CONST::ONE_FOURTH * centerz - cz;

        mx = centery * sz0 - centerz * sy0;
        my = centerz * sx0 - centerx * sz0;
        mz = centerx * sy0 - centery * sx0;

        resx += f * mx;
        resy += f * my;
        resz += f * mz;
      }
    }
    else
    {
      fprintf(stderr, "Error: in K_POST::integMomentNormUnstructNodeCenter.\n");
      fprintf(stderr, "Unsupported type of element, %s.\n", eltTypes[ic]);
    }
  }

  result[0] = resx;
  result[1] = resy;
  result[2] = resz;

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}
