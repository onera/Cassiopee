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
// Compute surface integral of the moment M (OM^F), coordinates 
// and field have the same size
// I(ABCD) = Aire(ABCD)*(F(A)+F(B)+F(C)+F(D))/4
// Aire(ABCD) = ||AB^AC||/2 + ||DB^DC||/2
// ============================================================================
void K_POST::integMomentUnstruct(
  FldArrayI& cn, const char* eltType,
  const E_Float cx, const E_Float cy, const E_Float cz, const E_Float* ratio,
  const E_Float* xt, const E_Float* yt, const E_Float* zt, const E_Float* surf,
  const E_Float* vx, const E_Float* vy, const E_Float* vz, E_Float* result
)
{
  E_Float res1, res2, res3;
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  res1 = 0.0;
  res2 = 0.0;
  res3 = 0.0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    
    E_Int ind1, ind2, ind3;
    E_Float fx, fy, fz;
    E_Float m1x, m2x, m3x;
    E_Float m1y, m2y, m3y;
    E_Float m1z, m2z, m3z;
    E_Float dx, dy, dz, r1, r2, r3, si;

    if (strcmp(eltTypes[ic], "TRI") == 0)
    {
      for (E_Int i = 0; i < nelts; i++)
      {
        ind1 = cm(i, 1) - 1;
        ind2 = cm(i, 2) - 1;
        ind3 = cm(i, 3) - 1;

        r1 = ratio[ind1];
        dx = xt[ind1] - cx;
        dy = yt[ind1] - cy;
        dz = zt[ind1] - cz;
        fx = vx[ind1];
        fy = vy[ind1];
        fz = vz[ind1];

        m1x = dy * fz - dz * fy;
        m1y = dz * fx - dx * fz;
        m1z = dx * fy - dy * fx;

        r2 = ratio[ind2];
        dx = xt[ind2] - cx;
        dy = yt[ind2] - cy;
        dz = zt[ind2] - cz;
        fx = vx[ind2];
        fy = vy[ind2];
        fz = vz[ind2];

        m2x = dy * fz - dz * fy;
        m2y = dz * fx - dx * fz;
        m2z = dx * fy - dy * fx;

        r3 = ratio[ind3];
        dx = xt[ind3] - cx;
        dy = yt[ind3] - cy;
        dz = zt[ind3] - cz;
        fx = vx[ind3];
        fy = vy[ind3];
        fz = vz[ind3];

        m3x = dy * fz - dz * fy;
        m3y = dz * fx - dx * fz;
        m3z = dx * fy - dy * fx;

        si = surf[i];
        res1 += si * (r1 * m1x + r2 * m2x + r3 * m3x);
        res2 += si * (r1 * m1y + r2 * m2y + r3 * m3y);
        res3 += si * (r1 * m1z + r2 * m2z + r3 * m3z);
      }
    }
    else
    {
      fprintf(stderr, "Error: in K_POST::integMomentUnstruct.\n");
      fprintf(stderr, "Unsupported type of element, %s.\n", eltTypes[ic]);
    }
  }

  result[0] = res1 * K_CONST::ONE_THIRD; // TODO
  result[1] = res2 * K_CONST::ONE_THIRD;
  result[2] = res3 * K_CONST::ONE_THIRD;

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

// ============================================================================
// Compute linear integral of the moment M (OM^ F), coordinates
// and field have the same size
// ============================================================================
void K_POST::integMomentUnstruct1D(
  FldArrayI& cn, const char* eltType,
  const E_Float cx, const E_Float cy, const E_Float cz, const E_Float* ratio,
  const E_Float* xt, const E_Float* yt, const E_Float* zt, const E_Float* length,
  const E_Float* vx, const E_Float* vy, const E_Float* vz, E_Float* result
)
{
  E_Float res1, res2, res3;
  E_Int nc = cn.getNConnect();

  res1 = 0.0;
  res2 = 0.0;
  res3 = 0.0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    
    E_Int ind1, ind2;
    E_Float fx, fy, fz;
    E_Float m1x, m2x;
    E_Float m1y, m2y;
    E_Float m1z, m2z;
    E_Float li;
    E_Float dx, dy, dz;

    for (E_Int i = 0; i < nelts; i++)
    {
      ind1 = cm(i, 1) - 1;
      ind2 = cm(i, 2) - 1;
      li = length[i];

      dx = xt[ind1] - cx;
      dy = yt[ind1] - cy;
      dz = zt[ind1] - cz;
      fx = ratio[ind1] * vx[ind1];
      fy = ratio[ind1] * vy[ind1];
      fz = ratio[ind1] * vz[ind1];

      m1x = dy * fz - dz * fy;
      m1y = dz * fx - dx * fz;
      m1z = dx * fy - dy * fx;

      dx = xt[ind2] - cx;
      dy = yt[ind2] - cy;
      dz = zt[ind2] - cz;
      fx = ratio[ind2] * vx[ind2];
      fy = ratio[ind2] * vy[ind2];
      fz = ratio[ind2] * vz[ind2];

      m2x = dy * fz - dz * fy;
      m2y = dz * fx - dx * fz;
      m2z = dx * fy - dy * fx;

      res1 += li * (m1x + m2x);
      res2 += li * (m1y + m2y);
      res3 += li * (m1z + m2z);
    }
  }

  result[0] = res1 * K_CONST::ONE_HALF;
  result[1] = res2 * K_CONST::ONE_HALF;
  result[2] = res3 * K_CONST::ONE_HALF;
}

// ============================================================================
// Compute surface integral of the moment M (OM^F)
// coordinates are defined in nodes and F is defined in center (unstructured)
// ============================================================================
void K_POST::integMomentUnstructNodeCenter(
  FldArrayI& cn, const char* eltType,
  const E_Float cx, const E_Float cy, const E_Float cz,
  const E_Float* ratio, const E_Float* xt, const E_Float* yt,
  const E_Float* zt, const E_Float* surf,
  const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float* result
)
{
  E_Float res1, res2, res3;
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  res1 = 0.0;
  res2 = 0.0;
  res3 = 0.0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    E_Int ind1, ind2, ind3;
    E_Float mx, my, mz, sri;
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

        centerx = K_CONST::ONE_THIRD * (xt[ind1] + xt[ind2] + xt[ind3]) - cx;
        centery = K_CONST::ONE_THIRD * (yt[ind1] + yt[ind2] + yt[ind3]) - cy;
        centerz = K_CONST::ONE_THIRD * (zt[ind1] + zt[ind2] + zt[ind3]) - cz;

        mx = centery * vz[i] - centerz * vy[i];
        my = centerz * vx[i] - centerx * vz[i];
        mz = centerx * vy[i] - centery * vx[i];

        sri = surf[i] * ratio[i];

        res1 += sri * mx;
        res2 += sri * my;
        res3 += sri * mz;
      }
    }
    else
    {
      fprintf(stderr, "Error: in K_POST::integMomentNormUnstructNodeCenter.\n");
      fprintf(stderr, "Unsupported type of element, %s.\n", eltTypes[ic]);
    }
  }

  result[0] = res1;
  result[1] = res2;
  result[2] = res3;

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

// ============================================================================
// Compute surface integral of the moment M (OM^F)
// coordinates are defined in nodes and F is defined in center (1D case)
// ============================================================================
void K_POST::integMomentUnstructNodeCenter1D(
  FldArrayI& cn, const char* eltType,
  const E_Float cx, const E_Float cy, const E_Float cz, const E_Float* ratio,
  const E_Float* xt, const E_Float* yt, const E_Float* zt, const E_Float* surf,
  const E_Float* vx, const E_Float* vy, const E_Float* vz,
  E_Float* result
)
{
  E_Float res1, res2, res3;
  E_Int nc = cn.getNConnect();

  res1 = 0.0;
  res2 = 0.0;
  res3 = 0.0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    E_Int ind1, ind2;
    E_Float mx, my, mz;
    E_Float centerx, centery, centerz;
    E_Float sri;
  
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();

    for (E_Int i = 0; i < nelts; i++)
    {
      ind1 = cm(i, 1) - 1;
      ind2 = cm(i, 2) - 1;

      centerx = K_CONST::ONE_HALF * (xt[ind1] + xt[ind2]) - cx;
      centery = K_CONST::ONE_HALF * (yt[ind1] + yt[ind2]) - cy;
      centerz = K_CONST::ONE_HALF * (zt[ind1] + zt[ind2]) - cz;

      mx = centery * vz[i] - centerz * vy[i];
      my = centerz * vx[i] - centerx * vz[i];
      mz = centerx * vy[i] - centery * vx[i];

      sri = surf[i] * ratio[i];

      res1 += sri * mx;
      res2 += sri * my;
      res3 += sri * mz;
    }
  }

  result[0] = res1;
  result[1] = res2;
  result[2] = res3;
}
