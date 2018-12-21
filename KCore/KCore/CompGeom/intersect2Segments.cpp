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
# include "CompGeom/compGeom.h"

using namespace K_FLD;
using namespace K_FUNC;
using namespace std;

#define SMALL_NUM 0.00000001
#define dot(u,v) ( u[0]*v[0] + u[1]*v[1] + u[2]*v[2] )
#define perp(u,v) (u[0]*v[1] - u[1]*v[0])
#define dist(P1,P2) ((P1[0]-P2[0])*(P1[0]-P2[0])+(P1[1]-P2[1])*(P1[1]-P2[1])+(P1[2]-P2[2])*(P1[2]-P2[2]))


//=============================================================================
E_Int K_COMPGEOM::intersect2Segments(E_Float* ps1a, E_Float* ps1b,
                                     E_Float* ps2a, E_Float* ps2b,
                                     E_Float* pi0, E_Float* pi1)
{
  E_Float u[3]; E_Float v[3]; E_Float w[3];
  u[0] = ps1b[0] - ps1a[0];
  u[1] = ps1b[1] - ps1a[1];
  u[2] = ps1b[2] - ps1a[2];
  v[0] = ps2b[0] - ps2a[0];
  v[1] = ps2b[1] - ps2a[1];
  v[2] = ps2b[2] - ps2a[2];
  w[0] = ps1a[0] - ps2a[0];
  w[1] = ps1a[1] - ps2a[1];
  w[2] = ps1a[2] - ps2a[2];
  E_Float D = perp(u, v);

  // Test si les segments sont paralleles
  if (E_abs(D) < SMALL_NUM) // segments paralleles
  {
    if (fEqualZero(perp(u,w)) == false ||
        fEqualZero(perp(v,w)) == false) return 0;
    E_Float du = dot(u,u);
    E_Float dv = dot(v,v);
    if (fEqualZero(du) == true && 
        fEqualZero(dv) == true)
    {
      if (fEqualZero(dist(ps1a,ps2a)) == false) return 0;
      pi0[0] = ps1a[0]; pi0[1] = ps1a[1]; pi0[2] = ps1a[2];
      return 1;
    }
    if (fEqualZero(du) == true)
    {
      if (inSegment(ps1a, ps2a, ps2b) == 0) return 0;
      pi0[0] = ps1a[0]; pi0[1] = ps1a[1]; pi0[2] = ps1a[2];
      return 1;
    }
    if (fEqualZero(dv) == true)
    {
      if (inSegment(ps2a, ps1a, ps1b) == 0) return 0;
      pi0[0] = ps2a[0]; pi0[1] = ps2a[1]; pi0[2] = ps2a[2];
      return 1;
    }
    
    // Ils sont dans le meme plan
    E_Float t0, t1;
    E_Float w2[3];
    w2[0] = ps1b[0] - ps2a[0];
    w2[1] = ps1b[1] - ps2a[1];
    w2[2] = ps1b[2] - ps2a[2];
    if (fEqualZero(v[0]) == false)
    {
      t0 = w[0] / v[0];
      t1 = w2[0] / v[0];
    }
    else if (fEqualZero(v[1]) == false)
    {
      t0 = w[1] / v[1];
      t1 = w2[1] / v[1];
    }
    else
    {
      t0 = w[2] / v[2];
      t1 = w2[2] / v[2];
    }
    if (t0 > t1)
    {
      E_Float t = t0; t0 = t1; t1 = t;
    }
    if (t0 > 1 || t1 < 0) return 0; // No overlap
    t0 = t0 < 0 ? 0 : t0;
    t1 = t1 > 1 ? 1 : t1;
    if (fEqualZero(t0-t1) == true)
    {
      pi0[0] = ps2a[0] + t0 * v[0];
      pi0[1] = ps2a[1] + t0 * v[1];
      pi0[2] = ps2a[2] + t0 * v[2];
      return 1;
    }
    pi0[0] = ps2a[0] + t0*v[0];
    pi0[1] = ps2a[1] + t0*v[1];
    pi0[2] = ps2a[2] + t0*v[2];
    pi1[0] = ps2a[0] + t1*v[0];
    pi1[1] = ps2a[1] + t1*v[1];
    pi1[2] = ps2a[2] + t1*v[2];
    return 2;
  }

  // Les segments sont skew
  E_Float si = perp(v,w) / D;
  if (si < 0 || si > 1) return 0;
  
  E_Float ti = perp(u,w) / D;
  if (ti < 0 || ti > 1) return 0;
  pi0[0] = ps1a[0] + si * u[0];
  pi0[1] = ps1a[1] + si * u[1];
  pi0[2] = ps1a[2] + si * u[2];
  return 1;
}

//=============================================================================
E_Int K_COMPGEOM::inSegment(E_Float* p, E_Float* psa, E_Float* psb)
{
  if (fEqualZero(psa[0]-psb[0]) == false)
  {
    if (psa[0] <= p[0] && p[0] <= psb[0]) return 1;
    if (psa[0] >= p[0] && p[0] >= psb[0]) return 1;
  }
  else
  {
    if (psa[1] <= p[1] && p[1] <= psb[1]) return 1;
    if (psa[1] >= p[1] && p[1] >= psb[1]) return 1;
  }
  return 0;
}
