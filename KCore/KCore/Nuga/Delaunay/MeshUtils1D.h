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

#ifndef __DELAUNAY_MESH_UTILS_1D_H__
#define __DELAUNAY_MESH_UTILS_1D_H__

#include "Def/DefTypes.h"
#include "Fld/DynArray.h"
#include <vector>

namespace DELAUNAY
{


class MeshUtils1D
{
public:
  MeshUtils1D(void);
  ~MeshUtils1D(void);

  template <E_Int DIM>
  static void mesh_line (const E_Float* P0, const E_Float* P1, E_Int N,
    K_FLD::FloatArray& pos, K_FLD::IntArray& connect);

  template <E_Int DIM>
  static void mesh_arc (E_Float* C, E_Float* axis, E_Float* P0, E_Float alpha, E_Int N,
    K_FLD::FloatArray& pos, K_FLD::IntArray& connect);

  static void
    compute_iso_metric(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
    const std::vector<E_Int>& hard_nodes,
    std::vector<E_Float> & metric, E_Float &hmax, E_Float& hmin);
};

template <E_Int DIM>
void
MeshUtils1D::mesh_line
(const E_Float* P0, const E_Float* P1, E_Int N,
 K_FLD::FloatArray& pos, K_FLD::IntArray& connect)
{
  E_Float E[DIM];
  E_Float Pi[DIM];
  E_Int   Ei[2];

  K_FUNC::diff<DIM>(P1, P0, E);
  E_Float L = K_FUNC::normalize<DIM>(E);
  E_Float ti = L/N;

  pos.pushBack(P0, P0+DIM);
  E_Int I = pos.cols()-1;

  for (E_Int i = 0; i < N; ++i)
  {
    for (E_Int j = 0; j < DIM; ++j)
      Pi[j] = pos(j, I) + ti * E[j];

    pos.pushBack(Pi, Pi+DIM);
    Ei[0] = ++I-1;
    Ei[1] = I;
    connect.pushBack(Ei, Ei+2);
  }
}

template <E_Int DIM>
void
MeshUtils1D::mesh_arc (E_Float* C, E_Float* axis, E_Float* P0, E_Float alpha, E_Int N,
                       K_FLD::FloatArray& pos, K_FLD::IntArray& connect)
{
  E_Float R = ::sqrt(K_FUNC::sqrDistance(C, P0, DIM));
  E_Float ai = alpha*R/N;
  E_Float Pj[DIM], CPi[DIM], PiPj[DIM];
  const E_Float* Pi;
  E_Int Ei[2], I;

  pos.pushBack(P0, P0+DIM);
  I = pos.cols()-1;

  K_FUNC::normalize<3>(axis);

  for (E_Int i = 0; i < N; ++i)
  {
    Pi = pos.col(I);
    K_FUNC::diff<DIM>(Pi, C, CPi);
    K_FUNC::crossProduct<3>(axis, CPi, PiPj);
    K_FUNC::normalize<3>(PiPj);

    for (E_Int j = 0; j < DIM; ++j)
      Pj[j] = *(Pi+j) + ai * PiPj[j];

    pos.pushBack(Pj, Pj+DIM);

    Ei[0] = ++I-1;
    Ei[1] = I;

    connect.pushBack(Ei, Ei+2);
  }
}

}


#endif

