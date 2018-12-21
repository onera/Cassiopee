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

#include "MeshUtils1D.h"
#include "macros.h"
#include <algorithm>
#include <set>

namespace DELAUNAY{

MeshUtils1D::MeshUtils1D(void)
{
}

MeshUtils1D::~MeshUtils1D(void)
{
}

void
MeshUtils1D::compute_iso_metric
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
 const std::vector<E_Int>& hard_nodes,
 std::vector<E_Float> & metric, E_Float &hmin, E_Float& hmax)
{
  /// Compute the metric at connect nodes and assigne the mean to other hard nodes.

  std::vector<E_Int> nodes;
  std::vector<E_Int> count;
  K_FLD::IntArray::const_iterator pS;
  E_Int Ni, Nj, DIM(pos.rows()), NBOUND(connect.cols()), bnods;
  E_Float L;

  hmin = K_CONST::E_MAX_FLOAT;
  hmax = -hmin;

  connect.uniqueVals(nodes);
  bnods = nodes.size();
  E_Int max = *std::max_element(ALL(nodes));
  if (!hard_nodes.empty())
    max = std::max(max, *std::max_element(ALL(hard_nodes)));
  
  metric.clear();
  metric.resize(max+1, 0.);
  count.resize(max+1, 0);

  for (E_Int i = 0; i < NBOUND; ++i)
  {
    pS = connect.col(i);
    Ni = *pS;
    Nj = *(pS+1);
    L = ::sqrt(K_FUNC::sqrDistance(pos.col(Ni), pos.col(Nj), DIM));
    metric[Ni] += L;
    metric[Nj] += L;
    ++count[Ni];
    ++count[Nj];
  }

  for (E_Int i = 0; i < bnods; ++i)
  {
    Ni = nodes[i];
    if (count[Ni] != 0)metric[Ni] /= count[Ni];
    if (metric[Ni] != 0.)
      hmin = std::min(metric[Ni], hmin);
    hmax = std::max(metric[Ni], hmax);
  }

  //hmin = *std::min_element(ALL(metric));
  //hmax = *std::max_element(ALL(metric));

  // Fixme : set the mean value at hard nodes that are not boundary nodes.
  E_Float hmean = 0.5*(hmin+hmax);
  for (size_t i = 0; i < metric.size(); ++i)
    metric[i] = (metric[i] == 0. ) ? hmean : metric[i];

}
}
