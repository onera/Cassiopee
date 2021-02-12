/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#include "Nuga/include/MeshUtils1D.h"
#include "Nuga/include/macros.h"
#include "Nuga/include/defs.h"
#include "Nuga/include/maths.hxx"
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

  hmin = NUGA::FLOAT_MAX;
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
    L = ::sqrt(NUGA::sqrDistance(pos.col(Ni), pos.col(Nj), DIM));
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
