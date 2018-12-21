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

#ifndef __DELAUNAY_REFINER_H__
#define __DELAUNAY_REFINER_H__

#include "Fld/DynArray.h"
#include "Def/DefContainers.h"
#include "Search/KdTree.h"
#include "MeshData.h"
#include "macros.h"

#ifdef DEBUG_METRIC
#include "Linear/DelaunayMath.h"
#include "iodata.h"
#include "IO/io.h"
#endif
#include "GeomMetric.h"

namespace DELAUNAY
{

  template <typename MetricType>
  class Refiner
  {
  public:
    typedef K_CONT_DEF::size_type                  size_type;
    typedef K_CONT_DEF::int_vector_type            int_vector_type;
    typedef K_CONT_DEF::int_set_type               int_set_type;
    typedef K_SEARCH::KdTree<>                     tree_type;
    typedef K_MESH::Triangle                       element_type;
    typedef K_CONT_DEF::non_oriented_edge_set_type non_oriented_edge_set_type;

  public:

    Refiner(MetricType& metric, E_Float growth_ratio, E_Int nb_smooth_iter, E_Bool symmetrize):_metric(metric), _threshold(0.5), 
            _gr(growth_ratio), _nb_smooth_iter(nb_smooth_iter), _symmetrize(symmetrize), _debug(false){}

    ~Refiner(void){}

    void computeRefinePoints(E_Int iter, MeshData& data, const int_set_type& box_nodes,
                             const non_oriented_edge_set_type& hard_edges,
                             int_vector_type& refine_nodes, E_Int N0/*for metrics changes when smoothing*/);

    void filterRefinePoints(MeshData& data, const int_set_type& box_nodes,
                            int_vector_type& refine_nodes,
                            tree_type& filter_tree);

  private:
    MetricType&             _metric;
    std::vector<size_type>  _tmpNodes;
    E_Float                 _threshold;
    E_Float                 _gr;
    E_Int                   _nb_smooth_iter;
    E_Bool                  _symmetrize;
  public:
    E_Bool                  _debug;    
  };

  ///
  template <typename MetricType>
  void
    Refiner<MetricType>::computeRefinePoints
    (E_Int iter, MeshData& data, const int_set_type& box_nodes,
     const non_oriented_edge_set_type& hard_edges,
     int_vector_type& refine_nodes, E_Int N0/*for metrics changes when smoothing*/)
  {
    refine_nodes.clear();

    std::set<K_MESH::NO_Edge> all_edges;

    K_FLD::IntArray::const_iterator pS;

    // Get all the inner edges (non-hard).
    size_type cols = data.connectM.cols();
    for (size_type j = 0; j < cols; ++j)
    {
      if (data.colors[j] == 0) continue;

      pS = data.connectM.col(j);

      for (size_type i = 0; i < element_type::NB_NODES; ++i)
      {
        const size_type& Ni = *(pS+i);
        const size_type& Nj = *(pS + (i+1) % element_type::NB_NODES);

        if (IS_IN(box_nodes, Ni) || IS_IN(box_nodes, Nj))
          continue;

        K_MESH::NO_Edge Ei(Ni, Nj);

        if (IS_IN(hard_edges, Ei))
          continue;

        all_edges.insert(Ei);
      }
    }

    //smooth the metric at each nodes.
    // so do it at leat once whatever the user ask for to impovre the overall
    // mesh quality by setting the right metrics at bone nodes (__init_refine_points)
    _metric._N0 = N0;

    bool do_smooth  = (_nb_smooth_iter > 0) && (iter > 1);
         do_smooth |= (_symmetrize && (iter == 1)); // T3 Mesher use : always smooth at first iter for skeleton nodes.

    if(do_smooth)
    {
    
#ifdef DEBUG_METRIC
      {
        std::ostringstream o;
        o << "ellipse_beforesmooth_iter_" << iter << ".mesh";
        _metric.draw_ellipse_field(o.str().c_str(), *data.pos, data.connectM, &data.mask);
      }
#endif
      
      _metric.smoothing_loop(all_edges, _gr, _nb_smooth_iter, N0);

#ifdef DEBUG_METRIC
      {
        std::ostringstream o;
        o << "ellipse_aftersmooth_iter_" << iter << ".mesh";
        _metric.draw_ellipse_field(o.str().c_str(), *data.pos, data.connectM, &data.mask);
      }
#endif
    
    }

    std::vector<std::pair<E_Float, size_type> > length_to_points;
    if (_symmetrize && iter== 0)
      // Compute the bone mesh (for a 2D mesh, not a Geom Mesh)
      for (const auto& Ei : all_edges)
        _metric.__init_refine_points(*data.pos, Ei.node(0), Ei.node(1), _threshold, length_to_points, _tmpNodes);
    else
      for (const auto& Ei : all_edges)
        _metric.__compute_refine_points(*data.pos, Ei.node(0), Ei.node(1), _threshold, length_to_points, _tmpNodes);

    std::sort(ALL(length_to_points));

    size_type sz = (size_type)length_to_points.size();
    for (size_type i = 0; i < sz; ++i)
      refine_nodes.push_back(length_to_points[i].second);

  }

  ///
  template <typename MetricType>
  void
    Refiner<MetricType>::filterRefinePoints
    (MeshData& data, const int_set_type& box_nodes, 
    int_vector_type& refine_nodes, tree_type& filter_tree)
  {
    size_type Ni, sz, nb_nodes;
    E_Float mBox[2], MBox[2];
    E_Float coeff(0.5*::sqrt(2.)/*fixme sqrt...*/), Ri;
    int_vector_type nodes;
    int_vector_type tmp(refine_nodes);
    refine_nodes.clear();

    sz = (size_type) tmp.size();
    for (size_type i = 0; i < sz; ++i)
    {
      Ni = tmp[i];
      const typename MetricType::value_type& Mi = _metric[Ni];
      Ri = coeff * _metric.getRadius(Ni);//fixme

      for (size_type j = 0; j < 2; ++j)
      {
        mBox[j] = (*data.pos)(j, Ni) - Ri;
        MBox[j] = (*data.pos)(j, Ni) + Ri;
      }

      nodes.clear();
      filter_tree.getInBox(mBox, MBox, nodes);

      bool discard = false;
      nb_nodes = (size_type)nodes.size();
      for (size_type n = 0; (n < nb_nodes) && !discard; ++n)
      {
        if (IS_IN(box_nodes, nodes[n]))// skip box nodes
          continue;

        const typename MetricType::value_type& Mn = _metric[nodes[n]];

        E_Float dmi = _metric.lengthEval(Ni, Mi, nodes[n], Mi);
        E_Float dmn = _metric.lengthEval(Ni, Mn, nodes[n], Mn);
        discard = (dmi < coeff) || (dmn < coeff); // Regarding Mi or Mn metric

        
#ifdef DEBUG_METRIC
        /*if (discard && Ni == 19551)
        {
          K_FLD::FloatArray crdo;
          K_FLD::IntArray cnto;
          
          E_Int E[] = {0,1};
          crdo.pushBack(data.pos->col(Ni), data.pos->col(Ni) + 2);
          crdo.pushBack(data.pos->col(nodes[n]), data.pos->col(nodes[n]) + 2);
          
          if (dmi < coeff)
            _metric.append_unity_ellipse(*data.pos, Ni, crdo, cnto);
          if (dmn < coeff)
            _metric.append_unity_ellipse(*data.pos, nodes[n], crdo, cnto);
          
          cnto.pushBack(E, E+2);
          std::ostringstream o;
          o << "discarded_" << Ni << ".mesh";
          MIO::write(o.str().c_str(), crdo, cnto, "BAR");     
        }*/
#endif
      }

      if (!discard)
      {
        refine_nodes.push_back(Ni);
        filter_tree.insert(Ni);
      }
    }
  }

/*
  ///
  template <>
  void
  Refiner<E_Float>::__compute_refine_points
  (size_type Ni, size_type Nj, std::vector<std::pair<E_Int, size_type> >& length_to_points)
  {
  E_Float d;
  E_Float NiNj[2];
  K_FUNC::diff<2>(_pos.col(Nj), _pos.col(Ni), NiNj);
  // Geom repartition.fixme
  E_Float h0 = _metric[Ni], h1 = _metric[Nj];
  E_Float d0 = ::sqrt(K_FUNC::sqrDistance(_pos.col(Ni), _pos.col(Nj), _pos.rows()));
  NiNj[0] /= d0;
  NiNj[1] /= d0;//Normalize
  E_Float x(h1 / h0), r1 = x;
  E_Float  r = 1;
  size_type n;
  if (::abs(r1 - 1.) < 0.01)
  d = d0 / std::max(h0, h1);
  else
  {
  r = (d0 + h1) / (d0 + h0);
  d = ::log(r1)/::log(r);
  //n = size_type(d);
  //r1 = ::pow (r1, 1./n);
  //d =  (h1 - r1 * h0)/(r1 - 1.);
  //r = r1;
  }
  n = std::max(size_type(d), 1);

  if ((n * (n + 1)) < (d * d))
  ++n;

  if (n == 1) // Saturated
  return;

  E_Float alpha = r*h0;
  size_type Nstart = Ni;
  E_Float newP[2], s, h;

  for (size_type i = 0; i < n-1; ++i)
  {
  newP[0] = _pos(0, Nstart) + alpha * NiNj[0];
  newP[1] = _pos(1, Nstart) + alpha * NiNj[1];

  _pos.pushBack (newP, newP+2);
  Nstart = _pos.cols()-1;
  s = (i+1.)/(n+1.);
  h = h0 * ::pow (x, s);//(1 - s) * h0 + s * h1;
  _metric.push_back (h);
  alpha *= r;
  length_to_points.push_back(std::make_pair(-n, Nstart));
  }
  }
  */

}

#endif

