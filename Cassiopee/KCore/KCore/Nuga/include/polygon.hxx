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
//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef NUGA_POLYGON_HXX
#define NUGA_POLYGON_HXX

#include "Nuga/include/Polygon.h"

namespace NUGA
{
///autonomous (i.e. self contained) polygon
struct aPolygon : public K_MESH::Polygon
{
  using parent_type = K_MESH::Polygon;
  
  std::vector<E_Int> m_nodes;
  K_FLD::FloatArray  m_crd;
  std::vector<long>  m_poids;
  E_Float            m_Lref2;
  mutable E_Float    m_normal[3];
  mutable E_Float    m_centroid[3];
    
  aPolygon():parent_type(nullptr, 0)
  {
    /*m_nodes.resize(1); //hack to have an allocated place and hence an address to pass to _nodes
    parent_type::_nodes = &m_nodes[0];
    m_nodes.clear(); //end of hack : clear this dummy value*/

    // init remaining stuff
    _triangles = nullptr;
    _shift = 0;
    m_normal[0] = m_centroid[0] = NUGA::FLOAT_MAX;
  }
  
  // from "mesh" to autonmous
  aPolygon(const E_Int* nodes, E_Int nb_nodes, E_Int idx_start, const K_FLD::FloatArray& crd):parent_type(nodes, nb_nodes, -idx_start)
  {
    _nb_nodes = nb_nodes;
 
    // keeping track of node history before compressing for autonomy
    m_poids.clear();
    m_poids.resize(_nb_nodes, IDX_NONE);
    for (E_Int n = 0; n < _nb_nodes; ++n)m_poids[n] = nodes[n] - idx_start;

    //compress
    NUGA::MeshTool::compact_to_mesh(crd, nodes, nb_nodes, idx_start, m_crd);

    //autonomous values & plugging
    m_nodes.clear();
    K_CONNECT::IdTool::init_inc(m_nodes, nb_nodes, 0);

    _nodes = &m_nodes[0];

    // init remaining stuff
    _triangles = nullptr;
    _shift = 0;
    m_normal[0] = m_centroid[0] = NUGA::FLOAT_MAX;

  }

  // from "mesh" to autonmous (2)
  aPolygon(const parent_type& pg, const K_FLD::FloatArray& crd):parent_type(pg.begin(), pg.nb_nodes(), pg.shift()), m_Lref2(-1.)
  { 
    // keeping track of node history before compressing for autonomy
    m_poids.clear();
    m_poids.resize(_nb_nodes, IDX_NONE);
    for (E_Int n = 0; n < _nb_nodes; ++n)m_poids[n] = _nodes[n] + _shift;
    
    //compress
    NUGA::MeshTool::compact_to_mesh(crd, _nodes, _nb_nodes, -_shift, m_crd);
    
    //autonomous values & plugging
    m_nodes.clear();
    K_CONNECT::IdTool::init_inc(m_nodes, pg.nb_nodes(), 0);
    
    parent_type::_nb_nodes = (E_Int)m_nodes.size();
    parent_type::_nodes = &m_nodes[0];
    
    // init remaining stuff
    _triangles = nullptr;
    parent_type::_shift = 0;
    m_normal[0] = m_centroid[0] = NUGA::FLOAT_MAX;
  }
  
  aPolygon(const parent_type& pg, const K_FLD::FloatArray& crd, E_Float L2r) :aPolygon(pg, crd) { m_Lref2 = L2r; m_normal[0] = m_centroid[0] = NUGA::FLOAT_MAX; }

  aPolygon(K_FLD::FloatArray && crd):parent_type(nullptr, 0), m_crd(std::move(crd))
  {
    //autonomous values & plugging
    m_nodes.clear();
    K_CONNECT::IdTool::init_inc(m_nodes, m_crd.cols(), 0);

    parent_type::_nb_nodes = m_nodes.size();
    parent_type::_nodes = &m_nodes[0];

    // init remaining stuff
    _triangles = nullptr;
    parent_type::_shift = 0;
    m_normal[0] = m_centroid[0] = NUGA::FLOAT_MAX;
  }

  aPolygon(const ngon_unit& ngu, E_Int ith, const K_FLD::FloatArray& crd): parent_type(nullptr, 0)
  {
    parent_type::_nb_nodes = ngu.stride(ith);
    const E_Int* p = ngu.get_facets_ptr(ith);

    // keeping track of node history before compressing for autonomy
    m_poids.clear();
    m_poids.resize(_nb_nodes, IDX_NONE);
    for (E_Int n = 0; n < _nb_nodes; ++n)m_poids[n] = p[n] - 1;

    //compress
    NUGA::MeshTool::compact_to_mesh(crd, ngu.get_facets_ptr(ith), parent_type::_nb_nodes, 1, m_crd);
    
    //autonomous values & plugging
    m_nodes.clear();
    K_CONNECT::IdTool::init_inc(m_nodes, parent_type::_nb_nodes, 0);

    parent_type::_nodes = &m_nodes[0];

    // init remianing stuff
    _triangles = nullptr;
    parent_type::_shift = 0;
    m_normal[0] = m_centroid[0] = NUGA::FLOAT_MAX;
  }

  aPolygon(const aPolygon& rhs) : parent_type(nullptr, 0)
  {
    *this = rhs;
  }
    
  void clear() { m_nodes.clear();}
  bool empty() { return m_nodes.empty(); }
  
  aPolygon& operator=(const parent_type& rhs) = delete;
  aPolygon& operator=(const aPolygon& rhs)
  {
    m_nodes = rhs.m_nodes;
    m_crd = rhs.m_crd;
    m_poids = rhs.m_poids;
    m_Lref2 = rhs.m_Lref2;

    parent_type::_nb_nodes = m_nodes.size();
    parent_type::_nodes = &m_nodes[0];
    _triangles = nullptr;
    parent_type::_shift = 0;

    m_normal[0] = rhs.m_normal[0];
    m_normal[1] = rhs.m_normal[1];
    m_normal[2] = rhs.m_normal[2];

    m_centroid[0] = rhs.m_centroid[0];
    m_centroid[1] = rhs.m_centroid[1];
    m_centroid[2] = rhs.m_centroid[2];

    return *this;
  }

  aPolygon(aPolygon&& rhs) :/* = default; rejected by old compiler intel (15)*/
  parent_type(rhs), m_nodes(std::move(rhs.m_nodes)), m_crd(std::move(rhs.m_crd)), m_poids(std::move(rhs.m_poids)), m_Lref2(rhs.m_Lref2)
  {
    m_normal[0] = rhs.m_normal[0];
    m_normal[1] = rhs.m_normal[1];
    m_normal[2] = rhs.m_normal[2];

    m_centroid[0] = rhs.m_centroid[0];
    m_centroid[1] = rhs.m_centroid[1];
    m_centroid[2] = rhs.m_centroid[2];
  }

  aPolygon& operator=(aPolygon&& rhs)
  {
    m_nodes = std::move(rhs.m_nodes);
    m_crd = std::move(rhs.m_crd);
    m_poids = std::move(rhs.m_poids);
    m_Lref2 = rhs.m_Lref2;

    parent_type::_nb_nodes = m_nodes.size();
    parent_type::_nodes = &m_nodes[0];
    _triangles = rhs._triangles;
    parent_type::_shift = rhs._shift;

    m_normal[0] = rhs.m_normal[0];
    m_normal[1] = rhs.m_normal[1];
    m_normal[2] = rhs.m_normal[2];

    m_centroid[0] = rhs.m_centroid[0];
    m_centroid[1] = rhs.m_centroid[1];
    m_centroid[2] = rhs.m_centroid[2];

    rhs._nb_nodes = 0;
    rhs._shift = 0;
    rhs._triangles = nullptr;
    rhs._nodes = nullptr;
    rhs.m_normal[0] = NUGA::FLOAT_MAX;//fixme: necessary ?

    return *this;
  }

  bool operator==(const aPolygon& rhs)
  {
    int nnodes = nb_nodes();
    if (nnodes != rhs.nb_nodes())     return false;
    if (m_crd.cols() != rhs.m_crd.cols()) return false;

    if ((m_normal[0] != NUGA::FLOAT_MAX) && (rhs.m_normal[0] != NUGA::FLOAT_MAX))
    {
      if (::fabs((m_normal[0] - rhs.m_normal[0])) > EPSILON) return false;
      if (::fabs((m_normal[1] - rhs.m_normal[1])) > EPSILON) return false;
      if (::fabs((m_normal[2] - rhs.m_normal[2])) > EPSILON) return false;
    }

    if ((m_centroid[0] != NUGA::FLOAT_MAX) && (rhs.m_centroid[0] != NUGA::FLOAT_MAX))
    {
      if (::fabs((m_centroid[0] - rhs.m_centroid[0])) > EPSILON) return false;
      if (::fabs((m_centroid[1] - rhs.m_centroid[1])) > EPSILON) return false;
      if (::fabs((m_centroid[2] - rhs.m_centroid[2])) > EPSILON) return false;
    }

    // check for vertices coincidence
    const double* pt0 = m_crd.col(0);
    E_Int n0 = IDX_NONE;
    double dmin2 = NUGA::FLOAT_MAX;
    for (E_Int n = 0; n < nnodes; ++n)
    {
      const double* rptn = rhs.m_crd.col(n);

      double d2 = (pt0[0] - rptn[0])*(pt0[0] - rptn[0]) + (pt0[1] - rptn[1])*(pt0[1] - rptn[1]) + (pt0[2] - rptn[2])*(pt0[2] - rptn[2]);
      if (d2 < dmin2)
      {
        n0 = n;
        dmin2 = d2;
      }
    }

    if (dmin2 > EPSILON*EPSILON) return false; // none coincident node

    for (E_Int n = 1; n < nnodes; ++n)
    {
      const double* pt  = m_crd.col(n);
      const double* rpt = rhs.m_crd.col((n0 + n) % nnodes);

      double d2 = (pt[0] - rpt[0])*(pt[0] - rpt[0]) + (pt[1] - rpt[1])*(pt[1] - rpt[1]) + (pt[2] - rpt[2])*(pt[2] - rpt[2]);
      if (d2 > EPSILON*EPSILON) return false;
    }

    // each node has exactly one coincident node in rhs
    return true;
  }
  
  void reverse_orient()
  {
    std::reverse(ALL(m_nodes));
    if (m_normal[0] != NUGA::FLOAT_MAX)
    {
      m_normal[0] *= -1.;
      m_normal[1] *= -1.;
      m_normal[2] *= -1.;
    }
  }
  
  template <short DIM> void normal(double* norm) const 
  {
    parent_type::normal<K_FLD::FloatArray, DIM>(m_crd, norm);
  }

  template <short DIM> void centroid(double* G) const
  {
    parent_type::centroid<DIM>(m_crd, _nodes, _nb_nodes, 0, G);
  }
  
  template <typename TriangulatorType>
  E_Int triangulate(const TriangulatorType& dt) const
  {
    return parent_type::triangulate(dt, m_crd);
  }

  double extent() const
  {
    return parent_type::surface<K_FLD::FloatArray, 3>(m_crd, parent_type::_nodes, parent_type::_nb_nodes, parent_type::_shift);
  }

  double metrics() const
  {
    normal<3>(m_normal);
    return extent();
  }

  const double* get_normal() const { 
    if (m_normal[0] == NUGA::FLOAT_MAX)
      normal<3>(m_normal);
    return m_normal;
  }

  const double* get_centroid() const {
    if (m_centroid[0] == NUGA::FLOAT_MAX)
      parent_type::centroid<3>(m_crd, parent_type::_nodes, parent_type::_nb_nodes, parent_type::_shift, m_centroid);
    return m_centroid;
  }

  double Lref2() const { return (m_Lref2 > 0.) ? m_Lref2 : parent_type::Lref2(m_crd);} // if passed by mesh_t, return it, otherwise compute it first
};

}

#endif /* POLYGON_HXX */

