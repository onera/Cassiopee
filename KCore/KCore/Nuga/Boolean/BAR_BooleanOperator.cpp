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

#include "Nuga/Boolean/BAR_BooleanOperator.h"
#include "Nuga/Boolean/BAR_Conformizer.h"
#include "Nuga/Delaunay/T3Mesher.h"
#include "Connect/EltAlgo.h"
#include "Nuga/GapFixer/FittingBox.h"

#ifdef DEBUG_BOOLEANBAR
#include "IO/io.h"
#endif

namespace NUGA
{

BAR_BooleanOperator::~BAR_BooleanOperator(void)
{
  delete _dT3;
}

BAR_BooleanOperator::BAR_BooleanOperator
(const K_FLD::FloatArray& coord1, const K_FLD::IntArray& cB1,
 const K_FLD::FloatArray& coord2, const K_FLD::IntArray& cB2,
 E_Float tolerance):parent_type(coord1, cB1, coord2, cB2, tolerance, new BAR_Conformizer<2>()), _dT3(0)
{
  _normal[0]=_normal[1]=_normal[2]=0.;
}

///B1 -  Intersection(B1, B2).
E_Int
BAR_BooleanOperator::get_1_minus_2(K_FLD::FloatArray& coord, K_FLD::IntArray& connect, std::vector<E_Int>& colors)
{
  if (!initialized()) return 1;
  coord = _coord;
  connect = _connect_1_out_2;
  //todo colors
  return 0;
}

///B2 -  Intersection(B1, B2).
E_Int
BAR_BooleanOperator::get_2_minus_1(K_FLD::FloatArray& coord, K_FLD::IntArray& connect, std::vector<E_Int>& colors)
{
  if (!initialized()) return 1;
  coord = _coord;
  connect = _connect_2_out_1;
  //todo colors
  return 0;
}

///  Sum(B1, B2) - Intersection(B1, B2).
E_Int
BAR_BooleanOperator::getUnion(K_FLD::FloatArray& coord, K_FLD::IntArray& connect, std::vector<E_Int>& colors)
{
  if (!initialized()) return 1;
  coord = _coord;
  connect = _connectUnion;
  //todo colors
  return 0;
}

/// Intersection between B1 & B2.
E_Int
BAR_BooleanOperator::getIntersection(K_FLD::FloatArray& coord, K_FLD::IntArray& connect, std::vector<E_Int>& colors)
{
  if (!initialized()) return 1;
  coord = _coord;
  connect = _connectInter;
  //todo colors
  return 0;
}

/// Intersection nodes
E_Int
BAR_BooleanOperator::getIntersectionBorder
(K_FLD::FloatArray& coord, K_FLD::IntArray& connect)
{
  connect.clear();
  coord.clear();
  
  K_FLD::IntArray conn = _connects[0];
  conn.pushBack(_connects[1]);
  
  std::vector<E_Int> nids;
  K_CONNECT::MeshTool::removeDuplicated(conn, nids, false/*strict orient*/);
  
  // Get the bound_to_elts map (to detect external boundaries
  typedef K_MESH::Edge::boundary_type node_type;
  typedef std::map<node_type, K_CONT_DEF::int_vector_type > BoundToEltType;
  BoundToEltType bound_to_elts;
  K_FLD::ArrayAccessor<K_FLD::IntArray> actv(conn);
  K_CONNECT::EltAlgo<K_MESH::Edge>::getBoundToElements(actv, bound_to_elts);
  
  E_Int DIM = _coord.rows();
  for (BoundToEltType::const_iterator it = bound_to_elts.begin();
                                      it != bound_to_elts.end(); ++it)
  {
    if (it->second.size() > 2)
    {
      const E_Int& Ni = it->first;
      coord.pushBack(_coord.col(Ni), _coord.col(Ni)+DIM);
    }
  }

  return 0;
  
}

///
E_Int
BAR_BooleanOperator::check_sanity()
{
  // Compute the normal to contours
  E_Float normal2[] = {0.,0.,0.};
  FittingBox::computeNormalToContour(_coord, _connects[0], _normal);
  FittingBox::computeNormalToContour(_coord, _connects[1], normal2);
  
  //check normals equality
  E_Float x = K_FUNC::dot<3>(_normal, normal2);
  
  if (::fabs(x + 1.) < E_EPSILON) // opposite orientation, one needs a reverse
  {
    for (E_Int i=0; i < _connects[1].cols(); ++i)
      std::swap(_connects[1](0,i), _connects[1](1,i));
    x = 1.;
  }
  
  if (::fabs(x - 1.) > E_EPSILON) return 1; //not on the same plane
  
  // check planarity
  E_Float Pt[3];
  for (E_Int i=0; i < _connects[0].cols(); ++i)
  {
    K_FUNC::diff<3>(_coord.col(_connects[0](1,i)), _coord.col(_connects[0](0,i)), Pt);
    E_Float x = K_FUNC::dot<3>(_normal, Pt);
    if ((x < -E_EPSILON)||(x > E_EPSILON)) return 1;
  }
  for (E_Int i=0; i < _connects[1].cols(); ++i)
  {
    K_FUNC::diff<3>(_coord.col(_connects[1](1,i)), _coord.col(_connects[1](0,i)), Pt);
    E_Float x = K_FUNC::dot<3>(_normal, Pt);
    if ((x < -E_EPSILON)||(x > E_EPSILON)) return 1;
  }
  
  return 0;
}

///
E_Int
BAR_BooleanOperator::compute_zones()
{
  K_FLD::FloatArray coord2D(_coord);
  
  // Go to the plane
  K_FLD::FloatArray P, iP;
  FittingBox::computeAFrame(_normal, P);
  iP = P;
  K_FLD::FloatArray::inverse3(iP);
  FittingBox::transform(coord2D, iP);// Transform to computed frame.
  
  coord2D.resize(2, coord2D.cols());
  
  //Triangulate the whole thing
  K_FLD::IntArray connect(_connects[0]);
  connect.pushBack(_connects[1]); 
  
  _dT3 = new DELAUNAY::MeshData(coord2D, connect);
  DELAUNAY::MesherMode mode;
  mode.mesh_mode = DELAUNAY::MesherMode::TRIANGULATION_MODE;
  //mode.remove_holes = false;
  DELAUNAY::T3Mesher<E_Float> mesher(mode);
  
  if (mesher.run(*_dT3))
    return 1;
  
#ifdef DEBUG_BOOLEANBAR
  MIO::write("triangulation.mesh", _dT3->pos, _dT3->connectM, "TRI", 0, &_dT3->colors);
#endif
  
  // Split the sub domains.
  const K_FLD::IntArray& connectT3 = _dT3->connectM;
  const std::vector<E_Int> &colors = _dT3->colors;
 
  std::map<E_Int, K_FLD::IntArray> connects;//use map instead of vector because colors can be arbitrary values (when the Mesher2D calls clean_data, some colors are removed)
  K_FLD::IntArray::const_iterator pS;
  for (E_Int i = 0; i < connectT3.cols(); ++i)
  {
    pS = connectT3.col(i);
    connects[colors[i]].pushBack(pS, pS+3);
  }
    
  // Store the edges in h sets : hB global outer boundary, h[0] BAR1 boundary, h[1] for BAR2
  typedef K_MESH::Edge edge_type;
  std::set<edge_type> hB, hBu, h[2];
  std::vector<edge_type> vB12;
  for (size_t j = 0; j < 2; ++j)
  {
    for (E_Int i = 0; i < _connects[j].cols(); ++i)
    {
      pS =  _connects[j].col(i);
      h[j].insert(edge_type(pS));
    }
  }
  
  K_CONNECT::MeshTool::getBoundary(connectT3, _connectUnion);
  
#ifdef DEBUG_BOOLEANBAR
  MIO::write("union.mesh", _coord, _connectUnion, "BAR");
#endif
  
  convert_to_hset(_connectUnion, hBu);
  
  // Loop through sub domains and dispatch edges to the appropriate container
  K_FLD::IntArray cB;
  for (std::map<E_Int, K_FLD::IntArray>::const_iterator i = connects.begin(); i != connects.end(); ++i)
  {
    const K_FLD::IntArray& connect = i->second;
    K_CONNECT::MeshTool::getBoundary(connect, cB);
    convert_to_hset(cB, hB);
    vB12.clear();
    std::set_difference(hB.begin(), hB.end(), hBu.begin(), hBu.end(), std::back_inserter(vB12));
    
    if (vB12.size() > 0 && vB12.size() < hB.size())
    {
      edge_type ne(vB12[0].node(1), vB12[0].node(0));
      if (h[0].find(ne) != h[0].end()) // opposite belongs to 1 => out 2
        _connect_2_out_1.pushBack(cB);
      else if (h[1].find(ne) != h[1].end()) // opposite belongs to 2 => out 1 
        _connect_1_out_2.pushBack(cB);
      else
        _connectInter.pushBack(cB);
    }
    else // cB is a totally enclosed contour.
    {
      //assert (vB12.size() == hB.size());
      _connectInter.pushBack(cB);
    }
  }
  return 0;
}

///
void
BAR_BooleanOperator::convert_to_hset(const K_FLD::IntArray& cB, std::set<K_MESH::Edge>& hB)
{
  hB.clear();
  K_FLD::IntArray::const_iterator   pS;
  
  for (E_Int i = 0; i < cB.cols(); ++i)
  {
    pS =  cB.col(i);
    hB.insert(K_MESH::Edge(pS));
  }
}
  
  
}

