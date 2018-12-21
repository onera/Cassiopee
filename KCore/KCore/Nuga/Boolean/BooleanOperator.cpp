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
//Author : Sâm Landier (sam.landier@onera.fr)

#include "BooleanOperator.h"
#include "Connect/MeshTool.h"

#ifdef DEBUG_BOOLEAN
#include "IO/io.h"
#endif

namespace NUGA
{
  
//
BooleanOperator::BooleanOperator
(const K_FLD::FloatArray& coord1, const K_FLD::IntArray& connect1,
 const K_FLD::FloatArray& coord2, const K_FLD::IntArray& connect2, E_Float tolerance,
 ConformizerRoot* c):_conformizer(c), _initialized(false)
{
  E_Int                             ret(1), ROWS(connect1.rows());
  std::vector<E_Int>                ancestors;
  K_FLD::FloatArray                 coord(coord1);
  K_FLD::IntArray                   connect(connect1);
  K_FLD::IntArray::const_iterator   pS;
  std::vector<E_Int>                color(connect1.cols(), 0);
  
  //Fast return.
  if (coord1.cols() * connect1.cols() * coord2.cols() * connect2.cols()== 0)
    return; // One input container is empty.
  
#ifdef DEBUG_BOOLEAN
  std::string EType = "TRI";
  if (connect1.rows() == 2)
    EType = "BAR";
#endif

#ifdef DEBUG_BOOLEAN
  {
    E_Int mini, maxi;
    K_CONNECT::MeshTool::computeMinMaxIndices(connect1, mini, maxi);
    assert (mini >= 0);
    assert (maxi < coord1.cols() || maxi == E_IDX_NONE);
  }
    {
    E_Int mini, maxi;
    K_CONNECT::MeshTool::computeMinMaxIndices(connect2, mini, maxi);
    assert (mini >= 0);
    assert (maxi < coord2.cols() || maxi == E_IDX_NONE);
  }
#endif
  
  {
    K_FLD::IntArray connect_tmp(connect2);
    connect_tmp.shift(coord1.cols());
    connect.pushBack(connect_tmp);
    color.resize(connect.cols(), 1);
  }

  coord.pushBack(coord2);
  
  //E_Float zero = 0.;
  //coord.resize(DIM, coord.cols(), &zero);
    
  std::vector<E_Int> nids;
  K_CONNECT::MeshTool::compact_to_mesh(coord, connect, nids);
 
#ifdef DEBUG_BOOLEAN
  {
    E_Int mini, maxi;
    K_CONNECT::MeshTool::computeMinMaxIndices(connect, mini, maxi);
    assert (mini >= 0);
    assert (maxi < coord.cols() || maxi == E_IDX_NONE);
  }
#endif
 
  ret = _conformizer->run(coord, connect, ancestors, 0, tolerance, connect1.cols()/*X0*/);
  if (ret) return;
  
#ifdef DEBUG_BOOLEAN
  MIO::write("conformized.mesh", coord, connect, EType.c_str(), 0, &ancestors);
#endif

  _coord = coord;

  // Split the 2 parts.
  for (E_Int i = 0; i < connect.cols(); ++i)
  {
    pS = connect.col(i);
    _connects[color[ancestors[i]]].pushBack(pS, pS+ROWS);
  }
  
  if (_connects[0].cols()*_connects[1].cols() == 0)
    return; // Errror : at least one is empty.
  
#ifdef DEBUG_BOOLEAN
  //MIO::write("connect1split0.mesh", _coord, _connects[0], EType.c_str());
  //MIO::write("connect2split0.mesh", _coord, _connects[1], EType.c_str());
#endif
  
  // Remove eventual duplicates in each part individually
  {
    std::vector<E_Int> dumIds;
    K_CONNECT::MeshTool::removeDuplicated(_connects[0], dumIds, false/*strict orient*/);
    K_CONNECT::MeshTool::removeDuplicated(_connects[1], dumIds, false/*strict orient*/);
  }

#ifdef DEBUG_BOOLEAN
  MIO::write("connect1split.mesh", _coord, _connects[0], EType.c_str());
  MIO::write("connect2split.mesh", _coord, _connects[1], EType.c_str());
#endif
}

///
E_Int
BooleanOperator::initialized()
{
  if (_initialized)
    return 1;
  
  if (this->check_sanity())
    return 0;
  
  if (this->compute_zones())
    return 0;
  
  _initialized = true;
  return 1;
}

///
E_Int
BooleanOperator::getSum
(K_FLD::FloatArray& coord, K_FLD::IntArray& connect, std::vector<E_Int>& colors)
{
  if (!initialized()) return 1;

  coord = _coord;
  connect = _connects[0];
  connect.pushBack(_connects[1]);

  colors.clear();
  colors.resize(_connects[0].cols(), 0);
  colors.resize(colors.size() + _connects[1].cols(), 1);

  return 0;
}





}


