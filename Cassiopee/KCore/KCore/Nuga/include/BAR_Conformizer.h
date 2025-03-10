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

#ifndef __BAR_CONFORMIZER_H__
#define	__BAR_CONFORMIZER_H__

#include "Conformizer.h"

#include "Nuga/include/defs.h"
#include <vector>

namespace NUGA
{

template <short DIM>
class BAR_Conformizer : public Conformizer<DIM, K_MESH::Edge> {

public :
  typedef Conformizer<DIM, K_MESH::Edge>  parent_type;
  
public:
  BAR_Conformizer(bool whisto = false) :parent_type(whisto) {}

  virtual ~BAR_Conformizer(){}
  
  std::vector<std::pair<E_Int, E_Int>> get_x_history();

  // Overridden Methods ///////////////////////////////////////////////////////

protected:
  
  ///
  void __set_tolerances(E_Float Lmin, E_Float Lmax, E_Float  user_tolerance);

  ///
  E_Int __intersect(K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
                     E2& e1, E2& e2, E_Float tol);
  ///
  void __update_data(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const std::vector<E_Int>& newIDs);
  /// Splits the triangles by triangulation.
  E_Int __split_Elements(const K_FLD::FloatArray& pos, K_FLD::IntArray & connect,
                        NUGA::bool_vector_type& xc,
                        NUGA::int_vector_type& ancestors);

#ifdef DEBUG_CONFORMIZER
  ///
  void drawElements(const char* fname, const char* filefmt, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connect, const std::vector<E2> & elts, bool localid = false, std::vector<E_Int>* colors = 0)
  {
    K_FLD::IntArray cB;
    for (size_t i = 0; i < elts.size(); ++i)
    {
      cB.pushBack(connect.col(elts[i].id), connect.col(elts[i].id)+2);
    }
    K_FLD::FloatArray crd = coord;
    crd.resize(3, coord.cols(), 0.);
    medith::write(fname, crd, cB, "BAR");

  }
#endif
    
private:
  
  void __reorder_nodes_on_edge(const K_FLD::FloatArray& pos, std::vector<E_Int>& nodes, const E_Float *P0);
  
  BAR_Conformizer(const BAR_Conformizer& orig){}
  
private:

};
}

#include "BAR_Conformizer.cxx"

#endif	/* BAR_CONFORMIZER_H */

