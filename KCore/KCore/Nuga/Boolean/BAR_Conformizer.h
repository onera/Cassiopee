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

#ifndef __BAR_CONFORMIZER_H__
#define	__BAR_CONFORMIZER_H__

#include "Conformizer.h"

#include "Def/DefTypes.h"
#include <vector>

namespace NUGA
{

template <int DIM>
class BAR_Conformizer : public Conformizer<DIM, K_MESH::Edge> {

public :
  typedef Conformizer<DIM, K_MESH::Edge>  parent_type;
  
public:
  BAR_Conformizer();
  virtual ~BAR_Conformizer(){}
  
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
                        K_CONT_DEF::bool_vector_type& xc,
                        K_CONT_DEF::int_vector_type& ancestors);

#ifdef DEBUG_CONFORMIZER
  ///
  void drawElements(const char* fname, const char* filefmt, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connect, const std::vector<E2> & elts, bool localid = false, std::vector<E_Int>* colors = 0){/*todo*/ }
#endif
    
private:
  
  void __reorder_nodes_on_edge(const K_FLD::FloatArray& pos, std::vector<E_Int>& nodes, const E_Float *P0);
  
  BAR_Conformizer(const BAR_Conformizer& orig){}
  
private:

};
}

#include "BAR_Conformizer.cxx"

#endif	/* BAR_CONFORMIZER_H */

