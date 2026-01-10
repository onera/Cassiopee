/*    
    Copyright 2013-2026 ONERA.

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

#ifndef __BOOLEANOPERATOR_H__
#define	__BOOLEANOPERATOR_H__

#include "Conformizer.h"
#include "Nuga/include/DynArray.h"
#include <vector>

namespace NUGA
{
  
class BooleanOperator
{ 

public:
  /// Conformized union of the input surfaces.
  E_Int getSum(K_FLD::FloatArray& coord, K_FLD::IntArray& connect, std::vector<E_Int>& colors);
  ///S1 -  Intersection(S1, S2).
  virtual E_Int get_1_minus_2(K_FLD::FloatArray& coord, K_FLD::IntArray& connect, std::vector<E_Int>& colors) = 0;
  ///S2 -  Intersection(S1, S2).
  virtual E_Int get_2_minus_1(K_FLD::FloatArray& coord, K_FLD::IntArray& connect, std::vector<E_Int>& colors) = 0;
  ///  Sum(S1, S2) - Intersection(S1, S2).
  virtual E_Int getUnion(K_FLD::FloatArray& coord, K_FLD::IntArray& connect, std::vector<E_Int>& colors) = 0;
  /// Intersection between S1 & S2.
  virtual E_Int getIntersection(K_FLD::FloatArray& coord, K_FLD::IntArray& connect, std::vector<E_Int>& colors) = 0;
  /// Intersection border (BARs for TRI, POINTs for BAR) 
  virtual E_Int getIntersectionBorder(K_FLD::FloatArray& coord, K_FLD::IntArray& connect) = 0;
  /// Destructor.
  virtual ~BooleanOperator(void){ delete _conformizer;}
  
protected:
  /// Constructor with the 2 input surfaces S1 & S2.
  BooleanOperator(const K_FLD::FloatArray& pos1, const K_FLD::IntArray& connect1,
                  const K_FLD::FloatArray& pos2, const K_FLD::IntArray& connect2,
                  E_Float tolerance, ConformizerRoot* c, int itermax=10);
  ///
  E_Int initialized();
  
private:
  ///
  virtual E_Int check_sanity() = 0;
  ///
  virtual E_Int compute_zones() = 0;
  
protected:
  ///
  ConformizerRoot*   _conformizer;
  ///
  K_FLD::FloatArray   _coord;
  ///
  K_FLD::IntArray     _connects[2];
  ///
  std::vector<E_Int>  _colors[2];
  ///
  bool                _initialized;
  ///
  K_FLD::IntArray     _connect_1_out_2;
  ///
  K_FLD::IntArray     _connect_2_out_1;
};

}

#endif	/* __BOOLEANOPERATOR_H__ */

