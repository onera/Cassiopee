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

#pragma once
#include "Fld/DynArray.h"
#include "Search/BbTree.h"
#include <vector>

struct XPair{

  XPair(E_Int s0, E_Int c0, E_Int s1, E_Int c1):S0(s0), C0(c0), S1(s1), C1(c1){}

  E_Int S0;
  E_Int C0;
  E_Int S1;
  E_Int C1;
};

class Intersector
{
  typedef K_SEARCH::BBox3D BBoxType;
  typedef K_SEARCH::BbTree3D TreeType;

public:
  static void getXPairs(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray* >& Cs,
                 std::vector<XPair>& pairs);
private:

  ///
  static void __getXPairs(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& Cs,
                          E_Int s1, E_Int s2,
                          const std::vector<TreeType* >& trees, 
                          const std::vector< std::vector<BBoxType*> >& boxes,
                          const std::vector<std::vector<E_Float> >& Lengths,
                          E_Float tol, bool tol_is_absolute,
                          std::vector<XPair>& pairs);
  
  ///
  static void __create_boxes (const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt, 
                              std::vector<K_SEARCH::BBox3D*>& boxes, K_SEARCH::BBox3D*& pool, K_SEARCH::BBox3D& GBbox, E_Float enlarge_factor=0.);

  ///
  inline static bool __overlap(const K_FLD::FloatArray& pos, K_FLD::IntArray::const_iterator pS1,
                               K_FLD::IntArray::const_iterator pS2,
                               E_Float tolerance, E_Bool tol_is_absolute, E_Float Lc);

  static bool __pointIsInside(const K_FLD::FloatArray& pos, K_FLD::IntArray::const_iterator pS1,
                              K_FLD::IntArray::const_iterator pS2, E_Float tolerance, E_Bool tol_is_absolute,
                              E_Float Lc);

  static bool __edgeOverlaps(const K_FLD::FloatArray& pos, K_FLD::IntArray::const_iterator pS1,
                             K_FLD::IntArray::const_iterator pS2, E_Float tolerance, E_Bool tol_is_absolute,
                             E_Float Lc);


  Intersector(void){}
  ~Intersector(void){}
};
