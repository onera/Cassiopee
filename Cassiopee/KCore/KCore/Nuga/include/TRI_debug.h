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

#ifndef __TRI_DEBUG_H__
#define	__TRI_DEBUG_H__

# include "Nuga/include/DynArray.h"
#include <sstream>


#define COL_DEFAULT 0
#define COL_RED     1
#define COL_GREEN   2
#define COL_NEXT    5
#define COL_NEXT2   7

class TRI_debug
{
//    
public: 
  
  static void draw_connected_to_node_T3s(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, E_Int Ni);
  //static void draw_connected_to_edge_T3s(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, E_Int Ni);
  
  static void get_connected_to_T3(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, E_Int Ti, std::vector<E_Int> & cIds);
  static void draw_connected_to_T3(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, E_Int Ti);
  
  static void draw_connected_to_T3(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, E_Int Ni, E_Int Nj, E_Int Nk/*, bool shrink*/);
  static void connected_to_T3(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, E_Int Ni, E_Int Nj, E_Int Nk/*, bool shrink*/,
                              K_FLD::FloatArray& coordOut, K_FLD::IntArray& connOut, std::vector<E_Int>& colors);
  static void draw_connected_to_E2(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, E_Int Ni, E_Int Nj);
  
  static void remove_T3 (const K_FLD::FloatArray& coord, K_FLD::IntArray& connectT3, E_Int Ti);
  static void remove_T3 (const K_FLD::FloatArray& coord, K_FLD::IntArray& connectT3, E_Int Ni, E_Int Nj, E_Int Nk);
  static E_Int get_T3_index(const K_FLD::IntArray& connect, E_Int Ni, E_Int Nj, E_Int Nk);
  
  static void get_T3_neighbors(const char* fname, E_Int Ti, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, K_FLD::IntArray& neighbors, bool both_orient=false);
  
  static void coloring_frames(const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, K_FLD::IntArray& neighbors, E_Int Kseed/*one of the target faulty PHT3s triangles*/);
  
  static void get_same_ancestor_T3s (E_Int PGi, const std::vector<E_Int>& ancestors, std::vector<E_Int>& oids);
  static void draw_same_ancestor_T3s(E_Int PGi, E_Int iter, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, const std::vector<E_Int>& ancestors);
  
  static void write_wired(const char* fname, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, bool normal, const std::vector<E_Int>*T3colors=0, const std::vector<bool>* mask=0);
  static void write_wired(const char* fname, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connectT3, const K_FLD::FloatArray& normals, const std::vector<E_Int>* T3colors=0, const std::vector<bool>* keep=0, bool localid=false);

  /// Gives info for debugging a faulty PHT3
  static bool analyze_T3_set(E_Int setid, const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt, const std::vector<E_Int>& T3s, bool& is_manifold, bool& is_closed, std::vector<E_Int>* colors=0);
};

#endif