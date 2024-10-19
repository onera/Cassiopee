/*    
    Copyright 2013-2024 Onera.

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

#include "Nuga/include/IdTool.h"
#include <algorithm>
#include <map>
#include "Nuga/include/macros.h"

//#include <iostream>
namespace K_CONNECT
{
  /// bijective
  void IdTool::reverse_indirection
  (const Vector_t<E_Int> & assoc, Vector_t<E_Int>& reverse_assoc)
  {
    E_Int maxval = -1;
  
    for (size_t i = 0; i < assoc.size(); ++i)
    {
      if (assoc[i] != IDX_NONE)
        maxval = (maxval < assoc[i]) ? assoc[i] : maxval;
    }
  
    if (maxval == IDX_NONE)
      return;
  
    reverse_assoc.clear();
    reverse_assoc.resize(maxval+1, IDX_NONE);
  
    for (size_t i = 0; i < assoc.size(); ++i)
    {
      const E_Int & Ni = assoc[i];
      if (Ni != IDX_NONE)
        reverse_assoc[Ni]=i;
    }
  }

  ///
  void IdTool::reverse_indirection(const E_Int*oids, E_Int sz, std::map<E_Int, std::vector<E_Int>>& split_map)
  {
    split_map.clear();

#ifdef DEBUG_NGON_UNIT
    E_Int max_id = -1;
#endif

    for (E_Int i = 0; i < sz; ++i)
    {
      E_Int id = oids[i]; //+ shift;
      if (id == IDX_NONE) continue; //created entity
      split_map[id].push_back(i);

#ifdef DEBUG_NGON_UNIT
      max_id = (max_id < id) ? id : max_id;
#endif

    }

#ifdef DEBUG_NGON_UNIT
    assert(max_id < nb_pgs);
#endif
  }
  
  /// n-to-1 => 1-to-n
  void IdTool::reverse_indirection(E_Int nb_pgs, const E_Int*oids, E_Int sz, ngon_unit& split_graph)
  {
    std::map<E_Int, std::vector<E_Int> > molecules;

    reverse_indirection(oids, sz, molecules);

    E_Int empty = IDX_NONE;
    std::map<E_Int, std::vector<E_Int> >::const_iterator it;

    split_graph.clear();
    split_graph._NGON.resize(2, 0);

    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      it = molecules.find(i);
      if (it == molecules.end())
        split_graph.add(1, &empty); // has disappeared
      else
        split_graph.add(it->second.size(), &it->second[0]);
    }

    split_graph.updateFacets();
  }
  
  ///
  void IdTool::negative(Vector_t<bool>& flag)
  {
    for (size_t i = 0; i < flag.size(); ++i)flag[i]=!flag[i];
  }
  void IdTool::negative(Vector_t<E_Int>& flag)
  {
    for (size_t i = 0; i < flag.size(); ++i)flag[i] = -flag[i];
  }
  
  ///
  void IdTool::propagate(const Vector_t<E_Int>& nids, Vector_t<E_Int>& oids)
  {
    size_t sz = oids.size();
    for (size_t i = 0; i < sz; ++i)
    {
      E_Int& n = oids[i];
      if (n != IDX_NONE)
        n = nids[n];
    }
  }
  
  ///
  void IdTool::build_indir(const std::vector<bool>& keep, std::vector<E_Int> & nids)
  {
    size_t sz = keep.size();    
    bool carry_on(false);
    size_t i1(0), i2(sz-1);

    nids.clear();
    nids.resize(sz, IDX_NONE);
    
    size_t osz = std::count(keep.begin(), keep.end(), true);
    if (!osz) // keep nothing
      return;
        
    do{

      while ((i1 <= i2) && keep[i1]){nids[i1] = i1; ++i1;}  // Get the first empty column.
      while ((i1 <= i2) && !keep[i2]){--i2;} // Get the first column to move from the tail.

      carry_on = (i1 < i2);
      if (carry_on)
        nids[i2] = i1;
    }
    while (carry_on);
  }
  
  ///
  void IdTool::build_indir(const std::vector<bool>& keep, std::vector<E_Int> & nids, std::vector<E_Int> & oids)
  {
    size_t sz = keep.size();    
    bool carry_on(false);
    size_t i1(0), i2(sz-1);

    nids.clear();
    nids.resize(sz, IDX_NONE);
    
    oids.clear();
    size_t osz = std::count(keep.begin(), keep.end(), true);
    oids.resize(osz, IDX_NONE);
    
    if (!osz) // keep nothing
      return;
    
    do{

      while ((i1 <= i2) && keep[i1]){nids[i1] = i1; oids[i1]=i1; ++i1;}  // Get the first empty column.
      while ((i1 <= i2) && !keep[i2]){--i2;} // Get the first column to move from the tail.

      if ( (carry_on = (i1 < i2)) )
      {
        nids[i2] = i1;
        oids[i1] = i2;
        ++i1; --i2;
      }
    }
    while (carry_on);
  }
  
  ///
  E_Int IdTool::max(const K_FLD::IntArray& connect)
  {
    E_Int mid = -1, id;
    for (E_Int i = 0; i < connect.cols(); ++i)
    {
      for (E_Int j = 0; j < connect.rows(); ++j)
      {
        id = connect(j,i);
        if (id != IDX_NONE) mid=std::max(mid, id);
      }
    }
    return mid;
  }
  
  ///
  E_Int IdTool::max(const Vector_t<E_Int>& vec)
  {
    E_Int mid = -1;
    for (size_t i=0; i < vec.size(); ++i)
    {
      const E_Int& id = vec[i];
      if (id != IDX_NONE)
          mid=std::max(mid, id);
    }
    return mid;
  }
  
  ///
  E_Int IdTool::min(const K_FLD::IntArray& connect)
  {
    E_Int mid = IDX_NONE;
    for (E_Int i = 0; i < connect.cols(); ++i)
    {
      for (E_Int j = 0; j < connect.rows(); ++j)
        mid=std::min(mid, connect(j,i));
    }
    return mid;
  }
  
  ///
  E_Int IdTool::min(const Vector_t<E_Int>& vec)
  {
    E_Int mid = IDX_NONE;
    for (size_t i=0; i < vec.size(); ++i)
      mid=std::min(mid, vec[i]);
    return mid;
  }
  
  ///
  void IdTool::compact(std::vector<E_Int>& vec, const std::vector<bool> & keep)
  {
    std::vector<E_Int> tmp;
    //tmp.reserve(vec.size());
    
    for (size_t i = 0; i < vec.size(); ++i)
    {
      if (keep[i])
        tmp.push_back(vec[i]);
    }
    
    vec=tmp;
  }
  
  ///
  void IdTool::reverse_sorting(Vector_t<E_Int> & vec)
  {
    std::reverse(vec.begin(), vec.end());
  }
  
  ///
  bool IdTool::equal(const E_Int* p, const E_Int* q, E_Int n, bool permut_accepted, bool strict_orient)
  {
    
    /*for (size_t i=0; i < n; ++i)
    {
      std::cout << p[i] << " vs " << q[i] << std::endl;
    }*/
    
    E_Int start=0;
    if (permut_accepted) //find the rank in q of p[0]
    {
      for (; start < n; ++start)
        if (p[0] == q[start])break;
      if (start == n)
        return false;
    }
   
    E_Int r=1;
    if (!strict_orient)
    {
      if (p[1] == q[((start-1)%n + n)%n])
        r=-1; //go reverse
    }
   
    bool equal=true;
    for (E_Int i=0; (i<n) && equal; ++i)
      equal &= (p[i] == q[((start+r*i)%n + n)%n]);
   
    return equal;
 }

  void IdTool::right_shift(E_Int* list, E_Int sz, int shift)
  {
    if (shift == 0) return;

    STACK_ARRAY(E_Int, sz, tmp);

    for (int i = 0; i < sz; ++i)
      tmp[i] = list[(i + shift) % sz];

    for (int i = 0; i < sz; ++i)
      list[i] = tmp[i];
  }
}
