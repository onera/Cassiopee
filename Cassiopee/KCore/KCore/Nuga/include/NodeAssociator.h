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

#pragma once
#include "Nuga/include/DynArray.h"
#include "Intersector.h"
#include <vector>

class NodeAssociator
{
public:
  ///
  NodeAssociator(void);

  ///
  virtual ~NodeAssociator(void);

  ///
  virtual void make_pairs(const K_FLD::FloatArray& pos, const std::vector< K_FLD::IntArray* >& components,
                          std::vector<E_Int> &mates, K_FLD::IntArray& OneSurface);

protected:
  ///
  void __removeOverlaps(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& components,
                        std::vector<E_Int>& nmates, std::vector<K_FLD::IntArray>& componentsOut);

  ///
  void __make_pairs(const K_FLD::FloatArray& pos, const std::vector< std::vector<E_Int> >& sorted_nodes,
                            const std::vector<E_Int>& surface_colors, std::vector<E_Int> &pairs);

  ///
  static void __reorientComponent(const K_FLD::FloatArray& pos, K_FLD::IntArray& component, const std::vector<E_Int>& nmates); 

private:
  ///
  static void __priorize_X_zones(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& components, const std::vector<XPair>& pairs, 
                                 std::vector<std::vector<E_Int> >& componentsKeep, std::vector<std::vector<E_Int> >& componentsRemove,
                                 std::vector<std::vector<E_Int> >& componentsUnchanged);

  ///
  static void __priorize_X_zones2(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& components, const std::vector<XPair>& pairs, 
                                  std::vector<std::vector<E_Int> >& componentsKeep, std::vector<std::vector<E_Int> >& componentsRemove);

  ///
  static void __getRelativeOrient(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray>& surfaces, const std::vector<E_Int>& nmates,
                                  std::map<K_MESH::NO_Edge, E_Int>& surfpair_to_orient); 

  ///
  static void __getRelativeOrient(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray>& comps, const std::vector<XPair>& pairs,
                                  std::map<K_MESH::NO_Edge, E_Int>& comppair_to_orient);

  ///
  static void __reorientComponents(const K_FLD::FloatArray& pos, std::vector<K_FLD::IntArray>& component, const std::vector<XPair>& pairs);

public:
  static bool reorient;

};


