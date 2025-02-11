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

#ifndef NUGA_V1_SMOOTHER_HXX
#define NUGA_V1_SMOOTHER_HXX

#include "smoother.hxx"


namespace NUGA
{

template <typename mesh_t>
struct V1_smoother : public smoother<mesh_t>
{
  using output_t = adap_incr_type<mesh_t::SUBTYPE>;
  
  V1_smoother() = default;

  virtual void smooth(const mesh_t& hmesh, output_t& adap_incr) ;

  ~V1_smoother(){}
};

template <typename mesh_t>
void V1_smoother<mesh_t>::smooth(const mesh_t& hmesh, output_t& adap_incr)  {

  const auto& ng = hmesh._ng;
  const auto& PHtree = hmesh._PHtree;

  if (adap_incr.cell_adap_incr.size() == 0) return;
  
  std::stack<E_Int> stck;

  for (E_Int i = 0; i < ng.PHs.size(); i++)
    if (adap_incr.cell_adap_incr[i] != 0) stck.push(i);

  std::vector<E_Int> neighbours;

  while (!stck.empty()) {

    E_Int ind_PHi = stck.top(); // index of ith PH
    stck.pop();

    //fixme : dangerous : operator[] returns a ref
    // so doing :
    //     auto incr = adap_incr.cell_adap_incr[ind_PHi] + PHtree.get_level(ind_PHi);
    // updates adap_incr.cell_adap_incr[ind_PHi] !!
    // so do it now in 2 steps (next 2 lines) but other problems could occur elsewhere
    auto incr = adap_incr.cell_adap_incr[ind_PHi];
    incr += PHtree.get_level(ind_PHi);

    E_Int s = ng.PHs.stride(ind_PHi);
    neighbours.clear();
    neighbours.resize(4*s);//fixme 4s
    E_Int nb_neighbours = 0;
    
    hmesh.get_enabled_neighbours(ind_PHi, neighbours.begin(), nb_neighbours); // get the effective neighbours (enabled ones) to assert the 2:1 rule

    for (int i = 0; i < nb_neighbours; ++i)
    {
      //fixme : dangerous : operator[] returns a ref
      // so doing :
      //    auto incr_neigh = adap_incr.cell_adap_incr[neighbours[i]] + PHtree.get_level(neighbours[i])
      // updates adap_incr.cell_adap_incr[neighbours[i]] !!
      // so do it now in 2 steps (next 2 lines) but other problems could occur elsewhere
      auto incr_neigh = adap_incr.cell_adap_incr[neighbours[i]];
      incr_neigh += PHtree.get_level(neighbours[i]);

      if (NUGA::abs(incr - incr_neigh) <= 1) // 2:1 rule respected
        continue;

      // not respected : increase by 1 the adap incr of the lowest incr (ind_PHi or neighbours[i])
      E_Int PH_to_mod = (incr > incr_neigh) ? neighbours[i] : ind_PHi;

      if (adap_incr.cell_adap_incr[PH_to_mod] >= 0) // increase by 1 to refine : +=1 the adap incr and add the cell in the stack
      {
        adap_incr.cell_adap_incr[PH_to_mod] += 1;
        stck.push(PH_to_mod);
      }
      else // increase by 1 an agglomeration : +=1 the adap incr for the cell and its siblings, and add all of them in the stack
      {
        E_Int father = PHtree.parent(PH_to_mod);
        const E_Int* p = PHtree.children(father);
        E_Int nb_children = PHtree.nb_children(father);
        for (int j = 0; j < nb_children; ++j)
        {
          adap_incr.cell_adap_incr[p[j]] += 1;
          stck.push(p[j]);
        }
      }
    }
  }
}

}
#endif
