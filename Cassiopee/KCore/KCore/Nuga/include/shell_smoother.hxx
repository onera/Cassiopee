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

#ifndef NUGA_SHELL_SMOOTHER_HXX
#define NUGA_SHELL_SMOOTHER_HXX

#include "smoother.hxx"

namespace NUGA
{

template <typename mesh_t>
struct shell_smoother : public smoother<mesh_t>
{
  using output_t = adap_incr_type<mesh_t::SUBTYPE>;
  
  Vector_t<E_Int> _Ln;

  shell_smoother() = default;

  void smooth(const mesh_t& hmesh, output_t& adap_incr);

  void __update_nodal_levels(const mesh_t& hmesh, E_Int PHi, E_Int aincr);

  ~shell_smoother() = default;
};

template <typename mesh_t>
void shell_smoother<mesh_t>::smooth(const mesh_t& hmesh, output_t& adap_incr)
{
  if (adap_incr.cell_adap_incr.size() == 0) return;

  E_Int n_nodes = hmesh._crd.cols();
  _Ln.resize(n_nodes, 0);

  for (size_t i = 0; i<adap_incr.cell_adap_incr.size(); i++)
  {
    if (adap_incr.cell_adap_incr[i] > 0) __update_nodal_levels(hmesh, i, adap_incr.cmax(i));// hack CLEF cell_adap_incr[i]);
  }

  // NODAL SMOOTHING LOOP
  E_Int Lc;
  E_Int n(0);
  bool carry_on = true;
  while (carry_on)
  {
    n = n + 1;
    carry_on = false;
    for (int i = 0; i<hmesh._ng.PHs.size(); i++)
    {  
      E_Int Mnodes(0);

      if (!hmesh._PHtree.is_enabled(i)) continue;

      Lc = hmesh._PHtree.get_level(i);
      E_Int Lcincr = Lc + adap_incr.cmax(i);// cell_adap_incr[i];

      const E_Int* pPHi = hmesh._ng.PHs.get_facets_ptr(i);

      // Mnodes definition
      for (int j = 0; j<hmesh._ng.PHs.stride(i); j++) {
        E_Int PGi = *(pPHi + j) - 1;
        E_Int PHn = NEIGHBOR(i, hmesh._F2E, PGi);
        // if the neighbor exist and its children are enabled then travers PGi children nodes
        bool traverse{false};
        if ((PHn != IDX_NONE) && !(hmesh._PHtree.is_enabled(PHn)))
        {
          // check that its children (if exist) are enabled
          E_Int nbc = hmesh._PHtree.nb_children(PHn);
          if (nbc)
          {
            // if any child is enabled, all of the are, so check on first one
            const E_Int* children = hmesh._PHtree.children(PHn);
            traverse = hmesh._PHtree.is_enabled(children[0]);
          }
        }

        if (traverse)
        {
          E_Int nb_ch = hmesh._PGtree.nb_children(PGi);
          assert (nb_ch != 0);
          const E_Int* enf = hmesh._PGtree.children(PGi);
          for (int k = 0; k< nb_ch; k++)
          {
            E_Int stride = hmesh._ng.PGs.stride(*(enf + k));
            const E_Int* ind = hmesh._ng.PGs.get_facets_ptr(*(enf + k));
            for (int l = 0; l<stride; l++)
              Mnodes = std::max(_Ln[*(ind + l) - 1], Mnodes);
          }
        }
        else // get the max among PGi nodes
        {
          E_Int stride = hmesh._ng.PGs.stride(PGi);
          const E_Int *ind = hmesh._ng.PGs.get_facets_ptr(PGi);
          for (int l = 0; l< stride; l++)
            Mnodes = std::max(_Ln[*(ind + l) - 1], Mnodes);
        }
      }
      // fin Mnodes definition
      if (Lcincr < Mnodes - 1)
      {
        const E_Int* faces = hmesh._ng.PHs.get_facets_ptr(i);
        E_Int nb_faces = hmesh._ng.PHs.stride(i);
        bool admissible_elt = K_MESH::Polyhedron<0>::is_basic(hmesh._ng.PGs, faces, nb_faces);

        if (!admissible_elt)
          continue;

        adap_incr.cell_adap_incr[i] = Mnodes - 1 - Lc;
        __update_nodal_levels(hmesh, i, adap_incr.cmax(i));// hack CLEF cell_adap_incr[i]);
        carry_on = true;
      }
    }
  }
}

template <typename mesh_t>
void shell_smoother<mesh_t>::__update_nodal_levels(const mesh_t& hmesh, E_Int PHi, E_Int aincr)
{
  E_Int Lc = hmesh._PHtree.get_level(PHi);
  E_Int nb_faces = hmesh._ng.PHs.stride(PHi);

  for (int j = 0; j< nb_faces; j++)
  {
    const E_Int* PGj = hmesh._ng.PHs.get_facets_ptr(PHi);
    E_Int PG = PGj[j] - 1;
    const E_Int* pN = hmesh._ng.PGs.get_facets_ptr(PG);
    E_Int nb_nodes = hmesh._ng.PGs.stride(PG);
    for (int l = 0; l< nb_nodes; l++)
      _Ln[*(pN + l) - 1] = std::max(_Ln[*(pN + l) - 1], Lc + aincr);
  }
}

}
#endif
