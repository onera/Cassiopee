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

#ifndef NUGA_BALANCING_SENSOR_HXX
#define NUGA_BALANCING_SENSOR_HXX

#include "Nuga/include/sensor.hxx"
#include "Nuga/include/V1_smoother.hxx"
#include "Nuga/include/join_plan.hxx"


namespace NUGA
{

//
template <typename mesh_t>
class join_sensor : public sensor<mesh_t, std::map<E_Int, K_FLD::IntArray>> 
{
  public:
    using pg_arr_t = typename mesh_t::pg_arr_t;
    using input_t = std::map<E_Int, K_FLD::IntArray>;
    using parent_t = sensor<mesh_t, input_t>;
    using output_t = typename mesh_t::output_t;

    join_sensor(mesh_t& mesh) : parent_t(mesh, new V1_smoother<mesh_t>()) {}

    bool fill_adap_incr(output_t& adap_incr, bool do_agglo) override;
    bool update() override;
    bool stop() override { return parent_t::_data.empty(); }
};

/// WARNING : only for ISO currently
template <typename mesh_t>
bool join_sensor<mesh_t>::fill_adap_incr(output_t& adap_incr, bool do_agglo)
{
  bool filled{ false };

  adap_incr.face_adap_incr.clear();
  adap_incr.cell_adap_incr.clear();

  //
  E_Int nb_elts = parent_t::_hmesh._ng.PHs.size();
  E_Int nb_faces = parent_t::_hmesh._ng.PGs.size();

  using cell_incr_t = typename output_t::cell_incr_t;
  using face_incr_t = typename output_t::face_incr_t;

  adap_incr.cell_adap_incr.resize(nb_elts, cell_incr_t(0));
  adap_incr.face_adap_incr.resize(nb_faces, face_incr_t(0));

  // 1. filling adap_incr.face_adap_incr : select if face plan ask for and not already subdivided
  for (auto& dat : parent_t::_data)
  {
    E_Int PGi = dat.first;
    auto& plan = dat.second;
    if (plan.getSize() == 0) continue;                            // nothing planned
    //std::cout << plan << std::endl;

    if (parent_t::_hmesh._PGtree.nb_children(PGi) != 0) continue; // already subdivided

    adap_incr.face_adap_incr[PGi] = 1;
    filled = true;
  }

  // 2. filling adap_incr.cell_adap_incr : select if face plan is going beyond next generation
  for (int i = 0; i < nb_elts; ++i)
  {
    const E_Int* faces = parent_t::_hmesh._ng.PHs.get_facets_ptr(i);
    E_Int nb_faces = parent_t::_hmesh._ng.PHs.stride(i);
    bool admissible_elt = K_MESH::Polyhedron<0>::is_basic(parent_t::_hmesh._ng.PGs, faces, nb_faces);//fixme : not true for ISO_HEX
    if (!admissible_elt) continue;

    bool require = false;
    for (E_Int j = 0; (j < nb_faces) && !require; ++j)
    {
      E_Int PGi = faces[j] - 1;
      auto idat = parent_t::_data.find(PGi);
      if (idat == parent_t::_data.end()) continue;

      auto& plan = idat->second;
      require = join_plan<pg_arr_t>::one_child_requires(plan);
    }

    if (require) {
      adap_incr.cell_adap_incr[i] = 1;
      filled = true;
    }
  }
  return filled;
}

///
template <typename mesh_t>
bool join_sensor<mesh_t>::update()
{
  // compute plans for next generation and enable joins PG for that generation
  input_t new_data;

  //std::cout << "SENSOR ENABLING UPDATE" << std::endl;
  
  for (auto& dat : parent_t::_data)
  {
    E_Int PGi = dat.first;
    auto& plan = dat.second;

    if (plan.getSize() == 0) continue;                            // nothing planned
    //std::cout << plan << std::endl;
  
    // FIXME : ENABLING FACES IS NOT WORKING PROPERLY! (see also hierarchical_mesh::enable_PGs)
    // the OLD MODE approach (maintained only for DIR mode now) is
    // to enable at hierarchical_mesh::enable_PGs stage any natural face of enabled cells
    // so at the end of this stage some join face can be wrongly disbaled
    // (those which ensure conformal join with the other side)
    // indeed, for some join, children must be enabled instead of natural face,
    // these are SUPPOSED TO BE fixed here at this stage
    // but some cases with a complex enabling state of the leaves is not handled.
    // The correct behaviour should be to enable on each side the merged highest enabled leaves
    // meanwhile, a NEW MODE is introduced : this stage handles entirely the enabling (nothing done in join_sensor)
    // by blindly and systematicaly enabling all the leaves. Since each side have the same face hierarchies
    // this should work, by providing an overdefined state (more subdivisions than required).
    // NEED TO BE DOUBLE CHECKED though when enabling agglomeration.
    bool OLD_MODE = (mesh_t::SUBTYPE == DIR);

    if (OLD_MODE)
    {
      if (parent_t::_hmesh._PGtree.is_enabled(PGi))
      {
        // enable its children instead : if they are in plan, it means they are enabled on the other join side
        E_Int nbc{ parent_t::_hmesh._PGtree.nb_children(PGi) };
        if (nbc != 0)
        {
          const E_Int* children = parent_t::_hmesh._PGtree.children(PGi);
          for (E_Int i = 0; i < nbc; ++i)
            parent_t::_hmesh._PGtree.enable(children[i]);
        }
      }
    }

    // get plans for children for the next adaptation pass
    join_plan<pg_arr_t>::extract_sub_plans(parent_t::_hmesh._PGtree, PGi, plan, new_data);
  }

  parent_t::_data.clear();
  if (!new_data.empty())
  {
    parent_t::_data = new_data;
    return true;
  }
  return false;
}

}
#endif
