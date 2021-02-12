/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_BALANCING_SENSOR_HXX
#define NUGA_BALANCING_SENSOR_HXX

#include "Nuga/include/sensor.hxx"
#include "Nuga/include/V1_smoother.hxx"
#include "Nuga/include/join_plan.hxx"


namespace NUGA
{

//
template <typename mesh_t>
class join_sensor : public sensor<mesh_t, std::map<E_Int, typename mesh_t::pg_arr_t>> // Vector_t might be templatized to have ISO mode with a metric field
{
  public:
    using pg_arr_t = typename mesh_t::pg_arr_t;
    using input_t = std::map<E_Int, pg_arr_t>;
    using parent_t = sensor<mesh_t, input_t>;
    using output_t = typename mesh_t::output_t; //fixme: static assert to add : must be ISO => IntVec

    join_sensor(mesh_t& mesh) : parent_t(mesh, new V1_smoother<mesh_t>()) {}

    void fill_adap_incr(output_t& adap_incr, bool do_agglo) override;
    bool update() override;
    bool stop() override { return parent_t::_data.empty(); }
};

/// WARNING : only for ISO currently
template <typename mesh_t>
void join_sensor<mesh_t>::fill_adap_incr(output_t& adap_incr, bool do_agglo)
{
  adap_incr.face_adap_incr.clear();
  adap_incr.cell_adap_incr.clear();

  //
  E_Int nb_elts = parent_t::_hmesh._ng.PHs.size();
  E_Int nb_faces = parent_t::_hmesh._ng.PGs.size();

  adap_incr.cell_adap_incr.resize(nb_elts, 0);
  adap_incr.face_adap_incr.resize(nb_faces, 0);

  // 1. filling adap_incr.face_adap_incr : select if face plan ask for and not already subdivided
  for (auto& dat : parent_t::_data)
  {
    E_Int PGi = dat.first;
    auto& plan = dat.second;
    if (plan.getSize() == 0) continue;                            // nothing planned
    //std::cout << plan << std::endl;

    if (parent_t::_hmesh._PGtree.nb_children(PGi) != 0) continue; // already subdivided

    adap_incr.face_adap_incr[PGi] = 1;
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

    if (require) adap_incr.cell_adap_incr[i] = 1;
  }
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
    
    // enable children : if they are in plan, it means they are enabled on the other join side
    E_Int nbc{parent_t::_hmesh._PGtree.nb_children(PGi)};
    if (nbc != 0)
    {
      const E_Int* children = parent_t::_hmesh._PGtree.children(PGi);
      for (E_Int i=0; i < nbc; ++i)
        parent_t::_hmesh._PGtree.enable(children[i]);
    }

    // get plans for children fort next adaptation pass
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
