/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_BALANCING_SENSOR_HXX
#define NUGA_BALANCING_SENSOR_HXX

#include "Nuga/include/sensor.hxx"
#include "Nuga/include/V1_smoother.hxx"


namespace NUGA
{

  using input_type = std::map<E_Int, K_FLD::IntArray>;

template <typename mesh_t>
class join_sensor : public sensor<mesh_t, input_type> // Vector_t might be templatized to have ISO mode with a metric field
{
  public:
    using input_t = input_type;
    using parent_t = sensor<mesh_t, input_t>;
    using output_t = typename mesh_t::output_t; //fixme: static assert to add : must be ISO => IntVec

    join_sensor(mesh_t& mesh) : parent_t(mesh, new V1_smoother<mesh_t>()) {}

    void fill_adap_incr(output_t& adap_incr, bool do_agglo) override;
    bool update() override;
    bool stop() override { return parent_t::_data.empty(); }
  
  private:
    ///
    void extract_subplan(E_Int i, const K_FLD::IntArray& plan, K_FLD::IntArray& sub_plan);
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
  bool has_face_adap{ false };
  for (auto& dat : parent_t::_data)
  {
    E_Int PGi = dat.first;
    K_FLD::IntArray& plan = dat.second;
    if (plan.cols() == 0) continue;                               // nothing planned
    //std::cout << plan << std::endl;

    if (parent_t::_hmesh._PGtree.nb_children(PGi) != 0) continue; // already subdivided

    adap_incr.face_adap_incr[PGi] = 1;
    has_face_adap = true;
  }

  if (!has_face_adap) return;

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

      K_FLD::IntArray& plan = idat->second;
      if (plan.cols() == 0) continue; // nothing planned
      
      for (size_t k = 0; (k < plan.rows()) && !require; ++k)
      {
        const E_Int& plan_for_childk = plan(k, 0);
        //require if grand children are planned : at least one of PGi's children is planned for having children
        require |= (plan_for_childk != NO_CHIDREN);
      }
    }

    if (require) adap_incr.cell_adap_incr[i] = 1;
  }
}

///
template <typename mesh_t>
bool join_sensor<mesh_t>::update()
{
  // compute plans for next generation
  input_t new_data;
  K_FLD::IntArray sub_plan;

  for (auto& dat : parent_t::_data)
  {
    E_Int PGi = dat.first;
    const E_Int* children = parent_t::_hmesh._PGtree.children(PGi);
    E_Int nbchildren = parent_t::_hmesh._PGtree.nb_children(PGi);

    assert(nbchildren != 0); //must have been subdivided or had children already
    assert(nbchildren == 4);           // WARNING : ISO ASSUMPTION  : NBC = 4

    K_FLD::IntArray& plan = dat.second;
    if (plan.cols() == 0) // no plan
      continue;
    if (plan.rows() == 1) // no grand children
      continue;

    // at least one column (sized as NBC)
    E_Int nnodes = plan.cols();
    assert(nnodes >= 1);
    assert(plan.rows() == nbchildren);

    for (size_t n= 0; n <plan.rows(); ++n)
    {
      extract_subplan(plan(n,0), plan, sub_plan);
      if (sub_plan.cols() != 0) new_data[children[n]] = sub_plan;
    }
  }

  parent_t::_data.clear();
  if (!new_data.empty())
  {
    parent_t::_data = new_data;
    return true;
  }
  return false;
}

///
template <typename mesh_t>
void join_sensor<mesh_t>::extract_subplan(E_Int i, const K_FLD::IntArray& plan, K_FLD::IntArray& sub_plan)
{
  //
  sub_plan.clear();

  if (i == NO_CHIDREN) return;

  if (i == NO_GRAND_CHILDREN)
  {
    sub_plan.resize(plan.rows(), 1, NO_CHIDREN);
    return;
  }

  std::vector<E_Int> nodepool;
  nodepool.push_back(i);

  std::vector<bool> keep(plan.cols(), false);
  keep[i] = true;
  
  while (!nodepool.empty())
  {
    E_Int c = nodepool.back();
    nodepool.pop_back();

    for (E_Int j = 0; j < plan.rows(); ++j)
    {
      E_Int v = plan(j, c);
      if (v <= 0) continue;
      nodepool.push_back(v);
      keep[v] = true;
    }
  }
  //std::cout << plan << std::endl;
  sub_plan = plan;
  std::vector<E_Int> nids;
  K_CONNECT::keep<bool> pred_keep(keep);
  K_CONNECT::IdTool::compress(sub_plan, pred_keep, nids);
  //std::cout << sub_plan << std::endl;
  // now update "pointers" (column references)
  for (E_Int i = 0; i < sub_plan.cols(); ++i)
  {
    for (E_Int j = 0; j < sub_plan.rows(); ++j)
    {
      E_Int& p = sub_plan(j, i);
      if (p <= 0) continue;
      assert(p > 1); // first column is gone and is the only one to point on 2nd so must not encounter p==1
      p = nids[p];
      assert(p > 0); // must at least point on the second column
    }
  }
  //std::cout << sub_plan << std::endl;
}

}
#endif
