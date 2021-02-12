/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_CELL_SENSOR_HXX
#define NUGA_CELL_SENSOR_HXX

#include "Nuga/include/sensor.hxx"
#include "Nuga/include/shell_smoother.hxx"
#include "Nuga/include/V1_smoother.hxx"


namespace NUGA
{

template <typename mesh_t>
class cell_sensor : public sensor<mesh_t, Vector_t<E_Int>>
{
  public:
    bool _single_pass_done;

    using sensor_input_t = Vector_t<E_Int>;
    using parent_t = sensor<mesh_t, sensor_input_t>;
    using output_t = typename mesh_t::output_t; //fixme: static assert to add : must be ISO => IntVec

    cell_sensor(mesh_t& mesh, eSmoother smoo_type): parent_t(mesh, nullptr), _single_pass_done(false)
    {
      if (smoo_type == eSmoother::V1_NEIGH)
        parent_t::_smoother = new V1_smoother<mesh_t>();
      else if (smoo_type == eSmoother::SHELL)
        parent_t::_smoother = new shell_smoother<mesh_t>();
    }

    virtual E_Int assign_data(const sensor_input_t& data) override;
    
    void fill_adap_incr(output_t& adap_incr, bool do_agglo) override ;

    virtual bool stop() override { return _single_pass_done; }
};

/// 
template <typename mesh_t>
E_Int cell_sensor<mesh_t>::assign_data(const sensor_input_t& data)
{
  E_Int nphs = parent_t::_hmesh._ng.PHs.size();
  parent_t::_data.clear();
  parent_t::_data.resize(nphs, 0.); // resize anyway to ensure same size as cells

  // now tranfer input data sized as enabled cells

  size_t pos{0};
  for (E_Int i = 0; i < nphs; ++i)
  {
    if (pos >= data.size()) break; //input was too small (shoudl not happen)
    if (!parent_t::_hmesh._PHtree.is_enabled(i)) continue;

    parent_t::_data[i] = data[pos++];
  }

  return 0;
}


///
template <typename mesh_t>
void cell_sensor<mesh_t>::fill_adap_incr(output_t& adap_incr, bool do_agglo)
{
  // almost nothing to do, just pass the data as cell_adap_incr
  adap_incr.face_adap_incr.clear();
  E_Int nb_faces = parent_t::_hmesh._ng.PGs.size();
  adap_incr.face_adap_incr.resize(nb_faces, 0);

  adap_incr.cell_adap_incr = parent_t::_data;

  _single_pass_done = true; // to exit after first outer loop (based on sensor update)

  //std::cout << "cell adapt incr sz : " << adap_incr.cell_adap_incr.size() << std::endl;
  //E_Int minv = *std::min_element(ALL(adap_incr.cell_adap_incr));
  //E_Int maxv = *std::max_element(ALL(adap_incr.cell_adap_incr));
  //std::cout << "cell adapt incr min/max : " << minv << "/" << maxv << std::endl;
}

}
#endif
