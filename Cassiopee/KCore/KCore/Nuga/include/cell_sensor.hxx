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
class cell_sensor : public sensor<mesh_t, std::vector<typename mesh_t::cell_incr_t>>
{
  public:
    bool _single_pass_done;

    using cell_incr_t = typename mesh_t::cell_incr_t;
    using sensor_input_t = std::vector<cell_incr_t>; // either vector<int> (ISO) or vector(incr_t<3>) (DIR_PROTO)
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
    
    bool fill_adap_incr(output_t& adap_incr, bool do_agglo) override ;

    virtual bool stop() override { return _single_pass_done; }
};

/// 
template <typename mesh_t>
E_Int cell_sensor<mesh_t>::assign_data(const sensor_input_t& data)
{
  E_Int nphs = parent_t::_hmesh._ng.PHs.size();
  parent_t::_data.clear();
  parent_t::_data.resize(nphs, cell_incr_t(0)); // resize anyway to ensure same size as cells

  // now tranfer input data sized as enabled cells

  size_t pos{0};
  for (E_Int i = 0; i < nphs; ++i)
  {
    if (pos >= data.size()) break; //input was too small (should not happen)
    if (!parent_t::_hmesh._PHtree.is_enabled(i)) continue;

    parent_t::_data[i] = data[pos++];
  }

  _single_pass_done = false;//reinit to enable
  
  return 0;
}


///
template <typename mesh_t>
bool cell_sensor<mesh_t>::fill_adap_incr(output_t& adap_incr, bool do_agglo)
{
  using cell_incr_t = typename output_t::cell_incr_t;
  using face_incr_t = typename output_t::face_incr_t;
  // almost nothing to do, just pass the data as cell_adap_incr
  adap_incr.face_adap_incr.clear();
  E_Int nb_faces = parent_t::_hmesh._ng.PGs.size();
  adap_incr.face_adap_incr.resize(nb_faces, face_incr_t(0));

  //adap_incr.cell_adap_incr = parent_t::_data;
  adap_incr.cell_adap_incr.clear();
  adap_incr.cell_adap_incr.resize(parent_t::_data.size(), cell_incr_t(0));
  for (size_t k = 0; k < adap_incr.cell_adap_incr.size(); ++k)
    adap_incr.cell_adap_incr[k] = parent_t::_data[k];


  _single_pass_done = true; // to exit after first outer loop (based on sensor update)

  E_Int cmax = 0;
  for (size_t k = 0; k < adap_incr.cell_adap_incr.size(); ++k)
    cmax = std::max(cmax, adap_incr.cmax(k));

  E_Int cmin = IDX_NONE; // max int32
  for (size_t k = 0; k < adap_incr.cell_adap_incr.size(); ++k)
    cmin = std::min(cmin, adap_incr.cmin(k));

  //std::cout << "cell adapt incr min/max : " << cmin << "/" << cmax << std::endl;

  return (cmin != 0 || cmax != 0);
}

}
#endif
