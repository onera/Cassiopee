/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr); Imad Hammani (imad.hammani@onera.fr)

#ifndef NUGA_METRIC_SENSOR_HXX
#define NUGA_METRIC_SENSOR_HXX

#include "Nuga/include/sensor.hxx"
#include "Nuga/include/Metric.h"
//#include "Nuga/include/shell_smoother.hxx"
//#include "Nuga/include/V1_smoother.hxx"


namespace NUGA
{

template <typename mesh_t>
class metric_sensor : public sensor<mesh_t, DELAUNAY::VarMetric<DELAUNAY::Aniso3D>>
{
  public:
    bool _single_pass_done;
    
    //using cell_incr_t = DELAUNAY::Aniso3D;
    using sensor_input_t = DELAUNAY::VarMetric<DELAUNAY::Aniso3D>; // either vector<int> (ISO) or vector(incr_t<3>) (DIR)
    using parent_t = sensor<mesh_t, sensor_input_t>;
    using output_t = typename mesh_t::output_t; //fixme: static assert to add : must be ISO => IntVec

    metric_sensor(mesh_t& mesh) : parent_t(mesh, nullptr) { parent_t::_hmesh._mfield = &this->_data; /*in order to have an access to the field values from hmesh*/}

    E_Int assign_data(const sensor_input_t& data) override;
    
    bool fill_adap_incr(output_t& adap_incr, bool do_agglo) override ;

    void fix_adap_incr(mesh_t& hmesh, output_t& adap_incr) override {};

    bool stop() override { return _single_pass_done; }

    
};

/// 
template <typename mesh_t>
E_Int metric_sensor<mesh_t>::assign_data(const sensor_input_t& data)
{
  E_Int npts = parent_t::_hmesh._crd.cols();

  assert(npts == data.size()); // input must be sized as number of pts (if decided to do differently,  need here to resize and fill missing field values)

  parent_t::_data = data;

  _single_pass_done = false;//reinit to enable

  return 0;
}


///
template <typename mesh_t>
bool metric_sensor<mesh_t>::fill_adap_incr(output_t& adap_incr, bool do_agglo)
{
  using cell_incr_t = typename output_t::cell_incr_t; //int_tuple<3>
  using face_incr_t = typename output_t::face_incr_t; // int_tuple<2>
  
  // init face_adap_incr
  adap_incr.face_adap_incr.clear();
  E_Int nb_faces = parent_t::_hmesh._ng.PGs.size();
  adap_incr.face_adap_incr.resize(nb_faces, face_incr_t(0));

  // init cell_adap_incr
  adap_incr.cell_adap_incr.clear();
  E_Int nb_cells = parent_t::_hmesh._ng.PHs.size();
  adap_incr.cell_adap_incr.resize(nb_cells, cell_incr_t(0));

  _single_pass_done = true; // to exit after first outer loop (based on sensor update)

  //todo Imad

  //1. compute face_adap_incr from nodal ellipses
  for (int i = 0; i < nb_faces; ++i)
  {
    face_incr_t f_inc;

    //dummy code :
    

    //todo
    f_inc.n[0] = f_inc.n[1] = parent_t::_data.getRadius(0); //dummy : 

    adap_incr.face_adap_incr[i] = f_inc;
  }

  //2. compute cell_adap_incr from face_adap_incr
  
  for (int i = 0; i < nb_cells; ++i)
  {
    cell_incr_t c_inc;

    //todo
    c_inc.n[0] = c_inc.n[1] = c_inc.n[2] = parent_t::_data.getRadius(0);; //dummy


    adap_incr.cell_adap_incr[i] = c_inc;
  }

 

  E_Int cmax = 0;
  for (size_t k = 0; (k < adap_incr.cell_adap_incr.size()) && (cmax==0); ++k)
    cmax = std::max(cmax, adap_incr.cmax(k));

  E_Int cmin = IDX_NONE; // max int32
  for (size_t k = 0; (k < adap_incr.cell_adap_incr.size()) && (cmin != 0); ++k)
    cmin = std::min(cmin, adap_incr.cmin(k));

  //std::cout << "cell adapt incr min/max : " << cmin << "/" << cmax << std::endl;

  return (cmin != 0 || cmax != 0);
}

}
#endif
