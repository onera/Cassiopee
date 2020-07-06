/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_NODAL_SENSOR_HXX
#define NUGA_NODAL_SENSOR_HXX

#include "Nuga/include/sensor.hxx"


namespace NUGA
{

template <typename mesh_t>
class nodal_sensor : public sensor<mesh_t, Vector_t<E_Int>> // Vector_t might be templatized to have ISO mode with a metric field
{
  public:
  
    using sensor_input_t = Vector_t<E_Int>;
    using parent_t = sensor<mesh_t, sensor_input_t>;
    using sensor_output_t = typename mesh_t::sensor_output_t; //fixme: static assert to add : must be ISO => IntVec

    nodal_sensor(mesh_t& mesh): parent_t(mesh, new V1_smoother<mesh_t>()){}

    E_Int assign_data(sensor_input_t& data);
    
    void fill_adap_incr(sensor_output_t& adap_incr, bool do_agglo) override;
    bool update() override;
};

/// 
template <typename mesh_t>
E_Int nodal_sensor<mesh_t>::assign_data(sensor_input_t& data)
{

  parent_t::assign_data(data);

  int ncrd = parent_t::_hmesh._crd.cols();
  if (ncrd > data.size())
    parent_t::_data.resize(ncrd, 0.);
}

///
template <typename mesh_t>
void nodal_sensor<mesh_t>::fill_adap_incr(sensor_output_t& adap_incr, bool do_agglo)
{
  //
  sensor_input_t& Ln = parent_t::_data;
  // E_Int n_nodes= parent_t::_hmesh._crd.size();
  // Ln.resize(n_nodes, 0);
  E_Int nb_faces = parent_t::_hmesh._ng.PGs.size();
  E_Int nb_elt= parent_t::_hmesh._ng.PHs.size();
  bool flag(false);
  adap_incr.cell_adap_incr.clear();
  adap_incr.cell_adap_incr.resize(nb_elt, 0);
  adap_incr.face_adap_incr.clear();
  adap_incr.face_adap_incr.resize(nb_faces, 0);

  for (int i=0; i< nb_elt; i++){
    if (parent_t::_hmesh._PHtree.is_enabled(i)){
      const E_Int* faces= parent_t::_hmesh._ng.PHs.get_facets_ptr(i);
      E_Int n_faces= parent_t::_hmesh._ng.PHs.stride(i);
      for (int k=0; k< n_faces; k++){
        E_Int PGk= *(faces+k)-1;
        const E_Int* pN= parent_t::_hmesh._ng.PGs.get_facets_ptr(PGk);
        E_Int n_face_nodes = parent_t::_hmesh._ng.PGs.stride(PGk);
        
        for (int l=0; l< n_face_nodes; l++){
          E_Int nodes_faces= *(pN+l)-1;
          if (Ln[nodes_faces]>0){
            adap_incr.cell_adap_incr[i]= 1;
            flag=true;
            break;
          }
        }
        if (flag==true) break;
      }
    }
    flag=false;
  }
}


///
template <typename mesh_t>
bool nodal_sensor<mesh_t>::update()
{
  sensor_input_t& Ln = parent_t::_data;
  E_Int nb_pts = parent_t::_hmesh._crd.size();
  Ln.resize(nb_pts, 0);
  
  for (int i=0; i< nb_pts; i++)
  {
    if (Ln[i]>0)
      Ln[i]--;
  }
  return false;
}

}
#endif
