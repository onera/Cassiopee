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
  
    using sensor_data_t = Vector_t<E_Int>;
    using parent_t = sensor<mesh_t, sensor_data_t>;
    using sensor_output_t = typename mesh_t::sensor_output_t; //fixme: static assert to add : must be ISO => IntVec

 
    nodal_sensor(mesh_t& mesh): parent_t(mesh){}
    //E_Int assign_data(sensor_data_t& level);
    bool compute(sensor_output_t& adap_incr, bool do_agglo) override;
    void update() override;
};

//template <typename mesh_t>
//E_Int nodal_sensor<mesh_t>::assign_data(sensor_data_t& level)
//{
//  _Ln =level;
//  return 0;  
//}      

template <typename mesh_t>
bool nodal_sensor<mesh_t>::compute(sensor_output_t& adap_incr, bool do_agglo)
{
  sensor_data_t& _Ln = *parent_t::_data;
  // E_Int n_nodes= parent_t::_hmesh._crd->size();
  // _Ln.resize(n_nodes, 0);
  E_Int nb_elt= parent_t::_hmesh._ng->PHs.size();
  bool flag(false);
  adap_incr.clear();
  adap_incr.resize(nb_elt, 0);

  for (int i=0; i< nb_elt; i++){
    if (parent_t::_hmesh._PHtree.is_enabled(i)){
      E_Int* faces= parent_t::_hmesh._ng->PHs.get_facets_ptr(i);
      E_Int n_faces= parent_t::_hmesh._ng->PHs.stride(i);
      for (int k=0; k< n_faces; k++){
        E_Int PGk= *(faces+k)-1;
        E_Int* pN= parent_t::_hmesh._ng->PGs.get_facets_ptr(PGk);
        E_Int n_face_nodes = parent_t::_hmesh._ng->PGs.stride(PGk);
        
        for (int l=0; l< n_face_nodes; l++){
          E_Int nodes_faces= *(pN+l)-1;
          if (_Ln[nodes_faces]>0){
            adap_incr[i]= 1;
            flag=true;
            break;
          }
        }
        if (flag==true) break;
      }
    }
    flag=false;
  }

  //apply the 2:1 rule
  parent_t::_hmesh.smooth(adap_incr);

  // fix inconsistencies
  fix_adap_incr(parent_t::_hmesh, adap_incr);
    
  for (int i=0; i< adap_incr.size(); i++){
    if (adap_incr[i]==1){
      return true;
    }
  }  
  return false;
}


///
template <typename mesh_t>
void nodal_sensor<mesh_t>::update()
{
  sensor_data_t& _Ln = *parent_t::_data;
  E_Int n_nodes= parent_t::_hmesh._crd->size();
  _Ln.resize(n_nodes, 0);
  
  for (int i=0; i< n_nodes; i++)
  {
    if (_Ln[i]>0)
      _Ln[i]--;
  }
}

}
#endif
