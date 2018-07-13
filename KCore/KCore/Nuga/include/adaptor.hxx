/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr), Alexis Gay (alexis.gay@onera.fr)

#ifndef NUGA_ADAPTOR_HXX
#define NUGA_ADAPTOR_HXX

namespace NUGA
{
  
template <typename mesh_t, typename sensor_t>
class adaptor
{
  public:
    
    using elt_type = typename mesh_t::elt_type;
    
    static E_Int run(mesh_t& hmesh, sensor_t& sensor, typename sensor_t::data_type& data);
};

}


template <typename mesh_t, typename sensor_t>
E_Int NUGA::adaptor<mesh_t, sensor_t>::run(mesh_t& hmesh, sensor_t& sensor, typename sensor_t::data_type & data)
{
  E_Int err(0);
  
  Vector_t<E_Int> adap_incr;
  
  err = hmesh.init();
  if (err) return err;

  err = sensor.init(data);
  if (err) return err;

  while (!err && sensor.compute(data, adap_incr))
    err = hmesh.adapt(adap_incr);

  return err;
}

#endif
