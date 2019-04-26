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
    
    static E_Int run(mesh_t& hmesh, sensor_t& sensor, typename sensor_t::data_type& data, bool do_agglo = false);
};

}


template <typename mesh_t, typename sensor_t>
E_Int NUGA::adaptor<mesh_t, sensor_t>::run(mesh_t& hmesh, sensor_t& sensor, typename sensor_t::data_type & data, bool do_agglo)
{
  E_Int err(0);
  
  Vector_t<E_Int> adap_incr;

  err = hmesh.init();
  if (err) return err;

  err = sensor.init(data);
  if (err) return err;

#ifdef DEBUG_2019
  std::cout << "avant" << std::endl;
  hmesh.quality_measure();
  
  std::cout << "///////////////////////////////////////////////" << std::endl;
  std::cout << "///////////////////////////////////////////////" << std::endl;
  E_Int k(0);
#endif

  while (!err && sensor.compute(data, adap_incr, do_agglo))
  {
#ifdef DEBUG_2019
    k=k+1;
#endif

    err = hmesh.adapt(adap_incr, do_agglo);

#ifdef DEBUG_2019
    std::cout << "k= " << k << std::endl;
    
    hmesh.quality_measure();
    std::cout << "///////////////////////////////////////////////" << std::endl;
    std::cout << "///////////////////////////////////////////////" << std::endl;
#endif
  
  }

#ifdef DEBUG_2019
  sensor.verif();
#endif

  return err;
}

#endif
