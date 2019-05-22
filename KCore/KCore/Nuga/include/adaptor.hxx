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
  #ifdef DEBUG_2019  
  auto start0 = std::chrono::system_clock::now();
  #endif

  
  E_Int err(0);
  
  Vector_t<E_Int> adap_incr;

  err = hmesh.init();
  if (err) return err;

  err = sensor.init(data);
  if (err) return err;

#ifdef DEBUG_2019
  //  E_Float Vol_init = hmesh.vol_init();
  std::cout << "avant" << std::endl;
  hmesh.quality_measure();
// 
  //std::cout << "_ng.PHs size= " <<hmesh._ng.PHs.size() << std::endl;
  //hmesh.check_vol(Vol_init, false);
  std::cout << "///////////////////////////////////////////////" << std::endl;
  std::cout << "///////////////////////////////////////////////" << std::endl;
  E_Int k(0);
  auto end0 = std::chrono::system_clock::now();
  std::chrono::duration<double> t0 = end0-start0;
  std::cout << "tps init= " << t0.count() << "s"<< std::endl;
#endif

  while (!err && sensor.compute(data, adap_incr, do_agglo))
  {
#ifdef DEBUG_2019
    //std::cout << "nouveau adapt" << std::endl;
    k=k+1;
    auto start = std::chrono::system_clock::now();
#endif

    err = hmesh.adapt(adap_incr, do_agglo);

#ifdef DEBUG_2019
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> t1 = end-start;
    std::cout << "k= " << k << " tps adapt= " << t1.count() << "s"<< std::endl;
    
    //std::cout << "k= " << k << std::endl;
    //E_Float tol= 1e-6;
    //hmesh.check_vol(Vol_init, false);

    //hmesh.quality_measure();
//    std::cout << "///////////////////////////////////////////////" << std::endl;
//    std::cout << "///////////////////////////////////////////////" << std::endl;
#endif
  
  }

#ifdef DEBUG_2019
    hmesh.quality_measure();
    std::cout << "///////////////////////////////////////////////" << std::endl;
    std::cout << "///////////////////////////////////////////////" << std::endl;

//  auto start = std::chrono::system_clock::now();  
//  E_Int v=sensor.verif();
//  v++;
//  hmesh.verif3();
//  hmesh.verif4();
//  auto end = std::chrono::system_clock::now();
//  std::chrono::duration<double> t1 = end-start;
//  std::cout << "tps verif= " << t1.count() << "s"<< std::endl;
#endif

  return err;
}

#endif
