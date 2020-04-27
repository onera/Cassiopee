/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr), Alexis Gay (alexis.gay@onera.fr)

#include<chrono>

#ifndef NUGA_ADAPTOR_HXX
#define NUGA_ADAPTOR_HXX

namespace NUGA
{
  
template <typename mesh_t, typename sensor_t>
class adaptor
{
  public:
      
    static E_Int run(mesh_t& hmesh, sensor_t& sensor, bool do_agglo = false);
  
};

}


template <typename mesh_t, typename sensor_t>
E_Int NUGA::adaptor<mesh_t, sensor_t>::run(mesh_t& hmesh, sensor_t& sensor, bool do_agglo)
{

  E_Int err(0);
  
  typename mesh_t::sensor_output_t adap_incr;

  hmesh.init();  

  //int iter = 0;
  while (!err)
  {

#ifdef ADAPT_STEPS
    std::cout << "nuga/adapt::sensor.compute iter " << iter++ << "..." << std::endl;
    start0 = std::chrono::system_clock::now();
#endif
    bool require = sensor.compute(adap_incr, do_agglo);
#ifdef ADAPT_STEPS
    end0 = std::chrono::system_clock::now();
    t0 = end0 - start0;
    std::cout << "nuga/adapt::sensor.compute iter " << iter << " : " << t0.count() << "s" << std::endl;
#endif

    if (!require)
    {
#ifdef ADAPT_STEPS
      std::cout << "nuga/adapt : no more adaptation required." << std::endl;
#endif
      break;
    }

#ifdef ADAPT_STEPS
    start0 = std::chrono::system_clock::now();
    std::cout << "nuga/adapt::hmesh.adapt    iter " << iter << "..." << std::endl;
#endif

    err = hmesh.adapt(adap_incr, do_agglo);

    sensor.update(); 

#ifdef ADAPT_STEPS
    end0 = std::chrono::system_clock::now();
    t0 = end0 - start0;
    std::cout << "nuga/adapt::hmesh.adapt    iter " << iter++ << " : " << t0.count() << "s" << std::endl;
#endif

  }

#ifdef ADAPT_STEPS
  end0 = std::chrono::system_clock::now();
  t0 = end0 - startG;
  std::cout << "nuga/adapt : total time : " << t0.count() << "s" << std::endl;
#endif

  return err;
}

#endif
