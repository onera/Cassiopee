/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr), Alexis Gay (alexis.gay@onera.fr)

#include<chrono>

#ifndef NUGA_ADAPTOR_HXX
#define NUGA_ADAPTOR_HXX

//#define ADAPT_STEPS

#ifdef ADAPT_STEPS
#include <chrono>
#endif

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

#ifdef ADAPT_STEPS
  int iter = -1;
#endif

  while (!err)
  {
    E_Int nbphs = hmesh._ng.PHs.size();
#ifdef ADAPT_STEPS
    std::cout << "nuga/adapt::sensor.compute iter " << ++iter << "... with agglo ? " << do_agglo << std::endl;
    auto start0 = std::chrono::system_clock::now();
#endif

    bool require = sensor.compute(adap_incr, do_agglo);

#ifdef ADAPT_STEPS
    auto end0 = std::chrono::system_clock::now();
    auto t0 = end0 - start0;
    std::cout << "nuga/adapt::sensor.compute iter " << iter << " : " << t0.count() << "s" << std::endl;
#endif

    if (require)
    {
#ifdef ADAPT_STEPS
      start0 = std::chrono::system_clock::now();
      std::cout << "nuga/adapt::hmesh.adapt    iter " << iter << "..." << std::endl;
#endif
      err = hmesh.adapt(adap_incr, do_agglo);
    }

    bool updated = sensor.update(); 

    if (!require && !updated)
    {
#ifdef ADAPT_STEPS
      std::cout << "nuga/adapt : no more adaptation required." << std::endl;
#endif
      break;
    }

#ifdef ADAPT_STEPS
    end0 = std::chrono::system_clock::now();
    t0 = end0 - start0;
    std::cout << "nuga/adapt::hmesh.adapt    iter " << iter << " : " << t0.count() << "s" << std::endl;
#endif

  }

  return err;
}

#endif
