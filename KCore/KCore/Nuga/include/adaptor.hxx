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
    
    using output_type = typename mesh_t::output_type;
    
    static E_Int run(mesh_t& hmesh, sensor_t& sensor, typename sensor_t::data_type& data, bool do_agglo = false);
};

}


template <typename mesh_t, typename sensor_t>
E_Int NUGA::adaptor<mesh_t, sensor_t>::run(mesh_t& hmesh, sensor_t& sensor, typename sensor_t::data_type & data, bool do_agglo)
{

  E_Int err(0);
  
  output_type adap_incr;

#ifdef ADAPT_STEPS
  std::cout << "nuga/adapt::hmesh.init..." << std::endl;
  auto start0 = std::chrono::system_clock::now();
  auto startG = start0;
#endif
  
  err = hmesh.init();
  if (err) return err;

#ifdef ADAPT_STEPS
  auto end0 = std::chrono::system_clock::now();
  std::chrono::duration<double> t0 = end0 - start0;
  std::cout << "nuga/adapt::hmesh.init : " << t0.count() << "s" << std::endl;
  start0 = std::chrono::system_clock::now();
  std::cout << "nuga/adapt::sensor.init..." << std::endl;
#endif

  err = sensor.init(data);
  if (err) return err;

#ifdef ADAPT_STEPS
  end0 = std::chrono::system_clock::now();
  t0 = end0-start0;
  std::cout << "nuga/adapt::sensor.init : " << t0.count() << "s"<< std::endl;
  int iter = 1;
#endif

  while (!err)
  {

#ifdef ADAPT_STEPS
    std::cout << "nuga/adapt::sensor.compute iter " << iter << "..." << std::endl;
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
