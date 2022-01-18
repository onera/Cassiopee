/*
 
 
 
 --------- NUGA v1.0 
 
 
 
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

  enum ePara { SEQ, FINE_OMP, COARSE_OMP /*, DISTRIB*/ };

///
template <typename hmesh_t, typename sensor_t>
class adaptor
{
  public:

    static E_Int run(hmesh_t& hmesh, sensor_t& sensor, bool do_agglo = false);

    template <typename communicator_t>
    static E_Int run(std::vector<hmesh_t*>& hzones, std::vector<sensor_t*>& sensors, bool do_agglo, ePara PARA, communicator_t* com);
  
};

}

///
template <typename hmesh_t, typename sensor_t>
E_Int NUGA::adaptor<hmesh_t, sensor_t>::run(hmesh_t& hmesh, sensor_t& sensor, bool do_agglo)
{

  E_Int err(0);
  
  typename hmesh_t::output_t adap_incr;

  hmesh.init();  

#ifdef ADAPT_STEPS
  int iter = -1;
#endif

  while (!err)
  {
    //E_Int nbphs = hmesh._ng.PHs.size();
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

///
template <typename mesh_t, typename sensor_t>
template <typename communicator_t>
E_Int NUGA::adaptor<mesh_t, sensor_t>::run(std::vector<mesh_t*>& hzones, std::vector<sensor_t*>& sensors, bool do_agglo, ePara PARA, communicator_t* COM)
{
  E_Int err(0);
  E_Int NBZ{ E_Int(hzones.size()) };

  assert ( NBZ == sensors.size() );

#pragma omp parallel for if(PARA == COARSE_OMP)
  for (E_Int i = 0; i < NBZ; ++i)
  {
    //std::cout << "running zone :" << i << std::endl;
    if (hzones[i] == nullptr || sensors[i] == nullptr) continue;

    NUGA::adaptor<mesh_t, sensor_t>::run(*hzones[i], *sensors[i], do_agglo);
  }

  if (COM == nullptr) return 0;

  // ADAPT PARA :: A POSTERIORI JOIN SMOOTHING
  //std::cout << "post smoothing" << std::endl;
  bool has_changes{ true };
  //E_Int iter{ 0 };
  while (has_changes)
  {
    has_changes = false;

    //std::cout << "smoothing iter : " << iter++ << std::endl;

    // prepare & send
    for (size_t i = 0; i < COM->agents.size(); ++i)
    {
      //std::cout << "agent : " << i << " : " << COM->agents[i] << std::endl;
      if (COM->agents[i] == nullptr) continue;
      //std::cout << "pack" << std::endl;
      bool has_packs = COM->agents[i]->pack();
      //if (has_packs) std::cout << "isend" << std::endl;
      if (has_packs) COM->isend(i);
    }

    // receive and process
    for (size_t i = 0; i < COM->agents.size(); ++i)
    {
      //std::cout << "agent : " << i << " : " << COM->agents[i] << std::endl;
      if (COM->agents[i] == nullptr) continue;
      //std::cout << "unpack : " << i << std::endl;
      bool has_packs = (COM->agents[i])->unpack(); //should wait untill all data are received
      //if (has_packs) std::cout << "process : " << i << std::endl;
      if (has_packs) has_changes |= (COM->agents[i])->process();
      //if (has_packs) std::cout << "process done" << std::endl;
    }
  }

  // done : all zones are sync (2:1 smoothing)
  //std::cout << "adapt multiseq done" << std::endl;
  return err;
}

#endif
