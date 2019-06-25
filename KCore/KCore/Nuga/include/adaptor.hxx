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
    
    static constexpr eSUBDIV_TYPE sub_type = mesh_t::sub_type;
    using output_type = typename mesh_t::output_type;
    
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
  
  output_type adap_incr;

  err = hmesh.init();
  if (err) return err;

  err = sensor.init(data);
  if (err) return err;

#ifdef DEBUG_2019
  //  E_Float Vol_init = hmesh.vol_init();
  //std::cout << "avant" << std::endl;
  //hmesh.quality_measure();
// 
  //std::cout << "_ng.PHs size= " <<hmesh._ng.PHs.size() << std::endl;
  //hmesh.check_vol(Vol_init, false);
  //std::cout << "///////////////////////////////////////////////" << std::endl;
  //std::cout << "///////////////////////////////////////////////" << std::endl;
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
//  sensor.verif();
//  hmesh.verif3();
//  hmesh.verif4();
//  auto end = std::chrono::system_clock::now();
//  std::chrono::duration<double> t1 = end-start;
//  std::cout << "tps verif= " << t1.count() << "s"<< std::endl;
    Vector_t<E_Int> s(10,0);
  for (int i=0; i< hmesh._ng.PHs.size(); i++){
      if (hmesh._PHtree.is_enabled(i)){
          E_Int level= hmesh._PHtree.get_level(i);
          for (int j=0; j<10; j++){
              if (level==j){
                  s[j] += 1 ;
              }
          }
      }
  }
  E_Int sT(0);
  for (int i=0; i<10; i++){
    std::cout << "s[" << i << "]= " << s[i] << std::endl;
    sT += s[i];
  }
  std::cout << "sT= " << sT << std::endl;
  std::cout << "///////////////////////"<< std::endl;
  hmesh.verif5();

#endif

  return err;
}

#endif
