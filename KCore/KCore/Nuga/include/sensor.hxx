/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_SENSOR_HXX
#define NUGA_SENSOR_HXX

#include<vector>
#include "Nuga/include/macros.h"
#include "subdivision.hxx"

namespace NUGA
{

template <eSUBDIV_TYPE STYPE>
struct sensor_output_data;

template<>
struct sensor_output_data<ISO>
{
  using type = Vector_t<E_Int>;
};

template<>
struct sensor_output_data<DIR>
{
  using type = dir_incr_type;
};

//template<>
//struct sensor_output_data<ISO2>
//{
//  using type = Vector_t<E_Float>;
//};

/// Geometric sensor
template <typename mesh_t, typename sensor_data_t>
class sensor
{
  public:
 
    using sensor_output_t = typename sensor_output_data<mesh_t::SUBTYPE>::type;
    
    sensor(mesh_t& mesh) :_hmesh(mesh) {};
    
    virtual E_Int assign_data(sensor_data_t& data) { _data = &data; return 0; }
    
    virtual bool compute(sensor_output_t& adap_incr, bool do_agglo) = 0;
  
    virtual void update() {};

    virtual ~sensor() = default;

  protected:
    mesh_t const & _hmesh;
    sensor_data_t* _data;
};

}
#endif
