/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_SMOOTHER_HXX
#define NUGA_SMOOTHER_HXX


namespace NUGA
{

/// Base smoother class
template <typename mesh_t>
struct smoother
{
  using sensor_output_t = typename sensor_output_data<mesh_t::SUBTYPE>::type;
     
  smoother() = default;

  virtual void smooth(const mesh_t& hmesh, sensor_output_t& adap_incr)  = 0;
  
  virtual ~smoother() {}
};

}
#endif
