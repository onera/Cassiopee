/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_SMOOTHER_HXX
#define NUGA_SMOOTHER_HXX

#include "Nuga/include/subdivision.hxx"


namespace NUGA
{
  enum eSmoother { V1_NEIGH = 0, SHELL = 1 };

/// Base smoother class
template <typename mesh_t>
struct smoother
{
  using output_t = adap_incr_type<mesh_t::SUBTYPE>;
     
  smoother() = default;

  virtual void smooth(const mesh_t& hmesh, output_t& adap_incr)  = 0;
  
  virtual ~smoother() {}
};

}
#endif
