/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_MASKER_HXX
#define NUGA_MASKER_HXX

#include "Nuga/include/classifyer.hxx"

namespace NUGA
{
  ///
  template<typename zmesh_t, typename bound_mesh_t = typename NUGA::boundary_t<zmesh_t>>
  class masker : public classifyer<COLLISION, zmesh_t, bound_mesh_t>
  {
  public:
    using parent_t = classifyer<COLLISION, zmesh_t, bound_mesh_t>;
    using wdata_t = typename parent_t::wdata_t;
    using outdata_t = typename parent_t::outdata_t;

    masker(double RTOL) : parent_t(RTOL), _col_X(0.5) {}
    
    outdata_t __process_X_cells(zmesh_t const & z_mesh, std::vector< bound_mesh_t*> const & mask_bits, wdata_t & wdata)
    {
      for (size_t i = 0; i < wdata.size(); ++i)
        wdata[i] = (wdata[i] == X) ? _col_X : wdata[i];

      return wdata;
    };

  private:
    double _col_X;

    
  };
}
#endif
