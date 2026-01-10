/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
//Authors : Sam Landier (sam.landier@onera.fr)

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

    masker(E_Float RTOL) : parent_t(RTOL), _col_X(0.5) {}
    
    outdata_t __process_X_cells(zmesh_t const & z_mesh, std::vector< bound_mesh_t*> const & mask_bits, wdata_t & wdata)
    {
      for (size_t i = 0; i < wdata.size(); ++i)
        wdata[i] = (wdata[i] == X) ? _col_X : wdata[i];

      return wdata;
    };

  private:
    E_Float _col_X;

  };
}
#endif
