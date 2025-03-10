/*    
    Copyright 2013-2025 Onera.

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
