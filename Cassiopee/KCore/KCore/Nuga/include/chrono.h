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

#ifndef __DELAUNAY_CHRONO_H__
#define __DELAUNAY_CHRONO_H__

#include "time.h"
#include "Nuga/include/defs.h"

namespace NUGA
{
  class chrono
  {
    public:
      chrono():_t0(0){}    
      void start(){_t0 = ::clock();}
      E_Float elapsed(){return E_Float(::clock() - _t0) / CLOCKS_PER_SEC;}
      
    private:
      clock_t _t0;
  }; 
}

#endif
