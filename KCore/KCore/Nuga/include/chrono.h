/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

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
