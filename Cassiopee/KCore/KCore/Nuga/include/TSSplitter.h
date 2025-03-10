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

#ifndef _TSSPLITTER_H_
#define _TSSPLITTER_H_

#include "Nuga/include/DynArray.h"
#include <vector>

class TSSplitter
{
  public:
    /* Splits the connectivity (only) with the edge N0N1 and stores the 
       bits into separate containers.
       WARNING: it is the responsibility of the caller to destroy connectBout.
    */
    static void split(
      const K_FLD::IntArray& connectBin, const K_FLD::IntArray& polyLine,
      std::vector<K_FLD::IntArray> & connectBout);

  /* Splits the connectivity with the edge N0N1 and stores the connectivity 
     bits and corresponding coordinates into separate containers.
     WARNING: it is the responsibility of the caller to destroy posOut and 
     connectBout. */
    static void split(
      const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectBin, 
      const K_FLD::IntArray& polyLine,
      std::vector<K_FLD::FloatArray> & posOut, 
      std::vector<K_FLD::IntArray> & connectBout);

  private:
    TSSplitter(void){}
    ~TSSplitter(void){}
};

#endif
