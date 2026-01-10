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

//=============================================================================
// Keep information on the transcient software state for DL lists
//=============================================================================
#ifndef _CPLOT_STATEDL_H_
# define _CPLOT_STATEDL_H_
# include <vector>
# include "CPlotState.h"
# include "Def/DefTypes.h"

struct CPlotStateDL : public CPlotState
{
  // Display lists  
  std::vector<unsigned int> freeDLList; // list des DL a liberer

  virtual void freeGPUResources() {
    for (std::vector<unsigned int>::iterator itr = freeDLList.begin(); itr != freeDLList.end(); itr++)
      glDeleteLists(*itr, 1);
    freeDLList.clear();
    freeGPURes = 0; // nothing to clean any more
  }
};

#endif
