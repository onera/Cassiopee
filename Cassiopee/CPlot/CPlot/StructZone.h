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

#ifndef _CPLOT_STRUCTZONE_H_
#define _CPLOT_STRUCTZONE_H_

#include "Zone.h"
//#include "kcore.h"

/* Define a structured zone */
class StructZone : public Zone
{
  public:
    StructZone(CPlotState* states, ZoneImpl* impl);
    virtual ~StructZone();
    void compNorm();
      
  public:
    E_Int ni, nj, nk;                  // number of points
 
    E_Int activePlane; // not used, i=0, j=1, k=2
    // In 2D or 3D mode, CPlot displays 3 2D planes (i=cste, j=cste, k=cste)
    E_Int iPlane; // from 0 to ni-1, -1 means min-max planes
    E_Int jPlane;
    E_Int kPlane;
    // In 1D mode, CPlot displays only one line:
    // In 1D i mode: line j = jPlane, k = kLine
    // In 1D j mode: line i = iLine, k = kPlane
    // In 1D k mode: line i = iPlane, j = jLine
    E_Int iLine;
    E_Int jLine;
    E_Int kLine;
};
#endif
