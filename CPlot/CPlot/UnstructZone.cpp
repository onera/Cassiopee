/*    
    Copyright 2013-2019 Onera.

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

#include "UnstructZone.h"
#include <stdlib.h>

//=============================================================================
UnstructZone::UnstructZone( CPlotState* states, ZoneImpl* impl ) : Zone(states, impl)
{
  connect = NULL; 
  posFaces = NULL; 
  nelts1D = 0; nelts2D = 0;
  posElts1D = NULL; posElts2D = NULL;
}

//=============================================================================
UnstructZone::~UnstructZone()
{
  if (connect != NULL) delete [] connect;
  if (posFaces != NULL) delete [] posFaces;
  if (posElts1D != NULL) delete [] posElts1D;
  if (posElts2D != NULL) delete [] posElts2D;
}
