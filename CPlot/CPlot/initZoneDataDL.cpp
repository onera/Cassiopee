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
# include "kcore.h"
# include "DataDL.h"
# include "cplot.h"
# include "CPlotStateDL.h"
# include "ZoneImplDL.h"
#include <time.h>
#if defined(_WIN32) || defined(_WIN64)
#  include <winsock.h>
#endif

using namespace K_FLD;
using namespace std;

//=============================================================================
// Initialise les donnees des zones a partir des arrays.
// si retourne 0: echec
// si retourne 1: OK
//=============================================================================
int
DataDL::initZoneData(vector<FldArrayF*>& structF,
			vector<char*>& structVarString,
			vector<E_Int>& nit,
			vector<E_Int>& njt,
			vector<E_Int>& nkt,
			vector<FldArrayF*>& unstrF,
			vector<char*>& unstrVarString,
			vector<FldArrayI*>& cnt,
			vector<char*>& eltType,
			vector<char*>& zoneNames,
			vector<char*>& zoneTags,
      E_Int referenceNfield,
      char** referenceVarNames)
{
  // Dit a display de liberer les DL des zones
  ptrState->syncGPURes();
  for (int i = 0; i < _numberOfZones; i++) {
    ZoneImplDL* z;
    z = static_cast<ZoneImplDL*>(_zones[i]->ptr_impl);
    z->freeGPURes(ptrState); z->_GPUResUse = 0; 
  }

  return Data::initZoneData(structF, structVarString, nit, njt, nkt,
			     unstrF, unstrVarString, cnt, eltType, zoneNames, zoneTags,
           referenceNfield, referenceVarNames);
}
