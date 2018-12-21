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
#ifndef _CPLOT_ZONEIMPLDL_H_
#define _CPLOT_ZONEIMPLDL_H_
# include "ZoneImpl.h"

class ZoneImplDL : public ZoneImpl
{
 public:
  ZoneImplDL();
  unsigned int _DLmesh; // Display list for mesh lines
  unsigned int _DLsolid; // Display list for solid
  unsigned int _DLiso; // Display list pour les isos
  int _DLisoField; // Champ actuellement contenu dans la DLiso
  int _DLisoField2; // champ 2 si champ vectoriel dans DLiso
  int _DLisoField3; // champ 3 si champ vectoriel dans DLiso

  virtual void freeGPURes( CPlotState* pt_state, bool freeIso = true );
  virtual void destroyIsoField();
};

#endif
