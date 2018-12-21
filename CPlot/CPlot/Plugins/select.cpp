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
#include "../Data.h"

//=============================================================================
// Selection plugins
//=============================================================================
// selectNextZone:
// modify _state.selectedZone, and z.selected, no other modification.

// selectPreviousZone:
// modify _state.selectedZone, and z.selected, no other modification.

//=============================================================================
/* 
   Incremental selection
   Select zones by number.
*/
//=============================================================================
void selectNextZoneIncr(Data* d)
{
  CPlotState* s = d->ptrState;
  if (s->selectedZone != 0) d->_zones[s->selectedZone-1]->selected = 0;
  s->selectedZone++;
  if (s->selectedZone == d->_numberOfZones+1)
    s->selectedZone = 0;
  if (s->selectedZone != 0) d->_zones[s->selectedZone-1]->selected = 1;
}

void selectPreviousZoneIncr(Data* d)
{
  CPlotState* s = d->ptrState;
  if (s->selectedZone != 0) d->_zones[s->selectedZone-1]->selected = 0;
  s->selectedZone--;
  if (s->selectedZone == -1)
    s->selectedZone = d->_numberOfZones;
  if (s->selectedZone != 0) d->_zones[s->selectedZone-1]->selected = 1;
}

//=============================================================================
/*
  Nearest neighbour selection
  Select the nearest zone of selected zone.
*/
//=============================================================================
void selectNextZoneNearest(Data* d)
{

}

void selectPreviousZoneNearest(Data* d)
{

}
