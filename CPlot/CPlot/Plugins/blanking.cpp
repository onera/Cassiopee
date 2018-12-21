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
// Blanking plugins
//=============================================================================

// Parameters :
// IN: p1 -> point index
//     blank -> no of the blanking variable
//     zone -> no of zone the point belongs to.
// OUT: return 0 or 1 depending if p1 is blanked or not.

//=============================================================================
/* 
   blankCellN.
   Tell if p1 is blanked following cellN variable.
   This variable is: blanked (0), interpolated or discretized (1)
*/
//=============================================================================
int blankCellN(Data* d, int p1, int blank, int zone)
{
  Zone* z = d->_zones[zone];
  double* cellN = z->f[blank-1]; 
  double val = cellN[p1];
  if (val >= 0 && val < 1) return 0; // blanked
  else return 1;
}

//=============================================================================
/* 
   blankCellNF.
   Tell if p1 point is blanked following cellNF variable.
   This variable is : blanked (0), discretized (1), orphan (-1000),
   interpolated (-noOfInterpolatedDomain)
*/
//=============================================================================
int blankCellNF(Data* d, int p1, int blank, int zone)
{
  Zone* z = d->_zones[zone];
  double* cellN = z->f[blank-1];
  double val = cellN[p1];
  if (val == 0) return 0; // blanked
  else return 1;
}

//=============================================================================
/* 
   blankStatus.
   Tell if p1 point is blanked following status variable.
   This variable is: blanked (0), discretized (1), orphan (-1000),
   interpolated (-noOfInterpolatedDomain)
*/
//=============================================================================
int blankStatus(Data* d, int p1, int blank, int zone)
{
  Zone* z = d->_zones[zone];
  double* cellN = z->f[blank-1];
  double val = cellN[p1];
  if (val >= -0.1 && val <= 0.1) return 0; // blanked
  else return 1;
}
