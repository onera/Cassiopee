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
#include "../Data.h"

//=============================================================================
// Compute the distance to a bounding box
//=============================================================================
double Data::dist2BB(double x, double y, double z,
double xmin, double ymin, double zmin,
double xmax, double ymax, double zmax)
{
  double xp, yp, zp;

  if (x <= xmin) xp = xmin;
  else if (x >= xmax) xp = xmax;
  else xp = x;

  if (y <= ymin) yp = ymin;
  else if (y >= ymax) yp = ymax;
  else yp = y;

  if (z <= zmin) zp = zmin;
  else if (z >= zmax) zp = zmax;
  else zp = z;

  return (xp-x)*(xp-x)+(yp-y)*(yp-y)+(zp-z)*(zp-z);
}
