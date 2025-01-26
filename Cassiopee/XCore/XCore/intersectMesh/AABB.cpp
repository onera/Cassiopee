/*    
    Copyright 2013-2024 Onera.

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
#include "AABB.h"

void AABB_clamp(AABB &box, const AABB &parent)
{
    box.xmin = std::max(parent.xmin, box.xmin);
    box.ymin = std::max(parent.ymin, box.ymin);
    box.zmin = std::max(parent.zmin, box.zmin);
    box.xmax = std::min(parent.xmax, box.xmax);
    box.ymax = std::min(parent.ymax, box.ymax);
    box.zmax = std::min(parent.zmax, box.zmax);
}