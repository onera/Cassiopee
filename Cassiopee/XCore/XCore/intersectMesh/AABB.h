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
#pragma once

#include "xcore.h"

struct AABB {
    E_Float xmin;
    E_Float ymin;
    E_Float zmin;
    E_Float xmax;
    E_Float ymax;
    E_Float zmax;
    E_Float dx;
    E_Float dy;
    E_Float dz;

    bool is_point_inside(const E_Float p[3]) const
    {
        return (p[0] >= xmin && p[0] <= xmax && p[1] >= ymin && p[1] <= ymax &&
            p[2] >= zmin && p[2] <= zmax);
    }

    void print() const
    {
        printf("[%f %f %f] - [%f %f %f]\n", xmin, ymin, zmin, xmax, ymax, zmax);
    }
};

const AABB AABB_HUGE = {EFLOATMIN, EFLOATMIN, EFLOATMIN,
                        EFLOATMAX, EFLOATMAX, EFLOATMAX};

void AABB_clamp(AABB &box, const AABB &parent);
