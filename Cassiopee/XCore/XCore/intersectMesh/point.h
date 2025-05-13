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

struct Point {
    E_Float x, y, z;
    bool operator<(const Point &p) const
    {
        return (x < p.x) || (x == p.x && y < p.y) ||
            (x == p.x && y == p.y && z < p.z);
    }
    Point operator-(const Point &p) const
    {
        return {x-p.x, y-p.y, z-p.z};
    }
    E_Float &operator[](E_Int idx)
    {
        return (idx == 0) ? x : ((idx == 1) ? y : z);
    }
    E_Float operator[](E_Int idx) const
    {
        return (idx == 0) ? x : ((idx == 1) ? y : z);
    }
};

typedef Point Point3D;

struct Point2D
{
    E_Float x, y;
};

struct PointLoc {
    E_Int fid = -1;
    E_Int v_idx = -1;
    E_Int e_idx = -1;
    E_Float bcrd[3] = {EFLOATMAX, EFLOATMAX, EFLOATMAX};
    E_Float t = EFLOATMAX;
    E_Float x = EFLOATMAX;
    E_Float y = EFLOATMAX;
    E_Float z = EFLOATMAX;
    E_Int sub = -1;
};

struct PointData {
    E_Int pid; // my id in spatch
    E_Int fid; // id of the mface I am on
    E_Int sid; // if of the smesh I belong to
    E_Float x, y, z;
};

struct pointFace {
    E_Int F;

    pointFace(E_Int face);
};
