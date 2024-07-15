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
#include "ray.h"
#include "primitives.h"

E_Int MollerTrumbore(E_Float px, E_Float py, E_Float pz, E_Float dx, E_Float dy,
    E_Float dz, E_Float ax, E_Float ay, E_Float az, E_Float bx, E_Float by, E_Float bz,
    E_Float cx, E_Float cy, E_Float cz, TriangleIntersection &TI)
{
    E_Float e1x = bx - ax;
    E_Float e1y = by - ay;
    E_Float e1z = bz - az;

    E_Float e2x = cx - ax;
    E_Float e2y = cy - ay;
    E_Float e2z = cz - az;

    E_Float pvecx = dy * e2z - dz * e2y;
    E_Float pvecy = dz * e2x - dx * e2z;
    E_Float pvecz = dx * e2y - dy * e2x;

    E_Float det = e1x * pvecx + e1y * pvecy + e1z * pvecz;

    if (det > -TOL && det < TOL)
        return 0;

    E_Float inv_det = 1.0 / det;

    E_Float tvecx = px - ax;
    E_Float tvecy = py - ay;
    E_Float tvecz = pz - az;

    TI.u = inv_det * (tvecx * pvecx + tvecy * pvecy + tvecz * pvecz);

    if (TI.u < 0 || TI.u > 1)
        return 0;

    E_Float qvecx = tvecy * e1z - tvecz * e1y;
    E_Float qvecy = tvecz * e1x - tvecx * e1z;
    E_Float qvecz = tvecx * e1y - tvecy * e1x;

    TI.v = inv_det * (dx * qvecx + dy * qvecy + dz * qvecz);

    if (TI.v < 0 || TI.u + TI.v > 1)
        return 0;
    
    TI.t = inv_det * (e2x * qvecx + e2y * qvecy + e2z * qvecz);

    if (TI.t > TOL) {
        TI.x = px + TI.t * dx;
        TI.y = py + TI.t * dy;
        TI.z = pz + TI.t * dz;
        return 1;
    }
    
    return 0;
}