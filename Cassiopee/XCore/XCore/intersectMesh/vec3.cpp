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
#include "vec3.h"

Vec3::Vec3()
: x(), y(), z()
{}

Vec3::Vec3(E_Float X, E_Float Y, E_Float Z)
: x(X), y(Y), z(Z)
{}

Vec3 Vec3::operator-(const Vec3 &p) const
{
    return { x - p.x, y - p.y, z - p.z };
}

Vec3 Vec3::operator+(const Vec3 &p) const
{
    return { x + p.x, y + p.y, z + p.z };
}

Vec3 Vec3::operator*(E_Float a) const
{
    return { a * x, a * y, a * z };
}

Vec3 operator*(E_Float a, const Vec3 & v)
{
    return { a * v.x, a * v.y, a * v.z };
}