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
{}

Vec3::Vec3(E_Float X, E_Float Y, E_Float Z)
{
    ptr[0] = X;
    ptr[1] = Y;
    ptr[2] = Z;
}

Vec3 Vec3::operator-(const Vec3 &p) const
{
    Vec3 ret;
    for (E_Int i = 0; i < 3; i++) ret[i] = ptr[i] - p[i];
    return ret;
}

Vec3 Vec3::operator+(const Vec3 &p) const
{
    Vec3 ret;
    for (E_Int i = 0; i < 3; i++) ret[i] = ptr[i] + p[i];
    return ret;
}

Vec3 Vec3::operator*(E_Float a) const
{
    Vec3 ret;
    for (E_Int i = 0; i < 3; i++) ret[i] = a * ptr[i];
    return ret;
}

Vec3 operator*(E_Float a, const Vec3 & v)
{
    Vec3 ret;
    for (E_Int i = 0; i < 3; i++) ret[i] = a * v[i];
    return ret;
}

E_Float &Vec3::operator[](E_Int idx)
{
    return ptr[idx];
}

E_Float Vec3::operator[](E_Int idx) const
{
    return ptr[idx];
}
