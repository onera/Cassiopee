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