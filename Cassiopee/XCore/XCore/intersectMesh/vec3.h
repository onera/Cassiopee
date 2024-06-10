#pragma once

#include "xcore.h"

struct Vec3 {
    E_Float x, y, z;

    Vec3();

    Vec3(E_Float X, E_Float Y, E_Float Z);

    Vec3 operator-(const Vec3 &v) const;

    Vec3 operator+(const Vec3 &v) const;

    Vec3 operator*(E_Float a) const;
};

Vec3 operator*(E_Float a, const Vec3 &v);