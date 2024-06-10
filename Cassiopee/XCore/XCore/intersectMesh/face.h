#pragma once

#include "xcore.h"

struct Hedge;

struct Face {
    Hedge *rep;
    E_Int oid[2];

    Face();
};
