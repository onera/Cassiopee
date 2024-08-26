#pragma once

struct AABB {
    E_Float xmin = FLOATMAX;
    E_Float xmax = FLOATMIN;
    E_Float ymin = FLOATMAX;
    E_Float ymax = FLOATMIN;
    E_Float zmin = FLOATMAX;
    E_Float zmax = FLOATMIN;
};