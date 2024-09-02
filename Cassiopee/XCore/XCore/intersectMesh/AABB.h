#pragma once

struct AABB {
    E_Float xmin = EFLOATMAX;
    E_Float xmax = EFLOATMIN;
    E_Float ymin = EFLOATMAX;
    E_Float ymax = EFLOATMIN;
    E_Float zmin = EFLOATMAX;
    E_Float zmax = EFLOATMIN;
};