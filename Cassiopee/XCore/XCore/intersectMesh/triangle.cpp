#include <limits>
#include <cmath>

#include "triangle.h"
#include "primitives.h"

static
E_Int orient(E_Float ax, E_Float ay, E_Float bx, E_Float by, E_Float cx, E_Float cy)
{
    E_Float det = (bx - ax)*(cy - ay) - (by - ay)*(cx - ax);
    return Sign(det);
}

E_Int Triangle::ispointInside(E_Float px, E_Float py, E_Float ax, E_Float ay,
    E_Float bx, E_Float by, E_Float cx, E_Float cy)
{
    if (orient(ax, ay, bx, by, px, py) < 0) return 0;
    if (orient(bx, by, cx, cy, px, py) < 0) return 0;
    if (orient(cx, cy, ax, ay, px, py) < 0) return 0;
    
    return 1;
}