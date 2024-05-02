#include "proto.h"

int sign(double x)
{
    if (x > TOL)
        return 1;
    if (x < -TOL)
        return -1;
    return 0;
}

int sign(int i)
{
    if (i > 0)
        return 1;
    if (i < 0)
        return -1;
    return 0;
}

int xy_cmp(double x0, double y0, double x1, double y1)
{
    int cmp = sign(x0 - x1);
    if (cmp)
        return cmp;
    cmp = sign(y0 - y1);
    return cmp;
}
