#include "AABB.h"
#include "mesh.h"

AABB::AABB()
: xmin(EFLOATMAX), ymin(EFLOATMAX), zmin(EFLOATMAX),
  xmax(EFLOATMIN), ymax(EFLOATMIN), zmax(EFLOATMIN)
{}

AABB::AABB(const IMesh &M, E_Int *ids, E_Int count)
{
    xmin = ymin = zmin = EFLOATMAX;
    xmax = ymax = zmax = EFLOATMIN;

    for (E_Int i = 0; i < count; i++) {
        E_Int fid = ids[i];

        const auto &pn = M.F[fid];

        for (E_Int p : pn) {
            E_Float x = M.X[p];
            E_Float y = M.Y[p];
            E_Float z = M.Z[p];
            if (x < xmin) xmin = x;
            if (y < ymin) ymin = y;
            if (z < zmin) zmin = z;
            if (x > xmax) xmax = x;
            if (y > ymax) ymax = y;
            if (z > zmax) zmax = z;
        }
    }
}