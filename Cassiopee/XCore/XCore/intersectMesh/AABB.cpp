#include "AABB.h"

void AABB_clamp(AABB &box, const AABB &parent)
{
    box.xmin = std::max(parent.xmin, box.xmin);
    box.ymin = std::max(parent.ymin, box.ymin);
    box.zmin = std::max(parent.zmin, box.zmin);
    box.xmax = std::min(parent.xmax, box.xmax);
    box.ymax = std::min(parent.ymax, box.ymax);
    box.zmax = std::min(parent.zmax, box.zmax);
}