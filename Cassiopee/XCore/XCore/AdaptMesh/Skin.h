#pragma once

#include "xcore.h"
#include "Vec.h"

struct SkinGraph {
    E_Int nf;
    E_Int *skin;
    E_Int *xadj;
    E_Int *fpts;
    E_Int *fnei;
    Vec3f *fnml;
};

void SkinGraph_free(SkinGraph *skin_graph);
