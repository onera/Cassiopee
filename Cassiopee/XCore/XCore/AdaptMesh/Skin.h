#pragma once

#include "xcore.h"
#include "Vec.h"

struct Mesh;
struct DynMesh;

struct SkinGraph {
    E_Int nf;
    E_Int *skin;
    E_Int *xadj;
    E_Int *fpts;
    E_Int *fnei;
    Vec3f *fnml;
};

void SkinGraph_free(SkinGraph *skin_graph);

void SkinGraph_make_skin_neighbours(SkinGraph *skin_graph);

void SkinGraph_smooth_ref_data(const SkinGraph *skin_graph,
    E_Int *fdat, const Mesh *M);
