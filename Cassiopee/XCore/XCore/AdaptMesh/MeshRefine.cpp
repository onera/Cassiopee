/*    
    Copyright 2013-2024 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "Mesh.h"
#include "common/mem.h"
#include "Hexa.h"
#include "Quad.h"
#include "Edge.h"

void Mesh_refine(Mesh *M, std::vector<Int> &ref_cells,
    std::vector<Int> &ref_faces, std::set<UEdge> &ref_edges)
{
    if (M->mode_2D) Mesh_refine_dir(M, ref_cells, ref_faces, ref_edges);
    else Mesh_refine_iso(M, ref_cells, ref_faces, ref_edges);
}

void Mesh_get_ref_entities(Mesh *M, std::vector<Int> &ref_cells,
    std::vector<Int> &ref_faces, std::set<UEdge> &ref_edges)
{
    assert(ref_cells.empty());
    assert(ref_faces.empty());
    assert(ref_edges.empty());

    for (Int i = 0; i < M->nc; i++) {
        if (M->cref[i] > 0) ref_cells.push_back(i);
    }

    for (Int i = 0; i < M->nf; i++) {
        if (M->fref[i] > 0) ref_faces.push_back(i);
    }

    for (Int fid : ref_faces) {
        Int *face = Mesh_get_face(M, fid);
        Int *frange = Mesh_get_frange(M, fid);

        for (Int i = 0; i < M->fstride[fid]; i++) {
            if (frange[i] == 2) continue;

            Int *pn = face + 2*i;

            Int p = pn[0];
            assert(pn[1] == -1);
            Int q = pn[2];

            UEdge E(p, q);
            auto it = ref_edges.find(E);

            if (it == ref_edges.end()) {
                ref_edges.insert(E);
            }
        }
    }
}

void Mesh_resize_for_refinement(Mesh *M, const std::vector<Int> &ref_cells,
    const std::vector<Int> &ref_faces, const std::set<UEdge> &ref_edges)
{
    // 7 new cells per refined cells

    Int cell_incr = ref_cells.size() * 7;
    
    // 3 new faces per refined faces + 13 new faces per refined cell
    
    Int face_incr = ref_faces.size() * 3 + ref_cells.size() * 12;
    
    // New points:
    // + 1 per refined edge,
    // + 1 per refined face,
    // + 1 centroid per refined cell.
    
    //Int point_incr = 5 * ref_faces.size() + ref_cells.size();
    Int point_incr = ref_edges.size() + ref_faces.size() + ref_cells.size();

    Int new_nc = M->nc + cell_incr;
    Int new_nf = M->nf + face_incr;
    Int new_np = M->np + point_incr;

    //printf("%d -> cells before: %d - cells after: %d\n", M->pid, M->nc, new_nc);

    M->X = (Float *)XRESIZE(M->X, new_np * sizeof(Float));
    M->Y = (Float *)XRESIZE(M->Y, new_np * sizeof(Float));
    M->Z = (Float *)XRESIZE(M->Z, new_np * sizeof(Float));

    M->faces   = (Int *)XRESIZE(M->faces,   (8 * new_nf) * sizeof(Int));
    M->frange  = (Int *)XRESIZE(M->frange,  (4 * new_nf) * sizeof(Int));
    M->fstride = (Int *)XRESIZE(M->fstride, new_nf       * sizeof(Int));

    for (Int i = 8*M->nf; i < 8*new_nf; i++) M->faces[i] = -1;
    for (Int i = 4*M->nf; i < 4*new_nf; i++) M->frange[i] = -1;
    for (Int i = M->nf; i < new_nf; i++) M->fstride[i] = -1;

    M->cells   = (Int *)XRESIZE(M->cells,   (24 * new_nc) * sizeof(Int));
    M->crange  = (Int *)XRESIZE(M->crange,  (6 * new_nc)  * sizeof(Int));
    M->cstride = (Int *)XRESIZE(M->cstride, new_nc        * sizeof(Int));

    for (Int i = 24*M->nc; i < 24*new_nc; i++) M->cells[i] = -1;
    for (Int i = 6*M->nc; i < 6*new_nc; i++) M->crange[i] = -1;
    for (Int i = M->nc; i < new_nc; i++) M->cstride[i] = -1;

    M->owner = (Int *)XRESIZE(M->owner, new_nf * sizeof(Int));
    M->neigh = (Int *)XRESIZE(M->neigh, new_nf * sizeof(Int));

    for (Int i = M->nf; i < new_nf; i++) {
        M->owner[i] = -1;
        M->neigh[i] = -1;
    }

    M->ftype  = (Int *)XRESIZE(M->ftype,  new_nf * sizeof(Int));
    M->flevel = (Int *)XRESIZE(M->flevel, new_nf * sizeof(Int));

    M->ctype  = (Int *)XRESIZE(M->ctype,  new_nc * sizeof(Int));
    M->clevel = (Int *)XRESIZE(M->clevel, new_nc * sizeof(Int));

    M->fref = (Int *)XRESIZE(M->fref, new_nf * sizeof(Int));
    for (Int i = M->nf; i < new_nf; i++) M->fref[i] = 0;

    M->fparent = (Int *)XRESIZE(M->fparent, new_nf * sizeof(Int));
    for (Int i = M->nf; i < new_nf; i++) M->fparent[i] = -1;
}

void Mesh_sort_ref_entities_by_level(Mesh *M,
    std::vector<Int> &ref_cells, std::vector<Int> &ref_faces,
    std::set<UEdge> &ref_edges)
{
    std::sort(ref_faces.begin(), ref_faces.end(),
        [&] (Int i, Int j) { return M->flevel[i] < M->flevel[j]; });
    
    std::sort(ref_cells.begin(), ref_cells.end(),
        [&] (Int i, Int j) { return M->clevel[i] < M->clevel[j]; });
}