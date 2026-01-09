/*    
    Copyright 2013-2026 ONERA.

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

void Mesh_refine(Mesh *M, std::vector<E_Int> &ref_cells,
    std::vector<E_Int> &ref_faces, std::set<UEdge> &ref_edges)
{
    if (M->mode_2D) Mesh_refine_dir(M, ref_cells, ref_faces, ref_edges);
    else Mesh_refine_iso(M, ref_cells, ref_faces, ref_edges);
}

void Mesh_get_ref_entities(Mesh *M, std::vector<E_Int> &ref_cells,
    std::vector<E_Int> &ref_faces, std::set<UEdge> &ref_edges)
{
    assert(ref_cells.empty());
    assert(ref_faces.empty());
    assert(ref_edges.empty());

    for (E_Int i = 0; i < M->nc; i++) {
        if (M->cref[i] > 0) ref_cells.push_back(i);
    }

    for (E_Int i = 0; i < M->nf; i++) {
        if (M->fref[i] > 0) ref_faces.push_back(i);
    }

    for (E_Int fid : ref_faces) {
        E_Int *face = Mesh_get_face(M, fid);
        E_Int *frange = Mesh_get_frange(M, fid);
        E_Int size = 2 * M->fstride[fid];
        assert(size == 8);

        for (E_Int i = 0; i < size; i += 2) {
            if (frange[i/2] == 2) continue;

            E_Int p = face[i];
            assert(face[i+1] == -1);
            E_Int q = face[(i+2)%size];

            UEdge E(p, q);
            auto it = ref_edges.find(E);

            if (it == ref_edges.end()) {
                ref_edges.insert(E);
            }
        }
    }
}

void Mesh_resize_face_data(Mesh *M, E_Int new_nf)
{
    M->faces   = (E_Int *)XRESIZE(M->faces,   (8 * new_nf) * sizeof(E_Int));
    M->frange  = (E_Int *)XRESIZE(M->frange,  (4 * new_nf) * sizeof(E_Int));
    M->fstride = (E_Int *)XRESIZE(M->fstride, new_nf       * sizeof(E_Int));

    for (E_Int i = 8*M->nf; i < 8*new_nf; i++) M->faces[i] = -1;
    for (E_Int i = 4*M->nf; i < 4*new_nf; i++) M->frange[i] = -1;
    for (E_Int i = M->nf; i < new_nf; i++) M->fstride[i] = -1;

    M->owner = (E_Int *)XRESIZE(M->owner, new_nf * sizeof(E_Int));
    M->neigh = (E_Int *)XRESIZE(M->neigh, new_nf * sizeof(E_Int));

    for (E_Int i = M->nf; i < new_nf; i++) {
        M->owner[i] = -1;
        M->neigh[i] = -1;
    }

    M->ftype  = (E_Int *)XRESIZE(M->ftype,  new_nf * sizeof(E_Int));
    M->flevel = (E_Int *)XRESIZE(M->flevel, new_nf * sizeof(E_Int));

    M->fref = (E_Int *)XRESIZE(M->fref, new_nf * sizeof(E_Int));
    for (E_Int i = M->nf; i < new_nf; i++) M->fref[i] = 0;

    M->fparent = (E_Int *)XRESIZE(M->fparent, new_nf * sizeof(E_Int));
    for (E_Int i = M->nf; i < new_nf; i++) M->fparent[i] = -1;

    M->ftag = (E_Int *)XRESIZE(M->ftag, new_nf * sizeof(E_Int));
    for (E_Int i = M->nf; i < new_nf; i++) M->ftag[i] = -1;

    M->fpattern = (E_Int *)XRESIZE(M->fpattern, new_nf * sizeof(E_Int));
    for (E_Int i = M->nf; i < new_nf; i++) M->fpattern[i] = -1;
}

void Mesh_resize_cell_data(Mesh *M, E_Int new_nc)
{
    M->cells   = (E_Int *)XRESIZE(M->cells,   (24 * new_nc) * sizeof(E_Int));
    M->crange  = (E_Int *)XRESIZE(M->crange,  (6 * new_nc)  * sizeof(E_Int));
    M->cstride = (E_Int *)XRESIZE(M->cstride, new_nc        * sizeof(E_Int));

    for (E_Int i = 24*M->nc; i < 24*new_nc; i++) M->cells[i] = -1;
    for (E_Int i = 6*M->nc; i < 6*new_nc; i++) M->crange[i] = -1;
    for (E_Int i = M->nc; i < new_nc; i++) M->cstride[i] = -1;

    M->ctype  = (E_Int *)XRESIZE(M->ctype,  new_nc * sizeof(E_Int));
    M->clevel = (E_Int *)XRESIZE(M->clevel, new_nc * sizeof(E_Int));
}

void Mesh_resize_for_refinement(Mesh *M, const std::vector<E_Int> &ref_cells,
    const std::vector<E_Int> &ref_faces, const std::set<UEdge> &ref_edges)
{
    // 7 new cells per refined cells

    E_Int cell_incr = ref_cells.size() * 7;
    
    // 3 new faces per refined face + 12 new faces per refined cell
    
    E_Int face_incr = ref_faces.size() * 3 + ref_cells.size() * 12;
    
    // New points:
    // + 1 per refined edge,
    // + 1 per refined face,
    // + 1 centroid per refined cell.
    
    //E_Int point_incr = 5 * ref_faces.size() + ref_cells.size();
    E_Int point_incr = ref_edges.size() + ref_faces.size() + ref_cells.size();

    E_Int new_nc = M->nc + cell_incr;
    E_Int new_nf = M->nf + face_incr;
    E_Int new_np = M->np + point_incr;

    M->X = (E_Float *)XRESIZE(M->X, new_np * sizeof(E_Float));
    M->Y = (E_Float *)XRESIZE(M->Y, new_np * sizeof(E_Float));
    M->Z = (E_Float *)XRESIZE(M->Z, new_np * sizeof(E_Float));

    Mesh_resize_face_data(M, new_nf);
    Mesh_resize_cell_data(M, new_nc);
}

void Mesh_sort_ref_entities_by_level(Mesh *M,
    std::vector<E_Int> &ref_cells, std::vector<E_Int> &ref_faces,
    std::set<UEdge> &ref_edges)
{
    std::sort(ref_faces.begin(), ref_faces.end(),
        [&] (E_Int i, E_Int j) { return M->flevel[i] < M->flevel[j]; });
    
    std::sort(ref_cells.begin(), ref_cells.end(),
        [&] (E_Int i, E_Int j) { return M->clevel[i] < M->clevel[j]; });
}
