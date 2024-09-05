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
#pragma once

#include <vector>
#include <cstddef>
#include <map>

#include "common/common.h"
#include "queue.h"
#include "point.h"

struct Vertex;
struct Hedge;
struct Face;
struct Segment;
struct Smesh;
struct Cycle;

struct Dcel {
    std::vector<Vertex *> V;
    std::vector<Hedge *> H;
    std::vector<Face *> F;
    std::vector<Cycle *> C;

    Queue Q; // Filters out duplicate vertices

    Face *f_unbounded[2];

    static E_Int RED;
    static E_Int BLACK;
    static E_Int NO_IDEA;

    std::map<Vertex *, std::vector<Hedge *>> Up;
    std::map<Vertex *, std::vector<Hedge *>> Cp;
    std::map<Vertex *, std::vector<Hedge *>> Lp;

    Dcel(Smesh &M0, Smesh &M1);

    ~Dcel();
    
    void init_vertices(const Smesh &M0, const Smesh &M1);

    void init_hedges_and_faces(Smesh &M, E_Int color);

    static E_Int check_hedges(const std::vector<Hedge *> &H);

    static E_Int check_faces(const std::vector<Hedge *> &H,
        const std::vector<Face *> &F);
    
    void make_cycles();

    void set_face_labels(std::vector<Face *> &F);

    Hedge *get_hedge_of_color(Face *f, E_Int color);

    std::vector<Face *> make_cycle_faces(const std::vector<Cycle *> &C);

    void update_hedge_faces(const std::vector<Face *> &F);

    void set_cycles_inout(const Smesh &M, const Smesh &S);

    std::vector<E_Int> extract_indices_of_type(E_Int inout);
    
    std::vector<Face *> extract_faces_of_indices(
        const std::vector<E_Int> &indices);

    void write_ngon(const char *fname, const std::vector<Face *> &faces) const;

    void write_degen_faces(const char *fname);
    
    void write_outer_faces(const char *fname);
    
    void write_inner_faces(const char *fname);

    static std::vector<Vertex *> get_face_vertices(Face *f);

    void locate_spoints(const Smesh &M, const Smesh &S);

    void find_intersections_3D(const Smesh &M, const Smesh &S);

    void cut_hedge_at_vertex(Hedge *h, Vertex *v);

    void resolve_hedges(const Smesh &M, const Smesh &S);

    void reconstruct(const Smesh &M, const Smesh &S);

    E_Int get_next_face(const Smesh &M, E_Float px, E_Float py, E_Float pz,
        const std::vector<E_Int> &pf, E_Float dir[3]);

    void handle_intersecting_endpoint(Vertex *v, const Smesh &M);

    void trace_hedge(Hedge *sh, const Smesh &M, const Smesh &S, E_Int hid);

    void sort_leaving_hedges(std::vector<Hedge *> &leaving, const E_Float N[3],
        const Smesh &M) const;
};
