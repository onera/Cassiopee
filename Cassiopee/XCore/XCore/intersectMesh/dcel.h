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
#include "point.h"
#include "primitives.h"

struct Smesh;

struct Dcel {

    static const int RED = 0;
    static const int BLACK = 1;
    static const int NO_IDEA = -1;

    struct Vertex {
        E_Float x, y, z;
        E_Int oids[2] = {-1, -1};
        E_Int id = -1;

        Vertex(E_Float x_, E_Float y_, E_Float z_)
        {
            x = x_;
            y = y_;
            z = z_;
            oids[0] = -1;
            oids[1] = -1;
        }
    };

    int cmp_vtx(const Vertex *p, const Vertex *q) const
    {
        return cmp_points(p->x, p->y, p->z, q->x, q->y, q->z);
    }

    struct Face;
    struct Cycle;

    struct Hedge {
        Vertex *orig = NULL;
        Hedge *twin = NULL;
        Hedge *next = NULL;
        Hedge *prev = NULL;
        Face *left = NULL;
        int color = NO_IDEA;
        Cycle *cycle = NULL;

        Hedge(Vertex *V, int color_)
        {
            orig = V;
            color = color_;
        }
    };

    struct Face {
        Hedge *rep = NULL;
        E_Int oids[2] = {-1, -1};
    };

    struct Cycle {
        Hedge *rep = NULL;
        int inout = 0;

        static const int HOLE = -1;
        static const int DEGEN = 0;
        static const int INNER = 1;
        static const int OUTER = 2;

        Cycle(Hedge *h)
        : rep(h)
        {}
    };

    std::vector<Vertex *> V;
    std::vector<Hedge *> H;
    std::vector<Face *> F;
    std::vector<Cycle *> C;

    E_Int inner = 0;
    E_Int outer = 0;
    E_Int degen = 0;
    E_Int hole = 0;

    E_Int dup_x = 0; // Number of duplicate intersections
    std::set<Vertex *> vertices_crossed; // M vertices crossed by S hedges
    std::map<Hedge *, std::vector<Vertex *>> hedge_intersections;

    std::map<E_Int, std::vector<E_Int>> grid;

    Dcel(const Smesh &M, E_Int color);

    Dcel(Smesh &M0, Smesh &M1);

    ~Dcel();

    void init_mh_sv_intersections(const Smesh &M);
    
    void init_vertices(const Smesh &M0, const Smesh &M1);

    void init_hedges_and_faces(const Smesh &M, E_Int color);

    static E_Int check_hedges(const std::vector<Hedge *> &H);

    static E_Int check_faces(const std::vector<Hedge *> &H,
        const std::vector<Face *> &F);
    
    void make_cycles();

    void set_face_labels(std::vector<Face *> &F);

    Hedge *get_hedge_of_color(Face *f, E_Int color);

    std::vector<Face *> make_cycle_faces(const std::vector<Cycle *> &C);

    void update_hedge_faces(const std::vector<Face *> &F);

    void set_cycles_inout(const Smesh &M, E_Int color);

    

    static std::vector<Vertex *> get_face_vertices(Face *f);

    void locate_spoints(const Smesh &M, const Smesh &S);

    void find_intersections_3D(const Smesh &M, const Smesh &S);

    void cut_hedge_at_vertex(Hedge *h, Vertex *v);

    void resolve_hedges(const Smesh &M, const Smesh &S);

    void reconstruct(const Smesh &M, const Smesh &S);

    E_Int get_next_face(const Smesh &M, E_Float px, E_Float py, E_Float pz,
        const std::vector<E_Int> &pf, E_Float dir[3], E_Int hid);
    
    E_Int trace_hedge_2(Hedge *sh, const Smesh &M, const Smesh &S, E_Int hid);

    void cut_hedges();

    void sort_leaving_hedges(std::vector<Hedge *> &leaving, const E_Float N[3]) const;
    
    Smesh export_smesh(bool check_Euler=true) const;

    // Extract

    std::vector<E_Int> extract_indices_of_type(int inout) const;
    std::vector<Cycle *> extract_cycles_of_indices(
        const std::vector<E_Int> &indices) const;

    // IO

    void write_face(const char *fname, const Face *face) const;
    void write_faces(const char *fname, const std::vector<Face *> &faces,
        E_Float scale = 1.0) const;
    void write_hedge(const char *fname, const Hedge *h) const;
    void write_point(const char *fname, const Vertex *v) const;
    void write_point(const char *fname, const std::vector<Dcel::Vertex *> &I) const;
    void write_ngon(const char *fname, const std::vector<Cycle *> &cycles) const;
    void write_degen_cycles(const char *fname) const;
    void write_inner_cycles(const char *fname) const;
    void write_hole_cycles(const char *fname) const;
    void write_outer_cycles(const char *fname) const;
    void write_cycles_of_type(const char *fname, E_Int type) const;
};
