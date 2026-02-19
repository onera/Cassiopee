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
#pragma once

#include <vector>
#include <cstddef>
#include <map>

#include "common/common.h"
#include "point.h"
#include "primitives.h"
#include "u_edge.h"

struct Smesh;

struct Dcel {

    static const int RED = 0;
    static const int BLACK = 1;
    static const int NO_IDEA = -1;

    struct Vertex {
        E_Float x, y, z;
        E_Int oids[2] = {-1, -1};
        E_Int id = -1;

        PointLoc ploc;
        E_Float d2;

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

    struct cmp_vertex {
        bool operator()(const Vertex *p, const Vertex *q) const
        {
            return cmp_points(p->x, p->y, p->z, q->x, q->y, q->z) < 0;
        }
    };

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

    std::vector<std::vector<Vertex *>> Fv;

    std::set<Vertex *, cmp_vertex> vertex_set;

    E_Int inner = 0;
    E_Int outer = 0;
    E_Int degen = 0;
    E_Int hole = 0;

    std::map<Hedge *, std::vector<Vertex *>> hedge_intersections;
    std::map<std::pair<Vertex *, Vertex *>, Vertex *> vcenter[2];

    std::map<E_Int, E_Int> l2gp;
    std::map<E_Int, E_Int> g2lp;
    std::map<E_Int, E_Int> l2gf;
    std::map<E_Int, E_Int> g2lf;

    Dcel();
    Dcel(const Smesh &Mf, const Smesh &Sf, const std::vector<PointLoc> &plocs);
    Dcel(const Smesh &Mf, int color);
    ~Dcel();
    
    // Intersection

    void init_vertices(const Smesh &Mf, const Smesh &Sf,
        const std::vector<PointLoc> &plocs);
    void init_hedges_and_faces(const Smesh &Mf, int color);
    void sort_leaving_hedges(std::vector<Hedge *> &leaving,
        const E_Float N[3]) const;
    void make_cycles();
    void set_cycles_inout();
    void set_face_labels(std::vector<Face *> &F);
    std::vector<Face *> make_cycle_faces(const std::vector<Cycle *> &C);
    void reconstruct(const Smesh &Mf, const Smesh &Sf);
    static Dcel intersect(const Smesh &Mf, const Smesh &Sf,
        const std::vector<PointLoc> &plocs);
    void triangulate(const Smesh &Mf, const Smesh &Sf);

    // Checks

    static E_Int check_hedges(const std::vector<Hedge *> &H);
    static E_Int check_faces(const std::vector<Hedge *> &H,
        const std::vector<Face *> &F);
    
    // Helpers

    Hedge *get_hedge_of_color(Face *f, int color);
    void update_hedge_faces(std::vector<Face *> &F);
    void get_face_vertices(const Face *f, std::vector<Vertex *> &vertices);
    //void get_vertex_normal(Vertex *q, const Smesh &Mf, E_Float N[3]);
    //bool is_vertex_in_triangle(Vertex *v, Vertex *a, Vertex *b, Vertex *c);
    //bool vertex_list_is_convex(const Vertex *a, const Vertex *b,
    //    const Vertex *c, const Smesh &Mf);

    // Export

    Smesh export_smesh(bool check_Euler=true) const;
    void reconstruct(Smesh &Mf, int color) const;

    // Extract

    std::vector<E_Int> extract_indices_of_type(int inout) const;
    std::vector<Cycle *> extract_cycles_of_indices(
        const std::vector<E_Int> &indices) const;

    // IO

    void write_face(const char *fname, const Face *face) const;
    void write_faces(const char *fname, const std::vector<Face *> &faces,
        E_Float scale = 1.0) const;
    void write_hedge(const char *fname, const Hedge *h) const;
    void write_vertex(const char *fname, const Vertex *v) const;
    void write_point(const char *fname, const std::vector<Vertex *> &I) const;
    void write_ngon(const char *fname, const std::vector<Cycle *> &cycles) const;
    void write_ngon(const char *fname, const std::vector<Face *> &faces) const;
    void write_ngon(const char *fname, const std::vector<E_Int> &fids) const;
    void write_degen_cycles(const char *fname) const;
    void write_inner_cycles(const char *fname) const;
    void write_hole_cycles(const char *fname) const;
    void write_outer_cycles(const char *fname) const;
    void write_cycles_of_type(const char *fname, int type) const;
};
