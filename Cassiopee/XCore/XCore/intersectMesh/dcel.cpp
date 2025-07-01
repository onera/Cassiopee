/*    
    Copyright 2013-2025 Onera.

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
#include <cstdio>
#include <cassert>
#include <cmath>
#include <algorithm>

#include "dcel.h"
#include "primitives.h"
#include "io.h"
#include "smesh.h"
#include "triangle.h"

Dcel::Dcel() 
{}

void Dcel::init_vertices(const Smesh &Mf, const Smesh &Sf,
    const std::vector<PointLoc> &plocs)
{
    E_Int duplicate_vertices = 0;

    for (E_Int i = 0; i < Mf.np; i++) {
        Vertex *v = new Vertex(Mf.X[i], Mf.Y[i], Mf.Z[i]);
        v->oids[Dcel::RED] = i;
        vertex_set.insert(v);
    }

    for (E_Int i = 0; i < Sf.np; i++) {
        Vertex tmp(Sf.X[i], Sf.Y[i], Sf.Z[i]);
        auto it = vertex_set.find(&tmp);
        Vertex *v = NULL;
        if (it == vertex_set.end()) {
            v = new Vertex(Sf.X[i], Sf.Y[i], Sf.Z[i]);
            v->oids[Dcel::BLACK] = i;
            vertex_set.insert(v);
        } else {
            v = *it;
            v->oids[Dcel::BLACK] = i;
            duplicate_vertices++;
        }
        v->ploc = plocs[i];
    }

    V.reserve(vertex_set.size());
    for (Vertex *v : vertex_set) {
        v->id = V.size();
        V.push_back(v);
    }

    printf("Duplicate vertices: %d\n", duplicate_vertices);
}

void Dcel::get_face_vertices(const Face *f, std::vector<Vertex *> &ret) 
{
    assert(f);
    ret.clear();
    Hedge *h = f->rep;
    ret.push_back(h->orig);
    Hedge *w = h->next;
    while (w != h) {
        ret.push_back(w->orig);
        w = w->next;
    }
}

struct HitData {
    E_Float t, s;
};

Dcel Dcel::intersect(const Smesh &Mf, const Smesh &Sf,
    const std::vector<PointLoc> &plocs)
{
    Dcel D;
    D.init_vertices(Mf, Sf, plocs);

    D.init_hedges_and_faces(Mf, RED);
    D.init_hedges_and_faces(Sf, BLACK);

    if (D.check_hedges(D.H) != 0) {
        fprintf(stderr, "Dcel: Inconsistent half-edge records!\n");
        abort();
    }
    
    D.make_cycles();
    D.set_cycles_inout();

    if (check_faces(D.H, D.F) != 0) {
        fprintf(stderr, "Dcel: Inconsistent face records!\n");
        abort();
    }

    //std::vector<Point> xpoints;
    //std::vector<Point> dpoints;

    // Register the intersections between Sf points and Mf edges
    E_Int v_on_e = 0;
    for (Vertex *v : D.V) {
        // Strictly the Sf points
        if (v->oids[0] == -1 && v->oids[1] != -1) {
            const auto &ploc = v->ploc;
            if (ploc.e_idx == -1) continue;

            // Set the reference edge for normal computation
            E_Int fid = ploc.fid;
            const auto &pe = Mf.F2E[fid];
            E_Int e = pe[ploc.e_idx];

            // Add the vertex to hedge intersections
            Hedge *h = D.H[2*e];
            assert(h->color == RED);
            Vertex *O = h->orig;
            Vertex *T = h->twin->orig;
            if (D.cmp_vtx(O, T) > 0) h = h->twin;
            D.hedge_intersections[h].push_back(v);

            //dpoints.push_back({v->x,v->y,v->z});

            v_on_e++;
        }
    }
    printf("Sf points on Mf edges: %d\n", v_on_e);
    //point_write("dpoints.im", dpoints);

    // Trace
    for (E_Int eid_s = 0; eid_s < Sf.ne; eid_s++) {
        Hedge *hs = D.H[2*(Mf.ne + eid_s)];
        assert(hs->color == BLACK);
        Vertex *O = hs->orig;
        Vertex *T = hs->twin->orig;
        if (D.cmp_vtx(O, T) > 0) hs = hs->twin;

        E_Int p = hs->orig->oids[1];
        E_Int q = hs->twin->orig->oids[1];

        E_Float spx = Sf.X[p], spy = Sf.Y[p], spz = Sf.Z[p];
        E_Float sqx = Sf.X[q], sqy = Sf.Y[q], sqz = Sf.Z[q];

        E_Float DIR[3] = {sqx-spx, sqy-spy, sqz-spz};
        E_Float NORM = K_MATH::norm(DIR, 3);
        DIR[0] /= NORM, DIR[1] /= NORM, DIR[2] /= NORM;

        std::vector<E_Int> orig_faces;
        std::vector<E_Int> tail_faces;

        E_Int last_vertex = -1, last_edge = -1, dummy;

        Mf.get_shared_faces(plocs[p], orig_faces, last_vertex, last_edge); 
        Mf.get_shared_faces(plocs[q], tail_faces, dummy, dummy);

        E_Int starting_face = Mf.deduce_face(orig_faces, spx, spy, spz,
            DIR, last_vertex, last_edge, eid_s);
        assert(starting_face != -1);

        bool found_tail = false;
        E_Int cur_fid = starting_face;
        E_Float cur_pos[3] = {spx, spy, spz};

        E_Int walk = 0;
        E_Int max_walks = 20;

        std::vector<E_Int> path;

        while (!found_tail && walk <= max_walks) {

            path.push_back(cur_fid);

            for (auto fid : tail_faces) {
                if (fid == cur_fid) {
                    found_tail = true;
                    break;
                }
            }

            if (found_tail) break;

            // Update the direction
            DIR[0] = sqx-cur_pos[0];
            DIR[1] = sqy-cur_pos[1];
            DIR[2] = sqz-cur_pos[2];
            E_Float NORM = K_MATH::norm(DIR, 3);
            DIR[0] /= NORM, DIR[1] /= NORM, DIR[2] /= NORM;

            // Project
            E_Float proj[3];
            Mf.get_unit_projected_direction(cur_fid, DIR, proj);

            const auto &pn = Mf.Fc[cur_fid];
            const auto &pe = Mf.F2E[cur_fid];
            assert(pe.size() == pn.size());
            //const E_Float *fN = &Mf.fnormals[3*cur_fid];

            E_Int next_fid = -1;
            E_Float next_pos[3] = {EFLOATMAX, EFLOATMAX, EFLOATMAX};

            bool hit = false;

            std::vector<HitData> hitData;

            for (size_t i = 0; i < pn.size(); i++) {
                E_Int p = pn[i];
                E_Int q = pn[(i+1)%pn.size()];
                E_Int e = pe[i];

                if (p == last_vertex || q == last_vertex || e == last_edge)
                    continue;
                
                E_Float px = Mf.X[p], py = Mf.Y[p], pz = Mf.Z[p];
                E_Float qx = Mf.X[q], qy = Mf.Y[q], qz = Mf.Z[q];
            
                E_Float t = -1.0, s = -1.0;
                hit = ray_edge_intersect(
                    cur_pos[0], cur_pos[1], cur_pos[2],
                    proj[0], proj[1], proj[2],
                    px, py, pz, qx, qy, qz,
                    t, s
                );

                hitData.push_back({t, s});

                if (hit) {
                    if (s > RAY_EDGE_TOL && s < 1 - RAY_EDGE_TOL) {
                        // Hit edge middle
                        const auto &pe = Mf.F2E[cur_fid];
                        E_Int eid_m = pe[i];
                        last_edge = eid_m;
                        last_vertex = -1;
                        assert(Mf.E2F[eid_m][0] == cur_fid || Mf.E2F[eid_m][1] == cur_fid);
                        if (Mf.E2F[eid_m][0] == cur_fid) next_fid = Mf.E2F[eid_m][1];
                        else next_fid = Mf.E2F[eid_m][0];

                        next_pos[0] = cur_pos[0] + t * proj[0];
                        next_pos[1] = cur_pos[1] + t * proj[1];
                        next_pos[2] = cur_pos[2] + t * proj[2];

                        //if (eid_s == 80) {
                        //    point_write("x.im", next_pos[0], next_pos[1], next_pos[2]);
                        //}

                        // Create a new intersection vertex
                        Vertex tmp(next_pos[0], next_pos[1], next_pos[2]);
                        assert(D.vertex_set.find(&tmp) == D.vertex_set.end());
                        Vertex *x = new Vertex(next_pos[0], next_pos[1],
                            next_pos[2]);
                        
                        x->ploc.fid = cur_fid;
                        x->ploc.e_idx = i;

                        x->id = D.V.size();
                        D.V.push_back(x);
                        D.vertex_set.insert(x);

                        // Register the intersection
                        Hedge *hm = D.H[2*eid_m];
                        assert(hm->color == RED);
                        Vertex *O = hm->orig;
                        Vertex *T = hm->twin->orig;
                        if (D.cmp_vtx(O, T) > 0) hm = hm->twin;

                        D.hedge_intersections[hm].push_back(x);
                        D.hedge_intersections[hs].push_back(x);

                    } else {
                        // Hit an edge endpoint
                        bool hit_p = (s <= RAY_EDGE_TOL);
                        bool hit_q = (s >= 1 - RAY_EDGE_TOL);
                        assert(!(hit_p && hit_q));
                        last_edge = -1;
                        if (hit_p) last_vertex = p;
                        else last_vertex = q;
                        next_pos[0] = Mf.X[last_vertex];
                        next_pos[1] = Mf.Y[last_vertex];
                        next_pos[2] = Mf.Z[last_vertex];
                        const auto &pf = Mf.P2F[last_vertex];

                        //if (eid_s == 80) {
                        //    Mf.write_ngon("pf.im", pf);
                        //    point_write("x.im", next_pos[0], next_pos[1], next_pos[2]);
                        //}

                        next_fid = Mf.deduce_face(pf,
                            next_pos[0], next_pos[1], next_pos[2],
                            DIR, last_vertex, last_edge, eid_s
                        );
                        assert(next_fid != -1);

                        // Find Vertex corresponding to hit vertex
                        Vertex tmp(Mf.X[last_vertex], Mf.Y[last_vertex], Mf.Z[last_vertex]);
                        auto it = D.vertex_set.find(&tmp);
                        assert(it != D.vertex_set.end());
                        Vertex *x = *it;

                        // Register the intersection
                        D.hedge_intersections[hs].push_back(x);

                        //xpoints.push_back({x->x, x->y, x->z});
                    }
                    break;
                }
            }
            if (!hit) {
                fprintf(stderr, "Failed to hit. Hit data:\n");
                for (const auto &hd : hitData) {
                    printf("t = %.12e | s = %.12e\n", hd.t, hd.s);
                }
                fflush(stdout);
                std::vector<E_Int> fids;
                fids.push_back(cur_fid);
                Mf.write_ngon("cur_fid.im", fids);
                edge_write("cur_sedge.im", cur_pos[0], cur_pos[1], cur_pos[2],
                    proj[0], proj[1], proj[2]);
            }
            assert(hit);
            assert(next_fid != cur_fid);
            cur_fid = next_fid;
            cur_pos[0] = next_pos[0];
            cur_pos[1] = next_pos[1];
            cur_pos[2] = next_pos[2];
            walk++;
        }

        if (!found_tail) {
            edge_write("lost_edge.im", spx, spy, spz, sqx, sqy, sqz);
            Mf.write_ngon("orig_faces.im", orig_faces);
            Mf.write_ngon("tail_faces.im", tail_faces);
            Mf.write_ngon("path.im", path);
        }

        assert(found_tail);
        assert(walk <= max_walks);
    }

    //point_write("xpoints.im", xpoints);

    // Cut
    for (auto &h2x : D.hedge_intersections) {
        Hedge *h = h2x.first;
        auto it = D.hedge_intersections.find(h->twin);
        assert(it == D.hedge_intersections.end());

        auto &xs = h2x.second;

        Vertex *o = h->orig;
        Vertex *tail = h->twin->orig;
        assert(D.cmp_vtx(o, tail) < 0);

        // TODO(Imad): check that the intersections are 'sufficiently' spaced out
        for (Vertex *x : xs) {
            E_Float d[3] = {x->x-o->x, x->y-o->y, x->z-o->z};
            x->d2 = K_MATH::dot(d, d, 3);
        }

        std::sort(xs.begin(), xs.end(), [&] (const Vertex *a, const Vertex *b)
        {
            assert(Sign(a->d2-b->d2) != 0);
            return a->d2 < b->d2;
        });

        // Before cutting, cache the intersections for volume mesh reconstruction
        {
            Vertex *start = o;
            Vertex *end = tail;
            for (size_t i = 0; i < xs.size(); i++) {
                D.vcenter[h->color][{start, end}] = xs[i];
                start = xs[i];
            }
        }

        // Cut
        Hedge *current_h = h;
        Hedge *t = h->twin;

        for (Vertex *x : xs) {
            Hedge *e1 = new Hedge(x, current_h->color);
            Hedge *e2 = new Hedge(x, t->color);

            D.H.push_back(e1);
            D.H.push_back(e2);

            e1->left = current_h->left;
            e2->left = t->left;

            current_h->twin = e2;
            e2->twin = current_h;
            t->twin = e1;
            e1->twin = t;

            current_h->next = e1;
            e1->prev = current_h;
            t->next = e2;
            e2->prev = t;

            current_h = e1;
        }
    }

    // Resolve
    std::vector<std::vector<Hedge *>> list(D.V.size());
    for (Hedge *h : D.H) {
        Vertex *o = h->orig;
        list[o->id].push_back(h);
    }

    for (size_t vid = 0; vid < D.V.size(); vid++) {
        Vertex *v = D.V[vid];

        E_Float N[3] = {0, 0, 0};
        
        if (v->oids[0] != -1) {
            const E_Float *pN = &Mf.pnormals[3*v->oids[0]];
            for (E_Int i = 0; i < 3; i++) N[i] = pN[i];
        } else if (v->oids[1] != -1) {
            const E_Float *fN = &Mf.fnormals[3*v->ploc.fid];
            for (E_Int i = 0; i < 3; i++) N[i] = fN[i];
        } else {
            E_Int fid_m = v->ploc.fid;
            E_Int eid_m = Mf.F2E[fid_m][v->ploc.e_idx];
            const auto &pf = Mf.E2F[eid_m];
            assert(pf[0] == fid_m || pf[1] == fid_m);
            const E_Float *fN1 = &Mf.fnormals[3*pf[0]];
            const E_Float *fN2 = &Mf.fnormals[3*pf[1]];
            for (int i = 0; i < 3; i++) N[i] += fN1[i] + fN2[i];
            E_Float NORM = K_MATH::norm(N, 3);
            for (int i = 0; i < 3; i++) N[i] /= NORM;
        }

        auto &leaving = list[vid];

        D.sort_leaving_hedges(leaving, N);

        for (size_t i = 0; i < leaving.size(); i++) {
            Hedge *h = leaving[i];
            Hedge *w = leaving[(i+1)%leaving.size()];
            h->twin->next = w;
            w->prev = h->twin;
        }
    }

    if (D.check_hedges(D.H) != 0) {
        fprintf(stderr, "Dcel: Inconsistent half-edge records!\n");
        assert(0);
        abort();
    }

    D.make_cycles();
    D.set_cycles_inout();

    //write_hole_cycles("hole.im");
    //write_degen_cycles("degen.im");
    //write_inner_cycles("inner.im");

    auto new_F = D.make_cycle_faces(D.C);
    D.set_face_labels(new_F);
    D.update_hedge_faces(new_F);
    D.check_faces(D.H, new_F);

    for (Face *f : D.F) delete f;
    D.F = new_F;

    D.Fv.resize(D.F.size());
    for (size_t fid = 0; fid < D.F.size(); fid++) {
        D.get_face_vertices(D.F[fid], D.Fv[fid]);
    }

    //D.triangulate(Mf, Sf);

    return D;
}

static
void get_vertex_normal(const Dcel::Vertex *q, const Smesh &Mf, E_Float N[3])
{
    N[0] = N[1] = N[2] = 0;
    if (q->oids[0] != -1) {
        const E_Float *pN = &Mf.pnormals[3*q->oids[0]];
        for (E_Int i = 0; i < 3; i++) N[i] = pN[i];
    } else if (q->oids[1] != -1) {
        const E_Float *fN = &Mf.fnormals[3*q->ploc.fid];
        for (E_Int i = 0; i < 3; i++) N[i] = fN[i];
    } else {
        E_Int fid_m = q->ploc.fid;
        E_Int eid_m = Mf.F2E[fid_m][q->ploc.e_idx];
        const auto &pf = Mf.E2F[eid_m];
        assert(pf[0] == fid_m || pf[1] == fid_m);
        const E_Float *fN1 = &Mf.fnormals[3*pf[0]];
        const E_Float *fN2 = &Mf.fnormals[3*pf[1]];
        for (int i = 0; i < 3; i++) N[i] += fN1[i] + fN2[i];
        E_Float NORM = K_MATH::norm(N, 3);
        for (int i = 0; i < 3; i++) N[i] /= NORM;
    }
}

// Circular doubly linked list

struct VNode {
    Dcel::Vertex *v;
    VNode *next;
    VNode *prev;

    VNode(Dcel::Vertex *V)
    {
        v = V;
        next = this;
        prev = this;
    }
};

void VNode_push_front(VNode **head, Dcel::Vertex *v)
{
    VNode *node = new VNode(v);
    if (*head == NULL) {
        *head = node;
    } else {
        VNode *last = (*head)->prev;
        node->next = *head;
        node->prev = last;
        last->next = node;
        (*head)->prev = node;
        *head = node;
    }
}

void VNode_push_back(VNode **head, Dcel::Vertex *v)
{
    VNode *node = new VNode(v);
    if (*head == NULL) {
        *head = node;
    } else {
        VNode *last = (*head)->prev;
        node->next = *head;
        node->prev = last;
        last->next = node;
        (*head)->prev = node;
    }
}

bool VNode_erase(VNode **head, Dcel::Vertex *v)
{
    if (*head == NULL) return false;

    VNode *current = *head;
    do {
        if (current->v == v) {
            if (current->next == current) {
                // Single node in the list
                *head = NULL;
            } else {
                VNode *prev = current->prev;
                VNode *next = current->next;
                prev->next = next;
                next->prev = prev;
                if (current == *head) {
                    // If head to be deleted, set it to the next element
                    *head = next;
                }
            }
            delete current;
            return true;
        }
        current = current->next;
    } while (current != *head);

    assert(0);
    return false;
}

VNode *VNode_find(const VNode *head, const Dcel::Vertex *v)
{
    if (!head) return NULL;

    VNode *current = (VNode *)head;
    do {
        if (current->v == v) return current;
        current = current->next;
    } while (current != head);
    return NULL;
}

static
bool vertex_is_in_triangle(const Dcel::Vertex *v, const Dcel::Vertex *a,
    const Dcel::Vertex *b, const Dcel::Vertex *c)
{
    return Triangle::is_point_inside(v->x, v->y, v->z,
        a->x, a->y, a->z, b->x, b->y, b->z, c->x, c->y, c->z);
}

static
bool vertex_is_ear(const Dcel::Vertex *b, const VNode *polygon,
    const VNode *convex, const VNode *reflex)
{
    // Polygon empty
    if (!polygon) return false;

    // Vertex not in polygon
    VNode *node = VNode_find(polygon, b);
    if (!node) return false;

    // Vertex not convex
    if (!VNode_find(convex, b)) return false;

    // No reflex vertices
    if (!reflex) return true;
    
    const Dcel::Vertex *a = node->prev->v;
    const Dcel::Vertex *c = node->next->v;

    // Test the inclusion of all reflex vertices within triangle {a, b, c}
    VNode *current = (VNode *)reflex;
    bool is_ear = true;

    do {
        Dcel::Vertex *v = current->v;

        if (v != a && v != b && v != c &&
            vertex_is_in_triangle(v, a, b, c)) {
            is_ear = false;
            break;
        }

        current = current->next;
        
    } while (current != reflex);

    return is_ear;
}

void VNode_free_list(VNode *head)
{
    if (!head) return;
    VNode *current = head;
    do {
        VNode *next = current->next;
        delete current;
        current = next;
    } while (current != head);
}

void VNode_print_list(const VNode *head)
{
    if (!head) return;

    VNode *current = (VNode *)head;

    do {
        printf("%d ", current->v->id);
        current = current->next;
    } while (current != head);

    printf("\n");
    fflush(stdout);
}

struct Vertex_triple
{
    Dcel::Vertex *a, *b, *c;
};

static
bool vertex_list_is_convex(const Dcel::Vertex *p, const Dcel::Vertex *q,
    const Dcel::Vertex *r, const Smesh &Mf)
{
    E_Float A[3] = {q->x-p->x, q->y-p->y, q->z-p->z};
    E_Float B[3] = {r->x-q->x, r->y-q->y, r->z-q->z};
    E_Float C[3];
    K_MATH::cross(A, B, C);

    E_Float N[3];
    get_vertex_normal(q, Mf, N);

    E_Float dp = K_MATH::dot(C, N, 3);
    if (dp > TOL) return true;
    return false;
}

static
bool vertex_is_convex(const Dcel::Vertex *b, const VNode *polygon,
    const Smesh &Mf)
{
    VNode *node = VNode_find(polygon, b);
    if (!node) return false;

    const Dcel::Vertex *a = node->prev->v;
    const Dcel::Vertex *c = node->next->v;

    return vertex_list_is_convex(a, b, c, Mf);
}

void Dcel::triangulate(const Smesh &Mf, const Smesh &Sf)
{
    E_Int non_convex_count = 0;
    std::vector<E_Int> non_convex_faces;

    for (size_t fid = 0; fid < F.size(); fid++) {
        // TODO(Imad): skip single color faces
        
        const auto &vertices = Fv[fid];
        if (vertices.size() == 3) continue;

        for (size_t i = 0; i < vertices.size(); i++) {
            Vertex *p = vertices[i];
            Vertex *q = vertices[(i+1)%vertices.size()];
            Vertex *r = vertices[(i+2)%vertices.size()];
            
            E_Float A[3] = {q->x-p->x, q->y-p->y, q->z-p->z};
            E_Float B[3] = {r->x-q->x, r->y-q->y, r->z-q->z};
            E_Float C[3];
            K_MATH::cross(A, B, C);

            E_Float N[3];
            get_vertex_normal(q, Mf, N);

            E_Float dp = K_MATH::dot(C, N, 3);
            if (dp < 0) {
                /*
                write_vertex("p.im", p);
                write_vertex("q.im", q);
                write_vertex("r.im", r);
                */

                std::vector<E_Int> face;
                face.push_back(fid);
                char fname[16] = {0};
                sprintf(fname, "fid%d.im", non_convex_count);
                write_ngon(fname, face);

                non_convex_faces.push_back(fid);
                non_convex_count++;
                break;
            }
        }
    }

    printf("Total faces: %zu\n", F.size());
    printf("Non-convex count: %d\n", non_convex_count);
    write_ngon("non_convex.im", non_convex_faces);

    std::vector<Face *> new_faces;

    for (size_t i = 0; i < non_convex_faces.size(); i++) {
        E_Int fid = non_convex_faces[i];

        const auto &vertices = Fv[fid];
        assert(vertices.size() > 3);

        // Store the polygon

        VNode *polygon = NULL;
        for (Vertex *v : vertices) VNode_push_back(&polygon, v);

        {
            VNode *current = polygon;
            E_Int vid = 0;
            do {
                Vertex *v = current->v;
                char fname[128] = {0};
                sprintf(fname, "vertex%d.im", vid);
                point_write(fname, v->x, v->y, v->z);
                current = current->next;
                vid++;
            } while (current != polygon);
        }

        // Find the convex/reflex vertices

        VNode *convex = NULL, *reflex = NULL;
        VNode *current = polygon;
        do {
            if (vertex_is_convex(current->v, polygon, Mf)) {
                VNode_push_back(&convex, current->v);
            } else {
                VNode_push_back(&reflex, current->v);
            }

            current = current->next;
        } while (current != polygon);

        // Store the ears

        VNode *ears = NULL;
        assert(current == polygon);
        do {
            if (vertex_is_ear(current->v, polygon, convex, reflex))
                VNode_push_back(&ears, current->v);
            current = current->next;
        } while (current != polygon);

        // Ear-clipping algorithm

        size_t polygon_size = vertices.size();

        std::vector<Vertex_triple> tris;

        if (i == 0) {
            point_write("polygon_head.im", polygon->v->x, polygon->v->y, polygon->v->z);
            printf("polygon before: ");
            VNode_print_list(polygon);
            printf("ears before: ");
            VNode_print_list(ears);
            printf("convex before: ");
            VNode_print_list(convex);
            printf("reflex before: ");
            VNode_print_list(reflex);
        }

        while (polygon_size != 3) {
            // Current ear is one of the resulting triangles
            Vertex *b = ears->v;
            VNode *node = VNode_find(polygon, b);
            Vertex *a = node->prev->v;
            Vertex *c = node->next->v;
            tris.push_back({a, b, c});

            if (i == 0) {
                point_write("a.im", a->x, a->y, a->z);
                point_write("b.im", b->x, b->y, b->z);
                point_write("c.im", c->x, c->y, c->z);
            }

            // Delete current ear tip from ear tip list
            VNode_erase(&ears, b);

            // Delete current ear tip from polygon
            VNode_erase(&polygon, b);
            polygon_size--;

            // Delete current ear tip from convex list
            VNode_erase(&convex, b);

            // Rules after ear tip deletion:
            // - if an adjacent vertex was convex, it remains convex, and may become an ear.
            // - if an adjacent vertex was an ear, it does not necessarily remains an ear.
            // - if an adjacent vertex was reflex, it may become convex and possibly and ear.

            // Update prev

            bool was_convex = (VNode_find(convex, a) != NULL);
            if (was_convex) {
                if (!VNode_find(ears, a)) {
                    if (vertex_is_ear(a, polygon, convex, reflex)) {
                        VNode_push_back(&ears, a);
                    }
                }
            } else {
                assert(VNode_find(reflex, a));
                if (vertex_is_convex(a, polygon, Mf)) {
                    VNode_erase(&reflex, a);
                    VNode_push_back(&convex, a);

                    assert(!VNode_find(ears, a));
                    if (vertex_is_ear(a, polygon, convex, reflex)) {
                        VNode_push_back(&ears, a);
                    }
                }
            }

            // Update next

            was_convex = (VNode_find(convex, c) != NULL);
            if (was_convex) {
                if (!VNode_find(ears, c)) {
                    if (vertex_is_ear(c, polygon, convex, reflex)) {
                        VNode_push_back(&ears, c);
                    }
                }
            } else {
                assert(VNode_find(reflex, c));
                if (vertex_is_convex(c, polygon, Mf)) {
                    VNode_erase(&reflex, c);
                    VNode_push_back(&convex, c);

                    assert(!VNode_find(ears, c));
                    if (vertex_is_ear(c, polygon, convex, reflex)) {
                        VNode_push_back(&ears, c);
                    }
                }
            }
        }

        if (i == 0) {
            printf("polygon after: ");
            VNode_print_list(polygon);
            printf("ears after: ");
            VNode_print_list(ears);
            printf("convex after: ");
            VNode_print_list(convex);
            printf("reflex after: ");
            VNode_print_list(reflex);
            puts("");
        }

        tris.push_back({polygon->prev->v, polygon->v, polygon->next->v});

        assert(tris.size() == vertices.size()-2);

        // From the triangles, create new face records.
        // These face records inherit the color of the parent face.
        if (i == 0) {
            for (size_t j = 0; j < tris.size(); j++) {
                write_vertex("a.im", tris[j].a);
                write_vertex("b.im", tris[j].b);
                write_vertex("c.im", tris[j].c);
                printf("bleu");
            }
        }

        // Replace fid by the first

        for (size_t j = 0; j < tris.size(); j++) {
            const auto &tri = tris[j];
            if (j == 0) {
                Fv[fid] = {tri.a, tri.b, tri.c};
            } else {
                Face *new_f = new Face;
                new_f->oids[0] = F[fid]->oids[0];
                new_f->oids[1] = F[fid]->oids[1];
                Fv.push_back({tri.a, tri.b, tri.c});
                F.push_back(new_f);
            }   
        }

        VNode_free_list(polygon);
        VNode_free_list(reflex);
        VNode_free_list(ears);
        VNode_free_list(convex);
    }

    //printf("Total faces: %zu\n", F.size());
}

void Dcel::update_hedge_faces(std::vector<Face *> &new_F)
{
    for (Hedge *h : H) h->left = NULL;

    for (Face *f : new_F) 
    {
        assert(f);

        Hedge *h = f->rep;
        assert(h);
        Cycle *c = h->cycle;
        assert(c->inout == Cycle::INNER);

        h->left = f;
        assert(f->rep->left == f);

        Hedge *w = h->next;
        while (w != h) 
        {
            assert(w->left == NULL);
            w->left = f;
            w = w->next;
        }
    }

    for (Face *f : new_F) {
        assert(f->rep->left == f);
    }

    for (Hedge *h : H) {
        Face *f = h->left;
        Hedge *w = h->next;
        while (w != h) {
            assert(w->left == f);
            w = w->next;
        }
    }
}

E_Int Dcel::check_faces(const std::vector<Hedge *> &H,
    const std::vector<Face *> &F)
{
    for (size_t i = 0; i < H.size(); i++) {
        Hedge *h = H[i];
        if (h->prev->left != h->left) { assert(0); return 1; }
        if (h->next->left != h->left) { assert(0); return 1; }
    }

    for (size_t i = 0; i < F.size(); i++) {
        Face *f = F[i];
        //Hedge *h = f->rep;
        //Cycle *c = h->cycle;
        //assert(c->inout == Cycle::DEGEN || c->inout == Cycle::INNER);
        if (f->rep->left != f) { assert(0); return 1; }
    }

    puts("CHECK: FACES OK.");

    return 0;
}



void Dcel::set_face_labels(std::vector<Face *> &new_F)
{
    // Label each face with the ids of the original faces containing it

    for (size_t i = 0; i < new_F.size(); i++) {
        Face *f = new_F[i];

        // Get the first RED and BLACK half-edges in the face cycle
        Hedge *h = f->rep;

        Cycle *c = h->cycle;
        assert(c->inout == Cycle::INNER);

        Hedge *R = NULL;
        Hedge *B = NULL;
        E_Int RB = 0;

        if (h->color == Dcel::RED) {
            R = h;
            B = get_hedge_of_color(f, Dcel::BLACK);
            if (B) RB = 1;
        } else if (h->color == Dcel::BLACK) {
            B = h;
            R = get_hedge_of_color(f, Dcel::RED);
            if (R) RB = 1;
        } else {
            assert(0);
        }

        if (RB) {
            // RED face should always exist
            assert(R->left);
            f->oids[Dcel::RED] = R->left->oids[Dcel::RED];
            assert(f->oids[Dcel::RED] != -1);

            // BLACK face might not exist
            f->oids[Dcel::BLACK] = (B->left) ? B->left->oids[Dcel::BLACK] : -1;
        } else {
            // Only single color possible is RED
            if (!R && B) {
                write_hedge("black.im", B);
                point_write("orig.im", B->orig->x, B->orig->y, B->orig->z);
            }
            assert(R && !B);
            //Hedge *REP = (R != NULL) ? R : B;
            //if (REP != R) {
            //    hedge_write("black", REP);
            //}
            f->oids[Dcel::RED] = R->left->oids[Dcel::RED];
        }
    }
}

Dcel::Hedge *Dcel::get_hedge_of_color(Face *f, int color)
{
    Hedge *h = f->rep;
    if (h->color == color) return h;
    Hedge *w = h->next;
    while (w != h) {
        if (w->color == color) return w;
        w = w->next;
    }
    return NULL;
}

std::vector<Dcel::Face *> Dcel::make_cycle_faces(const std::vector<Cycle *> &C)
{
    std::vector<Face *> new_F;

    for (Cycle *c : C) {
        //if (c->inout == Cycle::HOLE) continue;
        if (c->inout != Cycle::INNER) continue;
        
        // Create a face record
        Face *f = new Face;

        // Set its rep hedge to some edge of the cycle
        Hedge *h = c->rep;
        f->rep = h;

        new_F.push_back(f);
    }

    return new_F;
}

void Dcel::init_hedges_and_faces(const Smesh &Mf, int color)
{
    printf("Doing color %d\n", color);

    size_t current_nh = H.size();

    H.reserve(current_nh + Mf.ne*2);

    std::vector<std::vector<Hedge *>> list(Mf.np);

    // Create hedge records

    const auto &E = Mf.E;

    for (E_Int i = 0; i < Mf.ne; i++) {
        const auto &e = E[i];

        E_Int p = e.p;
        E_Int q = e.q;

        Hedge *h = NULL;
        Hedge *t = NULL;

        {
            Vertex tmp(Mf.X[p], Mf.Y[p], Mf.Z[p]);
            auto it = vertex_set.find(&tmp);
            assert(it != vertex_set.end());
            Vertex *P = *it;
            h = new Hedge(P, color);
            list[p].push_back(h);
        }

        {
            Vertex tmp(Mf.X[q], Mf.Y[q], Mf.Z[q]);
            auto it = vertex_set.find(&tmp);
            assert(it != vertex_set.end());
            Vertex *Q = *it;
            t = new Hedge(Q, color);
            list[q].push_back(t);
        }

        h->twin = t;
        t->twin = h;

        //h->oid = i;
        //t->oid = i;

        H.push_back(h);
        H.push_back(t);
    }
    
    // Pair-up hedges

    const auto &pnormals = Mf.pnormals;

    for (E_Int pid = 0; pid < Mf.np; pid++) {
        auto &hedges = list[pid];

        const E_Float *N = &pnormals[3*pid];
        sort_leaving_hedges(hedges, N);

        for (size_t i = 0; i < hedges.size(); i++) {
            Hedge *h = hedges[i];
            Hedge *w = hedges[(i+1)%hedges.size()];
            h->twin->next = w;
            w->prev = h->twin;
        }
    }

    // Create face records

    const auto &F2E = Mf.F2E;
    const auto &E2F = Mf.E2F;

    size_t current_nf = F.size();

    F.reserve(current_nf + Mf.nf);

    for (E_Int fid = 0; fid < Mf.nf; fid++) {
        const auto &pe = F2E[fid];
        E_Int first_edge = pe[0];
        E_Int where = 2*first_edge;
        Hedge *h = H[current_nh + where];
        Hedge *t = H[current_nh + where + 1];
        assert(h->twin == t);
        assert(t->twin == h);

        Face *f = new Face;
        f->oids[color] = fid;

        // Face must lie to the left of hedge
        assert(E2F[first_edge][0] == fid || E2F[first_edge][1] == fid);
        Hedge *REP = (E2F[first_edge][0] == fid) ? h : t;

        assert(REP->left == NULL);

        f->rep = REP;
        REP->left = f;
        Hedge *w = REP->next;
        while (w != REP) { w->left = f; w = w->next; }
        
        F.push_back(f);
    }
}

void Dcel::sort_leaving_hedges(std::vector<Hedge *> &leaving,
    const E_Float N[3]) const
{
    // Choose a vector that is not parallel to N

    E_Float ref_vec[3] = {0, N[2], -N[1]};
    E_Float NORM = K_MATH::norm(ref_vec, 3);
    if (Sign(NORM) == 0) {
        ref_vec[0] = -N[2];
        ref_vec[1] = 0;
        ref_vec[2] =  N[0];
        NORM = K_MATH::norm(ref_vec, 3);
        assert(Sign(NORM) != 0);
    }

    E_Float dp = K_MATH::dot(ref_vec, N, 3);
    for (E_Int i = 0; i < 3; i++) ref_vec[i] = ref_vec[i] - dp * N[i];
    NORM = K_MATH::norm(ref_vec, 3);
    for (E_Int i = 0; i < 3; i++) ref_vec[i] /= NORM;

    std::vector<E_Float> angles(leaving.size());

    for (size_t i = 0; i < leaving.size(); i++) {
        Hedge *h = leaving[i];
        Hedge *t = h->twin;

        Vertex *P = h->orig;
        Vertex *Q = t->orig;
        assert(P != Q);

        // Project the hedge onto the plane (pid, N)
        E_Float PQ[3] = {Q->x-P->x, Q->y-P->y, Q->z-P->z};
        E_Float dp = K_MATH::dot(PQ, N, 3);
        E_Float PQ_proj[3];
        for (E_Int j = 0; j < 3; j++) PQ_proj[j] = PQ[j] - dp * N[j];

        E_Float costheta = K_MATH::dot(ref_vec, PQ_proj, 3) / K_MATH::norm(PQ_proj, 3);
        costheta = std::min(costheta, 1.0);
        costheta = std::max(costheta, -1.0);
        E_Float angle = acos(costheta);
        
        // Determine the direction of the angle
        E_Float C[3];
        K_MATH::cross(ref_vec, PQ_proj, C);
        
        if (K_MATH::dot(N, C, 3) > 0) angle = 2*K_MATH::PI - angle;
        
        angles[i] = angle;
    }

    std::vector<E_Int> indices(leaving.size());
    for (size_t i = 0; i < leaving.size(); i++) indices[i] = i;

    std::sort(indices.begin(), indices.end(), [&](E_Int i, E_Int j)
    {
        if (angles[i] < angles[j]) return true;
        
        if (angles[i] > angles[j]) return false;
        
        Hedge *h = leaving[i];
        Hedge *w = leaving[j];

        if (h->color == w->color) {
            write_hedge("h.im", h);
            write_hedge("w.im", w);
        }

        assert(h->color != w->color);

        Vertex *O_h = h->orig;
        Vertex *T_h = h->twin->orig;

        Vertex *O_w = w->orig;
        Vertex *T_w = w->twin->orig;

        assert(O_w == O_h);
        assert(T_w == T_h);

        // If this is the origin: red comes first
        if (cmp_vtx(O_w, T_w) < 0) {
            if (h->color == RED) return true;
            return false;
        } else {
            if (h->color == RED) return false;
            return true;
        }
    });
    
    std::vector<Hedge *> tmp(leaving);
    for (size_t i = 0; i < leaving.size(); i++) leaving[i] = tmp[indices[i]];
}

E_Int Dcel::check_hedges(const std::vector<Hedge *> &H)
{
    for (size_t i = 0; i < H.size(); i++) {
        Hedge *h = H[i];
        if (h->prev->next != h) { assert(0); return 1; }
        if (h->next->prev != h) { assert(0); return 1; }
        if (h->twin->twin != h) { assert(0); return 1; }
        if (h->twin->next->orig != h->orig) { assert(0); return 1; }
        if (h->prev->twin->orig != h->orig) { assert(0); return 1; }
    }

    puts("CHECK: EDGES OK.");

    return 0;
}

void Dcel::make_cycles()
{
    C.clear();
    for (Hedge *h : H) h->cycle = NULL;

    for (size_t i = 0; i < H.size(); i++) {
        Hedge *h = H[i];

        if (h->cycle) continue;

        Cycle *c = new Cycle(h);
        C.push_back(c);

        h->cycle = c;

        Hedge *w = h->next;
        while (w != h) {
            w->cycle = c;
            w = w->next;
        }
    }
}

void Dcel::set_cycles_inout()
{
    inner = 0;
    degen = 0;
    hole = 0;

    for (Cycle *cycle : C) {

        // Hole cycle is a cycle where all hedges have a null left face
        bool is_hole = true;

        Hedge *h = cycle->rep;
        if (!h->left) {
            Hedge *w = h->next;
            while (w != h) {
                if (w->left) {
                    is_hole = false;
                    break;
                }
                w = w->next;
            }
        } else {
            is_hole = false;
        }

        if (is_hole) {
            cycle->inout = Cycle::HOLE;
            hole++;
            continue;
        }

        // Get the leftmost vertex in the cycle
        Vertex *b = h->orig;

        Hedge *e2 = h; // Half-edge starting at leftmost vertex

        Hedge *w = h->next;
        while (w != h) {
            Vertex *p = w->orig;
            E_Int cmp = cmp_vtx(p, b);
            if (cmp < 0) {
                b = p;
                e2 = w;
            }

            w = w->next;
        }

        Hedge *e1 = e2->prev;

        assert(e2->orig == b);
        assert(e1->twin->orig == b);

        Vertex *a = e1->orig;
        Vertex *c = e2->twin->orig;

        // Vectors ab and bc
        E_Float ab[3] = {b->x-a->x, b->y-a->y, b->z-a->z};
        E_Float bc[3] = {c->x-b->x, c->y-b->y, c->z-b->z};

        E_Float N[3];
        K_MATH::cross(ab, bc, N);

        if (Sign(K_MATH::norm(N, 3) == 0)) {
            cycle->inout = Cycle::DEGEN;
            degen++;
        } else {
            cycle->inout = Cycle::INNER;
            inner++;
        }
    }

    printf("Inner cycles: " SF_D_ "\n", inner);
    printf("Degen cycles: " SF_D_ "\n", degen);
    printf("Hole cycles: " SF_D_ "\n", hole);
}

Dcel::~Dcel()
{
    for (Vertex *v : V) delete v;
    for (Hedge *h : H) delete h;
    for (Face *f : F) delete f;
    for (Cycle *c : C) delete c;
}

Dcel::Dcel(const Smesh &Mf, int color)
{
    const auto &X = Mf.X;
    const auto &Y = Mf.Y;
    const auto &Z = Mf.Z;

    V.reserve(Mf.np);
    for (E_Int pid = 0; pid < Mf.np; pid++) {
        Vertex *v = new Vertex(X[pid], Y[pid], Z[pid]);
        v->oids[color] = pid;
        vertex_set.insert(v);
        V.push_back(v);
    }

    init_hedges_and_faces(Mf, color);

    if (check_hedges(H) != 0) {
        fprintf(stderr, "Dcel: Inconsistent half-edge records!\n");
        abort();
    }

    if (check_faces(H, F) != 0) {
        fprintf(stderr, "Dcel: Inconsistent face records!\n");
        abort();
    }

    make_cycles();
    set_cycles_inout();
}


