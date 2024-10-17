#include "dcel.h"
#include "smesh.h"
#include "primitives.h"
#include "io.h"





E_Int Dcel::trace_hedge_2(Hedge *sh, const Smesh &M, const Smesh &S, E_Int hid)
{
    assert(0);
    /*
    Vertex *O = sh->orig;
    Vertex *T = sh->twin->orig;
    E_Float D[3] = {T->x-O->x, T->y-O->y, T->z-O->z};
    E_Float NORM = K_MATH::norm(D, 3);
    D[0] /= NORM, D[1] /= NORM, D[2] /= NORM;

    // The origin could be inside a face, on an edge or on a vertex
    // If it is inside a face, no problem
    // If it is on an edge/vertex, which face is intersected first by sh?

    std::vector<E_Int> orig_faces, tail_faces;
    E_Int last_vertex = -1, last_edge = -1;
    E_Int dummy;
    M.get_shared_faces(O->loc, orig_faces, last_vertex, last_edge);
    M.get_shared_faces(T->loc, tail_faces, dummy, dummy);

    // If O is inside fid, we could skip this check
    // We keep it for consistency
    E_Int starting_face = M.deduce_face(
        orig_faces, O->x, O->y, O->z, D,
        last_vertex, last_edge
    );
    assert(starting_face != -1);

    bool found_tail = false;
    E_Int current_fid = starting_face;
    E_Float current_pos[3] = {O->x, O->y, O->z};

    E_Int walk = 0;
    E_Int max_walks = 20;

    //if (hid == hid) hedge_write("sh", sh);

    while (!found_tail && walk <= max_walks) {

        //if (hid == hid) face_write("current_face", F[current_fid]);

        // We are on current_face

        // If the current_face shares the tail, stop
        for (auto fid : tail_faces) {
            if (fid == current_fid) {
                found_tail = true;
                break;
            }
        }

        if (found_tail) break;

        // The tail is not within current_face, nor is it on one of its
        // edges, nor is it one of its vertices
        // So we must:
        // 1 - project the traced hedge direction on the current face
        // 2 - shoot a ray {current position, projected dir} and intersect xedge
        // 3 - if intersection is on an xedge endpoint, deduce the next face
        //     else, travel to the neighbour face of current face wrt to xedge

        // Project D onto current face
        E_Float proj[3];
        M.get_unit_projected_direction(current_fid, D, proj);

        E_Int next_fid = -1;
        E_Float next_pos[3] = {EFLOATMAX, EFLOATMAX, EFLOATMAX};

        // Ray-edge intersection
        const auto &pn = M.F[current_fid];
        const auto &pe = M.F2E[current_fid];
        bool hit = false;
        for (size_t i = 0; i < pn.size(); i++) {
            E_Int p = pn[i];
            E_Int q = pn[(i+1)%pn.size()];
            E_Int e = pe[i];
            if (p == last_vertex || q == last_vertex || e == last_edge)
                continue;
            
            E_Float t, s;
            hit = ray_edge_intersect(
                current_pos[0], current_pos[1], current_pos[2],
                proj[0], proj[1], proj[2],
                M.X[p], M.Y[p], M.Z[p],
                M.X[q], M.Y[q], M.Z[q],
                t, s
            );

            if (hit) {

                //if (hid == hid) point_write("hit", current_pos[0] + t * proj[0], current_pos[1] + t * proj[1], current_pos[2] + t * proj[2]);

                // Intersection within the edge
                if (s > TOL && s < 1 - TOL) {

                    // Simply move to the neighbour face
                    const auto &pe = M.F2E[current_fid];
                    E_Int eid = pe[i];
                    last_edge = eid;
                    last_vertex = -1;
                    const auto &e = M.E[eid];
                    if (e.p == p) assert(e.q == q);
                    else {
                        assert(e.p == q);
                        assert(e.q == p);
                    }
                    assert(M.E2F[eid][0] == current_fid ||
                           M.E2F[eid][1] == current_fid);
                    if (M.E2F[eid][0] == current_fid)
                        next_fid = M.E2F[eid][1];
                    else
                        next_fid = M.E2F[eid][0];
                    // We should be within M still
                    assert(next_fid != -1);

                    // Compute the next_pos
                    next_pos[0] = current_pos[0] + t * proj[0];
                    next_pos[1] = current_pos[1] + t * proj[1];
                    next_pos[2] = current_pos[2] + t * proj[2];

                    // Register intersection
                    Event *xit = Q.lookup(next_pos[0], next_pos[1], next_pos[2]);
                    Vertex *x = NULL;
                    if (xit == NULL) {
                        xit = Q.insert(next_pos[0], next_pos[1], next_pos[2]);
                        x = xit->key;
                        x->id = V.size();
                        V.push_back(x);
                        x->meid = last_edge;
                    } else {
                        x = xit->key;
                        dup_x++;
                    }

                    hedge_intersections[sh].push_back(x);
                    // TODO(Imad): H[eid] or its twin?
                    Hedge *mh = H[2*eid];
                    if (cmp_vtx(mh->orig, mh->twin->orig) > 0)
                        mh = mh->twin;
                    hedge_intersections[mh].push_back(x);
                }

                // Intersection on an endpoint
                else {
                    bool hit_p = (s <= TOL);
                    bool hit_q = (s >= 1 - TOL);


                    last_edge = -1;

                    if (hit_p) {
                        assert(hit_q == false);

                        last_vertex = p;
 
                    } else if (hit_q) {

                        last_vertex = q;
                    }

                    next_pos[0] = M.X[last_vertex];
                    next_pos[1] = M.Y[last_vertex];
                    next_pos[2] = M.Z[last_vertex];
                    
                    const auto &pf = M.P2F[last_vertex];
                    
                    next_fid = M.deduce_face(pf,
                        next_pos[0], next_pos[1], next_pos[2], D,
                        last_vertex, last_edge
                    );
                    assert(next_fid != -1);

                    // Register intersection

                    Event *xit = Q.lookup(M.X[last_vertex],
                        M.Y[last_vertex], M.Z[last_vertex]);
                    assert(xit != NULL);
                    Vertex *x = xit->key;

                    // Edges from S do not cross
                    assert(x->oid[1] == -1);

                    assert(vertices_crossed.find(x) == vertices_crossed.end());
                    vertices_crossed.insert(x);

                    hedge_intersections[sh].push_back(x);

                }

                break;
            }
        }

        if (!hit) {
            hedge_write("sh", sh);
            face_write("current_face", F[current_fid]);
            assert(0);
        }

        assert(next_fid != current_fid);
        current_fid = next_fid;
        current_pos[0] = next_pos[0];
        current_pos[1] = next_pos[1];
        current_pos[2] = next_pos[2];
        walk++;
    }

    if (walk > max_walks) {
        fprintf(stderr, "Warning : Could not reach the tail of edge %d after %d max walks!", hid, max_walks);
        assert(0);
        return 1;
    }

    assert(found_tail);
    */

    return 0;
}
