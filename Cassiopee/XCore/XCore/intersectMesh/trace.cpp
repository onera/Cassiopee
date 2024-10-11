#include "dcel.h"
#include "smesh.h"
#include "primitives.h"
#include "hedge.h"
#include "io.h"
#include "event.h"

static
bool ray_edge_intersect(E_Float ox, E_Float oy, E_Float oz,
    E_Float dx, E_Float dy, E_Float dz,
    E_Float px, E_Float py, E_Float pz,
    E_Float qx, E_Float qy, E_Float qz,
    E_Float &t, E_Float &u)
{
    E_Float v[3]= {px-ox, py-oy, pz-oz};

    E_Float dl[3] = {qx-px, qy-py, qz-pz};

    E_Float dr[3] = {dx, dy, dz};

    /*
    E_Float w[3];
    K_MATH::cross(v, dl, w);
    E_Float det = K_MATH::dot(w, dr, 3);

    // ray and edge must be coplanar
    if (Sign(det) != 0) return false;
    */

    // ray and edge must not be parallel
    E_Float n[3];
    K_MATH::cross(dr, dl, n);
    E_Float denom = K_MATH::dot(n, n, 3);
    if (Sign(denom) == 0) return false;

    E_Float tmp[3];
    K_MATH::cross(v, dl, tmp);

    t = K_MATH::dot(tmp, n, 3) / denom;

    if (t < TOL) return false;

    K_MATH::cross(v, dr, tmp);

    u = K_MATH::dot(tmp, n, 3) / denom;

    if (u < -TOL || u > 1 + TOL) return false;

    return true;
}

static
void get_unit_projected_direction(E_Int fid, const Smesh &M, const E_Float D[3],
    E_Float proj[3])
{
    assert(fid >= 0);
    assert(fid < M.nf);

    // Unit normal
    const E_Float *fN = &M.fnormals[3*fid];

    E_Float dp = K_MATH::dot(D, fN, 3); 

    proj[0] = D[0] - dp * fN[0];
    proj[1] = D[1] - dp * fN[1];
    proj[2] = D[2] - dp * fN[2];
    E_Float NORM = K_MATH::norm(proj, 3);
    proj[0] /= NORM, proj[1] /= NORM, proj[2] /= NORM;
}

static
E_Int deduce_face(const std::vector<E_Int> &pf, const Smesh &M,
    E_Float ox, E_Float oy, E_Float oz, E_Float D[3], 
    E_Int last_vertex, E_Int last_edge)
{
    // Intersect the projection of D with all the faces in pf
    // At least one intersection must exist
    // Return the face with the earliest intersection

    // For debugging
    E_Int faces_hit = 0;

    E_Float t_min = EFLOATMAX;
    E_Int ret_face = -1;

    for (auto fid : pf) {

        // Compute the unit projection of D on this face

        E_Float proj[3];
        get_unit_projected_direction(fid, M, D, proj);

        const auto &pn = M.F[fid];
        const auto &pe = M.F2E[fid];
        assert(pn.size() == pe.size());

        for (size_t i = 0; i < pn.size(); i++) {

            E_Int p = pn[i];
            E_Int q = pn[(i+1)%pn.size()];
            E_Int e = pe[i];

            if (p == last_vertex || q == last_vertex || e == last_edge)
                continue;

            E_Float t, s;

            bool hit = ray_edge_intersect(ox, oy, oz,
                proj[0], proj[1], proj[2],
                M.X[p], M.Y[p], M.Z[p],
                M.X[q], M.Y[q], M.Z[q],
                t, s
            );

            if (hit) {
                faces_hit += 1;

                if (t < t_min) {
                    t_min = t;
                    ret_face = fid;
                }

                // Hit an edge of the face, stop
                break;
            }
        }
    }

    // We must have hit a face
    assert(faces_hit > 0);

    return ret_face;
}

static
void get_shared_faces(const Vertex *v, const Smesh &M, std::vector<E_Int> &ret,
    E_Int &pid, E_Int &eid)
{
    ret.clear();

    const auto &loc = v->loc;

    E_Int fid = loc.fid;
    assert(fid != -1);

    if (loc.e_idx != -1) {
        assert(loc.v_idx == -1);
        const auto &pe = M.F2E[fid];
        eid = pe[loc.e_idx];
        const auto &pf = M.E2F[eid];
        assert(pf[0] == fid || pf[1] == fid);
        ret.push_back(pf[0]);
        // O could be on a boundary edge
        if (pf[1] != -1) ret.push_back(pf[1]);
    }
    else if (loc.v_idx != -1) {
        assert(loc.e_idx == -1);
        const auto &pn = M.F[fid];
        pid = pn[loc.v_idx];
        const auto &pf = M.P2F[pid];
        // For consistency
        bool found_fid = false;
        for (auto face : pf) {
            ret.push_back(face);
            if (face == fid) {
                found_fid = true;
            }
        }
        assert(found_fid == true);
    }
    else {
        ret.push_back(fid);
    }
}

E_Int Dcel::trace_hedge_2(Hedge *sh, const Smesh &M, const Smesh &S, E_Int hid)
{
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
    get_shared_faces(O, M, orig_faces, last_vertex, last_edge);
    get_shared_faces(T, M, tail_faces, dummy, dummy);

    // If O is inside fid, we could skip this check
    // We keep it for consistency
    E_Int starting_face = deduce_face(
        orig_faces, M, O->x, O->y, O->z, D,
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
        get_unit_projected_direction(current_fid, M, D, proj);

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
                    
                    next_fid = deduce_face(pf, M,
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

    return 0;
}
