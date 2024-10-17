#include "smesh.h"
#include "primitives.h"

void Smesh::get_unit_projected_direction(E_Int fid, const E_Float D[3],
    E_Float proj[3]) const
{
    assert(fid >= 0);
    assert(fid < nf);

    // Unit normal
    const E_Float *fN = &fnormals[3*fid];

    E_Float dp = K_MATH::dot(D, fN, 3); 

    proj[0] = D[0] - dp * fN[0];
    proj[1] = D[1] - dp * fN[1];
    proj[2] = D[2] - dp * fN[2];
    E_Float NORM = K_MATH::norm(proj, 3);
    proj[0] /= NORM, proj[1] /= NORM, proj[2] /= NORM;
}

E_Int Smesh::deduce_face(const std::vector<E_Int> &pf,
    E_Float ox, E_Float oy, E_Float oz, E_Float D[3], 
    E_Int last_vertex, E_Int last_edge) const
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
        get_unit_projected_direction(fid, D, proj);

        const auto &pn = Fc[fid];
        const auto &pe = F2E[fid];
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
                X[p], Y[p], Z[p],
                X[q], Y[q], Z[q],
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