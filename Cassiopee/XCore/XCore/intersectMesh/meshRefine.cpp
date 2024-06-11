#include <stack>

#include "mesh.h"
#include "../common/common.h"

int meshes_mutual_refinement(Mesh &M, Mesh &S)
{
    size_t refM, refS;
    int iter = 0;

    do {
        iter++;
        /*
        if (iter > 1) {
            RAISE("Cannot refine surface meshes more than once, for now.\n");
            return 1;
        }
        */
        refM = M.refine(M.patch, S.patch, S);
        refS = S.refine(S.patch, M.patch, M);
    } while (refM > 0 || refS > 0);

    return 0;
}

size_t Mesh::refine(std::set<E_Int> &mpatch, std::set<E_Int> &spatch,
    Mesh &S)
{
    // Isolate spatch points
    std::set<E_Int> spoints;
    for (E_Int sface : spatch) {
        assert(S.face_is_active(sface));
        const auto &spn = S.F[sface];
        for (E_Int sp : spn)
            spoints.insert(sp);
    }

    // Locate spatch points within mpatch
    std::map<E_Int, std::vector<pointFace>> spoints_to_mfaces;

    for (E_Int spt : spoints) {
        auto pf = locate(S.X[spt], S.Y[spt], mpatch);
        spoints_to_mfaces[spt] = pf;
    }

    // Mface to spoints located within it
    std::map<E_Int, std::vector<E_Int>> mfpoints;

    for (const auto &spt_to_mfs : spoints_to_mfaces) {
        E_Int spt = spt_to_mfs.first;
        const auto &mptfaces = spt_to_mfs.second;

        for (const auto &mptface : mptfaces) {
            E_Int mface = mptface.F;
            mfpoints[mface].push_back(spt);
        }
    }

    // Keep the mfaces which contain 3 points or more
    std::map<E_Int, std::vector<E_Int>> filtered_mfaces_map;

    for (const auto &mfpts : mfpoints) {
        if (mfpts.second.size() >= 3) {
            filtered_mfaces_map.insert({mfpts.first, mfpts.second});
        }
    }

    // Make point-to-face connectivity for S
    S.make_point_faces();

    // Sensor
    std::map<E_Int, std::vector<E_Int>> ref_mfaces_to_sfaces;

    for (const auto &mfdata : filtered_mfaces_map) {
        E_Int mface = mfdata.first;
        const auto &spoints = mfdata.second;

        // For every sface, how many points are contained within mface?
        std::map<E_Int, size_t> scontained;

        for (E_Int sp : spoints) {
            // Sfaces sharing sp
            const auto &sfaces = S.P2F[sp];
 
            for (E_Int sface : sfaces) {
                assert(S.face_is_active(sface));

                // Sface should belong to spatch
                if (spatch.find(sface) == spatch.end()) continue;

                // mface contains another point of sface
                scontained[sface]++;
            }
        }

        for (const auto &sfdat : scontained) {
            E_Int sface = sfdat.first;

            // The number of points of sface contained within mface...
            size_t npts = sfdat.second;

            // ... should be at most equal to the stride of the sface!
            assert(npts <= S.F[sface].size());

            // Discard sface if it is not completely enclosed by mface
            if (npts < S.F[sface].size()) continue;

            // Discard sface if it is a duplicate of mface
            if (faces_are_dups(mface, sface, S)) continue;

            // Add sface to the set of faces enclosed by mface
            ref_mfaces_to_sfaces[mface].push_back(sface);
        }
    }

    printf("Faces to refine: %lu\n", ref_mfaces_to_sfaces.size());

    E_Int iter = 0;
    size_t ret = 0;

    while (!ref_mfaces_to_sfaces.empty()) {
        iter++;

        auto ref_data = smooth_ref_data(ref_mfaces_to_sfaces);

        auto ref_faces = prepare_for_refinement(ref_data);

        ret += ref_faces.size();

        refine_faces(ref_faces);

        std::map<E_Int, std::vector<E_Int>> new_sensor;

        for (const auto &fdata : ref_mfaces_to_sfaces) {
            E_Int parent = fdata.first;
            const auto &sfaces = fdata.second;

            const auto &children = fchildren[parent];

            for (E_Int child : children) {
                for (E_Int sface : sfaces) {
                    if (face_contains_sface(child, sface, S)) {
                        assert(!faces_are_dups(child, sface, S));
                        new_sensor[child].push_back(sface);
                    }
                }
            }

            // Update mpatch
            patch.erase(parent);
            for (E_Int child : children) patch.insert(child);
        }

        ref_mfaces_to_sfaces = new_sensor;
    }

    return ret;
}


std::vector<E_Int> Mesh::smooth_ref_data(
    const std::map<E_Int, std::vector<E_Int>> &sensor)
{
    // TODO(Imad): shouldn't the size be factive.size() instead of nf?
    std::vector<E_Int> ref_data(nf, 0);
    std::stack<E_Int> stk;

    for (const auto &fdata : sensor) {
        assert(face_is_active(fdata.first));
        ref_data[fdata.first] = 1; // Warning: size must be nf for this to be ok
        stk.push(fdata.first);
    }

    /*
    while (!stk.empty()) {
        E_Int face = stk.top();
        stk.pop();

        E_Int face_incr = ref_data[face] + flevel[face];

        auto neis = get_active_neighbours(face);

        for (auto nei : neis) {
            assert(nei) != -1;

            E_Int nei_incr = ref_data[nei] + flevel[nei];

            E_Int diff = abs(face_incr - nei_incr);

            if (diff <= 1) continue;

            E_Int fid = face_incr > nei_incr ? nei : face;

            ref_data[fid] += 1;

            stk.push(fid);
        }
    }
    */

    return ref_data;
}

std::vector<E_Int> Mesh::prepare_for_refinement(const std::vector<E_Int> &ref_data)
{
    std::vector<E_Int> ref_faces;
    
    for (size_t i = 0; i < ref_data.size(); i++) {
        if (ref_data[i] > 0) ref_faces.push_back(i);
    }

    // Refine the lower-level faces first
    std::sort(ref_faces.begin(), ref_faces.end(),
        [&] (E_Int i, E_Int j) { return flevel[i] < flevel[j]; });
    
    // Resize data structures
    resize_point_data(ref_faces.size());
    resize_face_data(ref_faces.size());

    return ref_faces;
}

/*
std::vector<E_Int> Mesh::get_active_neighbours(E_Int face)
{
    return std::vector<E_Int>();
}
*/

void Mesh::refine_faces(const std::vector<E_Int> &ref_faces)
{
    for (E_Int ref_face : ref_faces) {
        if (face_is_tri(ref_face)) refine_tri(ref_face);
        else refine_quad(ref_face);
    }
}

void Mesh::resize_point_data(size_t nref_faces)
{
    size_t nnew_points = np + nref_faces * 5;
    X.resize(nnew_points);
    Y.resize(nnew_points);
    Z.resize(nnew_points);
}

void Mesh::resize_face_data(size_t nref_faces)
{
    size_t nnew_faces = nf + nref_faces * 4;
    F.resize(nnew_faces);
    flevel.resize(nnew_faces, -1);
}

void Mesh::refine_quad(E_Int quad)
{
    // Refine the edges
    const auto &pn = F[quad];
    E_Int ec[4];

    for (size_t i = 0; i < pn.size(); i++) {
        E_Int p = pn[i];
        E_Int q = pn[(i+1)%pn.size()];
        UEdge edge(p, q);

        auto it = ecenter.find(edge);

        if (it == ecenter.end()) {
            X[np] = 0.5 * (X[p] + X[q]);
            Y[np] = 0.5 * (Y[p] + Y[q]);
            Z[np] = 0.5 * (Z[p] + Z[q]);
            ecenter[edge] = np;
            ec[i] = np;
            np++;
        } else {
            ec[i] = it->second;
        }
    }

    // Face centroid

    X[np] = Y[np] = Z[np] = 0.0;
    for (E_Int i = 0; i < 4; i++) {
        E_Int p = pn[i];
        X[np] += X[p];
        Y[np] += Y[p];
        Z[np] += Z[p];
    }
    X[np] *= 0.25;
    Y[np] *= 0.25;
    Z[np] *= 0.25;

    // New face points
    E_Int nf0 = nf, nf1 = nf+1, nf2 = nf+2, nf3 = nf+3;

    F[nf0] = { pn[0], ec[0], np,    ec[3] };
    F[nf1] = { ec[0], pn[1], ec[1], np    };
    F[nf2] = { np   , ec[1], pn[2], ec[2] };
    F[nf3] = { ec[3], np   , ec[2], pn[3] };

    // Disable quad and enable its children
    factive.erase(quad);
    factive.insert(nf0);
    factive.insert(nf1);
    factive.insert(nf2);
    factive.insert(nf3);

    // Set quad children pointers
    fchildren[quad] = { nf0, nf1, nf2, nf3 };
    flevel[nf0] = flevel[nf1] = flevel[nf2] = flevel[nf3] = flevel[quad] + 1;

    np += 1;
    nf += 4;
}

void Mesh::refine_tri(E_Int tri)
{
    printf(SF_D_ "\n", tri);
    assert(0);
}