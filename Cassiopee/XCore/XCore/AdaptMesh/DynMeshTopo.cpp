#include "DynMesh.h"

#define OUT 0
#define IN 1

#define INTERNAL 0
#define EXTERNAL 1

#define DSMALL 1e-14

E_Int DynMesh::orient_skin(E_Int normal_direction)
{
    // flag external cells and faces
    std::vector<E_Int> fflags, efaces;
    flag_and_get_external_faces(fflags, efaces);

    // extract external faces connectivity
    std::vector<E_Int> fadj;
    std::vector<E_Int> xadj(1, 0);
    for (E_Int i = 0; i < nf; i++) {
        if (fflags[i] == EXTERNAL) {
            const auto &pn = F[i];
            xadj.push_back(pn.size());
            for (size_t j = 0; j < pn.size(); j++) fadj.push_back(pn[j]);
        }
    }

    E_Int nefaces = (E_Int)efaces.size();

    for (E_Int i = 0; i < nefaces; i++) xadj[i+1] += xadj[i];

    // build skin neighbourhood
    std::vector<E_Int> fneighbours;
    K_CONNECT::build_face_neighbourhood(fadj, xadj, fneighbours);

    // color the faces by connex part
    std::vector<E_Int> colors(xadj.size()-1);
    E_Int nconnex = K_CONNECT::colorConnexParts(&fneighbours[0], &xadj[0],
        nefaces, &colors[0]);

    //printf("orient_boundary(): connex parts: " SF_D_ "\n", nconnex);

    assert(efaces.size() == xadj.size()-1);
    std::vector<E_Int> forient(nefaces, 0);
    std::vector<E_Int> cflags;
    E_Int ret = 0;

    if (nconnex > 1) {
        // extract nconnex nface-ngon for separate orientation
        for (E_Int color = 0; color < nconnex; color++) {
            std::vector<bool> keep_pgs(nf, false);
            for (E_Int i = 0; i < nefaces; i++) {
                keep_pgs[efaces[i]] = (colors[i] == color);
            }
            
            // extract nface corresponding to kept faces
            std::vector<E_Int> NFACE, cxadj(1, 0), cells;
            extract_nface_of_kept_pgs(keep_pgs, NFACE, cxadj, cells);

            std::vector<E_Int> cflags;
            flag_marked_external_cells(cells, fflags, cflags);

            ret |= orient_boundary((E_Int)cells.size(),
                &fadj[0], &xadj[0], nefaces, &fneighbours[0], &efaces[0],
                forient, cflags, fflags, &cells[0], normal_direction);
        }
    } else {
        std::vector<E_Int> cflags;
        flag_all_external_cells(fflags, cflags);
        ret = orient_boundary(nc, &fadj[0], &xadj[0], nefaces, &fneighbours[0],
            &efaces[0], forient, cflags, fflags, NULL, normal_direction);
    }

    // Apply orientation
    E_Int nrev = 0;
    for (E_Int i = 0; i < nefaces; i++) {
        if (forient[i] == -1) {
            E_Int face = efaces[i]; // 0-based
            auto &pn = F[face];
            std::reverse(pn.begin(), pn.end());
            nrev++;
        }
    }

    //printf("orient_boundary(): reversed " SF_D_ " faces\n", nrev);

    return ret;
}

void DynMesh::flag_and_get_external_faces(std::vector<E_Int> &fflags,
    std::vector<E_Int> &efaces)
{
    std::vector<E_Int> face_count(nf, 0);

    // Loop through the elements and increment face_count
    for (E_Int i = 0; i < nc; i++) {
        const auto &pf = C[i];
        for (E_Int face : pf) face_count[face]++;
    }

    // External faces are those with a count equal to 1
    fflags.resize(nf);

    for (E_Int i = 0; i < nf; i++) {
        if (face_count[i] == 1) {
            fflags[i] = EXTERNAL;
            efaces.push_back(i);
        } else {
            fflags[i] = INTERNAL;
        }
    }
}

void DynMesh::extract_nface_of_kept_pgs(const std::vector<bool> &kept_pgs,
    std::vector<E_Int> &NFACE, std::vector<E_Int> &xadj,
    std::vector<E_Int> &cells)
{
    NFACE.clear();
    xadj.resize(1, 0);
    cells.clear();

    for (E_Int i = 0; i < nc; i++) {
        const auto &pf = C[i];
        bool keep = false;
        for (size_t j = 0; j < pf.size() && !keep; j++) keep = kept_pgs[pf[j]];
        if (keep) {
            cells.push_back(i);
            xadj.push_back(pf.size());
            for (size_t j = 0; j < pf.size(); j++) NFACE.push_back(pf[j]);
        }
    }

    for (size_t i = 0; i < xadj.size(); i++) xadj[i+1] += xadj[i];
}

void DynMesh::flag_marked_external_cells(const std::vector<E_Int> &cells,
    const std::vector<E_Int> &fflags, std::vector<E_Int> &cflags)
{
    // External cells are those with at least one external face
    cflags.resize(cells.size(), INTERNAL);

    for (size_t i = 0; i < cells.size(); i++) {
        E_Int cell = cells[i];
        const auto &pf = C[cell];
        for (size_t j = 0; j < pf.size(); j++) {
            E_Int face = pf[j];
            if (fflags[face] == EXTERNAL) {
                cflags[i] = EXTERNAL;
                break;
            }
        }
    }
}

E_Int DynMesh::orient_boundary(E_Int ncells, E_Int *efadj, E_Int *efxadj,
    E_Int nefaces, E_Int *fneis, E_Int *efaces, std::vector<E_Int> &forient,
    const std::vector<E_Int> &cflags, const std::vector<E_Int> &fflags,
    E_Int *cells, E_Int normal_direction)
{
    // Look for a cell whose volume is unambiguously computed
    E_Float cvol = 0.0;
    E_Int seed = -1;
    E_Int refPG = -1;
    E_Int refIdx = -1;

    while (++seed < ncells) {
        if (cflags[seed] != EXTERNAL) continue;

        E_Int cid = (cells != NULL) ? cells[seed] : seed;

        const auto &pf = C[cid];
        refPG = -1;
        E_Int local_idx = -1;
        for (size_t j = 0; j < pf.size(); j++) {
            E_Int face = pf[j];
            if (fflags[face] == EXTERNAL) {
                refPG = face;
                local_idx = j;
                break;
            }
        }

        if (refPG == -1) {
            fprintf(stderr, "orient_boundary(): couldn't find an external face "
                            "within external cell " SF_D_ "\n", cid);
            return 1;
        }

        // Look for index of refPG in efaces (0-based)
        refIdx = -1;
        for (E_Int i = 0; i < nefaces; i++) {
            if (efaces[i] == refPG) {
                refIdx = i;
                break;
            }
        }

        if (refIdx == -1) {
            fprintf(stderr, "orient_boundary(): couldn't find reference face "
                            SF_D_ " in external faces list\n", refPG);
            return 1;
        }

        // Set orientation of refPG to +1.
        // Reorient seed's faces based on orientation of refPG.
        // Compute cvol, the volume of seed.
        // If cvol > 0, orientation of all faces including refPG, is outwards
        // Otherwise, set orientation of refPG to -1.

        compute_cell_volume(cid, cvol, local_idx);

        if (fabs(cvol) < DSMALL) continue;

        // set reference orientation of refPG and exit
        forient[refIdx] = (cvol > 0.0) ? 1 : -1;

        if (normal_direction == IN) forient[refIdx] = -forient[refIdx];

        break;
    }

    if (seed >= ncells) {
        fprintf(stderr, "orient_boundary_ngon(): couldn't find reference "
                        "polyhedron\n");
        assert(0);
        return 1;
    }

    // propagate
    K_CONNECT::reversi_connex(efadj, efxadj, nefaces, fneis, refIdx, forient);

    return 0;
}

void DynMesh::compute_cell_volume(E_Int cell, E_Float &vol, E_Int refIdx)
{
    // Orient the faces coherently
    std::vector<E_Int> NGON;
    std::vector<E_Int> INDPG(1, 0);
    const auto &pf = C[cell];
    E_Int stride = (E_Int)pf.size();

    for (E_Int i = 0; i < stride; i++) {
        E_Int face = pf[i];
        const auto &pn = F[face];
        E_Int np = pn.size();
        INDPG.push_back(np);
        for (E_Int j = 0; j < np; j++) NGON.push_back(pn[j]);
    }

    for (E_Int i = 0; i < stride; i++) INDPG[i+1] += INDPG[i];

    // Fix orientation of first face
    std::vector<E_Int> orient(stride);
    orient[refIdx] = 1;
    std::vector<E_Int> neis(NGON.size());
    K_CONNECT::build_face_neighbourhood(NGON, INDPG, neis);
    K_CONNECT::reversi_connex(&NGON[0], &INDPG[0], stride, &neis[0], refIdx,
        orient);

    // Apply orientation in local NGON
    for (E_Int i = 0; i < stride; i++) {
        if (orient[i] == -1) {
            E_Int start = INDPG[i];
            E_Int np = INDPG[i+1] - start;
            E_Int *pn = &NGON[start];
            std::reverse(pn+1, pn+np);
        }
    }

    // Compute faces area and center
    std::vector<E_Float> faceAreas(3*stride, 0.0);
    std::vector<E_Float> faceCenters(3*stride, 0.0);

    for (E_Int i = 0; i < stride; i++) {
        E_Int face = pf[i];
        E_Float *fa = &faceAreas[3*i];
        E_Float *fc = &faceCenters[3*i];
        E_Int np = INDPG[i+1]-INDPG[i];
        E_Int *pn = &NGON[INDPG[i]];
        for (E_Int j = 0; j < np; j++) pn[j] += 1;
        K_METRIC::compute_face_center_and_area(face, np, pn, X.data(), Y.data(),
            Z.data(), fc, fa);
    }
  
    // Estimate cell centroid as average of face centers
    E_Float cc[3] = {0,0,0};
    for (E_Int i = 0; i < stride; i++) {
        E_Float *fc = &faceCenters[3*i];
        for (E_Int j = 0; j < 3; j++) cc[j] += fc[j];
    }
    for (E_Int i = 0; i < 3; i++) cc[i] /= stride;

    // Compute cell volume
    vol = 0.0;

    for (E_Int i = 0; i < stride; i++) {
        E_Float *fa = &faceAreas[3*i];
        E_Float *fc = &faceCenters[3*i];

        // Compute 3*face-pyramid volume contribution
        E_Float d[3] = {fc[0]-cc[0], fc[1]-cc[1], fc[2]-cc[2]};
        //E_Float pyr3Vol = K_MATH::dot(fa, fc, 3);
        E_Float pyr3Vol = K_MATH::dot(fa, d, 3);

        vol += pyr3Vol;
    }

    vol *= K_MATH::ONE_THIRD;
}

void DynMesh::flag_all_external_cells(const std::vector<E_Int> &fflags,
    std::vector<E_Int> &cflags)
{
    // External cells are those with at least one external face
    cflags.resize(nc, INTERNAL);
    for (E_Int i = 0; i < nc; i++) {
        const auto &pf = C[i];
        E_Int stride = (E_Int)pf.size();
        for (E_Int j = 0; j < stride; j++) {
            E_Int face = pf[j];
            if (fflags[face] == EXTERNAL) {
                cflags[i] = EXTERNAL;
                break;
            }
        }
    }
}
