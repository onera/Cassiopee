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
#include "mesh.h"

#define OUT 0
#define IN 1

#define INTERNAL 0
#define EXTERNAL 1

#define DSMALL 1e-14

Int IMesh::orient_skin(Int normal_direction)
{
    // flag external cells and faces
    std::vector<Int> fflags, efaces;
    flag_and_get_external_faces(fflags, efaces);

    // extract external faces connectivity
    std::vector<Int> fadj;
    std::vector<Int> xadj(1, 0);
    for (Int i = 0; i < nf; i++) {
        if (fflags[i] == EXTERNAL) {
            const auto &pn = F[i];
            xadj.push_back(pn.size());
            for (size_t j = 0; j < pn.size(); j++) fadj.push_back(pn[j]);
        }
    }

    Int nefaces = (Int)efaces.size();

    for (Int i = 0; i < nefaces; i++) xadj[i+1] += xadj[i];

    // build skin neighbourhood
    std::vector<Int> fneighbours;
    K_CONNECT::build_face_neighbourhood(fadj, xadj, fneighbours);

    // color the faces by connex part
    std::vector<Int> colors(xadj.size()-1);
    Int nconnex = K_CONNECT::colorConnexParts(&fneighbours[0], &xadj[0],
        nefaces, &colors[0]);

    printf("orient_boundary(): connex parts: " SF_D_ "\n", nconnex);

    assert(efaces.size() == xadj.size()-1);
    std::vector<Int> forient(nefaces, 0);
    std::vector<Int> cflags;
    Int ret = 0;

    if (nconnex > 1) {
        // extract nconnex nface-ngon for separate orientation
        for (Int color = 0; color < nconnex; color++) {
            std::vector<bool> keep_pgs(nf, false);
            for (Int i = 0; i < nefaces; i++) {
                keep_pgs[efaces[i]] = (colors[i] == color);
            }
            
            // extract nface corresponding to kept faces
            std::vector<Int> NFACE, cxadj(1, 0), cells;
            extract_nface_of_kept_pgs(keep_pgs, NFACE, cxadj, cells);

            std::vector<Int> cflags;
            flag_marked_external_cells(cells, fflags, cflags);

            ret |= orient_boundary((Int)cells.size(),
                &fadj[0], &xadj[0], nefaces, &fneighbours[0], &efaces[0],
                forient, cflags, fflags, &cells[0], normal_direction);
        }
    } else {
        std::vector<Int> cflags;
        flag_all_external_cells(fflags, cflags);
        ret = orient_boundary(nc, &fadj[0], &xadj[0], nefaces, &fneighbours[0],
            &efaces[0], forient, cflags, fflags, NULL, normal_direction);
    }

    // Apply orientation
    Int nrev = 0;
    for (Int i = 0; i < nefaces; i++) {
        if (forient[i] == -1) {
            Int face = efaces[i]; // 0-based
            auto &pn = F[face];
            std::reverse(pn.begin(), pn.end());
            nrev++;
        }
    }

    printf("orient_boundary(): reversed " SF_D_ " faces\n", nrev);

    return ret;
}

void IMesh::flag_and_get_external_faces(std::vector<Int> &fflags,
    std::vector<Int> &efaces)
{
    std::vector<Int> face_count(nf, 0);

    // Loop through the elements and increment face_count
    for (Int i = 0; i < nc; i++) {
        const auto &pf = C[i];
        for (Int face : pf) face_count[face]++;
    }

    // External faces are those with a count equal to 1
    fflags.resize(nf);

    for (Int i = 0; i < nf; i++) {
        if (face_count[i] == 1) {
            fflags[i] = EXTERNAL;
            efaces.push_back(i);
        } else {
            fflags[i] = INTERNAL;
        }
    }
}

void IMesh::extract_nface_of_kept_pgs(const std::vector<bool> &kept_pgs,
    std::vector<Int> &NFACE, std::vector<Int> &xadj,
    std::vector<Int> &cells)
{
    NFACE.clear();
    xadj.resize(1, 0);
    cells.clear();

    for (Int i = 0; i < nc; i++) {
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

void IMesh::flag_marked_external_cells(const std::vector<Int> &cells,
    const std::vector<Int> &fflags, std::vector<Int> &cflags)
{
    // External cells are those with at least one external face
    cflags.resize(cells.size(), INTERNAL);

    for (size_t i = 0; i < cells.size(); i++) {
        Int cell = cells[i];
        const auto &pf = C[cell];
        for (size_t j = 0; j < pf.size(); j++) {
            Int face = pf[j];
            if (fflags[face] == EXTERNAL) {
                cflags[i] = EXTERNAL;
                break;
            }
        }
    }
}

Int IMesh::orient_boundary(Int ncells, Int *efadj, Int *efxadj,
    Int nefaces, Int *fneis, Int *efaces, std::vector<Int> &forient,
    const std::vector<Int> &cflags, const std::vector<Int> &fflags,
    Int *cells, Int normal_direction)
{
    // Look for a cell whose volume is unambiguously computed
    Float cvol = 0.0;
    Int seed = -1;
    Int refPG = -1;
    Int refIdx = -1;

    while (++seed < ncells) {
        if (cflags[seed] != EXTERNAL) continue;

        Int cid = (cells != NULL) ? cells[seed] : seed;

        const auto &pf = C[cid];
        refPG = -1;
        Int local_idx = -1;
        for (size_t j = 0; j < pf.size(); j++) {
            Int face = pf[j];
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
        for (Int i = 0; i < nefaces; i++) {
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

        if (normal_direction == OUT) forient[refIdx] = -forient[refIdx];

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

void IMesh::compute_cell_volume(Int cell, Float &vol, Int refIdx)
{
    // Orient the faces coherently
    std::vector<Int> NGON;
    std::vector<Int> INDPG(1, 0);
    const auto &pf = C[cell];
    Int stride = (Int)pf.size();

    for (Int i = 0; i < stride; i++) {
        Int face = pf[i];
        const auto &pn = F[face];
        Int np = pn.size();
        INDPG.push_back(np);
        for (Int j = 0; j < np; j++) NGON.push_back(pn[j]);
    }

    for (Int i = 0; i < stride; i++) INDPG[i+1] += INDPG[i];

    // Fix orientation of first face
    std::vector<Int> orient(stride);
    orient[refIdx] = 1;
    std::vector<Int> neis(NGON.size());
    K_CONNECT::build_face_neighbourhood(NGON, INDPG, neis);
    K_CONNECT::reversi_connex(&NGON[0], &INDPG[0], stride, &neis[0], refIdx,
        orient);

    // Apply orientation in local NGON
    for (Int i = 0; i < stride; i++) {
        if (orient[i] == -1) {
            Int start = INDPG[i];
            Int np = INDPG[i+1] - start;
            Int *pn = &NGON[start];
            std::reverse(pn+1, pn+np);
        }
    }

    // Compute faces area and center
    std::vector<Float> faceAreas(3*stride, 0.0);
    std::vector<Float> faceCenters(3*stride, 0.0);

    for (Int i = 0; i < stride; i++) {
        Int face = pf[i];
        Float *fa = &faceAreas[3*i];
        Float *fc = &faceCenters[3*i];
        Int np = INDPG[i+1]-INDPG[i];
        Int *pn = &NGON[INDPG[i]];
        for (Int j = 0; j < np; j++) pn[j] += 1;
        K_METRIC::compute_face_center_and_area(face, np, pn, X.data(), Y.data(),
            Z.data(), fc, fa);
    }
  
    // Estimate cell centroid as average of face centers
    Float cc[3] = {0,0,0};
    for (Int i = 0; i < stride; i++) {
        Float *fc = &faceCenters[3*i];
        for (Int j = 0; j < 3; j++) cc[j] += fc[j];
    }
    for (Int i = 0; i < 3; i++) cc[i] /= stride;

    // Compute cell volume
    vol = 0.0;

    for (Int i = 0; i < stride; i++) {
        Float *fa = &faceAreas[3*i];
        Float *fc = &faceCenters[3*i];

        // Compute 3*face-pyramid volume contribution
        Float d[3] = {fc[0]-cc[0], fc[1]-cc[1], fc[2]-cc[2]};
        //Float pyr3Vol = K_MATH::dot(fa, fc, 3);
        Float pyr3Vol = K_MATH::dot(fa, d, 3);

        vol += pyr3Vol;
    }

    vol *= K_MATH::ONE_THIRD;
}

void IMesh::flag_all_external_cells(const std::vector<Int> &fflags,
    std::vector<Int> &cflags)
{
    // External cells are those with at least one external face
    cflags.resize(nc, INTERNAL);
    for (Int i = 0; i < nc; i++) {
        const auto &pf = C[i];
        Int stride = (Int)pf.size();
        for (Int j = 0; j < stride; j++) {
            Int face = pf[j];
            if (fflags[face] == EXTERNAL) {
                cflags[i] = EXTERNAL;
                break;
            }
        }
    }
}