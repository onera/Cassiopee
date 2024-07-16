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
#include <vector>
#include <algorithm>
#include <stack>

#include "Mesh.h"
#include "../common/mem.h"
#include "Edge.h"

#define INTERNAL 0
#define EXTERNAL 1

static
void flag_and_get_external_faces(Mesh *M, std::vector<Int> &fflags,
    std::vector<Int> &efaces)
{
    fflags.clear();
    fflags.resize(M->nf, 0);

    for (Int i = 0; i < M->nc; i++) {
        Int *cell = Mesh_get_cell(M, i);
        Int *crange = Mesh_get_crange(M, i);
        Int cstride = M->cstride[i];

        for (Int j = 0; j < cstride; j++) {
            Int *pf = cell + 4*j;

            for (Int k = 0; k < crange[j]; k++) {
                fflags[pf[k]]++;
            }
        }
    }

    efaces.clear();

    for (Int i = 0; i < M->nf; i++) {
        assert(fflags[i] == 1 || fflags[i] == 2);
        if (fflags[i] == 1) {
            fflags[i] = EXTERNAL;
            efaces.push_back(i);
        } else {
            fflags[i] = INTERNAL;
        }
    }
}

void Mesh_get_faces_connectivity(Mesh *M, const std::vector<Int> &faces, 
    std::vector<Int> &xadj, std::vector<Int> &fadj)
{
    xadj.clear();
    fadj.clear();

    xadj.resize(faces.size()+1);
    xadj[0] = 0;

    for (size_t i = 0; i < faces.size(); i++) {
        Int fid = faces[i];

        Int *face = Mesh_get_face(M, fid);
        Int *frange = Mesh_get_frange(M, fid);
        Int fstride = M->fstride[fid];

        Int np = 0;
    
        for (Int j = 0; j < fstride; j++) {
            Int *pn = face + 2*j;

            for (Int k = 0; k < frange[j]; k++) {
                Int point = pn[k];
                fadj.push_back(point);
                np++;
            }
        }

        xadj[i+1] = xadj[i] + np;
    }

    assert((size_t)xadj[faces.size()] == fadj.size());
}

typedef std::pair<Int, Int> IntPair;

void build_faces_neighbourhood(const std::vector<Int> &xadj,
    const std::vector<Int> &fpts, std::vector<Int> &fneis)
{
    fneis.clear();
    fneis.resize(fpts.size());

    std::map<UEdge, std::pair<IntPair, IntPair>> EM;

    size_t nf = xadj.size() - 1;

    for (size_t i = 0; i < nf; i++) {
        Int start = xadj[i];
        Int end = xadj[i+1];
        Int stride = end - start;

        const Int *pn = &fpts[start];

        for (Int j = 0; j < stride; j++) {
            Int p = pn[j];
            Int q = pn[(j+1)%stride];

            UEdge E(p, q);

            auto it = EM.find(E);

            if (it == EM.end()) {
                // First time encountering this edge
                EM[E].first = std::make_pair(i, j);
                EM[E].second = std::make_pair(-1, -1);
            } else {
                if (it->second.second.first == -1)
                    it->second.second = std::make_pair(i, j);
                else
                    it->second.second = std::make_pair(E_IDX_NONE, E_IDX_NONE);
            }
        }
    }

    for (const auto &edata : EM) {
        Int pg0 = edata.second.first.first;
        Int n0 = edata.second.first.second;
        Int pg1 = edata.second.second.first;
        Int n1 = edata.second.second.second;

        if (pg1 == -1 || pg1 == E_IDX_NONE)
            continue;

        Int s0 = xadj[pg0];
        Int s1 = xadj[pg1];

        fneis[s0 + n0] = pg1;
        fneis[s1 + n1] = pg0;
    }

    std::map<UEdge, Int> edge_to_count;
    for (size_t i = 0; i < nf; i++) {
        Int start = xadj[i];
        Int end = xadj[i+1];
        Int stride = end-start;
        const Int *pn = &fpts[start];
        for (Int j = 0; j < stride; j++) {
            Int ni = pn[j];
            Int nj = pn[(j+1)%stride];
            UEdge E(ni, nj);
            auto it = edge_to_count.find(E);
            if (it == edge_to_count.end())
                edge_to_count.insert(std::make_pair(E, 1));
            else
                it->second++;
        }
    }

    for (size_t i = 0; i < nf; i++) {
        Int start = xadj[i];
        Int end = xadj[i+1];
        Int stride = end-start;
        const Int *pn = &fpts[start];
        Int *pk = &fneis[start];
        for (Int j = 0; j < stride; j++) {
            Int ni = pn[j];
            Int nj = pn[(j+1)%stride];
            UEdge E(ni, nj);
            if (edge_to_count[E] != 2)
                pk[j] = -1;
        }
    }
}

/*
static
void flag_all_external_cells(Mesh *M, const std::vector<Int> &fflags,
    std::vector<Int> &cflags)
{
    cflags.resize(M->nc, INTERNAL);

    for (Int i = 0; i < M->nc; i++) {
        Int *cell = Mesh_get_cell(M, i);
        Int *crange = Mesh_get_crange(M, i);
        Int cstride = M->cstride[i];

        Int found = 0;

        for (Int j = 0; j < cstride && !found; j++) {
            Int *pf = cell + 4*j;

            for (Int k = 0; k < crange[j]; k++) {
                Int face = pf[k];
                if (fflags[face] == EXTERNAL) {
                    cflags[i] = EXTERNAL;
                    found = 1;
                    break;
                }
            }
        }
    }
}
*/

static
void compute_cell_volume_for_orient(Mesh *M, Int seed, Int ref_face,
    Int ref_idx, Float &cvol)
{
    // Orient the cell faces coherently
    std::vector<Int> NGON, INDPG(1, 0);
    Int *cell = Mesh_get_cell(M, seed);
    Int *crange = Mesh_get_crange(M, seed);

    Int nf = 0;
    

    for (Int i = 0; i < M->cstride[seed]; i++) {
        Int *pf = cell + 4*i;

        for (Int j = 0; j < crange[i]; j++) {
            Int fid = pf[j];
            if (nf == ref_idx) assert(fid == ref_face);

            Int *face = Mesh_get_face(M, fid);
            Int *frange = Mesh_get_frange(M, fid);

            Int np = 0;

            for (Int k = 0; k < M->fstride[fid]; k++) {
                Int *pn = face + 2*k;

                for (Int l = 0; l < frange[k]; l++) {
                    Int point = pn[l]+1;
                    NGON.push_back(point);
                    np++;
                }
            }

            INDPG.push_back(np);
            
            nf++;
        }
    }

    for (Int i = 0; i < nf; i++) INDPG[i+1] += INDPG[i];

    // Make the rest of the faces follow the orientation of ref_face
    std::vector<Int> orient(nf, 1);
    orient[ref_idx] = 1;
    std::vector<Int> neis(NGON.size());
    build_faces_neighbourhood(INDPG, NGON, neis);
    K_CONNECT::reversi_connex(NGON.data(), INDPG.data(), nf, neis.data(),
        ref_idx, orient);

    // Apply orientation in local NGON
    for (Int i = 0; i < nf; i++) {
        Int start = INDPG[i];
        Int np = INDPG[i+1] - start;
        Int *pn = &NGON[start];
        if (orient[i] == -1) {
            std::reverse(pn+1, pn+np);
        }
    }

    std::vector<Float> fareas(3*nf, 0);
    std::vector<Float> fcenters(3*nf, 0);

    Int j = 0;
    for (Int i = 0; i < 24; i++) {
        Int face = cell[i];
        if (face == -1) continue;
        Float *fc = &fcenters[3*j];
        Float *fa = &fareas[3*j];
        Int np = INDPG[j+1] - INDPG[j];
        Int *pn = &NGON[INDPG[j]];
        K_METRIC::compute_face_center_and_area(face, np, pn, M->X, M->Y, M->Z,
            fc, fa);
        j++;
    }

    assert(j == nf);

    // Estimate cell centroid as average of face centers
    Float cc[3] = {0, 0, 0};

    j = 0;
    for (Int i = 0; i < 24; i++) {
        Int face = cell[i];
        if (face == -1) continue;
        Float *fc = &fcenters[3*j];

        for (Int k = 0; k < 3; k++) cc[k] += fc[k];
        j++;
    }
    assert(j == nf);

    for (Int i = 0; i < 3; i++) cc[i] /= nf;

    // Compute cell volume
    cvol = 0;

    for (Int i = 0; i < nf; i++) {
        Float *fa = &fareas[3*i];
        Float *fc = &fcenters[3*i];

        Float d[3] = {fc[0]-cc[0], fc[1]-cc[1], fc[2]-cc[2]};

        Float contrib = K_MATH::dot(fa, d, 3);

        cvol += contrib;
    }

    cvol /= 3.0;
}   

static
Int orient_skin(Mesh *M, Int seed, Int ref_face, Int ref_idx, Int *xadj,
    Int *fpts, Int nefaces, Int *fneis, Int *efaces, std::vector<Int> &forient)
{
    Float cvol = 0;
    compute_cell_volume_for_orient(M, seed, ref_face, ref_idx, cvol);

    if (cvol == 0) {
        assert(0);
        return 1;
    }

    // Find face index in efaces
    Int face_idx = -1;

    for (Int i = 0; i < nefaces; i++) {
        if (efaces[i] == ref_face) {
            face_idx = i;
            break;
        }
    }

    if (face_idx == -1) {
        assert(0);
        return 1;
    }

    forient[face_idx] = (cvol > 0) ? 1 : -1;

    // Propagate orientation
    K_CONNECT::reversi_connex(fpts, xadj, nefaces, fneis, face_idx, forient);

    return 0;
}

static
Int Mesh_orient_skin(Mesh *M)
{
    // Extract external faces
    std::vector<Int> fflags, efaces;
    flag_and_get_external_faces(M, fflags, efaces);

    // Extract external faces connectivity
    std::vector<Int> xadj, fpts;
    Mesh_get_faces_connectivity(M, efaces, xadj, fpts);

    // Build skin neighbourhood
    std::vector<Int> fneis;
    build_faces_neighbourhood(xadj, fpts, fneis);

    // Color the faces by connex part
    std::vector<Int> colors(efaces.size());

    Int nefaces = efaces.size();

    Int nconnex = K_CONNECT::colorConnexParts(fneis.data(), xadj.data(),
        nefaces, colors.data());
    
    M->nconnex = nconnex;

    std::vector<Int> forient(nefaces, 1);
    
    Int ret = 0;

    if (nconnex > 1) {
        for (Int color = 0; color < nconnex; color++) {

            std::vector<Int> kept_faces(M->nf, 0), EFACES;
            
            for (Int i = 0; i < nefaces; i++) {
                if (colors[i] == color) {
                    kept_faces[efaces[i]] = 1;
                }
            }

            // Get the cell seed
            Int cseed = -1;
            Int ref_face = -1;
            Int ref_idx = 0;

            for (Int i = 0; (i < M->nc) && (cseed == -1); i++) {
                Int *cell = Mesh_get_cell(M, i);
                Int *crange = Mesh_get_crange(M, i);

                ref_idx = 0;

                for (Int j = 0; (j < M->cstride[i]) && (cseed == -1); j++) {
                    Int *pf = cell + 4*j;

                    for (Int k = 0; k < crange[j]; k++) {
                        Int face = pf[k];
                        if (kept_faces[face] == 1) {
                            cseed = i;
                            ref_face = face;
                            break;
                        }
                        ref_idx++;
                    }
                }
            }

            if (cseed == -1) {
                merr("Couldn't find cseed.");
                assert(0);
                return 1;
            }

            ret |= orient_skin(M, cseed, ref_face, ref_idx, xadj.data(), fpts.data(),
            nefaces, fneis.data(), efaces.data(), forient);
        }
    } else {
        Int cseed = -1;
        Int ref_face = -1;
        Int ref_idx = 0;

        for (Int i = 0; (i < M->nc) && (cseed == -1); i++) {
            Int *cell = Mesh_get_cell(M, i);
            Int *crange = Mesh_get_crange(M, i);

            ref_idx = 0;

            for (Int j = 0; (j < M->cstride[i]) && (cseed == -1); j++) {
                Int *pf = cell + 4*j;

                for (Int k = 0; k < crange[j]; k++) {
                    Int face = pf[k];
                    if (fflags[face] == EXTERNAL) {
                        cseed = i;
                        ref_face = face;
                        break;
                    }
                    ref_idx++;
                }
            }
        }

        if (cseed == -1) {
            merr("Couldn't find cseed.");
            assert(0);
            return 1;
        }

        ret |= orient_skin(M, cseed, ref_face, ref_idx, xadj.data(), fpts.data(),
            nefaces, fneis.data(), efaces.data(), forient);
    }

    if (ret != 0) {
        assert(0);
        return 1;
    }

    for (Int i = 0; i < nefaces; i++) {
        if (forient[i] == -1) {
            Int fid = efaces[i];
            Mesh_reverse_face_points(M, fid);
        }
    }

    return 0;
}

static
void build_cells_neighbourhood(Mesh *M, std::vector<Int> &xadj,
    std::vector<Int> &cneis)
{
    cneis.resize(24*M->nc, -1);

    std::vector<Int> neigh(M->nf, -1);

    for (Int count = 0; count < 2; count++) {
        for (Int i = 0; i < M->nc; i++) {
            Int *cell = Mesh_get_cell(M, i);
            Int *pn = &cneis[24*i];

            for (Int j = 0; j < 24; j++) {
                Int fid = cell[j];
                if (fid == -1) continue;

                Int &nei = neigh[fid];
                Int &Kn = pn[j]; 

                if (nei != -1 && nei != i) Kn = nei;

                neigh[fid] = i;
            }
        }
    }
}

Int Mesh_build_pe(Mesh *M)
{
    std::vector<Int> xadj, cneis;
    build_cells_neighbourhood(M, xadj, cneis);

    std::vector<Int> exPH(M->nc, E_IDX_NONE);

    //assert(M->owner == NULL);
    //assert(M->neigh == NULL);

    XFREE(M->owner);
    XFREE(M->neigh);

    M->owner = IntArray(M->nf);
    M->neigh = IntArray(M->nf);
    
    memset(M->owner, -1, M->nf * sizeof(Int));
    memset(M->neigh, -1, M->nf * sizeof(Int));

    for (Int i = 0; i < M->nc; i++) {
        Int *cell = Mesh_get_cell(M, i);
        Int *neis = &cneis[24*i];
        
        for (Int j = 0; j < 24; j++) {
            Int fid = cell[j];
            if (fid == -1) continue;
            Int nei = neis[j];
            if (nei == -1) {
                exPH[i] = fid+1;
                break;
            }
        }
    }

    std::vector<Int> processed(M->nc, 0);
    Int nconnex = 0;
    Int seed = 0;
    std::stack<Int> cpool;

    while (1) {
        while ((seed < M->nc) &&
              ((processed[seed] == 1) || (exPH[seed] == E_IDX_NONE)))
            seed++;

        if (seed >= M->nc) break;

        nconnex++;

        cpool.push(seed);

        while (!cpool.empty()) {
            Int cid = cpool.top();
            assert(cid != -1);
            cpool.pop();

            if (processed[cid]) continue;

            processed[cid] = 1;

            std::vector<Int> oids, xpgs(1, 0), pgs;
            
            Int *cell = Mesh_get_cell(M, cid);
            Int *crange = Mesh_get_crange(M, cid);
            Int cstride = M->cstride[cid];
            Int nf = 0;

            for (Int i = 0; i < cstride; i++) {
                Int *pf = cell + 4*i;

                for (Int j = 0; j < crange[i]; j++) {
                    Int fid = pf[j];

                    Int *face = Mesh_get_face(M, fid);
                    Int *frange = Mesh_get_frange(M, fid);
                    Int fstride = M->fstride[fid];

                    Int np = 0;

                    for (Int k = 0; k < fstride; k++) {
                        Int *pn = face + 2*k;

                        for (Int l = 0; l < frange[k]; l++) {
                            Int point = pn[l];
                            pgs.push_back(point);
                            np++;
                        }
                    }

                    xpgs.push_back(np);
                    
                    nf++;
                    oids.push_back(fid+1); // need the sign
                }
            }

            for (Int i = 0; i < nf; i++) xpgs[i+1] += xpgs[i];
            assert((size_t)xpgs[nf] == pgs.size());

            std::vector<Int> fneis;
            build_faces_neighbourhood(xpgs, pgs, fneis);
            
            Int revers = 0;

            // Reference face is the external face
            Int ref_face = exPH[cid];
            assert(ref_face != E_IDX_NONE);
            assert(ref_face != -E_IDX_NONE);
            
            //assert(ref_face != -1);

            if (ref_face < 0) {
                revers = 1;
                ref_face = -ref_face;
            }

            // Find reference face index in oids
            Int ref_idx = -1;
            for (size_t i = 0; i < oids.size(); i++) {
                if (ref_face == oids[i]) {
                    ref_idx = i;
                    break;
                }
            }

            if (ref_idx == -1) {
                assert(0);
                return 1;
            }

            std::vector<Int> orient(nf, 1);
            if (revers) orient[ref_idx] = -1;

            K_CONNECT::reversi_connex(pgs.data(), xpgs.data(), nf,
                fneis.data(), ref_idx, orient);
            
            nf = 0;
            Int *neis = &cneis[24*cid];

            for (Int i = 0; i < 24; i++) {
                Int fid = cell[i];
                if (fid == -1) continue;
                Int nei = neis[i];

                M->owner[fid] = cid;
                M->neigh[fid] = nei;

                Int forient = orient[nf++];

                if (nei == -1) {
                    assert(forient == 1);
                    continue;
                }

                exPH[nei] = -(fid+1);

                if (forient == -1) {
                    std::swap(M->owner[fid], M->neigh[fid]);
                    exPH[nei] = fid+1;
                }

                if (!processed[nei]) cpool.push(nei);
            }
        }
    }

    assert(nconnex == M->nconnex);

    return 0;
}

Int Mesh_set_orientation(Mesh *M)
{
    Int ret = Mesh_orient_skin(M);

    if (ret != 0) return ret;
    
    return Mesh_build_pe(M);
}
