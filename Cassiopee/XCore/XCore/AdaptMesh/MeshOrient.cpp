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
#include <vector>
#include <algorithm>
#include <stack>

#include "Mesh.h"
#include "common/mem.h"
#include "Edge.h"

#define INTERNAL 0
#define EXTERNAL 1

static
void flag_and_get_external_faces(Mesh *M, std::vector<E_Int> &fflags,
    std::vector<E_Int> &efaces)
{
    fflags.clear();
    fflags.resize(M->nf, 0);

    for (E_Int i = 0; i < M->nc; i++) {
        E_Int *cell = Mesh_get_cell(M, i);
        E_Int *crange = Mesh_get_crange(M, i);
        E_Int cstride = M->cstride[i];

        for (E_Int j = 0; j < cstride; j++) {
            E_Int *pf = cell + 4*j;

            for (E_Int k = 0; k < crange[j]; k++) {
                fflags[pf[k]]++;
            }
        }
    }

    efaces.clear();

    for (E_Int i = 0; i < M->nf; i++) {
        assert(fflags[i] == 1 || fflags[i] == 2);
        if (fflags[i] == 1) {
            fflags[i] = EXTERNAL;
            efaces.push_back(i);
        } else {
            fflags[i] = INTERNAL;
        }
    }
}

void Mesh_get_faces_connectivity(Mesh *M, const std::vector<E_Int> &faces, 
    std::vector<E_Int> &xadj, std::vector<E_Int> &fadj)
{
    xadj.clear();
    fadj.clear();

    xadj.resize(faces.size()+1);
    xadj[0] = 0;

    for (size_t i = 0; i < faces.size(); i++) {
        E_Int fid = faces[i];

        E_Int *face = Mesh_get_face(M, fid);
        E_Int *frange = Mesh_get_frange(M, fid);
        E_Int fstride = M->fstride[fid];

        E_Int np = 0;
    
        for (E_Int j = 0; j < fstride; j++) {
            E_Int *pn = face + 2*j;

            for (E_Int k = 0; k < frange[j]; k++) {
                E_Int point = pn[k];
                fadj.push_back(point);
                np++;
            }
        }

        xadj[i+1] = xadj[i] + np;
    }

    assert((size_t)xadj[faces.size()] == fadj.size());
}

typedef std::pair<E_Int, E_Int> E_IntPair;

void build_faces_neighbourhood(const std::vector<E_Int> &xadj,
    const std::vector<E_Int> &fpts, std::vector<E_Int> &fneis)
{
    fneis.clear();
    fneis.resize(fpts.size());

    std::map<UEdge, std::pair<E_IntPair, E_IntPair>> EM;

    size_t nf = xadj.size() - 1;

    for (size_t i = 0; i < nf; i++) {
        E_Int start = xadj[i];
        E_Int end = xadj[i+1];
        E_Int stride = end - start;

        const E_Int *pn = &fpts[start];

        for (E_Int j = 0; j < stride; j++) {
            E_Int p = pn[j];
            E_Int q = pn[(j+1)%stride];

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
        E_Int pg0 = edata.second.first.first;
        E_Int n0 = edata.second.first.second;
        E_Int pg1 = edata.second.second.first;
        E_Int n1 = edata.second.second.second;

        if (pg1 == -1 || pg1 == E_IDX_NONE)
            continue;

        E_Int s0 = xadj[pg0];
        E_Int s1 = xadj[pg1];

        fneis[s0 + n0] = pg1;
        fneis[s1 + n1] = pg0;
    }

    std::map<UEdge, E_Int> edge_to_count;
    for (size_t i = 0; i < nf; i++) {
        E_Int start = xadj[i];
        E_Int end = xadj[i+1];
        E_Int stride = end-start;
        const E_Int *pn = &fpts[start];
        for (E_Int j = 0; j < stride; j++) {
            E_Int ni = pn[j];
            E_Int nj = pn[(j+1)%stride];
            UEdge E(ni, nj);
            auto it = edge_to_count.find(E);
            if (it == edge_to_count.end())
                edge_to_count.insert(std::make_pair(E, 1));
            else
                it->second++;
        }
    }

    for (size_t i = 0; i < nf; i++) {
        E_Int start = xadj[i];
        E_Int end = xadj[i+1];
        E_Int stride = end-start;
        const E_Int *pn = &fpts[start];
        E_Int *pk = &fneis[start];
        for (E_Int j = 0; j < stride; j++) {
            E_Int ni = pn[j];
            E_Int nj = pn[(j+1)%stride];
            UEdge E(ni, nj);
            if (edge_to_count[E] != 2)
                pk[j] = -1;
        }
    }
}

/*
static
void flag_all_external_cells(Mesh *M, const std::vector<E_Int> &fflags,
    std::vector<E_Int> &cflags)
{
    cflags.resize(M->nc, INTERNAL);

    for (E_Int i = 0; i < M->nc; i++) {
        E_Int *cell = Mesh_get_cell(M, i);
        E_Int *crange = Mesh_get_crange(M, i);
        E_Int cstride = M->cstride[i];

        E_Int found = 0;

        for (E_Int j = 0; j < cstride && !found; j++) {
            E_Int *pf = cell + 4*j;

            for (E_Int k = 0; k < crange[j]; k++) {
                E_Int face = pf[k];
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
void compute_cell_volume_for_orient(Mesh *M, E_Int seed, E_Int ref_face,
    E_Int ref_idx, E_Float &cvol)
{
    // Orient the cell faces coherently
    std::vector<E_Int> NGON, INDPG(1, 0);
    E_Int *cell = Mesh_get_cell(M, seed);
    E_Int *crange = Mesh_get_crange(M, seed);

    E_Int nf = 0;
    

    for (E_Int i = 0; i < M->cstride[seed]; i++) {
        E_Int *pf = cell + 4*i;

        for (E_Int j = 0; j < crange[i]; j++) {
            E_Int fid = pf[j];
            if (nf == ref_idx) assert(fid == ref_face);

            E_Int *face = Mesh_get_face(M, fid);
            E_Int *frange = Mesh_get_frange(M, fid);

            E_Int np = 0;

            for (E_Int k = 0; k < M->fstride[fid]; k++) {
                E_Int *pn = face + 2*k;

                for (E_Int l = 0; l < frange[k]; l++) {
                    E_Int point = pn[l]+1;
                    NGON.push_back(point);
                    np++;
                }
            }

            INDPG.push_back(np);
            
            nf++;
        }
    }

    for (E_Int i = 0; i < nf; i++) INDPG[i+1] += INDPG[i];

    // Make the rest of the faces follow the orientation of ref_face
    std::vector<E_Int> orient(nf, 1);
    orient[ref_idx] = 1;
    std::vector<E_Int> neis(NGON.size());
    build_faces_neighbourhood(INDPG, NGON, neis);
    K_CONNECT::reversi_connex(NGON.data(), INDPG.data(), nf, neis.data(),
        ref_idx, orient);

    // Apply orientation in local NGON
    for (E_Int i = 0; i < nf; i++) {
        E_Int start = INDPG[i];
        E_Int np = INDPG[i+1] - start;
        E_Int *pn = &NGON[start];
        if (orient[i] == -1) {
            std::reverse(pn+1, pn+np);
        }
    }

    std::vector<E_Float> fareas(3*nf, 0);
    std::vector<E_Float> fcenters(3*nf, 0);

    E_Int j = 0;
    for (E_Int i = 0; i < 24; i++) {
        E_Int face = cell[i];
        if (face == -1) continue;
        E_Float *fc = &fcenters[3*j];
        E_Float *fa = &fareas[3*j];
        E_Int np = INDPG[j+1] - INDPG[j];
        E_Int *pn = &NGON[INDPG[j]];
        K_METRIC::compute_face_center_and_area(face, np, pn, M->X, M->Y, M->Z,
            fc, fa);
        j++;
    }

    assert(j == nf);

    // Estimate cell centroid as average of face centers
    E_Float cc[3] = {0, 0, 0};

    j = 0;
    for (E_Int i = 0; i < 24; i++) {
        E_Int face = cell[i];
        if (face == -1) continue;
        E_Float *fc = &fcenters[3*j];

        for (E_Int k = 0; k < 3; k++) cc[k] += fc[k];
        j++;
    }
    assert(j == nf);

    for (E_Int i = 0; i < 3; i++) cc[i] /= nf;

    // Compute cell volume
    cvol = 0;

    for (E_Int i = 0; i < nf; i++) {
        E_Float *fa = &fareas[3*i];
        E_Float *fc = &fcenters[3*i];

        E_Float d[3] = {fc[0]-cc[0], fc[1]-cc[1], fc[2]-cc[2]};

        E_Float contrib = K_MATH::dot(fa, d, 3);

        cvol += contrib;
    }

    cvol /= 3.0;
}   

static
E_Int orient_skin(Mesh *M, E_Int seed, E_Int ref_face, E_Int ref_idx, E_Int *xadj,
    E_Int *fpts, E_Int nefaces, E_Int *fneis, E_Int *efaces, std::vector<E_Int> &forient)
{
    E_Float cvol = 0;
    compute_cell_volume_for_orient(M, seed, ref_face, ref_idx, cvol);

    if (cvol == 0) {
        assert(0);
        return 1;
    }

    // Find face index in efaces
    E_Int face_idx = -1;

    for (E_Int i = 0; i < nefaces; i++) {
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
E_Int Mesh_orient_skin(Mesh *M)
{
    // Extract external faces
    std::vector<E_Int> fflags, efaces;
    flag_and_get_external_faces(M, fflags, efaces);

    // Extract external faces connectivity
    std::vector<E_Int> xadj, fpts;
    Mesh_get_faces_connectivity(M, efaces, xadj, fpts);

    // Build skin neighbourhood
    std::vector<E_Int> fneis;
    build_faces_neighbourhood(xadj, fpts, fneis);

    // Color the faces by connex part
    std::vector<E_Int> colors(efaces.size());

    E_Int nefaces = efaces.size();

    E_Int nconnex = K_CONNECT::colorConnexParts(fneis.data(), xadj.data(),
        nefaces, colors.data());
    
    M->nconnex = nconnex;

    std::vector<E_Int> forient(nefaces, 1);
    
    E_Int ret = 0;

    if (nconnex > 1) {
        for (E_Int color = 0; color < nconnex; color++) {

            std::vector<E_Int> kept_faces(M->nf, 0), EFACES;
            
            for (E_Int i = 0; i < nefaces; i++) {
                if (colors[i] == color) {
                    kept_faces[efaces[i]] = 1;
                }
            }

            // Get the cell seed
            E_Int cseed = -1;
            E_Int ref_face = -1;
            E_Int ref_idx = 0;

            for (E_Int i = 0; (i < M->nc) && (cseed == -1); i++) {
                E_Int *cell = Mesh_get_cell(M, i);
                E_Int *crange = Mesh_get_crange(M, i);

                ref_idx = 0;

                for (E_Int j = 0; (j < M->cstride[i]) && (cseed == -1); j++) {
                    E_Int *pf = cell + 4*j;

                    for (E_Int k = 0; k < crange[j]; k++) {
                        E_Int face = pf[k];
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
        E_Int cseed = -1;
        E_Int ref_face = -1;
        E_Int ref_idx = 0;

        for (E_Int i = 0; (i < M->nc) && (cseed == -1); i++) {
            E_Int *cell = Mesh_get_cell(M, i);
            E_Int *crange = Mesh_get_crange(M, i);

            ref_idx = 0;

            for (E_Int j = 0; (j < M->cstride[i]) && (cseed == -1); j++) {
                E_Int *pf = cell + 4*j;

                for (E_Int k = 0; k < crange[j]; k++) {
                    E_Int face = pf[k];
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

    for (E_Int i = 0; i < nefaces; i++) {
        if (forient[i] == -1) {
            E_Int fid = efaces[i];
            Mesh_reverse_face_points(M, fid);
        }
    }

    return 0;
}

static
void build_cells_neighbourhood(Mesh *M, std::vector<E_Int> &xadj,
    std::vector<E_Int> &cneis)
{
    cneis.resize(24*M->nc, -1);

    std::vector<E_Int> neigh(M->nf, -1);

    for (E_Int count = 0; count < 2; count++) {
        for (E_Int i = 0; i < M->nc; i++) {
            E_Int *cell = Mesh_get_cell(M, i);
            E_Int *pn = &cneis[24*i];

            for (E_Int j = 0; j < 24; j++) {
                E_Int fid = cell[j];
                if (fid == -1) continue;

                E_Int &nei = neigh[fid];
                E_Int &Kn = pn[j]; 

                if (nei != -1 && nei != i) Kn = nei;

                neigh[fid] = i;
            }
        }
    }
}

E_Int Mesh_build_pe(Mesh *M)
{
    std::vector<E_Int> xadj, cneis;
    build_cells_neighbourhood(M, xadj, cneis);

    std::vector<E_Int> exPH(M->nc, E_IDX_NONE);

    //assert(M->owner == NULL);
    //assert(M->neigh == NULL);

    XFREE(M->owner);
    XFREE(M->neigh);

    M->owner = IntArray(M->nf);
    M->neigh = IntArray(M->nf);
    
    memset(M->owner, -1, M->nf * sizeof(E_Int));
    memset(M->neigh, -1, M->nf * sizeof(E_Int));

    for (E_Int i = 0; i < M->nc; i++) {
        E_Int *cell = Mesh_get_cell(M, i);
        E_Int *neis = &cneis[24*i];
        
        for (E_Int j = 0; j < 24; j++) {
            E_Int fid = cell[j];
            if (fid == -1) continue;
            E_Int nei = neis[j];
            if (nei == -1) {
                exPH[i] = fid+1;
                break;
            }
        }
    }

    std::vector<E_Int> processed(M->nc, 0);
    E_Int nconnex = 0;
    E_Int seed = 0;
    std::stack<E_Int> cpool;

    while (1) {
        while ((seed < M->nc) &&
              ((processed[seed] == 1) || (exPH[seed] == E_IDX_NONE)))
            seed++;

        if (seed >= M->nc) break;

        nconnex++;

        cpool.push(seed);

        while (!cpool.empty()) {
            E_Int cid = cpool.top();
            assert(cid != -1);
            cpool.pop();

            if (processed[cid]) continue;

            processed[cid] = 1;

            std::vector<E_Int> oids, xpgs(1, 0), pgs;
            
            E_Int *cell = Mesh_get_cell(M, cid);
            E_Int *crange = Mesh_get_crange(M, cid);
            E_Int cstride = M->cstride[cid];
            E_Int nf = 0;

            for (E_Int i = 0; i < cstride; i++) {
                E_Int *pf = cell + 4*i;

                for (E_Int j = 0; j < crange[i]; j++) {
                    E_Int fid = pf[j];

                    E_Int *face = Mesh_get_face(M, fid);
                    E_Int *frange = Mesh_get_frange(M, fid);
                    E_Int fstride = M->fstride[fid];

                    E_Int np = 0;

                    for (E_Int k = 0; k < fstride; k++) {
                        E_Int *pn = face + 2*k;

                        for (E_Int l = 0; l < frange[k]; l++) {
                            E_Int point = pn[l];
                            pgs.push_back(point);
                            np++;
                        }
                    }

                    xpgs.push_back(np);
                    
                    nf++;
                    oids.push_back(fid+1); // need the sign
                }
            }

            for (E_Int i = 0; i < nf; i++) xpgs[i+1] += xpgs[i];
            assert((size_t)xpgs[nf] == pgs.size());

            std::vector<E_Int> fneis;
            build_faces_neighbourhood(xpgs, pgs, fneis);
            
            E_Int revers = 0;

            // Reference face is the external face
            E_Int ref_face = exPH[cid];
            assert(ref_face != E_IDX_NONE);
            assert(ref_face != -E_IDX_NONE);
            
            //assert(ref_face != -1);

            if (ref_face < 0) {
                revers = 1;
                ref_face = -ref_face;
            }

            // Find reference face index in oids
            E_Int ref_idx = -1;
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

            std::vector<E_Int> orient(nf, 1);
            if (revers) orient[ref_idx] = -1;

            K_CONNECT::reversi_connex(pgs.data(), xpgs.data(), nf,
                fneis.data(), ref_idx, orient);
            
            nf = 0;
            E_Int *neis = &cneis[24*cid];

            for (E_Int i = 0; i < 24; i++) {
                E_Int fid = cell[i];
                if (fid == -1) continue;
                E_Int nei = neis[i];

                M->owner[fid] = cid;
                M->neigh[fid] = nei;

                E_Int forient = orient[nf++];

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

E_Int Mesh_set_orientation(Mesh *M)
{
    E_Int ret = Mesh_orient_skin(M);

    if (ret != 0) return ret;
    
    return Mesh_build_pe(M);
}
