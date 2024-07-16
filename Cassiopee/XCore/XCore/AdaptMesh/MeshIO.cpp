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
#include "Mesh.h"
#include "xcore.h"
#include "Hexa.h"

// Assumes HEXA only for now
static
PyObject *export_BE_mesh(Mesh *M)
{
    const char *var_string = "CoordinateX,CoordinateY,CoordinateZ";

    PyObject *karray = K_ARRAY::buildArray3(3, var_string, M->np, M->nc,
        "HEXA", false, 3);
    
    FldArrayF *f;
    FldArrayI *cn;

    K_ARRAY::getFromArray3(karray, f, cn);

    Int *pcn = cn->begin();
    Int *ptr = pcn;

    Int local[9];
    Int NODES[26];
    Int i0;

    for (Int cid = 0; cid < M->nc; cid++) {
        Int *cell = Mesh_get_cell(M, cid);
        Int *crange = Mesh_get_crange(M, cid);

        assert(M->ctype[cid] == HEXA);

        for (Int i = 0; i < 26; i++) NODES[i] = -1;

        // BOT
        Int *BOT = cell;

        if (crange[0] == 4) {
            Int fid = BOT[0];
            Int *pn = Mesh_get_face(M, fid);
            for (Int i = 0; i < 4; i++) local[i] = pn[2*i];
            Int reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[0]);
            if (reorient) std::swap(local[1], local[3]);
            NODES[0]  = local[0];
            NODES[8]  = local[1];
            NODES[12] = local[2];
            NODES[11] = local[3];

            // Setup second face
            fid = BOT[1];
            pn = Mesh_get_face(M, fid);
            for (Int i = 0; i < 4; i++) local[i] = pn[2*i];
            i0 = Get_pos(NODES[8], local, 4);
            Right_shift(local, i0, 4);
            reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[0]);
            if (reorient) std::swap(local[1], local[3]);
            assert(local[0] == NODES[8]);
            assert(local[3] == NODES[12]);
            NODES[1] = local[1];
            NODES[9] = local[2];

            // Setup third face
            fid = BOT[2];
            pn = Mesh_get_face(M, fid);
            for (Int i = 0; i < 4; i++) local[i] = pn[2*i];
            i0 = Get_pos(NODES[12], local, 4);
            Right_shift(local, i0, 4);
            reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[0]);
            if (reorient) std::swap(local[1], local[3]);
            assert(local[0] == NODES[12]);
            assert(local[1] == NODES[9]);
            NODES[2]  = local[2];
            NODES[10] = local[3];

            // Setup fourth face
            fid = BOT[3];
            pn = Mesh_get_face(M, fid);
            for (Int i = 0; i < 4; i++) local[i] = pn[2*i];
            i0 = Get_pos(NODES[11], local, 4);
            Right_shift(local, i0, 4);
            reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[0]);
            if (reorient) std::swap(local[1], local[3]);
            assert(local[0] == NODES[11]);
            assert(local[1] == NODES[12]);
            assert(local[2] == NODES[10]);
            NODES[3] = local[3];
        } else {
            assert(crange[0] == 1);

            assert(BOT[1] == -1);
            assert(BOT[2] == -1);
            assert(BOT[3] == -1);
            Int *pn = Mesh_get_face(M, BOT[0]);
            for (Int i = 0; i < 8; i++) local[i] = pn[i];
            Int reorient = Mesh_get_reorient(M, BOT[0], cid, normalIn_H[0]);
            if (reorient) std::reverse(local+1, local+8);
            NODES[0]  = local[0];
            NODES[1]  = local[2];
            NODES[2]  = local[4];
            NODES[3]  = local[6];
            NODES[8]  = local[1];
            NODES[9]  = local[3];
            NODES[10] = local[5];
            NODES[11] = local[7];
        }

        // Find LFT
        Int top, lft, rgt, fro, bck;
        top = lft = rgt = fro = bck = 0;

        for (Int i = 1; i < 6 && (top == 0 || lft == 0); i++) {
            Int *SIDE = cell + 4*i;

            std::set<Int> points;

            for (Int j = 0; j < crange[i]; j++) {
                Int fid = SIDE[j];
                assert(fid != -1);
                Int *pn = Mesh_get_face(M, fid);
                for (Int k = 0; k < 8; k += 2) points.insert(pn[k]);
            }

            Int common[4] = {0, 0, 0, 0};

            for (Int point : points) {
                for (Int j = 0; j < 4; j++) {
                    if (NODES[j] == point) {
                        common[j] = 1;
                        break;
                    }
                }
            }

            if      (common[0] && common[3]) {
                lft = i;
            }
            else if (common[1] && common[2]) {
                rgt = i;
            }
            else if (common[0] && common[1]) {
                fro = i;
            }
            else if (common[2] && common[3]) {
                bck = i;
            }
            else                             {
                top = i;
            }
        }

        assert(top != 0);
        assert(lft != 0);

        // Rearrange LFT
        Int *LFT = cell + 4*lft;

        if (crange[lft] == 1) {
            assert(LFT[1] == -1);
            assert(LFT[2] == -1);
            assert(LFT[3] == -1);
            Int *pn = Mesh_get_face(M, LFT[0]);
            for (Int i = 0; i < 8; i++) local[i] = pn[i];
            i0 = Get_pos(NODES[0], local, 8);
            Right_shift(local, i0, 8);
            Int reorient = Mesh_get_reorient(M, LFT[0], cid, normalIn_H[lft]);
            if (reorient) std::reverse(local+1, local+8);
            assert(local[0] == NODES[0]);
            assert(local[1] == NODES[11]);
            assert(local[2] == NODES[3]);
            NODES[13] = local[3];
            NODES[7]  = local[4];
            NODES[14] = local[5];
            NODES[4]  = local[6];
            NODES[15] = local[7];
        } else if (crange[lft] == 2) {
            Int first = -1;

            for (Int i = 0; i < 2 && first == -1; i++) {
                Int fid = LFT[i];
                Int *pn = Mesh_get_face(M, fid);
                for (Int j = 0; j < 4; j++) {
                    Int point = pn[2*j];
                    if (point == NODES[0]) {
                        first = i;
                        break;
                    }
                }
            }

            assert(first != -1);

            Int second = (first+1)%2;

            // Setup first face
            Int fid = LFT[first];
            Int *pn = Mesh_get_face(M, fid);
            for (Int i = 0; i < 4; i++) local[i] = pn[2*i];
            Int i0 = Get_pos(NODES[0], local, 4);
            Right_shift(local, i0, 4);
            Int reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[lft]);
            if (reorient) std::swap(local[1], local[3]);
            assert(local[0] == NODES[0]);
            assert(local[1] == NODES[11]);
            NODES[14] = local[2];
            NODES[4] = local[3];

            // Setup second face
            fid = LFT[second];
            pn = Mesh_get_face(M, fid);
            for (Int i = 0; i < 4; i++) local[i] = pn[2*i];
            i0 = Get_pos(NODES[11], local, 4);
            Right_shift(local, i0, 4);
            reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[lft]);
            if (reorient) std::swap(local[1], local[3]);
            assert(local[0] == NODES[11]);
            assert(local[1] == NODES[3]);
            assert(local[3] == NODES[14]);
            NODES[7] = local[2];
        } else {
            assert(crange[lft] == 4);
            Int first = -1;

            for (Int i = 0; i < 4 && first == -1; i++) {
                Int fid = LFT[i];
                Int *pn = Mesh_get_face(M, fid);
                for (Int j = 0; j < 4; j++) {
                    Int point = pn[2*j];
                    if (point == NODES[0]) {
                        first = i;
                        break;
                    }
                }
            }

            assert(first != -1);

            // Setup first face
            Int fid = LFT[first];
            Int *pn = Mesh_get_face(M, fid);
            for (Int i = 0; i < 4; i++) local[i] = pn[2*i];
            Int i0 = Get_pos(NODES[0], local, 4);
            Right_shift(local, i0, 4);
            Int reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[lft]);
            if (reorient) std::swap(local[1], local[3]);
            assert(local[0] == NODES[0]);
            assert(local[1] == NODES[11]);
            NODES[16] = local[2];
            NODES[15] = local[3];

            // Get second, third and fourth sides
            Int second, third, fourth;
            second = third = fourth = -1;

            for (Int i = 0; i < 4; i++) {
                if (i == first) continue;

                Int fid = LFT[i];
                Int *pn = Mesh_get_face(M, fid);
                Int common[4] = {0, 0, 0, 0};

                for (Int j = 0; j < 4; j++) {
                    Int point = pn[2*j];
                    for (Int k = 0; k < 4; k++) {
                        if (local[k] == point) {
                            common[k] = 1;
                            break;
                        }
                    }
                }

                if (common[1] && common[2]) second = i;
                else if (common[2] && common[3]) fourth = i;
                else third = i;
            }

            assert(second != -1);
            assert(third != -1);
            assert(fourth != -1);

            // Setup second face
            fid = LFT[second];
            pn = Mesh_get_face(M, fid);
            for (Int i = 0; i < 4; i++) local[i] = pn[2*i];
            i0 = Get_pos(NODES[11], local, 4);
            Right_shift(local, i0, 4);
            reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[lft]);
            if (reorient) std::swap(local[1], local[3]);
            assert(local[0] == NODES[11]);
            assert(local[1] == NODES[3]);
            assert(local[3] == NODES[16]);
            NODES[13] = local[2];

            // Setup third face
            fid = LFT[third];
            pn = Mesh_get_face(M, fid);
            for (Int i = 0; i < 4; i++) local[i] = pn[2*i];
            i0 = Get_pos(NODES[16], local, 4);
            Right_shift(local, i0, 4);
            reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[lft]);
            if (reorient) std::swap(local[1], local[3]);
            assert(local[0] == NODES[16]);
            assert(local[1] == NODES[13]);
            NODES[7] = local[2];
            NODES[14] = local[3];

            // Setup fourth face
            fid = LFT[fourth];
            pn = Mesh_get_face(M, fid);
            for (Int i = 0; i < 4; i++) local[i] = pn[2*i];
            i0 = Get_pos(NODES[15], local, 4);
            Right_shift(local, i0, 4);
            reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[lft]);
            if (reorient) std::swap(local[1], local[3]);
            assert(local[0] == NODES[15]);
            assert(local[1] == NODES[16]);
            assert(local[2] == NODES[14]);
            NODES[4] = local[3];
        }

        // Rearrange TOP
        Int *TOP = cell + 4*top;

        if (crange[top] == 4) {
            // First face must share NODES[4]

            Int first = -1;

            for (Int i = 0; i < 4 && first == -1; i++) {
                Int fid = TOP[i];
                Int *pn = Mesh_get_face(M, fid);
                for (Int j = 0; j < 4; j++) {
                    Int point = pn[2*j];
                    if (point == NODES[4]) {
                        first = i;
                        break;
                    }
                }
            }

            assert(first != -1);

            // Setup first face
            Int fid = TOP[first];
            Int *pn = Mesh_get_face(M, fid);
            for (Int i = 0; i < 4; i++) local[i] = pn[2*i];
            Int i0 = Get_pos(NODES[4], local, 4);
            Right_shift(local, i0, 4);
            Int reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[top]);
            if (reorient) std::swap(local[1], local[3]);
            assert(local[0] == NODES[4]);
            NODES[21] = local[1];
            NODES[14] = local[3];
            NODES[25] = local[2];

            // Get second, third and fourth sides
            Int second, third, fourth;
            second = third = fourth = -1;

            for (Int i = 0; i < 4; i++) {
                if (i == first) continue;

                Int fid = TOP[i];
                Int *pn = Mesh_get_face(M, fid);

                Int common[4] = {0, 0, 0, 0};

                for (Int j = 0; j < 4; j++) {
                    Int point = pn[2*j];
                    for (Int k = 0; k < 4; k++) {
                        if (local[k] == point) {
                            common[k] = 1;
                            break;
                        }
                    }
                }

                if (common[1] && common[2]) second = i;
                else if (common[2] && common[3]) fourth = i;
                else third = i;
            }

            assert(second != -1);
            assert(third != -1);
            assert(fourth != -1);

            // Setup second face
            fid = TOP[second];
            pn = Mesh_get_face(M, fid);
            for (Int i = 0; i < 4; i++) local[i] = pn[2*i];
            i0 = Get_pos(NODES[21], local, 4);
            Right_shift(local, i0, 4);
            reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[top]);
            if (reorient) std::swap(local[1], local[3]);
            assert(local[0] == NODES[21]);
            assert(local[1] == NODES[5]);
            assert(local[2] == NODES[18]);
            assert(local[3] == NODES[25]);

            // Setup third face
            fid = TOP[third];
            pn = Mesh_get_face(M, fid);
            for (Int i = 0; i < 4; i++) local[i] = pn[2*i];
            i0 = Get_pos(NODES[25], local, 4);
            Right_shift(local, i0, 4);
            reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[top]);
            if (reorient) std::swap(local[1], local[3]);
            assert(local[0] == NODES[25]);
            assert(local[1] == NODES[18]);
            assert(local[2] == NODES[6]);
            assert(local[3] == NODES[23]);

            // Setup fourth face
            fid = TOP[fourth];
            pn = Mesh_get_face(M, fid);
            for (Int i = 0; i < 4; i++) local[i] = pn[2*i];
            i0 = Get_pos(NODES[14], local, 4);
            Right_shift(local, i0, 4);
            reorient = Mesh_get_reorient(M, fid, cid, normalIn_H[top]);
            if (reorient) std::swap(local[1], local[3]);
            assert(local[0] == NODES[14]);
            assert(local[1] == NODES[25]);
            assert(local[2] == NODES[23]);
            assert(local[3] == NODES[7]);
        } else {
            assert(crange[top] == 1);
            assert(TOP[1] == -1);
            assert(TOP[2] == -1);
            assert(TOP[3] == -1);
            Int *pn = Mesh_get_face(M, TOP[0]);
            for (Int i = 0; i < 8; i++) local[i] = pn[i];
            i0 = Get_pos(NODES[4], local, 8);
            Right_shift(local, i0, 8);
            Int reorient = Mesh_get_reorient(M, TOP[0], cid, normalIn_H[top]);
            if (reorient) std::reverse(local+1, local+8);
            assert(local[0] == NODES[4]);
            assert(local[6] == NODES[7]);
            assert(local[7] == NODES[14]);
            NODES[5] = local[2];
            NODES[6] = local[4];
        }

        // One-based
        for (Int j = 0; j < 8; j++) *ptr++ = NODES[j] + 1;
    }

    Float *px = f->begin(1);
    Float *py = f->begin(2);
    Float *pz = f->begin(3);

    memcpy(px, M->X, M->np * sizeof(Float));
    memcpy(py, M->Y, M->np * sizeof(Float));
    memcpy(pz, M->Z, M->np * sizeof(Float));

    PyObject *out = PyList_New(0);

    PyList_Append(out, karray);
    Py_DECREF(karray);

    // BCs
    PyList_Append(out, Py_None);


    // Comm is undefined
    PyList_Append(out, Py_None);

    return out;
}

static
PyObject *export_conformal_mesh(Mesh *M)
{
    Int nfld = 3;
    const char *var_string = "CoordinateX,CoordinateY,CoordinateZ";
    const char *elt_type = "NGON";
    Int ngon_type = 3;
    bool center = false;
    Int api = 3;

    Int sizeNFace = Mesh_get_sizeNFace(M);
    Int sizeNGon = Mesh_get_sizeNGon(M);

    PyObject *karray = K_ARRAY::buildArray3(nfld, var_string, M->np, M->nc,
        M->nf, elt_type, sizeNGon, sizeNFace, ngon_type, center,
        api);

    K_FLD::FldArrayF *f;
    K_FLD::FldArrayI *cn;
    K_ARRAY::getFromArray3(karray, f, cn);

    E_Float *px = f->begin(1);
    E_Float *py = f->begin(2);
    E_Float *pz = f->begin(3);

    memcpy(px, M->X, M->np * sizeof(Float));
    memcpy(py, M->Y, M->np * sizeof(Float));
    memcpy(pz, M->Z, M->np * sizeof(Float));

    Int *indpg = cn->getIndPG();
    Int *ngon = cn->getNGon();
    Int *indph = cn->getIndPH();
    Int *nface = cn->getNFace();

    // NGON

    memset(indpg, 0, (M->nf + 1) * sizeof(Int));

    for (Int i = 0; i < M->nf; i++) {
        Int *frange = Mesh_get_frange(M, i);
        for (Int j = 0; j < M->fstride[i]; j++) {
            indpg[i+1] += frange[j];
        }
    }

    indpg[0] = 0;
    for (Int i = 0; i < M->nf; i++) indpg[i+1] += indpg[i];

    assert(indpg[M->nf] == sizeNGon);
    
    Int *ptr = ngon;

    for (Int i = 0; i < M->nf; i++) {
        Int *face = Mesh_get_face(M, i);
        Int *frange = Mesh_get_frange(M, i);

        for (Int j = 0; j < M->fstride[i]; j++) {
            Int *pn = face + 2*j;

            for (Int k = 0; k < frange[j]; k++) *ptr++ = pn[k] + 1;
        }
    }


    // NFACE

    memset(indph, 0, (M->nc + 1) * sizeof(Int));

    for (Int i = 0; i < M->nc; i++) {
        Int *crange = Mesh_get_crange(M, i);
        for (Int j = 0; j < M->cstride[i]; j++) {
            indph[i+1] += crange[j];
        }
    }

    indph[0] = 0;
    for (Int i = 0; i < M->nc; i++) indph[i+1] += indph[i];

    assert(indph[M->nc] == sizeNFace);
    
    ptr = nface;

    for (Int i = 0; i < M->nc; i++) {
        Int *cell = Mesh_get_cell(M, i);
        Int *crange = Mesh_get_crange(M, i);

        for (Int j = 0; j < M->cstride[i]; j++) {
            Int *pf = cell + 4*j;

            for (Int k = 0; k < crange[j]; k++) *ptr++ = pf[k] + 1;
        }
    }

    PyObject *out = PyList_New(0);

    PyList_Append(out, karray);
    Py_DECREF(karray);

    // BCs

    PyObject *bcs = PyList_New(0);

    npy_intp dims[2];
    dims[1] = 1;

    for (Int i = 0; i < M->nbp; i++) {
        BPatch *P = &M->bps[i];

        dims[0] = P->nf;

        PyArrayObject *facelist = (PyArrayObject *)
            PyArray_SimpleNew(1, dims, E_NPY_INT);
        Int *pf = (Int *)PyArray_DATA(facelist);
        
        for (Int j = 0; j < P->nf; j++) pf[j] = P->pf[j] + 1;
        
        PyObject *bc = Py_BuildValue("[Olss]", facelist, P->gid, P->name, P->type);
        Py_DECREF(facelist);
        
        PyList_Append(bcs, bc);
        Py_DECREF(bc);
    }

    PyList_Append(out, bcs);
    Py_DECREF(bcs);

    // Comm patches

    PyObject *pps = PyList_New(0);

    for (Int i = 0; i < M->npp; i++) {
        PPatch *P = &M->pps[i];

        dims[0] = P->nf;

        PyArrayObject *facelist = (PyArrayObject *)
            PyArray_SimpleNew(1, dims, E_NPY_INT);
        Int *pf = (Int *)PyArray_DATA(facelist); 
    
        PyArrayObject *neilist = (PyArrayObject *)
            PyArray_SimpleNew(1, dims, E_NPY_INT);
        Int *pn = (Int *)PyArray_DATA(neilist);

        for (Int j = 0; j < P->nf; j++) pf[j] = P->pf[j] + 1;
        
        for (Int j = 0; j < P->nf; j++) pn[j] = P->pn[j];

        PyObject *patch = Py_BuildValue("[lOO]", P->nei, facelist, neilist);
        Py_DECREF(facelist);
        Py_DECREF(neilist);

        PyList_Append(pps, patch);
        Py_DECREF(patch);
    }

    PyList_Append(out, pps);
    Py_DECREF(pps);

    return out;
}

PyObject *Mesh_export_karray(Mesh *M, int conformize)
{
    if (conformize) return export_conformal_mesh(M);

    return export_BE_mesh(M);
}