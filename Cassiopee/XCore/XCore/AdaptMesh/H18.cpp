/*    
    Copyright 2013-2026 ONERA.

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
#include "Hexa.h"
#include "Mesh.h"

void H18_refine(E_Int hexa, Mesh *M)
{
    // This should be replaced by get_ordered_data...
    H18_reorder(hexa, M);

    E_Int *cell = Mesh_get_cell(M, hexa);

    E_Int FACES[24];
    memcpy(FACES, cell, 24 * sizeof(E_Int));

    E_Int *BOT = FACES;
    E_Int *TOP = FACES + 4;
    E_Int *LFT = FACES + 8;
    E_Int *RGT = FACES + 12;
    E_Int *FRO = FACES + 16;
    E_Int *BCK = FACES + 20;

    E_Int NODES[18];
    for (E_Int i = 0; i < 18; i++) NODES[i] = -1;

    // Local variables
    E_Int fid, *pn, i0, reorient, local[9];

    // BOT

    fid = BOT[0];
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[0]);
    if (reorient) std::swap(local[1], local[3]);
    NODES[0]  = local[0];
    NODES[8]  = local[1];
    NODES[12] = local[2];
    NODES[11] = local[3];

    fid = BOT[1];
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[8], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[0]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[8]);
    assert(local[3] == NODES[12]);
    NODES[1] = local[1];
    NODES[9] = local[2];

    fid = BOT[2];
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[12], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[0]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[12]);
    assert(local[1] == NODES[9]);
    NODES[2] = local[2];
    NODES[10] = local[3];

    fid = BOT[3];
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[11], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[0]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[11]);
    assert(local[1] == NODES[12]);
    assert(local[2] == NODES[10]);
    NODES[3] = local[3];

    // LFT

    fid = LFT[0];
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[0], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[2]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[0]);
    assert(local[1] == NODES[11]);
    NODES[13] = local[2];
    NODES[4] = local[3];

    fid = LFT[1];
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[11], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[2]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[11]);
    assert(local[1] == NODES[3]);
    assert(local[3] == NODES[13]);
    NODES[7] = local[2];

    // RGT

    fid = RGT[0];
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[1], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[3]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[1]);
    assert(local[1] == NODES[9]);
    NODES[14] = local[2];
    NODES[5] = local[3];

    fid = RGT[1];
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[9], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[3]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[9]);
    assert(local[1] == NODES[2]);
    assert(local[3] == NODES[14]);
    NODES[6] = local[2];

    // FRO

    fid = FRO[0];
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[1], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[4]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[1]);
    assert(local[1] == NODES[8]);
    assert(local[3] == NODES[5]);
    NODES[15] = local[2];

    fid = FRO[1];
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[8], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[4]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[8]);
    assert(local[1] == NODES[0]);
    assert(local[2] == NODES[4]);
    assert(local[3] == NODES[15]);

    // BCK

    fid = BCK[0];
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[2], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[5]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[2]);
    assert(local[1] == NODES[10]);
    assert(local[3] == NODES[6]);
    NODES[16] = local[2];

    fid = BCK[1];
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[10], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[5]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[10]);
    assert(local[1] == NODES[3]);
    assert(local[2] == NODES[7]);
    assert(local[3] == NODES[16]);

    // TOP

    fid = TOP[0];
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[1]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[4]);
    assert(local[1] == NODES[15]);
    assert(local[3] == NODES[13]);
    NODES[17] = local[2];

    fid = TOP[1];
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[15], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[1]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[15]);
    assert(local[1] == NODES[5]);
    assert(local[2] == NODES[14]);
    assert(local[3] == NODES[17]);

    fid = TOP[2];
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[17], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[1]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[17]);
    assert(local[1] == NODES[14]);
    assert(local[2] == NODES[6]);
    assert(local[3] == NODES[16]);

    fid = TOP[3];
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[13], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[1]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[13]);
    assert(local[1] == NODES[17]);
    assert(local[2] == NODES[16]);
    assert(local[3] == NODES[7]);

    // Set internal faces in ngon
    E_Int *face = NULL;

    face = Mesh_get_face(M, M->nf);
    face[0] = NODES[8];  face[2] = NODES[12];
    face[4] = NODES[17]; face[6] = NODES[15];

    face = Mesh_get_face(M, M->nf+1);
    face[0] = NODES[9];  face[2] = NODES[12];
    face[4] = NODES[17]; face[6] = NODES[14];

    face = Mesh_get_face(M, M->nf+2);
    face[0] = NODES[12]; face[2] = NODES[10];
    face[4] = NODES[16]; face[6] = NODES[17];

    face = Mesh_get_face(M, M->nf+3);
    face[0] = NODES[12]; face[2] = NODES[11];
    face[4] = NODES[13]; face[6] = NODES[17];

    // Update internal face strides, ranges and states
    for (E_Int i = 0; i < 4; i++) {
        E_Int fid = M->nf + i;
        E_Int *frange = Mesh_get_frange(M, fid);
        for (E_Int j = 0; j < 4; j++) frange[j] = 1;
        M->fstride[fid] = 4;
        M->fref[fid] = FACE_NEW;
    }

    // Update patterns
    M->fpattern[M->nf]   = DIR_X;
    M->fpattern[M->nf+1] = DIR_X;
    M->fpattern[M->nf+2] = DIR_X;
    M->fpattern[M->nf+3] = DIR_X;

    // Assemble children
    E_Int *child = NULL;

    // First child replaces hexa
    child = Mesh_get_cell(M, hexa);
    memset(child, -1, 24*sizeof(E_Int));
    child[0]  = BOT[0]; child[4]  = TOP[0];
    child[8]  = LFT[0]; child[12] = M->nf;
    child[16] = FRO[1]; child[20] = M->nf+3;

    // nc
    child = Mesh_get_cell(M, M->nc);
    child[0]  = BOT[1]; child[4]  = TOP[1];
    child[8]  = M->nf;  child[12] = RGT[0];
    child[16] = FRO[0]; child[20] = M->nf+1;

    // nc+1
    child = Mesh_get_cell(M, M->nc+1);
    child[0]  = BOT[2];  child[4]  = TOP[2];
    child[8]  = M->nf+2; child[12] = RGT[1];
    child[16] = M->nf+1; child[20] = BCK[0];

    // nc+2
    child = Mesh_get_cell(M, M->nc+2);
    child[0]  = BOT[3];  child[4]  = TOP[3];
    child[8]  = LFT[1];  child[12] = M->nf+2;
    child[16] = M->nf+3; child[20] = BCK[1];

    // Fix range and strides
    update_range_and_stride(M, hexa, M->nc, 3);

    // Update adaptation info
    M->clevel[hexa]++;

    for (E_Int i = 0; i < 3; i++) {
        M->clevel[M->nc+i] = M->clevel[hexa];
        M->ctype[M->nc+i] = M->ctype[hexa];
    }

    M->cchildren[hexa] = {hexa, M->nc, M->nc+1, M->nc+2};

    // Set shell faces owns and neis
    update_shell_pe(hexa, M);

    // Set owns and neis of internal faces
    M->owner[M->nf] = hexa;
    M->neigh[M->nf] = M->nc;

    M->owner[M->nf+1] = M->nc;
    M->neigh[M->nf+1] = M->nc+1;

    M->owner[M->nf+2] = M->nc+2;
    M->neigh[M->nf+2] = M->nc+1;

    M->owner[M->nf+3] = hexa;
    M->neigh[M->nf+3] = M->nc+2;

    // Update level/type of internal faces
    for (E_Int i = 0; i < 4; i++) {
        M->flevel[M->nf+i] = M->clevel[hexa];
        M->ftype[M->nf+i] = QUAD;
    }

    assert(check_canon_hexa(hexa, M) == 0);
    for (E_Int i = 0; i < 3; i++) assert(check_canon_hexa(M->nc+i, M) == 0);

    // Increment face/hexa count
    M->nf += 4;
    M->nc += 3;
}

void H18_reorder(E_Int hexa, Mesh *M)
{
    E_Int NODES[18];
    for (E_Int i = 0; i < 18; i++) NODES[i] = -1;

    E_Int local[8], i0;

    E_Int *cell = Mesh_get_cell(M, hexa);
    E_Int *crange = Mesh_get_crange(M, hexa);

    E_Int FACES[24];
    for (E_Int i = 0; i < 24; i++) FACES[i] = cell[i];

    E_Int *BOT = FACES;

    if (crange[0] == 4) {
        E_Int first = 0;

        E_Int fid = BOT[first];

        E_Int *pn = Mesh_get_face(M, fid);

        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        E_Int reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[0]);
        if (reorient) std::swap(local[1], local[3]);

        // Find second, third and fourth sides of BOT

        E_Int second, third, fourth;
        second = third = fourth = first;

        for (E_Int i = 1; i < 4; i++) {
            E_Int side = cell[i];
            E_Int *pn = Mesh_get_face(M, side);
            E_Int common[4] = {0, 0, 0, 0};
            for (E_Int j = 0; j < 4; j++) {
                E_Int point = pn[2*j];
                for (E_Int k = 0; k < 4; k++) {
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

        assert(second != first);
        assert(third != first);
        assert(fourth != first);

        // Fill bot nodes
        NODES[0]  = local[0];
        NODES[8]  = local[1];
        NODES[12] = local[2];
        NODES[11] = local[3];

        // Setup second face
        fid = BOT[second];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[8], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[0]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[8]);
        assert(local[3] == NODES[12]);
        NODES[1] = local[1];
        NODES[9] = local[2];

        // Setup third face
        fid = BOT[third];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[12], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[0]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[12]);
        assert(local[1] == NODES[9]);
        NODES[2]  = local[2];
        NODES[10] = local[3];

        // Setup fourth face
        fid = BOT[fourth];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[11], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[0]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[11]);
        assert(local[1] == NODES[12]);
        assert(local[2] == NODES[10]);
        NODES[3] = local[3];
        
        E_Int tmp[4] = {BOT[first], BOT[second], BOT[third], BOT[fourth]};
        for (E_Int i = 0; i < 4; i++) BOT[i] = tmp[i];
    } else {
        assert(crange[0] == 1);
        assert(BOT[1] == -1);
        assert(BOT[2] == -1);
        assert(BOT[3] == -1);
        E_Int *pn = Mesh_get_face(M, BOT[0]);
        for (E_Int i = 0; i < 8; i++) local[i] = pn[i];
        E_Int reorient = Mesh_get_reorient(M, BOT[0], hexa, normalIn_H[0]);
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

    for (E_Int i = 0; i < 4; i++) cell[i] = BOT[i];

    BOT = cell;
    E_Int *TOP = cell + 4;
    E_Int *LFT = cell + 8;
    E_Int *RGT = cell + 12;
    E_Int *FRO = cell + 16;
    E_Int *BCK = cell + 20;
    
    // Find TOP, LFT, RGT, FRO, BCK

    E_Int tmp_crange[6] = {crange[0], -1, -1, -1, -1, -1};

    for (E_Int i = 1; i < 6; i++) {
        E_Int *SIDE = FACES + 4*i;

        std::set<E_Int> points;

        for (E_Int j = 0; j < crange[i]; j++) {
            E_Int fid = SIDE[j];
            assert(fid != -1);
            E_Int *pn = Mesh_get_face(M, fid);
            for (E_Int k = 0; k < 8; k += 2) points.insert(pn[k]);
        }

        E_Int common[4] = {0, 0, 0, 0};

        for (E_Int point : points) {
            for (E_Int j = 0; j < 4; j++) {
                if (NODES[j] == point) {
                    common[j] = 1;
                    break;
                }
            }
        }

        if      (common[0] && common[3]) {
            tmp_crange[2] = crange[i];
            for (E_Int j = 0; j < 4; j++) LFT[j] = SIDE[j];
        }
        else if (common[1] && common[2]) {
            tmp_crange[3] = crange[i];
            for (E_Int j = 0; j < 4; j++) RGT[j] = SIDE[j];
        }
        else if (common[0] && common[1]) {
            tmp_crange[4] = crange[i];
            for (E_Int j = 0; j < 4; j++) FRO[j] = SIDE[j];
        }
        else if (common[2] && common[3]) {
            tmp_crange[5] = crange[i];
            for (E_Int j = 0; j < 4; j++) BCK[j] = SIDE[j];
        }
        else                             {
            tmp_crange[1] = crange[i];
            for (E_Int j = 0; j < 4; j++) TOP[j] = SIDE[j];
        }
    }

    for (E_Int i = 0; i < 6; i++) assert(tmp_crange[i] != -1);

    for (E_Int i = 0; i < 6; i++) crange[i] = tmp_crange[i];

    // Reorder LFT sides

    if (crange[2] == 2) {
        // First face must share NODES[0]
        
        E_Int first = -1;

        for (E_Int i = 0; i < 2 && first == -1; i++) {
            E_Int fid = LFT[i];
            E_Int *pn = Mesh_get_face(M, fid);
            for (E_Int j = 0; j < 4; j++) {
                E_Int point = pn[2*j];
                if (point == NODES[0]) {
                    first = i;
                    break;
                }
            }
        }

        assert(first != -1);

        E_Int second = (first+1)%2;

        // Setup first face
        E_Int fid = LFT[first];
        E_Int *pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        E_Int i0 = Get_pos(NODES[0], local, 4);
        Right_shift(local, i0, 4);
        E_Int reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[2]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[0]);
        assert(local[1] == NODES[11]);
        NODES[13] = local[2];
        NODES[4] = local[3];

        // Setup second face
        fid = LFT[second];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[11], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[2]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[11]);
        assert(local[1] == NODES[3]);
        assert(local[3] == NODES[13]);
        NODES[7] = local[2];

        E_Int tmp[2] = {LFT[first], LFT[second]};
        for (E_Int i = 0; i < 2; i++) LFT[i] = tmp[i];
    } else {
        assert(crange[2] == 1);
        assert(LFT[1] == -1);
        assert(LFT[2] == -1);
        assert(LFT[3] == -1);
        E_Int *pn = Mesh_get_face(M, LFT[0]);
        for (E_Int i = 0; i < 8; i++) local[i] = pn[i];
        i0 = Get_pos(NODES[0], local, 8);
        Right_shift(local, i0, 8);
        E_Int reorient = Mesh_get_reorient(M, LFT[0], hexa, normalIn_H[2]);
        if (reorient) std::reverse(local+1, local+8);
        assert(local[0] == NODES[0]);
        assert(local[1] == NODES[11]);
        assert(local[2] == NODES[3]);
        NODES[7]  = local[4];
        NODES[13] = local[5];
        NODES[4]  = local[6];
    }

    // Reorder RGT sides

    if (crange[3] == 2) {
        // First face must share NODES[1]
        
        E_Int first = -1;

        for (E_Int i = 0; i < 2 && first == -1; i++) {
            E_Int fid = RGT[i];
            E_Int *pn = Mesh_get_face(M, fid);
            for (E_Int j = 0; j < 4; j++) {
                E_Int point = pn[2*j];
                if (point == NODES[1]) {
                    first = i;
                    break;
                }
            }
        }

        assert(first != -1);

        E_Int second = (first+1)%2;

        // Setup first face
        E_Int fid = RGT[first];
        E_Int *pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        E_Int i0 = Get_pos(NODES[1], local, 4);
        Right_shift(local, i0, 4);
        E_Int reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[3]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[1]);
        assert(local[1] == NODES[9]);
        NODES[14] = local[2];
        NODES[5] = local[3];

        // Setup second face
        fid = RGT[second];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[9], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[3]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[9]);
        assert(local[1] == NODES[2]);
        assert(local[3] == NODES[14]);
        NODES[6] = local[2];

        E_Int tmp[2] = {RGT[first], RGT[second]};
        for (E_Int i = 0; i < 2; i++) RGT[i] = tmp[i];
    } else {
        assert(crange[3] == 1);
        assert(RGT[1] == -1);
        assert(RGT[2] == -1);
        assert(RGT[3] == -1);
        E_Int *pn = Mesh_get_face(M, RGT[0]);
        for (E_Int i = 0; i < 8; i++) local[i] = pn[i];
        i0 = Get_pos(NODES[1], local, 8);
        Right_shift(local, i0, 8);
        E_Int reorient = Mesh_get_reorient(M, RGT[0], hexa, normalIn_H[3]);
        if (reorient) std::reverse(local+1, local+8);
        assert(local[0] == NODES[1]);
        assert(local[1] == NODES[9]);
        assert(local[2] == NODES[2]);
        NODES[6]  = local[4];
        NODES[14] = local[5];
        NODES[5]  = local[6];
    }

    // Reorder FRO sides

    if (crange[4] == 2) {
        // First face must share NODES[1]
        
        E_Int first = -1;

        for (E_Int i = 0; i < 2 && first == -1; i++) {
            E_Int fid = FRO[i];
            E_Int *pn = Mesh_get_face(M, fid);
            for (E_Int j = 0; j < 4; j++) {
                E_Int point = pn[2*j];
                if (point == NODES[1]) {
                    first = i;
                    break;
                }
            }
        }

        assert(first != -1);

        E_Int second = (first+1)%2;

        // Setup first face
        E_Int fid = FRO[first];
        E_Int *pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        E_Int i0 = Get_pos(NODES[1], local, 4);
        Right_shift(local, i0, 4);
        E_Int reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[4]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[1]);
        assert(local[1] == NODES[8]);
        assert(local[3] == NODES[5]);
        NODES[15] = local[2];

        // Setup second face
        fid = FRO[second];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[8], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[4]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[8]);
        assert(local[1] == NODES[0]);
        assert(local[2] == NODES[4]);
        assert(local[3] == NODES[15]);

        E_Int tmp[2] = {FRO[first], FRO[second]};
        for (E_Int i = 0; i < 2; i++) FRO[i] = tmp[i];
    } else {
        assert(crange[4] == 1);
        assert(FRO[1] == -1);
        assert(FRO[2] == -1);
        assert(FRO[3] == -1);
        E_Int *pn = Mesh_get_face(M, FRO[0]);
        for (E_Int i = 0; i < 8; i++) local[i] = pn[i];
        i0 = Get_pos(NODES[1], local, 8);
        Right_shift(local, i0, 8);
        E_Int reorient = Mesh_get_reorient(M, FRO[0], hexa, normalIn_H[4]);
        if (reorient) std::reverse(local+1, local+8);
        assert(local[0] == NODES[1]);
        assert(local[1] == NODES[8]);
        assert(local[2] == NODES[0]);
        assert(local[4] == NODES[4]);
        assert(local[6] == NODES[5]);
        NODES[15] = local[5];
    }

    // Reorder BCK sides

    if (crange[5] == 2) {
        // First face must share NODES[2]

        E_Int nn0, nn1, pn0[8], pn1[8];
        Mesh_get_fpoints(M, BCK[0], nn0, pn0);
        Mesh_get_fpoints(M, BCK[1], nn1, pn1);
        
        E_Int first = -1;

        for (E_Int i = 0; i < 2 && first == -1; i++) {
            E_Int fid = BCK[i];
            E_Int *pn = Mesh_get_face(M, fid);
            for (E_Int j = 0; j < 4; j++) {
                E_Int point = pn[2*j];
                if (point == NODES[2]) {
                    first = i;
                    break;
                }
            }
        }

        assert(first != -1);

        E_Int second = (first+1)%2;

        // Setup first face
        E_Int fid = BCK[first];
        E_Int *pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        E_Int i0 = Get_pos(NODES[2], local, 4);
        Right_shift(local, i0, 4);
        E_Int reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[5]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[2]);
        assert(local[1] == NODES[10]);
        assert(local[3] == NODES[6]);
        NODES[16] = local[2];
        assert(NODES[16] != -1);

        // Setup second face
        fid = BCK[second];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[10], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[5]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[10]);
        assert(local[1] == NODES[3]);
        assert(local[2] == NODES[7]);
        assert(local[3] == NODES[16]);

        E_Int tmp[2] = {BCK[first], BCK[second]};
        for (E_Int i = 0; i < 2; i++) BCK[i] = tmp[i];
    } else {
        assert(crange[5] == 1);
        assert(BCK[1] == -1);
        assert(BCK[2] == -1);
        assert(BCK[3] == -1);
        E_Int *pn = Mesh_get_face(M, BCK[0]);
        for (E_Int i = 0; i < 8; i++) local[i] = pn[i];
        i0 = Get_pos(NODES[2], local, 8);
        Right_shift(local, i0, 8);
        E_Int reorient = Mesh_get_reorient(M, BCK[0], hexa, normalIn_H[5]);
        if (reorient) std::reverse(local+1, local+8);
        assert(local[0] == NODES[2]);
        assert(local[1] == NODES[10]);
        assert(local[2] == NODES[3]);
        assert(local[4] == NODES[7]);
        assert(local[6] == NODES[6]);
        NODES[16] = local[5];
    }

    // Reorder TOP sides

    if (crange[1] == 4) {
        // First face must share NODES[4]

        E_Int first = -1;

        for (E_Int i = 0; i < 4 && first == -1; i++) {
            E_Int fid = TOP[i];
            E_Int *pn = Mesh_get_face(M, fid);
            for (E_Int j = 0; j < 4; j++) {
                E_Int point = pn[2*j];
                if (point == NODES[4]) {
                    first = i;
                    break;
                }
            }
        }

        assert(first != -1);

        // Setup first face
        E_Int fid = TOP[first];
        E_Int *pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        E_Int i0 = Get_pos(NODES[4], local, 4);
        Right_shift(local, i0, 4);
        E_Int reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[1]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[4]);
        assert(local[1] == NODES[15]);
        assert(local[3] == NODES[13]);
        NODES[17] = local[2];

        // Get second, third and fourth sides
        E_Int second, third, fourth;
        second = third = fourth = -1;

        for (E_Int i = 0; i < 4; i++) {
            if (i == first) continue;
            
            E_Int fid = TOP[i];
            E_Int *pn = Mesh_get_face(M, fid);

            E_Int common[4] = {0, 0, 0, 0};

            for (E_Int j = 0; j < 4; j++) {
                E_Int point = pn[2*j];
                for (E_Int k = 0; k < 4; k++) {
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
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[15], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[1]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[15]);
        assert(local[1] == NODES[5]);
        assert(local[2] == NODES[14]);
        assert(local[3] == NODES[17]);

        // Setup third face
        fid = TOP[third];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[17], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[1]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[17]);
        assert(local[1] == NODES[14]);
        assert(local[2] == NODES[6]);
        assert(local[3] == NODES[16]);

        // Setup fourth face
        fid = TOP[fourth];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[13], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[1]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[13]);
        assert(local[1] == NODES[17]);
        assert(local[2] == NODES[16]);
        assert(local[3] == NODES[7]);

        E_Int tmp[4] = {TOP[first], TOP[second], TOP[third], TOP[fourth]};
        for (E_Int i = 0; i < 4; i++) TOP[i] = tmp[i];
    } else {
        assert(crange[1] == 1);
        assert(TOP[1] == -1);
        assert(TOP[2] == -1);
        assert(TOP[3] == -1);
        E_Int *pn = Mesh_get_face(M, TOP[0]);
        for (E_Int i = 0; i < 8; i++) local[i] = pn[i];
        i0 = Get_pos(NODES[4], local, 8);
        Right_shift(local, i0, 8);
        E_Int reorient = Mesh_get_reorient(M, TOP[0], hexa, normalIn_H[1]);
        if (reorient) std::reverse(local+1, local+8);
        assert(local[0] == NODES[4]);
        assert(local[1] == NODES[15]);
        assert(local[2] == NODES[5]);
        assert(local[3] == NODES[14]);
        assert(local[4] == NODES[6]);
        assert(local[5] == NODES[16]);
        assert(local[6] == NODES[7]);
        assert(local[7] == NODES[13]);
    }
}
