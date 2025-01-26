#include "Hexa.h"
#include "Quad.h"
#include "Mesh.h"

void H27_refine(E_Int hexa, Mesh *M)
{
    H27_reorder(hexa, M);

    E_Int *cell = Mesh_get_cell(M, hexa);

    E_Int FACES[24];

    memcpy(FACES, cell, 24 * sizeof(E_Int));

    E_Int *BOT = FACES;
    E_Int *TOP = FACES + 4;
    E_Int *LFT = FACES + 8;
    E_Int *RGT = FACES + 12;
    E_Int *FRO = FACES + 16;
    E_Int *BCK = FACES + 20;

    E_Int NODES[27];
    for (E_Int i = 0; i < 27; i++) NODES[i] = -1;

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
    NODES[16] = local[2];
    NODES[15] = local[3];

    fid = LFT[1]; 
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[11], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[2]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[11]);
    assert(local[1] == NODES[3]);
    assert(local[3] == NODES[16]);
    NODES[13] = local[2];

    fid = LFT[2]; 
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[16], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[2]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[16]);
    assert(local[1] == NODES[13]);
    NODES[7] = local[2];
    NODES[14] = local[3];

    fid = LFT[3]; 
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[15], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[2]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[15]);
    assert(local[1] == NODES[16]);
    assert(local[2] == NODES[14]);
    NODES[4] = local[3];

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
    NODES[20] = local[2];
    NODES[19] = local[3];

    fid = RGT[1]; 
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[9], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[3]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[9]);
    assert(local[1] == NODES[2]);
    assert(local[3] == NODES[20]);
    NODES[17] = local[2];

    fid = RGT[2]; 
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[20], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[3]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[20]);
    assert(local[1] == NODES[17]);
    NODES[6] = local[2];
    NODES[18] = local[3];

    fid = RGT[3]; 
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[19], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[3]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[19]);
    assert(local[1] == NODES[20]);
    assert(local[2] == NODES[18]);
    NODES[5] = local[3];

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
    assert(local[3] == NODES[19]);
    NODES[22] = local[2];

    fid = FRO[1]; 
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[8], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[4]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[8]);
    assert(local[1] == NODES[0]);
    assert(local[2] == NODES[15]);
    assert(local[3] == NODES[22]);

    fid = FRO[2]; 
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[22], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[4]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[22]);
    assert(local[1] == NODES[15]);
    assert(local[2] == NODES[4]);
    NODES[21] = local[3];

    fid = FRO[3]; 
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[19], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[4]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[19]);
    assert(local[1] == NODES[22]);
    assert(local[2] == NODES[21]);
    assert(local[3] == NODES[5]);

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
    assert(local[3] == NODES[17]);
    NODES[24] = local[2];

    fid = BCK[1]; 
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[10], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[5]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[10]);
    assert(local[1] == NODES[3]);
    assert(local[2] == NODES[13]);
    assert(local[3] == NODES[24]);

    fid = BCK[2]; 
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[24], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[5]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[24]);
    assert(local[1] == NODES[13]);
    assert(local[2] == NODES[7]);
    NODES[23] = local[3];

    fid = BCK[3]; 
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[17], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[5]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[17]);
    assert(local[1] == NODES[24]);
    assert(local[2] == NODES[23]);
    assert(local[3] == NODES[6]);

    // TOP

    fid = TOP[0]; 
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[4], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[1]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[4]);
    assert(local[1] == NODES[21]);
    assert(local[3] == NODES[14]);
    NODES[25] = local[2];

    fid = TOP[1]; 
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[21], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[1]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[21]);
    assert(local[1] == NODES[5]);
    assert(local[2] == NODES[18]);
    assert(local[3] == NODES[25]);

    fid = TOP[2]; 
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[25], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[1]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[25]);
    assert(local[1] == NODES[18]);
    assert(local[2] == NODES[6]);
    assert(local[3] == NODES[23]);

    fid = TOP[3]; 
    pn = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    i0 = Get_pos(NODES[14], local, 4);
    Right_shift(local, i0, 4);
    reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[1]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[0] == NODES[14]);
    assert(local[1] == NODES[25]);
    assert(local[2] == NODES[23]);
    assert(local[3] == NODES[7]);

 
    // Add hexa centroid
    // TODO(Imad): for now, mean of hexa points xyz
    NODES[26] = M->np;
    M->X[M->np] = M->Y[M->np] = M->Z[M->np] = 0;
    for (E_Int i = 0; i < 8; i++) {
        M->X[M->np] += M->X[NODES[i]];
        M->Y[M->np] += M->Y[NODES[i]];
        M->Z[M->np] += M->Z[NODES[i]];
    }
    M->X[M->np] *= 0.125;
    M->Y[M->np] *= 0.125;
    M->Z[M->np] *= 0.125;
    M->np++;

    // Set internal faces in ngon
    E_Int *face = NULL;

    // nfaces (hexa TOP)
    face = Mesh_get_face(M, M->nf);
    face[0] = NODES[15]; face[2] = NODES[22];
    face[4] = NODES[26]; face[6] = NODES[16];

    // nfaces+1 (hexa RGT)
    face = Mesh_get_face(M, M->nf+1);
    face[0] = NODES[8];  face[2] = NODES[12];
    face[4] = NODES[26]; face[6] = NODES[22];

    // nfaces+2 (hexa BCK)
    face = Mesh_get_face(M, M->nf+2);
    face[0] = NODES[12]; face[2] = NODES[11];
    face[4] = NODES[16]; face[6] = NODES[26];

    // nfaces+3 (nhexas TOP)
    face = Mesh_get_face(M, M->nf+3);
    face[0] = NODES[22]; face[2] = NODES[19];
    face[4] = NODES[20]; face[6] = NODES[26];

    // nfaces+4 (nhexas BCK)
    face = Mesh_get_face(M, M->nf+4);
    face[0] = NODES[9];  face[2] = NODES[12];
    face[4] = NODES[26]; face[6] = NODES[20];

    // nfaces+5 (nhexas+1 TOP)
    face = Mesh_get_face(M, M->nf+5);
    face[0] = NODES[26]; face[2] = NODES[20];
    face[4] = NODES[17]; face[6] = NODES[24];

    // nfaces+6 (nhexas+1 LFT)
    face = Mesh_get_face(M, M->nf+6);
    face[0] = NODES[12]; face[2] = NODES[10];
    face[4] = NODES[24]; face[6] = NODES[26];

    // nfaces+7 (nhexas+2 TOP)
    face = Mesh_get_face(M, M->nf+7);
    face[0] = NODES[16]; face[2] = NODES[26];
    face[4] = NODES[24]; face[6] = NODES[13];

    /*************/

    // nfaces+8 (nhexas+3 RGT)
    face = Mesh_get_face(M, M->nf+8);
    face[0] = NODES[22]; face[2] = NODES[26];
    face[4] = NODES[25]; face[6] = NODES[21];

    // nfaces+9 (nhexas+3 BCK)
    face = Mesh_get_face(M, M->nf+9);
    face[0] = NODES[26]; face[2] = NODES[16];
    face[4] = NODES[14]; face[6] = NODES[25];

    // nfaces+10 (nhexas+4 BCK)
    face = Mesh_get_face(M, M->nf+10);
    face[0] = NODES[20]; face[2] = NODES[26];
    face[4] = NODES[25]; face[6] = NODES[18];

    // nfaces+11 (nhexas+5 LFT)
    face = Mesh_get_face(M, M->nf+11);
    face[0] = NODES[26]; face[2] = NODES[24];
    face[4] = NODES[23]; face[6] = NODES[25];

    // Update internal face strides, ranges and states
    for (E_Int i = 0; i < 12; i++) {
        E_Int fid = M->nf + i;
        E_Int *frange = Mesh_get_frange(M, fid);
        for (E_Int j = 0; j < 4; j++) frange[j] = 1;
        M->fstride[fid] = 4;
        M->fref[fid] = FACE_NEW;
    }

    // Assemble children

    E_Int *child = NULL;

    // First child replaces hexa
    child = Mesh_get_cell(M, hexa);
    memset(child, -1, 24 * sizeof(E_Int));
    child[0]  = BOT[0]; child[4]  = M->nf;
    child[8]  = LFT[0]; child[12] = M->nf+1;
    child[16] = FRO[1]; child[20] = M->nf+2;

    // nhexas
    child = Mesh_get_cell(M, M->nc);
    child[0]  = BOT[1];  child[4]  = M->nf+3;
    child[8]  = M->nf+1; child[12] = RGT[0];
    child[16] = FRO[0];  child[20] = M->nf+4;

    // nhexas+1
    child = Mesh_get_cell(M, M->nc+1);
    child[0]  = BOT[2];  child[4]  = M->nf+5;
    child[8]  = M->nf+6; child[12] = RGT[1];
    child[16] = M->nf+4; child[20] = BCK[0];

    // nhexas+2
    child = Mesh_get_cell(M, M->nc+2);
    child[0] = BOT[3];   child[4]  = M->nf+7;
    child[8] = LFT[1];   child[12] = M->nf+6;
    child[16] = M->nf+2; child[20] = BCK[1];

    /*********/

    // nhexas+3
    child = Mesh_get_cell(M, M->nc+3);
    child[0] = M->nf;   child[4]  = TOP[0];
    child[8] = LFT[3];  child[12] = M->nf+8;
    child[16] = FRO[2]; child[20] = M->nf+9;

    // nhexas+4
    child = Mesh_get_cell(M, M->nc+4);
    child[0]  = M->nf+3;  child[4]  = TOP[1];
    child[8]  = M->nf+8;  child[12] = RGT[3];
    child[16] = FRO[3];   child[20] = M->nf+10;

    // nhexas+5
    child = Mesh_get_cell(M, M->nc+5);
    child[0]  = M->nf+5;  child[4]  = TOP[2];
    child[8]  = M->nf+11; child[12] = RGT[2];
    child[16] = M->nf+10; child[20] = BCK[3];

    // nhexas+6
    child = Mesh_get_cell(M, M->nc+6);
    child[0]  = M->nf+7; child[4]  = TOP[3];
    child[8]  = LFT[2];  child[12] = M->nf+11;
    child[16] = M->nf+9; child[20] = BCK[2];

    // Fix range and strides
    update_range_and_stride(M, hexa, M->nc, 7);

    // Update adaptation info
    M->clevel[hexa]++;

    for (E_Int i = 0; i < 7; i++) {
        M->clevel[M->nc+i] = M->clevel[hexa];
        M->ctype[M->nc+i] = M->ctype[hexa];
    }

    M->cchildren[hexa] = {hexa, M->nc, M->nc+1, M->nc+2, M->nc+3, M->nc+4,
                          M->nc+5, M->nc+6};
    
    // Set shell faces owns and neis
    update_shell_pe(hexa, M);

    // Set owns and neis of internal faces
    M->owner[M->nf]    = hexa;
    M->neigh[M->nf]    = M->nc+3;

    M->owner[M->nf+1]  = hexa;
    M->neigh[M->nf+1]  = M->nc;

    M->owner[M->nf+2]  = hexa;
    M->neigh[M->nf+2]  = M->nc+2; 

    M->owner[M->nf+3]  = M->nc;
    M->neigh[M->nf+3]  = M->nc+4;

    M->owner[M->nf+4]  = M->nc;
    M->neigh[M->nf+4]  = M->nc+1;

    M->owner[M->nf+5]  = M->nc+1;
    M->neigh[M->nf+5]  = M->nc+5;

    M->owner[M->nf+6]  = M->nc+2;        
    M->neigh[M->nf+6]  = M->nc+1; 

    M->owner[M->nf+7]  = M->nc+2;
    M->neigh[M->nf+7]  = M->nc+6;

    M->owner[M->nf+8]  = M->nc+3;
    M->neigh[M->nf+8]  = M->nc+4; 

    M->owner[M->nf+9]  = M->nc+3;
    M->neigh[M->nf+9]  = M->nc+6;

    M->owner[M->nf+10] = M->nc+4;
    M->neigh[M->nf+10] = M->nc+5; 

    M->owner[M->nf+11] = M->nc+6;
    M->neigh[M->nf+11] = M->nc+5;

    // Update level/type of internal faces
    for (E_Int i = 0; i < 12; i++) {
        M->flevel[M->nf+i] = M->clevel[hexa];
        M->ftype[M->nf+i] = QUAD;
    }

    assert(check_canon_hexa(hexa, M) == 0);
    for (E_Int i = 0; i < 7; i++) assert(check_canon_hexa(M->nc+i, M) == 0);

    // Increment face/hexa count
    M->nf += 12;
    M->nc += 7;
}

void H27_reorder(E_Int hexa, Mesh *M)
{
    E_Int NODES[26];
    memset(NODES, -1, 26*sizeof(E_Int));

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

    if (crange[2] == 4) {
        // First face must share NODES[0]

        E_Int first = -1;

        for (E_Int i = 0; i < 4 && first == -1; i++) {
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
        NODES[16] = local[2];
        NODES[15] = local[3];

        // Get second, third and fourth sides
        E_Int second, third, fourth;
        second = third = fourth = -1;

        for (E_Int i = 0; i < 4; i++) {
            if (i == first) continue;

            E_Int fid = LFT[i];
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
        fid = LFT[second];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[11], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[2]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[11]);
        assert(local[1] == NODES[3]);
        assert(local[3] == NODES[16]);
        NODES[13] = local[2];

        // Setup third face
        fid = LFT[third];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[16], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[2]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[16]);
        assert(local[1] == NODES[13]);
        NODES[7] = local[2];
        NODES[14] = local[3];

        // Setup fourth face
        fid = LFT[fourth];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[15], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[2]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[15]);
        assert(local[1] == NODES[16]);
        assert(local[2] == NODES[14]);
        NODES[4] = local[3];

        E_Int tmp[4] = {LFT[first], LFT[second], LFT[third], LFT[fourth]};
        for (E_Int i = 0; i < 4; i++) LFT[i] = tmp[i];

    } else {
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
        NODES[13] = local[3];
        NODES[7]  = local[4];
        NODES[14] = local[5];
        NODES[4]  = local[6];
        NODES[15] = local[7];
    }

    // Reorder RGT sides

    if (crange[3] == 4) {
        // First face must share NODES[1]

        E_Int first = -1;

        for (E_Int i = 0; i < 4 && first == -1; i++) {
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
        NODES[20] = local[2];
        NODES[19] = local[3];

        // Get second, third and fourth sides
        E_Int second, third, fourth;
        second = third = fourth = -1;

        for (E_Int i = 0; i < 4; i++) {
            if (i == first) continue;

            E_Int fid = RGT[i];
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
        fid = RGT[second];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[9], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[3]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[9]);
        assert(local[1] == NODES[2]);
        assert(local[3] == NODES[20]);
        NODES[17] = local[2];

        // Setup third face
        fid = RGT[third];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[20], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[3]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[20]);
        assert(local[1] == NODES[17]);
        NODES[6] = local[2];
        NODES[18] = local[3];

        // Setup fourth face
        fid = RGT[fourth];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[19], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[3]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[19]);
        assert(local[1] == NODES[20]);
        assert(local[2] == NODES[18]);
        NODES[5] = local[3];

        E_Int tmp[4] = {RGT[first], RGT[second], RGT[third], RGT[fourth]};
        for (E_Int i = 0; i < 4; i++) RGT[i] = tmp[i];

    } else {
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
        NODES[17] = local[3];
        NODES[6]  = local[4];
        NODES[18] = local[5];
        NODES[5]  = local[6];
        NODES[19] = local[7];
    }

    // Reorder FRO sides

    if (crange[4] == 4) {
        // First face must share NODES[1]

        E_Int first = -1;

        for (E_Int i = 0; i < 4 && first == -1; i++) {
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
        assert(local[3] == NODES[19]);
        NODES[22] = local[2];

        // Get second, third and fourth sides
        E_Int second, third, fourth;
        second = third = fourth = -1;

        for (E_Int i = 0; i < 4; i++) {
            if (i == first) continue;

            E_Int fid = FRO[i];
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
        fid = FRO[second];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[8], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[4]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[8]);
        assert(local[1] == NODES[0]);
        assert(local[2] == NODES[15]);
        assert(local[3] == NODES[22]);

        // Setup third face
        fid = FRO[third];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[22], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[4]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[22]);
        assert(local[1] == NODES[15]);
        assert(local[2] == NODES[4]);
        NODES[21] = local[3];

        // Setup fourth face
        fid = FRO[fourth];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[19], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[4]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[19]);
        assert(local[1] == NODES[22]);
        assert(local[2] == NODES[21]);
        assert(local[3] == NODES[5]);

        E_Int tmp[4] = {FRO[first], FRO[second], FRO[third], FRO[fourth]};
        for (E_Int i = 0; i < 4; i++) FRO[i] = tmp[i];

    } else {
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
        assert(local[3] == NODES[15]);
        assert(local[4] == NODES[4]);
        assert(local[6] == NODES[5]);
        assert(local[7] == NODES[19]);
        NODES[21] = local[5];
    }

    // Reorder BCK sides

    if (crange[5] == 4) {
        // First face must share NODES[2]

        E_Int first = -1;

        for (E_Int i = 0; i < 4 && first == -1; i++) {
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
        assert(local[3] == NODES[17]);
        NODES[24] = local[2];

        // Get second, third and fourth sides
        E_Int second, third, fourth;
        second = third = fourth = -1;

        for (E_Int i = 0; i < 4; i++) {
            if (i == first) continue;
            
            E_Int fid = BCK[i];
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
        fid = BCK[second];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[10], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[5]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[10]);
        assert(local[1] == NODES[3]);
        assert(local[2] == NODES[13]);
        assert(local[3] == NODES[24]);

        // Setup third face
        fid = BCK[third];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[24], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[5]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[24]);
        assert(local[1] == NODES[13]);
        assert(local[2] == NODES[7]);
        NODES[23] = local[3];

        // Setup fourth face
        fid = BCK[fourth];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[17], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[5]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[17]);
        assert(local[1] == NODES[24]);
        assert(local[2] == NODES[23]);
        assert(local[3] == NODES[6]);

        E_Int tmp[4] = {BCK[first], BCK[second], BCK[third], BCK[fourth]};
        for (E_Int i = 0; i < 4; i++) BCK[i] = tmp[i];

    } else {
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
        assert(local[3] == NODES[13]);
        assert(local[4] == NODES[7]);
        assert(local[6] == NODES[6]);
        assert(local[7] == NODES[17]);
        NODES[23] = local[5];
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
        assert(local[1] == NODES[21]);
        assert(local[3] == NODES[14]);
        NODES[25] = local[2];

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
        i0 = Get_pos(NODES[21], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[1]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[21]);
        assert(local[1] == NODES[5]);
        assert(local[2] == NODES[18]);
        assert(local[3] == NODES[25]);

        // Setup third face
        fid = TOP[third];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[25], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[1]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[25]);
        assert(local[1] == NODES[18]);
        assert(local[2] == NODES[6]);
        assert(local[3] == NODES[23]);

        // Setup fourth face
        fid = TOP[fourth];
        pn = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
        i0 = Get_pos(NODES[14], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[1]);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[14]);
        assert(local[1] == NODES[25]);
        assert(local[2] == NODES[23]);
        assert(local[3] == NODES[7]);

        E_Int tmp[4] = {TOP[first], TOP[second], TOP[third], TOP[fourth]};
        for (E_Int i = 0; i < 4; i++) TOP[i] = tmp[i];

    } else {
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
        assert(local[1] == NODES[21]);
        assert(local[2] == NODES[5]);
        assert(local[3] == NODES[18]);
        assert(local[4] == NODES[6]);
        assert(local[5] == NODES[23]);
        assert(local[6] == NODES[7]);
        assert(local[7] == NODES[14]);
    }
}

void reconstruct_quad(Mesh *M, E_Int hexa, E_Int *fids, E_Int crange, E_Int normalIn,
    E_Int NODE, E_Int pn[4])
{
    E_Int fid = fids[0];

    E_Int local[4];
    E_Int NODES[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
    E_Int i0 = 0;

    // Setup first face
    fid = fids[0];
    E_Int *face = Mesh_get_face(M, fid);
    for (E_Int i = 0; i < 4; i++) local[i] = face[2*i];
    if (NODE != -1) {
        i0 = Get_pos(NODE, local, 4);
        Right_shift(local, i0, 4);
    }
    E_Int reorient = Mesh_get_reorient(M, fid, hexa, normalIn);
    if (reorient) std::swap(local[1], local[3]);

    if (crange == 1) {
        for (E_Int i = 0; i < 4; i++) pn[i] = local[i];
    } else if (crange == 4) {
        NODES[0] = local[0];
        NODES[4] = local[1];
        NODES[8] = local[2];
        NODES[7] = local[3];

        E_Int map[4];
        for (E_Int i = 0; i < 4; i++) map[i] = local[i];

        // Find face order

        E_Int second, third, fourth;
        second = third = fourth = 0;

        for (E_Int i = 1; i < 4; i++) {
            E_Int fid = fids[i];
            E_Int *pn = Mesh_get_face(M, fid);
            for (E_Int j = 0; j < 4; j++) local[j] = pn[2*j];
            E_Int common[4] = {0, 0, 0, 0};
            for (E_Int j = 0; j < 4; j++) {
                E_Int point = local[j];
                for (E_Int k = 0; k < 4; k++) {
                    if (map[k] == point) {
                        common[k] = 1;
                        break;
                    }
                }
            }
            
            assert(common[0] == 0);
            if (common[1] && common[2]) second = i;
            else if (common[2] && common[3]) fourth = i;
            else third = i;
        }

        assert(second == 1);
        assert(third == 2);
        assert(fourth == 3);

        // Second face
        fid = fids[second];
        face = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = face[2*i];
        i0 = Get_pos(NODES[4], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[4]);
        assert(local[3] == NODES[8]);
        NODES[1] = local[1];
        NODES[5] = local[2];

        // Third face
        fid = fids[third];
        face = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = face[2*i];
        i0 = Get_pos(NODES[8], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[8]);
        assert(local[1] == NODES[5]);
        NODES[2] = local[2];
        NODES[6] = local[3];

        // Fourth face
        fid = fids[fourth];
        face = Mesh_get_face(M, fid);
        for (E_Int i = 0; i < 4; i++) local[i] = face[2*i];
        i0 = Get_pos(NODES[7], local, 4);
        Right_shift(local, i0, 4);
        reorient = Mesh_get_reorient(M, fid, hexa, normalIn);
        if (reorient) std::swap(local[1], local[3]);
        assert(local[0] == NODES[7]);
        assert(local[1] == NODES[8]);
        assert(local[2] == NODES[6]);
        NODES[3] = local[3];

        for (E_Int i = 0; i < 4; i++) pn[i] = NODES[i];
    } else {
        assert(0);
    }
}


E_Int check_canon_hexa(E_Int hexa, Mesh *M)
{
    E_Int NODES[8] = {-1, -1, -1, -1, -1, -1, -1, -1};

    E_Int *cell = Mesh_get_cell(M, hexa);
    E_Int *crange = Mesh_get_crange(M, hexa);

    E_Int local[4];

    E_Int *BOT = cell;
    E_Int *TOP = cell+4;
    E_Int *LFT = cell+8;
    E_Int *RGT = cell+12;
    E_Int *FRO = cell+16;
    E_Int *BCK = cell+20;

    NODES[0] = Mesh_get_face(M, BOT[0])[0];

    // BOT
    reconstruct_quad(M, hexa, BOT, crange[0], normalIn_H[0], NODES[0], local);
    for (E_Int i = 0; i < 4; i++) NODES[i] = local[i];

    // LFT
    reconstruct_quad(M, hexa, LFT, crange[2], normalIn_H[2], NODES[0], local);
    assert(local[0] == NODES[0]);
    assert(local[1] == NODES[3]);
    NODES[7] = local[2];
    NODES[4] = local[3];

    // RGT
    reconstruct_quad(M, hexa, RGT, crange[3], normalIn_H[3], NODES[1], local);
    assert(local[0] == NODES[1]);
    assert(local[1] == NODES[2]);
    NODES[6] = local[2];
    NODES[5] = local[3];

    // FRO
    reconstruct_quad(M, hexa, FRO, crange[4], normalIn_H[4], NODES[1], local);
    assert(local[0] == NODES[1]);
    assert(local[1] == NODES[0]);
    assert(local[2] == NODES[4]);
    assert(local[3] == NODES[5]);

    // BCK
    reconstruct_quad(M, hexa, BCK, crange[5], normalIn_H[5], NODES[2], local);
    assert(local[0] == NODES[2]);
    assert(local[1] == NODES[3]);
    assert(local[2] == NODES[7]);
    assert(local[3] == NODES[6]);

    // TOP
    reconstruct_quad(M, hexa, TOP, crange[1], normalIn_H[1], NODES[4], local);
    assert(local[0] == NODES[4]);
    assert(local[1] == NODES[5]);
    assert(local[2] == NODES[6]);
    assert(local[3] == NODES[7]);

    return 0;
}

void update_range_and_stride(Mesh *M, E_Int hexa, E_Int cpos, E_Int nchildren)
{
    E_Int *crange = Mesh_get_crange(M, hexa);
    for (E_Int i = 0; i < M->cstride[hexa]; i++) {
        crange[i] = 1;
    }

    for (E_Int i = 0; i < nchildren; i++) {
        E_Int child = cpos + i;

        M->cstride[child] = M->cstride[hexa];

        crange = Mesh_get_crange(M, child);
        for (E_Int j = 0; j < M->cstride[child]; j++) {
            crange[j] = 1;
        }
    }
}
