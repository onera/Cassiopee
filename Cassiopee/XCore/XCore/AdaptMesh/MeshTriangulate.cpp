#include "Mesh.h"
#include "Hexa.h"

void Mesh_face_to_prism(Mesh *M, E_Int fid)
{
    E_Int hexa = M->owner[fid];
    assert(M->ctype[hexa] == HEXA);

    // Make fid the bottom face
    E_Int *cell = Mesh_get_cell(M, hexa);
    E_Int pos = Get_pos(fid, cell, 24);
    assert(pos != -1);
    assert(pos % 4 == 0);

    Right_shift(cell, pos, 24);
    assert(cell[0] == fid);
    E_Int *crange = Mesh_get_crange(M, hexa);
    Right_shift(crange, pos/4, 6);

    // Get the top and left face
    // Top face shares no point, left face shares p0 and p3

    E_Int *quad = Mesh_get_face(M, fid);
    
    E_Int map[4];
    for (E_Int i = 0; i < 4; i++)
        map[i] = quad[2*i];
    E_Int reorient = Mesh_get_reorient(M, fid, hexa, normalIn_H[0]);
    E_Int BOT[2] = {fid, M->nf};
    if (reorient) {
        std::swap(map[1], map[3]);
        std::swap(BOT[0], BOT[1]);
    }

    E_Int top=-1, lft=-1, rgt=-1, fro=-1, bck=-1;

    for (E_Int i = 1; i < 6; i++) {
        E_Int f = cell[4*i];
        E_Int *face = Mesh_get_face(M, f);
        E_Int common[4] = {0, 0, 0, 0};

        for (E_Int j = 0; j < 8; j += 2) {
            E_Int p = face[j];
            for (E_Int k = 0; k < 4; k++) {
                if (p == map[k]) {
                    common[k] = 1;
                    break;
                }
            }
        }

        if (common[0] && common[3]) lft = f;
        else if (common[1] && common[2]) rgt = f;
        else if (common[1] && common[0]) fro = f;
        else if (common[2] && common[3]) bck = f;
        else top = f;

    }

    assert(top != -1);
    assert(lft != -1);
    assert(rgt != -1);
    assert(fro != -1);
    assert(bck != -1);
    
    // Setup the canon config

    E_Int NODES[8] = {};
    for (E_Int i = 0; i < 4; i++)
        NODES[i] = map[i];
    
    E_Int local[4];

    // Align left with bottom

    E_Int *pn = Mesh_get_face(M, lft);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    pos = Get_pos(NODES[0], local, 4);
    assert(pos != -1);
    assert(pos % 2 == 0);
    Right_shift(local, pos, 4); 
    assert(local[0] == NODES[0]);
    reorient = Mesh_get_reorient(M, lft, hexa, normalIn_H[2]);
    if (reorient) std::swap(local[1], local[3]);
    assert(local[1] == NODES[3]);
    NODES[4] = local[3];
    NODES[7] = local[2];

    // Align top with left
    pn = Mesh_get_face(M, top);
    for (E_Int i = 0; i < 4; i++) local[i] = pn[2*i];
    pos = Get_pos(NODES[4], local, 4);
    assert(pos != -1);
    assert(pos % 2 == 0);
    Right_shift(local, pos, 4); 
    assert(local[0] == NODES[4]);
    E_Int TOP[2] = {top, M->nf+1};
    reorient = Mesh_get_reorient(M, top, hexa, normalIn_H[1]);
    if (reorient) {
        std::swap(local[1], local[3]);
        if (pos == 0 || pos == 3) std::swap(TOP[0], TOP[1]);
    } else {
        if (pos == 2 || pos == 3) std::swap(TOP[0], TOP[1]);
    }
    assert(local[3] == NODES[7]);
    NODES[5] = local[1];
    NODES[6] = local[2];

    // Triangulate bottom
    M->fstride[fid] = 3;
    quad = Mesh_get_face(M, fid);
    E_Int *qrange = Mesh_get_frange(M, fid);
    qrange[2] = 1;

    M->fstride[M->nf] = 3;
    E_Int *tri = Mesh_get_face(M, M->nf);
    tri[0] = quad[6]; tri[2] = quad[0]; tri[4] = quad[4];
    E_Int *trange = Mesh_get_frange(M, M->nf);
    trange[0] = qrange[3]; trange[1] = 1; trange[2] = qrange[2];
    
    // Triangulate top
    M->fstride[top] = 3;
    quad = Mesh_get_face(M, top);
    qrange = Mesh_get_frange(M, top);
    qrange[2] = 1;

    M->fstride[M->nf+1] = 3;
    tri = Mesh_get_face(M, M->nf+1);
    tri[0] = quad[6]; tri[2] = quad[0]; tri[4] = quad[4];
    trange = Mesh_get_frange(M, M->nf+1);
    trange[0] = qrange[3]; trange[1] = 1; trange[2] = qrange[2];

    // Update face info
    M->fref[fid] = M->fref[top] = FACE_REFINED;
    M->fref[M->nf] = M->fref[M->nf+1] = FACE_NEW;
    M->flevel[M->nf] = M->flevel[fid];
    M->flevel[M->nf+1] = M->flevel[top];
    M->ftag[M->nf] = M->ftag[fid];
    M->ftag[M->nf+1] = M->ftag[top];
    M->fparent[M->nf] = M->nf;
    M->fparent[M->nf+1] = M->nf+1;
    M->owner[M->nf] = M->owner[fid];
    M->neigh[M->nf] = M->neigh[fid];
    M->owner[M->nf+1] = M->owner[top];
    M->neigh[M->nf+1] = M->neigh[top];
    M->ftype[fid] = M->ftype[M->nf] = TRI;
    M->ftype[top] = M->ftype[M->nf+1] = TRI;

    
    // Conformize parent cells
    E_Int ret, own, nei;

    M->fchildren[fid] = {fid, M->nf};
    M->fchildren[top] = {top, M->nf+1};

    // fid
    ret = Mesh_conformize_cell_face(M, hexa, fid, M->nf, 2);
    assert(ret == 0);
    nei = M->neigh[fid];
    assert(nei == -1);
    /*
    if (nei != -1) {
        ret = Mesh_conformize_cell_face(M, nei, hexa, M->nf, 2);
        assert(ret == 0);
    }
    */

    // top
    own = M->owner[top];
    ret = Mesh_conformize_cell_face(M, own, top, M->nf+1, 2);
    assert(ret == 0);
    nei = M->neigh[top];
    if (nei != -1) {
        ret = Mesh_conformize_cell_face(M, nei, top, M->nf+1, 2);
        assert(ret == 0);
    }

    // Set internal face in ngon
    
    E_Int *iquad = Mesh_get_face(M, M->nf+2);
    iquad[0] = NODES[0]; iquad[2] = NODES[2];
    iquad[4] = NODES[6]; iquad[6] = NODES[4];

    // Set internal face stride, range and state
    E_Int *irange = Mesh_get_frange(M, M->nf+2);
    for (E_Int i = 0; i < 4; i++) irange[i] = 1;
    M->fstride[M->nf+2] = 4;
    M->fref[M->nf+2] = FACE_NEW;

    // Make the first prism
    E_Int *prism = Mesh_get_cell(M, hexa);
    prism[0]  = BOT[0];  prism[4]  = TOP[0];
    prism[8]  = M->nf+2; prism[12] = rgt;
    prism[16] = fro;
    
    // Make the second prism
    prism = Mesh_get_cell(M, M->nc);
    prism[0]  = BOT[1]; prism[4]  = TOP[1];
    prism[8]  = bck;    prism[12] = M->nf+2;
    prism[16] = lft;

    // Update ranges and strides
    
    M->cstride[hexa] = 5;
    M->cstride[M->nc] = 5;
    
    E_Int *hrange = Mesh_get_crange(M, hexa);
    
    E_Int *prange = Mesh_get_crange(M, M->nc);
    prange[0] = 1; prange[1] = 1; prange[2] = hrange[5];
    prange[3] = 1; prange[4] = hrange[2];

    hrange[0] = 1; hrange[1] = 1; hrange[2] = 1;

    // Update adaptation info

    M->clevel[M->nc] = M->clevel[hexa];
    M->ctype[M->nc] = M->ctype[hexa] = PENTA;

    M->cchildren[hexa] = {hexa, M->nc};
    
    // Set shell faces owns and neis

    update_shell_pe(hexa, M);

    // Set PE pf internal face

    M->owner[M->nf+2] = M->nc;
    M->neigh[M->nf+2] = hexa; 
    
    M->nc++;
    M->nf += 3;
}

void Mesh_generate_prisms(Mesh *M, E_Int *faces, E_Int nf)
{
    E_Int new_nf = M->nf + 3*nf;
    E_Int new_nc = M->nc + nf;

    Mesh_resize_face_data(M, new_nf);
    Mesh_resize_cell_data(M, new_nc);

    for (E_Int i = 0; i < nf; i++) {
        E_Int fid = faces[i];
        Mesh_face_to_prism(M, fid);
    }

    Mesh_update_bpatches(M);

    // TODO(Imad): parallel consequences
}

void Mesh_triangulate_face(Mesh *M, E_Int fid)
{ 
    // Get quad data
    E_Int *quad = Mesh_get_face(M, fid);
    E_Int *qrange = Mesh_get_frange(M, fid);

    // Create the second triangle
    E_Int *tri = Mesh_get_face(M, M->nf);
    E_Int *trange = Mesh_get_frange(M, M->nf);

    M->fstride[M->nf] = 3;
    tri[0] = quad[4]; tri[1] = quad[5];
    tri[2] = quad[6]; tri[3] = quad[7];
    tri[4] = quad[0]; assert(tri[5] == -1); 
    trange[0] = qrange[2];
    trange[1] = qrange[4];
    trange[2] = 1;

    // Update first quad -> triangle
    M->fstride[fid] = 3;
    qrange[2] = 1;

    // Conformize parent elements
    E_Int own = M->owner[fid];
    E_Int ret = Mesh_conformize_cell_face(M, own, fid, M->nf, 2);
    assert(ret == 0);

    E_Int nei = M->neigh[fid];
    assert(nei == -1);
    if (nei != -1) {
        ret = Mesh_conformize_cell_face(M, nei, fid, M->nf, 2);
        assert(ret == 0);
    }

    M->owner[M->nf] = own;
    M->neigh[M->nf] = nei; 

    M->flevel[M->nf] = M->flevel[fid];
    M->ftype[M->nf] = M->ftype[fid] = TRI;
    M->fparent[M->nf] = M->nf;
    
    M->fref[fid] = FACE_REFINED;
    M->fref[M->nf] = FACE_NEW;
    M->fchildren[fid] = {fid, M->nf};
    
    M->ftag[M->nf] = M->ftag[fid];

    M->nf++;
}

void Mesh_triangulate_faces(Mesh *M, E_Int *faces, E_Int nf)
{
    E_Int new_nf = M->nf + nf;

    Mesh_resize_face_data(M, new_nf);

    for (E_Int i = 0; i < nf; i++) {
        E_Int fid = faces[i];
        assert(M->ftag[fid] == 1);
        assert(M->fstride[fid] == 4);
        Mesh_triangulate_face(M, fid);
    }

    Mesh_update_bpatches(M);

    // TODO(Imad): parallel consequences
}
