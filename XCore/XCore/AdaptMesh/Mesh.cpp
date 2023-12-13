#include "Proto.h"
#include <cassert>

const E_Int normalIn_T[4] = {1, 0, 1, 0};
const E_Int normalIn_P[5] = {1, 0, 1, 0, 1};
const E_Int normalIn_H[6] = {1, 0, 1, 0, 1, 0};

Edge::Edge()
{}

Edge::Edge(E_Int p0, E_Int p1) :
  p0_(std::min(p0, p1)), p1_(std::max(p0, p1))
{}

void Edge::set(E_Int p0, E_Int p1)
{
  p0_ = std::min(p0, p1);
  p1_ = std::max(p0, p1);
}

bool Edge::operator<(const Edge &e) const
{
  return (p0_ < e.p0_) || (p0_ == e.p0_ && p1_ < e.p1_);
}

AMesh::AMesh() :
  ncells(-1), nfaces(-1), npoints(-1),
  x(NULL), y(NULL), z(NULL),
  owner(NULL), neigh(NULL),
  nface(NULL), indPH(NULL), ngon(NULL), indPG(NULL),
  ecenter(),
  patches(NULL), npatches(-1),
  pid(-1), npc(-1), nreq(-1), req(NULL),
  cellTree(NULL), faceTree(NULL),
  fc(NULL), fa(NULL), cx(NULL), cy(NULL), cz(NULL),
  xn(NULL), yn(NULL), zn(NULL), fn(NULL),
  lsqG(NULL), lsqGG(NULL), lsqH(NULL), lsqHH(NULL),
  ref_data(NULL), ref_Tr(-1.0), unref_Tr(-1.0),
  nref_hexa(-1), nref_tetra(-1), nref_penta(-1), nref_pyra(-1),
  nunref_hexa(-1), nunref_tetra(-1), nunref_penta(-1), nunref_pyra(-1)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &npc);
  nreq = 0;
  req = (MPI_Request *)XMALLOC(2*npc * sizeof(MPI_Request));
}

E_Int *get_face(E_Int i, E_Int &np, E_Int *ngon, E_Int *indPG)
{
  np = indPG[i+1] - indPG[i];
  return &ngon[indPG[i]];
}

E_Int *get_cell(E_Int i, E_Int &nf, E_Int *nface, E_Int *indPH)
{
  nf = indPH[i+1] - indPH[i];
  return &nface[indPH[i]];
}

E_Int get_stride(E_Int i, E_Int *indir)
{
  return indir[i+1] - indir[i];
}

E_Int *get_facets(E_Int i, E_Int *cn, E_Int *indir)
{
  return &cn[indir[i]];
}

void compute_face_center_and_area(E_Int id, E_Int stride,
  E_Int *pn, E_Float *x, E_Float *y, E_Float *z, E_Float *fc, E_Float *fa)
{
  // Init
  fa[0] = fa[1] = fa[2] = 0.0;
  fc[0] = fc[1] = fc[2] = 0.0;

  // Approximate face center
  E_Float fcenter[3] = {0,0,0};

  for (E_Int i = 0; i < stride; i++) {
    E_Int point = pn[i];
    fcenter[0] += x[point];
    fcenter[1] += y[point];
    fcenter[2] += z[point];
  }

  for (E_Int i = 0; i < 3; i++) fcenter[i] /= stride;

  // Sum of triangle area vectors
  E_Float sumN[3] = {0,0,0};
  // Sum of triangle areas
  E_Float sumA = 0;
  // Sum of area-weighted triangle centers
  E_Float sumAc[3] = {0,0,0};

  // Compute area vector and center of stride-2 triangles
  // formed by p0p1p2, p0p2p3, ... p0p(stride-2)p(stride-1)
  E_Int p0 = pn[0];
  for (E_Int i = 1; i < stride-1; i++) {
    E_Int p1 = pn[i];
    E_Int p2 = pn[i+1];

    // Triangle center
    E_Float tc[3];
    tc[0] = x[p0] + x[p1] + x[p2];
    tc[1] = y[p0] + y[p1] + y[p2];
    tc[2] = z[p0] + z[p1] + z[p2];

    // Area vector
    E_Float n[3];
    E_Float v10[3] = {x[p1]-x[p0], y[p1]-y[p0], z[p1]-z[p0]};
    E_Float v20[3] = {x[p2]-x[p0], y[p2]-y[p0], z[p2]-z[p0]};
    K_MATH::cross(v10, v20, n);

    // Area
    E_Float a = K_MATH::norm(n, 3);

    for (E_Int j = 0; j < 3; j++) {
      sumN[j] += n[j];
      sumAc[j] += a*tc[j];
    }
    sumA += a;
  }

  // Deal with zero-area faces
  if (sumA < K_MATH::SMALL) {
    fprintf(stderr, "compute_face_area_and_center(): "
      "Warning: Face: %d - Area: %f - Tol: %.2e\n", id, sumA, K_MATH::SMALL);
    for (E_Int i = 0; i < 3; i++) {
      fc[i] = fcenter[i];
      fa[i] = 0.0;
    }
  } else {
    for (E_Int i = 0; i < 3; i++) {
      fc[i] = sumAc[i]/(3.0*sumA);
      fa[i] = 0.5*sumN[i];
    }
  }
}

void make_cell_centers(AMesh *M)
{
  // Make face centers and areas
  M->fc = (E_Float *)XRESIZE(M->fc, M->nfaces*3 * sizeof(E_Float));
  M->fa = (E_Float *)XRESIZE(M->fa, M->nfaces*3 * sizeof(E_Float));
  
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int np = -1;
    E_Int *pn = get_face(i, np, M->ngon, M->indPG);
    compute_face_center_and_area(i, np, pn, M->x, M->y, M->z,
      &M->fc[3*i], &M->fa[3*i]);
  }

  // Estimate cell centers as average of face centers
  E_Float *cEst = (E_Float *)XCALLOC(M->ncells*3, sizeof(E_Float));
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Float *fc = &M->fc[3*i];
    E_Float *cE;

    E_Int own = M->owner[i];
    assert(own < M->ncells);
    cE = &cEst[3*own];

    for (E_Int j = 0; j < 3; j++) cE[j] += fc[j];

    E_Int nei = M->neigh[i];
    assert(nei < M->ncells);
    if (nei == -1) continue;

    cE = &cEst[3*nei];

    for (E_Int j = 0; j < 3; j++) cE[j] += fc[j];
  }

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int nf = -1;
    get_cell(i, nf, M->nface, M->indPH);
    E_Float *cE = &cEst[3*i];
    for (E_Int j = 0; j < 3; j++)
      cE[j] /= nf;
  }

  E_Float *vols = (E_Float *)XMALLOC(M->ncells * sizeof(E_Float));
  M->cx = (E_Float *)XRESIZE(M->cx, M->ncells * sizeof(E_Float));
  M->cy = (E_Float *)XRESIZE(M->cy, M->ncells * sizeof(E_Float));
  M->cz = (E_Float *)XRESIZE(M->cz, M->ncells * sizeof(E_Float));
   
  memset(vols,  0, M->ncells*sizeof(E_Float));
  memset(M->cx, 0, M->ncells*sizeof(E_Float));
  memset(M->cy, 0, M->ncells*sizeof(E_Float));
  memset(M->cz, 0, M->ncells*sizeof(E_Float));

  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Float *fa = &M->fa[3*i];
    E_Float *fc = &M->fc[3*i];
    E_Float pyr3vol, pc[3], *cE, d[3];

    E_Int own = M->owner[i];
    cE = &cEst[3*own];
    for (E_Int j = 0; j < 3; j++) d[j] = fc[j]-cE[j];
    pyr3vol = K_MATH::dot(fa, d, 3);
    for (E_Int j = 0; j < 3; j++) pc[j]= 0.75*fc[j] + 0.25*cE[j];
    M->cx[own] += pyr3vol*pc[0];
    M->cy[own] += pyr3vol*pc[1];
    M->cz[own] += pyr3vol*pc[2];
    vols[own] += pyr3vol;

    E_Int nei = M->neigh[i];
    if (nei == -1) continue;

    cE = &cEst[3*nei];
    for (E_Int j = 0; j < 3; j++) d[j] = cE[j]-fc[j];
    pyr3vol = K_MATH::dot(fa, d, 3);
    for (E_Int j = 0; j < 3; j++) pc[j]= 0.75*fc[j] + 0.25*cE[j];
    M->cx[nei] += pyr3vol*pc[0];
    M->cy[nei] += pyr3vol*pc[1];
    M->cz[nei] += pyr3vol*pc[2];
    vols[nei] += pyr3vol;
  }

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Float coeff = 1.0/vols[i];
    M->cx[i] *= coeff;
    M->cy[i] *= coeff;
    M->cz[i] *= coeff;
  }

  XFREE(cEst);
  XFREE(vols);
}

E_Int set_faces_type(AMesh *M)
{
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int np = get_stride(i, M->indPG);
    switch (np) {
      case 3:
        M->faceTree[i]->type = TRI;
        break;
      case 4:
        M->faceTree[i]->type = QUAD;
        break;
      default:
        return 1;
    }
  }

  return 0;
}

E_Int set_cells_type(AMesh *M)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int nf = get_stride(i, M->indPH);
    switch (nf) {
      case 4:
        M->cellTree[i]->type = TETRA;
        break;
      case 5: {
        // Prism: 2 TRI + 3 QUAD
        // Pyra: 4 TRI + 1 QUAD
        E_Int ntri, nquad;
        ntri = nquad = 0;
        E_Int *pf = &M->nface[M->indPH[i]];
        for (E_Int j = 0; j < nf; j++)
          M->faceTree[pf[j]]->type == TRI ? ntri++ : nquad++;
        E_Int is_penta = ntri == 2 && nquad == 3;
        E_Int is_pyra = ntri == 4 && nquad == 1;
        if (is_penta) M->cellTree[i]->type = PENTA;
        else M->cellTree[i]->type = PYRA;
        break;
      }
      case 6:
        M->cellTree[i]->type = HEXA;
        break;
      default:
        assert(0);
        return 1;
    }
  }

  return 0;
}

E_Int get_reorient(E_Int face, E_Int cell, E_Int normalIn, AMesh *M)
{
  assert(M->owner[face] == cell || M->neigh[face] == cell);
  if (M->neigh[face] == cell && normalIn == 1) return 0;
  else if (M->owner[face] == cell && normalIn == 0) return 0;
  else return 1;
}

static
void reorder_tetra(E_Int i, E_Int nf, E_Int *pf, AMesh *M)
{
  E_Int common[3], map[3];
  E_Int bot = pf[0];
  E_Int *pn = get_facets(bot, M->ngon, M->indPG);

  for (E_Int j = 0; j < 3; j++) map[j] = pn[j];
  E_Int reorient = get_reorient(bot, i, normalIn_T[0], M);
  if (reorient) std::swap(map[1], map[2]);

  E_Int lft, rgt, fro;
  lft = rgt = fro = -1;

  for (E_Int j = 1; j < 4; j++) {
    E_Int face = pf[j];

    for (E_Int k = 0; k < 3; k++) common[k] = 0;

    E_Int *pnn = get_facets(face, M->ngon, M->indPG);

    for (E_Int k = 0; k < 3; k++) {
      E_Int point = pnn[k];

      for (E_Int l = 0; l < 3; l++) {
        if (map[l] == point) {
          common[l] = 1;
          break;
        }
      }
    }

    if      (common[0] && common[2]) lft = face;
    else if (common[1] && common[2]) rgt = face;
    else if (common[1] && common[0]) fro = face;
    else assert(0);
  }
  assert(lft != -1 && rgt != -1 && fro != -1);

  pf[1] = lft;
  pf[2] = rgt;
  pf[3] = fro;
}

static
void reorder_penta(E_Int i, E_Int nf, E_Int *pf, AMesh *M)
{
  E_Int common[4], map[3];
  E_Int bot = pf[0];
  E_Int *pn = get_facets(bot, M->ngon, M->indPG);
  
  for (E_Int i = 0; i < 3; i++) map[i] = pn[i];
  E_Int reorient = get_reorient(bot, i, normalIn_P[0], M);
  if (reorient) std::swap(map[1], map[2]);

  E_Int top, lft, rgt, bck;
  top = lft = rgt = bck = -1;

  for (E_Int j = 1; j < 5; j++) {
    E_Int face = pf[j];

    E_Int stride = get_stride(face, M->indPG);

    for (E_Int k = 0; k < stride; k++) common[k] = 0;

    E_Int *pn = get_facets(face, M->ngon, M->indPG);

    for (E_Int k = 0; k < stride; k++) {
      E_Int point = pn[k];

      for (E_Int l = 0; l < 3; l++) {
        if (map[l] == point) {
          common[l] = 1;
          break;
        }
      }
    }

    if      (common[0] && common[2]) lft = face;
    else if (common[0] && common[1]) rgt = face;
    else if (common[1] && common[2]) bck = face;
    else                             top = face;
  }
  assert(top != -1);
  assert(lft != -1);
  assert(rgt != -1);
  assert(bck != -1);

  pf[1] = top;
  pf[2] = lft;
  pf[3] = rgt;
  pf[4] = bck;
}

static
void reorder_pyra(E_Int i, E_Int nf, E_Int *pf, AMesh *M)
{}

static
void reorder_hexa(E_Int i, E_Int nf, E_Int *pf, AMesh *M)
{
  E_Int common[4], map[4];
  E_Int bot = pf[0];
  E_Int *pn = get_facets(bot, M->ngon, M->indPG);
  
  for (E_Int i = 0; i < 4; i++) map[i] = pn[i];
  E_Int reorient = get_reorient(bot, i, normalIn_H[0], M);
  if (reorient) std::swap(map[1], map[3]);

  E_Int top, lft, rgt, fro, bck;
  top = lft = rgt = fro = bck = -1;

  for (E_Int j = 1; j < 6; j++) {
    E_Int face = pf[j];

    for (E_Int k = 0; k < 4; k++) common[k] = 0;

    E_Int *pn = get_facets(face, M->ngon, M->indPG);

    for (E_Int k = 0; k < 4; k++) {
      E_Int point = pn[k];

      for (E_Int l = 0; l < 4; l++) {
        if (map[l] == point) {
          common[l] = 1;
          break;
        }
      }
    }

    if      (common[0] && common[3]) lft = face;
    else if (common[0] && common[1]) fro = face;
    else if (common[1] && common[2]) rgt = face;
    else if (common[3] && common[2]) bck = face;
    else                             top = face;
  }
  assert(top != -1);
  assert(lft != -1);
  assert(rgt != -1);
  assert(fro != -1);
  assert(bck != -1);

  pf[1] = top;
  pf[2] = lft;
  pf[3] = rgt;
  pf[4] = fro;
  pf[5] = bck;
}

void reorder_cells(AMesh *M)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int nf = -1;
    E_Int *pf = get_cell(i, nf, M->nface, M->indPH);
    E_Int ctype = M->cellTree[i]->type;

    switch (ctype) {
      case TETRA:
        reorder_tetra(i, nf, pf, M);
        break;
      case PENTA:
        reorder_penta(i, nf, pf, M);
        break;
      case PYRA:
        reorder_pyra(i, nf, pf, M);
        break;
      case HEXA:
        reorder_hexa(i, nf, pf, M);
        break;
      default:
        assert(0);
    }
  }
}

void check_canon_cells(AMesh *M)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int type = M->cellTree[i]->type;
    switch (type) {
      case HEXA:
        check_canon_hexa(i, M);
        break;
      case TETRA:
        check_canon_tetra(i, M);
        break;
      case PENTA:
        check_canon_penta(i, M);
        break;
      default:
        assert(0);
        break;
    }
  }
}

void Order_tri(E_Int *local, E_Int *pn, E_Int reorient, E_Int i0)
{
  for (E_Int i = 0; i < 3; i++) local[i] = pn[i];
  Right_shift(local, i0, 3);
  if (reorient) std::swap(local[1], local[2]);
}

void Order_quad(E_Int *local, E_Int *pn, E_Int reorient, E_Int i0)
{
  for (E_Int i = 0; i < 4; i++) local[i] = pn[i];
  Right_shift(local, i0, 4);
  if (reorient) std::swap(local[1], local[3]);
}

void check_canon_tetra(E_Int cell, AMesh *M)
{
  E_Int NODES[4] = {-1, -1, -1, -1};

  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, i0, reorient, *pn, local[3];

  // BOT
  face = pf[0];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn_T[0], M);
  Order_tri(local, pn, reorient, i0);
  for (E_Int i = 0; i < 3; i++) NODES[i] = local[i];

  // LFT
  face = pf[1];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 3);
  reorient = get_reorient(face, cell, normalIn_T[1], M);
  Order_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[2] == NODES[2]);
  NODES[3] = local[1];

  // RGT
  face = pf[2];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[1], pn, 3);
  reorient = get_reorient(face, cell, normalIn_T[2], M);
  Order_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[3]);
  assert(local[2] == NODES[2]);

  // FRO
  face = pf[3];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 3);
  reorient = get_reorient(face, cell, normalIn_T[3], M);
  Order_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[1]);
  assert(local[2] == NODES[3]);
}

void check_canon_hexa(E_Int cell, AMesh *M)
{
  E_Int NODES[8] = {-1, -1, -1, -1, -1, -1, -1, -1};

  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, i0, reorient, *pn, local[4];

  // BOT
  face = pf[0];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn_H[0], M);
  Order_quad(local, pn, reorient, i0);
  for (E_Int i = 0; i < 4; i++) NODES[i] = local[i];

  // LFT
  face = pf[2];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 4);
  reorient = get_reorient(face, cell, normalIn_H[2], M);
  Order_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[3]);
  NODES[7] = local[2];
  NODES[4] = local[3];

  // RGT
  face = pf[3];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[1], pn, 4);
  reorient = get_reorient(face, cell, normalIn_H[3], M);
  Order_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[2]);
  NODES[6] = local[2];
  NODES[5] = local[3];

  // FRO
  face = pf[4];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[1], pn, 4);
  reorient = get_reorient(face, cell, normalIn_H[4], M);
  Order_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[0]);
  assert(local[2] == NODES[4]);
  assert(local[3] == NODES[5]);

  // BCK
  face = pf[5];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[2], pn, 4);
  reorient = get_reorient(face, cell, normalIn_H[5], M);
  Order_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[3]);
  assert(local[2] == NODES[7]);
  assert(local[3] == NODES[6]);

  // TOP
  face = pf[1];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[4], pn, 4);
  reorient = get_reorient(face, cell, normalIn_H[1], M);
  Order_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[4]);
  assert(local[1] == NODES[5]);
  assert(local[2] == NODES[6]);
  assert(local[3] == NODES[7]);
}

void check_canon_penta(E_Int cell, AMesh *M)
{
  E_Int NODES[6] = {-1, -1, -1, -1, -1, -1};

  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, i0, reorient, *pn, local[4];

  // BOT
  face = pf[0];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn_P[0], M);
  Order_tri(local, pn, reorient, i0);
  for (E_Int i = 0; i < 3; i++) NODES[i] = local[i];

  // LFT (in)
  face = pf[2];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 4);
  assert(i0 != -1);
  reorient = get_reorient(face, cell, normalIn_P[2], M);
  Order_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[2]);
  NODES[5] = local[2];
  NODES[3] = local[3];

  // RGT (out)
  face = pf[3];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 4);
  reorient = get_reorient(face, cell, normalIn_P[3], M);
  Order_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[1]);
  assert(local[3] == NODES[3]);
  NODES[4] = local[2];

  // BCK (in)
  face = pf[4];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[2], pn, 4);
  reorient = get_reorient(face, cell, normalIn_P[4], M);
  Order_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[1]);
  assert(local[2] == NODES[4]);
  assert(local[3] == NODES[5]);

  // TOP
  face = pf[1];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[3], pn, 3);
  reorient = get_reorient(face, cell, normalIn_P[1], M);
  Order_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[3]);
  assert(local[1] == NODES[4]);
  assert(local[2] == NODES[5]);
}


/*
static
void compute_tetra_vol(E_Int cell, AMesh *M)
{
  // Make nodes
  E_Int *pf = &M->nface[M->indPH[cell]];
  E_Int bot = pf[0];
  E_Int *pn = &M->ngon[M->indPG[bot]];
  E_Int n[4] = {pn[0], pn[1], pn[2], -1};

  E_Int lft = pf[1];
  pn = &M->ngon[M->indPG[lft]];

  for (E_Int i = 0; i < 3; i++) {
    E_Int pt = pn[i];

    E_Int found = 0;
    // pt must not be in n
    for (E_Int j = 0; j < 3; j++) {
      if (pt == n[j]) {
        found = 1;
        break;
      }
    }

    if (!found) {
      n[3] = pt;
      break;
    }
  }

  E_Float xa = M->x[n[0]]; 
  E_Float ya = M->y[n[0]]; 
  E_Float za = M->z[n[0]]; 
  E_Float xb = M->x[n[1]]; 
  E_Float yb = M->y[n[1]]; 
  E_Float zb = M->z[n[1]]; 
  E_Float xc = M->x[n[2]]; 
  E_Float yc = M->y[n[2]]; 
  E_Float zc = M->z[n[2]]; 
  E_Float xd = M->x[n[3]]; 
  E_Float yd = M->y[n[3]]; 
  E_Float zd = M->z[n[3]]; 
  E_Float ad[3] = {xa-xd, ya-yd, za-zd};
  E_Float bd[3] = {xb-xd, yb-yd, zb-zd};
  E_Float cd[3] = {xc-xd, yc-yd, zc-zd};
  E_Float tmp[3];
  K_MATH::cross(bd, cd, tmp);
  E_Float V = fabs(K_MATH::dot(ad, tmp, 3))/6.0;
  printf("vol %d: %.4e\n", cell, V);
}

static
void compute_tetra_vols(AMesh *M)
{
  for (E_Int i = 0; i < M->ncells; i++)
    compute_tetra_vol(i, M);
}
*/
