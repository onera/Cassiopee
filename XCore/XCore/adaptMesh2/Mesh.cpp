#include "Proto.h"
#include <cassert>

E_Int master_face(E_Int face, AMesh *M)
{
  return M->faceTree->children(face) ? face : M->faceTree->parent(face);
}

//E_Int master_cell(E_Int cell, AMesh *M)
//{
//  return cell == -1 ? cell : 
//    (M->cellTree->children(cell) ? cell : M->cellTree->parent(cell)); 
//}

E_Int master_cell(E_Int cell, AMesh *M)
{
  if (cell == -1) return -1;
  return M->cellTree->children(cell) ? cell : M->cellTree->parent(cell);
}

const char *type_to_string(E_Int type)
{
  switch (type) {
    case HEXA: return "HEXA";
    case TETRA: return "TETRA";
    case PENTA: return "PENTA";
    case PYRA: return "PYRA";
    default:
      assert(0);
      return NULL;
  }
}

void ngon_print(AMesh *M)
{
  puts("");
  for (E_Int i = 0; i < M->nfaces; i++) {
    for (E_Int j = M->indPG[i]; j < M->indPG[i+1]; j++)
      printf("%d ", M->ngon[j]);
    puts("");
  }
  puts("");
}

void nface_print(AMesh *M)
{
  puts("");
  for (E_Int i = 0; i < M->ncells; i++) {
    for (E_Int j = M->indPH[i]; j < M->indPH[i+1]; j++)
      printf("%d ", M->nface[j]);
    puts("");
  }
  puts("");
}

const E_Int normalIn_T[4] = {1, 0, 1, 0};
const E_Int normalIn_Pe[5] = {1, 0, 1, 0, 1};
const E_Int normalIn_Py[5] = {1, 1, 0, 1, 0};
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
  ncells(-1), nfaces(-1), npoints(-1), nif(-1), nbf(-1),
  x(NULL), y(NULL), z(NULL),
  nface(NULL), indPH(NULL), ngon(NULL), indPG(NULL),
  owner(NULL), neigh(NULL),
  nbc(-1), ptlists(NULL), bcsizes(NULL), bcnames(NULL),
  patches(NULL), npatches(-1),
  pid(-1), npc(-1), nreq(-1), req(NULL),
  ecenter(NULL),
  cellTree(NULL), faceTree(NULL),
  ref_data(NULL), Tr(-1.0), Tu(-1.0),
  closed_indPG(NULL), closed_ngon(NULL),
  fc(NULL), cx(NULL), cy(NULL), cz(NULL)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  MPI_Comm_size(MPI_COMM_WORLD, &npc);
  nreq = 0;
  req = (MPI_Request *)XMALLOC(2*npc * sizeof(MPI_Request));
}

E_Int get_neighbour(E_Int cell, E_Int face, AMesh *M)
{
  assert(cell == M->owner[face] || cell == M->neigh[face]);
  if (cell == M->owner[face]) return M->neigh[face];
  return M->owner[face];
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

E_Int set_faces_type(AMesh *M)
{
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int np = get_stride(i, M->indPG);
    switch (np) {
      case 3:
        M->faceTree->type_[i] = TRI;
        break;
      case 4:
        M->faceTree->type_[i] = QUAD;
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
        M->cellTree->type_[i] = TETRA;
        break;
      case 5: {
        // Prism: 2 TRI + 3 QUAD
        // Pyra: 4 TRI + 1 QUAD
        E_Int ntri, nquad;
        ntri = nquad = 0;
        E_Int *pf = &M->nface[M->indPH[i]];
        for (E_Int j = 0; j < nf; j++)
          M->faceTree->type_[pf[j]] == TRI ? ntri++ : nquad++;
        E_Int is_penta = ntri == 2 && nquad == 3;
        E_Int is_pyra = ntri == 4 && nquad == 1;
        if (is_penta) M->cellTree->type_[i] = PENTA;
        else if (is_pyra) M->cellTree->type_[i] = PYRA;
        else assert(0);
        break;
      }
      case 6:
        M->cellTree->type_[i] = HEXA;
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
  // First tri is bottom
  E_Int bot = -1;
  for (E_Int i = 0; i < nf; i++) {
    if (get_stride(pf[i], M->indPG) == 3) {
      bot = pf[i];
      Right_shift(pf, i, 5);
      break;
    }
  }
  assert(bot != -1);

  E_Int common[3], map[3];
  E_Int *pn = get_facets(bot, M->ngon, M->indPG);
  
  for (E_Int i = 0; i < 3; i++) map[i] = pn[i];
  E_Int reorient = get_reorient(bot, i, normalIn_Pe[0], M);
  if (reorient) std::swap(map[1], map[2]);

  E_Int top, lft, rgt, bck;
  top = lft = rgt = bck = -1;

  for (E_Int j = 1; j < 5; j++) {
    for (E_Int k = 0; k < 3; k++) common[k] = 0;

    E_Int face = pf[j];
    E_Int stride = get_stride(face, M->indPG);
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
{
  E_Int common[4], map[4];
  E_Int bot = pf[0];
  E_Int *pn = get_facets(bot, M->ngon, M->indPG);
  
  for (E_Int i = 0; i < 4; i++) map[i] = pn[i];
  E_Int reorient = get_reorient(bot, i, normalIn_Py[0], M);
  if (reorient) std::swap(map[1], map[3]);

  E_Int lft, rgt, fro, bck;
  lft = rgt = fro = bck = -1;

  for (E_Int j = 1; j < 5; j++) {
    for (E_Int k = 0; k < 4; k++) common[k] = 0;

    E_Int face = pf[j];
    E_Int *pn = get_facets(face, M->ngon, M->indPG);

    for (E_Int k = 0; k < 3; k++) {
      E_Int point = pn[k];

      for (E_Int l = 0; l < 4; l++) {
        if (map[l] == point) {
          common[l] = 1;
          break;
        }
      }
    }

    if      (common[0] && common[3]) lft = face;
    else if (common[1] && common[2]) rgt = face;
    else if (common[1] && common[0]) fro = face;
    else                             bck = face;
  }

  assert(lft != -1);
  assert(rgt != -1);
  assert(fro != -1);
  assert(bck != -1);

  pf[1] = lft;
  pf[2] = rgt;
  pf[3] = fro;
  pf[4] = bck;
}

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
    E_Int ctype = M->cellTree->type_[i];

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
        break;
    }
  }
}

E_Int check_canon_cells(AMesh *M)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int type = M->cellTree->type_[i];
    switch (type) {
      case HEXA:
        assert(check_canon_hexa(i, M));
        break;
      case TETRA:
        assert(check_canon_tetra(i, M));
        break;
      case PENTA:
        assert(check_canon_penta(i, M));
        break;
      case PYRA:
        assert(check_canon_pyra(i, M));
        break;
      default:
        assert(0);
        break;
    }
  }

  return 1;
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

E_Int check_canon_tetra(E_Int cell, AMesh *M)
{
  E_Int NODES[4] = {-1, -1, -1, -1};

  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, i0, reorient, *pn, local[3];

  // BOT (In)
  face = pf[0];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn_T[0], M);
  Order_tri(local, pn, reorient, i0);
  for (E_Int i = 0; i < 3; i++) NODES[i] = local[i];

  // LFT (Out)
  face = pf[1];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 3);
  reorient = get_reorient(face, cell, normalIn_T[1], M);
  Order_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[2] == NODES[2]);
  NODES[3] = local[1];

  // RGT (In)
  face = pf[2];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[1], pn, 3);
  reorient = get_reorient(face, cell, normalIn_T[2], M);
  Order_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[3]);
  assert(local[2] == NODES[2]);

  // FRO (Out)
  face = pf[3];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 3);
  reorient = get_reorient(face, cell, normalIn_T[3], M);
  Order_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[1]);
  assert(local[2] == NODES[3]);

  return 1;
}

E_Int check_canon_hexa(E_Int cell, AMesh *M)
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

  return 1;
}

E_Int check_canon_penta(E_Int cell, AMesh *M)
{
  E_Int NODES[6] = {-1, -1, -1, -1, -1, -1};

  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, i0, reorient, *pn, local[4];

  // BOT
  face = pf[0];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn_Pe[0], M);
  Order_tri(local, pn, reorient, i0);
  for (E_Int i = 0; i < 3; i++) NODES[i] = local[i];

  // LFT (in)
  face = pf[2];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 4);
  assert(i0 != -1);
  reorient = get_reorient(face, cell, normalIn_Pe[2], M);
  Order_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[2]);
  NODES[5] = local[2];
  NODES[3] = local[3];

  // RGT (out)
  face = pf[3];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 4);
  reorient = get_reorient(face, cell, normalIn_Pe[3], M);
  Order_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[1]);
  assert(local[3] == NODES[3]);
  NODES[4] = local[2];

  // BCK (in)
  face = pf[4];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[2], pn, 4);
  reorient = get_reorient(face, cell, normalIn_Pe[4], M);
  Order_quad(local, pn, reorient, i0);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[1]);
  assert(local[2] == NODES[4]);
  assert(local[3] == NODES[5]);

  // TOP
  face = pf[1];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[3], pn, 3);
  reorient = get_reorient(face, cell, normalIn_Pe[1], M);
  Order_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[3]);
  assert(local[1] == NODES[4]);
  assert(local[2] == NODES[5]);

  return 1;
}

E_Int check_canon_pyra(E_Int cell, AMesh *M)
{
  E_Int NODES[5] = {-1, -1, -1, -1, -1};

  E_Int *pf = get_facets(cell, M->nface, M->indPH);

  E_Int face, i0, reorient, *pn, local[4];

  // BOT (in)
  face = pf[0];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = 0;
  reorient = get_reorient(face, cell, normalIn_Py[0], M);
  Order_quad(local, pn, reorient, i0);
  for (E_Int i = 0; i < 4; i++) NODES[i] = local[i];

  // LFT (in)
  face = pf[1];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[0], pn, 3);
  reorient = get_reorient(face, cell, normalIn_Py[1], M);
  Order_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[0]);
  assert(local[1] == NODES[3]);
  NODES[4] = local[2];

  // RGT (out)
  face = pf[2];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[1], pn, 3);
  reorient = get_reorient(face, cell, normalIn_Py[2], M);
  Order_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[2]);
  assert(local[2] == NODES[4]);

  // FRO (in)
  face = pf[3];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[1], pn, 3);
  reorient = get_reorient(face, cell, normalIn_Py[3], M);
  Order_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[1]);
  assert(local[1] == NODES[0]);
  assert(local[2] == NODES[4]);

  // BCK (out)
  face = pf[4];
  pn = get_facets(face, M->ngon, M->indPG);
  i0 = Get_pos(NODES[2], pn, 3);
  reorient = get_reorient(face, cell, normalIn_Py[4], M);
  Order_tri(local, pn, reorient, i0);
  assert(local[0] == NODES[2]);
  assert(local[1] == NODES[3]);
  assert(local[2] == NODES[4]);

  return 1;
}

E_Int is_internal_face(E_Int face, AMesh *M)
{
  return M->neigh[face] != -1;
}

AMesh *init_mesh(K_FLD::FldArrayI &cn, E_Float *px, E_Float *py,
  E_Float *pz, E_Int npts)
{
  E_Int *ngon = cn.getNGon();
  E_Int *indPG = cn.getIndPG();
  E_Int *nface = cn.getNFace();
  E_Int *indPH = cn.getIndPH();
  E_Int ncells = cn.getNElts();
  E_Int nfaces = cn.getNFaces();

  AMesh *M = new AMesh;

  M->npoints = npts;
  M->nfaces = nfaces;
  M->ncells = ncells;

  M->x = (E_Float *)XMALLOC(M->npoints * sizeof(E_Float));
  M->y = (E_Float *)XMALLOC(M->npoints * sizeof(E_Float));
  M->z = (E_Float *)XMALLOC(M->npoints * sizeof(E_Float));

  memcpy(M->x, px, M->npoints * sizeof(E_Float));
  memcpy(M->y, py, M->npoints * sizeof(E_Float));
  memcpy(M->z, pz, M->npoints * sizeof(E_Float));

  M->indPH = (E_Int *)XMALLOC((M->ncells+1) * sizeof(E_Int));
  M->indPG = (E_Int *)XMALLOC((M->nfaces+1) * sizeof(E_Int));

  memcpy(M->indPH, indPH, (M->ncells+1) * sizeof(E_Int));
  memcpy(M->indPG, indPG, (M->nfaces+1) * sizeof(E_Int));

  M->nface = (E_Int *)XMALLOC(M->indPH[M->ncells] * sizeof(E_Int));
  M->ngon  = (E_Int *)XMALLOC(M->indPG[M->nfaces] * sizeof(E_Int));

  E_Int *ptr = M->nface;

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int nf = -1;
    E_Int *pf = cn.getElt(i, nf, nface, indPH);
    for (E_Int j = 0; j < nf; j++)
      *ptr++ = pf[j]-1;
  }

  ptr = M->ngon;

  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int np = -1;
    E_Int *pn = cn.getFace(i, np, ngon, indPG);
    for (E_Int j = 0; j < np; j++)
      *ptr++ = pn[j]-1;
  }

  M->ecenter = new std::map<Edge, E_Int>;

  return M;
}

void get_full_cell(E_Int cell, AMesh *M, E_Int &nf, E_Int *pf)
{
  E_Int stride = -1;
  E_Int *FACES = get_cell(cell, stride, M->nface, M->indPH);
  nf = 0;

  for (E_Int i = 0; i < stride; i++) {
    E_Int face = FACES[i];

    E_Int flvl = M->faceTree->level(face);
    E_Int clvl = M->cellTree->level(cell);

    assert(flvl >= clvl);

    if (flvl == clvl) {
      pf[nf++] = face;
    } else {
      assert(master_face(face, M) == face);
      Children *children = M->faceTree->children(face);
      for (E_Int j = 0; j < children->n; j++)
        pf[nf++] = children->pc[j];
    }

    /*
    if (flvl == clvl) {
      pf[nf++] = face;
    } else {
      assert(master_face(face, M) == face);
      E_Int master = master_face(face, M);
      if (M->faceTree->state(face) == UNREFINED) {
        pf[nf++] = face;
      } else {
        Children *children = M->faceTree->children(master);
        if (children == NULL) {
          assert(0);
          pf[nf++] = master;
        } else {
          for (E_Int j = 0; j < children->n; j++)
            pf[nf++] = children->pc[j];
        }
      }
    }
    */
  }
}

// TODO(Imad): update comm patch faces

// Keep refined faces and their children, untouched faces and unrefined
// faces

void update_boundary_faces(AMesh *M)
{
  M->nbf = 0;
  for (E_Int i = 0; i < M->nbc; i++) {
    E_Int *ptlist = M->ptlists[i];
    
    // How many faces on current boundary have been refined/unrefined?
    E_Int new_bcsize = 0;
    for (E_Int j = 0; j < M->bcsizes[i]; j++) {
      E_Int face = ptlist[j];
      E_Int state = M->faceTree->state(face);

      if (state == UNTOUCHED)
        new_bcsize += 1;
      else if (state == REFINED)
        new_bcsize += M->faceTree->children(face)->n;
      else if (state == UNREFINED)
        new_bcsize += 1;
    }

    M->nbf += new_bcsize;

    E_Int *new_ptlist = (E_Int *)XMALLOC(new_bcsize * sizeof(E_Int));
    E_Int *ptr = new_ptlist;

    for (E_Int j = 0; j < M->bcsizes[i]; j++) {
      E_Int face = ptlist[j];
      E_Int state = M->faceTree->state(face);
      
      if (state == UNTOUCHED || state == UNREFINED) {
        *ptr++ = face;
      } else if (state == REFINED) {
        Children *children = M->faceTree->children(face);
        for (E_Int k = 0; k < children->n; k++) *ptr++ = children->pc[k];
      }
    }

    XFREE(M->ptlists[i]);

    M->ptlists[i] = new_ptlist;
    M->bcsizes[i] = new_bcsize;
  }

  M->nif = M->nfaces - M->nbf;
}

void renumber_mesh(AMesh *M, const std::vector<E_Int> &new_cells,
  const std::vector<E_Int> &new_faces, E_Int new_ncells, E_Int new_nfaces,
  E_Int sizeNFace, E_Int sizeNGon)
{
  // ngon
  E_Int *indPG = (E_Int *)XMALLOC((new_nfaces+1) * sizeof(E_Int));
  indPG[0] = 0;
  
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int new_face = new_faces[i];
    if (new_face == -1) continue;

    indPG[new_face+1] = get_stride(i, M->indPG);
  }

  for (E_Int i = 0; i < new_nfaces; i++) indPG[i+1] += indPG[i];
  assert(indPG[new_nfaces] == sizeNGon);

  E_Int *ngon = (E_Int *)XMALLOC(sizeNGon * sizeof(E_Int));
  
  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int new_face = new_faces[i];
    if (new_face == -1) continue;

    E_Int np = -1;
    E_Int *pn = get_face(i, np, M->ngon, M->indPG);

    E_Int *ptr = &ngon[indPG[new_face]];

    for (E_Int j = 0; j < np; j++)
      *ptr++ = pn[j];
  }



  // nface
  E_Int *indPH = (E_Int *)XMALLOC((new_ncells+1) * sizeof(E_Int));
  indPH[0] = 0;
  
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int new_cell = new_cells[i];
    if (new_cell == -1) continue;

    indPH[new_cell+1] = get_stride(i, M->indPH);
  }

  for (E_Int i = 0; i < new_ncells; i++) indPH[i+1] += indPH[i];
  assert(indPH[new_ncells] == sizeNFace);

  E_Int *nface = (E_Int *)XMALLOC(sizeNFace * sizeof(E_Int));
  
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int new_cell = new_cells[i];
    if (new_cell == -1) {
      continue;
    }

    E_Int nf = -1;
    E_Int *pf = get_cell(i, nf, M->nface, M->indPH);

    E_Int *ptr = &nface[indPH[new_cell]];

    for (E_Int j = 0; j < nf; j++) {
      *ptr++ = new_faces[pf[j]];
    }
  }

  // points
  E_Int new_npoints = 0;
  std::vector<E_Int> new_points(M->npoints, -1);
  for (E_Int i = 0; i < sizeNGon; i++) {
    E_Int point = ngon[i];
    if (new_points[point] == -1)
      new_points[point] = new_npoints++;
  }

  // ecenter
  Edge E;
  std::map<Edge, E_Int> *ecenter = new std::map<Edge, E_Int>;

  for (auto& e : *(M->ecenter)) {
    E_Int ec = e.second;
    
    // Delete GONE edge centers
    if (new_points[ec] == -1)
      continue;
    
    E_Int n0 = e.first.p0_;
    E_Int n1 = e.first.p1_;
    assert(new_points[n0] != -1);
    assert(new_points[n1] != -1);
    
    E.set(new_points[n0], new_points[n1]);

    ecenter->insert({E, new_points[ec]});
  }

  for (E_Int i = 0; i < sizeNGon; i++) {
    ngon[i] = new_points[ngon[i]];
  }

  // xyz
  E_Float *x = (E_Float *)XMALLOC(new_npoints * sizeof(E_Float));
  E_Float *y = (E_Float *)XMALLOC(new_npoints * sizeof(E_Float));
  E_Float *z = (E_Float *)XMALLOC(new_npoints * sizeof(E_Float));

  for (E_Int i = 0; i < M->npoints; i++) {
    if (new_points[i] == -1) continue;
    x[new_points[i]] = M->x[i];
    y[new_points[i]] = M->y[i];
    z[new_points[i]] = M->z[i];
  }

  // owner and neigh
  E_Int *owner = (E_Int *)XMALLOC(new_nfaces * sizeof(E_Int));
  E_Int *neigh = (E_Int *)XMALLOC(new_nfaces * sizeof(E_Int));
  
  for (E_Int i = 0; i < M->nfaces; i++) {
    if (new_faces[i] == -1) continue;

    owner[new_faces[i]] = new_cells[M->owner[i]];

    if (M->neigh[i] == -1) neigh[new_faces[i]] = -1;
    else {
      neigh[new_faces[i]] = new_cells[M->neigh[i]];
    }
  }

  // patches
  for (E_Int i = 0; i < M->nbc; i++) {
    E_Int *ptlist = M->ptlists[i];

    // Replace
    for (E_Int j = 0; j < M->bcsizes[i]; j++) {
      ptlist[j] = new_faces[ptlist[j]];
    }
  }

  // ref_data
  //puts("Warning: ref data is not renumbered");
  int *ref_data = (E_Int *)XMALLOC(new_ncells * sizeof(int));
  for (E_Int i = 0; i < M->ncells; i++) {
    if (new_cells[i] != -1)
      ref_data[new_cells[i]] = M->ref_data[i];
  }

  /*
  for (E_Int i = 0; i < new_ncells; i++) {
    for (E_Int j = indPH[i]; j < indPH[i+1]; j++)
      printf("%d ", nface[j]);
    printf("| ");
  }
  puts("");
  */

  // face/cell centers
  E_Float *fc = (E_Float *)XMALLOC(3*new_nfaces * sizeof(E_Float));
  E_Float *cx = (E_Float *)XMALLOC(new_ncells * sizeof(E_Float));
  E_Float *cy = (E_Float *)XMALLOC(new_ncells * sizeof(E_Float));
  E_Float *cz = (E_Float *)XMALLOC(new_ncells * sizeof(E_Float));

  for (E_Int i = 0; i < M->nfaces; i++) {
    if (new_faces[i] == -1) continue;

    E_Int new_face = new_faces[i];
    E_Float *np = &fc[3*new_face];
    E_Float *op = &M->fc[3*i];
    for (E_Int j = 0; j < 3; j++)
      np[j] = op[j];
  }

  for (E_Int i = 0; i < M->ncells; i++) {
    if (new_cells[i] == -1) continue;
    E_Int new_cell = new_cells[i];
    cx[new_cell] = M->cx[i];
    cy[new_cell] = M->cy[i];
    cz[new_cell] = M->cz[i];
  }

  // Free and replace
  XFREE(M->x);
  XFREE(M->y);
  XFREE(M->z);
  XFREE(M->nface);
  XFREE(M->indPH);
  XFREE(M->ngon);
  XFREE(M->indPG);
  XFREE(M->owner);
  XFREE(M->neigh);
  XFREE(M->ref_data); 
  delete M->ecenter;
  XFREE(M->fc);
  XFREE(M->cx);
  XFREE(M->cy);
  XFREE(M->cz);

  M->ncells = new_ncells;
  M->nfaces = new_nfaces;
  M->npoints = new_npoints;
  M->x = x;
  M->y = y;
  M->z = z;
  M->nface = nface;
  M->indPH = indPH;
  M->ngon = ngon;
  M->indPG = indPG;
  M->owner = owner;
  M->neigh = neigh;
  M->ref_data = ref_data;
  M->ecenter = ecenter; 
  M->fc = fc;
  M->cx = cx;
  M->cy = cy;
  M->cz = cz;
}
