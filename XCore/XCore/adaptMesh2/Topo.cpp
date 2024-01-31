#include "Proto.h"
#include <stack>

#define INTERNAL 0
#define EXTERNAL 1

#define DSMALL 1e-14

static
void set_PENTA_for_2D(E_Int cell, AMesh *M)
{
  E_Int nf = -1;
  E_Int *pf = get_cell(cell, nf, M->nface, M->indPH);
  assert(nf == 5);

  E_Int start = 0;
  for (; start < nf; start++) {
    E_Int *pn = get_face(cell, nf, M->ngon, M->indPG);

    // Compute normal to the triangle formed by the first 3 points
    E_Float e0[3] = {M->x[pn[1]] - M->x[pn[0]],
                     M->y[pn[1]] - M->y[pn[0]],
                     M->z[pn[1]] - M->z[pn[0]]};

    E_Float e1[3] = {M->x[pn[2]] - M->x[pn[0]],
                     M->y[pn[2]] - M->y[pn[0]],
                     M->z[pn[2]] - M->z[pn[0]]};

    E_Float n[3];
    K_MATH::cross(e0, e1, n);
    
    if (fabs(K_MATH::dot(n, M->mode_2D, 3) - 1.0) <= K_MATH::SMALL) {
      break;
    }
  }

  if (start == nf) {
    fprintf(stderr, "Couldn't align PENTA %d with prescribed 2D normal.\n",
      cell);
    exit(1);
  }

  if (M->faceTree->type(pf[start]) != TRI) {
    fprintf(stderr, "PENTA %d is aligned with prescribed 2D normal at"
      "face %d but it is not a triangle. This is weird.\n", cell, pf[start]);
    exit(1);
  }

  Right_shift(pf, start, nf);
}

E_Int check_face_aligned_with_vector(E_Int face, E_Float vec[3], AMesh *M)
{
  E_Int np = -1;
  E_Int *pn = get_face(face, np, M->ngon, M->indPG);

  // Compute normal to the triangle formed by the first 3 points
  E_Float e0[3] = {M->x[pn[1]] - M->x[pn[0]],
                   M->y[pn[1]] - M->y[pn[0]],
                   M->z[pn[1]] - M->z[pn[0]]};

  E_Float e1[3] = {M->x[pn[2]] - M->x[pn[0]],
                   M->y[pn[2]] - M->y[pn[0]],
                   M->z[pn[2]] - M->z[pn[0]]};

  E_Float n[3];
  K_MATH::cross(e0, e1, n);

  E_Float dp = fabs(K_MATH::dot(n, vec, 3)) / K_MATH::norm(n, 3);

  return fabs(dp - 1.0) <= 1.0e-6;
}

static
void set_HEXA_for_2D(E_Int cell, AMesh *M)
{
  E_Int nf = -1;
  E_Int *pf = get_cell(cell, nf, M->nface, M->indPH);

  E_Int start = 0;
  for (; start < nf; start++) {
    E_Int face = pf[start];
    if (check_face_aligned_with_vector(face, M->mode_2D, M))
      break;
  }

  if (start == nf) {
    fprintf(stderr, "Couldn't align HEXA %d with prescribed 2D normal.\n",
      cell);
    assert(0);
    exit(1);
  }

  Right_shift(pf, start, nf);
}

void set_cells_for_2D(AMesh *M)
{
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int type = M->cellTree->type(i);
    switch (type) {
      case HEXA:
        set_HEXA_for_2D(i, M);
        break;
      case PENTA:
        set_PENTA_for_2D(i, M);
        break;
      default:
        assert(0);
        break;
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

  if (M->mode_2D) {
    for (E_Int i = 0; i < M->ncells; i++) {
      E_Int type = M->cellTree->type(i);
      if (type != HEXA && type != PENTA) {
        fprintf(stderr, "AdaptMesh: 2D mode incompatible with cell type %s.\n",
          cell_type_to_string(type));
        exit(1);
      }
    }
  }

  return 0;
}


E_Int Get_pos(E_Int e, E_Int *pn, E_Int size)
{
  for (E_Int i = 0; i < size; i++) {
    if (pn[i] == e)
      return i;
  }
  return -1;
}

void Right_shift(E_Int *pn, E_Int pos, E_Int size)
{
  E_Int tmp[10];
  assert(size <= 10);
  for (E_Int i = 0; i < size; i++) tmp[i] = pn[i];
  for (E_Int i = 0; i < size; i++)
    pn[i] = tmp[(i+pos)%size];
}

void close_mesh(AMesh *M)
{
  Edge e;

  E_Int closed_ngon_size = 0;

  for (E_Int i = 0; i < M->nfaces; i++) {
    assert(M->faceTree->state(i) != GONE);

    E_Int np = -1;
    E_Int *pn = get_face(i, np, M->ngon, M->indPG);
    
    for (E_Int j = 0; j < np; j++) {
      e.set(pn[j], pn[(j+1)%np]);

      // Was this edge refined?
      auto search = M->ecenter->find(e);

      if (search != M->ecenter->end()) {
        // Yes
        closed_ngon_size += 2;
      } else {
        closed_ngon_size += 1;
      }
    }
  }

  M->closed_ngon = (E_Int *)XRESIZE(M->closed_ngon, closed_ngon_size * sizeof(E_Int));
  M->closed_indPG = (E_Int *)XRESIZE(M->closed_indPG, (M->nfaces+1) * sizeof(E_Int));
  memset(M->closed_ngon,  -1, closed_ngon_size * sizeof(E_Int));
  memset(M->closed_indPG, -1, (M->nfaces+1) * sizeof(E_Int));

  E_Int *closed_ngon = M->closed_ngon;
  E_Int *closed_indPG = M->closed_indPG;

  closed_indPG[0] = 0;

  E_Int *ptr = closed_ngon;

  for (E_Int i = 0; i < M->nfaces; i++) {
    E_Int np = -1;
    E_Int *pn = get_face(i, np, M->ngon, M->indPG);

    E_Int stride = 0;

    for (E_Int j = 0; j < np; j++) {
      E_Int p0 = pn[j];
      E_Int p1 = pn[(j+1)%np];

      e.set(p0, p1);

      auto search = M->ecenter->find(e);

      if (search != M->ecenter->end()) {
        *ptr++ = p0;
        *ptr++ = search->second;
        stride += 2;
      } else {
        *ptr++ = p0;
        stride += 1;
      }
    }

    closed_indPG[i+1] = stride;
  }

  for (E_Int i = 0; i < M->nfaces; i++)
    closed_indPG[i+1] += closed_indPG[i];
  
  assert(closed_indPG[M->nfaces] == closed_ngon_size);
}


static
void Compute_cell_volume(E_Int cell, AMesh *M, E_Float *x, E_Float *y,
  E_Float *z, E_Float &vol, E_Int refIdx)
{
  E_Int *nface = M->nface;
  E_Int *indPH = M->indPH;
  E_Int *ngon = M->ngon;
  E_Int *indPG = M->indPG;

  // Orient the faces coherently
  std::vector<E_Int> NGON;
  std::vector<E_Int> INDPG(1, 0);
  E_Int stride = -1;
  E_Int *pf = get_cell(cell, stride, nface, indPH);

  for (E_Int i = 0; i < stride; i++) {
    E_Int face = pf[i];
    E_Int np = -1;
    E_Int *pn = get_face(face, np, ngon, indPG);
    INDPG.push_back(np);
    for (E_Int j = 0; j < np; j++)
      NGON.push_back(pn[j]);
  }

  for (E_Int i = 0; i < stride; i++)
    INDPG[i+1] += INDPG[i];

  // Fix orientation of first face
  std::vector<E_Int> orient(stride);
  orient[refIdx] = 1;
  std::vector<E_Int> neis(NGON.size());
  K_CONNECT::build_face_neighbourhood(NGON, INDPG, neis);
  K_CONNECT::reversi_connex(&NGON[0], &INDPG[0], stride, &neis[0], refIdx, orient);

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
    K_METRIC::compute_face_center_and_area(face, np, pn, x, y, z, fc, fa);
  }
  
  // Estimate cell centroid as average of face centers
  E_Float cc[3] = {0,0,0};
  for (E_Int i = 0; i < stride; i++) {
    E_Float *fc = &faceCenters[3*i];
    for (E_Int j = 0; j < 3; j++)
      cc[j] += fc[j];
  }
  for (E_Int i = 0; i < 3; i++)
    cc[i] /= stride;

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

static
void Flag_and_get_external_faces(AMesh *M, std::vector<E_Int> &fflags,
  std::vector<E_Int> &efaces)
{
  std::vector<E_Int> face_count(M->nfaces, 0);

  // Loop through the elements and increment face_count
  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int stride = -1;
    E_Int *pf = get_cell(i, stride, M->nface, M->indPH);
    for (E_Int j = 0; j < stride; ++j)
      face_count[pf[j]]++;
  }

  // External faces are those with a count equal to 1
  fflags.resize(M->nfaces);
  for (E_Int i = 0; i < M->nfaces; i++) {
    if (face_count[i] == 1) {
      fflags[i] = EXTERNAL;
      efaces.push_back(i);
    } else {
      fflags[i] = INTERNAL;
    }
  }
}

static
E_Int _Orient_boundary
(
  AMesh *M,
  E_Float *x, E_Float *y, E_Float *z,
  E_Int ncells,
  E_Int *efadj, E_Int *efxadj, E_Int nefaces,
  E_Int *fneis, E_Int *efaces, std::vector<E_Int> &forient,
  const std::vector<E_Int> &cflags, const std::vector<E_Int> &fflags,
  E_Int *cells)
{
  // Look for a cell whose volume is "definitely well computed"
  E_Float cvol = 0.0;
  E_Int seed = -1;
  E_Int refPG = -1;
  E_Int refIdx = -1;

  while (++seed < ncells) {
    if (cflags[seed] != EXTERNAL) continue;

    E_Int cid = (cells != NULL) ? cells[seed] : seed;

    E_Int stride = -1;
    E_Int *pf = get_cell(cid, stride, M->nface, M->indPH);
    refPG = -1;
    E_Int local_idx = -1;
    for (E_Int j = 0; j < stride; j++) {
      E_Int face = pf[j];
      if (fflags[face] == EXTERNAL) {
        refPG = face;
        local_idx = j;
        break;
      }
    }

    if (refPG == -1) {
      fprintf(stderr, "orient_boundary(): couldn't find an external face within external cell %d\n", cid);
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
      fprintf(stderr, "orient_boundary(): couldn't find reference face %d in external faces list\n", refPG);
      return 1;
    }

    // Set orientation of refPG to +1.
    // Reorient seed's faces based on orientation of refPG.
    // Compute cvol, the volume of seed.
    // If cvol > 0, orientation of all faces including refPG, is outwards
    // Otherwise, set orientation of refPG to -1.

    Compute_cell_volume(cid, M, x, y, z, cvol, local_idx);

    if (fabs(cvol) < DSMALL) continue;

    // set reference orientation of refPG and exit
    forient[refIdx] = (cvol > 0.0) ? 1 : -1;

    break;
  }

   if (seed >= ncells) {
    fprintf(stderr, "orient_boundary_ngon(): couldn't find reference polyhedron\n");
    assert(0);
    return 1;
  }

  // propagate
  K_CONNECT::reversi_connex(efadj, efxadj, nefaces, fneis, refIdx, forient);

  return 0;
}

static
void Extract_nface_of_kept_pgs(AMesh *M, const std::vector<bool> &kept_pgs,
  std::vector<E_Int> &NFACE, std::vector<E_Int> &xadj,
  std::vector<E_Int> &cells)
{
  E_Int *nface = M->nface;
  E_Int *indPH = M->indPH;
  E_Int ncells = M->ncells;

  NFACE.clear();
  xadj.resize(1, 0);
  cells.clear();

  for (E_Int i = 0; i < ncells; i++) {
    E_Int stride = -1;
    E_Int *pf = get_cell(i, stride, nface, indPH);
    bool keep = false;
    for (E_Int j = 0; j < stride && !keep; j++)
      keep = kept_pgs[pf[j]];
    if (keep) {
      cells.push_back(i);
      xadj.push_back(stride);
      for (E_Int j = 0; j < stride; j++)
        NFACE.push_back(pf[j]);
    }
  }

  for (size_t i = 0; i < xadj.size(); i++)
    xadj[i+1] += xadj[i];
}

static
void Flag_marked_external_cells(AMesh *M, const std::vector<E_Int> &cells,
  const std::vector<E_Int> &fflags, std::vector<E_Int> &cflags)
{
  E_Int *nface = M->nface;
  E_Int *indPH = M->indPH;

  // External cells are those with at least one external face
  cflags.resize(cells.size(), INTERNAL);
  for (size_t i = 0; i < cells.size(); i++) {
    E_Int cell = cells[i];
    E_Int stride = -1;
    E_Int *pf = get_cell(cell, stride, nface, indPH);
    for (E_Int j = 0; j < stride; j++) {
      E_Int face = pf[j];
      if (fflags[face] == EXTERNAL) {
        cflags[i] = EXTERNAL;
        break;
      }
    }
  }
}

static
void Flag_all_external_cells(AMesh *M, const std::vector<E_Int> &fflags,
  std::vector<E_Int> &cflags)
{
  E_Int *nface = M->nface;
  E_Int *indPH = M->indPH;
  E_Int ncells = M->ncells;

  // External cells are those with at least one external face
  cflags.resize(ncells, INTERNAL);
  for (E_Int i = 0; i < ncells; i++) {
    E_Int stride = -1;
    E_Int *pf = get_cell(i, stride, nface, indPH);
    for (E_Int j = 0; j < stride; j++) {
      E_Int face = pf[j];
      if (fflags[face] == EXTERNAL) {
        cflags[i] = EXTERNAL;
        break;
      }
    }
  }
}

E_Int Orient_boundary(AMesh *M)
{
  E_Int nfaces = M->nfaces;
  E_Int *indPG = M->indPG;
  E_Int *ngon = M->ngon;
  E_Int ncells = M->ncells;
  E_Float *x = M->x;
  E_Float *y = M->y;
  E_Float *z = M->z;

  // flag external cells and faces
  std::vector<E_Int> fflags, efaces;
  Flag_and_get_external_faces(M, fflags, efaces);

  // extract external faces connectivity
  std::vector<E_Int> fadj;
  std::vector<E_Int> xadj(1, 0);
  for (E_Int i = 0; i < nfaces; i++) {
    if (fflags[i] == EXTERNAL) {
      E_Int stride = -1;;
      E_Int *pn = get_face(i, stride, ngon, indPG);
      xadj.push_back(stride);
      for (E_Int j = 0; j < stride; j++)
        fadj.push_back(pn[j]);
    }
  }

  E_Int nefaces = (E_Int)efaces.size();

  for (E_Int i = 0; i < nefaces; i++)
    xadj[i+1] += xadj[i];

  // build skin neighbourhood
  std::vector<E_Int> fneighbours;
  K_CONNECT::build_face_neighbourhood(fadj, xadj, fneighbours);

  // color the faces by connex part
  std::vector<E_Int> colors(xadj.size()-1);
  E_Int nconnex = K_CONNECT::colorConnexParts(&fneighbours[0], &xadj[0],
    nefaces, &colors[0]);

  //printf("orient_boundary(): connex parts: %d\n", nconnex);

  assert(efaces.size() == xadj.size()-1);
  std::vector<E_Int> forient(nefaces, 0);
  std::vector<E_Int> cflags;
  E_Int ret = 0;
  if (nconnex > 1) {
    // extract nconnex nface-ngon for separate orientation
    for (E_Int color = 0; color < nconnex; color++) {
      std::vector<bool> keep_pgs(nfaces, false);
      for (E_Int i = 0; i < nefaces; i++) {
        keep_pgs[efaces[i]] = (colors[i] == color);
      }
      // extract nface corresponding to kept faces
      std::vector<E_Int> NFACE, cxadj(1, 0), cells;
      Extract_nface_of_kept_pgs(M, keep_pgs, NFACE, cxadj, cells);

      std::vector<E_Int> cflags;
      Flag_marked_external_cells(M, cells, fflags, cflags);

      ret |= _Orient_boundary(M, x, y, z, (E_Int)cells.size(),
        &fadj[0], &xadj[0], nefaces, &fneighbours[0], &efaces[0], forient,
        cflags, fflags, &cells[0]);
    }
  } else {
    std::vector<E_Int> cflags;
    Flag_all_external_cells(M, fflags, cflags);
    ret = _Orient_boundary(M, x, y, z, ncells, &fadj[0], &xadj[0],
      nefaces, &fneighbours[0], &efaces[0], forient, cflags, fflags, NULL);
  }

  // Apply orientation
  //E_Int nrev = 0;
  for (E_Int i = 0; i < nefaces; i++) {
    if (forient[i] == -1) {
      E_Int face = efaces[i]; // 0-based
      E_Int stride = -1;
      E_Int *pn = get_face(face, stride, ngon, indPG);
      std::reverse(pn+1, pn+stride);
      //nrev++;
    }
  }
  //printf("orient_boundary(): reversed %d faces\n", nrev);

  return ret;
}

static
void Build_cell_neighbourhood(AMesh *M, std::vector<E_Int>& neighbours,
  std::vector<E_Int> &xadj)
{
  E_Int nfaces = M->nfaces;
  E_Int ncells = M->ncells;
  E_Int *indPH = M->indPH;
  E_Int *nface = M->nface;

  xadj.resize(ncells+1);
  xadj[0] = 0;

  // TODO(Imad): this is a safe resize
  E_Int *ptr = &xadj[0]+1;
  for (E_Int i = 0; i < ncells; i++) {
    E_Int stride = -1;
    get_cell(i, stride, nface, indPH);
    *ptr++ = stride;
  }

  for (E_Int i = 0; i < ncells; i++) xadj[i+1] += xadj[i];

  E_Int sz = xadj[ncells];
  neighbours.resize(sz, -1);

  std::vector<E_Int> neigh(nfaces, -1);

  E_Int count = 0;
  while (count++ != 2) {
    for (E_Int i = 0; i < ncells; i++) {
      E_Int stride = -1;
      E_Int *pf = get_cell(i, stride, nface, indPH);
      E_Int *pn = &neighbours[xadj[i]];
      for (E_Int j = 0; j < stride; j++) {
        E_Int face = pf[j];
        E_Int &nei = neigh[face];
        E_Int &Kn = pn[j];
        if (nei != -1 && nei != i)
          Kn = nei;
        neigh[face] = i;
      }
    }
  }
}

E_Int Build_own_nei(AMesh *M)
{
  E_Int ncells = M->ncells;
  E_Int *nface = M->nface;
  E_Int *ngon = M->ngon;
  E_Int *indPH = M->indPH;
  E_Int *indPG = M->indPG;
  E_Int *owner = M->owner;
  E_Int *neigh = M->neigh;

  std::vector<E_Int> neighbours, xadj;
  Build_cell_neighbourhood(M, neighbours, xadj);

  std::vector<E_Int> exPH(ncells, -1);
  for (E_Int i = 0; i < ncells; i++) {
    E_Int stride = -1;
    E_Int *pf = get_cell(i, stride, nface, indPH);
    E_Int *pn = &neighbours[xadj[i]];
    for (E_Int j = 0; j < stride; j++) {
      if (pn[j] == -1) {
        owner[pf[j]] = i;
        exPH[i] = pf[j]+1;
        break;
      }
    }
  }

  // look for first external cell
  std::vector<E_Int> processed(ncells, 0);

  E_Int nconnex = 0;

  E_Int seed = 0;

  std::stack<E_Int> cpool;

   while (1) {
    while ((seed < ncells) && (processed[seed] || exPH[seed] == -1))
      ++seed;

    if (seed >= ncells)
      break;

    nconnex++;

    cpool.push(seed);

    while (!cpool.empty()) {
      E_Int cell = cpool.top();
      cpool.pop();

      if (processed[cell])
        continue;

      processed[cell] = 1;

      // build faces neighbourhood based on shared edges' nodes order
      E_Int stride = -1;
      E_Int *pf = get_cell(cell, stride, nface, indPH);
      std::vector<E_Int> oids;
      std::vector<E_Int> orient(stride, 1);
      std::vector<E_Int> pgs;
      std::vector<E_Int> xpgs(1, 0);
      for (E_Int i = 0; i < stride; i++) {
        E_Int face = pf[i];
        E_Int np = -1;
        E_Int *pn = get_face(face, np, ngon, indPG);
        for (E_Int j = 0; j < np; j++)
          pgs.push_back(pn[j]);
        xpgs.push_back(np);
        oids.push_back(face+1);
      }

      for (E_Int i = 0; i < stride; i++)
        xpgs[i+1] += xpgs[i];

      std::vector<E_Int> PGneighbours(pgs.size());
      K_CONNECT::build_face_neighbourhood(pgs, xpgs, PGneighbours);

      E_Int revers = 0;

      // reference face is the external face
      E_Int PGref = exPH[cell];
      // face can be negative
      if (PGref < 0) {
        revers = 1;
        PGref = -PGref;
      }

      // find reference face index in oids
      E_Int iref = -1;
      for (size_t i = 0; i < oids.size(); i++) {
        if (PGref == oids[i]) {
          iref = i;
          break;
        }
      }
      assert(iref != -1);

      // set orientation of face if prescribed
      if (revers)
        orient[iref] = -1;

      // all connected faces must follow the orientation of the reference face
      K_CONNECT::reversi_connex(&pgs[0], &xpgs[0], stride, &PGneighbours[0],
        iref, orient);

       // set the owner and neighbour of the faces
      E_Int *pn = &neighbours[xadj[cell]];
      for (E_Int i = 0; i < stride; i++) {
        E_Int face = pf[i];
        E_Int nei = pn[i];
        assert(nei < ncells && nei >= -1);

        owner[face] = cell;
        neigh[face] = nei;

        if (nei == -1)
          continue;

        // set the reference face for neighbour
        exPH[nei] = -(face+1);

        if (orient[i] == -1) {
          std::swap(owner[face], neigh[face]);
          exPH[nei] = face+1;
        }

        if (!processed[nei])
          cpool.push(nei);
      }
    }
  }

  //printf("build_parent_elements(): connex parts: %d\n", nconnex);

  return 0;
}

E_Int check_closed_cell(E_Int cell, E_Int *pf, E_Int stride, AMesh *M)
{
  // Closed cell: edge count = 2 for all edges
  Edge E;
  std::map<Edge, E_Int> edgeCount;
  for (E_Int i = 0; i < stride; i++) {
    E_Int face = pf[i];
    E_Int np = -1;
    E_Int *pn = get_face(face, np, M->ngon, M->indPG);
    for (E_Int j = 0; j < np; j++) {
      E_Int p0 = pn[j];
      E_Int p1 = pn[(j+1)%np];
      E.set(p0, p1);
      edgeCount[E]++;
    }
  }

  for (auto& ec : edgeCount) {
    if (ec.second != 2) {
      fprintf(stderr, "Warning: Cell %d is not closed.\n", cell); 
      return 0;
    }
  }

  return 1;
}
