#include "xcore.h"
#include <mpi.h>
#include <unordered_map>

#include "scotch/ptscotch.h"

#define EXIT \
  do { \
    MPI_Finalize(); \
    exit(0); \
  } while (0);

static
E_Int get_proc(E_Int element, const std::vector<E_Int> &distribution,
  E_Int nproc)
{
  for (E_Int j = 0; j < nproc; j++) {
    if (element >= distribution[j] && element < distribution[j+1])
      return j;
  }
  printf("\nWarning: could not find distribution of element %d\n", element);
  assert(0);
  return -1;
}

struct Poly {
  E_Int n[8];
  E_Int stride;

  Poly(E_Int s)
  : stride(s)
  {}

  void set(E_Int *p)
  {
    std::sort(p, p+stride);
    for (E_Int i = 0; i < stride; i++)
      n[i] = p[i];
  }

  bool operator==(const Poly &P) const
  {
    for (E_Int i = 0; i < stride; i++) {
      if (n[i] != P.n[i])
        return false;
    }
    
    return true;
  }
};

struct Poly_hash {
  uint32_t hash(uint32_t val, uint32_t seed) const
  {
    uint32_t HASH = seed;
    HASH += val;
    HASH += HASH << 10;
    HASH ^= HASH >> 6;
    return HASH;
  }

  uint32_t operator()(const Poly &P) const
  {
    uint32_t res = 0;
    for (E_Int i = 0; i < P.stride; i++)
      res = hash(P.n[i], res);
    
    res += res << 3;
    res ^= res >> 11;
    res += res << 15;

    return res;
  }
};

static
void add_facet(const Poly &P, E_Int &nfaces,
  std::unordered_map<Poly, E_Int, Poly_hash> &faces, std::vector<E_Int> &NFACE,
  std::vector<Poly> &face2Poly)
{
  auto search = faces.find(P);
  if (search == faces.end()) {
    faces[P] = nfaces;
    NFACE.push_back(nfaces);
    face2Poly.push_back(P);
    nfaces += 1;
  } else {
    NFACE.push_back(search->second);
  }
}

static
void make_facets_HEXA(E_Int *cn, E_Int ncells, E_Int &nfaces,
  std::vector<E_Int> &NFACE, std::unordered_map<Poly, E_Int, Poly_hash> &faces,
  std::vector<Poly> &f2Q)
{
  nfaces = 1;
  NFACE.clear();
  faces.clear();

  Poly Quad(4);
  E_Int n[4];

  for (E_Int i = 0; i < ncells; i++) {
    E_Int *pn = &cn[8*i];
    
    // bot
    n[0] = pn[3]; n[1] = pn[0]; n[2] = pn[1]; n[3] = pn[2];
    Quad.set(n);
    add_facet(Quad, nfaces, faces, NFACE, f2Q);
    
    // top
    n[0] = pn[7]; n[1] = pn[4]; n[2] = pn[5]; n[3] = pn[6];
    Quad.set(n);
    add_facet(Quad, nfaces, faces, NFACE, f2Q);
    
    // left
    n[0] = pn[3]; n[1] = pn[2]; n[2] = pn[6]; n[3] = pn[7];
    Quad.set(n);
    add_facet(Quad, nfaces, faces, NFACE, f2Q);
    
    // right
    n[0] = pn[0]; n[1] = pn[1]; n[2] = pn[5]; n[3] = pn[4];
    Quad.set(n);
    add_facet(Quad, nfaces, faces, NFACE, f2Q);
    
    // front
    n[0] = pn[0]; n[1] = pn[3]; n[2] = pn[7]; n[3] = pn[4];
    Quad.set(n);
    add_facet(Quad, nfaces, faces, NFACE, f2Q);
    
    // back
    n[0] = pn[1]; n[1] = pn[2]; n[2] = pn[6]; n[3] = pn[5];
    Quad.set(n);
    add_facet(Quad, nfaces, faces, NFACE, f2Q);
  }

  assert((E_Int)NFACE.size() == 6*ncells);

  nfaces--;

  assert(nfaces == (E_Int)faces.size());
}

static
void make_facets_TETRA(E_Int *cn, E_Int ncells, E_Int &nfaces,
  std::vector<E_Int> &NFACE, std::unordered_map<Poly, E_Int, Poly_hash> &faces,
  std::vector<Poly> &f2Q)
{
  nfaces = 1;
  NFACE.clear();
  faces.clear();

  Poly Tri(3);
  E_Int n[3];

  for (E_Int i = 0; i < ncells; i++) {
    E_Int *pn = &cn[4*i];
    
    // bot
    n[0] = pn[0];
    n[1] = pn[1];
    n[2] = pn[2];
    Tri.set(n);
    add_facet(Tri, nfaces, faces, NFACE, f2Q);

    // right
    n[0] = pn[2];
    n[1] = pn[1];
    n[2] = pn[3];
    Tri.set(n);
    add_facet(Tri, nfaces, faces, NFACE, f2Q);

    // front
    n[0] = pn[1];
    n[1] = pn[0];
    n[2] = pn[3];
    Tri.set(n);
    add_facet(Tri, nfaces, faces, NFACE, f2Q);

    // left
    n[0] = pn[0];
    n[1] = pn[2];
    n[2] = pn[3];
    Tri.set(n);
    add_facet(Tri, nfaces, faces, NFACE, f2Q);
  }

  assert((E_Int)NFACE.size() == 4*ncells);

  nfaces--;

  assert(nfaces == (E_Int)faces.size());
}

static
void make_facets_PENTA(E_Int *cn, E_Int ncells, E_Int &nfaces,
  std::vector<E_Int> &NFACE, std::unordered_map<Poly, E_Int, Poly_hash> &faces,
  std::vector<Poly> &f2Q)
{
  nfaces = 1;
  NFACE.clear();
  faces.clear();

  Poly Tri(3), Quad(4);
  E_Int n[4];

  for (E_Int i = 0; i < ncells; i++) {
    E_Int *pn = &cn[6*i];
    
    // bot
    n[0] = pn[0]; n[1] = pn[1]; n[2] = pn[2];
    Tri.set(n);
    add_facet(Tri, nfaces, faces, NFACE, f2Q);

    // top
    n[0] = pn[3]; n[1] = pn[4]; n[2] = pn[5];
    Tri.set(n);
    add_facet(Tri, nfaces, faces, NFACE, f2Q);
    
    // left
    n[0] = pn[0]; n[1] = pn[3]; n[2] = pn[5]; n[3] = pn[2];
    Quad.set(n);
    add_facet(Quad, nfaces, faces, NFACE, f2Q);

    // right
    n[0] = pn[0]; n[1] = pn[1]; n[2] = pn[4]; n[3] = pn[3];
    Quad.set(n);
    add_facet(Quad, nfaces, faces, NFACE, f2Q);

    // back
    n[0] = pn[1]; n[1] = pn[2]; n[2] = pn[5]; n[3] = pn[4];
    Quad.set(n);
    add_facet(Quad, nfaces, faces, NFACE, f2Q);
  }

  assert((E_Int)NFACE.size() == 5*ncells);

  nfaces--;

  assert(nfaces == (E_Int)faces.size());
}

static
void make_facets_PYRA(E_Int *cn, E_Int ncells, E_Int &nfaces,
  std::vector<E_Int> &NFACE, std::unordered_map<Poly, E_Int, Poly_hash> &faces,
  std::vector<Poly> &f2Q)
{
  nfaces = 1;
  NFACE.clear();
  faces.clear();

  Poly Tri(3), Quad(4);
  E_Int n[5];

  for (E_Int i = 0; i < ncells; i++) {
    E_Int *pn = &cn[5*i];
    
    // bot
    n[0] = pn[0];
    n[1] = pn[1];
    n[2] = pn[2];
    n[3] = pn[3];
    Quad.set(n);
    add_facet(Quad, nfaces, faces, NFACE, f2Q);

    // left
    n[0] = pn[0];
    n[1] = pn[3];
    n[2] = pn[4];
    Tri.set(n);
    add_facet(Tri, nfaces, faces, NFACE, f2Q);

    // right
    n[0] = pn[2];
    n[1] = pn[1];
    n[2] = pn[4];
    Tri.set(n);
    add_facet(Tri, nfaces, faces, NFACE, f2Q);

    // front
    n[0] = pn[1];
    n[1] = pn[0];
    n[2] = pn[4];
    Tri.set(n);
    add_facet(Tri, nfaces, faces, NFACE, f2Q);

    // back
    n[0] = pn[2];
    n[1] = pn[3];
    n[2] = pn[4];
    Tri.set(n);
    add_facet(Tri, nfaces, faces, NFACE, f2Q);
  }

  assert((E_Int)NFACE.size() == 5*ncells);

  nfaces--;

  assert(nfaces == (E_Int)faces.size());
}

static
E_Int make_facets(E_Int *cn, E_Int ncells, E_Int stride,
  const char *elemName, E_Int &nfaces, std::vector<E_Int> &NFACE,
  std::unordered_map<Poly, E_Int, Poly_hash> &faces,
  std::vector<Poly> &f2Q)
{
  if (stride == 5) {
    make_facets_PYRA(cn, ncells, nfaces, NFACE, faces, f2Q);
  } else if (stride == 4 && strcmp(elemName, "TETRA") == 0) {
    make_facets_TETRA(cn, ncells, nfaces, NFACE, faces, f2Q);
  } else if (stride == 6) {
    make_facets_PENTA(cn, ncells, nfaces, NFACE, faces, f2Q);
  } else if (stride == 8) {
    make_facets_HEXA(cn, ncells, nfaces, NFACE, faces, f2Q);
  } else {
    fprintf(stderr, "chunk2part_elt(): implemented for TETRA, PYRA, PENTA and HEXA\n");
    return 1;
  }

  return 0;
}

struct comm_patch {
  E_Int nei; // nei proc
  std::vector<E_Int> local;
  std::vector<E_Int> remote;

  comm_patch(E_Int proc)
  : nei(proc)
  {}

  void add(E_Int lf, E_Int rf)
  {
    local.push_back(lf);
    remote.push_back(rf);
  }
};

PyObject *K_XCORE::chunk2part_elt(PyObject *self, PyObject *args)
{
  PyObject *XYZ, *STRIDE, *CN;
  char *eltName;

  if (!PyArg_ParseTuple(args, "OsOO", &XYZ, &eltName, &STRIDE, &CN)) {
    PyErr_SetString(PyExc_ValueError, "chunk2part_elt(): bad input");
    return NULL;
  }

  // deduce fstride: number of faces per elements
  E_Int fstride = -1;
  if (strcmp(eltName, "HEXA") == 0)
    fstride = 6;
  else if (strcmp(eltName, "PENTA") == 0)
    fstride = 5;
  else if (strcmp(eltName, "TETRA") == 0)
    fstride = 4;
  else if (strcmp(eltName, "PYRA") == 0)
    fstride = 5;
  else {
    PyErr_SetString(PyExc_TypeError, "chunk2part_elt(): only for HEXA, PENTA and TETRA elements.");
    return NULL;
  }

  E_Int rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  PyObject *Xi, *Yi, *Zi;
  E_Int res;

  // coordinates
  E_Float *X, *Y, *Z;
  E_Int npoints, nfld;

  Xi = PyList_GetItem(XYZ, 0);
  res = K_NUMPY::getFromNumpyArray(Xi, X, npoints, nfld, true);
  assert(res == 1);

  Yi = PyList_GetItem(XYZ, 1);
  res = K_NUMPY::getFromNumpyArray(Yi, Y, npoints, nfld, true);
  assert(res == 1);

  Zi = PyList_GetItem(XYZ, 2);
  res = K_NUMPY::getFromNumpyArray(Zi, Z, npoints, nfld, true);
  assert(res == 1);

  // element stride
  if (!PyLong_Check(STRIDE)) {
    PyErr_SetString(PyExc_TypeError, "chunk2part_elt(): bad stride.");
    return NULL;
  }
  E_Int stride;
  stride = PyLong_AsLong(STRIDE);

  // connectivity array
  E_Int *cn;
  E_Int cnsize;
  res = K_NUMPY::getFromNumpyArray(CN, cn, cnsize, nfld, true);
  assert(res == 1);
  E_Int ncells = cnsize / stride;

  // build cell and points distributions
  std::vector<E_Int> cdist((nproc+1));
  std::vector<E_Int> pdist((nproc+1));
  cdist[0] = pdist[0] = 0;

  MPI_Allgather(&npoints, 1, MPI_INT, &pdist[0]+1, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&ncells, 1, MPI_INT, &cdist[0]+1, 1, MPI_INT, MPI_COMM_WORLD);
  for (E_Int i = 0; i < nproc; i++) {
    cdist[i+1] += cdist[i];
    pdist[i+1] += pdist[i];
  }

  E_Int gncells = 0;
  MPI_Allreduce(&ncells, &gncells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  assert(gncells == cdist[nproc]);
  if (rank == 0)
    printf("Total number of %s: %d (Average: %d)\n", eltName, gncells, gncells/nproc);


  // make cellFaces connectivity
  E_Int nfaces = 0;
  std::vector<E_Int> NFACE;
  std::unordered_map<Poly, E_Int, Poly_hash> faces;
  std::vector<Poly> face2Poly;
  res = make_facets(cn, ncells, stride, eltName, nfaces, NFACE, faces, face2Poly);
  if (res == 1)
    return NULL;

  E_Int firstFace;
  MPI_Scan(&nfaces, &firstFace, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  firstFace -= nfaces;

  // adjust global index of faces
  for (auto &face : faces)
    face.second = face.second + firstFace;
  
  // adjust NFACE
  for (auto &face : NFACE)
    face = face + firstFace;

  // make local ParentElement
  E_Int firstCell = cdist[rank];

  std::unordered_map<E_Int, std::vector<E_Int>> lPE;

  for (E_Int i = 0; i < ncells; i++) {
    E_Int *pf = &NFACE[fstride*i];
    for (E_Int j = 0; j < fstride; j++) {
      lPE[pf[j]].push_back(i+firstCell);
    }
  }

  // potential pface data to exchange with every proc
  // Note(Imad): is it sufficient to send the facet hash only?
  std::vector<E_Int> SEND;

  for (auto &face : lPE) {
    if (face.second.size() == 1) {
      assert(face.first > 0);
      SEND.push_back(face.first);
      const auto &P = face2Poly[face.first - firstFace - 1];
      SEND.push_back(P.stride);
      for (E_Int j = 0; j < P.stride; j++)
        SEND.push_back(P.n[j]);
    }
  }

  E_Int ssize = SEND.size();
  std::vector<E_Int> RCOUNT(nproc);
  MPI_Allgather(&ssize, 1, MPI_INT, &RCOUNT[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> RDIST(nproc+1);
  RDIST[0] = 0;
  for (E_Int i = 0; i < nproc; i++)
    RDIST[i+1] = RDIST[i] + RCOUNT[i];
  
  std::vector<E_Int> RECV(RDIST[nproc], -1);

  MPI_Allgatherv(&SEND[0], (E_Int)SEND.size(), MPI_INT,
                 &RECV[0], &RCOUNT[0], &RDIST[0], MPI_INT,
                 MPI_COMM_WORLD);


  // look for received faces in face table
  std::vector<E_Int> spfaces;
  std::vector<E_Int> spcount(nproc, 0), rpcount(nproc);

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &RECV[RDIST[i]];
    if (i == rank) continue;
    for (E_Int j = 0; j < RCOUNT[i];) {
      E_Int face = ptr[j++];
      E_Int stride = ptr[j++];
      Poly P(stride);
      E_Int n[8];
      for (E_Int k = 0; k < stride; k++)
        n[k] = ptr[j++];
      P.set(n);
      auto search = faces.find(P);

      if (search != faces.end()) {
        spfaces.push_back(face);
        spfaces.push_back(search->second);
        spcount[i] += 2;
      }
    }
  }

  MPI_Alltoall(&spcount[0], 1, MPI_INT, &rpcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> spdist(nproc+1), rpdist(nproc+1);
  spdist[0] = rpdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    spdist[i+1] = spdist[i] + spcount[i];
    rpdist[i+1] = rpdist[i] + rpcount[i];
  }

  std::vector<E_Int> rpfaces(rpdist[nproc]);
  MPI_Alltoallv(&spfaces[0], &spcount[0], &spdist[0], MPI_INT,
                &rpfaces[0], &rpcount[0], &rpdist[0], MPI_INT,
                MPI_COMM_WORLD);

  // comm graph
  std::vector<comm_patch> cpatches;
  E_Int ndup = 0;
  std::set<E_Int> dups;

  for (E_Int i = 0; i < nproc; i++) {
    if (rank == i) continue;
    E_Int *ptr = &rpfaces[rpdist[i]];

    comm_patch cpatch(i);
    for (E_Int j = 0; j < rpcount[i];) {
      E_Int lf = ptr[j++];
      E_Int rf = ptr[j++];

      // rank and i communicate through duplicate faces lf/rf
      cpatch.add(lf, rf);

      dups.insert(lf);

      ndup++;
    }

    if (rpcount[i] > 0)
      cpatches.push_back(cpatch);
  }

  // sort cpatches faces following master rank
  for (auto& p : cpatches) {
    std::vector<E_Int> indices(p.local.size());

    for (size_t i = 0; i < indices.size(); i++)
      indices[i] = i;

    if (rank < p.nei) {
      // local rules
      std::sort(indices.begin(), indices.end(),
      [&](E_Int a, E_Int b)
      {
        return p.local[a] < p.local[b];
      });

    } else {
      // remote rules
      std::sort(indices.begin(), indices.end(),
      [&](E_Int a, E_Int b)
      {
        return p.remote[a] < p.remote[b];
      });
    }

    std::vector<E_Int> tmplocal(p.local);
    std::vector<E_Int> tmpremote(p.remote);

    for (size_t i = 0; i < indices.size(); i++) {
      p.local[i] = tmplocal[indices[i]];
      p.remote[i] = tmpremote[indices[i]];
    }
  }

  // free faces
  E_Int nfree = nfaces - ndup;
  E_Int gstart;
  MPI_Scan(&nfree, &gstart, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  gstart -= nfree;

  std::unordered_map<E_Int, E_Int> old2new;

  E_Int nf = gstart+1;

  for (auto &face : NFACE) {
    if (dups.find(face) != dups.end()) {
      face = -face;
    } else {
      auto search = old2new.find(face);
      if (search == old2new.end()) {
        old2new[face] = nf;
        face = nf++;
      } else{
        face = search->second;
      }
    }
  }

  E_Int npatches = cpatches.size();

  // how many faces do i control?
  E_Int ncontrol = 0;
  for (E_Int i = 0; i < npatches; i++) {
    const auto &patch = cpatches[i];
    E_Int nei = patch.nei;
    if (rank > nei) continue;

    ncontrol += patch.local.size();
  }

  std::vector<E_Int> control_dist(nproc+1);
  control_dist[0] = 0;
  MPI_Allgather(&ncontrol, 1, MPI_INT, &control_dist[0]+1, 1, MPI_INT, MPI_COMM_WORLD);
  for (E_Int i = 0; i < nproc; i++)
    control_dist[i+1] += control_dist[i];

  E_Int gnfree;
  MPI_Allreduce(&nfree, &gnfree, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  E_Int my_start = gnfree + control_dist[rank] + 1;

  std::vector<MPI_Request> req(npatches);
  E_Int nreq = 0;
  
  for (E_Int i = 0; i < npatches; i++) {
    auto& patch = cpatches[i];
    E_Int nei = patch.nei;
    auto& local = patch.local;
    auto& remote = patch.remote;

    if (rank < nei) {
      // renumber faces i control
      for (auto& face : local) {
        old2new[-face] = my_start;
        face = my_start++;
      }

      // update neighbour
      MPI_Isend(&local[0], local.size(), MPI_INT, nei, 0, MPI_COMM_WORLD, &req[nreq++]); 
    } else {
      MPI_Irecv(&remote[0], remote.size(), MPI_INT, nei, 0, MPI_COMM_WORLD, &req[nreq++]);
    }

    assert(nreq <= npatches);
  }

  MPI_Waitall(nreq, &req[0], MPI_STATUSES_IGNORE);
  nreq = 0;

  for (E_Int i = 0; i < npatches; i++) {
    auto& patch = cpatches[i];
    E_Int nei = patch.nei;
    auto& local = patch.local;
    auto& remote = patch.remote;

    if (rank < nei) continue;
    
    // update
    for (size_t i = 0; i < local.size(); i++) {
      E_Int rf = remote[i];
      E_Int lf = local[i];
      old2new[-lf] = rf;
    }
  }

  for (auto &face : NFACE) {
    if (face < 0)
      face = old2new[face];
  }

  //EXIT;

  // make face distribution
  std::vector<E_Int> fdist(nproc+1);
  fdist[0] = 0;
  E_Int gnfaces;
  E_Int nlf = ncontrol + nfree;
  MPI_Allreduce(&nlf, &gnfaces, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0)
    printf("total number of unique faces: %d\n", gnfaces);
  E_Int k = gnfaces;
  for (E_Int i = 0; i < nproc; i++) {
    E_Int l = k / (nproc-i);
    fdist[i+1] = fdist[i] + l;
    k -= l;
  }
  nf = fdist[rank+1] - fdist[rank]; // how many and which faces i will get

  // update lPE
  // Note(Imad): for now we create a new map
  // WARNING(Imad): check that this is ok!!!
  std::unordered_map<E_Int, std::vector<E_Int>> PE;
  for (auto i = 0; i < ncells; i++) {
    E_Int *pf = &NFACE[fstride*i];
    for (E_Int j = 0; j < fstride; j++)
      PE[pf[j]].push_back(i + firstCell);
  }

  // we know PE of every face
  std::vector<E_Int> f_scount(nproc, 0);
  for (auto& face : PE) {
    E_Int gf = face.first;
    assert(gf > 0);
    E_Int target = get_proc(gf-1, fdist, nproc);
    f_scount[target] += 1 + 1 + face.second.size(); // id + size + own/nei
  }

  std::vector<E_Int> f_rcount(nproc, 0);
  MPI_Alltoall(&f_scount[0], 1, MPI_INT, &f_rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> f_sdist(nproc+1), f_rdist(nproc+1);
  f_sdist[0] = f_rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    f_sdist[i+1] = f_sdist[i] + f_scount[i];
    f_rdist[i+1] = f_rdist[i] + f_rcount[i];
  }

  std::vector<E_Int> idx(nproc), f_sdata(f_sdist[nproc]), f_rdata(f_rdist[nproc]);
  for (E_Int i = 0; i < nproc; i++)
    idx[i] = f_sdist[i];

  for (const auto& face : PE) {
    E_Int target = get_proc(face.first-1, fdist, nproc);
    f_sdata[idx[target]++] = face.first;
    f_sdata[idx[target]++] = face.second.size();
    for (const auto& elem : face.second)
      f_sdata[idx[target]++] = elem;
  }

  MPI_Alltoallv(&f_sdata[0], &f_scount[0], &f_sdist[0], MPI_INT,
                &f_rdata[0], &f_rcount[0], &f_rdist[0], MPI_INT,
                MPI_COMM_WORLD);


  std::vector<E_Int> A(2*nf, -1);
  std::vector<E_Int> c(nf, 0);

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &f_rdata[f_rdist[i]];
    for (E_Int j = 0; j < f_rcount[i]; ) {
      assert(rank == get_proc(ptr[j]-1, fdist, nproc));
      E_Int face = ptr[j++] - 1 - fdist[rank];
      E_Int size = ptr[j++];
      for (E_Int k = 0; k < size; k++) {
        A[2*face + c[face]++] = ptr[j++];
        assert(c[face] <= 2);
      }
    }
  }

  if (rank == 0)
    puts("ParentElements OK");

  // dual graph
  std::unordered_map<E_Int, std::vector<E_Int>> CADJ;
  for (E_Int i = 0; i < nf; i++) {
    E_Int nei = A[2*i+1];
    if (nei >= 0) {
      E_Int own = A[2*i];
      CADJ[own].push_back(nei);
      CADJ[nei].push_back(own);
    }
  }

  std::vector<E_Int> scount(nproc, 0), rcount(nproc);
  
  for (const auto &cell : CADJ) {
    E_Int target = get_proc(cell.first, cdist, nproc);
    scount[target] += 1 + 1 + cell.second.size(); // id + size + neis
  }

  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> sdist(nproc+1), rdist(nproc+1);
  sdist[0] = rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  std::vector<E_Int> sdata(sdist[nproc]);
  std::vector<E_Int> rdata(rdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];

  for (const auto& cell : CADJ) {
    E_Int target = get_proc(cell.first, cdist, nproc);
    sdata[idx[target]++] = cell.first;
    sdata[idx[target]++] = cell.second.size();
    for (const auto& elem : cell.second)
      sdata[idx[target]++] = elem;
  }

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], MPI_INT,
                &rdata[0], &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  std::vector<std::vector<E_Int>> cadj(ncells);
  E_Int nedges = 0;
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &rdata[rdist[i]];
    for (E_Int j = 0; j < rcount[i];) {
      assert(get_proc(ptr[j], cdist, nproc) == rank);
      E_Int cell = ptr[j++] - firstCell;
      E_Int size = ptr[j++];
      nedges += size;
      for (E_Int k = 0; k < size; k++) {
        cadj[cell].push_back(ptr[j++]);
      }
    }
  }

  std::vector<E_Int> ADJ;
  ADJ.reserve(nedges);
  std::vector<E_Int> xadj(ncells+1);
  xadj[0] = 0;
  E_Int count = 0;
  for (const auto& adj : cadj) {
    ADJ.insert(ADJ.end(), adj.begin(), adj.end());
    xadj[count+1] = xadj[count] + adj.size();
    count++;
  }
  assert(count == ncells);
  assert(xadj[count] == nedges);

  if (rank == 0)
    printf("Dual graph OK\n");

  SCOTCH_Dgraph graph;
  SCOTCH_dgraphInit(&graph, MPI_COMM_WORLD);

  E_Int ret;
  ret = SCOTCH_dgraphBuild(
    &graph,       // grafptr
    0,            // baseval
    ncells,       // vertlocnbr
    ncells,       // vertlocmax
    &xadj[0],     // vertloctab
    NULL,         // vendloctab
    NULL,         // veloloctab
    NULL,         // vlblocltab
    nedges,       // edgelocnbr
    nedges,       // edgelocsiz
    &ADJ[0],      // edgeloctab
    NULL,         // edgegsttab
    NULL);        // edloloctab

  if (ret != 0) {
    fprintf(stderr, "SCOTCH_dgraphBuild(): Failed to load graph\n");
  }

  ret = SCOTCH_dgraphCheck(&graph);
  if (ret != 0) {
    fprintf(stderr, "SCOTCH_dgraphCheck(): Failed graph check\n");
  }

  SCOTCH_Strat strat;
  ret = SCOTCH_stratInit(&strat);
  if (ret != 0) {
    fprintf(stderr, "SCOTCH_startInit(): Failed to init strat\n");
  }

  std::vector<E_Int> part(ncells);
  ret = SCOTCH_dgraphPart(&graph, nproc, &strat, &part[0]); 

  if (ret != 0) {
    fprintf(stderr, "SCOTCH_dgraphPart(): Failed to map graph\n");
  }

  if (rank == 0)
    printf("Graph map OK\n");


  // distribute
  std::vector<E_Int> c_scount(nproc, 0), c_rcount(nproc);
  
  for (E_Int i = 0; i < ncells; i++) {
    E_Int target = part[i];
    c_scount[target] += 1;
  }

  MPI_Alltoall(&c_scount[0], 1, MPI_INT, &c_rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> c_sdist(nproc+1, 0);
  std::vector<E_Int> c_rdist(nproc+1, 0);
  for (E_Int i = 0; i < nproc; i++) {
    c_sdist[i+1] = c_sdist[i] + c_scount[i];
    c_rdist[i+1] = c_rdist[i] + c_rcount[i];
  }

  E_Int nncells = c_rdist[nproc];
  std::vector<E_Int> scells(c_sdist[nproc]);
  std::vector<E_Int> rcells(c_rdist[nproc]);

  for (E_Int i = 0; i < nproc; i++) idx[i] = c_sdist[i];

  for (E_Int i = 0; i < ncells; i++) {
    E_Int target = part[i];
    assert(target >= 0 && target < nproc);
    scells[idx[target]] = i + cdist[rank];
    idx[target]++;
  }

  MPI_Alltoallv(&scells[0], &c_scount[0], &c_sdist[0], MPI_INT,
                &rcells[0], &c_rcount[0], &c_rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  // send connectivity
  for (E_Int i = 0; i < nproc; i++) {
    scount[i] = stride*c_scount[i];
    rcount[i] = stride*c_rcount[i];
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];
  
  sdata.resize(sdist[nproc]);
  std::vector<E_Int> ncn(rdist[nproc]);
  assert(rdist[nproc] == nncells*stride);
  
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pc = &scells[c_sdist[i]];
    E_Int *ptr = &sdata[sdist[i]];
    for (E_Int j = 0; j < c_scount[i]; j++) {
      E_Int cell = pc[j] - firstCell;
      E_Int *pn = &cn[stride*cell];
      for (E_Int k = 0; k < stride; k++)
        *ptr++ = pn[k];
    }
  }

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], MPI_INT,
                &ncn[0], &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);
  
  // hash cells
  std::unordered_map<E_Int, E_Int> CT;
  for (E_Int i = 0; i < nncells; i++)
    CT[rcells[i]] = i;

  // request points
  E_Int nnpoints = 0;
  std::unordered_map<E_Int, E_Int> PT;
  std::vector<E_Int> p_rcount(nproc, 0), p_scount(nproc);

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pc = &rcells[c_rdist[i]];
    for (E_Int j = 0; j < c_rcount[i]; j++) {
      E_Int lc = CT[pc[j]];
      E_Int *pn = &ncn[stride*lc];
      for (E_Int k = 0; k < stride; k++) {
        E_Int point = pn[k];
        if (PT.find(point) == PT.end()) {
          PT[point] = nnpoints++;
          E_Int src = get_proc(point-1, pdist, nproc);
          p_rcount[src]++;
        }
      }
    }
  }
  
  MPI_Alltoall(&p_rcount[0], 1, MPI_INT, &p_scount[0], 1, MPI_INT, MPI_COMM_WORLD);
  
  std::vector<E_Int> p_rdist(nproc+1);
  std::vector<E_Int> p_sdist(nproc+1);
  p_rdist[0] = p_sdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    p_rdist[i+1] = p_rdist[i] + p_rcount[i];
    p_sdist[i+1] = p_sdist[i] + p_scount[i];
  }

  assert(nnpoints == p_rdist[nproc]);

  std::vector<E_Int> rpoints(nnpoints);
  std::vector<E_Int> spoints(p_sdist[nproc]);

  for (E_Int i = 0; i < nproc; i++) idx[i] = p_rdist[i];

  std::vector<E_Int> vpoints(nnpoints, 0);
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pc = &rcells[c_rdist[i]];
    for (E_Int j = 0; j < c_rcount[i]; j++) {
      E_Int lc = CT[pc[j]];
      E_Int *pn = &ncn[stride*lc];
      for (E_Int k = 0; k < stride; k++) {
        E_Int point = pn[k];
        if (!vpoints[PT[point]]) {
          vpoints[PT[point]]++;
          E_Int src = get_proc(point-1, pdist, nproc);
          rpoints[idx[src]++] = point;
        }
      }
    }
  }

  MPI_Alltoallv(&rpoints[0], &p_rcount[0], &p_rdist[0], MPI_INT,
                &spoints[0], &p_scount[0], &p_sdist[0], MPI_INT,
                MPI_COMM_WORLD);

  PT.clear();
  nnpoints = 0; 
  for (const auto& point : rpoints)
    PT[point] = nnpoints++;

  // send coordinates
  for (E_Int i = 0; i < nproc; i++) {
    scount[i] = 3*p_scount[i];
    rcount[i] = 3*p_rcount[i];
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }
  
  std::vector<E_Float> sxyz(sdist[nproc]);
  std::vector<E_Float> rxyz(rdist[nproc]);

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pp = &spoints[p_sdist[i]];
    E_Float *px = &sxyz[sdist[i]];
    for (E_Int j = 0; j < p_scount[i]; j++) {
      E_Int point = pp[j] - 1 - pdist[rank];
      *px++ = X[point];
      *px++ = Y[point];
      *px++ = Z[point];
    }
  }

  MPI_Alltoallv(&sxyz[0], &scount[0], &sdist[0], MPI_DOUBLE,
                &rxyz[0], &rcount[0], &rdist[0], MPI_DOUBLE,
                MPI_COMM_WORLD);

  if (rank == 0)
    puts("Coordinates OK");

  // Build output
  const char *varString = "CoordinateX,CoordinateY,CoordinateZ";

  PyObject* m = K_ARRAY::buildArray3(3, varString, nnpoints, nncells, eltName, false, 3);

  FldArrayF *fo;
  FldArrayI *cno;
  K_ARRAY::getFromArray3(m, fo, cno);

  for (E_Int n = 0; n < 3; n++) { 
    E_Float *px = fo->begin(n+1);
    for (E_Int i = 0; i < nnpoints; i++)
      px[i] = rxyz[3*i+n];
  }

  E_Int *pn = cno->begin();
  for (E_Int i = 0; i < stride*nncells; i++)
    pn[i] = PT[ncn[i]]+1;

  return m;
}
