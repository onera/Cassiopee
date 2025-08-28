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
  printf("\nWarning: could not find distribution of element " SF_D_ "\n", element);
  assert(0);
  return -1;
}

static
E_Int get_chunk(E_Int element, const std::vector<E_Int> &dist)
{
  for (size_t i = 0; i < dist.size(); i++) {
    if (element >= dist[i] && element < dist[i+1])
      return i;
  }
  printf("\nWarning: could not find chunk of element " SF_D_ "\n", element);
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
    for (E_Int i = 0; i < stride; i++)
      n[i] = p[i];
    std::sort(n, n+stride);
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
  NFACE.clear();

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
}

static
void make_facets_TETRA(E_Int *cn, E_Int ncells, E_Int &nfaces,
  std::vector<E_Int> &NFACE, std::unordered_map<Poly, E_Int, Poly_hash> &faces,
  std::vector<Poly> &f2Q)
{
  NFACE.clear();

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
}

static
void make_facets_PENTA(E_Int *cn, E_Int ncells, E_Int &nfaces,
  std::vector<E_Int> &NFACE, std::unordered_map<Poly, E_Int, Poly_hash> &faces,
  std::vector<Poly> &f2Q)
{
  NFACE.clear();

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
}

static
void make_facets_PYRA(E_Int *cn, E_Int ncells, E_Int &nfaces,
  std::vector<E_Int> &NFACE, std::unordered_map<Poly, E_Int, Poly_hash> &faces,
  std::vector<Poly> &f2Q)
{
  NFACE.clear();

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
}

static
void make_facets_QUAD(E_Int *cn, E_Int ncells, E_Int &nfaces,
  std::vector<E_Int> &NFACE, std::unordered_map<Poly, E_Int, Poly_hash> &faces,
  std::vector<Poly> &f2Q)
{
  NFACE.clear();

  Poly Edge(2);
  E_Int n[4];

  for (E_Int i = 0; i < ncells; i++) {
    E_Int *pn = &cn[4*i];
    
    // bot
    n[0] = pn[0];
    n[1] = pn[1];
    Edge.set(n);
    add_facet(Edge, nfaces, faces, NFACE, f2Q);

    // right
    n[0] = pn[1];
    n[1] = pn[2];
    Edge.set(n);
    add_facet(Edge, nfaces, faces, NFACE, f2Q);

    // top
    n[0] = pn[2];
    n[1] = pn[3];
    Edge.set(n);
    add_facet(Edge, nfaces, faces, NFACE, f2Q);

    // left
    n[0] = pn[3];
    n[1] = pn[0];
    Edge.set(n);
    add_facet(Edge, nfaces, faces, NFACE, f2Q);
  }

  assert((E_Int)NFACE.size() == 4*ncells);
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
  } else if (stride == 4 && strcmp(elemName, "QUAD") == 0) {
    make_facets_QUAD(cn, ncells, nfaces, NFACE, faces, f2Q);
  } else if (stride == 6) {
    make_facets_PENTA(cn, ncells, nfaces, NFACE, faces, f2Q);
  } else if (stride == 8) {
    make_facets_HEXA(cn, ncells, nfaces, NFACE, faces, f2Q);
  } else {
    fprintf(stderr, "chunk2part_elt(): %s not supported.\n", elemName);
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

PyObject *K_XCORE::chunk2partElt(PyObject *self, PyObject *args)
{
  int rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  PyObject *XYZ, *CHUNKS;

  if (!PyArg_ParseTuple(args, "OO", &XYZ, &CHUNKS)) {
    PyErr_SetString(PyExc_ValueError, "chunk2partElt(): bad input");
    return NULL;
  }

  // coordinates
  PyObject *Xi, *Yi, *Zi;
  E_Float *X, *Y, *Z;
  E_Int npoints, nfld;

  Xi = PyList_GetItem(XYZ, 0);
  K_NUMPY::getFromNumpyArray(Xi, X, npoints, nfld);

  Yi = PyList_GetItem(XYZ, 1);
  K_NUMPY::getFromNumpyArray(Yi, Y, npoints, nfld);

  Zi = PyList_GetItem(XYZ, 2);
  K_NUMPY::getFromNumpyArray(Zi, Z, npoints, nfld);

  // chunks
  E_Int nchunks = PyList_Size(CHUNKS);

  std::vector<const char *> eltNames(nchunks);
  std::vector<E_Int> strides(nchunks);
  std::vector<E_Int *> cns(nchunks);
  std::vector<E_Int> local_cdist(nchunks+1, 0);
  std::vector<E_Int> csize(nchunks);
  std::vector<E_Int> fstrides(nchunks);


  E_Int ncells = 0;

  for (E_Int i = 0; i < nchunks; i++) {
    PyObject *o = PyList_GetItem(CHUNKS, i);

    // [name, stride, connectivity]
    PyObject *NAME = PyList_GetItem(o, 0);
#if PY_VERSION_HEX >= 0x03000000
    eltNames[i] = PyUnicode_AsUTF8(NAME);
#else
    eltNames[i] = PyString_AsString(NAME);
#endif

    PyObject *STRIDE = PyList_GetItem(o, 1);
    strides[i] = PyLong_AsLong(STRIDE);

    PyObject *CN = PyList_GetItem(o, 2);
    E_Int cnsize, nfld;
    K_NUMPY::getFromNumpyArray(CN, cns[i], cnsize, nfld);
    csize[i] = cnsize / strides[i];
    local_cdist[i+1] = csize[i];
    ncells += local_cdist[i+1];

    if (strcmp(eltNames[i], "HEXA") == 0)
      fstrides[i] = 6;
    else if (strcmp(eltNames[i], "PENTA") == 0)
      fstrides[i] = 5;
    else if (strcmp(eltNames[i], "TETRA") == 0)
      fstrides[i] = 4;
    else if (strcmp(eltNames[i], "PYRA") == 0)
      fstrides[i] = 5;
    //TODO(Imad): QUAD et TRI
    /*
    else if (strcmp(eltNames[i], "TRI") == 0)
      fstrides[i] = 1;
    else if (strcmp(eltNames[i], "QUAD") == 0)
      fstrides[i] = 1;
    */
  }

  for (E_Int i = 0; i < nchunks; i++)
    local_cdist[i+1] += local_cdist[i];
  
  assert(ncells == local_cdist[nchunks]);

  // build cell and points distributions
  std::vector<E_Int> cdist((nproc+1));
  std::vector<E_Int> pdist((nproc+1));
  cdist[0] = pdist[0] = 0;

  MPI_Allgather(&npoints, 1, XMPI_INT, &pdist[0]+1, 1, XMPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(&ncells, 1, XMPI_INT, &cdist[0]+1, 1, XMPI_INT, MPI_COMM_WORLD);
  for (E_Int i = 0; i < nproc; i++) {
    cdist[i+1] += cdist[i];
    pdist[i+1] += pdist[i];
  }

  E_Int gncells = 0;
  MPI_Allreduce(&ncells, &gncells, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);
  assert(gncells == cdist[nproc]);
  if (rank == 0)
    printf("Total number of cells: " SF_D_ " (Average: " SF_D_ ")\n", gncells, gncells/nproc);

  // make local cell labels
  std::vector<E_Int> cellLabels(nchunks);

  for (E_Int i = 0; i < nchunks; i++)
    cellLabels[i] = csize[i];
  
  for (E_Int i = 0; i < nchunks-1; i++)
    cellLabels[i+1] += cellLabels[i];
  
  for (E_Int i = 0; i < nchunks; i++) {
    cellLabels[i] -= csize[i];
    cellLabels[i] += cdist[rank];
  }

  // make cellFaces connectivity
  E_Int nfaces = 1;
  std::vector<std::vector<E_Int>> NFACES(nchunks);
  std::unordered_map<Poly, E_Int, Poly_hash> faces;
  std::vector<Poly> face2Poly;
  for (E_Int i = 0; i < nchunks; i++) {
    E_Int res = make_facets(cns[i], csize[i], strides[i], eltNames[i], nfaces,
      NFACES[i], faces, face2Poly);
    if (res == 1)
      return NULL;
  }

  nfaces--;
  E_Int firstFace;
  MPI_Scan(&nfaces, &firstFace, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);

  firstFace -= nfaces;

  // adjust global index of faces
  for (auto &face : faces)
    face.second = face.second + firstFace;
  
  // adjust NFACE
  for (auto& NFACE : NFACES)
    for (auto &face : NFACE)
      face = face + firstFace;
  
  // make local ParentElement

  std::unordered_map<E_Int, std::vector<E_Int>> lPE;

  for (E_Int chunk = 0; chunk < nchunks; chunk++) {
    auto& NFACE = NFACES[chunk];
    E_Int nc = csize[chunk];
    E_Int fstr = fstrides[chunk];
    E_Int cstart = cellLabels[chunk];

    for (E_Int i = 0; i < nc; i++) {
      E_Int *pf = &NFACE[fstr*i];
      for (E_Int j = 0; j < fstr; j++) {
        lPE[pf[j]].push_back(i+cstart);
      }
    }
  }

  // potential pface data to exchange with every proc
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

  int ssize = SEND.size();
  std::vector<int> RCOUNT(nproc);
  MPI_Allgather(&ssize, 1, MPI_INT, &RCOUNT[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> RDIST(nproc+1);
  RDIST[0] = 0;
  for (E_Int i = 0; i < nproc; i++)
    RDIST[i+1] = RDIST[i] + RCOUNT[i];
  
  std::vector<E_Int> RECV(RDIST[nproc], -1);

  MPI_Allgatherv(&SEND[0], ssize, XMPI_INT,
                 &RECV[0], &RCOUNT[0], &RDIST[0], XMPI_INT,
                 MPI_COMM_WORLD);


  // look for received faces in face table
  std::vector<E_Int> spfaces;
  std::vector<int> spcount(nproc, 0), rpcount(nproc);

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

  std::vector<int> spdist(nproc+1), rpdist(nproc+1);
  spdist[0] = rpdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    spdist[i+1] = spdist[i] + spcount[i];
    rpdist[i+1] = rpdist[i] + rpcount[i];
  }

  std::vector<E_Int> rpfaces(rpdist[nproc]);
  MPI_Alltoallv(&spfaces[0], &spcount[0], &spdist[0], XMPI_INT,
                &rpfaces[0], &rpcount[0], &rpdist[0], XMPI_INT,
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
  MPI_Scan(&nfree, &gstart, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);
  gstart -= nfree;

  std::unordered_map<E_Int, E_Int> old2new;

  E_Int nf = gstart+1;

  for (auto &NFACE : NFACES) {
    for (auto &face : NFACE) {
      if (dups.find(face) != dups.end()) {
        face = -face;
      } else {
        auto search = old2new.find(face);
        if (search == old2new.end()) {
          old2new[face] = nf;
          face = nf++;
        } else {
          face = search->second;
        }
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
  MPI_Allgather(&ncontrol, 1, XMPI_INT, &control_dist[0]+1, 1, XMPI_INT, MPI_COMM_WORLD);
  for (E_Int i = 0; i < nproc; i++)
    control_dist[i+1] += control_dist[i];

  E_Int gnfree;
  MPI_Allreduce(&nfree, &gnfree, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);
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
      MPI_Isend(&local[0], local.size(), XMPI_INT, nei, 0, MPI_COMM_WORLD, &req[nreq++]); 
    } else {
      MPI_Irecv(&remote[0], remote.size(), XMPI_INT, nei, 0, MPI_COMM_WORLD, &req[nreq++]);
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

  for (auto& NFACE : NFACES) {
    for (auto &face : NFACE) {
      if (face < 0)
        face = old2new[face];
    }
  }

  // make face distribution
  std::vector<E_Int> fdist(nproc+1);
  fdist[0] = 0;
  E_Int gnfaces;
  E_Int nlf = ncontrol + nfree;
  MPI_Allreduce(&nlf, &gnfaces, 1, XMPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (rank == 0)
    printf("total number of unique faces: " SF_D_ "\n", gnfaces);
  E_Int k = gnfaces;
  for (E_Int i = 0; i < nproc; i++) {
    E_Int l = k / (nproc-i);
    fdist[i+1] = fdist[i] + l;
    k -= l;
  }
  nf = fdist[rank+1] - fdist[rank]; // how many and which faces i will get

  // update lPE
  // Note(Imad): for now we create a new map
  std::unordered_map<E_Int, std::vector<E_Int>> PE;
  for (E_Int chunk = 0; chunk < nchunks; chunk++) {
    auto& NFACE = NFACES[chunk];
    E_Int nc = csize[chunk];
    E_Int fstr = fstrides[chunk];
    E_Int cstart = cellLabels[chunk];

    for (E_Int i = 0; i < nc; i++) {
      E_Int *pf = &NFACE[fstr*i];
      for (E_Int j = 0; j < fstr; j++) {
        PE[pf[j]].push_back(i+cstart);
      }
    }
  }

  // we know PE of every face
  std::vector<int> f_scount(nproc, 0);
  for (auto& face : PE) {
    E_Int gf = face.first;
    assert(gf > 0);
    E_Int target = get_proc(gf-1, fdist, nproc);
    f_scount[target] += 1 + 1 + face.second.size(); // id + size + own/nei
  }

  std::vector<int> f_rcount(nproc, 0);
  MPI_Alltoall(&f_scount[0], 1, MPI_INT, &f_rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> f_sdist(nproc+1), f_rdist(nproc+1);
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

  MPI_Alltoallv(&f_sdata[0], &f_scount[0], &f_sdist[0], XMPI_INT,
                &f_rdata[0], &f_rcount[0], &f_rdist[0], XMPI_INT,
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

  std::vector<int> scount(nproc, 0), rcount(nproc);
  
  for (const auto &cell : CADJ) {
    E_Int target = get_proc(cell.first, cdist, nproc);
    scount[target] += 1 + 1 + cell.second.size(); // id + size + neis
  }

  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> sdist(nproc+1), rdist(nproc+1);
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

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], XMPI_INT,
                &rdata[0], &rcount[0], &rdist[0], XMPI_INT,
                MPI_COMM_WORLD);

  std::vector<std::vector<E_Int>> cadj(ncells);
  E_Int nedges = 0;
  E_Int firstCell = cdist[rank];
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

  SCOTCH_Arch arch;
  ret = SCOTCH_archInit(&arch);
  if (ret != 0) fprintf(stderr, "Bad SCOTCH_archInit");

  ret = SCOTCH_archCmplt(&arch, nproc);
  if (ret != 0) fprintf(stderr, "Bad archCmplt");

  std::vector<E_Int> part(ncells);
  //ret = SCOTCH_dgraphPart(&graph, nproc, &strat, &part[0]); 
  ret = SCOTCH_dgraphMap(&graph, &arch, &strat, &part[0]);

  if (ret != 0) {
    fprintf(stderr, "SCOTCH_dgraphPart(): Failed to map graph\n");
  }

  SCOTCH_dgraphExit(&graph);
  SCOTCH_stratExit(&strat);
  SCOTCH_archExit(&arch);

  if (rank == 0)
    printf("Graph map OK\n");

  // distribute
  std::vector<int> c_scount(nproc, 0), c_rcount(nproc);
  
  for (E_Int i = 0; i < ncells; i++) {
    E_Int target = part[i];
    c_scount[target] += 1;
  }

  MPI_Alltoall(&c_scount[0], 1, MPI_INT, &c_rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> c_sdist(nproc+1, 0);
  std::vector<int> c_rdist(nproc+1, 0);
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

  MPI_Alltoallv(&scells[0], &c_scount[0], &c_sdist[0], XMPI_INT,
                &rcells[0], &c_rcount[0], &c_rdist[0], XMPI_INT,
                MPI_COMM_WORLD);


  // send connectivity
  for (E_Int i = 0; i < nproc; i++)
    scount[i] = 0;

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &scells[c_sdist[i]];
    for (E_Int j = 0; j < c_scount[i]; j++) {
      E_Int cell = ptr[j] - cdist[rank];
      E_Int chunk = get_chunk(cell, local_cdist);
      scount[i] += strides[chunk];
    }
  }

  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  sdist[0] = rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  sdata.resize(sdist[nproc]);
  rdata.resize(rdist[nproc]);

  std::vector<E_Int> sstride(c_sdist[nproc]);
  std::vector<E_Int> rstride(c_rdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &scells[c_sdist[i]];
    E_Int *sp = &sstride[c_sdist[i]];

    for (E_Int j = 0; j < c_scount[i]; j++) {
      E_Int cell = ptr[j];

      E_Int lc = cell - cdist[rank];

      E_Int chunk = get_chunk(lc, local_cdist);

      E_Int stride = strides[chunk];
      E_Int lsize = local_cdist[chunk];
      const auto& cn = cns[chunk];

      *sp++ = stride;
      lc -= lsize;

      E_Int *pn = &cn[lc * stride];
      for (E_Int k = 0; k < stride; k++)
        sdata[idx[i]++] = pn[k];
    }
  }


  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], XMPI_INT,
                &rdata[0], &rcount[0], &rdist[0], XMPI_INT,
                MPI_COMM_WORLD);
  
  MPI_Alltoallv(&sstride[0], &c_scount[0], &c_sdist[0], XMPI_INT,
                &rstride[0], &c_rcount[0], &c_rdist[0], XMPI_INT,
                MPI_COMM_WORLD);

  if (rank == 0)
    puts("Connectivity OK");

  // hash cells
  std::unordered_map<E_Int, E_Int> CT;
  for (E_Int i = 0; i < nncells; i++)
    CT[rcells[i]] = i;
  
  // request points
  E_Int nnpoints = 0;
  std::unordered_map<E_Int, E_Int> PT;
  std::vector<int> p_rcount(nproc, 0), p_scount(nproc);

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ps = &rstride[c_rdist[i]];
    E_Int *pr = &rdata[rdist[i]];
    for (E_Int j = 0; j < c_rcount[i]; j++) {
      E_Int stride = ps[j];
      for (E_Int k = 0; k < stride; k++) {
        E_Int point = *pr++;
        if (PT.find(point) == PT.end()) {
          PT[point] = nnpoints++;
          E_Int src = get_proc(point-1, pdist, nproc);
          p_rcount[src]++;
        }
      }
    }
  }

  MPI_Alltoall(&p_rcount[0], 1, MPI_INT, &p_scount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> p_rdist(nproc+1);
  std::vector<int> p_sdist(nproc+1);
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
    E_Int *ps = &rstride[c_rdist[i]];
    E_Int *pr = &rdata[rdist[i]];
    for (E_Int j = 0; j < c_rcount[i]; j++) {
      E_Int stride = ps[j];
      for (E_Int k = 0; k < stride; k++) {
        E_Int point = *pr++;
        if (!vpoints[PT[point]]) {
          vpoints[PT[point]]++;
          E_Int src = get_proc(point-1, pdist, nproc);
          rpoints[idx[src]++] = point;
        }
      }
    }
  }

  MPI_Alltoallv(&rpoints[0], &p_rcount[0], &p_rdist[0], XMPI_INT,
                &spoints[0], &p_scount[0], &p_sdist[0], XMPI_INT,
                MPI_COMM_WORLD);

  PT.clear();
  nnpoints = 0; 
  for (const auto& point : rpoints)
    PT[point] = nnpoints++;

  // send coordinates
  std::vector<int> sxcount(nproc);
  std::vector<int> rxcount(nproc);
  std::vector<int> sxdist(nproc+1);
  std::vector<int> rxdist(nproc+1);
  sxdist[0] = rxdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    sxcount[i] = 3*p_scount[i];
    rxcount[i] = 3*p_rcount[i];
    sxdist[i+1] = sxdist[i] + sxcount[i];
    rxdist[i+1] = rxdist[i] + rxcount[i];
  }

  std::vector<E_Float> sxyz(sxdist[nproc]);
  std::vector<E_Float> rxyz(rxdist[nproc]);

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pp = &spoints[p_sdist[i]];
    E_Float *px = &sxyz[sxdist[i]];
    for (E_Int j = 0; j < p_scount[i]; j++) {
      E_Int point = pp[j] - 1 - pdist[rank];
      *px++ = X[point];
      *px++ = Y[point];
      *px++ = Z[point];
    }
  }

  MPI_Alltoallv(&sxyz[0], &sxcount[0], &sxdist[0], MPI_DOUBLE,
                &rxyz[0], &rxcount[0], &rxdist[0], MPI_DOUBLE,
                MPI_COMM_WORLD);

  if (rank == 0)
    puts("Coordinates OK");
  
  // build output
  if (rank == 0)
    puts("Exporting...");

  PyObject *out = PyList_New(0);

  const char *varString = "CoordinateX,CoordinateY,CoordinateZ";

  // count number of cells and points per element type
  std::unordered_map<E_Int, E_Int> stride2nelem;
  std::unordered_map<E_Int, std::set<E_Int>> stride2points;

  for (E_Int i = 0; i < nchunks; i++)
    stride2nelem[strides[i]] = 0;

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ps = &rstride[c_rdist[i]];
    E_Int *pr = &rdata[rdist[i]];

    for (E_Int j = 0; j < c_rcount[i]; j++) {
      E_Int stride = ps[j];

      stride2nelem[stride]++;


      for (E_Int k = 0; k < stride; k++) {
        E_Int point = *pr++;
        stride2points[stride].insert(point);
      }
    }
  }

  std::unordered_map<E_Int, PyObject *> stride2arr;
  std::unordered_map<E_Int, FldArrayF *> stride2fo;
  std::unordered_map<E_Int, FldArrayI *> stride2cno;
  std::unordered_map<E_Int, E_Int *> begin_ptr;
  std::unordered_map<E_Int, E_Int> stride2count;
  std::unordered_map<E_Int, std::unordered_map<E_Int, E_Int>> stride2localp;

  for (auto& zone : stride2points) {
    E_Int stride = zone.first;
    const auto& points = zone.second;
    E_Int ncells = stride2nelem[stride];
    E_Int npoints = points.size();

    char eltName[10];

    // TODO(Imad): what if stride == 4 but elt is QUAD?
    if (stride == 4)
      strcpy(eltName, "TETRA");
    else if (stride == 5)
      strcpy(eltName, "PYRA");
    else if (stride == 6)
      strcpy(eltName, "PENTA");
    else if (stride == 8)
      strcpy(eltName, "HEXA");

    stride2arr[stride] = K_ARRAY::buildArray3(3, varString, npoints, ncells, eltName, false, 3);
    K_ARRAY::getFromArray3(stride2arr[stride], stride2fo[stride], stride2cno[stride]);
  
    begin_ptr[stride] = stride2cno[stride]->begin();

    stride2count[stride] = 0;

    // make local points table
    E_Int np = 0;
    auto& localp = stride2localp[stride];
    for (const auto &point : points)
      localp[point] = np++;
    
    // fill in xyz
    for (E_Int n = 0; n < 3; n++) {
      E_Float *px = stride2fo[stride]->begin(n+1);
      for (auto& point : points) {
        E_Int pt = localp[point];
        px[pt] = rxyz[3*PT[point]+n];
      }
    }
  }

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ps = &rstride[c_rdist[i]];
    E_Int *pn = &rdata[rdist[i]];

    for (E_Int j = 0; j < c_rcount[i]; j++) {
      E_Int stride = ps[j];
      E_Int& where = stride2count[stride];
      E_Int *ptr = begin_ptr[stride];
      auto& localp = stride2localp[stride];
      for (E_Int k = 0; k < stride; k++) {
        ptr[where] = localp[*pn]+1;
        where++;
        pn++;
      }
    }
  }

  for (auto& m : stride2arr) {
    //E_Int stride = m.first;
    //if (stride2nelem[stride] == 0) continue;
    PyList_Append(out, m.second);
    Py_DECREF(m.second);
  }

  return out;
}
