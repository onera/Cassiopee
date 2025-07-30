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
#include "common/mem.h"
#include "common/common.h"
#include <numeric>

#include "scotch/ptscotch.h"

#define EXIT \
  do { \
    MPI_Finalize(); \
    exit(0); \
  } while (0);

#define INTMAX E_IDX_NONE
#define INTMIN -(E_IDX_NONE-1)

static
E_Int get_proc(E_Int element, E_Int *distribution, E_Int nproc)
{
  for (E_Int j = 0; j < nproc; j++) {
    if (element >= distribution[j] && element < distribution[j+1])
      return j;
  }
  printf("\nWarning: could not find distribution of element " SF_D_ "\n", element);
  assert(0);
  return -1;
}

struct proc_patch {
  E_Int proc;
  std::vector<E_Int> faces;
  std::vector<E_Int> neis;

  proc_patch()
  {}

  proc_patch(E_Int pid) :
    proc(pid), faces(), neis()
  {}
};

void paraSort(E_Int *arr, int size, std::vector<E_Int> &plist_out,
  std::vector<E_Int> &pivots)
{
  int rank, nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // sort local segment
  std::sort(arr, arr+size);

  // select regularly spaced samples
  std::vector<E_Int> lsamples(nproc);
  E_Int jump = size/(nproc*nproc);
  for (E_Int i = 0; i < nproc; i++) {
    lsamples[i] = arr[i*jump];
  }

  // gather and sort all the samples
  std::vector<E_Int> gsamples(nproc*nproc);
  MPI_Allgather(&lsamples[0], nproc, XMPI_INT, &gsamples[0], nproc, XMPI_INT,
    MPI_COMM_WORLD);
  std::sort(gsamples.begin(), gsamples.end());

  // select nproc-1 pivots
  pivots.resize(nproc+1);
  for (E_Int i = 1; i < nproc; i++) {
    pivots[i] = gsamples[i*gsamples.size()/nproc];
  }
  pivots[0] = INTMIN;
  pivots[nproc] = INTMAX;

  std::vector<int> scount(nproc, 0), rcount(nproc);

  E_Int j = 0;
  for (E_Int i = 0; i < size; ) {
    if (arr[i] < pivots[j+1]) {
      scount[j]++;
      i++;
    } else {
      j++;
    }
  }
  
  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> sdist(nproc+1), rdist(nproc+1);
  sdist[0] = rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  plist_out.resize(rdist[nproc]);
  MPI_Alltoallv(arr, &scount[0], &sdist[0], XMPI_INT,
                &plist_out[0], &rcount[0], &rdist[0], XMPI_INT,
                MPI_COMM_WORLD);
  
  std::sort(plist_out.begin(), plist_out.end());

  if (plist_out.size() > 0)
    assert(plist_out[rdist[nproc]-1] < pivots[rank+1]);
}

PyObject* K_XCORE::chunk2partNGon(PyObject *self, PyObject *args)
{
  PyObject *array;
  if (!PyArg_ParseTuple(args, "O", &array)) {
    return NULL;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  clock_t tic = clock();
  
  PyObject *o, *l;
  E_Int nfld=1;
  E_Float *X, *Y, *Z;
  E_Int npoints, ncells, nfaces;
  E_Int faces_size, cells_size, *faces, *cells, *xfaces, *xcells;
  int rank, nproc;
  E_Int res;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  l = PyList_GetItem(array, 0);

  // 1 must be coordinateX chunk
  o = PyList_GetItem(l, 0);
  res = K_NUMPY::getFromNumpyArray(o, X, npoints, nfld);
  if (res != 1) { RAISE("Input error."); return NULL; };

  // 2 must be coordinateY chunk
  o = PyList_GetItem(l, 1);
  res = K_NUMPY::getFromNumpyArray(o, Y, npoints, nfld);
  if (res != 1) { RAISE("Input error."); return NULL; };
    
  // 3 must be coordinateZ chunk
  o = PyList_GetItem(l, 2);
  res = K_NUMPY::getFromNumpyArray(o, Z, npoints, nfld);
  if (res != 1) { RAISE("Input error."); return NULL; };
 
  // 4 must be ngon chunk
  o = PyList_GetItem(l, 3);
  res = K_NUMPY::getFromNumpyArray(o, faces, faces_size, nfld);
  if (res != 1) { RAISE("Input error."); return NULL; };
    
  // 5 must be ngon so chunk
  o = PyList_GetItem(l, 4);
  res = K_NUMPY::getFromNumpyArray(o, xfaces, nfaces, nfld);
  if (res != 1) { RAISE("Input error."); return NULL; };
    
  // 6 must be nface chunk
  o = PyList_GetItem(l, 5);
  res = K_NUMPY::getFromNumpyArray(o, cells, cells_size, nfld);
  if (res != 1) { RAISE("Input error."); return NULL; };
    
  // 7 must be nface so chunk
  o = PyList_GetItem(l, 6);
  res = K_NUMPY::getFromNumpyArray(o, xcells, ncells, nfld);
  if (res != 1) { RAISE("Input error."); return NULL; };

  ncells--;
  nfaces--;

  // construct cells distribution
  E_Int *cells_dist = (E_Int *)XCALLOC((nproc+1), sizeof(E_Int));
  cells_dist[0] = 0;
 
  MPI_Allgather(&ncells, 1, XMPI_INT, cells_dist+1, 1, XMPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++)
    cells_dist[i+1] += cells_dist[i];

  // construct faces distribution
  E_Int *faces_dist = (E_Int *)XCALLOC((nproc+1), sizeof(E_Int));
  faces_dist[0] = 0;
 
  MPI_Allgather(&nfaces, 1, XMPI_INT, faces_dist+1, 1, XMPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++)
    faces_dist[i+1] += faces_dist[i];

  // construct points distribution
  E_Int *points_dist = (E_Int *)XCALLOC((nproc+1), sizeof(E_Int));
  points_dist[0] = 0;
 
  MPI_Allgather(&npoints, 1, XMPI_INT, points_dist+1, 1, XMPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++)
      points_dist[i+1] += points_dist[i];

  // shift xcells and xfaces to start from zero
  E_Int cell_shift = xcells[0];
  for (E_Int i = 0; i < ncells+1; i++)
    xcells[i] -= cell_shift;

  E_Int face_shift = xfaces[0];
  for (E_Int i = 0; i < nfaces+1; i++)
    xfaces[i] -= face_shift;

  // global info
  if (rank == 0) {
    printf("Total number of cells: " SF_D_ "\n", cells_dist[nproc]);
    printf("Total number of faces: " SF_D_ "\n", faces_dist[nproc]);
    printf("Total number of points: " SF_D_ "\n", points_dist[nproc]);
  }

  E_Int sfaces_exist = 0; // switch for handling signed faces
  std::set<E_Int> SF;
  for (E_Int i = 0; i < ncells; i++) {
    for (E_Int j = xcells[i]; j < xcells[i+1]; j++) {
      E_Int face = cells[j];
      if (face < 0) {
        SF.insert(-face);
        cells[j] = -face;
      } 
    }
  }
  sfaces_exist = (SF.size() > 0);

  if (rank == 0) {
    if (sfaces_exist)
      printf("Found signed faces\n");
  }

  // make ParentElement
  std::unordered_map<E_Int, std::vector<E_Int>> PE;

  for (E_Int i = 0; i < ncells; i++) {
    for (E_Int j = xcells[i]; j < xcells[i+1]; j++) {
      E_Int face = cells[j];
      PE[face-1].push_back(i+cells_dist[rank]);
    }
  }

  // filter faces that need completing
  std::vector<int> scount(nproc, 0);
  for (const auto& face : PE) {
    E_Int target = get_proc(face.first, faces_dist, nproc);
    scount[target] += 1 + 1 + face.second.size(); // id + size + own/nei
  }

  std::vector<int> rcount(nproc);
  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> sdist(nproc+1);
  std::vector<int> rdist(nproc+1);
  sdist[0] = rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  std::vector<E_Int> sdata(sdist[nproc]);
  std::vector<E_Int> idx(nproc);
  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];

  for (const auto& face : PE) {
    E_Int target = get_proc(face.first, faces_dist, nproc);
    sdata[idx[target]++] = face.first;
    sdata[idx[target]++] = face.second.size();
    for (const auto& elem : face.second)
      sdata[idx[target]++] = elem;
  }
  
  std::vector<E_Int> rdata(rdist[nproc]);

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], XMPI_INT,
                &rdata[0], &rcount[0], &rdist[0], XMPI_INT,
                MPI_COMM_WORLD);
 
  PE.clear();

  std::vector<E_Int> A(2*nfaces, -1);
  std::vector<E_Int> c(nfaces, 0);

  E_Int first_face = faces_dist[rank]; 

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &rdata[rdist[i]];
    for (E_Int j = 0; j < rcount[i]; ) {
      E_Int face = ptr[j++] - first_face;
      E_Int size = ptr[j++];
      for (E_Int k = 0; k < size; k++) {
        A[2*face + c[face]++] = ptr[j++];
      }
    }
  }

  if (rank == 0)
    printf("ParentElements OK\n");

  // dual graph
  std::unordered_map<E_Int, std::vector<E_Int>> CADJ;
  for (E_Int i = 0; i < nfaces; i++) {
    E_Int nei = A[2*i+1];
    if (nei >= 0) {
      E_Int own = A[2*i];
      CADJ[own].push_back(nei);
      CADJ[nei].push_back(own);
    }
  }

  for (E_Int i = 0; i < nproc; i++)
    scount[i] = 0;

  for (const auto& cell : CADJ) {
    E_Int target = get_proc(cell.first, cells_dist, nproc);
    scount[target] += 1 + 1 + cell.second.size(); // cell + size + neis
  }
  
  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  sdata.resize(sdist[nproc]);
  rdata.resize(rdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];

  for (const auto& cell : CADJ) {
    E_Int target = get_proc(cell.first, cells_dist, nproc);
    sdata[idx[target]++] = cell.first;
    sdata[idx[target]++] = cell.second.size();
    for (const auto& elem : cell.second)
      sdata[idx[target]++] = elem;
  }

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], XMPI_INT,
                &rdata[0], &rcount[0], &rdist[0], XMPI_INT,
                MPI_COMM_WORLD);

  E_Int first_cell = cells_dist[rank];
  std::vector<std::vector<E_Int>> cadj(ncells);
  E_Int nedges = 0;
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &rdata[rdist[i]];
    for (E_Int j = 0; j < rcount[i]; ) {
      E_Int cell = ptr[j++];
      E_Int psize = ptr[j++];
      nedges += psize;
      for (E_Int k = 0; k < psize; k++)
        cadj[cell - first_cell].push_back(ptr[j++]);
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

  assert(sizeof(SCOTCH_Num) <= sizeof(SCOTCH_Idx));
  assert(sizeof(SCOTCH_Idx) >= sizeof(void *));
  assert(sizeof(SCOTCH_Num) == sizeof(E_Int));
  assert(SCOTCH_numSizeof() == sizeof(SCOTCH_Num));

  int ret;
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

  SCOTCH_dgraphExit(&graph);
  SCOTCH_stratExit(&strat);

  // Send cells to their target proc!
  std::vector<int> c_scount(nproc, 0);
  std::vector<int> c_rcount(nproc, 0);

  for (E_Int i = 0; i < ncells; i++) {
    E_Int where = part[i];
    c_scount[where] += 1;
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
  std::vector<E_Int> c_stride(c_sdist[nproc]);
  std::vector<E_Int> nxcells(nncells+1, 0);

  for (E_Int i = 0; i < nproc; i++) idx[i] = c_sdist[i];

  for (E_Int i = 0; i < ncells; i++) {
    E_Int where = part[i];
    scells[idx[where]] = i + cells_dist[rank];
    c_stride[idx[where]] = xcells[i+1] - xcells[i];
    idx[where]++;
  }

  MPI_Alltoallv(&scells[0], &c_scount[0], &c_sdist[0], XMPI_INT,
                &rcells[0], &c_rcount[0], &c_rdist[0], XMPI_INT,
                MPI_COMM_WORLD);
  
  MPI_Alltoallv(&c_stride[0], &c_scount[0], &c_sdist[0], XMPI_INT,
                &nxcells[0]+1, &c_rcount[0], &c_rdist[0], XMPI_INT,
                MPI_COMM_WORLD);
  
  nxcells[0] = 0;
  for (E_Int i = 0; i < nncells; i++)
    nxcells[i+1] += nxcells[i];

  // send NFACE
  for (E_Int i = 0; i < nproc; i++) {
    scount[i] = 0;
    rcount[i] = 0;
  }

  for (E_Int i = 0; i < ncells; i++) {
    E_Int where = part[i];
    scount[where] += xcells[i+1] - xcells[i];
  }

  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  sdist[0] = rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  sdata.resize(sdist[nproc]);
  std::vector<E_Int> NFACE(rdist[nproc]);

  for (E_Int i = 0; i < nproc; i++) idx[i] = sdist[i];

  if (sfaces_exist) {
    for (E_Int i = 0; i < ncells; i++) {
      E_Int where = part[i];
      for (E_Int j = xcells[i]; j < xcells[i+1]; j++) {
        E_Int face = cells[j];
        if (SF.find(face) != SF.end())
          sdata[idx[where]++] = -face;
        else
          sdata[idx[where]++] = face;
      }
    }
  } else {
    for (E_Int i = 0; i < ncells; i++) {
      E_Int where = part[i];
      for (E_Int j = xcells[i]; j < xcells[i+1]; j++)
        sdata[idx[where]++] = cells[j];
    }
  }

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], XMPI_INT,
                &NFACE[0], &rcount[0], &rdist[0], XMPI_INT,
                MPI_COMM_WORLD);

  if (rank == 0)
    puts("NFACE OK");

  // hash cells
  std::unordered_map<E_Int, E_Int> CT;
  for (E_Int i = 0; i < nncells; i++)
    CT[rcells[i]] = i;
  
  // store signed faces again
  if (sfaces_exist) {
    SF.clear();
    for (E_Int i = 0; i < nncells; i++) {
      for (E_Int j = nxcells[i]; j < nxcells[i+1]; j++) {
        if (NFACE[j] < 0) {
          SF.insert(-NFACE[j]);
          NFACE[j] = -NFACE[j];
        }
      }
    }
  }

  // hash and request faces
  E_Int nnfaces = 0;
  std::unordered_map<E_Int, E_Int> FT; // to avoid face duplication
  std::vector<int> f_rcount(nproc, 0), f_scount(nproc, 0);

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pc = &rcells[c_rdist[i]];
    for (E_Int j = 0; j < c_rcount[i]; j++) {
      E_Int lc = CT[pc[j]];
      for (E_Int k = nxcells[lc]; k < nxcells[lc+1]; k++) {
        E_Int face = NFACE[k];
        if (FT.find(face) == FT.end()) {
          FT[face] = nnfaces++;
          E_Int src = get_proc(face-1, faces_dist, nproc);
          f_rcount[src]++;
        }
      }
    }
  }

  MPI_Alltoall(&f_rcount[0], 1, MPI_INT, &f_scount[0], 1, MPI_INT, MPI_COMM_WORLD);
  
  std::vector<int> f_rdist(nproc+1);
  std::vector<int> f_sdist(nproc+1);
  f_rdist[0] = f_sdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    f_rdist[i+1] = f_rdist[i] + f_rcount[i];
    f_sdist[i+1] = f_sdist[i] + f_scount[i];
  }

  assert(nnfaces == f_rdist[nproc]);

  std::vector<E_Int> rfaces(nnfaces);
  std::vector<E_Int> sfaces(f_sdist[nproc]);

  for (E_Int i = 0; i < nproc; i++) idx[i] = f_rdist[i];

  std::vector<E_Int> vfaces(nnfaces, 0);
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pc = &rcells[c_rdist[i]];
    for (E_Int j = 0; j < c_rcount[i]; j++) {
      E_Int lc = CT[pc[j]];
      for (E_Int k = nxcells[lc]; k < nxcells[lc+1]; k++) {
        E_Int face = NFACE[k];
        if (!vfaces[FT[face]]) {
          vfaces[FT[face]]++;
          E_Int src = get_proc(face-1, faces_dist, nproc);
          rfaces[idx[src]++] = face;
        }
      }
    }
  }

  MPI_Alltoallv(&rfaces[0], &f_rcount[0], &f_rdist[0], XMPI_INT,
                &sfaces[0], &f_scount[0], &f_sdist[0], XMPI_INT,
                MPI_COMM_WORLD);

  FT.clear();
  nnfaces = 0; 
  for (const auto& face : rfaces)
    FT[face] = nnfaces++;
  
  // send face strides
  std::vector<E_Int> f_stride(f_sdist[nproc]);
  std::vector<E_Int> nxfaces(nnfaces+1, 0);
  for (E_Int i = 0; i < nproc; i++) scount[i] = 0;

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pf = &sfaces[f_sdist[i]];
    E_Int *ps = &f_stride[f_sdist[i]];
    for (E_Int j = 0; j < f_scount[i]; j++) {
      E_Int face = pf[j]-1-faces_dist[rank];
      E_Int stride = xfaces[face+1] - xfaces[face];
      ps[j] = stride;
      scount[i] += stride;
    }
  }

  MPI_Alltoallv(&f_stride[0], &f_scount[0], &f_sdist[0], XMPI_INT,
                &nxfaces[0]+1, &f_rcount[0], &f_rdist[0], XMPI_INT,
                MPI_COMM_WORLD);
  
  nxfaces[0] = 0;
  for (E_Int i = 0; i < nnfaces; i++)
    nxfaces[i+1] += nxfaces[i];
  
  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  // send NGON
  sdist[0] = rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  sdata.resize(sdist[nproc]);
  std::vector<E_Int> NGON(rdist[nproc]);

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pf = &sfaces[f_sdist[i]];
    E_Int *pn = &sdata[sdist[i]];
    for (E_Int j = 0; j < f_scount[i]; j++) {
      E_Int face = pf[j] - 1 - faces_dist[rank];
      for (E_Int k = xfaces[face]; k < xfaces[face+1]; k++)
        *pn++ = faces[k];
    }
  }

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], XMPI_INT,
                &NGON[0], &rcount[0], &rdist[0], XMPI_INT,
                MPI_COMM_WORLD);
  
  if (rank == 0)
    puts("NGON OK");

  // request points
  std::unordered_map<E_Int, E_Int> PT; // avoid point repetition
  E_Int nnpoints = 0;
  std::vector<int> p_rcount(nproc, 0);
  std::vector<int> p_scount(nproc);
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pf = &rfaces[f_rdist[i]];
    for (E_Int j = 0; j < f_rcount[i]; j++) {
      E_Int face = FT[pf[j]];
      for (E_Int k = nxfaces[face]; k < nxfaces[face+1]; k++) {
        E_Int point = NGON[k];
        if (PT.find(point) == PT.end()) {
          PT[point] = nnpoints++;
          E_Int src = get_proc(point-1, points_dist, nproc);
          p_rcount[src]++;
        }
      }
    }
  }

  MPI_Alltoall(&p_rcount[0], 1, MPI_INT, &p_scount[0], 1, MPI_INT,
    MPI_COMM_WORLD);

  std::vector<int> p_sdist(nproc+1);
  std::vector<int> p_rdist(nproc+1);
  p_sdist[0] = p_rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    p_sdist[i+1] = p_sdist[i] + p_scount[i];
    p_rdist[i+1] = p_rdist[i] + p_rcount[i];
  }
  std::vector<E_Int> spoints(p_sdist[nproc]);
  std::vector<E_Int> rpoints(p_rdist[nproc]);
  std::vector<E_Int> vpoints(nnpoints, 0);
  assert(nnpoints == p_rdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = p_rdist[i];

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pf = &rfaces[f_rdist[i]];
    for (E_Int j = 0; j < f_rcount[i]; j++) {
      E_Int face = FT[pf[j]];
      for (E_Int k = nxfaces[face]; k < nxfaces[face+1]; k++) {
        E_Int point = NGON[k];
        if (!vpoints[PT[point]]) {
          vpoints[PT[point]]++;
          E_Int src = get_proc(point-1, points_dist, nproc);
          rpoints[idx[src]++] = point;
        }
      }
    }
  }

  MPI_Alltoallv(&rpoints[0], &p_rcount[0], &p_rdist[0], XMPI_INT,
                &spoints[0], &p_scount[0], &p_sdist[0], XMPI_INT,
                MPI_COMM_WORLD);

  // renumber points
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
  assert(sdist[nproc] == 3*p_sdist[nproc]);
  assert(rdist[nproc] == 3*p_rdist[nproc]);
  
  std::vector<E_Float> sxyz(sdist[nproc]);
  std::vector<E_Float> rxyz(rdist[nproc]);

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pp = &spoints[p_sdist[i]];
    E_Float *px = &sxyz[sdist[i]];
    for (E_Int j = 0; j < p_scount[i]; j++) {
      E_Int point = pp[j] - 1 - points_dist[rank];
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

  // construct communication patches
  for (E_Int i = 0; i < nproc; i++) {
    scount[i] = 2*f_scount[i];
    rcount[i] = 2*f_rcount[i];
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  sdata.resize(sdist[nproc]);
  rdata.resize(rdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &sfaces[f_sdist[i]];
    for (E_Int j = 0; j < f_scount[i]; j++) {
      E_Int face = ptr[j];
      assert(rank == get_proc(face-1, faces_dist, nproc));
      face -= 1 + faces_dist[rank];
      E_Int own = A[2*face];
      E_Int nei = A[2*face+1];
      sdata[idx[i]++] = own;
      sdata[idx[i]++] = nei;
    }
  }

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], XMPI_INT,
                &rdata[0], &rcount[0], &rdist[0], XMPI_INT,
                MPI_COMM_WORLD);


  XFREE(faces_dist);

  std::vector<E_Int> pneis; 

  E_Int nif = 0;

  std::vector<int> rncount(nproc,0);
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &rdata[rdist[i]];
    E_Int j = 0;
    for (; j < rcount[i];) {
      E_Int own = ptr[j++];
      E_Int nei = ptr[j++];

      assert(own != -1);

      if (nei == -1) {
        assert(CT.find(own) != CT.end());
        continue; // boundary face
      }

      E_Int ho = CT.find(own) != CT.end();
      E_Int hn = CT.find(nei) != CT.end();

      assert(ho || hn);

      if (ho && hn) {
        nif++;
        continue; // internal face
      }

      if (ho && !hn) {
        // nei is remote neighbor
        E_Int src = get_proc(nei, cells_dist, nproc);
        rncount[src]++;
      } else if (!ho && hn) {
        E_Int src = get_proc(own, cells_dist, nproc);
        rncount[src]++;
      } else {
        assert(0);
      }
    }
    assert(j == rcount[i]);
  }

  std::vector<int> sncount(nproc);
  MPI_Alltoall(&rncount[0], 1, MPI_INT, &sncount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<int> sndist(nproc+1);
  std::vector<int> rndist(nproc+1);
  sndist[0] = rndist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    rndist[i+1] = rndist[i] + rncount[i];
    sndist[i+1] = sndist[i] + sncount[i];
  }

  std::vector<E_Int> sncells(sndist[nproc]);
  std::vector<E_Int> rncells(rndist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = rndist[i];

  // pfaces follow rncells dist
  std::vector<E_Int> pfaces(rndist[nproc]);

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &rdata[rdist[i]];
    E_Int *pf = &rfaces[f_rdist[i]];
    E_Int j = 0;
    E_Int count = 0;
    for (; j < rcount[i];) {
      E_Int face = pf[count++];
      E_Int own = ptr[j++];
      E_Int nei = ptr[j++];

      if (nei == -1) {
        continue; // boundary face
      }

      E_Int ho = CT.find(own) != CT.end();
      E_Int hn = CT.find(nei) != CT.end();


      if (ho && hn) {
        continue; // face is internal
      }

      if (ho && !hn) {
        // nei is remote neighbor
          E_Int src = get_proc(nei, cells_dist, nproc);
          rncells[idx[src]] = nei;
          pfaces[idx[src]] = face;
          idx[src]++;
      } else if (!ho && hn) {
        // own is remote neighbor
          E_Int src = get_proc(own, cells_dist, nproc);
          rncells[idx[src]] = own;
          pfaces[idx[src]] = face;
          idx[src]++;
      } else {
        assert(0);
      }
    }
  }

  MPI_Alltoallv(&rncells[0], &rncount[0], &rndist[0], XMPI_INT,
                &sncells[0], &sncount[0], &sndist[0], XMPI_INT,
                MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &sncells[sndist[i]];
    for (E_Int j = 0; j < sncount[i]; j++) {
      E_Int cell = ptr[j];
      assert(rank == get_proc(cell, cells_dist, nproc));
      cell -= cells_dist[rank];
      E_Int where = part[cell];
      ptr[j] = where;
    }
  }

  std::vector<E_Int> nei2proc(rndist[nproc]);
  MPI_Alltoallv(&sncells[0], &sncount[0], &sndist[0], XMPI_INT,
                &nei2proc[0], &rncount[0], &rndist[0], XMPI_INT,
                MPI_COMM_WORLD);

  // we have a one to one map between pfaces, pneis and nei2proc via rndist and rncount
  // construct proc_patches
  proc_patch **ppatches = new proc_patch *[nproc]; // nproc is max size;
  std::unordered_map<E_Int, E_Int> patch_id; // global proc to local proc
  E_Int npatches = 0;
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pf = &pfaces[rndist[i]];
    E_Int *pn = &rncells[rndist[i]];
    E_Int *pp = &nei2proc[rndist[i]];
    for (E_Int j = 0; j < rncount[i]; j++) {
      E_Int face = pf[j];
      E_Int nei = pn[j];
      E_Int proc = pp[j];

      auto ppatch = patch_id.find(proc);

      if (ppatch == patch_id.end()) {
        proc_patch *PP = new proc_patch(proc);
        PP->faces.push_back(face);
        PP->neis.push_back(nei);

        ppatches[npatches] = PP;

        patch_id[proc] = npatches++;
      } else {
        ppatches[ppatch->second]->faces.push_back(face);
        ppatches[ppatch->second]->neis.push_back(nei);
      }
    }
  }

  // Sort pfaces in ascending order, then replace them with local indices
  // Note(Imad): this is in order to get rid of huge FaceLoc2Glob array
  // TODO(Imad): is gneis ever useful?
  for (E_Int i = 0; i < npatches; i++) {
    auto &pfaces = ppatches[i]->faces;
    auto &gneis = ppatches[i]->neis;

    std::vector<E_Int> indices(pfaces.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
      [&](E_Int a, E_Int b) { return pfaces[a] < pfaces[b]; });
    
    std::vector<E_Int> sorted_pfaces(pfaces.size());
    std::vector<E_Int> sorted_gneis(pfaces.size());
    for (size_t j = 0; j < pfaces.size(); j++) {
      sorted_pfaces[j] = pfaces[indices[j]];
      sorted_gneis[j] = gneis[indices[j]];
    }

    for (size_t j = 0; j < pfaces.size(); j++) {
      pfaces[j] = FT[sorted_pfaces[j]] + 1;
      gneis[j] = sorted_gneis[j];

    }
  }

  if (rank == 0)
    puts("Comm patches OK");

  // Build output Python list
  PyObject* out = PyList_New(0);

  const char *varString = "CoordinateX,CoordinateY,CoordinateZ";
  // TODO(Imad): avoid copying
  PyObject* m = K_ARRAY::buildArray3(3, varString, nnpoints, nncells, nnfaces,
    "NGON", nxfaces[nnfaces], nxcells[nncells], 3, false, 3);

  K_FLD::FldArrayF *f; K_FLD::FldArrayI *cn;
  K_ARRAY::getFromArray3(m, f, cn); 

  for (E_Int n = 0; n < 3; n++) {
    E_Float *pt = f->begin(n+1);
    for (E_Int i = 0; i < nnpoints; i++)
      pt[i] = rxyz[3*i+n];
  }

  E_Int *ngon = cn->getNGon();
  E_Int *nface = cn->getNFace();
  E_Int *indPG = cn->getIndPG();
  E_Int *indPH = cn->getIndPH();

  for (E_Int i = 0; i <= nnfaces; i++) indPG[i] = nxfaces[i];
  for (E_Int i = 0; i <= nncells; i++) indPH[i] = nxcells[i];

  E_Int* ptr = ngon;

  for (E_Int i = 0; i < nnfaces; i++) {
    E_Int start = nxfaces[i];
    E_Int end = nxfaces[i+1];
    for (E_Int j = start; j < end; j++) 
      *ptr++ = PT[NGON[j]]+1;
  }

  ptr = nface;

  if (sfaces_exist) {
    for (E_Int i = 0; i < nncells; i++) {
      E_Int start = nxcells[i];
      E_Int end = nxcells[i+1];
      for (E_Int j = start; j < end; j++) {
        E_Int face = NFACE[j];
        if (SF.find(face) == SF.end())
          *ptr++ = FT[NFACE[j]]+1;
        else
          *ptr++ = -(FT[NFACE[j]]+1);
      }
    }
  } else {
    for (E_Int i = 0; i < nncells; i++) {
      E_Int start = nxcells[i];
      E_Int end = nxcells[i+1];
      for (E_Int j = start; j < end; j++) 
        *ptr++ = FT[NFACE[j]]+1;
    }
  }

  PyList_Append(out, m);  
  Py_DECREF(m);
  RELEASESHAREDU(m, f, cn);

  PyObject* comm_data = PyList_New(0);

  for (E_Int i = 0; i < npatches; i++) {
    proc_patch *PP = ppatches[i];
    E_Int proc = PP->proc;
    const auto& faces = PP->faces;
    const auto& neis = PP->neis;
  
    npy_intp dims[2];

    // faces array
    dims[1] = 1;
    dims[0] = (npy_intp)faces.size();
 
    PyArrayObject *f = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
    E_Int *pf = (E_Int *)PyArray_DATA(f);
    for (size_t j = 0; j < faces.size(); j++)
      pf[j] = faces[j];

    // neis array
    PyArrayObject *n = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
    E_Int *pn = (E_Int *)PyArray_DATA(n);
    for (size_t j = 0; j < neis.size(); j++)
      pn[j] = neis[j];
  
    PyObject *arr = Py_BuildValue("[lOO]", proc, f, n);
    Py_DECREF(f);
    Py_DECREF(n);
    PyList_Append(comm_data, arr);

    delete ppatches[i];
  }

  delete [] ppatches; 

  

  PyList_Append(out, comm_data);
  Py_DECREF(comm_data);

  npy_intp dims[2];

  // 8 must be an array of FlowSolutions#Centers chunks
  o = PyList_GetItem(l, 7);
  E_Int csize = PyList_Size(o);
  if (csize == 0) {
    PyList_Append(out, PyList_New(0));
  } else {
    E_Float **csols = (E_Float **)XCALLOC(csize, sizeof(E_Float *));

    for (E_Int i = 0; i < csize; i++) {
      PyObject *csol = PyList_GetItem(o, i);
      res = K_NUMPY::getFromNumpyArray(csol, csols[i], ncells, nfld);
      if (res != 1) { RAISE("Bad input."); return NULL; }
    }

    PyObject *clist = PyList_New(0);

    dims[1] = 1;
    dims[0] = (npy_intp)nncells;
    
    std::vector<E_Float> c_sdata(c_sdist[nproc]);

    for (E_Int k = 0; k < csize; k++) {
      PyArrayObject *ca = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);

      for (E_Int i = 0; i < nproc; i++)
        idx[i] = c_sdist[i];
      
      for (E_Int i = 0; i < nproc; i++) {
        E_Int *ptr = &scells[c_sdist[i]];
        for (E_Int j = 0; j < c_scount[i]; j++) {
          E_Int cell = ptr[j] - cells_dist[rank];
          assert(rank == get_proc(ptr[j], cells_dist, nproc));
          c_sdata[idx[i]++] = csols[k][cell];
        }
      }

      MPI_Alltoallv(&c_sdata[0], &c_scount[0], &c_sdist[0], MPI_DOUBLE,
                    PyArray_DATA(ca), &c_rcount[0], &c_rdist[0], MPI_DOUBLE,
                    MPI_COMM_WORLD);

      PyList_Append(clist, (PyObject *)ca);
      Py_DECREF(ca);
    }

    PyList_Append(out, clist);
    Py_DECREF(clist);
  }

  XFREE(cells_dist);

  // 9 must be an array of FlowSolutions chunks
  o = PyList_GetItem(l, 8);
  E_Int psize = PyList_Size(o);
  if (psize == 0) {
    PyList_Append(out, PyList_New(0));
  } else {
    E_Float **psols = (E_Float **)XCALLOC(psize, sizeof(E_Float *));

    for (E_Int i = 0; i < psize; i++) {
      PyObject *psol = PyList_GetItem(o, i);
      res = K_NUMPY::getFromNumpyArray(psol, psols[i], npoints, nfld);
      if (res != 1) { RAISE("Input error."); return NULL; };
    }

    PyObject *plist = PyList_New(0);

    dims[1] = 1;
    dims[0] = (npy_intp)nnpoints;
    
    std::vector<E_Float> p_sdata(p_sdist[nproc]);

    for (E_Int k = 0; k < psize; k++) {
      PyArrayObject *pa = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);

      for (E_Int i = 0; i < nproc; i++)
        idx[i] = p_sdist[i];
      
      for (E_Int i = 0; i < nproc; i++) {
        E_Int *ptr = &spoints[p_sdist[i]];
        for (E_Int j = 0; j < p_scount[i]; j++) {
          E_Int point = ptr[j] - points_dist[rank] - 1;
          assert(rank == get_proc(ptr[j]-1, points_dist, nproc));
          p_sdata[idx[i]++] = psols[k][point];
        }
      }

      MPI_Alltoallv(&p_sdata[0], &p_scount[0], &p_sdist[0], MPI_DOUBLE,
                    PyArray_DATA(pa), &p_rcount[0], &p_rdist[0], MPI_DOUBLE,
                    MPI_COMM_WORLD);

      PyList_Append(plist, (PyObject *)pa);
      Py_DECREF(pa);
    }

    PyList_Append(out, plist);
    Py_DECREF(plist);
  }

  XFREE(points_dist);

  // 10 must be an array of PointList chunks
  o = PyList_GetItem(l, 9);
  E_Int nbc = PyList_Size(o);

  if (nbc == 0) {
    PyList_Append(out, PyList_New(0));
  } else {
    E_Int **plists = (E_Int **)XCALLOC(nbc, sizeof(E_Int *));
    int *bcsize = (int *)XMALLOC(nbc * sizeof(int));
    std::vector<std::vector<E_Int>> myptlists(nbc);
    std::vector<std::vector<E_Int>> pivots(nbc);
    E_Int size;

    for (E_Int i = 0; i < nbc; i++) {
      PyObject *plist = PyList_GetItem(o, i);
      res = K_NUMPY::getFromNumpyArray(plist, plists[i], size, nfld);
      if (res != 1) { RAISE("Input error."); return NULL; };
      bcsize[i] = int(size);

      auto& list = plists[i];
      paraSort(&list[0], bcsize[i], myptlists[i], pivots[i]);

      // Note(Imad): we could free up plists[i] and bcsize[i] here
    }

    XFREE(plists);
    XFREE(bcsize);

    // make local PE
    std::unordered_map<E_Int, std::vector<E_Int>> lPE;
    for (E_Int i = 0; i < nncells; i++) {
      for (E_Int j = nxcells[i]; j < nxcells[i+1]; j++) {
        E_Int face = NFACE[j];
        lPE[face].push_back(i);
      }
    }

    // isolate bfaces
    std::set<E_Int> pfaces_set(pfaces.begin(), pfaces.end());
    std::vector<E_Int> bfaces;

    for (auto& f : lPE) {
      if (f.second.size() == 1) {
        if (pfaces_set.find(f.first) == pfaces_set.end()) {
          // not a pface
          bfaces.push_back(f.first);
        }
      }
    }

    // request bface info for each bc
    PyObject *bclist_out = PyList_New(0);
    for (E_Int bc = 0; bc < nbc; bc++) {
      std::vector<int> scount(nproc, 0), rcount(nproc, 0), sdist(nproc+1), rdist(nproc+1);
      
      auto& pivot = pivots[bc];
      for (auto bface : bfaces) {
        E_Int src = get_proc(bface, &pivot[0], nproc);
        scount[src]++;
      }
      
      MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);
      
      sdist[0] = rdist[0] = 0;
      for (E_Int i = 0; i < nproc; i++) {
        sdist[i+1] = sdist[i] + scount[i];
        rdist[i+1] = rdist[i] + rcount[i];
      }

      std::vector<E_Int> sdata(sdist[nproc]), rdata(rdist[nproc]);
      std::vector<int> idx(sdist);
      
      for (auto bface : bfaces) {
        E_Int src = get_proc(bface, &pivot[0], nproc);
        sdata[idx[src]++] = bface;
      }

      MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], XMPI_INT,
                    &rdata[0], &rcount[0], &rdist[0], XMPI_INT,
                    MPI_COMM_WORLD);
      
      auto& myptlist = myptlists[bc];
      std::set<E_Int> pointlist_set(myptlist.begin(), myptlist.end());

      std::vector<int> bscount(nproc, 0), brcount(nproc), bsdist(nproc+1), brdist(nproc+1);

      for (E_Int i = 0; i < nproc; i++) {
        E_Int *pf = &rdata[rdist[i]];
        for (E_Int j = 0; j < rcount[i]; j++) {
          E_Int bface = pf[j];
          // is bface in pointlist_set?
          if (pointlist_set.find(bface) != pointlist_set.end())
            bscount[i]++;
        }
      }

      MPI_Alltoall(&bscount[0], 1, MPI_INT, &brcount[0], 1, MPI_INT, MPI_COMM_WORLD);

      bsdist[0] = brdist[0] = 0;
      for (E_Int i = 0; i < nproc; i++) {
        bsdist[i+1] = bsdist[i] + bscount[i];
        brdist[i+1] = brdist[i] + brcount[i];
      }

      std::vector<E_Int> bsdata(bsdist[nproc]);

      for (E_Int i = 0; i < nproc; i++) idx[i] = bsdist[i];

      for (E_Int i = 0; i < nproc; i++) {
        E_Int *pf = &rdata[rdist[i]];
        for (E_Int j = 0; j < rcount[i]; j++) {
          E_Int bface = pf[j];
          if (pointlist_set.find(bface) != pointlist_set.end())
            bsdata[idx[i]++] = bface;
        }
      }

      int nrecv = brdist[nproc];
      dims[1] = 1;
      dims[0] = (npy_intp)nrecv;
      PyArrayObject *pa = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);

      MPI_Alltoallv(&bsdata[0], &bscount[0], &bsdist[0], XMPI_INT,
                    PyArray_DATA(pa), &brcount[0], &brdist[0], XMPI_INT,
                    MPI_COMM_WORLD);

      E_Int *ppa = (E_Int*)PyArray_DATA(pa);
      for (E_Int i = 0; i < nrecv; i++) {
        assert(FT.find(ppa[i]) != FT.end());
        ppa[i] = FT[ppa[i]]+1;
      }

      PyList_Append(bclist_out, (PyObject *)pa);
      Py_DECREF(pa);
    }

    PyList_Append(out, bclist_out);
    Py_DECREF(bclist_out);
  }

  // my global cells
  dims[1] = 1;
  dims[0] = (npy_intp)nncells;
  PyArrayObject *mycells = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
  E_Int *pc = (E_Int *)PyArray_DATA(mycells);
  for (E_Int i = 0; i < nncells; i++)
    pc[i] = rcells[i];

  // my global faces
  dims[0] = (npy_intp)nnfaces;
  PyArrayObject *myfaces = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
  E_Int *pf = (E_Int *)PyArray_DATA(myfaces);
  for (E_Int i = 0; i < nnfaces; i++)
    pf[i] = rfaces[i];
 
  // my global points
  dims[0] = (npy_intp)nnpoints;
  PyArrayObject *mypoints = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
  E_Int *pp = (E_Int *)PyArray_DATA(mypoints);
  for (E_Int i = 0; i < nnpoints; i++)
    pp[i] = rpoints[i];

  PyList_Append(out, (PyObject *)mycells); 
  PyList_Append(out, (PyObject *)myfaces); 
  PyList_Append(out, (PyObject *)mypoints); 

  Py_DECREF(mycells); 
  Py_DECREF(myfaces);
  Py_DECREF(mypoints);

  if (rank == 0)
    puts("Export OK");
  
  MPI_Barrier(MPI_COMM_WORLD);
  clock_t toc = clock();
  E_Float ptime = ((E_Float)(toc-tic)) / CLOCKS_PER_SEC;
  if (rank == 0) {
    printf("Partitioned mesh in %.2f s\n", ptime);
  }

  return out;
}
