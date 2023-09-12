/*    
    Copyright 2013-2023 Onera.

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
#include "../common/mem.h"

#include "scotch/ptscotch.h"

#define EXIT \
  do { \
    MPI_Finalize(); \
    exit(0); \
  } while (0);

static
E_Int get_proc(E_Int element, E_Int *distribution, E_Int nproc)
{
  for (E_Int j = 0; j < nproc; j++) {
    if (element >= distribution[j] && element < distribution[j+1])
      return j;
  }
  printf("\nWarning: could not find distribution of element %d\n", element);
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

PyObject* K_XCORE::chunk2part(PyObject *self, PyObject *args)
{
  PyObject *array;
  if (!PyArg_ParseTuple(args, "O", &array)) {
    return NULL;
  }

  E_Int nzones = PyList_Size(array);
  if (nzones != 1) {
    fprintf(stderr, "chunk2part(): should be one zone per chunk for now.\n");
    return NULL;
  }

  PyObject *o, *l;
  E_Int nfld;
  E_Float *X, *Y, *Z;
  E_Int npoints, ncells, nfaces;
  E_Int faces_size, cells_size, *faces, *cells, *xfaces, *xcells;
  E_Int rank, nproc, res;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  l = PyList_GetItem(array, 0);

  // 1 must be coordinateX chunk
  o = PyList_GetItem(l, 0);
  res = K_NUMPY::getFromNumpyArray(o, X, npoints, nfld, true);
  assert(res == 1);

  // 2 must be coordinateY chunk
  o = PyList_GetItem(l, 1);
  res = K_NUMPY::getFromNumpyArray(o, Y, npoints, nfld, true);
  assert(res == 1);
    
  // 3 must be coordinateZ chunk
  o = PyList_GetItem(l, 2);
  res = K_NUMPY::getFromNumpyArray(o, Z, npoints, nfld, true);
  assert(res == 1);
 
  // 4 must be ngon chunk
  o = PyList_GetItem(l, 3);
  res = K_NUMPY::getFromNumpyArray(o, faces, faces_size, nfld, true);
  assert(res == 1);  
    
  // 5 must be ngon so chunk
  o = PyList_GetItem(l, 4);
  res = K_NUMPY::getFromNumpyArray(o, xfaces, nfaces, nfld, true);
  assert(res == 1);
    
  // 6 must be nface chunk
  o = PyList_GetItem(l, 5);
  res = K_NUMPY::getFromNumpyArray(o, cells, cells_size, nfld, true);
  assert(res == 1);
    
  // 7 must be nface so chunk
  o = PyList_GetItem(l, 6);
  res = K_NUMPY::getFromNumpyArray(o, xcells, ncells, nfld, true);
  assert(res == 1);

  ncells--;
  nfaces--;

  // construct cells distribution
  E_Int *cells_dist = (E_Int *)XCALLOC((nproc+1), sizeof(E_Int));
  cells_dist[0] = 0;
 
  MPI_Allgather(&ncells, 1, MPI_INT, cells_dist+1, 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++)
    cells_dist[i+1] += cells_dist[i];

  // construct faces distribution
  E_Int *faces_dist = (E_Int *)XCALLOC((nproc+1), sizeof(E_Int));
  faces_dist[0] = 0;
 
  MPI_Allgather(&nfaces, 1, MPI_INT, faces_dist+1, 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++)
    faces_dist[i+1] += faces_dist[i];

  // construct points distribution
  E_Int *points_dist = (E_Int *)XCALLOC((nproc+1), sizeof(E_Int));
  points_dist[0] = 0;
 
  MPI_Allgather(&npoints, 1, MPI_INT, points_dist+1, 1, MPI_INT, MPI_COMM_WORLD);

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
    printf("total number of cells: %d\n", cells_dist[nproc]);
    printf("total number of faces: %d\n", faces_dist[nproc]);
    printf("total number of points: %d\n", points_dist[nproc]);
  }

  // check for signed faces and holes in face numbering
  E_Int sfaces_exist = 0; // switch for handling signed faces
  std::set<E_Int> SF;
  for (E_Int i = 0; i < ncells; i++) {
    for (E_Int j = xcells[i]; j < xcells[i+1]; j++) {
      E_Int face = cells[j];
      if (face < 0) {
        if (SF.find(-face) == SF.end())
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
  std::vector<E_Int> scount(nproc, 0);
  for (const auto& face : PE) {
    E_Int target = get_proc(face.first, faces_dist, nproc);
    scount[target] += 1 + 1 + face.second.size(); // id + size + own/nei
  }

  std::vector<E_Int> rcount(nproc);
  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> sdist(nproc+1);
  std::vector<E_Int> rdist(nproc+1);
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

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], MPI_INT,
                &rdata[0], &rcount[0], &rdist[0], MPI_INT,
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

  E_Int count = 0;
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

  count = 0;
  for (const auto& cell : CADJ) {
    E_Int target = get_proc(cell.first, cells_dist, nproc);
    sdata[idx[target]++] = cell.first;
    sdata[idx[target]++] = cell.second.size();
    for (const auto& elem : cell.second)
      sdata[idx[target]++] = elem;
  }

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], MPI_INT,
                &rdata[0], &rcount[0], &rdist[0], MPI_INT,
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
  count = 0;
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

  // Send cells to their target proc!
  std::vector<E_Int> c_scount(nproc, 0);
  std::vector<E_Int> c_rcount(nproc);

  // start with cell ids
  for (E_Int i = 0; i < ncells; i++) {
    E_Int where = part[i];
    c_scount[where]++;
  }

  MPI_Alltoall(&c_scount[0], 1, MPI_INT, &c_rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> c_sdist(nproc+1);
  std::vector<E_Int> c_rdist(nproc+1);
  
  c_sdist[0] = c_rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    c_sdist[i+1] = c_sdist[i] + c_scount[i];
    c_rdist[i+1] = c_rdist[i] + c_rcount[i];
  }

  // Note(Imad): new ncells may differ from cells_dist[rank+1] - cells_dist[rank]
  E_Int nncells = c_rdist[nproc];

  std::vector<E_Int> scells(c_sdist[nproc]);
  std::vector<E_Int> rcells(c_rdist[nproc]);
  std::vector<E_Int> sstride(c_sdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = c_sdist[i];

  for (E_Int i = 0; i < ncells; i++) {
    E_Int where = part[i];
    scells[idx[where]] = i + cells_dist[rank];
    sstride[idx[where]] = xcells[i+1] - xcells[i];
    idx[where]++;
  }

  MPI_Alltoallv(&scells[0], &c_scount[0], &c_sdist[0], MPI_INT,
                &rcells[0], &c_rcount[0], &c_rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  std::vector<E_Int> nxcells(nncells+1);

  MPI_Alltoallv(&sstride[0], &c_scount[0], &c_sdist[0], MPI_INT,
                &nxcells[0]+1, &c_rcount[0], &c_rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  nxcells[0] = 0;
  for (E_Int i = 0; i < nncells; i++)
    nxcells[i+1] += nxcells[i];

  // send NFACE connectivity
  for (E_Int i = 0; i < nproc; i++)
    scount[i] = 0;

  for (E_Int i = 0; i < ncells; i++) {
     E_Int where = part[i];
     scount[where] += xcells[i+1] - xcells[i];
  }

  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  sdata.resize(sdist[nproc]);
  std::vector<E_Int> NFACE(rdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];

  for (E_Int i = 0; i < ncells; i++) {
    E_Int where = part[i];
    for (E_Int j = xcells[i]; j < xcells[i+1]; j++)
      sdata[idx[where]++] = cells[j];
  }

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], MPI_INT,
                &NFACE[0], &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  
  if (rank == 0)
    printf("NFACE OK\n");

  // Hash and request faces
  std::unordered_map<E_Int, E_Int> FT;
  E_Int nnfaces = 0;
  for (E_Int i = 0; i < nproc; i++)
    rcount[i] = 0;

  for (E_Int i = 0; i < nncells; i++) {
    for (E_Int j = nxcells[i]; j < nxcells[i+1]; j++) {
      E_Int face = NFACE[j];
      if (FT.find(face) == FT.end()) {
        FT[face] = nnfaces++;
        E_Int source = get_proc(face-1, faces_dist, nproc);
        rcount[source]++;
      }
    }
  }
  
  // requesting, roles are inverted
  MPI_Alltoall(&rcount[0], 1, MPI_INT, &scount[0], 1, MPI_INT, MPI_COMM_WORLD);

  rdist[0] = sdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    rdist[i+1] = rdist[i] + rcount[i];
    sdist[i+1] = sdist[i] + scount[i];
  }

  assert(nnfaces == rdist[nproc]);
  std::vector<E_Int> sfaces(sdist[nproc]);
  std::vector<E_Int> rfaces(rdist[nproc]);
  std::vector<E_Int> nxfaces(nnfaces+1);
  nxfaces[0] = 0;

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = rdist[i];

  for (const auto& face : FT) {
    E_Int source = get_proc(face.first-1, faces_dist, nproc);
    rfaces[idx[source]++] = face.first;
  }

  MPI_Alltoallv(&rfaces[0], &rcount[0], &rdist[0], MPI_INT,
                &sfaces[0], &scount[0], &sdist[0], MPI_INT,
                MPI_COMM_WORLD);

  std::vector<E_Int> fstride(sdist[nproc]);
  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];
  
  std::vector<E_Int> sscount(nproc, 0);
  std::vector<E_Int> rrcount(nproc);
  
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &sfaces[sdist[i]];
    for (E_Int j = 0; j < scount[i]; j++) {
      E_Int face = ptr[j]-1;
      face -= faces_dist[rank];
      E_Int stride = xfaces[face+1] - xfaces[face];
      fstride[idx[i]++] = stride;
      sscount[i] += stride;
    }
  }

  MPI_Alltoallv(&fstride[0], &scount[0], &sdist[0], MPI_INT,
                &nxfaces[0]+1, &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  for (E_Int i = 0; i < nnfaces; i++)
    nxfaces[i+1] += nxfaces[i];

  MPI_Alltoall(&sscount[0], 1, MPI_INT, &rrcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> ssdist(nproc+1);
  std::vector<E_Int> rrdist(nproc+1);
  ssdist[0] = rrdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    ssdist[i+1] = ssdist[i] + sscount[i];
    rrdist[i+1] = rrdist[i] + rrcount[i];
  }

  std::vector<E_Int> sNGON(ssdist[nproc]);
  std::vector<E_Int> NGON(rrdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = ssdist[i];

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &sfaces[sdist[i]];
    for (E_Int j = 0; j < scount[i]; j++) {
      E_Int face = ptr[j] - 1 - faces_dist[rank];
      for (E_Int k = xfaces[face]; k < xfaces[face+1]; k++)
        sNGON[idx[i]++] = faces[k];
    }
  }

  MPI_Alltoallv(&sNGON[0], &sscount[0], &ssdist[0], MPI_INT,
                &NGON[0], &rrcount[0], &rrdist[0], MPI_INT,
                MPI_COMM_WORLD);
 
  if (rank == 0)
    printf("NGON OK\n");

  // renumber faces
  FT.clear();
  for (E_Int i = 0; i < nnfaces; i++) {
    E_Int face = rfaces[i];
    FT[face] = i;
  } 
  
  // Hash and request points
  std::unordered_map<E_Int, E_Int> PT;
  E_Int nnpoints = 0;
  std::vector<E_Int> p_rcount(nproc, 0);

  for (E_Int i = 0; i < nnfaces; i++) {
    for (E_Int j = nxfaces[i]; j < nxfaces[i+1]; j++) {
      E_Int point = NGON[j];
      if (PT.find(point) == PT.end()) {
        PT[point] = nnpoints++;
        E_Int source = get_proc(point-1, points_dist, nproc);
        p_rcount[source]++;
      }
    }
  }
  
  // requesting, roles are inverted
  std::vector<E_Int> p_scount(nproc);
  MPI_Alltoall(&p_rcount[0], 1, MPI_INT, &p_scount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> p_rdist(nproc+1);
  std::vector<E_Int> p_sdist(nproc+1);
  p_rdist[0] = p_sdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    p_rdist[i+1] = p_rdist[i] + p_rcount[i];
    p_sdist[i+1] = p_sdist[i] + p_scount[i];
  }

  assert(p_rdist[nproc] == nnpoints);
  std::vector<E_Int> rpoints(p_rdist[nproc]);
  std::vector<E_Int> spoints(p_sdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = p_rdist[i];

  for (const auto& point : PT) {
    E_Int source = get_proc(point.first-1, points_dist, nproc);
    rpoints[idx[source]++] = point.first;
  }
  
  MPI_Alltoallv(&rpoints[0], &p_rcount[0], &p_rdist[0], MPI_INT,
                &spoints[0], &p_scount[0], &p_sdist[0], MPI_INT,
                MPI_COMM_WORLD);

  std::vector<E_Int> sxcount(nproc,0);
  std::vector<E_Int> rxcount(nproc,0);
  std::vector<E_Int> sxdist(nproc+1);
  std::vector<E_Int> rxdist(nproc+1);
  sxdist[0] = rxdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    sxcount[i] = p_scount[i]*3;
    rxcount[i] = p_rcount[i]*3;
    sxdist[i+1] = sxdist[i] + sxcount[i];
    rxdist[i+1] = rxdist[i] + rxcount[i];
  }

  std::vector<E_Float> scrd(sxdist[nproc]);
  std::vector<E_Float> rcrd(rxdist[nproc]);

  E_Int first_point = points_dist[rank];

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sxdist[i];

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &spoints[p_sdist[i]];
    for (E_Int j = 0; j < p_scount[i]; j++) {
      E_Int point = ptr[j]-1;
      point -= first_point;
      scrd[idx[i]++] = X[point];
      scrd[idx[i]++] = Y[point];
      scrd[idx[i]++] = Z[point];
    }
  }

  MPI_Alltoallv(&scrd[0], &sxcount[0], &sxdist[0], MPI_DOUBLE, 
                &rcrd[0], &rxcount[0], &rxdist[0], MPI_DOUBLE,
                MPI_COMM_WORLD);

  if (rank == 0)
    printf("Points OK\n");

  // renumber points
  PT.clear();
  for (E_Int i = 0; i < nnpoints; i++) {
    E_Int point = rpoints[i];
    PT[point] = i;
  }

  FldArrayF local_crd;
  local_crd.malloc(nnpoints, 3);

  for (E_Int i = 0; i < nnpoints; i++) {
    E_Float *px = &rcrd[3*i];
    for (E_Int j = 0; j < 3; j++)
      local_crd(i, j+1) = px[j];
  }

  // construct communication patches
  for (E_Int i = 0; i < nproc; i++)
    rcount[i] = 0;

  for (E_Int i = 0; i < nnfaces; i++) {
      E_Int face = rfaces[i];
      E_Int src = get_proc(face-1, faces_dist, nproc);
      rcount[src]++;
  }

  MPI_Alltoall(&rcount[0], 1, MPI_INT, &scount[0], 1, MPI_INT, MPI_COMM_WORLD);

  rdist[0] = sdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    rdist[i+1] = rdist[i] + rcount[i];
    sdist[i+1] = sdist[i] + scount[i];
  }
  
  for (E_Int i = 0; i < nproc; i++)
    idx[i] = rdist[i];

  for (E_Int i = 0; i < nproc; i++) {
    sscount[i] = 2*scount[i];
    rrcount[i] = 2*rcount[i];
    ssdist[i+1] = ssdist[i] + sscount[i];
    rrdist[i+1] = rrdist[i] + rrcount[i];
  }

  std::vector<E_Int> ssdata(ssdist[nproc]);
  std::vector<E_Int> rrdata(rrdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = ssdist[i];

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &sfaces[sdist[i]];
    for (E_Int j = 0; j < scount[i]; j++) {
      E_Int face = ptr[j];
      assert(rank == get_proc(face-1, faces_dist, nproc));
      face -= 1 + faces_dist[rank];
      E_Int own = A[2*face];
      E_Int nei = A[2*face+1];
      ssdata[idx[i]++] = own;
      ssdata[idx[i]++] = nei;
    }
  }

  for (E_Int i = 0; i < nproc; i++)
    assert(idx[i] == ssdist[i+1]);

  assert(ssdist[nproc] == 2*sdist[nproc]);

  MPI_Alltoallv(&ssdata[0], &sscount[0], &ssdist[0], MPI_INT,
                &rrdata[0], &rrcount[0], &rrdist[0], MPI_INT,
                MPI_COMM_WORLD);

  std::unordered_map<E_Int, E_Int> CT;
  for (E_Int i = 0; i < nncells; i++) {
    CT[rcells[i]] = i;
  }

  std::vector<E_Int> pneis; 

  E_Int nif = 0;

  std::vector<E_Int> rncount(nproc,0);
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &rrdata[rrdist[i]];
    E_Int j = 0;
    for (; j < rrcount[i];) {
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
    assert(j == rrcount[i]);
  }

  std::vector<E_Int> sncount(nproc);
  MPI_Alltoall(&rncount[0], 1, MPI_INT, &sncount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> sndist(nproc+1);
  std::vector<E_Int> rndist(nproc+1);
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
    E_Int *ptr = &rrdata[rrdist[i]];
    E_Int *pf = &rfaces[rdist[i]];
    E_Int j = 0;
    E_Int count = 0;
    for (; j < rrcount[i];) {
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

  MPI_Alltoallv(&rncells[0], &rncount[0], &rndist[0], MPI_INT,
                &sncells[0], &sncount[0], &sndist[0], MPI_INT,
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
  MPI_Alltoallv(&sncells[0], &sncount[0], &sndist[0], MPI_INT,
                &nei2proc[0], &rncount[0], &rndist[0], MPI_INT,
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

  // request signed faces info if need be
  if (sfaces_exist) {
    for (E_Int i = 0; i < nproc; i++)
      rcount[i] = 0;

    // my faces are in rfaces
    for (E_Int i = 0; i < nnfaces; i++) {
      E_Int face = rfaces[i];
      E_Int source = get_proc(face-1, faces_dist, nproc);
      rcount[source]++;
    }

    MPI_Alltoall(&rcount[0], 1, MPI_INT, &scount[0], 1, MPI_INT, MPI_COMM_WORLD);

    sdist[0] = rdist[0] = 0;
    for (E_Int i = 0; i < nproc; i++) {
      rdist[i+1] = rdist[i] + rcount[i];
      sdist[i+1] = sdist[i] + scount[i];
    }

    rdata.resize(rdist[nproc]);
    sdata.resize(sdist[nproc]);

    for (E_Int i = 0; i < nproc; i++)
      idx[i] = rdist[i];

    for (E_Int i = 0; i < nnfaces; i++) {
      E_Int face = rfaces[i];
      E_Int source = get_proc(face-1, faces_dist, nproc);
      rdata[idx[source]++] = face;
    }

    MPI_Alltoallv(&rdata[0], &rcount[0], &rdist[0], MPI_INT,
                  &sdata[0], &scount[0], &sdist[0], MPI_INT,
                  MPI_COMM_WORLD);

    for (E_Int i = 0; i < nproc; i++)
      sscount[i] = 0;

    for (E_Int i = 0; i < nproc; i++) {
      E_Int *ptr = &sdata[sdist[i]];
      for (E_Int j = 0; j < scount[i]; j++) {
        E_Int face = ptr[j];
        if (SF.find(face) != SF.end())
          sscount[i]++;
      }
    }

    std::vector<E_Int> rscount(nproc); // signed rcount
    MPI_Alltoall(&sscount[0], 1, MPI_INT, &rscount[0], 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<E_Int> rsdist(nproc+1);
    ssdist[0] = rsdist[0] = 0;
    for (E_Int i = 0; i < nproc; i++) {
      ssdist[i+1] = ssdist[i] + sscount[i];
      rsdist[i+1] = rsdist[i] + rscount[i];
    }

    std::vector<E_Int> ssdata(ssdist[nproc]);
    std::vector<E_Int> rsdata(rsdist[nproc]);

    for (E_Int i = 0; i < nproc; i++)
      idx[i] = ssdist[i];

    for (E_Int i = 0; i < nproc; i++) {
      E_Int *ptr = &sdata[sdist[i]];
      for (E_Int j = 0; j < scount[i]; j++) {
        E_Int face = ptr[j];
        if (SF.find(face) != SF.end())
          ssdata[idx[i]++] = face;
      }
    }

    MPI_Alltoallv(&ssdata[0], &sscount[0], &ssdist[0], MPI_INT,
                  &rsdata[0], &rscount[0], &rsdist[0], MPI_INT,
                MPI_COMM_WORLD);
  
    // make local PE
    std::vector<E_Int> lPE(2*nnfaces, -1);
    std::vector<E_Int> cPE(nnfaces, 0);
    for (E_Int i = 0; i < nncells; i++) {
      for (E_Int j = nxcells[i]; j < nxcells[i+1]; j++) {
        E_Int face = NFACE[j];
        E_Int lf = FT[face];
        lPE[2*lf + cPE[lf]++] = i;
      }
    }

    // reset signed faces
    for (E_Int i = 0; i < nproc; i++) {
      E_Int *ptr = &rsdata[rsdist[i]];
      for (E_Int j = 0; j < rscount[i]; j++) {
        E_Int face = ptr[j];
        E_Int lf = FT[face];
        
        E_Int own = lPE[2*lf];
        for (E_Int k = nxcells[own]; k < nxcells[own+1]; k++) {
          if (face == NFACE[k]) {
            NFACE[k] = -face;
            break;
          }
        }
        
        E_Int nei = lPE[2*lf+1];
        if (nei == -1) continue;
        for (E_Int k = nxcells[nei]; k < nxcells[nei+1]; k++) {
          if (face == NFACE[k]) {
            NFACE[k] = -face;
            break;
          }
        }
      }
    }

    if (rank == 0)
      printf("Signed faces OK\n");
  }
  
  // TODO (Imad): use Cassiopee's data structures to avoid copying
  const char *varString = "CoordinateX,CoordinateY,CoordinateZ";
 
  /* export array3 / NGON v4 */
  // TODO(Imad): avoid copying
  PyObject* m = K_ARRAY::buildArray3(3, varString, nnpoints, nncells, nnfaces, "NGON",
                               nxfaces[nnfaces], nxcells[nncells], 3,  
                               false, 3);
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  K_ARRAY::getFromArray3(m, f, cn); 
  
  // copie des champs
  for (E_Int n = 0; n < 3; n++)
  {
    E_Float* pt = f->begin(n+1);
    for (E_Int i = 0; i < nnpoints; i++) pt[i] = local_crd(i, n+1);
  }
  // copie connect
  E_Int* ngon = cn->getNGon();
  E_Int* nface = cn->getNFace();
  E_Int* indPG = cn->getIndPG();
  E_Int* indPH = cn->getIndPH();
  for (E_Int i = 0; i <= nnfaces; i++) indPG[i] = nxfaces[i];
  for (E_Int i = 0; i <= nncells; i++) indPH[i] = nxcells[i];
  E_Int* ptr = ngon;
  E_Int start, end;
  for (E_Int i = 0; i < nnfaces; i++)
  {
    start = nxfaces[i];
    end = nxfaces[i+1];
    for (E_Int j = start; j < end; j++) 
    { *ptr = PT[NGON[j]]+1; ptr++; }
  }
  ptr = nface;
  for (E_Int i = 0; i < nncells; i++)
  {
    start = nxcells[i];
    end = nxcells[i+1];
    for (E_Int j = start; j < end; j++) 
    { *ptr = FT[NFACE[j]]+1; ptr++; }
  }

  // Build output Python list
  PyObject* out = PyList_New(0);

  PyList_Append(out, m);  
  Py_DECREF(m);

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
 
    PyArrayObject *f = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_INT);
    E_Int *pf = (E_Int *)PyArray_DATA(f);
    for (size_t j = 0; j < faces.size(); j++)
      pf[j] = faces[j];

    // neis array
    PyArrayObject *n = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_INT);
    E_Int *pn = (E_Int *)PyArray_DATA(n);
    for (size_t j = 0; j < neis.size(); j++)
      pn[j] = neis[j];
  
    PyObject *arr = Py_BuildValue("[lOO]", proc, f, n);
    Py_DECREF(f);
    Py_DECREF(n);
    PyList_Append(comm_data, arr);
  }

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
      res = K_NUMPY::getFromNumpyArray(csol, csols[i], ncells, nfld, true);
      assert(res == 1);
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

  // 9 must be an array of FlowSolutions chunks
  o = PyList_GetItem(l, 8);
  E_Int psize = PyList_Size(o);
  if (psize == 0) {
    PyList_Append(out, PyList_New(0));
  } else {
    E_Float **psols = (E_Float **)XCALLOC(psize, sizeof(E_Float *));

    for (E_Int i = 0; i < psize; i++) {
      PyObject *psol = PyList_GetItem(o, i);
      res = K_NUMPY::getFromNumpyArray(psol, psols[i], npoints, nfld, true);
      assert(res == 1);
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

  // my global cells
  dims[1] = 1;
  dims[0] = (npy_intp)nncells;
  PyArrayObject *mycells = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_INT);
  E_Int *pc = (E_Int *)PyArray_DATA(mycells);
  for (E_Int i = 0; i < nncells; i++)
    pc[i] = rcells[i];

  // my global faces
  dims[0] = (npy_intp)nnfaces;
  PyArrayObject *myfaces = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_INT);
  E_Int *pf = (E_Int *)PyArray_DATA(myfaces);
  for (E_Int i = 0; i < nnfaces; i++)
    pf[i] = rfaces[i];
 
  // my global points
  dims[0] = (npy_intp)nnpoints;
  PyArrayObject *mypoints = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_INT);
  E_Int *pp = (E_Int *)PyArray_DATA(mypoints);
  for (E_Int i = 0; i < nnpoints; i++)
    pp[i] = rpoints[i];

  PyList_Append(out, (PyObject *)mycells); 
  PyList_Append(out, (PyObject *)myfaces); 
  PyList_Append(out, (PyObject *)mypoints); 

  Py_DECREF(mycells); 
  Py_DECREF(myfaces);
  Py_DECREF(mypoints);
  
  return out;
}
