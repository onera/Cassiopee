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

PyObject* K_XCORE::chunk2part(PyObject *self, PyObject *args)
{
  PyObject *array;
  if (!PyArg_ParseTuple(args, "O", &array))
    return NULL;

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
  E_Int *cells_dist = (E_Int *)malloc((nproc+1) * sizeof(E_Int));
  cells_dist[0] = 0;
 
  MPI_Allgather(&ncells, 1, MPI_INT, cells_dist+1, 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++)
    cells_dist[i+1] += cells_dist[i];

  // construct faces distribution
  E_Int *faces_dist = (E_Int *)malloc((nproc+1) * sizeof(E_Int));
  faces_dist[0] = 0;
 
  MPI_Allgather(&nfaces, 1, MPI_INT, faces_dist+1, 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++)
    faces_dist[i+1] += faces_dist[i];

  // construct points distribution
  E_Int *points_dist = (E_Int *)malloc((nproc+1) * sizeof(E_Int));
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
  
  for (E_Int i = 0; i < nproc; i++)
    assert(idx[i] == sdist[i+1]);

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

  for (E_Int i = 0; i < nproc; i++)
    assert(idx[i] == sdist[i+1]);

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
      assert(get_proc(cell, cells_dist, nproc) == rank);
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

  // Send cells to their target proc!
  for (E_Int i = 0; i < nproc; i++)
    scount[i] = 0;

  // start with cell ids
  for (E_Int i = 0; i < ncells; i++) {
    E_Int where = part[i];
    scount[where]++;
  }

  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  // Note(Imad): new ncells may differ from cells_dist[rank+1] - cells_dist[rank]
  E_Int nncells = rdist[nproc];

  std::vector<E_Int> scells(sdist[nproc]);
  std::vector<E_Int> rcells(rdist[nproc]);
  std::vector<E_Int> sstride(sdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];

  for (E_Int i = 0; i < ncells; i++) {
    E_Int where = part[i];
    scells[idx[where]] = i + cells_dist[rank];
    sstride[idx[where]] = xcells[i+1] - xcells[i];
    idx[where]++;
  }

  MPI_Alltoallv(&scells[0], &scount[0], &sdist[0], MPI_INT,
                &rcells[0], &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  std::vector<E_Int> nxcells(nncells+1);

  MPI_Alltoallv(&sstride[0], &scount[0], &sdist[0], MPI_INT,
                &nxcells[0]+1, &rcount[0], &rdist[0], MPI_INT,
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

  // Hash and request faces
  std::vector<E_Int> target_proc;
  std::unordered_map<E_Int, E_Int> FT;
  std::vector<E_Int> l2gf;
  E_Int nnfaces = 0;
  for (E_Int i = 0; i < nproc; i++)
    rcount[i] = 0;

  for (E_Int i = 0; i < nncells; i++) {
    for (E_Int j = nxcells[i]; j < nxcells[i+1]; j++) {
      E_Int face = NFACE[j];
      if (FT.find(face) == FT.end()) {
        l2gf.push_back(face);
        FT[face] = nnfaces++;
        E_Int source = get_proc(face-1, faces_dist, nproc);
        target_proc.push_back(source);
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

  assert(rdist[nproc] == nnfaces);
  std::vector<E_Int> sfaces(sdist[nproc]);
  std::vector<E_Int> rfaces(rdist[nproc]);
  std::vector<E_Int> nxfaces(nnfaces+1);
  nxfaces[0] = 0;

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = rdist[i];

  for (const auto& face : FT) {
    E_Int source = target_proc[face.second];
    rfaces[idx[source]++] = face.first;
  }

  MPI_Alltoallv(&rfaces[0], &rcount[0], &rdist[0], MPI_INT,
                &sfaces[0], &scount[0], &sdist[0], MPI_INT,
                MPI_COMM_WORLD);

  // send face strides
  sstride.resize(sdist[nproc]);
  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &sfaces[sdist[i]];
    for (E_Int j = 0; j < scount[i]; j++) {
      E_Int face = ptr[j];
      assert(rank == get_proc(face-1, faces_dist, nproc));
      face = face - 1 - first_face;
      sstride[idx[i]++] = xfaces[face+1] - xfaces[face];
    }
  }

  std::vector<E_Int> rstride(rdist[nproc]);
  MPI_Alltoallv(&sstride[0], &scount[0], &sdist[0], MPI_INT,
                &rstride[0], &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  for (E_Int i = 0; i < rdist[nproc]; i++) {
    E_Int gf = rfaces[i];
    E_Int lf = FT[gf];
    nxfaces[lf+1] = rstride[lf];
  }

  for (E_Int i = 0; i < nnfaces; i++)
    nxfaces[i+1] += nxfaces[i];

  std::vector<E_Int> sfcount(nproc, 0);
  std::vector<E_Int> rfcount(nproc);
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pf = &sfaces[sdist[i]];
    for (E_Int j = 0; j < scount[i]; j++) {
      E_Int face = pf[j]-1-first_face;
      sfcount[i] += xfaces[face+1] - xfaces[face];
    }
  }

  MPI_Alltoall(&sfcount[0], 1, MPI_INT, &rfcount[0], 1, MPI_INT,
    MPI_COMM_WORLD);

  std::vector<E_Int> sfdist(nproc+1);
  std::vector<E_Int> rfdist(nproc+1);
  sfdist[0] = rfdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    sfdist[i+1] = sfdist[i] + sfcount[i];
    rfdist[i+1] = rfdist[i] + rfcount[i];
  }

  std::vector<E_Int> sNGON(sfdist[nproc]);
  std::vector<E_Int> rNGON(rfdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sfdist[i];

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *pf = &sfaces[sdist[i]];
    for (E_Int j = 0; j < scount[i]; j++) {
      E_Int face = pf[j];
      assert(rank == get_proc(face-1, faces_dist, nproc));
      face = face - 1 - first_face;
      for (E_Int k = xfaces[face]; k < xfaces[face+1]; k++)
        sNGON[idx[i]++] = faces[k];
    }
  }

  MPI_Alltoallv(&sNGON[0], &sfcount[0], &sfdist[0], MPI_INT,
                &rNGON[0], &rfcount[0], &rfdist[0], MPI_INT,
                MPI_COMM_WORLD);

  std::vector<E_Int> NGON(rfdist[nproc]);
  count = 0;
  for (E_Int i = 0; i < nnfaces; i++) {
    E_Int gf = rfaces[i];
    E_Int lf = FT[gf];
    for (E_Int j = nxfaces[lf]; j < nxfaces[lf+1]; j++)
      NGON[j] = rNGON[count++];
  }

  // Hash and request points
  target_proc.clear();
  std::unordered_map<E_Int, E_Int> PT;
  std::vector<E_Int> l2gp;
  E_Int nnpoints = 0;
  for (E_Int i = 0; i < nproc; i++)
    rcount[i] = 0;

  for (E_Int i = 0; i < nnfaces; i++) {
    for (E_Int j = nxfaces[i]; j < nxfaces[i+1]; j++) {
      E_Int point = NGON[j];
      if (PT.find(point) == PT.end()) {
        l2gp.push_back(point);
        PT[point] = nnpoints++;
        E_Int source = get_proc(point-1, points_dist, nproc);
        target_proc.push_back(source);
        rcount[source]++;
      }
    }
  }
  
  // requesting, roles are inverted
  MPI_Alltoall(&rcount[0], 1, MPI_INT, &scount[0], 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++) {
    rdist[i+1] = rdist[i] + rcount[i];
    sdist[i+1] = sdist[i] + scount[i];
  }

  assert(rdist[nproc] == nnpoints);
  std::vector<E_Int> rpoints(rdist[nproc]);
  std::vector<E_Int> spoints(sdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = rdist[i];

  for (const auto& point : PT) {
    E_Int source = target_proc[point.second];
    rpoints[idx[source]++] = point.first;
  }
  
  MPI_Alltoallv(&rpoints[0], &rcount[0], &rdist[0], MPI_INT,
                &spoints[0], &scount[0], &sdist[0], MPI_INT,
                MPI_COMM_WORLD);

  std::vector<E_Int> sxcount(nproc);
  std::vector<E_Int> rxcount(nproc);
  std::vector<E_Int> sxdist(nproc+1);
  std::vector<E_Int> rxdist(nproc+1);
  sxdist[0] = rxdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    sxcount[i] = scount[i]*3;
    rxcount[i] = rcount[i]*3;
    sxdist[i+1] = sxdist[i] + sxcount[i];
    rxdist[i+1] = rxdist[i] + rxcount[i];
  }

  std::vector<E_Float> scrd(sxdist[nproc]);
  std::vector<E_Float> rcrd(rxdist[nproc]);

  E_Int first_point = points_dist[rank];

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sxdist[i];

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &spoints[sdist[i]];
    for (E_Int j = 0; j < scount[i]; j++) {
      E_Int point = ptr[j]-1;
      assert(rank == get_proc(point, points_dist, nproc));
      point -= first_point;
      scrd[idx[i]++] = X[point];
      scrd[idx[i]++] = Y[point];
      scrd[idx[i]++] = Z[point];
    }
  }

  MPI_Alltoallv(&scrd[0], &sxcount[0], &sxdist[0], MPI_DOUBLE, 
                &rcrd[0], &rxcount[0], &rxdist[0], MPI_DOUBLE,
                MPI_COMM_WORLD);

  FldArrayF local_crd;
  local_crd.malloc(nnpoints, 3);

  count = 0;
  for (E_Int i = 0; i < nnpoints; i++) {
    E_Int gp = rpoints[i];
    E_Int lp = PT[gp];
    for (E_Int j = 0; j < 3; j++)
      local_crd(lp, j+1) = rcrd[count++];
  }

  // TODO (Imad): use Cassiopee's data structures to avoid copying
  const char *varString = "CoordinateX,CoordinateY,CoordinateZ";

  FldArrayI cn;
  cn.malloc(2+nnfaces+nxfaces[nnfaces] + 2+nncells+nxcells[nncells], 1);
  cn.setAllValuesAtNull();
  E_Int *ptr = cn.begin(1);
  ptr[0] = nnfaces;
  ptr[1] = nnfaces+nxfaces[nnfaces];
  ptr += 2;
  for (E_Int i = 0; i < nnfaces; i++) {
    E_Int start = nxfaces[i];
    E_Int end = nxfaces[i+1];
    E_Int stride = end - start;
    *ptr++ = stride;
    for (E_Int j = start; j < end; j++) {
      *ptr++ = PT[NGON[j]]+1;
    }
  }

  ptr = cn.begin(1) + 2 + nnfaces + nxfaces[nnfaces];
  ptr[0] = nncells;
  ptr[1] = nncells + nxcells[nncells];
  ptr += 2;
  for (E_Int i = 0; i < nncells; i++) {
    E_Int start = nxcells[i];
    E_Int end = nxcells[i+1];
    E_Int stride = end - start;
    *ptr++ = stride;
    for (E_Int j = start; j < end; j++) {
      *ptr++ = FT[NFACE[j]]+1;
    }
  }

    
  PyObject* m = K_ARRAY::buildArray(local_crd, varString, cn, 8, NULL, false);


  return m;
  return Py_None;
}
