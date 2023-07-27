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

#define ASSERT(expr) \
    if (!(expr)) { \
        printf("\nASSERTION failed on line %d of file %s: " #expr "\n", \
            __LINE__, __FILE__); \
        exit(1); \
    }

#define EXIT \
  do { \
    MPI_Finalize(); \
    exit(0); \
  } while (0);

struct ikv {
  E_Int key;
  E_Int val;

  ikv() : key(-1), val(-1)
  {}

  ikv(E_Int k, E_Int v) : key(k), val(v)
  {}
};

struct rkv {
  E_Float key;
  E_Int val;
};

bool compare_rkv(const rkv &a, const rkv &b)
{
  return (a.key < b.key) || ((a.key == b.key) && (a.val < b.val));
}

bool compare_ikv(const ikv &a, const ikv &b)
{
  return a.key < b.key;
}

bool compare_ikv_2(const ikv &a, const ikv &b)
{
  return (a.key < b.key) || ((a.key == b.key) && (a.val < b.val));
}


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


void sample_sort(E_Int ncells, std::vector<ikv> &elements, E_Int *cells_dist);
void bin_coords(E_Int ncells, E_Float *xyz, E_Int nbins, E_Int *bxyz);

PyObject *redistribute_cells(std::vector<ikv>& elements, E_Int ncells, const E_Int *cells,
  const E_Int *xcells, E_Int nfaces, const E_Int *faces, const E_Int *xfaces, const E_Int *l2gf,
  const E_Int *l2gp, E_Int *cells_dist, E_Int npoints, const FldArrayF& lcrd)
{
  E_Int rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  std::vector<E_Int> scount(nproc, 0);
  for (E_Int i = 0; i < ncells; i++)
    scount[elements[i].key]++;
 
  // we will first send the global ids and their stride

  std::vector<E_Int> rcount(nproc);
  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> sdist(nproc+1);
  std::vector<E_Int> rdist(nproc+1);
  sdist[0] = rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  std::vector<E_Int> scells(sdist[nproc]);
  std::vector<E_Int> rcells(rdist[nproc]);
  std::vector<E_Int> sstride(sdist[nproc]);

  std::vector<E_Int> idx(nproc);
  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];
  
  for (E_Int i = 0; i < ncells; i++) {
    scells[idx[elements[i].key]] = i + cells_dist[nproc];
    sstride[idx[elements[i].key]] = xcells[i+1] - xcells[i];
    idx[elements[i].key]++;
  }

  MPI_Alltoallv(&scells[0], &scount[0], &sdist[0], MPI_INT,
                &rcells[0], &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);
 
  std::vector<E_Int> nxcells(ncells+1);
  MPI_Alltoallv(&sstride[0], &scount[0], &sdist[0], MPI_INT,
                &nxcells[0]+1, &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  nxcells[0] = 0;
  for (E_Int i = 0; i < ncells; i++) {
    nxcells[i+1] += nxcells[i];
  }

  // now send the NFACE
  for (E_Int i = 0; i < nproc; i++)
    scount[i] = 0;

  for (E_Int i = 0; i < ncells; i++) {
    scount[elements[i].key] += xcells[i+1] - xcells[i];
  }

  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];

  std::vector<E_Int> sNFACE(sdist[nproc]);

  for (E_Int i = 0; i < ncells; i++) {
    for (E_Int j = xcells[i]; j < xcells[i+1]; j++)
      sNFACE[idx[elements[i].key]++] = l2gf[cells[j]];
  }

  std::vector<E_Int> NFACE(rdist[nproc]);
  MPI_Alltoallv(&sNFACE[0], &scount[0], &sdist[0], MPI_INT,
                &NFACE[0], &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  // make parent element
  std::vector<E_Int> PE(2*nfaces, -1);
  std::vector<E_Int> visited(nfaces, 0);
  for (E_Int i = 0; i < ncells; i++) {
    for (E_Int j = xcells[i]; j < xcells[i+1]; j++) {
      E_Int face = cells[j];
      PE[2*face + visited[face]++] = i;
      assert(visited[face] <= 2);
    }
  }

  // first send the global face ids and strides
  for (E_Int i = 0; i < nproc; i++)
    scount[i] = 0;

  for (E_Int i = 0; i < nfaces; i++) {
    for (E_Int j = 0; j < 2; j++) {
      E_Int elem = PE[2*i+j];
      E_Int dest = elements[elem].key;
      scount[dest]++;
    }
  }

  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }
  
  std::vector<E_Int> sfaces(sdist[nproc]);
  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];

  sstride.resize(sdist[nproc]);

  for (E_Int i = 0; i < nfaces; i++) {
    for (E_Int j = 0; j < 2; j++) {
      E_Int elem = PE[2*i+j];
      if (elem < 0) continue;
      E_Int dest = elements[elem].key;
      sfaces[idx[dest]] = l2gf[i];
      sstride[idx[dest]] = xfaces[i+1] - xfaces[i];
      idx[dest]++;
    }
  }
      
  std::vector<E_Int> rfaces(rdist[nproc]);
  MPI_Alltoallv(&sfaces[0], &scount[0], &sdist[0], MPI_INT,
                &rfaces[0], &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  E_Int nnfaces = rdist[nproc]; // number of new faces
  std::vector<E_Int> nxfaces(nnfaces+1);
  MPI_Alltoallv(&sstride[0], &scount[0], &sdist[0], MPI_INT,
                &nxfaces[0]+1, &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  nxfaces[0] = 0;
  for (E_Int i = 0; i < nnfaces; i++) {
    assert(nxfaces[i+1] == 4);
    nxfaces[i+1] += nxfaces[i];
  }

  // send ngon
  for (E_Int i = 0; i < nproc; i++)
    scount[i] = 0;

  for (E_Int i = 0; i < nfaces; i++) {
    for (E_Int j = 0; j < 2; j++) {
      E_Int elem = PE[2*i+j];
      if (elem < 0) continue;
      E_Int dest = elements[elem].key;
      scount[dest] += xfaces[i+1] - xfaces[i];
    }
  }

  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];
  
  std::vector<E_Int> sNGON(sdist[nproc]);

  for (E_Int i = 0; i < nfaces; i++) {
    for (E_Int j = 0; j < 2; j++) {
      E_Int elem = PE[2*i+j];
      if (elem < 0) continue;
      E_Int dest = elements[elem].key;
      for (E_Int k = xfaces[i]; k < xfaces[i+1]; k++)
        sNGON[idx[dest]++] = l2gp[faces[k]];
    }
  }

  std::vector<E_Int> NGON(rdist[nproc]);
  MPI_Alltoallv(&sNGON[0], &scount[0], &sdist[0], MPI_INT,
                &NGON[0], &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  
  // make pointCells
  std::vector<std::vector<E_Int>> pointCells(npoints);
  for (E_Int i = 0; i < nfaces; i++) {
    E_Int start = xfaces[i];
    E_Int end = xfaces[i+1];

    E_Int own = PE[2*i];
    E_Int nei = PE[2*i+1];

    for (E_Int j = start; j < end; j++) {
      pointCells[faces[j]].push_back(own);
      if (nei >= 0)
        pointCells[faces[j]].push_back(nei);
    }
  }

  // get rid of duplicate cells
  for (E_Int point = 0; point < npoints; point++) {
    std::sort(pointCells[point].begin(), pointCells[point].end());
    E_Int n = 1;
    for (size_t i = 1; i < pointCells[point].size(); i++) {
      if (pointCells[point][i-1] != pointCells[point][i])
        pointCells[point][n++] = pointCells[point][i];
    }

    pointCells[point].resize(n);
  }

  for (E_Int i = 0; i < nproc; i++)
    scount[i] = 0;

  for (E_Int i = 0; i < npoints; i++) {
    for (size_t j = 0; j < pointCells[i].size(); j++) {
      E_Int cell = pointCells[i][j];
      E_Int dest = elements[cell].key;
      scount[dest]++;
    }
  }

  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];

  std::vector<E_Int> spoints(sdist[nproc]);

  for (E_Int i = 0; i < npoints; i++) {
    for (size_t j = 0; j < pointCells[i].size(); j++) {
      E_Int cell = pointCells[i][j];
      E_Int dest = elements[cell].key;
      spoints[idx[dest]++] = l2gp[i];
    }
  }

  E_Int nnpoints = rdist[nproc];
  std::vector<E_Int> rpoints(rdist[nproc]);

  MPI_Alltoallv(&spoints[0], &scount[0], &sdist[0], MPI_INT,
                &rpoints[0], &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  // lastly, send coordinates
  for (E_Int i = 0; i < nproc; i++) {
    scount[i] *= 3; // xyz
    rcount[i] *= 3;
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];

  std::vector<E_Float> scrd(sdist[nproc]);

  for (E_Int i = 0; i < npoints; i++) {
    for (size_t j = 0; j < pointCells[i].size(); j++) {
      E_Int cell = pointCells[i][j];
      E_Int dest = elements[cell].key;
      for (E_Int k = 0; k < 3; k++)
        scrd[idx[dest]++] = lcrd(i, k+1);
    }
  }

  std::vector<E_Float> crd(rdist[nproc]);
  MPI_Alltoallv(&scrd[0], &scount[0], &sdist[0], MPI_DOUBLE,
                &crd[0], &rcount[0], &rdist[0], MPI_DOUBLE,
                MPI_COMM_WORLD);

  // hash points and get rid of dups
  E_Int np = 0;
  E_Int count = 0;
  std::unordered_map<E_Int, E_Int> PT;
  for (E_Int i = 0; i < nnpoints; i++) {
    if (PT.find(rpoints[i]) == PT.end()) {
      PT[rpoints[i]] = np;

      // point gid
      rpoints[np++] = rpoints[i];

      // xyz
      E_Float *px = &crd[3*i];
      for (E_Int j = 0; j < 3; j++)
        crd[count++] = px[j];
    }
  }

  assert(count == 3*np);
  
  nnpoints = np;

  FldArrayF local_crd;
  local_crd.malloc(nnpoints, 3);
  for (E_Int i = 0; i < nnpoints; i++) {
    E_Float *px = &crd[3*i];
    for (E_Int k = 0; k < 3; k++) {
      local_crd(i, k+1) = px[k];
    }
  }

  // hash faces and get rid of dups
  E_Int nf = 0;
  count = 0;
  std::unordered_map<E_Int, E_Int> FT;
  for (E_Int i = 0; i < nnfaces; i++) {
    if (FT.find(rfaces[i]) == FT.end()) {
      FT[rfaces[i]] = nf;

      // face gid
      rfaces[nf] = rfaces[i];

      // face stride
      E_Int stride = nxfaces[i+1] - nxfaces[i];
      nxfaces[nf+1] = nxfaces[nf] + stride;

      nf++;

      // NGON
      for (E_Int j = nxfaces[i]; j < nxfaces[i+1]; j++)
        NGON[count++] = NGON[j];
    }
  }

  // rcells
  // nxcells
  // NFACE
  // rfaces
  // nxfaces
  // NGON
  // rpoints
  // crd
  const char *varString = "CoordinateX,CoordinateY,CoordinateZ";

  FldArrayI cn;
  cn.malloc(2+nnfaces+nxfaces[nnfaces] + 2+ncells+nxcells[ncells], 1);
  cn.setAllValuesAtNull();
  E_Int *ptr = cn.begin(1);
  ptr[0] = nnfaces;
  ptr[1] = nnfaces+nxfaces[nnfaces];
  ptr += 2;
  for (E_Int i = 0; i < nnfaces; i++) {
    E_Int stride = nxfaces[i+1] - nxfaces[i];
    assert(stride == 4);
    *ptr++ = stride;
    for (E_Int j = nxfaces[i]; j < nxfaces[i+1]; j++) {
      assert(PT.find(NGON[j]) != PT.end());
      *ptr++ = PT[NGON[j]]+1;
    }
  }

  ptr = cn.begin(1) + nnfaces + nxfaces[nnfaces] + 2;
  ptr[0] = ncells;
  ptr[1] = ncells + nxcells[ncells];
  ptr += 2;
  for (E_Int i = 0; i < ncells; i++) {
    E_Int start = nxcells[i];
    E_Int end = nxcells[i+1];
    E_Int stride = end - start;
    assert(stride == 6);
    *ptr++ = stride;
    for (E_Int j = start; j < end; j++) {
      *ptr++ = FT[NFACE[j]]+1;
    }
  }

  PyObject* m = K_ARRAY::buildArray(local_crd, varString, cn, 8, NULL, false);

  return m;
}

PyObject* K_XCORE::chunk2part(PyObject *self, PyObject *args)
{
  PyObject *array;
  if (!PyArg_ParseTuple(args, "O", &array))
    return NULL;

  E_Int nzones = PyList_Size(array);
  assert(nzones == 1);

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
    
  // 6 must be nface so chunk
  o = PyList_GetItem(l, 6);
  res = K_NUMPY::getFromNumpyArray(o, xcells, ncells, nfld, true);
  assert(res == 1);
    
  ncells--;
  nfaces--;

  // construct cells distribution
  E_Int *cells_dist = (E_Int *)malloc((ncells+1) * sizeof(E_Int));
  cells_dist[0] = 0;
 
  MPI_Allgather(&ncells, 1, MPI_INT, cells_dist+1, 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++)
    cells_dist[i+1] += cells_dist[i];

  // construct faces distribution
  E_Int *faces_dist = (E_Int *)malloc((nfaces+1) * sizeof(E_Int));
  faces_dist[0] = 0;
 
  MPI_Allgather(&nfaces, 1, MPI_INT, faces_dist+1, 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++)
    faces_dist[i+1] += faces_dist[i];

  // construct points distribution
  E_Int *points_dist = (E_Int *)malloc((npoints+1) * sizeof(E_Int));
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
        assert(c[face] <= 2);
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
  // cells_dist should be recomputed?
  E_Int nncells = rdist[nproc];

  std::vector<E_Int> scells(sdist[nproc]);
  assert(nncells > 0);
  E_Int *rcells = (E_Int *)malloc(nncells * sizeof(E_Int));
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
                rcells, &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  std::vector<E_Int> nxcells(nncells+1);

  MPI_Alltoallv(&sstride[0], &scount[0], &sdist[0], MPI_INT,
                &nxcells[0]+1, &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  nxcells[0] = 0;
  for (E_Int i = 0; i < nncells; i++)
    nxcells[i+1] += nxcells[i];

  // construct new cells distribution
  std::vector<E_Int> ncells_dist(nproc);
  ncells_dist[0] = 0;
  MPI_Allgather(&nncells, 1, MPI_INT, &ncells_dist[0] + 1, 1, MPI_INT, MPI_COMM_WORLD);
  for (E_Int i = 0; i < nproc; i++)
    ncells_dist[i+1] += ncells_dist[i];

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

  for (E_Int i = 0; i < nproc; i++) {
    rdist[i+1] = rdist[i] + rcount[i];
    sdist[i+1] = sdist[i] + scount[i];
  }

  assert(rdist[nproc] == nnfaces);
  std::vector<E_Int> rfaces(rdist[nproc]);
  std::vector<E_Int> nxfaces(nnfaces+1);
  nxfaces[0] = 0;
  sdata.resize(sdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = rdist[i];

  for (const auto& face : FT) {
    E_Int source = target_proc[face.second];
    rfaces[idx[source]++] = face.first;
  }
  
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &rfaces[rdist[i]];
    for (E_Int j = 0; j < rcount[i]; j++)
      assert(i == get_proc(ptr[j]-1, faces_dist, nproc));
  }

  MPI_Alltoallv(&rfaces[0], &rcount[0], &rdist[0], MPI_INT,
                &sdata[0], &scount[0], &sdist[0], MPI_INT,
                MPI_COMM_WORLD);

  // send face strides
  sstride.resize(sdist[nproc]);
  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];

  std::vector<E_Int> sscount(nproc,0);
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &sdata[sdist[i]];
    for (E_Int j = 0; j < scount[i]; j++) {
      E_Int face = ptr[j];
      assert(rank == get_proc(face-1, faces_dist, nproc));
      face = face - 1 - first_face;
      sstride[idx[i]++] = xfaces[face+1] - xfaces[face];
      sscount[i] += xfaces[face+1] - xfaces[face];
    }
  }


  MPI_Alltoallv(&sstride[0], &scount[0], &sdist[0], MPI_INT,
                &nxfaces[0]+1, &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  for (E_Int i = 0; i < nnfaces; i++)
    nxfaces[i+1] += nxfaces[i];

  
  std::vector<E_Int> rrcount(nproc);
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
  assert(rrdist[nproc] == nxfaces[nnfaces]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = ssdist[i];

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &sdata[sdist[i]];
    for (E_Int j = 0; j < scount[i]; j++) {
      E_Int face = ptr[j];
      assert(rank == get_proc(face-1, faces_dist, nproc));
      face = face - 1 - first_face;
      for (E_Int k = xfaces[face]; k < xfaces[face+1]; k++)
        sNGON[idx[i]++] = faces[k];
    }
  }

  MPI_Alltoallv(&sNGON[0], &sscount[0], &ssdist[0], MPI_INT,
                &NGON[0], &rrcount[0], &rrdist[0], MPI_INT,
                MPI_COMM_WORLD);

  if (rank == 0) {
    for (E_Int i = 0; i < nnfaces; i++) {
      printf("%d -> ", rfaces[i]);
      for (E_Int j = nxfaces[i]; j < nxfaces[i+1]; j++)
        printf("%d ", NGON[j]);
      puts("");
    }

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
  sdata.resize(sdist[nproc]);

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = rdist[i];

  for (const auto& point : PT) {
    E_Int source = target_proc[point.second];
    rpoints[idx[source]++] = point.first;
  }
  
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &rpoints[rdist[i]];
    for (E_Int j = 0; j < rcount[i]; j++)
      assert(i == get_proc(ptr[j]-1, points_dist, nproc));
  }

  MPI_Alltoallv(&rpoints[0], &rcount[0], &rdist[0], MPI_INT,
                &sdata[0], &scount[0], &sdist[0], MPI_INT,
                MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++) {
    sscount[i] = scount[i]*3;
    rrcount[i] = rcount[i]*3;
    ssdist[i+1] = ssdist[i] + sscount[i];
    rrdist[i+1] = rrdist[i] + rrcount[i];
  }

  E_Float *sxyz = (E_Float *)malloc(ssdist[nproc] * sizeof(E_Float));
  std::vector<E_Float> crd(rrdist[nproc]);

  E_Int first_point = points_dist[rank];

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = ssdist[i];

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &sdata[sdist[i]];
    for (E_Int j = 0; j < scount[i]; j++) {
      E_Int point = ptr[j]-1;
      assert(rank == get_proc(point, points_dist, nproc));
      point -= first_point;
      sxyz[idx[i]++] = X[point];
      sxyz[idx[i]++] = Y[point];
      sxyz[idx[i]++] = Z[point];
    }
  }

  MPI_Alltoallv(&sxyz[0], &sscount[0], &ssdist[0], MPI_DOUBLE, 
                &crd[0], &rrcount[0], &rrdist[0], MPI_DOUBLE,
                MPI_COMM_WORLD);

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
    assert(stride == 4);
    *ptr++ = stride;
    for (E_Int j = start; j < end; j++) {
      assert(PT.find(NGON[j]) != PT.end());
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
    assert(stride == 6);
    *ptr++ = stride;
    for (E_Int j = start; j < end; j++) {
      assert(FT.find(NFACE[j]) != FT.end());
      *ptr++ = FT[NFACE[j]]+1;
    }
  }

  FldArrayF local_crd;
  local_crd.malloc(nnpoints, 3);

  for (E_Int i = 0; i < nnpoints; i++) {
    E_Int gp = rpoints[i];
    E_Int lp = PT[gp];
    E_Float *xyz = &crd[0] + 3*i;
    for (E_Int j = 0; j < 3; j++)
      local_crd(lp, j+1) = xyz[j];

  }
  
  PyObject* m = K_ARRAY::buildArray(local_crd, varString, cn, 8, NULL, false);


  return m;
  return Py_None;
}


/*

//============================================================================
// IN: Chunk of NGON2
//============================================================================
PyObject* K_XCORE::chunk2part(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;
  
  E_Int nzones = PyList_Size(array);
  assert(nzones == 1); // Note(Imad): for now...

  PyObject* o; PyObject* l;
  E_Int nfld;
  E_Float *X, *Y, *Z;
  E_Int npoints, ncells, nfaces;
  E_Int faces_size, cells_size;
  E_Int *faces, *cells, *xfaces, *xcells;

  E_Int rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  E_Int res;
  
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
    
  // PE a venir...

  // ..
  ncells--;
  nfaces--;

  // construct cells distribution
  E_Int *cells_dist = (E_Int *)malloc((ncells+1) * sizeof(E_Int));
  cells_dist[0] = 0;
 
  MPI_Allgather(&ncells, 1, MPI_INT, cells_dist+1, 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++)
    cells_dist[i+1] += cells_dist[i];

  // construct faces distribution
  E_Int *faces_dist = (E_Int *)malloc((nfaces+1) * sizeof(E_Int));
  faces_dist[0] = 0;
 
  MPI_Allgather(&nfaces, 1, MPI_INT, faces_dist+1, 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++)
    faces_dist[i+1] += faces_dist[i];

  // construct points distribution
  E_Int *points_dist = (E_Int *)malloc((npoints+1) * sizeof(E_Int));
  points_dist[0] = 0;
 
  MPI_Allgather(&npoints, 1, MPI_INT, points_dist+1, 1, MPI_INT, MPI_COMM_WORLD);

  for (E_Int i = 0; i < nproc; i++)
      points_dist[i+1] += points_dist[i];


  E_Int *scount = (E_Int *)calloc(nproc, sizeof(E_Int));
  E_Int *rcount = (E_Int *)calloc(nproc, sizeof(E_Int));

  // shift xcells, xfaces and xpoints to start from zero
  E_Int cell_shift = xcells[0];
  for (E_Int i = 0; i < ncells+1; i++)
    xcells[i] -= cell_shift;

  E_Int face_shift = xfaces[0];
  for (E_Int i = 0; i < nfaces+1; i++)
    xfaces[i] -= face_shift;


  // count how many faces are requested from each proc
  std::vector<std::set<E_Int>> face_tables(nproc);

  for (E_Int i = 0; i < ncells; i++) {
    E_Int start = xcells[i];
    E_Int end = xcells[i+1];

    for (E_Int j = start; j < end; j++) {
      E_Int face = cells[j];
      E_Int source = get_proc(face-1, faces_dist, nproc); 
      auto search = face_tables[source].find(face);

      if (search == face_tables[source].end()) {
        face_tables[source].insert(face);
        rcount[source]++;
      }
    }
  }

  MPI_Alltoall(rcount, 1, MPI_INT, scount, 1, MPI_INT, MPI_COMM_WORLD);

  E_Int *sdist = (E_Int *)malloc((nproc+1) * sizeof(E_Int));
  E_Int *rdist = (E_Int *)malloc((nproc+1) * sizeof(E_Int));
  
  sdist[0] = rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  // faces to be sent/received
  E_Int *sdata = (E_Int *)malloc(sdist[nproc] * sizeof(E_Int));
  E_Int *rdata = (E_Int *)malloc(rdist[nproc] * sizeof(E_Int));

  E_Int *idx = (E_Int *)malloc(nproc * sizeof(E_Int));
  for (E_Int i = 0; i < nproc; i++)
    idx[i] = rdist[i];

  for (E_Int i = 0; i < nproc; i++)
    face_tables[i].clear();

  for (E_Int i = 0; i < ncells; i++) {
    E_Int start = xcells[i];
    E_Int end = xcells[i+1];
    for (E_Int j = start; j < end; j++) {
      E_Int face = cells[j];
      E_Int source = get_proc(face-1, faces_dist, nproc);
      auto search = face_tables[source].find(face);
      if (search == face_tables[source].end()) {
        face_tables[source].insert(face);
        rdata[idx[source]++] = face;
      }
    }
  }

  for (E_Int i = 0; i < nproc; i++)
    assert(idx[i] == rdist[i+1]);

  for (E_Int i = 0; i < nproc; i++)
    face_tables[i].clear();

  // sdata contains the faces to be sent to each proc
  MPI_Alltoallv(rdata, rcount, rdist, MPI_INT,
                sdata, scount, sdist, MPI_INT,
                MPI_COMM_WORLD);

  // face data (stride + points)
  E_Int *fscount = (E_Int *)calloc(nproc, sizeof(E_Int));
  
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &sdata[sdist[i]];
    
    for (E_Int j = 0; j < scount[i]; j++) {
      E_Int face = ptr[j]-1;
      assert(get_proc(face, faces_dist, nproc) == rank);
      E_Int pos = face - faces_dist[rank];
      E_Int start = xfaces[pos];
      E_Int end = xfaces[pos+1];
      E_Int stride = end - start;
      fscount[i] += 1 + stride; // stride + points
    }
  }

  E_Int *frcount = (E_Int *)malloc(nproc * sizeof(E_Int));
  MPI_Alltoall(fscount, 1, MPI_INT, frcount, 1, MPI_INT, MPI_COMM_WORLD);

  E_Int *fsdist = (E_Int *)malloc((nproc+1) * sizeof(E_Int));
  E_Int *frdist = (E_Int *)malloc((nproc+1) * sizeof(E_Int));
  
  fsdist[0] = frdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    fsdist[i+1] = fsdist[i] + fscount[i];
    frdist[i+1] = frdist[i] + frcount[i];
  }

  E_Int *fsdata = (E_Int *)malloc(fsdist[nproc] * sizeof(E_Int));
  E_Int *frdata = (E_Int *)malloc(frdist[nproc] * sizeof(E_Int));

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = fsdist[i];

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &sdata[sdist[i]];

    for (E_Int j = 0; j < scount[i]; j++) {
      E_Int face = ptr[j]-1;
      E_Int pos = face - faces_dist[rank];
      E_Int start = xfaces[pos];
      E_Int end = xfaces[pos+1];
      E_Int stride = end - start;

      fsdata[idx[i]++] = stride;
      
      for (E_Int k = start; k < end; k++)
        fsdata[idx[i]++] = faces[k];
    }

    assert(idx[i] == fsdist[i+1]);
  }

  MPI_Alltoallv(fsdata, fscount, fsdist, MPI_INT,
                frdata, frcount, frdist, MPI_INT,
                MPI_COMM_WORLD);
 
  // frdata contains the stride and points of the faces requested in rdata

  std::unordered_map<E_Int, E_Int> face_table;
  E_Int n_local_faces = 0;
  E_Int *l2gf = (E_Int *)malloc(rdist[nproc] * sizeof(E_Int));
  // hash global faces -> local faces
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *req_faces = &rdata[rdist[i]];

    for (E_Int j = 0; j < rcount[i]; j++) {
      E_Int face = req_faces[j];
      // cannot encounter a face twice
      assert(face_table.find(face) == face_table.end());
      face_table[face] = n_local_faces;
      l2gf[n_local_faces++] = face;
    }
  }
  assert(n_local_faces == rdist[nproc]);

  // request point coordinates
  std::vector<std::set<E_Int>> point_tables(nproc);
  E_Int *prcount = (E_Int *)calloc(nproc, sizeof(E_Int));

  // loop over my local faces (rdata)
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *r_face_data = &frdata[frdist[i]];
    E_Int count = 0;

    for (E_Int j = 0; j < rcount[i]; j++) {
      E_Int stride = r_face_data[count++];
      assert(stride == 4);
      for (E_Int k = 0; k < stride; k++) {
        E_Int point = r_face_data[count++];
        E_Int source = get_proc(point - 1, points_dist, nproc);
        auto search = point_tables[source].find(point);
        if (search == point_tables[source].end()) {
          point_tables[source].insert(point);
          prcount[source]++;
        }
      }
    }
  }

  E_Int *prdist = (E_Int *)malloc((nproc+1) * sizeof(E_Int));
  prdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++)
    prdist[i+1] = prdist[i] + prcount[i];

  for (E_Int i = 0; i < nproc; i++)
    point_tables[i].clear();

  for (E_Int i = 0; i < nproc; i++)
    idx[i] = prdist[i];

  E_Int *prdata = (E_Int *)malloc(prdist[nproc]* sizeof(E_Int));

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *r_face_data = &frdata[frdist[i]];
    E_Int count = 0;

    for (E_Int j = 0; j < rcount[i]; j++) {
      E_Int stride = r_face_data[count++];
      for (E_Int k = 0; k < stride; k++) {
        E_Int point = r_face_data[count++];
        E_Int source = get_proc(point - 1, points_dist, nproc);
        auto search = point_tables[source].find(point);
        if (search == point_tables[source].end()) {
          point_tables[source].insert(point);
          prdata[idx[source]++] = point;
        }
      }
    }
  }

  for (E_Int i = 0; i < nproc; i++)
    point_tables[i].clear();

  for (E_Int i = 0; i < nproc; i++)
    assert(idx[i] == prdist[i+1]);

  E_Int *pscount = (E_Int *)malloc(nproc * sizeof(E_Int));
  MPI_Alltoall(prcount, 1, MPI_INT, pscount, 1, MPI_INT, MPI_COMM_WORLD);
  
  E_Int *psdist = (E_Int *)malloc((nproc+1) * sizeof(E_Int));
  psdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++)
    psdist[i+1] = psdist[i] + pscount[i];

  E_Int *psdata = (E_Int *)malloc(psdist[nproc] * sizeof(E_Int));
  MPI_Alltoallv(prdata, prcount, prdist, MPI_INT,
                psdata, pscount, psdist, MPI_INT,
                MPI_COMM_WORLD);

  // psdata contains the points that need to be sent to every proc
  // all points within psdata should belong to me

  E_Int *xsdist = (E_Int *)malloc((nproc+1) * sizeof(E_Int));
  E_Int *xrdist = (E_Int *)malloc((nproc+1) * sizeof(E_Int));
  E_Int *xscount = (E_Int *)malloc(nproc * sizeof(E_Int));
  E_Int *xrcount = (E_Int *)malloc(nproc * sizeof(E_Int));
  
  xsdist[0] = xrdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    xscount[i] = 3*pscount[i];
    xrcount[i] = 3*prcount[i];
    xsdist[i+1] = xsdist[i] + xscount[i];
    xrdist[i+1] = xrdist[i] + xrcount[i];
  }

  E_Float *xsdata = (E_Float *)malloc(xsdist[nproc] * sizeof(E_Float));
  E_Float *xrdata = (E_Float *)malloc(xrdist[nproc] * sizeof(E_Float));
  
  for (E_Int i = 0; i < nproc; i++)
    idx[i] = xsdist[i];

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &psdata[psdist[i]];
    for (E_Int j = 0; j < pscount[i]; j++) {
      E_Int point = ptr[j];
      point--;
      assert(get_proc(point, points_dist, nproc) == rank);
      point -= points_dist[rank];
      xsdata[idx[i]++] = X[point];
      xsdata[idx[i]++] = Y[point];
      xsdata[idx[i]++] = Z[point];
    }
  }

  for (E_Int i = 0; i < nproc; i++)
    assert(idx[i] == xsdist[i+1]);

  // xrdata contains xyz coordinates of the points requested in prdata
  MPI_Alltoallv(xsdata, xscount, xsdist, MPI_DOUBLE,
                xrdata, xrcount, xrdist, MPI_DOUBLE,
                MPI_COMM_WORLD);


  // hash received points
  std::unordered_map<E_Int, E_Int> point_table;
  E_Int n_local_points = prdist[nproc];
  FldArrayF local_crd;
  local_crd.malloc(n_local_points, 3);

  E_Int *l2gp = (E_Int *)malloc(n_local_points * sizeof(E_Int));
  E_Int count = 0;
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &prdata[prdist[i]];
    E_Float *X = &xrdata[xrdist[i]];
    for (E_Int j = 0; j < prcount[i]; j++) {
      E_Int point = ptr[j];
      
      // cannot encounter a point twice
      assert(point_table.find(point) == point_table.end());
      point_table[point] = count;
      l2gp[count] = point;

      // fill local coordinates
      E_Float *pX = &X[3*j];   
      local_crd(count, 1) = pX[0];
      local_crd(count, 2) = pX[1];
      local_crd(count, 3) = pX[2];

      count++;
    }
  }
  assert(count == n_local_points);

  // update cells connectivity with local faces
  for (E_Int i = 0; i < ncells; i++) {
    E_Int start = xcells[i];
    E_Int end = xcells[i+1];
    for (E_Int j = start; j < end; j++) {
      assert(face_table.find(cells[j]) != face_table.end());
      cells[j] = face_table[cells[j]];
    }
  }

  // local faces are the ones requested in rdata
  // stride and points are stored in frdata
  // replace points with their local id
  count = 0;
  for (E_Int i = 0; i < n_local_faces; i++) {
    E_Int stride = frdata[count++];
    assert(stride == 4);
    for (E_Int j = 0; j < stride; j++) {
      assert(point_table.find(frdata[count]) != point_table.end());
      frdata[count] = point_table[frdata[count]];
      count++;
    }
  }

  // compute xfaces and update faces
  E_Int *p = frdata;
  xfaces = (E_Int *)realloc(xfaces, (n_local_faces+1) * sizeof(E_Int));
  faces = (E_Int *)realloc(faces, (frdist[nproc]-n_local_faces)*sizeof(E_Int));
  xfaces[0] = 0;
  count = 0;
  for (E_Int i = 0; i < n_local_faces; i++) {
    E_Int stride = *p;
    assert(stride == 4);
    xfaces[i+1] = xfaces[i] + stride;
    p++;
    for (E_Int j = 0; j < stride; j++)
      faces[count++] = *p++;
  }

  assert(xfaces[n_local_faces] == frdist[nproc]-n_local_faces);

  // compute cell centers
  E_Float *fc = (E_Float *)calloc(n_local_faces*3, sizeof(E_Float));
  E_Float *cc = (E_Float *)calloc(ncells*3, sizeof(E_Float));

  count = 0;
  for (E_Int i = 0; i < n_local_faces; i++) {
    E_Float *fx = &fc[3*i];
    E_Int stride = frdata[count++];
    for (E_Int j = 0; j < stride; j++) {
      E_Int point = frdata[count++];
      fx[0] += local_crd(point, 1);
      fx[1] += local_crd(point, 2);
      fx[2] += local_crd(point, 3);
    }
    E_Float coeff = 1./stride;
    for (E_Int k = 0; k < 3; k++)
      fx[k] *= coeff;
  }

  for (E_Int i = 0; i < ncells; i++) {
    E_Float *cx = &cc[3*i];
    E_Int start = xcells[i];
    E_Int end = xcells[i+1];
    for (E_Int j = start; j < end; j++) {
      E_Int face = cells[j];
      E_Float *fx = &fc[3*face];
      for (E_Int k = 0; k < 3; k++)
        cx[k] += fx[k];
    }
    E_Float coeff = 1./(end - start);
    for (E_Int k = 0; k < 3; k++)
      cx[k] *= coeff;
  }

  // redistribute cells into contiguous parts
  const E_Int nbits = 9;
  const E_Int nbins = 1<<nbits;
  E_Int *bxyz = (E_Int *)malloc(3*ncells * sizeof(E_Int));
  memset(bxyz, -1, 3*ncells*sizeof(E_Int));
  bin_coords(ncells, cc, nbins, bxyz);

  std::vector<ikv> elements(ncells);
  E_Int first_cell = cells_dist[rank];
  for (E_Int i = 0; i < ncells; i++) {
    E_Int icoord = 0;
    // interleave the bits, starting from the most significant bit
    for (E_Int j = nbits-1; j >= 0; j--) {
      for (E_Int k = 0; k < 3; k++) {
        icoord = (icoord<<1) + (bxyz[i*3+k]&(1<<j) ? 1 : 0);
      }
    }

    elements[i].key = icoord;
    elements[i].val = i + first_cell;
  }

  sample_sort(ncells, elements, cells_dist);

  PyObject *m = redistribute_cells(elements, ncells, cells, xcells, n_local_faces, faces, xfaces, l2gf, l2gp, cells_dist, n_local_points, local_crd);


  // TODO (Imad): use Cassiopee's data structures to avoid copying
  const char *varString = "CoordinateX,CoordinateY,CoordinateZ";

  FldArrayI cn;
  cn.malloc(2+frdist[nproc] + 2+ncells+xcells[ncells], 1);
  cn.setAllValuesAtNull();
  E_Int *ptr = cn.begin(1);
  ptr[0] = n_local_faces;
  ptr[1] = frdist[nproc];
  ptr += 2;
  count = 0;
  for (E_Int i = 0; i < n_local_faces; i++) {
    E_Int stride = frdata[count++];
    assert(stride == 4);
    *ptr++ = stride;
    for (E_Int j = 0; j < stride; j++) {
      *ptr++ = frdata[count++]+1;
    }
  }

  ptr = cn.begin(1) + frdist[nproc] + 2;
  ptr[0] = ncells;
  ptr[1] = ncells + xcells[ncells];
  ptr += 2;
  for (E_Int i = 0; i < ncells; i++) {
    E_Int start = xcells[i];
    E_Int end = xcells[i+1];
    E_Int stride = end - start;
    assert(stride == 6);
    *ptr++ = stride;
    for (E_Int j = start; j < end; j++) {
      *ptr++ = cells[j]+1;
    }
  }

  PyObject* m = K_ARRAY::buildArray(local_crd, varString, cn, 8, NULL, false);

  // Release numpys
  l = PyList_GetItem(array, 0);
  Py_DECREF(PyList_GetItem(l, 0));
  Py_DECREF(PyList_GetItem(l, 1));
  Py_DECREF(PyList_GetItem(l, 2));
  Py_DECREF(PyList_GetItem(l, 3));
  Py_DECREF(PyList_GetItem(l, 4));
  Py_DECREF(PyList_GetItem(l, 5));
  Py_DECREF(PyList_GetItem(l, 6));
  
  //Py_INCREF(Py_None);
  //return Py_None;
  
  return m;    
}
*/

const E_Int decomp_vector[3] = {3,2,1};

static
void assign_to_processor_group(E_Int ncells, E_Int *processor_groups, E_Int n_proc_groups)
{
  E_Int size = ncells / n_proc_groups;
  E_Int sizeb = size + 1;
  E_Int fat_procs =  ncells - size*n_proc_groups;

  E_Int count = 0;
  E_Int j = 0;

  // procs with more elements than "size"
  for (j = 0; j < fat_procs; j++) {
    for (E_Int k = 0; k < sizeb; k++)
      processor_groups[count++] = j;
  }

  // the remaining procs
  for (; j < n_proc_groups; j++) {
    for (E_Int k = 0; k < size; k++)
      processor_groups[count++] = j;
  }
}

E_Int *simple_geom_decomp(E_Int ncells, E_Float *xyz)
{
  E_Int *final_decomp = (E_Int *)malloc(ncells * sizeof(E_Int));
  E_Int *processor_groups = (E_Int *)malloc(ncells * sizeof(E_Int));
  E_Int *cell_indices = (E_Int *)malloc(ncells * sizeof(E_Int));
  for (E_Int i = 0; i < ncells; i++)
    cell_indices[i] = i;

  // small angle rotation matrix to avoid jitters when mesh is aligned with XYZ
  E_Float angle = 0.001;
  E_Float c = 1. - 0.5*angle*angle; // taylor expansion of cos
  E_Float s = angle; // taylor expansion of sin
  E_Float c2 = c*c;
  E_Float s2 = s*s;

  E_Float rot[3][3] = {
    c2, c*s2 - s*c, s*c2 + s2,
    s*c, s2*s + c2, s2*c - c*s,
    -s, c*s, c2
  };

  E_Float *rotated_cells = (E_Float *)malloc(3*ncells * sizeof(E_Float));
  for (E_Int i = 0; i < ncells; i++) {
    E_Float *cx = &xyz[3*i];
    E_Float *rx = &rotated_cells[3*i];
    rx[0] = rot[0][0]*cx[0] + rot[0][1]*cx[1] + rot[0][2]*cx[2];
    rx[1] = rot[1][0]*cx[0] + rot[1][1]*cx[1] + rot[1][2]*cx[2];
    rx[2] = rot[2][0]*cx[0] + rot[2][1]*cx[1] + rot[2][2]*cx[2];
  }

  // sort by X
  std::sort(cell_indices, cell_indices+ncells, [&](E_Int a, E_Int b)
                                               {
                                                 return rotated_cells[3*a] < rotated_cells[3*b];
                                               });


  assign_to_processor_group(ncells, processor_groups, decomp_vector[0]);

  for (E_Int i = 0; i < ncells; i++)
    final_decomp[cell_indices[i]] = processor_groups[i];

  // sort by Y
  std::sort(cell_indices, cell_indices+ncells, [&](E_Int a, E_Int b)
                                               {
                                                 return rotated_cells[3*a+1] < rotated_cells[3*b+1];
                                               });

  assign_to_processor_group(ncells, processor_groups, decomp_vector[1]);

  for (E_Int i = 0; i < ncells; i++)
    final_decomp[cell_indices[i]] += decomp_vector[0] * processor_groups[i];

  // sort by Z
  std::sort(cell_indices, cell_indices+ncells, [&](E_Int a, E_Int b)
                                               {
                                                 return rotated_cells[3*a+2] < rotated_cells[3*b+2];
                                               });

  assign_to_processor_group(ncells, processor_groups, decomp_vector[2]);

  for (E_Int i = 0; i < ncells; i++)
    final_decomp[cell_indices[i]] += decomp_vector[0] * decomp_vector[1] * processor_groups[i];

  return final_decomp;
}

void sample_sort(E_Int ncells, std::vector<ikv> &elements, E_Int *cells_dist)
{
  E_Int nproc, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // number of local samples
  E_Int nlsamples = nproc;

  // sort the local elements
  std::sort(elements.begin(), elements.end(), compare_ikv_2);

  // select the local nlsamples equally spaced elements
  std::vector<ikv> mypicks(nproc+1, ikv());

  for (E_Int i = 0; i < nlsamples-1; i++) {
    E_Int k = ((i+1)*ncells/nlsamples)%ncells;
    mypicks[i].key = elements[k].key;
    mypicks[i].val = elements[k].val;
  }

  // gather all the picks
  std::vector<ikv> allpicks(nproc*nlsamples, ikv());

  MPI_Allgather(&mypicks[0], 2*(nlsamples-1), MPI_INT,
                &allpicks[0], 2*(nlsamples-1), MPI_INT,
                MPI_COMM_WORLD);

  // remove any bad samples
  E_Int ntsamples = 0;
  for (E_Int i = 0; i < nproc*nlsamples; i++) {
    if (allpicks[i].key == -1)
      break;
    ntsamples++;
  }

  std::sort(allpicks.begin(), allpicks.begin() + ntsamples, compare_ikv_2);

  // select the final splitters
  for (E_Int i = 1; i < nproc; i++)
    mypicks[i] = allpicks[i*ntsamples/nproc];
  mypicks[0].key = -1;
  mypicks[nproc].key = INT_MAX;

  // compute the number of elements that belong to each bucket
  std::vector<E_Int> scount(nproc, 0);
  E_Int j = 0;
  for (E_Int i = 0; i < ncells; i++) {
    if ((elements[i].key < mypicks[j+1].key) ||
       ((elements[i].key == mypicks[j+1].key) && (elements[i].val < mypicks[j+1].val)))
      scount[j]++;
    else
      scount[++j]++;
  }

  std::vector<E_Int> rcount(nproc);
  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> sdist(nproc+1);
  std::vector<E_Int> rdist(nproc+1);
  sdist[0] = rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    scount[i] *= 2; // ikv
    rcount[i] *= 2;
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  E_Int nrecv = rdist[nproc]/2; // number of received elements
  std::vector<ikv> relements(nrecv);
  
  MPI_Alltoallv(&elements[0], &scount[0], &sdist[0], MPI_INT,
                &relements[0], &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  std::sort(relements.begin(), relements.end(), compare_ikv);


  E_Int lastcell;
  MPI_Scan(&nrecv, &lastcell, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  E_Int firstcell = lastcell - nrecv;

  // shift excess cells wrt [cells_dist[rank], cells_dist[rank+1]] to rank-1 / rank+1
  j = 0;
  for (E_Int i = 0; i < nproc; i++) {
    if (cells_dist[i+1] > firstcell) {
      if (cells_dist[i+1] >= lastcell) {
        for (E_Int k = 0; k < lastcell-firstcell; k++)
          relements[j++].key = i;
      } else {
        for (E_Int k = 0; k < cells_dist[i+1]-firstcell; k++)
          relements[j++].key = i;

        firstcell = cells_dist[i+1];
      }
    }

    if (cells_dist[i+1] >= lastcell)
      break;
  }
  assert(j == nrecv);

  // send back
  MPI_Alltoallv(&relements[0], &rcount[0], &rdist[0], MPI_INT,
                &elements[0], &scount[0], &sdist[0], MPI_INT,
                MPI_COMM_WORLD);

}

void bin_coords(E_Int ncells, E_Float *xyz, E_Int nbins, E_Int *bxyz)
{
  E_Int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::vector<ikv> buckets(nbins);
  std::vector<rkv> cand(ncells);
  E_Float gmin, gmax;
  std::vector<E_Float> emarkers(nbins+1);
  std::vector<E_Float> nemarkers(nbins+1);
  std::vector<E_Int> lcounts(nbins);
  std::vector<E_Int> gcounts(nbins);
  std::vector<E_Int> lsums(nbins);
  std::vector<E_Int> gsums(nbins);
  E_Int cnbins;
  E_Float gsum;
  E_Int gncells;
  MPI_Allreduce(&ncells, &gncells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  for (E_Int dim = 0; dim < 3; dim++) {
    E_Float sum = 0.;
    for (E_Int i = 0; i < ncells; i++) {
      cand[i].key = xyz[3*i+dim];
      cand[i].val = i;
      sum += cand[i].key;
    }
    std::sort(cand.begin(), cand.end(), compare_rkv);

    // determine initial range
    MPI_Allreduce(&cand[0].key, &gmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&cand[ncells-1].key, &gmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&sum, &gsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    emarkers[0] = gmin;
    emarkers[1] = gsum/gncells;
    emarkers[2] = gmax*1.001;
    cnbins = 2;

    while (cnbins < nbins) {
      for (E_Int i = 0; i < cnbins; i++) {
        lcounts[i] = 0;
        lsums[i] = 0;
      }

      E_Int j = 0;
      for (E_Int i = 0; i < ncells; ) {
        if (cand[i].key < emarkers[j+1]) {
          lcounts[j]++;
          lsums[j] += cand[i].key;
          i++;
        } else {
          j++;
        }
      }

      MPI_Allreduce(&lcounts[0], &gcounts[0], cnbins, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&lsums[0], &gsums[0], cnbins, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      for (E_Int i = 0; i < cnbins; i++) {
        buckets[i].key = gcounts[i];
        buckets[i].val = i;
      }
      std::sort(buckets.begin(), buckets.end(), compare_ikv);

      j = 0;
      for (E_Int i = cnbins-1; i >= 0; i--, j++) {
        E_Int l = buckets[i].val;
        if (buckets[i].key > gncells/nbins && cnbins < nbins) {
          nemarkers[j++] = (emarkers[l]+emarkers[l+1])/2;
          cnbins++;
        }
        nemarkers[j] = emarkers[l];
      }
      assert(cnbins == j);
      
      std::sort(&nemarkers[0], &nemarkers[0]+cnbins);
      for (E_Int i = 0; i < cnbins; i++)
        emarkers[i] = nemarkers[i];
      emarkers[cnbins] = gmax*1.001;
    }

    // assign the coordinates to the appropriate bin
    E_Int j = 0;
    for (E_Int i = 0; i < ncells;) {
      if (cand[i].key < emarkers[j+1]) {
        bxyz[3*cand[i].val + dim] = j;
        i++;
      } else {
        j++;
      }
    }
  }
}
