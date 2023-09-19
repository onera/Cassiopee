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

PyObject *K_XCORE::chunk2part_elt(PyObject *self, PyObject *args)
{
  PyObject *XYZ, *NAME, *STRIDE, *CN;
  char *eltName;

  if (!PyArg_ParseTuple(args, "OsOO", &XYZ, &eltName, &STRIDE, &CN)) {
    PyErr_SetString(PyExc_ValueError, "chunk2part_elt(): bad input");
    return NULL;
  }

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
  E_Int rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  
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

  // make pointCells connectivity
  std::unordered_map<E_Int, std::vector<E_Int>> point2cells;
  E_Int firstCell = cdist[rank];
  for (E_Int i = 0; i < ncells; i++) {
    E_Int *pn = &cn[stride*i];
    for (E_Int j = 0; j < stride; j++)
      point2cells[pn[j]-1].push_back(firstCell+i);
  }

  // dispatch points
  std::vector<E_Int> scount(nproc, 0);
  for (const auto& point : point2cells) {
    E_Int target = get_proc(point.first, pdist, nproc);
    assert(target >= 0 && target < nproc);
    scount[target] += 1 + 1 + point.second.size(); // id + size + cells
  }

  std::vector<E_Int> rcount(nproc);
  MPI_Alltoall(&scount[0], 1, MPI_INT, &rcount[0], 1, MPI_INT, MPI_COMM_WORLD);

  std::vector<E_Int> sdist(nproc+1), rdist(nproc+1);
  sdist[0] = rdist[0] = 0;
  for (E_Int i = 0; i < nproc; i++) {
    sdist[i+1] = sdist[i] + scount[i];
    rdist[i+1] = rdist[i] + rcount[i];
  }

  std::vector<E_Int> sdata(sdist[nproc]);
  std::vector<E_Int> rdata(rdist[nproc]);

  std::vector<E_Int> idx(nproc);
  for (E_Int i = 0; i < nproc; i++)
    idx[i] = sdist[i];

  for (const auto& point : point2cells) {
    E_Int target = get_proc(point.first, pdist, nproc);
    assert(target >= 0 && target < nproc);
    sdata[idx[target]++] = point.first;
    sdata[idx[target]++] = point.second.size();
    for (const auto& cell : point.second)
      sdata[idx[target]++] = cell;
  }

  for (E_Int i = 0; i < nproc; i++)
    assert(idx[i] == sdist[i+1]);

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], MPI_INT,
                &rdata[0], &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  point2cells.clear();

  std::vector<std::vector<E_Int>> pointCells(npoints);
  E_Int firstPoint = pdist[rank];

  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &rdata[rdist[i]];
    for (E_Int j = 0; j < rcount[i];) {
      E_Int point = ptr[j++];
      assert(get_proc(point, pdist, nproc) == rank);
      E_Int size = ptr[j++];
      E_Int lp = point - firstPoint;
      for (E_Int k = 0; k < size; k++)
        pointCells[lp].push_back(ptr[j++]);
    }
  }

  if (rank == 0)
    puts("PointCells OK");
  

  std::unordered_map<E_Int, std::set<E_Int>> CADJ;

  // build graph
  // edge between two cells: shared point
  for (E_Int i = 0; i < npoints; i++) {
    const auto& neis = pointCells[i];
    for (size_t j = 0; j < neis.size(); j++) {
      E_Int cell = neis[j];
      for (size_t k = j+1; k < neis.size(); k++) {
        E_Int nei = neis[k];
        CADJ[cell].insert(nei);
        CADJ[nei].insert(cell);
      }
    }
  }

  for (E_Int i = 0; i < nproc; i++)
    scount[i] = 0;
  
  for (const auto &cell : CADJ) {
    E_Int target = get_proc(cell.first, cdist, nproc);
    scount[target] += 1 + 1 + cell.second.size(); // id + size + neis
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
    E_Int target = get_proc(cell.first, cdist, nproc);
    sdata[idx[target]++] = cell.first;
    sdata[idx[target]++] = cell.second.size();
    for (const auto& elem : cell.second)
      sdata[idx[target]++] = elem;
  }

  MPI_Alltoallv(&sdata[0], &scount[0], &sdist[0], MPI_INT,
                &rdata[0], &rcount[0], &rdist[0], MPI_INT,
                MPI_COMM_WORLD);

  CADJ.clear();

  std::vector<std::set<E_Int>> cadj(ncells);
  E_Int nedges = 0;
  for (E_Int i = 0; i < nproc; i++) {
    E_Int *ptr = &rdata[rdist[i]];
    for (E_Int j = 0; j < rcount[i];) {
      E_Int cell = ptr[j++] - firstCell;
      E_Int size = ptr[j++];
      for (E_Int k = 0; k < size; k++) {
        if (cadj[cell].find(ptr[j]) == cadj[cell].end()) {
          nedges += 1;
          cadj[cell].insert(ptr[j++]);
        } else {
          j++;
        }
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