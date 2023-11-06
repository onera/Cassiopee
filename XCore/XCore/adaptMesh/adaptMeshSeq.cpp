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
#include "proto.h"
#include <stack>

/*
 T = initial cgns tree
 H = createAdaptMesh(T)

 for i in range(niter):
  L = X.extractLeafMesh(H)
  fld = makeFieldSensor(L)
  H = X.adaptMeshSeq(H, fld)
*/

#define TAG 1.0
#define TOL 1e-12

static
E_Int FEQ(E_Float a, E_Float b)
{
  return fabs(a-b) < TOL;
}

static
void get_higher_level_neighbours(E_Int cell, E_Int face, mesh *M, std::vector<E_Int> &neis)
{
  auto& FT = M->ftree;
  E_Int nchild = tree_get_nchildren(&FT, face);
  //assert(nchild == 4);
  E_Int *children = tree_get_children(&FT, face);
  for (E_Int i = 0; i < nchild; i++) {
    E_Int fchild = children[i];
    E_Int nei = M->owner[face] == cell ? M->neigh[fchild] : M->owner[fchild];
    neis.push_back(nei);
  }
}

PyObject *K_XCORE::adaptMeshSeq(PyObject *self, PyObject *args)
{
  PyObject *MESH, *field, *FV=NULL;
  if (!PYPARSETUPLE_(args, OOO_, &MESH, &field, &FV)) {
    PyErr_SetString(PyExc_ValueError, "adaptMeshSeq: wrong input.");
    return NULL;
  }

  E_Int ret;

  // Unpack adaptMesh
  if (!PyCapsule_IsValid(MESH, "adaptMesh")) {
    PyErr_SetString(PyExc_TypeError, "adaptMeshSeq: bad adaptMesh hook.");
    return NULL;
  }
  mesh *M = (mesh *)PyCapsule_GetPointer(MESH, "adaptMesh");

  // Check field
  E_Float *fld;
  E_Int field_size, field_nfld;
  ret = K_NUMPY::getFromNumpyArray(field, fld, field_size, field_nfld, true);
  
  if (ret != 1) {
    PyErr_SetString(PyExc_TypeError, "adaptMeshSeq: field should be an array.");
    return NULL;
  }

  if (field_size != M->ctree.nleaves || field_nfld != 1) {
    PyErr_SetString(PyExc_TypeError, "adaptMeshSeq: bad field size.");
    return NULL;
  }

  // Freeze Vector
  E_Float *freezeVector = NULL;
  if (FV != Py_None) {
    ret = K_NUMPY::getFromNumpyArray(FV, freezeVector, field_size, field_nfld, true);
    if (ret == 0 || field_nfld != 1 || field_size != 3) {
      PyErr_SetString(PyExc_TypeError, "adaptMeshSeq: bad freeze vector.");
      return NULL;
    }
  }

  // Isolate ref cells
  auto &CT = M->ctree;
  auto &FT = M->ftree;
  E_Int nleaves = CT.nleaves;
  std::set<E_Int> rcells;
  for (E_Int i = 0; i < nleaves; i++) {
    if (!FEQ(fld[i], TAG)) continue;
    E_Int cell = CT.l2g[i];
    rcells.insert(cell);
  }

  // Smooth out ref data
  std::stack<E_Int> stk;
  for (auto cell : rcells) stk.push(cell);

  while (!stk.empty()) {
    E_Int cell = stk.top();
    stk.pop();

    E_Int incr_cell = 1;
    incr_cell += CT.level[cell];

    E_Int *pf = &M->NFACE[M->xcells[cell]];

    std::vector<E_Int> neis;

    for (E_Int i = 0; i < 6; i++) {
      E_Int face = pf[i];
      E_Int nei = get_neighbour(cell, face, M);
      if (nei == -1) continue;
      if (CT.enabled[nei])// || (CT.level[nei] > 0 && CT.enabled[CT.parent[nei]]))
        neis.push_back(nei);
      else
        get_higher_level_neighbours(cell, face, M, neis);
    }

    for (size_t i = 0; i < neis.size(); i++) {
      E_Int nei = neis[i];
      E_Int nei_in_rcells = rcells.find(nei) != rcells.end();
      E_Int incr_nei = nei_in_rcells ? 1 : 0;
      incr_nei += CT.level[nei];
      if (abs(incr_nei - incr_cell) <= 1) continue;

      E_Int cell_to_mod = incr_cell > incr_nei ? nei : cell;

      stk.push(cell_to_mod);
      rcells.insert(cell_to_mod);
    }

  }

  // Resize data structures
  E_Int nref_cells = rcells.size();
  resize_data_for_refinement(M, &CT, &FT, nref_cells);

  // Make ref dirs
  std::vector<E_Int> ref_cells;
  for (auto cell : rcells) ref_cells.push_back(cell);
  std::vector<E_Int> ref_data(3*nref_cells, 1);
  
  if (freezeVector) {
    std::vector<pDirs> Dirs(nref_cells);
    compute_cells_principal_vecs(M, ref_cells, Dirs);
    E_Int ok_vec = 1;
    for (E_Int i = 0; i < nref_cells; i++) {
      E_Int *pr = &ref_data[3*i];
      E_Float dp;

      dp = dot(Dirs[i].I, freezeVector, 3) / norm(Dirs[i].I, 3);
      if (fabs(fabs(dp)-1.0) < TOL) {
        pr[0] = 0;
        continue;
      }

      dp = dot(Dirs[i].J, freezeVector, 3) / norm(Dirs[i].J, 3);
      if (fabs(fabs(dp)-1.0) < TOL) {
        pr[1] = 0;
        continue;
      }

      dp = dot(Dirs[i].K, freezeVector, 3) / norm(Dirs[i].K, 3);
      if (fabs(fabs(dp)-1.0) < TOL) {
        pr[2] = 0;
        continue;
      }
      
      ok_vec = 0;
      fprintf(stderr, "Inconsistant freeze vector\n");
      break;
    }

    if (!ok_vec) {
      PyErr_SetString(PyExc_TypeError, "adaptMeshSeq: inconsistant freeze vector.");
      return NULL;
    }
  }

  // Adapt!
  for (E_Int i = 0; i < nref_cells; i++) {
    E_Int cell = ref_cells[i];
    E_Int *pr = &ref_data[3*i];
    if (pr[0] && pr[1] && pr[2]) {
      cut_cell_xyz(cell, M, &CT, &FT);
    } else if (pr[0] && pr[1] && !pr[2]) {
      cut_cell_xy(cell, M, &CT, &FT);
    } else if (pr[0] && !pr[1] && pr[2]) {
      cut_cell_xz(cell, M, &CT, &FT);
    } else if (!pr[0] && pr[1] && pr[2]) {
      cut_cell_yz(cell, M, &CT, &FT);
    } else if (pr[0] && !pr[1] && !pr[2]) {
      cut_cell_x(cell, M, &CT, &FT);
    } else if (!pr[0] && pr[1] && !pr[2]) {
      cut_cell_y(cell, M, &CT, &FT);
    } else if (!pr[0] && !pr[1] && pr[2]) {
      cut_cell_z(cell, M, &CT, &FT);
    } else {
      assert(0);
    }
  }

  Py_INCREF(Py_None);
  return Py_None;
}
