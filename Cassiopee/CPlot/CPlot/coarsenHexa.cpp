/*    
    Copyright 2013-2026 ONERA.

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
#include "kcore.h"
#include "cplot.h"
#include <unordered_map>

#define TARGET_NCELLS 100
#define FREE_CELL -1

typedef struct graph {
  E_Int ncells;
  E_Int nedges;
  E_Int *xadj;
  E_Int *cadj;
  E_Int *coarse_map;
  struct graph *coarser;
  struct graph *finer;
  E_Float *cell_centers;
} graph;

static
graph *setup_graph(E_Int ncells, E_Int *xadj, E_Int *cadj, E_Float *cc)
{
  graph *G = (graph *)malloc(sizeof(*G));

  G->ncells = ncells;
  G->nedges = xadj[ncells];
  G->xadj = xadj;
  G->cadj = cadj;
  G->coarse_map = NULL;
  G->coarser = NULL;
  G->finer = NULL;
  G->cell_centers = cc;

  return G;
}

static
graph *create_graph()
{
  graph *G = (graph *)malloc(sizeof(*G));
  
  G->ncells = -1;
  G->nedges = -1;
  G->xadj = NULL;
  G->cadj = NULL;
  G->coarse_map = NULL;
  G->coarser = NULL;
  G->finer = NULL;
  G->cell_centers = NULL;
  
  return G;
}

static
graph *setup_coarse_graph(graph *G, E_Int coarse_ncells)
{
  graph *cG = create_graph();

  cG->ncells = coarse_ncells;
  
  cG->xadj = (E_Int *)malloc((coarse_ncells+1) * sizeof(E_Int));
  cG->cadj = (E_Int *)malloc((G->nedges+1) * sizeof(E_Int));

  cG->finer = G;
  G->coarser = cG;

  cG->cell_centers = (E_Float *)malloc(3*coarse_ncells * sizeof(E_Float));

  return cG;
}


static
void create_coarse_graph(graph *G, E_Int coarse_ncells, const E_Int *match)
{
  graph *cG = setup_coarse_graph(G, coarse_ncells);
  auto &cxadj = cG->xadj;
  E_Int *ccadj = &cG->cadj[0];
  
  E_Int ncells = G->ncells;
  E_Int *xadj = G->xadj;
  E_Int *cadj = G->cadj;
  E_Int *coarse_map = G->coarse_map;

  E_Int coarse_nedges, nedges;
  cxadj[0] = coarse_ncells = coarse_nedges = 0;

  std::unordered_map<E_Int, E_Int> emap;

  E_Float *cc = G->cell_centers;
  E_Float *coarse_cc = cG->cell_centers;

  for (E_Int cell = 0; cell < ncells; cell++) {
    // treat edge only once
    E_Int nei = match[cell];
    if (nei < cell)
      continue;

    // compute coarse cell center
    E_Float *p0 = &cc[3*cell];
    E_Float *p1 = &cc[3*nei];
    E_Float *pt = &coarse_cc[3*coarse_ncells];
    
    pt[0] = pt[1] = pt[2] = 0;
    for (E_Int j = 0; j < 3; j++)
      pt[j] += p0[j] + p1[j];
    for (E_Int j = 0; j < 3; j++)
      pt[j] *= 0.5;

    nedges = 0;
    for (E_Int j = xadj[cell]; j < xadj[cell+1]; j++) {
      if (cadj[j] == nei) continue;
      E_Int k = coarse_map[cadj[j]];

      auto search = emap.find(k);
      if (search == emap.end()) {
        ccadj[nedges] = k;
        emap[k] = nedges++;
      }
    }
    
    if (nei != cell) {
      for (E_Int j = xadj[nei]; j < xadj[nei+1]; j++) {
        if (cadj[j] == cell) continue;
        E_Int k = coarse_map[cadj[j]];

        auto search = emap.find(k);
        if (search == emap.end()) {
          ccadj[nedges] = k;
          emap[k] = nedges++;
        }
      }
    }

    emap.clear();

    ccadj += nedges;
    coarse_nedges += nedges;
    cxadj[++coarse_ncells] = coarse_nedges;
  }

  cG->nedges = coarse_nedges;
}

// Computes matching of the graph by randomly selecting one of the unmatched adjacent cells
static
E_Int graph_random_match(graph *G)
{
  E_Int ncells = G->ncells;
  auto &xadj = G->xadj;
  auto &cadj = G->cadj;
  auto &coarse_map = G->coarse_map;
  auto &cc = G->cell_centers;

  E_Int *match = (E_Int *)malloc(ncells * sizeof(E_Int));
  memset(match, -1, ncells*sizeof(E_Int));
  
  E_Int n_free_cells = 0;
  E_Int coarse_ncells = 0;

  for (E_Int i = 0; i < ncells; i++) {
    if (match[i] == FREE_CELL) {
      E_Int idx = i;

      for (E_Int j = xadj[i]; j < xadj[i+1]; j++) {
        E_Int k = cadj[j];
        if (match[k] == FREE_CELL) {
          idx = k;
          break;
        }
      }

      if (idx == i) {
        n_free_cells++;
        idx = FREE_CELL;
      }

      if (idx != FREE_CELL) {
        coarse_map[i] = coarse_map[idx] = coarse_ncells++;
        match[i] = idx;
        match[idx] = i;
      }
    }
  }

  printf("coarse_ncells: %d\n", coarse_ncells);
  printf("n_free_cells: %d\n", n_free_cells);

  // match the unmatched vertices with themselves
  for (E_Int i = 0; i < ncells; i++) {
    if (match[i] == FREE_CELL) {
      match[i] = i;
      coarse_map[i] = coarse_ncells++;
    }
  }

  create_coarse_graph(G, coarse_ncells, match);

  free(match);

  return coarse_ncells;
}

static
graph *coarsen_graph(graph *G)
{
  do {
    if (!G->coarse_map)
      G->coarse_map = (E_Int *)malloc(G->ncells * sizeof(E_Int));
    
    graph_random_match(G);

    G = G->coarser;
  } while (G->ncells > TARGET_NCELLS);

  return G;
}

PyObject* K_CPLOT::coarsenHexa(PyObject *self, PyObject *args)
{
  PyObject *array;
  if (!PYPARSETUPLE_(args, O_, &array)) 
  {
    PyErr_SetString(PyExc_TypeError, "CoarsenHexa: couldn't read mesh.");
    return NULL;
  }

  char *varString, *eltType;
  FldArrayF *f;
  FldArrayI *cn;
  E_Int ni, nj, nk;
  E_Int posx, posy, posz;
  posx = posy = posz = -1;

  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, eltType);

  if (res != 1 && res != 2) 
  {
    PyErr_SetString(PyExc_TypeError, "CoarsenHexa: invalid array.");
    exit(1);
  }

  if (res == 1 || strcmp(eltType, "NGON")) {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError, "CoarsenHexa: works only for NGON meshes.");
    exit(1);
  }

  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1) {
    PyErr_SetString(PyExc_TypeError, "Coarsen_NGON: coordinates not found in array.");
    RELEASESHAREDB(res, array, f, cn);
    exit(1);
  }

  posx++; posy++; posz++;

  E_Int *cnp = cn->begin();
  E_Int nfaces = cnp[0];
  E_Int sizeFN = cnp[1];
  E_Int ncells = cnp[sizeFN+2];

  printf("nfaces: %d\n", nfaces);
  printf("ncells: %d\n", ncells);
  
  // create cell adjacency
  E_Int * cadj = (E_Int *)malloc(6*ncells * sizeof(E_Int));
  memset(cadj, -1, 6*ncells*sizeof(E_Int));
  
  FldArrayI cFE;
  K_CONNECT::connectNG2FE(*cn, cFE);

  E_Int *cFE1 = cFE.begin(1);	
  E_Int *cFE2 = cFE.begin(2);	

  E_Int *count_neis = (E_Int *)calloc(ncells, sizeof(E_Int));

  for (E_Int i = 0; i < nfaces; i++) {
    E_Int nei = cFE2[i]-1;

    if (nei != -1) {
      E_Int own = cFE1[i]-1;
      cadj[6*own + count_neis[own]++] = nei;
      cadj[6*nei + count_neis[nei]++] = own;
    }
  }

  E_Int * xadj = (E_Int *)malloc((ncells+1) * sizeof(E_Int));
  xadj[0] = 0;

  E_Int nc = 0;
  E_Int k = 0;
  for (E_Int i = 0; i < 6*ncells;) {
    if (cadj[i] != -1)
      cadj[k++] = cadj[i++];
    else {
      xadj[++nc] = k;
      while (cadj[i] == -1)
        i++;
    }
  }

  E_Int nedges = xadj[ncells]/2;

  // compute face centers
  FldArrayI posFaces;
  K_CONNECT::getPosFaces(*cn, posFaces);
  E_Float *px = f->begin(posx);
  E_Float *py = f->begin(posy);
  E_Float *pz = f->begin(posz);
  FldArrayF face_coords;
  face_coords.malloc(nfaces, 3);
  face_coords.setAllValuesAtNull();

  for (E_Int i = 0; i < nfaces; i++) {
    cnp = cn->begin();
    cnp += posFaces[i];
    E_Int stride = cnp[0];
    E_Float inv_stride = 1./stride;
    assert(stride == 4); // Hexa for now
    cnp++;
    for (E_Int j = 0; j < stride; j++) {
      E_Int point = cnp[j]-1;
      face_coords(i, 1) += px[point];
      face_coords(i, 2) += py[point];
      face_coords(i, 3) += pz[point];
    }
    for (E_Int j = 1; j <= 3; j++)
      face_coords(i, j) *= inv_stride;
  }

  FldArrayI posCells;
  K_CONNECT::getPosElts(*cn, posCells);
  E_Float *cell_centers = (E_Float *)calloc(3*ncells, sizeof(E_Float));
  
  const E_Float one_sixth = 1./6.;

  for (E_Int i = 0; i < ncells; i++) {
    cnp = cn->begin();
    cnp += posCells[i];
    E_Int stride = cnp[0];
    assert(stride == 6); // Hexa for now
    cnp++;
    E_Float *pt = &cell_centers[3*i];
    for (E_Int j = 0; j < stride; j++) {
      E_Int face = cnp[j]-1;
      pt[0] += face_coords(face, 1); 
      pt[1] += face_coords(face, 2); 
      pt[2] += face_coords(face, 3);
    }
    for (E_Int j = 0; j < 3; j++)
      pt[j] *= one_sixth;
  }
  
  // create dual graph
  graph *G = setup_graph(ncells, xadj, cadj, cell_centers);
  
  // coarsen
  graph *cG = coarsen_graph(G);

  PyObject *l;
  return l;
}
