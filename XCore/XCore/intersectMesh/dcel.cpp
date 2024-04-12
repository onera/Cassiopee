/*    
    Copyright 2013-2024 Onera.

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
#include <cassert>

// M is a 2D mesh
dcel::dcel(Mesh *M)
{
  points.resize(M->np, NULL);
  half_edges.resize(2*M->nf);
  faces.resize(M->nc);

  // Make PE for edges
  std::vector<E_Int> PE(2*M->nf, -1);
  for (E_Int i = 0; i < M->nc; i++) {
    for (E_Int j = M->indPH[i]; j < M->indPH[i+1]; j++) {
      E_Int e = M->nface[j];
      if (PE[2*e] == -1) PE[2*e] = i;
      else PE[2*e+1] = i;
    }
  }

  // Points
  for (E_Int i = 0; i < M->np; i++) {
    points[i] = new point(M->x[i], M->y[i], i, -1);
    // Incident Halfedge is set later
  }

  // Half-edges and faces
  
  for (E_Int i = 0; i < M->nf; i++) {
    E_Int left_face = PE[2*i];
    E_Int right_face = PE[2*i+1];
    assert(left_face != -1);

    E_Int np = -1;
    E_Int *pn = mesh_get_face(i, np, M);
    assert(np == 2);

    E_Int orig = pn[0];
    E_Int dest = pn[1];

    // Set points half-edges
    points[orig]->he = i;
    points[dest]->he = i+M->nf;

    // Set half-edge components
    auto &he = half_edges[i];
    auto &te = half_edges[i+M->nf];
    
    he.orig = orig;
    he.twin = i + M->nf;
    he.face = left_face;

    te.orig = dest;
    te.twin = i;
    te.face = right_face;

    // Set right and left faces outer components
    faces[left_face].outer_component = i;
    if (right_face != -1) faces[right_face].outer_component = i + M->nf;
  }

  // Set next and prev pointers
  for (E_Int i = 0; i < M->nc; i++) {
    E_Int ne = -1;
    E_Int *pe = mesh_get_cell(i, ne, M);

    for (E_Int j = 0; j < ne; j++) {
      E_Int e1 = pe[j];
      E_Int e2 = pe[(j+1)%ne];

      E_Int he1 = PE[2*e1] == i ? e1 : e1+M->nf;
      E_Int he2 = PE[2*e2] == i ? e2 : e2+M->nf;

      half_edges[he1].next = he2;
      half_edges[he2].prev = he1;
    }
  }

  assert(is_valid());
}

dcel::dcel(E_Int ni, E_Int nj,
  const std::vector<E_Float> &X, const std::vector<E_Float> &Y,
  const std::vector<E_Float> &Z)
{
  Mesh *M = mesh_make_surface_mesh_from_structured_points(ni, nj, X, Y, Z);

  mesh_write(M, "structured");

  *this = dcel(M);

  // TODO(Imad): delete M
  //mesh_drop(M);
}

int dcel::_check_duplicate_points()
{
  return 1;
}

int dcel::is_valid()
{
  // Tests
  for (size_t i = 0; i < points.size(); i++) {
    point *p = points[i];

    assert(p->he < half_edges.size());
    
    const auto& he = half_edges[p->he];

    E_Int orig = he.orig;
    E_Int dest = half_edges[he.twin].orig;

    if ((int)i != orig && (int)i != dest)
      return 0;
  }

  for (size_t i = 0; i < faces.size(); i++) {
    size_t he = faces[i].outer_component;
    if (he == -1) return 0;
    
    E_Int nwalks = 0;
    size_t next = -1;
    size_t curr = he;
    while (next != he) {
      if (++nwalks > 10) break;

      next = half_edges[curr].next;
      curr = next;
    }

    if (nwalks > 10) {
      printf("Warning: couldn't reach back incident edge after " SF_D_ " walks",
        nwalks);
      return 0;
    }
  }

  return 1;
}

dcel::dcel(const dcel &D0, const dcel &D1)
{
  points.resize(D0.points.size() + D1.points.size());
  half_edges.resize(D0.half_edges.size() + D1.half_edges.size());
  faces.resize(D0.faces.size() + D1.faces.size());

  // Copy points
  point **v_it = &points[0];
  for (size_t i = 0; i < D0.points.size(); i++) {
    *v_it++ = D0.points[i];
  }
  for (size_t i = 0; i < D1.points.size(); i++) {
    *v_it++ = D1.points[i];
  }

  // Copy half_edges
  half_edge *he_it = &half_edges[0];
  for (size_t i = 0; i < D0.half_edges.size(); i++) {
    *he_it++ = D0.half_edges[i];
  }
  for (size_t i = 0; i < D1.half_edges.size(); i++) {
    *he_it++ = D1.half_edges[i];
  }

  // Copy faces
  face *f_it = &faces[0];
  for (size_t i = 0; i < D0.faces.size(); i++) {
    *f_it++ = D0.faces[i];
  }
  for (size_t i = 0; i < D1.faces.size(); i++) {
    *f_it++ = D1.faces[i];
  }

  // Increment record ids from D1

  // Points
  for (size_t i = D0.points.size(); i < points.size(); i++) {
    points[i]->he += D0.half_edges.size();
    points[i]->id += D0.points.size();
  }

  // Half_edges
  for (size_t i = D0.half_edges.size(); i < half_edges.size(); i++) {
    half_edges[i].orig += D0.points.size();
    half_edges[i].twin += D0.half_edges.size();
    half_edges[i].next += D0.half_edges.size();
    half_edges[i].prev += D0.half_edges.size();
    half_edges[i].face += D0.faces.size();
  }

  // Faces
  for (size_t i = D0.faces.size(); i < faces.size(); i++) {
    faces[i].outer_component += D0.half_edges.size();
    if (faces[i].inner_components.size()) {
      fprintf(stderr, "Warning: face inner components unimplemented\n");
      assert(0);
    }
    //for (size_t j = 0; j < faces[i].inner_components.size(); j++)
    //  faces[i].inner_components[j] += D0.half_edges.size();
  }

  assert(is_valid());
}

void dcel::resolve()
{
  // Make the segments
  std::vector<segment *> segs;

  {
    std::vector<int> marked_edges(half_edges.size(), 0);

    for (size_t i = 0; i < half_edges.size(); i++) {
      // Skip if already processed
      if (marked_edges[i]) continue;

      const auto &HE = half_edges[i];
      size_t te = HE.twin;
      const auto &TE = half_edges[te];

      marked_edges[i] = 1;
      marked_edges[te] = 1;

      segment *s = new segment(points[HE.orig], points[TE.orig], i, i);
      assert(s->p != s->q);

      if (!segment_is_lexico(s)) {
          std::swap(s->p, s->q);
          s->id = i;
          s->he = te;
      }

      segs.push_back(s);
    }
  }

  // Point to list of intersecting segments
  std::unordered_map<point *, std::vector<segment *>> Xsegs;

  _sweep(segs, points, Xsegs);
}
