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

void mesh_write(Mesh *O, const char *file_name)
{
  FILE *fh = fopen(file_name, "w");
  assert(fh);

  // Points
  fprintf(fh, "Points\n");
  fprintf(fh, SF_D_ "\n", O->np);
  for (E_Int i = 0; i < O->np; i++)
    fprintf(fh, "%.12f %.12f %.12f\n", O->x[i], O->y[i], O->z[i]);

  // NGon
  fprintf(fh, "NGon\n");
  fprintf(fh, SF_D_ "\n", O->indPG[O->nf]);
  for (E_Int i = 0; i < O->nf; i++) {
    for (E_Int j = O->indPG[i]; j < O->indPG[i+1]; j++) {
      fprintf(fh, SF_D_ " ", O->ngon[j]);
    }
  }
  fprintf(fh, "\n");

  // NFace
  fprintf(fh, "NFace\n");
  fprintf(fh, SF_D_ "\n", O->indPH[O->nc]);
  for (E_Int i = 0; i < O->nc; i++) {
    for (E_Int j = O->indPH[i]; j < O->indPH[i+1]; j++) {
      fprintf(fh, SF_D_ " ", O->nface[j]);
    }
  }
  fprintf(fh, "\n");

  // IndPG
  fprintf(fh, "IndPG\n");
  fprintf(fh, SF_D_ "\n", O->nf+1);
  for (E_Int i = 0; i < O->nf+1; i++) {
    fprintf(fh, SF_D_ " ", O->indPG[i]);
  }
  fprintf(fh, "\n");

  // IndPH
  fprintf(fh, "IndPH\n");
  fprintf(fh, SF_D_ "\n", O->nc+1);
  for (E_Int i = 0; i < O->nc+1; i++) {
    fprintf(fh, SF_D_ " ", O->indPH[i]);
  }
  fprintf(fh, "\n");

  fclose(fh);
}

void point_set_write(const std::set<E_Int> &points, Mesh *M, const char *fname)
{
  FILE *fh = fopen(fname, "w");
  assert(fh);
  fprintf(fh, "Points\n");
  fprintf(fh, "%zu\n", points.size());
  for (auto point : points) {
    fprintf(fh, "%f %f %f\n", M->x[point], M->y[point], M->z[point]);
  }

  fclose(fh);
}

void edge_hits_write(const std::map<Edge_NO, Edge_Hit> &edge_hit_table,
  Mesh *M, const char *ename, const char *pname)
{
  FILE *eh = fopen(ename, "w");
  FILE *ph = fopen(pname, "w");
  assert(eh && ph);

  // Which points are present?
  std::map<E_Int, E_Int> PT;
  E_Int np = 0;
  std::vector<E_Int> points;

  for (auto edge : edge_hit_table) {
    E_Int p = edge.first.p;
    E_Int q = edge.first.q;

    if (PT.find(p) == PT.end()) {
      PT[p] = np++;
      points.push_back(p);
    }

    if (PT.find(q) == PT.end()) {
      PT[q] = np++;
      points.push_back(q);
    }
  }

  // Write xyz
  fprintf(eh, "Points\n");
  fprintf(eh, SF_D_ "\n", np);

  for (auto point : points)
    fprintf(eh, "%.12f %.12f %.12f\n", M->x[point], M->y[point], M->z[point]);

  // Write the edges
  fprintf(eh, "Edges\n");
  fprintf(eh, "%zu\n", edge_hit_table.size());

  fprintf(ph, "Points\n");
  fprintf(ph, "%zu\n", edge_hit_table.size());

  for (auto edge : edge_hit_table) {
    E_Int p = edge.first.p;
    E_Int q = edge.first.q;

    auto EH = edge.second;
    
    
    fprintf(ph, "%f %f %f\n", EH.x, EH.y, EH.z);

    fprintf(eh, SF_D2_ " ", PT[p], PT[q]);
  }

  fclose(eh);
  fclose(ph);
}

void point_hits_write(const std::map<E_Int, Edge_Hit> &point_hit_table,
  Mesh *M, const char *ename)
{
  FILE *eh = fopen(ename, "w");
  assert(eh);

  E_Int np = point_hit_table.size() * 2;

  fprintf(eh, "Points\n");
  fprintf(eh, SF_D_ "\n", np);
  
  for (auto edge : point_hit_table) {
    E_Int p = edge.first;
    const auto &EH = edge.second;

    fprintf(eh, "%f %f %f\n", M->x[p], M->y[p], M->z[p]);
    fprintf(eh, "%f %f %f\n", EH.x, EH.y, EH.z);
  }

  fprintf(eh, "Edges\n");
  fprintf(eh, SF_D_ "\n", np/2);
  for (E_Int i = 0; i < np; i++)
    fprintf(eh, SF_D_ " ", i);
  fprintf(eh, "\n");

  fclose(eh);
}

void point_hits_write(const std::unordered_map<E_Int, Edge_Hit> &ptable,
  const char *fname)
{
  FILE *fh = fopen(fname, "w");
  assert(fh);

  fprintf(fh, "POINTS\n");
  fprintf(fh, "%zu\n", ptable.size());

  for (auto P : ptable) {
    fprintf(fh, "%f %f %f\n", P.second.x, P.second.y, P.second.z);
  }

  fclose(fh);
}

void write_edge_path(Mesh *M, const Edge_NO &edge,
  const std::vector<E_Int> &path_faces, const char *fname)
{
  FILE *fh = fopen(fname, "w");
  assert(fh);

  fprintf(fh, "Points\n");
  fprintf(fh, "%zu\n", 2 + path_faces.size()*4);

  fprintf(fh, "%f %f %f\n", M->x[edge.p], M->y[edge.p], M->z[edge.p]);
  fprintf(fh, "%f %f %f\n", M->x[edge.q], M->y[edge.q], M->z[edge.q]);

  for (auto face : path_faces) {
    E_Int np = -1;
    E_Int *pn = mesh_get_face(face, np, M);
    for (E_Int i = 0; i < np; i++)
      fprintf(fh, "%f %f %f\n", M->x[pn[i]], M->y[pn[i]], M->z[pn[i]]);
  }

  fprintf(fh, "Quads\n");
  fprintf(fh, "%zu\n", path_faces.size());
  E_Int incr = 2;
  for (auto face : path_faces) {
    for (E_Int i = 0; i < 4; i++) {
      fprintf(fh, SF_D_ " ", incr);
      incr++;
    }
    fprintf(fh, "\n");
  }

  fclose(fh);
}

void write_point_list(E_Float *X, E_Float *Y, E_Float *Z, E_Int *list,
  E_Int n, const char *fname)
{
  FILE *fh = fopen(fname, "w");
  assert(fh);
  fprintf(fh, "POINTS\n");
  fprintf(fh, SF_D_ "\n", n);
  for (E_Int i = 0; i < n; i++) {
    E_Int ind = list[i];
    fprintf(fh, "%f %f %f\n", X[ind], Y[ind], Z[ind]);
  }
  fclose(fh);
}

void write_trimesh(const std::vector<Triangle> &TRIS, Mesh *M)
{
  FILE *fh = fopen("trimesh", "w");
  assert(fh);
  fprintf(fh, "Points\n");

  std::unordered_map<E_Int, E_Int> pmap;
  E_Int np = 0;

  // How many points
  for (const auto &TRI : TRIS) {
    E_Int A = TRI.A;
    E_Int B = TRI.B;
    E_Int C = TRI.C;

    if (pmap.find(A) == pmap.end()) {
      pmap[A] = np++;
    }

    if (pmap.find(B) == pmap.end()) {
      pmap[B] = np++;
    }

    if (pmap.find(C) == pmap.end()) {
      pmap[C] = np++;
    }
  }

  fprintf(fh, SF_D_ "\n", np);

  std::vector<E_Int> inv_map(np);
  for (auto P : pmap)
    inv_map[P.second] = P.first;

  for (auto P : inv_map) {
    fprintf(fh, "%f %f %f\n", M->x[P], M->y[P], M->z[P]);
  }

  fprintf(fh, "TRIS\n");
  fprintf(fh, "%zu\n", TRIS.size());

  for (const auto &TRI : TRIS) {
    E_Int A = TRI.A;
    E_Int B = TRI.B;
    E_Int C = TRI.C;

    fprintf(fh, SF_D3_ "\n", pmap[A], pmap[B], pmap[C]);
  }

  fclose(fh);
}
