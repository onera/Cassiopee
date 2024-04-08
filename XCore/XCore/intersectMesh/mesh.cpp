#include "proto.h"
#include <map>
#include <unordered_map>

#define MAXPOINTSPERCELL 8

Mesh::Mesh() :
  nc(0), nf(0), np(0), ne(0),
  nface(NULL), indPH(NULL), ngon(NULL), indPG(NULL),
  x(NULL), y(NULL), z(NULL),
  owner(NULL), neigh(NULL)
{}

Mesh *mesh_init(K_FLD::FldArrayI &cn, E_Float *X, E_Float *Y, E_Float *Z,
  E_Int np)
{
  Mesh *M = new Mesh;

  M->nc = cn.getNElts();
  M->nf = cn.getNFaces();
  M->np = np;

  // Coordinates
  M->x = (E_Float *)XMALLOC(M->np * sizeof(E_Float));
  M->y = (E_Float *)XMALLOC(M->np * sizeof(E_Float));
  M->z = (E_Float *)XMALLOC(M->np * sizeof(E_Float));

  memcpy(M->x, X, M->np * sizeof(E_Float));
  memcpy(M->y, Y, M->np * sizeof(E_Float));
  memcpy(M->z, Z, M->np * sizeof(E_Float));

  // Connectivities
  E_Int *indPH = cn.getIndPH();

  E_Int *indPG = cn.getIndPG();

  M->indPH = (E_Int *)XMALLOC((M->nc+1) * sizeof(E_Int));
  M->indPG = (E_Int *)XMALLOC((M->nf+1) * sizeof(E_Int));

  for (E_Int i = 0; i < M->nc+1; i++) M->indPH[i] = indPH[i];
  for (E_Int i = 0; i < M->nf+1; i++) M->indPG[i] = indPG[i];


  M->nface = (E_Int *)XMALLOC(cn.getSizeNFace() * sizeof(E_Int));
  M->ngon  = (E_Int *)XMALLOC(cn.getSizeNGon()  * sizeof(E_Int));

  E_Int *nface = cn.getNFace();  
  E_Int *ngon = cn.getNGon();
  

  for (E_Int i = 0; i < M->nc; i++) {
    E_Int nf = -1;
    E_Int *pf = cn.getElt(i, nf, nface, indPH);

    E_Int *ptr = &M->nface[M->indPH[i]];
    for (E_Int j = 0; j < nf; j++) ptr[j] = pf[j]-1;
  }

  for (E_Int i = 0; i < M->nf; i++) {
    E_Int np = -1;
    E_Int *pn = cn.getFace(i, np, ngon, indPG);

    E_Int *ptr = &M->ngon[M->indPG[i]];
    for (E_Int j = 0; j < np; j++) ptr[j] = pn[j]-1;
  }

  M->owner = (E_Int *)XMALLOC(M->nf * sizeof(E_Int));
  M->neigh = (E_Int *)XMALLOC(M->nf * sizeof(E_Int));

  //mesh_orient_boundary(M);
  //mesh_build_own_nei(M);

  return M;
}

E_Int *mesh_get_cell(E_Int cell, E_Int &stride, Mesh *M)
{
  stride = M->indPH[cell+1] - M->indPH[cell];
  return &M->nface[M->indPH[cell]]; 
}

E_Int *mesh_get_face(E_Int face, E_Int &stride, Mesh *M)
{
  stride = M->indPG[face+1] - M->indPG[face];
  return &M->ngon[M->indPG[face]];
}

std::vector<E_Int> mesh_make_cell_points(E_Int cell, Mesh *M)
{
  E_Int nf = -1;
  E_Int *pf = mesh_get_cell(cell, nf, M);

  // Max number of points is less than or equal to the number of faces
  std::vector<E_Int> cell_points(MAXPOINTSPERCELL);

  E_Int npoints = 0;

  // Add the points of the first face
  E_Int face = pf[0];
  E_Int np = -1;
  E_Int *pn = mesh_get_face(face, np, M);

  for (E_Int i = 0; i < np; i++) cell_points[npoints++] = pn[i];

  for (E_Int i = 1; i < nf; i++) {
    face = pf[i];
    np = -1;
    pn = mesh_get_face(face, np, M);

    // Check for duplicate points
    for (E_Int j = 0; j < np; j++) {
      E_Int cur_point = pn[j];

      E_Int found_point = 0;

      for (E_Int k = 0; k < npoints; k++) {
        if (cur_point == cell_points[k]) {
          found_point = 1;
          break;
        }
      }

      if (!found_point) cell_points[npoints++] = cur_point;
    }
  }

  assert(npoints <= MAXPOINTSPERCELL);

  cell_points.resize(npoints);

  return cell_points;
}

std::vector<std::vector<E_Int>> mesh_make_point_cells(Mesh *M)
{
  std::vector<std::vector<E_Int>> point_cells(M->np);

  // Count the number of cells per points
  std::vector<E_Int> npc(M->np, 0);
  
  for (E_Int i = 0; i < M->nc; i++) {
    auto cell_points = mesh_make_cell_points(i, M);
    for (auto point : cell_points) npc[point] += 1;
  }

  // Resize
  for (E_Int i = 0; i < M->np; i++) point_cells[i].resize(npc[i]);

  // Fill
  for (E_Int i = 0; i < M->np; i++) npc[i] = 0;

  for (E_Int i = 0; i < M->nc; i++) {
    auto cell_points = mesh_make_cell_points(i, M);

    for (auto point : cell_points) point_cells[point][npc[point]++] = i;
  }

  return point_cells;
}

E_Int mesh_get_stride(E_Int cell, E_Int *ind)
{
  return ind[cell+1] - ind[cell];
}

Mesh *mesh_extract_from_cell_set(Mesh *M, const std::set<E_Int> &cell_set)
{
  Mesh *O = new Mesh;

  O->nc = cell_set.size();
  O->indPH = (E_Int *)XMALLOC((O->nc+1) * sizeof(E_Int));

  // Count sizeNFace
  E_Int sizeNFace = 0;

  E_Int nc = 0;

  for (auto cell : cell_set) {
    E_Int stride = mesh_get_stride(cell, M->indPH);
    sizeNFace += stride;
    O->indPH[++nc] = stride;
  }

  O->indPH[0] = 0;
  for (E_Int i = 0; i < nc; i++) O->indPH[i+1] += O->indPH[i];

  O->nface = (E_Int *)XMALLOC(sizeNFace * sizeof(E_Int));

  nc = 0;
  E_Int *ptr = O->nface;

  for (auto cell : cell_set) {
    for (E_Int j = M->indPH[cell]; j < M->indPH[cell+1]; j++)
      *ptr++ = M->nface[j];
  }

  assert(ptr - O->nface == sizeNFace);

  // Hash kept faces and count sizeNGon
  std::map<E_Int, E_Int> face_table;
  E_Int sizeNGon = 0;

  O->nf = 0;

  for (E_Int i = 0; i < O->nc; i++) {
    for (E_Int j = O->indPH[i]; j < O->indPH[i+1]; j++) {
      E_Int face = O->nface[j];

      auto search = face_table.find(face);

      if (search == face_table.end()) {
        face_table[face] = O->nf++;
        sizeNGon += mesh_get_stride(face, M->indPG);
      }
    }
  }

  O->indPG = (E_Int *)XMALLOC((M->nf+1) * sizeof(E_Int));
  O->ngon = (E_Int *)XMALLOC(sizeNGon * sizeof(E_Int));

  O->indPG[0] = 0;
  for (auto F : face_table) {
    E_Int old_face = F.first;
    E_Int new_face = F.second;

    O->indPG[new_face+1] = mesh_get_stride(old_face, M->indPG);
  }

  for (E_Int i = 0; i < O->nf+1; i++) O->indPG[i+1] += O->indPG[i];

  // Fill in the points
  std::map<E_Int, E_Int> point_table;

  O->np = 0;

  for (E_Int i = 0; i < O->nc; i++) {
    for (E_Int j = O->indPH[i]; j < O->indPH[i+1]; j++) {
      E_Int old_face = O->nface[j];
      E_Int new_face = face_table[old_face];

      E_Int *ptr = &O->ngon[O->indPG[new_face]];

      E_Int np = -1;
      E_Int *pn = mesh_get_face(old_face, np, M);

      for (E_Int k = 0; k < np; k++) {
        E_Int point = pn[k];

        *ptr++ = point;

        auto search = point_table.find(point);

        if (search == point_table.end()) {
          point_table[point] = O->np;
          O->np++;
        }
      }
    }
  }

  // Get point coordinates
  O->x = (E_Float *)XMALLOC(O->np * sizeof(E_Float));
  O->y = (E_Float *)XMALLOC(O->np * sizeof(E_Float));
  O->z = (E_Float *)XMALLOC(O->np * sizeof(E_Float));


  for (auto P : point_table) {
    E_Int old_point = P.first;
    E_Int new_point = P.second;

    O->x[new_point] = M->x[old_point];
    O->y[new_point] = M->y[old_point];
    O->z[new_point] = M->z[old_point];
  }

  // Set new faces/points in nface/ngon
  for (E_Int i = 0; i < sizeNFace; i++) O->nface[i] = face_table[O->nface[i]];
  for (E_Int i = 0; i < sizeNGon; i++) O->ngon[i] = point_table[O->ngon[i]];

  return O;
}

std::vector<E_Int> mesh_get_external_faces_indices(Mesh *M)
{
  // External faces appear only once in nface
  std::vector<E_Int> face_count(M->nf, 0);

  for (E_Int i = 0; i < M->indPH[M->nc]; i++)
    face_count[M->nface[i]]++;

  std::vector<E_Int> external_faces;

  for (E_Int i = 0; i < M->nf; i++) {
    if (face_count[i] == 1)
      external_faces.push_back(i);
  }

  return external_faces;
}

Mesh *mesh_make_surface_mesh_from_face_list(Mesh *M,
  const std::vector<E_Int> &face_list)
{
  Mesh *E = new Mesh;

  E->nc = face_list.size();
  E->indPH = (E_Int *)XMALLOC((E->nc+1) * sizeof(E_Int));
  E->indPH[0] = 0;

  // Count sizeNFace
  E_Int sizeNFace = 0;

  for (auto i = 0; i < face_list.size(); i++) {
    E_Int face = face_list[i];

    // A face has as many edges as it has points
    sizeNFace += mesh_get_stride(face, M->indPG);
    E->indPH[i+1] = sizeNFace;
  }

  E->nface = (E_Int *)XMALLOC(sizeNFace * sizeof(E_Int));

  // Populate nface with edges
  E_Int *ptr = E->nface;

  std::map<Edge_NO, E_Int> edge_table;
  std::unordered_map<E_Int, E_Int> point_table;

  E->nf = 0;
  E->np = 0;

  for (auto i = 0; i < face_list.size(); i++) {
    E_Int face = face_list[i];

    E_Int np = -1;
    E_Int *pn = mesh_get_face(face, np, M);

    // First pass: hash the points
    for (E_Int j = 0; j < np; j++) {
      E_Int p = pn[j];
      
      auto search = point_table.find(p);

      if (search == point_table.end()) {
        point_table[p] = E->np++;
      }
    }

    // Second pass: hash the edges
    for (E_Int j = 0; j < np; j++) {
      E_Int p = pn[j];
      E_Int q = pn[(j+1)%np];

      Edge_NO e(p,q);

      auto search = edge_table.find(e);

      if (search == edge_table.end()) {
        *ptr++ = E->nf;
        edge_table[e] = E->nf++;
      } else {
        *ptr++ = search->second;
      }
    }
  }

  // Make ngon
  E->indPG = (E_Int *)XMALLOC((E->nf+1) * sizeof(E_Int));
  E->indPG[0] = 0;
  for (E_Int i = 1; i < E->nf+1; i++) E->indPG[i] = 2;
  for (E_Int i = 0; i < E->nf; i++) E->indPG[i+1] += E->indPG[i];

  E_Int sizeNGon = E->indPG[E->nf];
  E->ngon = (E_Int *)XMALLOC(sizeNGon * sizeof(E_Int));

  // Populate ngon with points
  for (auto e : edge_table) {
    E_Int pos = e.second;
    E_Int *ptr = &E->ngon[2*pos];
    *ptr++ = point_table[e.first.p];
    *ptr++ = point_table[e.first.q];
  }

  // Make xyz
  E->x = (E_Float *)XMALLOC(E->np * sizeof(E_Float));
  E->y = (E_Float *)XMALLOC(E->np * sizeof(E_Float));
  E->z = (E_Float *)XMALLOC(E->np * sizeof(E_Float));

  for (auto P : point_table) {
    E_Int old_point = P.first;
    E_Int new_point = P.second;

    E->x[new_point] = M->x[old_point];
    E->y[new_point] = M->y[old_point];
    E->z[new_point] = M->z[old_point];
  }

  return E;
}

E_Int mesh_get_neighbour(E_Int cell, E_Int face, Mesh *M)
{
  assert(M->owner[face] == cell || M->neigh[face] == cell);
  if (M->owner[face] == cell) return M->neigh[face];
  return M->owner[face];
}

Mesh *mesh_make_surface_mesh_from_structured_points(E_Int ni, E_Int nj,
  const std::vector<E_Float> &X, const std::vector<E_Float> &Y,
  const std::vector<E_Float> &Z)
{
  assert(ni*nj == (E_Int)X.size());

  Mesh *M = new Mesh();

  M->np = ni*nj;
  M->nf = ni*(nj-1) + nj*(ni-1);
  M->nc = 2 + M->nf - M->np - 1;

  // Points
  M->x = (E_Float *)XMALLOC(M->np * sizeof(E_Float));
  M->y = (E_Float *)XMALLOC(M->np * sizeof(E_Float));
  M->z = (E_Float *)XMALLOC(M->np * sizeof(E_Float));

  memcpy(M->x, &X[0], M->np * sizeof(E_Float));
  memcpy(M->y, &Y[0], M->np * sizeof(E_Float));
  memcpy(M->z, &Z[0], M->np * sizeof(E_Float));

  // Edges
  M->indPG = (E_Int *)XMALLOC((M->nf+1) * sizeof(E_Int));
  M->indPG[0] = 0;
  for (E_Int i = 1; i < M->nf+1; i++) M->indPG[i] = 2;
  for (E_Int i = 0; i < M->nf; i++) M->indPG[i+1] += M->indPG[i];

  M->ngon = (E_Int *)XMALLOC(M->indPG[M->nf] * sizeof(E_Int));
  E_Int *ptr = M->ngon;

  for (E_Int j = 0; j < nj; j++) {
    for (E_Int i = 0; i < ni-1; i++) {
      E_Int p = i + ni*j;
      E_Int q = p + 1;
      *ptr++ = p;
      *ptr++ = q;
    }
  }

  for (E_Int i = 0; i < ni; i++) {
    for (E_Int j = 0; j < nj-1; j++) {
      E_Int p = i + ni*j;
      E_Int q = p + ni;
      *ptr++ = p;
      *ptr++ = q;
    }
  }

  // Faces
  M->indPH = (E_Int *)XMALLOC((M->nc+1) * sizeof(E_Int));
  M->indPH[0] = 0;
  for (E_Int i = 1; i < M->nc+1; i++) M->indPH[i] = 4;
  for (E_Int i = 0; i < M->nc; i++) M->indPH[i+1] += M->indPH[i];

  M->nface = (E_Int *)XMALLOC(M->indPH[M->nc] * sizeof(E_Int));

  E_Int nie = nj*(ni-1);

  // Face is bounded from the:
  // Bottom by edge i+(ni-1)*j,
  // Top by edge i+(ni-1)*(j+1),
  // Left by edge nie + j + (nj-1)*i
  // Right by edge nie + j + (nj-1)*(i+1)

  ptr = M->nface;

  for (E_Int j = 0; j < nj-1; j++) {
    for (E_Int i = 0; i < ni-1; i++) {
      *ptr++ = i + (ni-1)*j;
      *ptr++ = i + (ni-1)*(j+1);
      *ptr++ = j + (nj-1)*i + nie;
      *ptr++ = j + (nj-1)*(i+1) + nie;
    }
  }

  return M;
}

/*
Mesh *mesh_make_face_point_mesh_from_face_list(Mesh *M,
  const std::vector<E_Int> &face_list)
{
  Mesh *O = new Mesh;

  return O;
}
*/
