#ifndef _INTERSECT_MESH_PROTO_H
#define _INTERSECT_MESH_PROTO_H

#include "struct.h"
#include <set>
#include <map>

// Mesh
Mesh *mesh_init(K_FLD::FldArrayI &cn, E_Float *X, E_Float *Y, E_Float *Z,
  E_Int np);
E_Int *mesh_get_cell(E_Int cell, E_Int &stride, Mesh *M);
E_Int *mesh_get_face(E_Int face, E_Int &stride, Mesh *M);
E_Int mesh_get_stride(E_Int elem, E_Int *ind);
Mesh *mesh_extract_from_cell_set(Mesh *M, const std::set<E_Int> &cell_set);
std::vector<E_Int> mesh_get_external_faces_indices(Mesh *M);
Mesh *mesh_make_surface_mesh_from_face_list(Mesh *M,
  const std::vector<E_Int> &face_list);
std::vector<E_Int> mesh_make_cell_points(E_Int cell, Mesh *M);
std::vector<std::vector<E_Int>> mesh_make_point_cells(Mesh *M);
E_Int mesh_get_neighbour(E_Int cell, E_Int face, Mesh *M);
E_Int mesh_orient_boundary(Mesh *M);
E_Int mesh_build_own_nei(Mesh *M);


// Io
void point_set_write(const std::set<E_Int> &points, Mesh *M, const char *fname);
void edge_hits_write(const std::map<Edge_NO, Edge_Hit> &edge_hit_table,
  Mesh *M, const char *ename, const char *pname);
void mesh_write(Mesh *M, const char *file_name);
void point_hits_write(const std::map<E_Int, Edge_Hit> &point_hit_table,
  Mesh *M, const char *ename);
void write_edge_path(Mesh *M, const Edge_NO &edge,
  const std::vector<E_Int> &path_faces, const char *fname);


// BBox
BBox bbox_make(Mesh *M);
void bbox_print(const BBox &bbox);
BBox bbox_union(const BBox &b1, const BBox &b2);
BBox bbox_from_two_meshes(Mesh *A, Mesh *B);
void bbox_normalize_mesh_xyz(const BBox &bbox, Mesh *M);
E_Int bbox_is_point_inside(const BBox &bbox, E_Float x, E_Float y, E_Float z);

// Geom
E_Int geom_ray_triangle_intersect(E_Float px, E_Float py, E_Float pz,
  E_Float qx, E_Float qy, E_Float qz,
  E_Float ax, E_Float ay, E_Float az,
  E_Float bx, E_Float by, E_Float bz,
  E_Float cx, E_Float cy, E_Float cz,
  Edge_Hit &EH);
void geom_compute_face_normal(E_Int face, Mesh *M, E_Float n[3]);
E_Int geom_is_right(E_Float ax, E_Float ay, E_Float az, E_Float bx, E_Float by,
  E_Float bz, E_Float cx, E_Float cy, E_Float cz);

#endif