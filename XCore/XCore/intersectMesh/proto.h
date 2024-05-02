#ifndef _INTERSECT_MESH_PROTO_H
#define _INTERSECT_MESH_PROTO_H

#include "struct.h"
#include <set>
#include <map>
#include <unordered_map>

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
Mesh *mesh_make_surface_mesh_from_structured_points(E_Int ni, E_Int nj,
  const std::vector<E_Float> &X, const std::vector<E_Float> &Y,
  const std::vector<E_Float> &Z);
void mesh_patch_intersect(Mesh *M, Mesh *S, E_Int *mpatch, E_Int mpatchc,
  E_Int *spatch, E_Int spatchc);

// Io
void point_set_write(const std::set<E_Int> &points, Mesh *M, const char *fname);
void edge_hits_write(const std::map<Edge_NO, Edge_Hit> &edge_hit_table,
  Mesh *M, const char *ename, const char *pname);
void mesh_write(Mesh *M, const char *file_name);
void point_hits_write(const std::map<E_Int, Edge_Hit> &point_hit_table,
  Mesh *M, const char *ename);
void write_edge_path(Mesh *M, const Edge_NO &edge,
  const std::vector<E_Int> &path_faces, const char *fname);
void write_point_list(E_Float *X, E_Float *Y, E_Float *Z, E_Int *list,
  E_Int n, const char *fname);
void point_hits_write(const std::unordered_map<E_Int, Edge_Hit> &ptable,
  const char *fname);
void write_trimesh(const std::vector<Triangle> &TRIS, Mesh *M);


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

void event_print(event *);
void snode_print(snode *);

int vertex_cmp(vertex *, vertex *);
int vertex_cmp_xy(vertex *, vertex *);
int vertex_orient(vertex *, vertex *, vertex *);
int vertex_cmp_cwise(vertex *, vertex *, vertex *);

int sign(double);
int sign(int);
int xy_cmp(double, double, double, double);

int segment_cmp_lexico(segment *, segment *);
void segment_sort(std::vector<segment *> &S, int (*cmp)(segment *, segment *),
    int, int);
int segments_are_colli(segment *, segment *);

std::vector<hedge *> dcel_get_incident_hedges(vertex *);
void dcel_resolve(vertex *, std::vector<segment *> &L,
    std::vector<segment *> &C, std::vector<segment *> &U,
    std::vector<hedge *> &);
int hedge_check_without_faces(const std::vector<hedge *> &);
void dcel_set_cycles_inout(std::vector<cycle *> &);
std::vector<cycle *> dcel_make_cycles(const std::vector<hedge *> &);
void dcel_read_two(const char *, dcel &, dcel &);
int dcel_check_faces(const std::vector<hedge *> &,
    const std::vector<face *> &);
vertex *dcel_get_leftmost_vertex_of_outer_cycle(face *);
face *dcel_get_face_from_cycle(cycle *);
cycle *dcel_get_outer_cycle_from_cycle(cycle *);
std::vector<hedge *> dcel_get_face_outer_loop(face *);
int dcel_faces_are_adjacent(face *, face *);
hedge *dcel_get_hedge_of_color_from_face(face *, int);
int dcel_check_hedges_without_faces(const std::vector<hedge *> &);


void hedge_sort_ccwise(std::vector<hedge *> &, int, int);
void hedge_sort_cwise(std::vector<hedge *> &, int, int);
int hedge_cmp_cwise(hedge *, hedge *);
int hedges_incident_overlap(hedge *h, hedge *w);

void sweep(std::vector<segment *> &S, std::vector<vertex *> &V,
    std::vector<hedge *> &H, queue &Q, status &T);

const char *color_to_str(int);





#endif
