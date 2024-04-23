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
void geom_compute_face_normal(E_Int face, Mesh *M, E_Float n[3]);
E_Int geom_is_right(E_Float ax, E_Float ay, E_Float az, E_Float bx, E_Float by,
  E_Float bz, E_Float cx, E_Float cy, E_Float cz);
void geom_barycentric(E_Float px, E_Float py, E_Float pz,
  E_Float ax, E_Float ay, E_Float az,
  E_Float bx, E_Float by, E_Float bz,
  E_Float cx, E_Float cy, E_Float cz,
  E_Float &u, E_Float &v, E_Float &w);
E_Int geom_segment_intersect_2D(
  E_Float Ax, E_Float Ay,
  E_Float Bx, E_Float By,
  E_Float Px, E_Float Py,
  E_Float Qx, E_Float Qy,
  Edge_Hit &EH1, Edge_Hit &EH2);
int geom_orient_2D(const point *a, const point *b, const point *c);
int geom_orient_2D(const point &a, const point &b, const point &c);
int geom_line_intersect_2D(const point &A, const point &B, const point &C,
  const point &D, point &I);
int geom_segment_intersect_2D(segment *S1, segment *S2, point &I, point &J);

// Math
int FEQ(E_Float a, E_Float b);

// DCEL
dcel *dcel_make_from_mesh(Mesh *M);
dcel *dcel_make_from_structured_points(E_Int ni, E_Int nj,
  const std::vector<E_Float> &X, const std::vector<E_Float> &Y,
  const std::vector<E_Float> &Z);
dcel *dcel_overlay(dcel *D1, dcel *D2);
E_Int dcel_check(dcel *D);

int segment_is_lexico(const segment *S);

int sign(E_Float x);
int f_eq(E_Float x, E_Float y);
int cmp_points(const point *p, const point *q);
int cmp_xyz(E_Float x1, E_Float y1, E_Float x2, E_Float y2);
int segments_are_colli(segment *s1, segment *s2);
int get_orient(point *a, point *b, point *c);
int cmp_segments_lexico(const segment *s1, const segment *s2);

// Sweep
void _sweep(std::vector<segment *> &segs, std::vector<point *> &points,
    std::unordered_map<point *, std::vector<segment *>> & Xsegs);


/* Event */

event *event_insert(event *&, point *);
event *event_lookup(event *, point *);
event *event_locate(event *, E_Float, E_Float);
event *event_min(event *);
event *event_delete(event *, point *);
void   event_show(event *);
void   event_print(event *);

/* Status */

status *status_insert(status *&, segment *, E_Float, E_Float);
status *status_lookup(status *, segment *, E_Float, E_Float);
status *status_locate(status *, segment *, E_Float, E_Float);
status *status_pred(status *, segment *, E_Float, E_Float);
status *status_succ(status *, segment *, E_Float, E_Float);
status *status_delete(status *, segment *, E_Float, E_Float);
void    status_show(status *);
void    status_print(status *);


#endif
