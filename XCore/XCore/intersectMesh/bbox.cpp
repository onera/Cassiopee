#include "proto.h"

BBox bbox_make(Mesh *M)
{
  BBox bbox;
  bbox.xmin = E_FLOAT_MAX;
  bbox.xmax = E_FLOAT_MIN;
  bbox.ymin = E_FLOAT_MAX;
  bbox.ymax = E_FLOAT_MIN;
  bbox.zmin = E_FLOAT_MAX;
  bbox.zmax = E_FLOAT_MIN;

  for (E_Int i = 0; i < M->np; i++) {
    if (M->x[i] < bbox.xmin) bbox.xmin = M->x[i];
    if (M->x[i] > bbox.xmax) bbox.xmax = M->x[i];
    if (M->y[i] < bbox.ymin) bbox.ymin = M->y[i];
    if (M->y[i] > bbox.ymax) bbox.ymax = M->y[i];
    if (M->z[i] < bbox.zmin) bbox.zmin = M->z[i];
    if (M->z[i] > bbox.zmin) bbox.zmax = M->z[i];
  }

  bbox.dmax = std::max(bbox.xmax-bbox.xmin, bbox.ymax-bbox.ymin);
  bbox.dmax = std::max(bbox.dmax, bbox.zmax-bbox.zmin);

  assert(bbox.dmax >= 0.0);

  return bbox;
}

void bbox_print(const BBox &bbox)
{
  printf("x: [%.4e %.4e] - y: [%.4e %.4e] - z: [%.4e %.4e]\n\n",
    bbox.xmin, bbox.xmax, bbox.ymin, bbox.ymax, bbox.zmin, bbox.zmax);
}

E_Int bbox_is_point_inside(const BBox &bbox, E_Float x, E_Float y, E_Float z)
{
  return (x >= bbox.xmin && x <= bbox.xmax) &&
         (y >= bbox.ymin && y <= bbox.ymax) &&
         (z >= bbox.zmin && z <= bbox.zmax);
}

BBox bbox_from_two_meshes(Mesh *A, Mesh *B)
{
  BBox bbox;
  bbox.xmin = E_FLOAT_MAX;
  bbox.xmax = E_FLOAT_MIN;
  bbox.ymin = E_FLOAT_MAX;
  bbox.ymax = E_FLOAT_MIN;
  bbox.zmin = E_FLOAT_MAX;
  bbox.zmax = E_FLOAT_MIN;

  for (E_Int i = 0; i < A->np; i++) {
    if (A->x[i] < bbox.xmin) bbox.xmin = A->x[i];
    if (A->x[i] > bbox.xmax) bbox.xmax = A->x[i];
    if (A->y[i] < bbox.ymin) bbox.ymin = A->y[i];
    if (A->y[i] > bbox.ymax) bbox.ymax = A->y[i];
    if (A->z[i] < bbox.zmin) bbox.zmin = A->z[i];
    if (A->z[i] > bbox.zmax) bbox.zmax = A->z[i];
  }

  for (E_Int i = 0; i < B->np; i++) {
    if (B->x[i] < bbox.xmin) bbox.xmin = B->x[i];
    if (B->x[i] > bbox.xmax) bbox.xmax = B->x[i];
    if (B->y[i] < bbox.ymin) bbox.ymin = B->y[i];
    if (B->y[i] > bbox.ymax) bbox.ymax = B->y[i];
    if (B->z[i] < bbox.zmin) bbox.zmin = B->z[i];
    if (B->z[i] > bbox.zmax) bbox.zmax = B->z[i];
  }

  bbox.dmax = std::max(bbox.xmax-bbox.xmin, bbox.ymax-bbox.ymin);
  bbox.dmax = std::max(bbox.dmax, bbox.zmax-bbox.zmin);

  return bbox;
}

BBox bbox_union(const BBox &b1, const BBox &b2)
{
  BBox bbox;

  bbox.xmin = std::min(b1.xmin, b2.xmin);
  bbox.xmax = std::max(b1.xmax, b2.xmax);
  bbox.ymin = std::min(b1.ymin, b2.ymin);
  bbox.ymax = std::max(b1.ymax, b2.ymax);
  bbox.zmin = std::min(b1.zmin, b2.zmin);
  bbox.zmax = std::max(b1.zmax, b2.zmax);

  bbox.dmax = std::max(bbox.xmax-bbox.xmin, bbox.ymax-bbox.ymin);
  bbox.dmax = std::max(bbox.dmax, bbox.zmax-bbox.zmin);

  return bbox;
}

void bbox_normalize_mesh_xyz(const BBox &bbox, Mesh *M)
{
  E_Float inv_dmax = 1. / bbox.dmax;

  for (E_Int i = 0; i < M->np; i++) {
    M->x[i] = (M->x[i] - bbox.xmin) * inv_dmax;
    M->y[i] = (M->y[i] - bbox.ymin) * inv_dmax;
    M->z[i] = (M->z[i] - bbox.zmin) * inv_dmax;
  }
}