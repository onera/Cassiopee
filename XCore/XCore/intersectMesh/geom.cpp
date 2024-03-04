#include "proto.h"

Edge_NO::Edge_NO(E_Int p_, E_Int q_)
{
  if (p_ < q_) {
    p = p_;
    q = q_;
  } else {
    p = q_;
    q = p_;
  }
}

bool Edge_NO::operator<(const Edge_NO &other) const
{
  return (p < other.p) || (p == other.p && q < other.q);
}

// Source: Moller-Trumbore algorithm
E_Int geom_ray_triangle_intersect(E_Float px, E_Float py, E_Float pz,
  E_Float qx, E_Float qy, E_Float qz,
  E_Float ax, E_Float ay, E_Float az,
  E_Float bx, E_Float by, E_Float bz,
  E_Float cx, E_Float cy, E_Float cz,
  Edge_Hit &EH)
{
  const E_Float TOL = 1.e-6;

  // Compute (normalized) edge direction
  E_Float dir[3] = {qx-px, qy-py, qz-pz};

  // Find vectors of two edges vertex a
  E_Float e1[3] = {bx-ax, by-ay, bz-az};
  E_Float e2[3] = {cx-ax, cy-ay, cz-az};

  // Begin calculating determinant
  E_Float pvec[3];
  K_MATH::cross(dir, e2, pvec);

  // If determinant is near zero, ray lies in plane of triangle
  E_Float det = K_MATH::dot(e1, pvec, 3);

  if (det > -TOL && det < TOL)
    return 0;
  
  E_Float inv_det = 1.0 / det;

  // Calculate distance from vertex a to ray origin
  E_Float tvec[3] = {px-ax, py-ay, pz-az};

  // Calculate u parameter and test bounds
  EH.u = K_MATH::dot(tvec, pvec, 3) * inv_det;

  if (EH.u < 0.0 || EH.u > 1.0)
    return 0;
  
  // Prepare to test v parameter
  K_MATH::cross(tvec, e1, pvec);

  // Calculate v parameter and test bounds
  EH.v = K_MATH::dot(dir, pvec, 3) * inv_det;

  if (EH.v < 0.0 || (EH.u + EH.v > 1.0))
    return 0;
  
  // Calculate t, ray intersects triangle
  EH.t = K_MATH::dot(e2, pvec, 3) * inv_det;

  if (EH.t > TOL) {
    EH.x = px + dir[0] * EH.t;
    EH.y = py + dir[1] * EH.t;
    EH.z = pz + dir[2] * EH.t;

    return 1;
  }

  return 0;
}

void geom_compute_face_normal(E_Int face, Mesh *M, E_Float n[3])
{
  assert(face < M->nf);

  E_Int np = -1;
  E_Int *pn = mesh_get_face(face, np, M);

  assert(np >= 3);

  E_Int p0 = pn[0];
  E_Int p1 = pn[1];
  E_Int p2 = pn[2];

  E_Float e0[3] = {M->x[p1]-M->x[p0], M->y[p1]-M->y[p0], M->z[p1]-M->z[p0]};
  E_Float e1[3] = {M->x[p2]-M->x[p0], M->y[p2]-M->y[p0], M->z[p2]-M->z[p0]};

  K_MATH::cross(e0, e1, n);

  E_Float inv_norm = 1.0 / K_MATH::norm(n, 3);

  for (E_Int i = 0; i < 3; i++) n[i] *= inv_norm;
}

E_Int geom_is_right(E_Float ax, E_Float ay, E_Float az, E_Float bx, E_Float by,
  E_Float bz, E_Float cx, E_Float cy, E_Float cz)
{
  return (bx-ax)*(cy-ay) - (by-ay)*(cx-ax) > 0.0;
}