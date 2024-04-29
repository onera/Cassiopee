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
  E_Float dx, E_Float dy, E_Float dz,
  E_Float ax, E_Float ay, E_Float az,
  E_Float bx, E_Float by, E_Float bz,
  E_Float cx, E_Float cy, E_Float cz,
  Edge_Hit &EH)
{
  E_Float ray_vector[3] = {dx, dy, dz};

  // Find vectors of two edges vertex a
  E_Float e1[3] = {bx-ax, by-ay, bz-az};
  E_Float e2[3] = {cx-ax, cy-ay, cz-az};
  E_Float ray_cross_e2[3];
  K_MATH::cross(ray_vector, e2, ray_cross_e2);
  E_Float det = K_MATH::dot(e1, ray_cross_e2, 3);

  if (det > -TOL && det < TOL)
    return 0; // This ray is parallel to this triangle

  E_Float inv_det = 1.0 / det;
  E_Float s[3] = {px-ax, py-ay, pz-az};
  EH.v = inv_det * K_MATH::dot(s, ray_cross_e2, 3);

  if (EH.v < -TOL || EH.v > 1+TOL)
    return 0;

  E_Float s_cross_e1[3];
  K_MATH::cross(s, e1, s_cross_e1);
  EH.w = inv_det * K_MATH::dot(ray_vector, s_cross_e1, 3);

  if (EH.w < -TOL || EH.w + EH.v > 1+TOL)
    return 0;

  EH.u = 1.0 - EH.w - EH.v;
  assert(EH.u > -TOL && EH.u < 1+TOL);

  // At this stage we can compute t to find out where the intersection point
  // is on the line
  EH.t = inv_det * K_MATH::dot(e2, s_cross_e1, 3);

  if (EH.t > TOL) {
    EH.x = px + dx * EH.t;
    EH.y = py + dy * EH.t;
    EH.z = pz + dz * EH.t;
    return 1;
  }

  return 0;
}