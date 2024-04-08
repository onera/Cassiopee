#include "proto.h"

int geom_sign(E_Float x)
{
  if (fabs(x) < TOL) return 0;
  else if (x >= TOL) return 1;
  else {
    assert(x <= -TOL);
    return -1;
  }
}

int FEQ(E_Float x, E_Float y)
{
  return fabs(x-y) < TOL;
}

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

void geom_barycentric(E_Float px, E_Float py, E_Float pz,
  E_Float ax, E_Float ay, E_Float az,
  E_Float bx, E_Float by, E_Float bz,
  E_Float cx, E_Float cy, E_Float cz,
  E_Float &u, E_Float &v, E_Float &w)
{
  E_Float v0[3] = {bx-ax, by-ay, bz-az};
  E_Float v1[3] = {cx-ax, cy-ay, cz-az};
  E_Float v2[3] = {px-ax, py-ay, pz-az};

  E_Float d00 = K_MATH::dot(v0, v0, 3);
  E_Float d01 = K_MATH::dot(v0, v1, 3);
  E_Float d11 = K_MATH::dot(v1, v1, 3);
  E_Float d20 = K_MATH::dot(v2, v0, 3);
  E_Float d21 = K_MATH::dot(v2, v1, 3);
  E_Float denom = d00 * d11 - d01 * d01;
  v = (d11 * d20 - d01 * d21) / denom;
  w = (d00 * d21 - d01 * d20) / denom;
  u = 1.0 - v - w;
}

int geom_orient_2D(const point &a, const point &b, const point &c)
{
  E_Float det = (b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x);
  if (det > TOL) return LEFT;
  if (det < -TOL) return RIGHT;
  return ON;
}

int geom_orient_2D(const point *a, const point *b, const point *c)
{
  E_Float det = (b->x-a->x)*(c->y-a->y) - (b->y-a->y)*(c->x-a->x);
  if (det > TOL) return LEFT;
  if (det < -TOL) return RIGHT;
  return ON;
}

int geom_line_intersect_2D(const point &A, const point &B, const point &C,
  const point &D, point &I)
{
  double s, t;
  double num, denom;

  denom = A.x * (D.y - C.y) +
          B.x * (C.y - D.y) +
          D.x * (B.y - A.y) +
          C.x * (A.y - B.y);
  
  assert(fabs(denom) >= TOL);

  num = A.x * (D.y - C.y) +
        C.x * (A.y - D.y) +
        D.x * (C.y - A.y);
  
  s = num / denom;

  num = A.x * (C.y - B.y) +
        B.x * (A.y - C.y) +
        C.x * (B.y - A.y);
  
  t = -num / denom;

  I.x = A.x + s * (B.x - A.x);
  I.y = A.y + s * (B.y - A.y);

  if (s <= -TOL || s >= 1.0+TOL || t <= -TOL || t >= 1.0+TOL)
    return 0;

  return 1;
}

int geom_segment_intersect_2D(segment *S1, segment *S2, point &I, point &J)
{
  const point &A = *(S1->p);
  const point &B = *(S1->q);
  const point &P = *(S2->p);
  const point &Q = *(S2->q);

  // Position of P and Q wrt [A,B]
  int sign_P = geom_orient_2D(A, B, P);
  int sign_Q = geom_orient_2D(A, B, Q);
  int mult_sign = sign_P * sign_Q;

  // P and Q on the same (strictly positive/negative) side of AB
  if (mult_sign == 1)
    return 0;

  // P and Q on different non-zero sides of AB
  else if (mult_sign == -1)
    return geom_line_intersect_2D(A, B, P, Q, I);

  // P not on AB and (A,B,Q) collinear: potential intersection at Q
  else if (abs(sign_P) == 1 && sign_Q == 0) {
    if (A.x <= Q.x && Q.x <= B.x) {
      I.x = Q.x;
      I.y = Q.y;
      return 1;
    } else {
      return 0;
    }
  }

  // Q not on AB and (A,B,P) collinear: potential intersection at P
  else if (abs(sign_Q) == 1 && sign_P == 0) {
    if (A.x <= P.x && P.x <= B.x) {
      I.x = P.x;
      I.y = P.y;
      return 1;
    } else {
      return 0;
    }
  }

  // AB and PQ collinear
  else if (sign_P == 0 && sign_Q == 0) {
    if (Q.x < A.x || P.x > B.x)
      return 0;
    
    if (P.x < A.x && B.x < Q.x) {
      I.x = A.x;
      I.y = A.y;
      J.x = B.x;
      J.y = B.y;
      return 2;
    }

    if (A.x < P.x && Q.x < B.x) {
      I.x = P.x;
      I.y = P.y;
      J.x = Q.x;
      J.y = Q.x;
      return 2;
    }

    if (P.x < A.x && (Q.x > A.x && Q.x <= B.x)) {
      I.x = Q.x;
      I.y = Q.y;
      J.x = A.x;
      J.y = A.y;
      return 2;
    }

    if (Q.x > B.x && (P.x >= A.x && P.x < B.x)) {
      I.x = P.x;
      I.y = P.y;
      J.x = B.x;
      J.y = B.y;
      return 2;
    }

    return 0;

  }

  // In case I missed a position...
  else {
    assert(0);
  }

  assert(0);
  return 0;
}
