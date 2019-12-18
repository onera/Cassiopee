/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <float.h>
#include <math.h>
#include <stdlib.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_triangle.h"
#include "pdm_line.h"
#include "pdm_plane.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/


/**
 * \brief Evaluates the position in a triangle
 *
 * \param [in]  x        Point coordinates to evaluate position
 * \param [in]  pts      Triangle vertices coordinates
 * \param [out] closest  Closest Point in Triangle or NULL
 * \param [out] minDist2 Square of the distance
 * \param [out] weights  Vertices weights or NULL
 *
 * \return      \ref PDM_TRIANGLE_INSIDE or \ref PDM_TRIANGLE_OUTSIDE
 *              if the projected is in the triangle or not
 *
 */

/*  This function is derived from VTK                                      */
/*  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen                */
/*  All rights reserved.                                                   */
/*  See Copyright.txt or http://www.kitware.com/Copyright.htm for details. */


PDM_triangle_status_t
PDM_triangle_evaluate_position
(
 const double  x[3],
 const double  pts[9],
       double *closestPoint,
       double *minDist2,
       double *weights
)
{
  double *pt1, *pt2, *pt3;
  double n[3], fabsn;
  double rhs[2], c1[2], c2[2];
  double det;
  double maxComponent;
  int idx=0, indices[2];
  double dist2Point, dist2Line1, dist2Line2;
  double *closest, closestPoint1[3], closestPoint2[3], cp[3];
  double pcoords[3];

  pcoords[2] = 0.0;

  double weights_local[3];
  double *_weights = weights_local;
  if (weights != NULL) {
    _weights = weights;
  }

  /*
   * Get normal for triangle, only the normal direction is needed, i.e. the
   * normal need not be normalized (unit length)
   */

  PDM_plane_normal (3, pts, n);

  pt1 = (double *) pts;
  pt2 = (double *) pts + 3;
  pt3 = (double *) pts + 6;

  /*
   * Project point to plane
   */

  PDM_plane_projection2 (x, pt1, n, cp);

  /*
   * Construct matrices.  Since we have over determined system, need to find
   * which 2 out of 3 equations to use to develop equations. (Any 2 should
   * work since we've projected point to plane.)
   */

  maxComponent = 0.0;
  for (int i = 0; i < 3; i++) {

    /*
     * Trying to avoid an expensive call to fabs()
     */

    if (n[i] < 0) {
      fabsn = -n[i];
    }
    else {
      fabsn = n[i];
    }
    if (fabsn > maxComponent) {
      maxComponent = fabsn;
      idx = i;
    }
  }

  for (int j = 0, i = 0; i < 3; i++) {
    if (i != idx) {
      indices[j++] = i;
    }
  }

  for (int i = 0; i < 2; i++) {
    rhs[i] = cp[indices[i]] - pt3[indices[i]];
    c1[i] = pt1[indices[i]] - pt3[indices[i]];
    c2[i] = pt2[indices[i]] - pt3[indices[i]];
  }

  if ((det = PDM_DETERMINANT2X2(c1,c2)) == 0.0) {
    pcoords[0] = 0.0;
    pcoords[1] = 0.0;
    return PDM_TRIANGLE_DEGENERATED;
  }

  pcoords[0] = PDM_DETERMINANT2X2(rhs,c2) / det;
  pcoords[1] = PDM_DETERMINANT2X2(c1,rhs) / det;

  /*
   * Okay, now find closest point to element
   */

  _weights[0] = 1 - (pcoords[0] + pcoords[1]);
  _weights[1] = pcoords[0];
  _weights[2] = pcoords[1];

  if ( _weights[0] >= 0.0 && _weights[0] <= 1.0 &&
       _weights[1] >= 0.0 && _weights[1] <= 1.0 &&
       _weights[2] >= 0.0 && _weights[2] <= 1.0 ) {

    /*
     * Projection distance
     */

    if (closestPoint) {
      double v_cp_x[3];
      for (int i = 0; i < 3; i++) {
        v_cp_x[i] = cp[i] - x[i];
      }

      *minDist2 = PDM_DOT_PRODUCT(v_cp_x, v_cp_x);
      closestPoint[0] = cp[0];
      closestPoint[1] = cp[1];
      closestPoint[2] = cp[2];
    }
    return PDM_TRIANGLE_INSIDE;
  }

  else {
    double t;
    if (closestPoint) {
      if (_weights[1] < 0.0 && _weights[2] < 0.0) {
        double v_pt3_x[3];
        for (int i = 0; i < 3; i++) {
          v_pt3_x[i] = pt3[i] - x[i];
        }
        dist2Point = PDM_DOT_PRODUCT (v_pt3_x, v_pt3_x);
        dist2Line1 = PDM_line_distance (x, pt1, pt3, &t, closestPoint1);
        dist2Line2 = PDM_line_distance (x, pt3, pt2, &t, closestPoint2);
        if (dist2Point < dist2Line1) {
          *minDist2 = dist2Point;
          closest = pt3;
        }
        else {
          *minDist2 = dist2Line1;
          closest = closestPoint1;
        }
        if (dist2Line2 < *minDist2) {
          *minDist2 = dist2Line2;
          closest = closestPoint2;
        }
        for (int i = 0; i < 3; i++) {
          closestPoint[i] = closest[i];
        }

      }
      else if (_weights[2] < 0.0 && _weights[0] < 0.0) {
        double v_pt1_x[3];
        for (int i = 0; i < 3; i++) {
          v_pt1_x[i] = pt1[i] - x[i];
        }
        dist2Point = PDM_DOT_PRODUCT(v_pt1_x, v_pt1_x);
        dist2Line1 = PDM_line_distance (x, pt1, pt3, &t, closestPoint1);
        dist2Line2 = PDM_line_distance (x, pt1, pt2, &t, closestPoint2);
        if (dist2Point < dist2Line1) {
          *minDist2 = dist2Point;
          closest = pt1;
        }
        else {
          *minDist2 = dist2Line1;
          closest = closestPoint1;
        }
        if (dist2Line2 < *minDist2) {
          *minDist2 = dist2Line2;
          closest = closestPoint2;
        }
        for (int i = 0; i < 3; i++) {
          closestPoint[i] = closest[i];
        }

      }
      else if ( _weights[1] < 0.0 && _weights[0] < 0.0 ) {
        double v_pt2_x[3];
        for (int i = 0; i < 3; i++) {
          v_pt2_x[i] = pt2[i] - x[i];
        }
        dist2Point = PDM_DOT_PRODUCT (v_pt2_x, v_pt2_x);
        dist2Line1 = PDM_line_distance (x, pt2, pt3, &t, closestPoint1);
        dist2Line2 = PDM_line_distance (x, pt1, pt2, &t, closestPoint2);
        if (dist2Point < dist2Line1) {
          *minDist2 = dist2Point;
          closest = pt2;
        }
        else {
          *minDist2 = dist2Line1;
          closest = closestPoint1;
        }
        if (dist2Line2 < *minDist2) {
          *minDist2 = dist2Line2;
          closest = closestPoint2;
        }
        for (int i = 0; i < 3; i++) {
          closestPoint[i] = closest[i];
        }

      }
      else if (_weights[0] < 0.0) {
        *minDist2 = PDM_line_distance (x, pt1, pt2, &t, closestPoint);

      }
      else if (_weights[1] < 0.0) {
        *minDist2 = PDM_line_distance (x, pt2, pt3, &t, closestPoint);

      }
      else if (_weights[2] < 0.0) {
        *minDist2 = PDM_line_distance (x, pt1, pt3, &t, closestPoint);

      }

    }
    return PDM_TRIANGLE_OUTSIDE;
  }
}


/**
 * \brief Computes polygon barycenter
 *
 * \param [in]   numPts  Number of polygon vertices
 * \param [in]   pts     Polygon vertices coordinates
 * \param [out]  bary    Barycenter
 *
 */

void
PDM_triangle_compute_barycenter
(
 const double pts[9],
       double bary[3]
)
{
  bary[0] = 0.;
  bary[1] = 0.;
  bary[2] = 0.;

  for (int i = 0; i < 3; i++) {
    for (int ipt = 0; ipt < 3; ipt++) {
      bary[i] += pts[3*ipt+i];
    }
    bary[i] /= 3;
  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

