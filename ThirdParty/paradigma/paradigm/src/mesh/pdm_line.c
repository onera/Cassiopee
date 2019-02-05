/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_line.h"

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

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 * \brief Compute 2x2 determinant
 *
 * \param [in]  a 
 * \param [in]  b
 * \param [in]  c
 * \param [in]  d
 *
 * \return      \ref PDM_TRUE or \ref PDM_FALSE
 *
 */

static double
_det_2x2
(
 double a,
 double b,
 double c,
 double d
)
{
  return (a * d - b * c);
}


/**
 * \brief Solve 2x2 system
 *
 * \param [in]     A  Matrix
 * \param [inout]  x  right hand side in, solution out
 *
 * \return   \ref PDM_FALSE if matrix is singular, \ref PDM_TRUE otherwise
 *
 */

static int
_solve_2x2
(
 double A[2][2],
 double x[2]
)
{
  double y[2];
  
  double det = _det_2x2 (A[0][0], A[0][1], A[1][0], A[1][1]);

  if (det == 0.0) {
    return PDM_FALSE;
  }
  
  y[0] = (A[1][1]*x[0] - A[0][1]*x[1]) / det;
  y[1] = (-A[1][0]*x[0] + A[0][0]*x[1]) / det;
  
  x[0] = y[0];
  x[1] = y[1];
  return PDM_TRUE;
}

/*=============================================================================
 * Public function prototypes 
 *============================================================================*/


/**
 * \brief Performs intersection of two finite 3D lines
 *
 *  An intersection is found if the projection of the two lines onto the plane 
 *  perpendicular to the cross product of the two lines intersect. 
 *  The parameters (u,v) are the parametric coordinates of the lines at the 
 *  position of closest approach 
 *
 * \param [in]  a1 Coordinates of the first line vertex of 'a'
 * \param [in]  a2 Coordinates of the second line vertex of 'a'
 * \param [in]  b1 Coordinates of the first line vertex of 'b'
 * \param [in]  b2 Coordinates of the second line vertex of 'b'
 * \param [out] u  Parameter of the intersection in line 'a' parametric coordinates 
 * \param [out] v  Parameter of the intersection in line 'b' parametric coordinates 
 *
 * \return      \ref PDM_TRUE or \ref PDM_FALSE
 *
 */

PDM_line_intersect_t
PDM_line_intersection
(
 const double a1[3],
 const double a2[3],
 const double b1[3],
 const double b2[3],
 double *u,
 double *v
 )
{
 
  double a21[3], b21[3], b1a1[3];
  double c[2];
  double A[2][2];

  *u = 0.;
  *v = 0.;

  /*
   * Determine vectors
   */

  a21[0] = a2[0] - a1[0];
  a21[1] = a2[1] - a1[1];
  a21[2] = a2[2] - a1[2];

  b21[0] = b2[0] - b1[0];
  b21[1] = b2[1] - b1[1];
  b21[2] = b2[2] - b1[2];
  
  b1a1[0] = b1[0] - a1[0];
  b1a1[1] = b1[1] - a1[1];
  b1a1[2] = b1[2] - a1[2];

  /*
   * Define least squares system matrix.
   */

  A[0][0] = PDM_DOT_PRODUCT ( a21, a21 );
  A[0][1] = - PDM_DOT_PRODUCT ( a21, b21 );
  A[1][0] = A[0][1];
  A[1][1] = PDM_DOT_PRODUCT( b21, b21 );

  /*
   * Compute the least squares system constant term.
   */

  c[0] = PDM_DOT_PRODUCT( a21, b1a1 );
  c[1] = - PDM_DOT_PRODUCT( b21, b1a1 );

  /*
   * Solve the system of equations
   */

  if ( _solve_2x2 (A, c) == 0 ) {
    return PDM_LINE_INTERSECT_ON_LINE;
  }
  else {
    *u = c[0];
    *v = c[1];
  }
  
  /*
   * Check parametric coordinates for intersection.
   */

  if ( (0.0 <= *u) && (*u <= 1.0) && (0.0 <= *v) && (*v <= 1.0) ) {
    return PDM_LINE_INTERSECT_YES;
  }
  else {
    return PDM_LINE_INTERSECT_NO;
  }
}


/**
 * \brief Computes point-line distance 
 *
 * \param [in]  x             Point coordinates
 * \param [in]  p1            First line vertex coordinates
 * \param [in]  p2            Second line vertex coordinates
 * \param [out] t             Parameter of the intersection in line parametric coordinates 
 * \param [out] closest_point Closest point
 *
 * \return   The square of the distance
 *
 */

double
PDM_line_distance
(
 const double x[3],
 const double p1[3],
 const double p2[3],
 double *t,
 double closestPoint[3]
 )
{

  const double _tol_dist = 1e-5;

  double p21[3], denom, num;
  double *closest;
  
  /*
   * Determine appropriate vectors
   */
  
  p21[0] = p2[0]- p1[0];
  p21[1] = p2[1]- p1[1];
  p21[2] = p2[2]- p1[2];

  /*
   *  Get parametric location
   */
  
  num = p21[0]*(x[0]-p1[0]) + p21[1]*(x[1]-p1[1]) + p21[2]*(x[2]-p1[2]);
  denom = PDM_DOT_PRODUCT(p21,p21);

  double tolerance = fabs (_tol_dist * num);
  if ( fabs(denom) < tolerance ) {
    closest = (double *) p1; 
  }

  /* 
   *  If parametric coordinate is within 0<=p<=1, then the point is closest to
   *  the line.  Otherwise, it's closest to a point at the end of the line.
   */
  
  else if ( denom <= 0.0 || (*t=num/denom) < 0.0 ) {
    closest = (double *) p1;
  }
  
  else if ( *t > 1.0 ) {
    closest = (double *) p2;
  }
  
  else {
    closest = p21;
    p21[0] = p1[0] + (*t)*p21[0];
    p21[1] = p1[1] + (*t)*p21[1];
    p21[2] = p1[2] + (*t)*p21[2];
  }

  closestPoint[0] = closest[0];
  closestPoint[1] = closest[1];
  closestPoint[2] = closest[2];

  double v[3] = {closest[0] - x[0],
                 closest[1] - x[1],
                 closest[2] - x[2]};

  return  PDM_DOT_PRODUCT(v,v);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */

