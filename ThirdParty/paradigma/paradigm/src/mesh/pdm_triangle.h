#ifndef __PDM_TRIANGLE_H__
#define __PDM_TRIANGLE_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"

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
 * \enum PDM_triangle_status_t
 * \brief Triangle status type
 *
 */

typedef enum {

  PDM_TRIANGLE_INSIDE      = 0,  /*!< Inside  */
  PDM_TRIANGLE_OUTSIDE     = 1,  /*!< Outside */               
  PDM_TRIANGLE_DEGENERATED = 2,  /*!< Degenerated */               

} PDM_triangle_status_t;


/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes 
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
 );


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
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_SURF_MESH_H__ */
