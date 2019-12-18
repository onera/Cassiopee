#ifndef __PDM_PLANE_H__
#define __PDM_PLANE_H__

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

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 * \brief Computes normal
 *
 * \param [in]   numPts  Number of polygon vertices
 * \param [in]   pts     Polygon vertices coordinates
 * \param [out]  n       Normal
 *
 */

void
PDM_plane_normal
(
 const int     numPts,
 const double *pts,
 double        n[3]
);


/**
 * \brief Performs plane projection
 *
 * \param [in]   x       Point to project
 * \param [in]   origin  Plane origin
 * \param [in]   n       Plane normal
 * \param [out]  cp      Projected point
 *
 */

void
PDM_plane_projection
(
const double x[3],
const double origin[3],
const double n[3],
      double cp[3]
);


/**
 * \brief Performs plane projection
 *
 * \param [in]   x       Point to project
 * \param [in]   pt      Point inside the plane
 * \param [in]   n       Plane normal
 * \param [out]  cp      Projected point
 *
 */

void
PDM_plane_projection2
(
const double x[3],
const double pt[3],
const double n[3],
      double cp[3]
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PLANE_H__ */
