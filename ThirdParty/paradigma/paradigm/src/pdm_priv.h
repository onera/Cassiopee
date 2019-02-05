#ifndef __PDM_PRIV_H__
#define __PDM_PRIV_H__

#include <stdio.h>
#include <math.h>

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
 * Macro definitions
 *============================================================================*/

enum {X, Y, Z};

/** 
 * Absolute value
 */

#define PDM_ABS(a)     ((a) <  0  ? -(a) : (a))

/** 
 * Minimum value
 */

#define PDM_MIN(a,b)   ((a) < (b) ?  (a) : (b))

/** 
 * Maximum value
 */

#define PDM_MAX(a,b)   ((a) > (b) ?  (a) : (b))

/** 
 * Dot product
 */

#define PDM_DOT_PRODUCT(vect1, vect2)                                   \
  (vect1[X] * vect2[X] + vect1[Y] * vect2[Y] + vect1[Z] * vect2[Z])

/** 
 * Module
 */

#define PDM_MODULE(vect)                                                \
  sqrt(vect[X] * vect[X] + vect[Y] * vect[Y] + vect[Z] * vect[Z])

#define PDM_CROSS_PRODUCT(prod_vect, vect1, vect2)  \
  (prod_vect[X] = vect1[Y] * vect2[Z] - vect2[Y] * vect1[Z], \
   prod_vect[Y] = vect2[X] * vect1[Z] - vect1[X] * vect2[Z], \
   prod_vect[Z] = vect1[X] * vect2[Y] - vect2[X] * vect1[Y])

#define PDM_DETERMINANT2X2(vect1, vect2) \
  (vect1[X] * vect2[Y] - vect2[X] * vect1[Y] )

#define PDM_DOT_PRODUCT_2D(vect1, vect2) \
  (vect1[X] * vect2[X] + vect1[Y] * vect2[Y])

#define PDM_PI 3.1415926535897931


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PRIV_H__ */
