/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "Def/DefTypes.h"
#include <math.h>
#include <limits.h>
#include <float.h>
#include "Def/DefCplusPlusConst.h"

// ==========================================================================

extern const E_Int K_CONST::E_EXIT_CODE = 13;

extern const E_Int K_CONST::E_ZERO_INT    = 0;
extern const E_Int K_CONST::E_BADVALUE_I  = -999;
extern const E_Int K_CONST::E_BADVALUE_B  = -999;

extern const E_Int K_CONST::X_AXIS_E      = 1;
extern const E_Int K_CONST::Y_AXIS_E      = 2;
extern const E_Int K_CONST::Z_AXIS_E      = 3;
extern const E_Int K_CONST::E_MAXITER_NWT = 10;


#ifdef E_DOUBLEREAL
extern const E_Float K_CONST::E_MAX_FLOAT    = DBL_MAX;
extern const E_Float K_CONST::E_MIN_FLOAT    = DBL_MIN;
extern const E_Float K_CONST::E_BADVALUE_F   = -999.999;
extern const E_Float K_CONST::E_ZERO_FLOAT   =  0.0;
extern const E_Float K_CONST::ONE_EIGHTH     =  0.125;
extern const E_Float K_CONST::ONE_FOURTH     =  0.25;
extern const E_Float K_CONST::ONE_THIRD      =  0.33333333333333333333;
extern const E_Float K_CONST::ONE_HALF       =  0.5;
extern const E_Float K_CONST::ONE            =  1.0;
extern const E_Float K_CONST::TWO            =  2.0;
extern const E_Float K_CONST::THREE          =  3.0;
extern const E_Float K_CONST::FOUR           =  4.0;
extern const E_Float K_CONST::FIVE           =  5.0;
extern const E_Float K_CONST::E_CUTOFF       = 1.e-15;
extern const E_Float K_CONST::E_ZERO_MACHINE = 1.e-15;
extern const E_Float K_CONST::E_GEOM_CUTOFF  = 1.e-8;
extern const E_Float K_CONST::E_MIN_SURFACE  = 1.e-30;
extern const E_Float K_CONST::E_MIN_VOL      = 1.e-30;
extern const E_Float K_CONST::E_PI_DEG       = 180.;
extern const E_Float K_CONST::E_TOLNWT       = 1.e-5;
extern const E_Float K_CONST::E_INFINITE     = 1.e20;

extern const E_Float  K_CONST::E_E           =  2.7182818284590452354 ;       /* e */
extern const E_Float  K_CONST::E_LOG2E       =  1.4426950408889634074 ;       /* log_2 e */
extern const E_Float  K_CONST::E_LOG10E      =  0.43429448190325182765;       /* log_10 e */
extern const E_Float  K_CONST::E_LN2         =  0.69314718055994530942;       /* log_e 2 */
extern const E_Float  K_CONST::E_LN10        =  2.30258509299404568402;       /* log_e 10 */
extern const E_Float  K_CONST::E_PI          =  3.14159265358979323846;       /* pi */
extern const E_Float  K_CONST::E_PI_2        =  1.57079632679489661923;       /* pi/2 */
extern const E_Float  K_CONST::E_PI_4        =  0.78539816339744830962;       /* pi/4 */
extern const E_Float  K_CONST::E_1_PI        =  0.31830988618379067154;       /* 1/pi */
extern const E_Float  K_CONST::E_2_PI        =  0.63661977236758134308;       /* 2/pi */
extern const E_Float  K_CONST::E_2_SQRTPI    =  1.12837916709551257390;       /* 2/sqrt(pi) */
extern const E_Float  K_CONST::E_SQRT2       =  1.41421356237309504880;       /* sqrt(2) */
extern const E_Float  K_CONST::E_SQRT1_2     =  0.70710678118654752440;       /* 1/sqrt(2) */
extern const E_Float  K_CONST::E_SQRT3       =  1.73205080756887729353;       /* sqrt(3) */

#else
extern const E_Float K_CONST::E_MAX_FLOAT    = FLT_MAX;
extern const E_Float K_CONST::E_MIN_FLOAT    = FLT_MIN;
extern const E_Float K_CONST::E_BADVALUE_F   = -999.999f;
extern const E_Float K_CONST::E_ZERO_FLOAT   =  0.0f;
extern const E_Float K_CONST::ONE_EIGHTH     =  0.125f;
extern const E_Float K_CONST::ONE_FOURTH     =  0.25f;
extern const E_Float K_CONST::ONE_THIRD      =  0.33333333333333333333f;
extern const E_Float K_CONST::ONE_HALF       =  0.5f;
extern const E_Float K_CONST::ONE            =  1.0f;
extern const E_Float K_CONST::TWO            =  2.0f;
extern const E_Float K_CONST::THREE          =  3.0f;
extern const E_Float K_CONST::FOUR           =  4.0f;
extern const E_Float K_CONST::FIVE           =  5.0f;
extern const E_Float K_CONST::E_CUTOFF       = 1.e-11f;
extern const E_Float K_CONST::E_ZERO_MACHINE = 1.e-15f;
extern const E_Float K_CONST::E_GEOM_CUTOFF  = 1.e-6f;
extern const E_Float K_CONST::E_MIN_SURFACE  = 1.e-30f;
extern const E_Float K_CONST::E_MIN_VOL      = 1.e-30f;
extern const E_Float K_CONST::E_PI_DEG       = 180.f;
extern const E_Float K_CONST::E_TOLNWT       = 1.e-5f;
extern const E_Float K_CONST::E_INFINITE     = 1.e8f;
 
extern const E_Float  K_CONST::E_E           =  2.7182818284590452354f ;      /* e */
extern const E_Float  K_CONST::E_LOG2E       =  1.4426950408889634074f ;      /* log_2 e */
extern const E_Float  K_CONST::E_LOG10E      =  0.43429448190325182765f;      /* log_10 e */
extern const E_Float  K_CONST::E_LN2         =  0.69314718055994530942f;      /* log_e 2 */
extern const E_Float  K_CONST::E_LN10        =  2.30258509299404568402f;      /* log_e 10 */
extern const E_Float  K_CONST::E_PI          =  3.14159265358979323846f;      /* pi */
extern const E_Float  K_CONST::E_PI_2        =  1.57079632679489661923f;      /* pi/2 */
extern const E_Float  K_CONST::E_PI_4        =  0.78539816339744830962f;      /* pi/4 */
extern const E_Float  K_CONST::E_1_PI        =  0.31830988618379067154f;      /* 1/pi */
extern const E_Float  K_CONST::E_2_PI        =  0.63661977236758134308f;      /* 2/pi */
extern const E_Float  K_CONST::E_2_SQRTPI    =  1.12837916709551257390f;      /* 2/sqrt(pi) */
extern const E_Float  K_CONST::E_SQRT2       =  1.41421356237309504880f;      /* sqrt(2) */
extern const E_Float  K_CONST::E_SQRT1_2     =  0.70710678118654752440f;      /* 1/sqrt(2) */
extern const E_Float  K_CONST::E_SQRT3       =  1.73205080756887729353f;      /* sqrt(3) */

#endif

// ===== KCore/Def/DefCplusPlusConst.cpp === Last line ===
