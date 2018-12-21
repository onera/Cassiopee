/*    
    Copyright 2013-2019 Onera.

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
#ifndef _KCORE_DEF_DEFCPLUSPLUSCONST_H_
#define _KCORE_DEF_DEFCPLUSPLUSCONST_H_

///+ C++ Constants
namespace K_CONST
{
/** exit code (13) */
  extern const E_Int E_EXIT_CODE;   // 13

/** 0 (integer) */
  extern const E_Int E_ZERO_INT;    // 0

/** */
  extern const E_Int E_BADVALUE_I;  // -999

/** */
  extern const E_Int E_BADVALUE_B;  // -999

/** 1 */
  extern const E_Int X_AXIS_E;

/** 2 */
  extern const E_Int Y_AXIS_E;

/** 3 */
  extern const E_Int Z_AXIS_E;

/** Maximum Default Newton Iteration Number (10) */
  extern const E_Int E_MAXITER_NWT;


// Real (float or double)
/** DBL_MAX */
  extern const E_Float E_MAX_FLOAT;

/** DBL_MIN */
  extern const E_Float E_MIN_FLOAT;

/** */
  extern const E_Float E_ZERO_FLOAT;

/** */
  extern const E_Float ONE_EIGHT;

/** */
  extern const E_Float ONE_FOURTH;

/** */
  extern const E_Float ONE_HALF;

/** */
  extern const E_Float ONE;

/** */
  extern const E_Float TWO;

/** */
  extern const E_Float THREE;

/** */
  extern const E_Float FOUR;

/** */
  extern const E_Float FIVE;

/** */
  extern const E_Float E_CUTOFF;               // 1.E-11

/** */
  extern const E_Float E_ZERO_MACHINE;         // 1.E-15

/** cut off for geometric operations    */ 
  extern const E_Float E_GEOM_CUTOFF;          // 1.E-8

/** */
  extern const E_Float E_MIN_SURFACE;

/** big (infinite) positive value       */
  extern const E_Float E_INFINITE;             //  1.E20

/** PI (in degrees)                     */
  extern const E_Float E_PI_DEG;

/** Newton convergence tolerance (1E-5) */
  extern const E_Float E_TOLNWT;

/** Should be NaN                       */
  extern const E_Float E_BADVALUE_F;

/** e                                   */
  extern const E_Float  E_E;
/**   log_2 e                           */
  extern const E_Float  E_LOG2E;
/** log_10 e                            */
  extern const E_Float  E_LOG10E;
/** log_e 2                             */
  extern const E_Float  E_LN2;
/** log_e 10                            */
  extern const E_Float  E_LN10;
/**  pi                                 */
  extern const E_Float  E_PI;
/** pi/2                                */
  extern const E_Float  E_PI_2;
/** pi/4                                */
  extern const E_Float  E_PI_4;
/** 1/pi                                */
  extern const E_Float  E_1_PI;
/** 2/pi                                */
  extern const E_Float  E_2_PI;
/** 2/sqrt(pi)                          */
  extern const E_Float  E_2_SQRTPI;
/** sqrt(2)                             */
  extern const E_Float  E_SQRT2;
/** 1/sqrt(2)                           */
  extern const E_Float  E_SQRT1_2;
/** sqrt(3)                             */
  extern const E_Float  E_SQRT3;
}
///-

#define E_MAXSTRINGSIZE 80

#endif

// ===== KCore/Def/DefCplusPlusConst.h === Last line ===
