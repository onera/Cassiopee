/*    
    Copyright 2013-2025 Onera.

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

//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef __NUGA_DEFS_H__
#define __NUGA_DEFS_H__

#define ZERO_M       1.e-15
#define EPSILON      1.e-12

#define Vector_t     std::vector
#define SQRT         std::sqrt

#define SIGN(a) ((a < -ZERO_M) ? -1 : ((a > ZERO_M) ? 1 : 0))  

// types
#ifndef NUGALIB

#include "Def/DefTypes.h"
#include "Def/DefCplusPlusConst.h"

#ifdef E_DOUBLEINT
  #define IDX_NONE E_Int(9223372036854775807)
#else
  #define IDX_NONE 2147483647
#endif

namespace NUGA
{
static const E_Float FLOAT_MAX = K_CONST::E_MAX_FLOAT;
static const E_Float SQRT3 = K_CONST::E_SQRT3;

static const E_Float PI   = K_CONST::E_PI;
static const E_Float PI_2 = K_CONST::E_PI_2;
static const E_Float PI_4 = K_CONST::E_PI_4;
}

#else

#define IDX_NONE     2147483647
using E_Float = double;
using E_Int   = int;
using E_Bool  = int;

#include <limits>
namespace NUGA
{
  static const E_Float FLOAT_MAX = std::numeric_limits<E_Float>::max();
  static const E_Float SQRT3     = 1.73205080756887729353;

  static const E_Float PI        = 3.14159265358979323846;
  static const E_Float PI_2      = 1.57079632679489661923;
  static const E_Float PI_4      = 0.78539816339744830962;
}

#endif

#endif
