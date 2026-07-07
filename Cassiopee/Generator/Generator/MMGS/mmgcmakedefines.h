/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

#ifndef MMGCMAKEDEFINE_H
#define MMGCMAKEDEFINE_H

#include "Def/DefTypes.h"

/* inttypes.h is needed to handle prints of MMG5_int using PRId32 and PRId64 macros */
#include <inttypes.h>

#ifndef _WIN32
#define MMG_POSIX
#define MMG_GNU
#endif

//@DEF_POSIX@
//@DEF_GNU@

//@DEF_MMG5_INT@ /*!< Integer type for C */
#ifdef E_DOUBLEINT
#define MMG5_int int64_t
#else
#define MMG5_int int32_t
#endif

//@DEF_MMG5_INTMAX@ /*!< INT_MAX or LONG_MAX depending on MMG5_INT size */
#define MMG5_INTMAX E_MAXINT


//@DEF_MMG5_PRId@ /*!< Printing format for MMG5_int type */
#if defined E_DOUBLEINT
#if defined _WIN32
#define MMG5_PRId "lld"
#define MMG5_abs llabs
#else
#define MMG5_PRId "ld"
#define MMG5_abs labs
#endif
#else
#define MMG5_PRId "d"
#define MMG5_abs abs
#endif

//@DEF_MMG_SWPBIN@ /*!< MMG5_swapbin function for MMG5_int */
#define MMG_SWPBIN 

//@DEF_MMG_ABS@ /*!< Abs function for MMG5_int */

//#cmakedefine   USE_POINTMAP /*!< Flag to enable and export the pointmap used */
//#cmakedefine01 MMG_DYN_LIB

#endif
