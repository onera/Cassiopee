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

#include <stdlib.h>
#include <limits.h>
#include <stdint.h>

// ==========================================================================
#ifndef _KCORE_DEF_DEFTYPES_H_
#define _KCORE_DEF_DEFTYPES_H_

#define E_DOUBLEREAL
//#define E_DOUBLEINT

// ==========================================================================
///+ Basic types

// Essai de trouve si long long existe (8 octets)
#ifdef LLONG_MAX
#define E_LONG long long
//int64_t
//long long
#else
#define E_LONG long
#endif

#ifdef E_DOUBLEREAL
  typedef double E_Float;
  #ifdef E_MPI
    #define E_PCM_FLOAT MPI_DOUBLE
  #else
    #define E_PCM_FLOAT sizeof(E_Float)
  #endif
#else
  typedef float E_Float;
  #ifdef E_MPI
    #define E_PCM_FLOAT MPI_FLOAT
  #else
    #define E_PCM_FLOAT sizeof(E_Float)
  #endif
#endif

// Int (int or long long)
#ifdef E_DOUBLEINT
  #if !defined(_ELSA_COMPILER_NEC_) && !defined(_ELSA_COMPILER_HP64_) && !defined(_ELSA_COMPILER_DEC_)
    #if defined(_ELSA_COMPILER_CRAY_)
      #define E_Int long 
      #define E_Boolean long
      #define E_Bool long
    #else
      typedef long long E_Int;
      typedef long long E_Boolean;
      typedef long long E_Bool;
    #endif
    #ifdef E_MPI
      #define E_PCM_INT MPI_LONG_LONG
    #else
      #define E_PCM_INT sizeof(E_Int)
    #endif
  #else
    typedef long E_Int;
    typedef long E_Boolean;
    typedef long E_Bool;
    #ifdef E_MPI
      #define E_PCM_INT MPI_LONG
    #else
      #define E_PCM_INT sizeof(E_Int)
    #endif
  #endif
#else
  typedef int  E_Int;
  typedef int  E_Boolean;
  typedef int  E_Bool;
  #ifdef E_MPI
    #define E_PCM_INT MPI_INT
  #else
    #define E_PCM_INT sizeof(E_Int)
  #endif
#endif

///-

/** None index in an array */
#define E_IDX_NONE        2147483647

#define E_EPSILON         1.e-12

///-

#endif
// ===== Def/DefTypes.h === Last line ===
