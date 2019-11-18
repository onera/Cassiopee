/*
 
 
 
              NUGA 
 
 
 
 */

#ifndef __NUGA_MACROS_H__
#define __NUGA_MACROS_H__

#include<memory>

#define ALL(c)  (c).begin(), (c).end()

#define IS_IN(c,v) (c).find(v) != (c).end()

#define zSIGN(a, TOL) ( (a < -TOL)? -1 : (a > TOL) ? 1 : 0 )

#define ZERO_M 1.e-15
#define EPSILON 1.e-12

#define Vector_t std::vector

#define STACK_ARRAY(T, n, name) std::unique_ptr<T[]> name(new T[n]);

#endif
