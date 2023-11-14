#ifndef _KCORE_MATH_H
#define _KCORE_MATH_H

#include "Def/DefTypes.h"

namespace K_MATH
{
  extern const E_Float ONE_THIRD;

  extern const E_Float PI;

  extern const E_Float SMALL;

  E_Float dot(const E_Float *, const E_Float *, const E_Int);
  
  void cross(const E_Float [3], const E_Float [3], E_Float [3]);
  
  E_Float norm(const E_Float *, const E_Int);
  
  E_Int feq(const E_Float, const E_Float, const E_Float tol = SMALL);
  
  void sqrmat_dot_vec(const E_Float *, const E_Float *, const E_Int, E_Float *);
  
  E_Int BiCGStab(const E_Float *, const E_Float *, const E_Int, E_Float *,
    const E_Float tol = 1e-6);
  
  void sym3mat_dot_vec(const E_Float *, const E_Float *, E_Float *);
  
  E_Float sign(const E_Float, const E_Float tol = 1e-15);
  
  void sym3mat_dot_sym3mat(const E_Float *, const E_Float *, E_Float *);
  
  E_Float sym3mat_det(const E_Float [6]);
  
  E_Float sym3mat_trace(const E_Float [6]);
  
  E_Float sym3mat_second_invariant(const E_Float *);
  
  E_Float sym3mat_third_invariant(const E_Float *);
  
  void sym3mat_eigen(const E_Float [6], E_Float L[3], E_Float v1[3],
    E_Float v2[3], E_Float v3[3], const E_Float tol = SMALL);
}

#endif