/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef __LINEAR_DELAUNAY_MATH_H__
#define __LINEAR_DELAUNAY_MATH_H__

#include "Nuga/include/defs.h"
#include "Nuga/include/DynArray.h"

namespace K_LINEAR
{
  class DelaunayMath
  {
  public:
    static void eigen_vectors 
      (E_Float a00, E_Float a11, E_Float a10, E_Float& lambda0, E_Float& lambda1, E_Float* v0, E_Float* v1) ;

    static void eigen_values 
      (E_Float a00, E_Float a11, E_Float a10, E_Float& lambda0, E_Float& lambda1) ;

    static void simultaneous_reduction
      (const K_FLD::FloatArray& M1, const K_FLD::FloatArray& M2,
       E_Float* V1, E_Float* V2);

    static void intersect(const K_FLD::FloatArray& Mi, const K_FLD::FloatArray& M2, K_FLD::FloatArray& I);

    static void resoLin(E_Float a11, E_Float a12, E_Float a21, E_Float a22,
                        E_Float d1, E_Float d2, E_Float&x, E_Float&y);
  };
}

#endif
