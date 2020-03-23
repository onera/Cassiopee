/*
 
 
 
              NUGA 
 
 
Author : Sam Landier (sam.landier@onera.fr) 
 */

#ifndef NUGA_VECTOR_HXX
#define NUGA_VECTOR_HXX

#include "Def/DefFunction.h"


namespace NUGA
{
  
  ///
  inline double project(double const * plane_pt, double const * plane_dir, double const * pt, double const * dir, double* proj_pt)
  {
    double PPt[3];
    K_FUNC::diff<3>(pt, plane_pt, PPt);

    double k = -K_FUNC::dot<3>(PPt, plane_dir);
    double c = K_FUNC::dot<3>(plane_dir, dir);

    //assert(SIGN(c, EPSILON) != 0); //fixme
    k /= c;

    K_FUNC::sum<3>(1., pt, k, dir, proj_pt); //project

    return k;
  }

} //NUGA

#endif