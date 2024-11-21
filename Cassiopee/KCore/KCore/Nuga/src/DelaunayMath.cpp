/*    
    Copyright 2013-2024 Onera.

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

#include "Nuga/include/DelaunayMath.h"
#include "Nuga/include/maths.hxx"

namespace K_LINEAR
{
//=============================================================================
void DelaunayMath::eigen_values 
(E_Float a00, E_Float a11, E_Float a10, E_Float& lambda0, E_Float& lambda1) 
{
  // eigen values are solution of: x^2 - tr(A)*x + det(A)
  // therefore x=(1/2)*(tr(A) +- sqrt(delta)) where delta=tr(A)^2-4*det(A)
  // Assumption: we always have 2 (not necessarily distinct) solutions.
  //             therefore delta is assumed positive or null.
   E_Float trA    = a00+a11;
   E_Float detA   = a00*a11 - a10*a10;
   E_Float delta  = trA*trA - 4.*detA;          
   if (delta > 0.) // Valid discriminant.
   {
     delta=::sqrt(delta);
     lambda0=0.5*(trA-delta);
     lambda1=lambda0+delta;
   }
   else //delta is set to 0. 
     lambda1=lambda0=0.5*trA;
}

//=============================================================================
void DelaunayMath::simultaneous_reduction
(const K_FLD::FloatArray& M1, const K_FLD::FloatArray& M2,
 E_Float* V1, E_Float* V2)
{

  K_FLD::FloatArray M1b(M1);
  K_FLD::FloatArray::inverse2(M1b);

  K_FLD::FloatArray NN(M1b * M2);

  if (::fabs(NN(1,0)) < 1.e-10)//fixme value epsilon?
  {
    V1[0] = 1.; V1[1] = 0.;
    V2[0] = 0.; V2[1] = 1.;
  }
  else
  {
    E_Float lambda1, lambda2;
    DelaunayMath::eigen_vectors(NN(0,0), NN(1,1), NN(1,0), lambda1, lambda2, V1, V2);
  }
}

//=============================================================================
void DelaunayMath::eigen_vectors 
(E_Float a00, E_Float a11, E_Float a10, 
 E_Float& lambda0, E_Float& lambda1, E_Float* v0, E_Float* v1) 
{
  eigen_values(a00, a11, a10, lambda0, lambda1);

  // Diag. matrix case
  if ((::fabs(a10) < EPSILON) && (::fabs(lambda0-a11)< EPSILON))
  {
    if (a00 <= a11)
    {
      *v0     = 1.;
      *(v0+1) = 0.;
    }
    else
    {
      *v0     = 0.;
      *(v0+1) = 1.;
    }
  }
  else //non diag matrix.
  { 
    *v0 = lambda0 - a11; 
    *(v0+1) = a10;
    NUGA::normalize<2>(v0);
  }

  *v1     = -*(v0+1);  // (vo, v1) is orthonormal
  *(v1+1) = *v0;
}

//=============================================================================
void DelaunayMath::intersect(const K_FLD::FloatArray& M1, 
                             const K_FLD::FloatArray& M2, 
                             K_FLD::FloatArray& I)
{

  E_Float lambda1, lambda2, mu1, mu2, V1[2], V2[2];
  
  eigen_vectors(M1(0,0), M1(1,1), M1(1,0), lambda1, lambda2, V1, V2);
  eigen_vectors(M2(0,0), M2(1,1), M2(1,0), mu1, mu2, V1, V2);

  simultaneous_reduction(M1, M2, V1, V2);

  E_Float l1 = (lambda1 < mu1) ? mu1 : lambda1;
  E_Float l2 = (lambda2 < mu2) ? mu2 : lambda2;

  E_Float b00 = V1[0];
  E_Float b10 = V1[1];
  E_Float b01 = V2[0];
  E_Float b11 = V2[1];

  E_Float det2i = (b00*b11) - (b10*b01);
  det2i = 1. / (det2i*det2i);

  I.resize(2,2);
  I(0,0) = det2i*((l1*b11*b11) + (l2*b01*b01));
  I(1,1) = det2i*((l1*b10*b10) + (l2*b00*b00));
  I(0,1) = I(1,0) = - det2i*(((l1*b10*b11) + (l2*b00*b01)));
}

//=============================================================================
void DelaunayMath::resoLin
(E_Float a11, E_Float a12, E_Float a21, E_Float a22,
 E_Float d1, E_Float d2, E_Float& x, E_Float& y)
{
  x = y = NUGA::FLOAT_MAX;
  E_Float detA = a11*a22 - a12*a21;

  if (detA <= 0.) return;
  detA = 1./detA;

  x = a22*d1 - a12*d2;
  y = a11*d2 - a21*d1;

  x *= detA;
  y *= detA;
}

}
