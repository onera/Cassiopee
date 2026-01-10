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
//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef __GENERATOR_FITTING_BOX_H__
#define __GENERATOR_FITTING_BOX_H__

#include "Nuga/include/DynArray.h"
#include "Nuga/include/maths.hxx"

class FittingBox
{
public:
  /// Computes the normal to the contour.
  static void computeNormalToContour(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, E_Float* Z);
  
  /// Computes the transformation matrix to the coordinate system (having W as 3rd axis) minimizing the bounding box.
  static void computeFittingFrame(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
                                  const E_Float* Z, K_FLD::FloatArray& F);
  /// Computes the transformation matrix to the coordinate system (having W as 3rd axis) minimizing the bounding box
  /** and optimizing the view over the contour */
  static E_Int computeOptimalViewFrame(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, K_FLD::FloatArray& F);
 

private:
  FittingBox(void){}
  ~FittingBox(void){}

  static E_Int __computeOptimalViewFrame(const K_FLD::FloatArray& posE2, const K_FLD::IntArray& connectE2,
                                         const E_Float* W0, const K_FLD::FloatArray& R, K_FLD::FloatArray& P);

  static void __fitByRotating(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
                              K_FLD::FloatArray& P);

#ifdef E_DEBUG
  static void drawBox(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const E_Float* mB, const E_Float* MB);
  static void make_box(const E_Float* minB, const E_Float* maxB, K_FLD::FloatArray& boxPs, K_FLD::IntArray& boxC);
#endif

private:
  static const E_Int _maxStep = 40;

};

#endif

