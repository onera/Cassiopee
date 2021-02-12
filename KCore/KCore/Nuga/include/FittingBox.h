/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

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

