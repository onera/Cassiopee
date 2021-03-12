/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef __GENERATOR_GAP_FIXER_H__
#define __GENERATOR_GAP_FIXER_H__

#include "Nuga/include/defs.h"
#include "Nuga/include/DefContainers.h"

class GapFixer
{
public:
  typedef NUGA::int_vector_type int_vector_type;
public:
  ///
  GapFixer();

  ///
  ~GapFixer(void);

  /// Main function to fix a gap.
  static E_Int run(const K_FLD::FloatArray& posC,
                   E_Int nj,
                   const K_FLD::FloatArray& posB0,
                   const K_FLD::IntArray& connectB0,
                   K_FLD::FloatArray& posG,
                   K_FLD::IntArray& connectG,
                   E_Bool refine,
                   const K_FLD::FloatArray* coordHP = 0); // hard points

private:
  /// Converts the parametrization to the quads global one.
  static void __convertToGlobalQInterp(const K_FLD::FloatArray& posB0, E_Int nj,
                                       const NUGA::int_vector_type& cell_indices,
                                       K_FLD::FloatArray& posUV);

};

#endif
