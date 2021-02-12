/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#pragma once
#include "Zipper.h"

class MergingZipper : public Zipper
{

public:

  ///
  MergingZipper(bool reorient);

  ///
  virtual ~MergingZipper(void);

  ///
  virtual void setMates(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
                        const K_FLD::IntArray& connectE2, std::vector<E_Int>& nmates);

  /// Merges a nodal mesh. Does nothing otherwise.
  virtual void merge(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
                     K_FLD::IntArray& connectE2, std::vector<E_Int>& nmates);

private:
  E_Float __computeTolerance(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
                             const std::vector<E_Int>& nodes);

  E_Int __removeDuplicates(K_FLD::IntArray& connect);

};
