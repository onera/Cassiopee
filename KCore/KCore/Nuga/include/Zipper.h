/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#pragma once
#include "Nuga/include/DynArray.h"
#include "Nuga/include/defs.h"

class Zipper
{
public:
  enum eMateType {FREE = -9, OVERLAP = -1, UNKNO = IDX_NONE};

public:

  ///
  Zipper(bool reorient);

  ///
  virtual ~Zipper(void);

  ///
  virtual void setMates(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
                        const K_FLD::IntArray& connectE2, std::vector<E_Int>& nmates);

  /// Merges a nodal mesh. Does nothing otherwise.
  virtual void merge(const K_FLD::FloatArray& pos, K_FLD::IntArray& connect,
                     K_FLD::IntArray& connectE2, std::vector<E_Int>& nmates){}

  ///
  static void zip(const K_FLD::IntArray& connectE2, const std::vector<E_Int>& nmates, K_FLD::IntArray& zipped);

protected:
  ///
  static void __setNodeMates(const K_FLD::FloatArray& vpos, const K_FLD::IntArray& vE2, const std::vector<E_Int>& realID,
                             std::vector<E_Int>& nmates);

private:
  ///
  static void __setBottoms(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectT3, const K_FLD::IntArray& connectE2, K_FLD::IntArray& bottomE2);

  ///
  static void __computeVirtualEdges(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectE2,
                                     const K_FLD::IntArray& bottomE2,
                                     K_FLD::FloatArray& vpos, K_FLD::IntArray& vE2,
                                     std::vector<E_Int>& realID);

/*
  ///
  static void __computeRadius(const K_FLD::FloatArray& vpos, const K_FLD::IntArray& vE2, const K_FLD::FloatArray& pos,
                              const std::vector<E_Int>& realID, std::vector<E_Float>& R2);
*/
protected:
  bool _reorient;

};
