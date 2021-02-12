/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#pragma once
#include "NodeAssociator.h"
#include "Zipper.h"


class PostNodeAssociator : public NodeAssociator
{
public:
  typedef NodeAssociator parent_type;

public:
  ///
  PostNodeAssociator(bool nodal);

  ///
  virtual ~PostNodeAssociator(void);

  ///
  void make_pairs(const K_FLD::FloatArray& pos, const std::vector< K_FLD::IntArray* >& components,
                  std::vector<E_Int> &mates, K_FLD::IntArray& OneSurface);

private:

  ///
  void __setZipMates(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& components,
                     std::vector<E_Int>& nmates);

private:
  Zipper* _zipper;
};


