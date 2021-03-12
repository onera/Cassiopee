/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#pragma once
#include "Nuga/include/DynArray.h"
#include "Nuga/include/BbTree.h"
#include <vector>

struct XPair{

  XPair(E_Int s0, E_Int c0, E_Int s1, E_Int c1):S0(s0), C0(c0), S1(s1), C1(c1){}

  E_Int S0;
  E_Int C0;
  E_Int S1;
  E_Int C1;
};

class Intersector
{
  typedef K_SEARCH::BBox3D BBoxType;
  typedef K_SEARCH::BbTree3D TreeType;

public:
  static void getXPairs(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray* >& Cs,
                 std::vector<XPair>& pairs);
private:

  ///
  static void __getXPairs(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& Cs,
                          E_Int s1, E_Int s2,
                          const std::vector<TreeType* >& trees, 
                          const std::vector< std::vector<BBoxType*> >& boxes,
                          const std::vector<std::vector<E_Float> >& Lengths,
                          E_Float tol, bool tol_is_absolute,
                          std::vector<XPair>& pairs);
  
  ///
  static void __create_boxes (const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt, 
                              std::vector<K_SEARCH::BBox3D*>& boxes, K_SEARCH::BBox3D*& pool, K_SEARCH::BBox3D& GBbox, E_Float enlarge_factor=0.);

  ///
  inline static bool __overlap(const K_FLD::FloatArray& pos, K_FLD::IntArray::const_iterator pS1,
                               K_FLD::IntArray::const_iterator pS2,
                               E_Float tolerance, E_Bool tol_is_absolute, E_Float Lc);

  static bool __pointIsInside(const K_FLD::FloatArray& pos, K_FLD::IntArray::const_iterator pS1,
                              K_FLD::IntArray::const_iterator pS2, E_Float tolerance, E_Bool tol_is_absolute,
                              E_Float Lc);

  static bool __edgeOverlaps(const K_FLD::FloatArray& pos, K_FLD::IntArray::const_iterator pS1,
                             K_FLD::IntArray::const_iterator pS2, E_Float tolerance, E_Bool tol_is_absolute,
                             E_Float Lc);


  Intersector(void){}
  ~Intersector(void){}
};
