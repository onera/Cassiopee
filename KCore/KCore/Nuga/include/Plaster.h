/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef __GENERATOR_PLASTER_H__
#define __GENERATOR_PLASTER_H__

#include "Nuga/include/DynArray.h"
#include "Nuga/include/defs.h"
#include "Nuga/include/DefContainers.h"

class Plaster
{
public:
  Plaster(void);
  ~Plaster(void);

  E_Int make(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectE2,
             K_FLD::FloatArray& posPlaster, E_Int& ni, E_Float bump_factor = 0.);

private:

  void __cartesian (const E_Float* minB, const E_Float* maxB, E_Int ni, E_Int nj, K_FLD::FloatArray& cart);

  void __smooth(std::vector<E_Float>& z, E_Int ni, E_Float bump_factor, E_Float tol = EPSILON);

  void __smooth_1(std::vector<E_Float>& z, E_Int ni, E_Float tol = EPSILON);

  void __smooth_2(std::vector<E_Float>& z, E_Int ni, E_Float tol = EPSILON);

#ifdef WIN32
#ifdef E_DEBUG
  void make_box(const E_Float* minB, const E_Float* maxB, K_FLD::FloatArray& pos, K_FLD::IntArray& connect);
  void drawBox(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const E_Float* mB, const E_Float* MB);
#endif
#endif

  void __initializePlaster(const K_FLD::FloatArray& plaster2D, E_Int ni, const K_FLD::FloatArray& pos2D,
                           const K_FLD::IntArray& connectE2/*, const std::vector<E_Int>& hard_nodes*/, const std::vector<E_Float>& zE2,
                           const K_FLD::IntArray& connectT3, std::vector<E_Float>& z, E_Float bump_factor = 0.);

  void __blockSurroundingNodes(const K_FLD::FloatArray& plaster2D, E_Int ni, 
                               const NUGA::bool_vector_type& mask,
                               const K_FLD::FloatArray& pos2D,
                               const K_FLD::IntArray& connectE2, const std::vector<E_Float>& zE2,
                               std::vector<E_Float>& z, NUGA::int_set_type& onodes);
/*
  void __blockInsideNodes(const K_FLD::FloatArray& plaster2D, E_Int ni, 
                          const K_FLD::FloatArray& pos2D, const std::vector<E_Int>& hard_nodes,
                          const std::vector<E_Float>& zE2, std::vector<E_Float>& z);
*/
  ///
  void __bumpPlaster(const K_FLD::FloatArray& plaster2D, E_Int ni,
                     const NUGA::bool_vector_type& mask,
                     E_Float bump_factor, const NUGA::int_set_type& onodes,
                     std::vector<E_Float>& z);

  void __mask(const K_FLD::FloatArray& pos2D, const K_FLD::FloatArray& plaster2D,
              const K_FLD::IntArray& connectT3, NUGA::bool_vector_type& mask);

  bool __IsStrictlyInT3(const E_Float* P, const E_Float* P0, const E_Float* P1, const E_Float* P2);

  void __getPlasterBoundary(const NUGA::bool_vector_type& mask, E_Int ni, bool outside,
                            NUGA::int_set_type& nodes);

  void __blockNodes(const K_FLD::FloatArray& pos2D, const K_FLD::FloatArray& plaster2D,
                    const K_FLD::IntArray& connectE2, const std::vector<E_Float>& zE2,
                    const NUGA::int_set_type& nodes, std::vector<E_Float>& z);

  E_Float __computeCharacteristicLength(const K_FLD::FloatArray& pos2D, const K_FLD::IntArray& connectE2);

#ifdef WIN32
#ifdef E_DEBUG
  static int _count;
#endif
#endif

};

#endif
