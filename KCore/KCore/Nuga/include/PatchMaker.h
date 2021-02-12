/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#pragma once
#include "Nuga/include/DynArray.h"
#include "Nuga/include/defs.h"
#include "Nuga/include/DefContainers.h"
#include "Nuga/include/KdTree.h"
#include <vector>

class PatchMaker
{
public:
  ///
  static void run (const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectS, 
                   const NUGA::int_vector_type & pairs, E_Float angle_tolerance,
                   std::vector<K_FLD::IntArray> & connectBout);

private:

  ///
  static E_Int __update_normals(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectB,
                                const std::vector<E_Int> &pairs, K_FLD::FloatArray& normals);

  ///
  static void
    __flag_critical_nodes(const K_FLD::FloatArray& normals, E_Float critical_angle,
                          const std::vector<E_Int> & pairs, const std::vector<E_Int> & colors,
                          std::vector< std::vector<E_Int> > & sorted_nodes,
                          std::vector<E_Int> & critical_nodes);

  ///
  static E_Int
    __flag_critical_nodes_on_contour(const std::vector<E_Int>& nodes, E_Int N0, const K_FLD::FloatArray& normals,
                                     E_Float critical_angle, const std::vector<E_Int> & pairs,
                                     const std::vector<E_Int>& colors, std::vector<E_Int> & critical_nodes);

  static E_Float __getAngle(const E_Float* n1, const E_Float* n2);

  static void __build_cutting_edges(const std::vector<E_Int>& pairs,
                                    const std::vector<E_Int>& colors,
                                    const std::vector<E_Int>& critical_nodes,
                                    K_FLD::IntArray& cuttingEdges,
                                    K_FLD::IntArray& periodicEdges);

  PatchMaker(void){}
  ~PatchMaker(void){}
};
