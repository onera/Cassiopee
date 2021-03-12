/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#pragma once
#include "Nuga/include/DynArray.h"
#include "Intersector.h"
#include <vector>

class NodeAssociator
{
public:
  ///
  NodeAssociator(void);

  ///
  virtual ~NodeAssociator(void);

  ///
  virtual void make_pairs(const K_FLD::FloatArray& pos, const std::vector< K_FLD::IntArray* >& components,
                          std::vector<E_Int> &mates, K_FLD::IntArray& OneSurface);

protected:
  ///
  void __removeOverlaps(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& components,
                        std::vector<E_Int>& nmates, std::vector<K_FLD::IntArray>& componentsOut);

  ///
  void __make_pairs(const K_FLD::FloatArray& pos, const std::vector< std::vector<E_Int> >& sorted_nodes,
                            const std::vector<E_Int>& surface_colors, std::vector<E_Int> &pairs);

  ///
  static void __reorientComponent(const K_FLD::FloatArray& pos, K_FLD::IntArray& component, const std::vector<E_Int>& nmates); 

private:
  ///
  static void __priorize_X_zones(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& components, const std::vector<XPair>& pairs, 
                                 std::vector<std::vector<E_Int> >& componentsKeep, std::vector<std::vector<E_Int> >& componentsRemove,
                                 std::vector<std::vector<E_Int> >& componentsUnchanged);

  ///
  static void __priorize_X_zones2(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& components, const std::vector<XPair>& pairs, 
                                  std::vector<std::vector<E_Int> >& componentsKeep, std::vector<std::vector<E_Int> >& componentsRemove);

  ///
  static void __getRelativeOrient(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray>& surfaces, const std::vector<E_Int>& nmates,
                                  std::map<K_MESH::NO_Edge, E_Int>& surfpair_to_orient); 

  ///
  static void __getRelativeOrient(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray>& comps, const std::vector<XPair>& pairs,
                                  std::map<K_MESH::NO_Edge, E_Int>& comppair_to_orient);

  ///
  static void __reorientComponents(const K_FLD::FloatArray& pos, std::vector<K_FLD::IntArray>& component, const std::vector<XPair>& pairs);

public:
  static bool reorient;

};


