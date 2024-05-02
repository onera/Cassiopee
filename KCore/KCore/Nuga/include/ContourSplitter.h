/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef _CONTOURSPLITTER_H_
#define _CONTOURSPLITTER_H_

#include "Nuga/include/EltAlgo.h"


template <typename ElementType, typename BoundaryType>
class ContourSplitter
{
public:
  /// Splits the connectivity (only) with the edge N0N1 and stores the bits into separate containers.
  /** WARNING : it is the responsibility of the caller to destroy connectBout.*/
  static void splitConnectivity(const K_FLD::IntArray& connectBin, const std::set<BoundaryType>& cuttingEntities,
                     std::vector<K_FLD::IntArray> & connectBout);

  ///
  static void splitConnectivity(const K_FLD::IntArray& connectBin, const std::set<BoundaryType>& cuttingEntities,
                                std::vector< std::vector<E_Int> > & scolors);

  /// Splits and compacts the coordinates for each connectivity bits.
  static void splitCoordinates(std::vector<K_FLD::IntArray>& connects, const K_FLD::FloatArray& pos,
                               std::vector<K_FLD::FloatArray>& poss);

private:

  ContourSplitter(void){}

  ~ContourSplitter(void){}
};

#include "ContourSplitter.cxx"

#endif
