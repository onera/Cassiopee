/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef _TSSPLITTER_H_
#define _TSSPLITTER_H_

#include "Nuga/include/DynArray.h"
#include <vector>

class TSSplitter
{
  public:
    /* Splits the connectivity (only) with the edge N0N1 and stores the 
       bits into separate containers.
       WARNING: it is the responsibility of the caller to destroy connectBout.
    */
    static void split(
      const K_FLD::IntArray& connectBin, const K_FLD::IntArray& polyLine,
      std::vector<K_FLD::IntArray> & connectBout);

  /* Splits the connectivity with the edge N0N1 and stores the connectivity 
     bits and corresponding coordinates into separate containers.
     WARNING: it is the responsibility of the caller to destroy posOut and 
     connectBout. */
    static void split(
      const K_FLD::FloatArray& pos, const K_FLD::IntArray& connectBin, 
      const K_FLD::IntArray& polyLine,
      std::vector<K_FLD::FloatArray> & posOut, 
      std::vector<K_FLD::IntArray> & connectBout);

  private:
    TSSplitter(void){}
    ~TSSplitter(void){}
};

#endif
