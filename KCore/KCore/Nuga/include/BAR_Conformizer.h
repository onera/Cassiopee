/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef __BAR_CONFORMIZER_H__
#define	__BAR_CONFORMIZER_H__

#include "Conformizer.h"

#include "Nuga/include/defs.h"
#include <vector>

namespace NUGA
{

template <E_Int DIM>
class BAR_Conformizer : public Conformizer<DIM, K_MESH::Edge> {

public :
  typedef Conformizer<DIM, K_MESH::Edge>  parent_type;
  
public:
  BAR_Conformizer(bool whisto = false) :parent_type(whisto) {}

  virtual ~BAR_Conformizer(){}
  
  std::vector<std::pair<int, int>> get_x_history();

  // Overridden Methods ///////////////////////////////////////////////////////

protected:
  
  ///
  void __set_tolerances(E_Float Lmin, E_Float Lmax, E_Float  user_tolerance);

  ///
  E_Int __intersect(K_FLD::FloatArray& pos, const K_FLD::IntArray& connect,
                     E2& e1, E2& e2, E_Float tol);
  ///
  void __update_data(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const std::vector<E_Int>& newIDs);
  /// Splits the triangles by triangulation.
  E_Int __split_Elements(const K_FLD::FloatArray& pos, K_FLD::IntArray & connect,
                        NUGA::bool_vector_type& xc,
                        NUGA::int_vector_type& ancestors);

#ifdef DEBUG_CONFORMIZER
  ///
  void drawElements(const char* fname, const char* filefmt, const K_FLD::FloatArray& coord, const K_FLD::IntArray& connect, const std::vector<E2> & elts, bool localid = false, std::vector<E_Int>* colors = 0)
  {
    K_FLD::IntArray cB;
    for (size_t i = 0; i < elts.size(); ++i)
    {
      cB.pushBack(connect.col(elts[i].id), connect.col(elts[i].id)+2);
    }
    K_FLD::FloatArray crd = coord;
    crd.resize(3, coord.cols(), 0.);
    medith::write(fname, crd, cB, "BAR");

  }
#endif
    
private:
  
  void __reorder_nodes_on_edge(const K_FLD::FloatArray& pos, std::vector<E_Int>& nodes, const E_Float *P0);
  
  BAR_Conformizer(const BAR_Conformizer& orig){}
  
private:

};
}

#include "BAR_Conformizer.cxx"

#endif	/* BAR_CONFORMIZER_H */

