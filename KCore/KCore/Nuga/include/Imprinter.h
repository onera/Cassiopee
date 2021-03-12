/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef __DELAUNAY_IMPRINTER_H__
#define __DELAUNAY_IMPRINTER_H__

#include "Nuga/include/defs.h"
#include "Nuga/include/Triangle.h"

namespace DELAUNAY
{

  class Imprinter
  {

  public:
    typedef   K_MESH::Triangle            element_type;
    typedef   NUGA::int_vector_type int_vector_type;
    typedef   NUGA::size_type       size_type;

  public:
    ///
    Imprinter(const K_FLD::FloatArray& posS,
              const K_FLD::IntArray& connectT3);
    
    ///
    ~Imprinter(){}

    ///
    void run(const K_FLD::FloatArray& posC, const K_FLD::IntArray& connectB0,
             K_FLD::FloatArray& posUV, int_vector_type& cell_indices, std::vector<E_Int>* hard_nodes = 0);

  private:
    ///
    const K_FLD::FloatArray& _posS;
    ///
    const K_FLD::IntArray& _connectT3; 
  };
}

#endif

