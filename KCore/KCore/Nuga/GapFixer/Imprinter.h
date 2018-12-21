/*    
    Copyright 2013-2019 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __DELAUNAY_IMPRINTER_H__
#define __DELAUNAY_IMPRINTER_H__

#include "Def/DefContainers.h"
#include "MeshElement/Triangle.h"

namespace DELAUNAY
{

  class Imprinter
  {

  public:
    typedef   K_MESH::Triangle            element_type;
    typedef   K_CONT_DEF::int_vector_type int_vector_type;
    typedef   K_CONT_DEF::size_type       size_type;

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

