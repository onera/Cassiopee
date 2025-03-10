/*    
    Copyright 2013-2025 Onera.

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

//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef __K_CONTAINERS_DEF_H__
#define __K_CONTAINERS_DEF_H__

#include "Nuga/include/DynArray.h"
#include "Nuga/include/Edge.h"
#include <vector>
#include <map>

namespace NUGA
{

typedef   E_Int                                         size_type;

typedef   std::vector<bool>                             bool_vector_type;
typedef   std::vector<size_type>                        int_vector_type;
typedef   std::set<size_type>                           int_set_type;
typedef   std::pair<size_type, size_type>               int_pair_type;

typedef   std::vector<int_pair_type>                    int_pair_vector_type;
typedef   std::set< int_pair_type >                     int_pair_set_type;

typedef   std::set<K_MESH::Edge>                        oriented_edge_set_type;
typedef   std::set<K_MESH::NO_Edge>                     non_oriented_edge_set_type;
typedef   std::set<K_MESH::NO_Edge>                     non_oriented_int_pair_set_type;

}

#endif
