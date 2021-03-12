/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

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
