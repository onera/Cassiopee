/*    
    Copyright 2013-2024 Onera.

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

#include "Nuga/include/Hexahedron.h"
//#include "Nuga/include/macros.h"
//#include "Nuga/include/Triangle.h"
//#include "../Def/DefCplusPlusConst.h"

namespace K_MESH
{
void Hexahedron::triangulate(E_Int* target)
{
  // WARNING: connectT3 is Apended (not cleared upon entry)
  
  K_MESH::Quadrangle q4;
  
  for (size_t f = 0; f < 6; ++f)
  {
    getBoundary(f, q4);
    q4.triangulate(q4.nodes(), target + f*6);
  }
}

void Hexahedron::get_edges(const E_Int* nodes, Vector_t<K_MESH::NO_Edge>& edges)
{
  // nodes length is 8 and there are 12 edges in total
  E_Int nb_nodes = 4;

  for (int i = 0; i < nb_nodes; ++i)
  {
	edges[i] = K_MESH::NO_Edge(nodes[i], nodes[(i+1)%nb_nodes]); // 0 to 3 : first PG : nb_nodes edges
	edges[nb_nodes+i] = K_MESH::NO_Edge(nodes[nb_nodes+i], nodes[nb_nodes+(i+1)%nb_nodes]); // 4 to 7 : second PG : nb_nodes edges
	edges[2*nb_nodes+i] = K_MESH::NO_Edge(nodes[i], nodes[nb_nodes+i]); // 8 to 11 : the nb_nodes remaining edges
  }
}

bool Hexahedron::cross(const ngon_type& ng, const K_FLD::FloatArray& crd, const E_Int* face, E_Int nb_faces, K_FLD::FloatArray& data, E_Float* P0, E_Float* P1, E_Float& lambda0, E_Float& lambda1, E_Float tolerance)
{
  // crossing points are defined by lambda coefficients, initialize them
  lambda0 = NUGA::FLOAT_MAX;
  lambda1 = NUGA::FLOAT_MAX;
  
  for (int i = 0; i < nb_faces; ++i)
  {
	E_Int PGi = face[i] - 1; // face comes from the NGON
	
	// 3 points of the PGi
	const E_Int* node = ng.PGs.get_facets_ptr(PGi);
	const E_Float* q0 = crd.col(node[0]-1);
	const E_Float* q1 = crd.col(node[1]-1);
	const E_Float* q2 = crd.col(node[2]-1);
	
	// intersect
	E_Bool overlap = false;
	E_Int tx = 0;
	E_Float l0, l1;
	bool x = K_MESH::Triangle::intersect<3>(q0, q1, q2, P0, P1, EPSILON, true, l0, l1, tx, overlap);
	
	if (!x) continue; // doesn't intersect
	
	if (overlap) continue; // both l0 & l1 need a value
	
	if ( (l0 <= tolerance) || (l0 >= 1 - tolerance) ) continue;
	
	if (lambda0 == NUGA::FLOAT_MAX) lambda0 = l0;
	else
	  lambda1 = l0;
	
	if ( (l1 <= tolerance) || (l1 >= 1 - tolerance) ) continue;
	
	if (lambda1 == NUGA::FLOAT_MAX) lambda1 = l1;
	else
	  lambda0 = l1;
  }
  
  // end 
  if ( (lambda0 != NUGA::FLOAT_MAX) && (lambda1 != NUGA::FLOAT_MAX) ) return true;
	
  return false;
}
  




}
