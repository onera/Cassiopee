/*    
    Copyright 2013-2026 ONERA.

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
// -- build edge map of unstructured zone --
#include <utility>
#include <unordered_map>
std::unordered_map<size_t, std::pair> edge;

E_Int np = zonep->np;
#define setEdge(n1,n2) \
    if (n1 < n2) edge[n1+n2*np] = {n1,n2}; \
    else edge[n2+n1*np] = {n1,n2};

for (size_t nc = 0; nc < zonep->connect.size(); nc++) {

  E_Int eltType = zonep->eltType[nc];
  E_Int* connect = zonep->connect[nc];
  E_Int ne = zonep->nec[nc];
  
  E_Int ne2 = 2*ne;
  E_Int ne3 = 3*ne;
  E_Int nd, l;

  switch (eltType)
  {    
    case 2: // TRI
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;       
        setEdge(n1,n2);

        n1 = connect[i+ne]-1;
        n2 = connect[i+ne2]-1;
        setEdge(n1,n2);
  
        n1 = connect[i+ne2]-1;
        n2 = connect[i]-1;
        setEdge(n1,n2);
      }
      break;
  }
}
