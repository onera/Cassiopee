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
#include "connect.h"
#include <stack>
#include <unordered_map>

E_Int K_CONNECT::colorConnexParts(E_Int *neis, E_Int *xadj, E_Int nelts, E_Int *colors)
{
  memset(colors, -1, nelts*sizeof(E_Int));
  std::stack<E_Int> pool;
  E_Int seed = 0, color = 0;

  while (1) {
    while (seed < nelts && colors[seed] != -1)
      seed++;

    if (seed >= nelts)
      return color;

    pool.push(seed);

    while (!pool.empty()) {
      E_Int elem = pool.top();
      pool.pop();
      
      if (colors[elem] != -1)
        continue;

      colors[elem] = color;

      E_Int start = xadj[elem];
      E_Int stride = xadj[elem+1] - start;
      E_Int *pn = &neis[start];
      for (E_Int i = 0; i < stride; i++) {
        E_Int nei = pn[i];
        if (nei != -1 && colors[nei] == -1)
          pool.push(nei);
      }
    }

    color++;
  }
}
