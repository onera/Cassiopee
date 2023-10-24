#include "metric.h"
#include <stack>
#include <unordered_map>

E_Int K_METRIC::colorConnexParts(E_Int *neis, E_Int *xadj, E_Int nelts, E_Int *colors)
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
