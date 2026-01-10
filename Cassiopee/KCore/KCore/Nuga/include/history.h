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
//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef NUGA_HISTORY_H
#define NUGA_HISTORY_H

#include<vector>
#include "Nuga/include/defs.h"

namespace NUGA
{
  struct morse_t
  {
    std::vector<E_Int> data, xr;
    
    morse_t() { xr.push_back(0); };
    
    void clear() { data.clear(); xr.clear(); xr.push_back(0); }
    void append(E_Int val) { data.push_back(val); xr.push_back(data.size());/*next pos*/}
    void append(std::vector<E_Int>& vals) { data.insert(data.end(), ALL(vals)); xr.push_back(data.size());}
  };

  struct history_t
  {
    std::vector<E_Int> ptoids, ptnids;
    morse_t pgoids, phoids, pgnids, phnids;

    void clear() { ptoids.clear();  pgoids.clear();  phoids.clear();  ptnids.clear(); pgnids.clear();  phnids.clear(); }

    template<typename T>
    void transfer_pg_colors(std::vector<E_Int>& src_ids, std::vector<T> & src_flags, std::vector<T>& target_flags)
    {
      //compute target size : find max id
      E_Int target_sz = 0;
      for (E_Int i=0; i < (E_Int)src_ids.size(); ++i)
      {
        E_Int srcid = src_ids[i];

        E_Int nb_target = pgnids.xr[srcid+1] - pgnids.xr[srcid];
        const E_Int* tgt_start = &pgnids.data[pgnids.xr[srcid]];
        for (E_Int t=0; t < nb_target; ++t)
        {
          if (tgt_start[t] == IDX_NONE) continue;
          if (tgt_start[t] < 0) target_sz = std::max(target_sz, -(tgt_start[t]+1));
          else target_sz = std::max(target_sz, tgt_start[t]);
        }
      }

      target_flags.resize(target_sz + 1, IDX_NONE);
      for (E_Int i=0; i < (E_Int)src_ids.size(); ++i)
      {
        E_Int srcid = src_ids[i];

        E_Int nb_target = pgnids.xr[srcid+1] - pgnids.xr[srcid];
        const E_Int* tgt_start = &pgnids.data[pgnids.xr[srcid]];

        if (nb_target == 1)
        {
          if (*tgt_start == IDX_NONE) // is gone
            continue;
          else if (*tgt_start < 0)    //  is agglomerated into target
          {
            E_Int tgtid = -(*tgt_start+1);
            if (target_flags[tgtid] == -1) //frozen because contrinutions have different flags
              continue;
            else if (target_flags[tgtid] == IDX_NONE)
              target_flags[tgtid] = src_flags[srcid];
            else if (target_flags[tgtid] != src_flags[srcid]) //freeze
              target_flags[tgtid] = -1;
          }
          else // just a move 
            target_flags[*tgt_start] = src_flags[srcid];
        }
        else // it's a split : all children inherits parent value
        {
          for (E_Int t=0; t < nb_target; ++t)
            target_flags[tgt_start[t]] = src_flags[srcid];
        }
      }
    }

  }; //history_t
}

#endif
