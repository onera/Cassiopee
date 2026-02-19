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

#ifndef NUGA_GHOST_H
#define NUGA_GHOST_H

#include "Nuga/include/ngon_t.hxx"

using ngon_type = ngon_t<K_FLD::IntArray>;

#define NEIGHBOR(PH0, F2E, shift, Fi) ((F2E[Fi] == PH0 ) ? F2E[Fi+shift] : F2E[Fi])

namespace NUGA
{
  
static E_Int append_block_with_ghosts(ngon_type& ng, K_FLD::FloatArray& crd, E_Int nb_neighblocks, const E_Int* sz, const E_Int* F2Esep, E_Int** pointList, const ngon_type* ngD, K_FLD::FloatArray** crdD, E_Int** F2EsepD, E_Int** pointListD, E_Int idx_start)
{ 
  // WARNING : expect a F2E formated as Cassiopee : non-interleaved, idx start at 1, NONE is 0
    
  std::vector<bool> flagFirstL, flagSecL;
  ngon_type ngFirstL, ngSecondL;
    
  ngon_type ngtmp;
  std::vector<E_Int> ngpids;

  for (size_t b = 0; b < nb_neighblocks; ++b)
  {
    //current block
    const ngon_type& ngb = ngD[b];
    const E_Int* f2eD    = F2EsepD[b];
    const E_Int* ptListD = pointListD[b];
      
    E_Int nb_phsb = ngb.PHs.size();
    E_Int nb_pgsb = ngb.PGs.size();
    
    //std::cout << "neighbor block size : " << nb_phsb << " elemnts" << std::endl;
    
    flagFirstL.clear();
    flagFirstL.resize(nb_phsb, false);
    flagSecL.clear();
    flagSecL.resize(nb_phsb, false);
    
    E_Int nb_pgf = sz[b];
    
    // first layer
    
    for (size_t i=0; i < nb_pgf; ++i)
    {
      const E_Int & PGiD = ptListD[i];
              
      E_Int PHi = f2eD[PGiD - idx_start] - idx_start;
      if (PHi == -idx_start) // meaning that orienttation is wrong , so consider the other one
        PHi = f2eD[nb_pgsb + PGiD - idx_start] - idx_start;

#ifdef DEBUG_GHOST
      if (PHi < 0 || ( PHi >= nb_phsb))
      {
        std::cout << "wrong F2E : " << std::endl;
        std::cout << "PGi : " << PGiD - idx_start << " (0-based)" << std::endl;
        std::cout << "one side : " << f2eD[PGiD - idx_start] - idx_start << " (0-based)" << std::endl;
        std::cout << "other side : " << f2eD[nb_pgsb + PGiD - idx_start] - idx_start << " (0-based)" << std::endl;
        
      }
#endif
        
      flagFirstL[PHi] = true;
    }
    
    // second layer
            
    //E_Int count(0);
    for (size_t i=0; i < nb_phsb; ++i)
    {
      if (flagFirstL[i] == false) continue; // consider only first layer element
      
      //++count;
      
      const E_Int * pPGi = ngb.PHs.get_facets_ptr(i);
      E_Int nb_faces = ngb.PHs.stride(i);
      //std::cout << "nb faces : " << nb_faces << std::endl;
      for (size_t n=0; n < nb_faces;++n)
      {
        E_Int PGi = *(pPGi + n) - idx_start;
        E_Int PHn = NEIGHBOR(i + idx_start, f2eD, nb_pgsb, PGi) - idx_start;
        
        //std::cout << "neighbor : " << PHn << std::endl; 
        
        if (PHn == -idx_start) continue; // NONE

#ifdef DEBUG_GHOST
         if (PHn < 0 || PHn >= nb_phsb){
            std::cout << "wrong F2E : " << std::endl;
            std::cout << "PGi : " << PGi - idx_start << std::endl;
            std::cout << "left : " << PHn << std::endl;
          }
#endif

        if (flagFirstL[PHn] == true) continue;
          
        flagSecL[PHn] = true;
      }
    }

    //
    ngtmp.clear();
    ngpids.clear();
    ngon_type::select_phs(ngb, flagFirstL, ngpids, ngtmp);

    {
      K_FLD::FloatArray crdt(*crdD[b]);
      ngon_type::compact_to_used_nodes(ngtmp.PGs, crdt);

      E_Int shft = crd.cols();
      ngtmp.PGs.shift(shft);

      crd.pushBack(crdt);
      ngFirstL.append(ngtmp);
    }

    //
    ngtmp.clear();
    ngpids.clear();
    ngon_type::select_phs(ngb, flagSecL, ngpids, ngtmp);

    {
      K_FLD::FloatArray crdt(*crdD[b]);
      ngon_type::compact_to_used_nodes(ngtmp.PGs, crdt);

      E_Int shft = crd.cols();
      ngtmp.PGs.shift(shft);

      crd.pushBack(crdt);
      ngSecondL.append(ngtmp);
    }

  } //block loop
   
#ifdef DEBUG_GHOST
  std::cout << "first lay has : " <<  ngFirstL.PGs.size() << " PGs" << std::endl;
  std::cout << "second lay has : " <<  ngSecondL.PGs.size() << " PGs" << std::endl;
  std::cout << "first lay has : " <<  ngFirstL.PHs.size() << " PHs" << std::endl;
  std::cout << "second lay has : " <<  ngSecondL.PHs.size() << " PHs" << std::endl;
#endif
  
  
  ng.append(ngFirstL);
  ng.append(ngSecondL);
  
} // end
  
#ifdef DEBUG_GHOST
  static E_Int check_inputs(ngon_type& ng, K_FLD::FloatArray& crd, E_Int nb_neighblocks, const E_Int* sz, const E_Int* F2Esep, E_Int** pointList, const ngon_type* ngD, K_FLD::FloatArray** crdD, E_Int** F2EsepD, E_Int** pointListD, E_Int idx_start)
  {
    //1. check F2E : good values (1-based, 0 is for NONE)
    E_Int nb_phs = ng.PHs.size();
    E_Int nb_pgs = ng.PGs.size();

    if (F2Esep == NULL) 
    {
      std::cout << "F2Esep is missing" << std::endl;
      return 1;
    }

    //
    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      if (F2Esep[i] < 0 || F2Esep[i] > nb_phs)
      {
        std::cout << "F2Esep contains a wrong value (" << F2Esep[i] << ") at pos " << i << std::endl;
        return 1;
      }
    }
    
    //2. check sz : good values : strictly positive
    for (E_Int i = 0; i < nb_neighblocks; ++i)
    {
      if (sz[i] <= 0)
      {
        std::cout << "sz contains a wrong value (" << sz[i] << ") at pos " << i << std::endl;
        return 1;
      }
    }
    
    // check block data
    for (size_t b = 0; b < nb_neighblocks; ++b)
    {
      //current block
      const ngon_type& ngb = ngD[b];
      
      //ngb.is_consistent(10.e6);
      
      const E_Int* f2eD    = F2EsepD[b];
      const E_Int* ptListD = pointListD[b];
      
      E_Int nb_phsb = ngb.PHs.size();
      E_Int nb_pgsb = ngb.PGs.size();
      std::cout << "current block pgs/phs : " << nb_pgsb << "/" << nb_phsb << std::endl;
      
      E_Int cur_sz = sz[b];
      
      // 1. check F2EsepD sanity : good values (1-based, 0 is for NONE)
      std::cout << "checking f2eD : should extend as twice the nb of pgs which is : " << nb_pgsb << std::endl;
      for (E_Int i = 0; i < nb_pgsb*2; ++i)
      {
        if (f2eD[i] < 0 || f2eD[i] > nb_phsb)
        {
          std::cout << "F2Esep of donnor " << b << " contains a wrong value (" << f2eD[i] << " with nb_phs=" << nb_phsb << ") at pos " << i << std::endl;
          return 1;
        }
      }
      
      // 2. check pointListD sanity
      for (E_Int i = 0; i < cur_sz; ++i)
      {
        if (ptListD[i] < 1 || ptListD[i] > nb_pgsb)
        {
          std::cout << "ptListD of donnor " << b << " contains a wrong value (" << ptListD[i] << ") at pos " << i << std::endl;
          return 1;
        }
      }
      
      // 3. check pointListD and F2EsepD consistency
      for (E_Int i = 0; i < cur_sz; ++i)
      {
        E_Int PGi = ptListD[i] - 1;
        
        if (f2eD[PGi] != 0 && f2eD[PGi + nb_pgsb] != 0)
        {
          std::cout << "ptListD/f2E inconsistency for donnor " << b << " at polygon " << PGi << "(0-based) : elements are :" << f2eD[PGi] << "/" << f2eD[PGi + nb_pgsb] << std::endl;
          return 1;
        }
        
      }
    }
    
    return 0;
    
  }
#endif  
}




#endif /* GHOST_H */

