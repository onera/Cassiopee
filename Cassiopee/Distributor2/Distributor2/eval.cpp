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

# include "distributor2.h"
# include "kcore.h"
using namespace std;
using namespace K_FLD;

//=============================================================================
// Evalue la distribution dis
//=============================================================================
E_Float K_DISTRIBUTOR2::eval(
  E_Int nb, E_Int NProc, E_Float meanPtsPerProc,
  vector<E_Float>& solver, vector<E_Float>& latence,
  vector<E_Float>& comSpeed, E_Int* com, E_Int* comd, E_Int sizeComd,
  FldArrayF& nbPtsPerProcs, vector<E_Float>& nbPts,
  E_Int* dis)
{
  nbPtsPerProcs.setAllValuesAtNull();
  E_Float* nbPtsPerProcsp = nbPtsPerProcs.begin();

  E_Int i, k, p, proci, prock;
  E_Float res, volcom, lati, speedi;

  for (i = 0; i < nb; i++)
  {
    nbPtsPerProcsp[dis[i]] += nbPts[i];
  }

  res = 0.;
  for (p = 0; p < NProc; p++)
    res += solver[p] * K_FUNC::E_abs(nbPtsPerProcsp[p]-meanPtsPerProc);

  // avec com
  if (com != NULL)
  {
    for (i = 0; i < nb; i++)
    {
      proci = dis[i];
      lati = latence[proci];
      speedi = comSpeed[proci];
      for (k = 0; k < nb; k++)
      {
        volcom = com[k + i*nb];
        if (volcom > 0)
        {
          prock = dis[k];
          // le voisin est-il sur le meme processeur?
          if (proci != prock)
          {
            res += lati + speedi*volcom;
          }
        }
      }
    }
  }

  // avec comd
  if (comd != NULL)
  {
    E_Int v1;
    for (E_Int v = 0; v < sizeComd/2; v++)
    {
      v1 = comd[2*v]; volcom = comd[2*v+1];
      k = E_Int(v1/nb);
      i = v1-k*nb;
      proci = dis[i];
      lati = latence[proci];
      speedi = comSpeed[proci];
      prock = dis[k];
      // le voisin est-il sur le meme processeur?
      if (proci != prock)
      {
        res += lati + speedi*volcom;
      }
    }
  }

  // Check for empty processors
  for (i = 0; i < NProc; i++)
  {
    if (K_FUNC::E_abs(nbPtsPerProcsp[i]) < 1.e-10) res += 1.e9;
  }
  return res;
}
