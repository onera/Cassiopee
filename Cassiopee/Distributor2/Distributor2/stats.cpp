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

# include "distributor2.h"
# include "kcore.h"
using namespace std;
using namespace K_FLD;

//============================================================================
// Remplit les stats
// IN: nbPts: nbre de points de chaque zone
// IN: out: affectation de zone sur procs
// OUT: empty, varMin, varMax, varRMS, volRatio
//============================================================================
void K_DISTRIBUTOR2::stats(vector<E_Float>& nbPts, E_Int NProc,
  E_Int* com, E_Int* comd, E_Int sizeComd,
  vector<E_Int>& out,
  E_Int& empty, E_Float& varMin, E_Float& varMax, E_Float& varRMS,
  E_Float& volRatio)
{
  E_Int proc, nptsCom;
  E_Int nb = nbPts.size();

  // Nb total de pts: nbTot = 0
  E_Float nbTot = 0;
  for (E_Int i = 0; i < nb; i++) nbTot += nbPts[i];

  // Nb de noeuds moyens devant etre contenus par chaque processeur
  E_Float meanPtsPerProc = nbTot*1./NProc;

  // Calcul du nombre de pts par processeurs
  FldArrayI nbNodePerProc(NProc);
  nbNodePerProc.setAllValuesAtNull();
  for (E_Int i = 0; i < nb; i++)
  {
    proc = out[i];
    nbNodePerProc[proc] += nbPts[i];
  }
  //printf("Nb de pts moyen par proc: " SF_D_ "\n", E_Int(meanPtsPerProc));

  //printf("Nb de pts par proc:\n");
  empty = 0;
  for (E_Int i = 0; i < NProc; i++)
  {
    //  printf("Proc " SF_D_ ": " SF_D_ " pts\n", i, nbNodePerProc[i]);
    if (nbNodePerProc[i] == 0)
    {
      empty = 1;
      printf("Warning: processor " SF_D_ " is empty!\n", i);
    }
  }

  // Variance
  varMin = 1.e6; varMax = 0.; varRMS = 0.;
  for (E_Int i = 0; i < NProc; i++)
  {
    E_Float v = K_FUNC::E_abs(nbNodePerProc[i] - meanPtsPerProc);
    varMin = K_FUNC::E_min(varMin, v);
    varMax = K_FUNC::E_max(varMax, v);
    varRMS = varRMS + v*v;
  }
  varMin = varMin / meanPtsPerProc;
  varMax = varMax / meanPtsPerProc;
  varRMS = sqrt(varRMS) / (NProc*meanPtsPerProc);
  printf("Info: varMin=" SF_F_ "%%, varMax=" SF_F_ "%%, varRMS=" SF_F_ "%%\n",
         100*varMin, 100*varMax, 100*varRMS);

  nptsCom = 0; E_Int volTot = 0;
  if (com != NULL)
  {
    for (E_Int i = 0; i < nb; i++)
    {
      E_Int proci = out[i];
      for (E_Int k = 0; k < nb; k++)
      {
        if (com[k+i*nb] > 0)
        {
          E_Int prock = out[k];
          volTot += com[k + i*nb];
          // le voisin est-il sur le meme processeur?
          if (proci != prock)
          {
            nptsCom += com[k + i*nb];
          }
        }
      }
    }
  }

  if (comd != NULL)
  {
    E_Int v1, volcom, proci, prock;
    E_Int i,k;
    for (E_Int v = 0; v < sizeComd/2; v++)
    {
      v1 = comd[2*v]; volcom = comd[2*v+1];
      k = E_Int(v1/nb);
      i = v1-k*nb;
      proci = out[i];
      prock = out[k];
      volTot += volcom;
      // le voisin est-il sur le meme processeur?
      if (proci != prock)
      {
        nptsCom += volcom;
      }
    }
  }

  //printf("Volume de communication=" SF_D_ "\n", nptsCom);
  if (volTot > 1.e-6) volRatio = E_Float(nptsCom)/E_Float(volTot);
  else volRatio = 0.;
  printf("Info: external com ratio=" SF_F_ "%%\n", volRatio*100);
  fflush(stdout);
}
