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
// IN: nbPts: pour chaque bloc, son nombre de pts
// IN: setBlocks: pour chaque bloc, -1 si ce bloc est a equilibrer,
// si c'est une autre valeur, on considere ce bloc comme deja
// place sur ce proc
// IN: NProc: nombre de processeurs
// OUT: out: pour chaque bloc, son no de proc
//=============================================================================
void K_DISTRIBUTOR2::gradient(
  vector<E_Float>& nbPts, vector<E_Int>& setBlocks,
  E_Int NProc, E_Int* com, E_Int* comd, E_Int sizeComd,
  vector<E_Float>& solver,
  vector<E_Float>& latence, vector<E_Float>& comSpeed, E_Int param,
  vector<E_Int>& out, E_Float& meanPtsPerProc, E_Float& varMin,
  E_Float& varMax, E_Float& varRMS, E_Int& nptsCom, E_Float& volRatio,
  E_Float& bestAdapt)
{
  // Nombre total de blocs: nb
  E_Int nb = nbPts.size();

  // Nombre de blocs a equilibrer: nbeq
  E_Int nbeq = 0;
  for (E_Int i = 0; i < nb; i++)
  {
    if (setBlocks[i] < 0) nbeq++;
  }

  // Nb de noeuds des maillages imposes par processeurs : nbNodePerProc
  FldArrayF nbNodePerProc(NProc);
  nbNodePerProc.setAllValuesAtNull();
  E_Int proc;
  for (E_Int i = 0; i < nb; i++)
  {
    proc = setBlocks[i];
    if (proc > 0) nbNodePerProc[proc] += nbPts[i];
  }

  // Nb total de pts: nbTot = 0
  E_Float nbTot = 0;
  for (E_Int i = 0; i < nb; i++) nbTot += nbPts[i];

  // Nb de noeuds moyens devant etre contenu par chaque processeur
  meanPtsPerProc = E_Int(E_Float(nbTot)/E_Float(NProc));

  // Bloc le plus petit
  E_Float nbPtsMin = K_CONST::E_MAX_FLOAT;
  for (E_Int i = 0; i < nb; i++) nbPtsMin = K_FUNC::E_min(nbPts[i], nbPtsMin);
  // Bloc le plus grand
  E_Float nbPtsMax = 0;
  for (E_Int i = 0; i < nb; i++) nbPtsMax = K_FUNC::E_max(nbPts[i], nbPtsMax);

  FldArrayI dis(nb);
  FldArrayI dis1(nb);

  // Init, on repartit les blocs par taille decroissante
  // Puis sur les procs les moins charges (bin-packing)
  FldArrayF nbPtsPerProcs(NProc);
  nbPtsPerProcs.setAllValuesAtNull();
  E_Float* nbPtsPerProcsp = nbPtsPerProcs.begin();

  E_Int* largest = new E_Int [nb];
  vector<E_Float> npts(nbPts);
  for (E_Int i = 0; i < nb; i++)
  {
    E_Float size = -1; E_Int jlarge = 0;
    for (E_Int j = 0; j < nb; j++)
    {
      if (npts[j] > size) {jlarge = j; size = npts[j];}
    }
    largest[i] = jlarge; npts[jlarge] = -1;
  }

  for (E_Int i = 0; i < nb; i++)
  {
    E_Int j = largest[i];
    if (setBlocks[j] < 0)
    {
      if (i < NProc)
      {
        dis1[j] = i;
        nbPtsPerProcsp[i] += nbPts[j];
      }
      else
      {
        E_Int kless = 0; E_Float minProc = K_CONST::E_MAX_FLOAT;
        for (E_Int k = 0; k < NProc; k++)
        {
          if (nbPtsPerProcsp[k] < minProc)
          { kless = k; minProc = nbPtsPerProcsp[k]; }
        }
        dis1[j] = kless;
        nbPtsPerProcsp[kless] += nbPts[j];
      }
    }
    else
    {
      dis1[j] = setBlocks[j];
      nbPtsPerProcs[setBlocks[j]] += nbPts[j];
    }
  }

  // Essai best fit avec coms
  /*
  nbPtsPerProcs.setAllValuesAtNull();
  E_Int* alreadySet = new E_Int [nb];
  for (E_Int i = 0; i < nb; i++) alreadySet[i] = 0;
  for (E_Int i = 0; i < nb; i++)
  {
    if (setBlocks[i] >= 0)
    {
      alreadySet[i] = 1; dis[i] = setBlocks[i];
      nbPtsPerProcsp[setBlocks[i]] += nbPts[i];
    }
  }

  for (E_Int i = 0; i < nb; i++)
  {
    E_Int j = largest[i];
    if (alreadySet[j] == 0)
    {
      // on le met sur le proc le moins charge
      E_Int kless = 0; E_Float minProc = K_CONST::E_MAX_FLOAT;
      for (E_Int k = 0; k < NProc; k++)
      {
        if (nbPtsPerProcs[k] < minProc)
        { kless = k; minProc = nbPtsPerProcs[k]; }
      }
      alreadySet[j] = 1;
      dis[j] = kless;
      nbPtsPerProcs[kless] += nbPts[j];

      // on ajoute les blocs en plus grosse com jusqu'a depasser la moyenne
      E_Int kmax = 0; E_Float volComMax; E_Float volCom;
      while (kmax >= 0)
      {
        volComMax = -1.; kmax = -1;
        for (E_Int k = 0; k < nb; k++)
        {
          if (alreadySet[k] == 0)
          {
            volCom = com[k + j*nb];
            if (volCom > volComMax &&
                nbPtsPerProcs[kless] + nbPts[k] <= 1.01*meanPtsPerProc)
            { volComMax = volCom; kmax = k; }
          }
        }
        if (kmax >= 0)
        {
          alreadySet[kmax] = 1;
          dis[kmax] = kless;
          nbPtsPerProcs[kless] += nbPts[kmax];
        }
      }
    }
  }
  delete [] alreadySet;
  */
  delete [] largest;

  E_Float evalp1 = K_DISTRIBUTOR2::eval(nb, NProc, meanPtsPerProc,
                                        solver, latence,
                                        comSpeed, com, comd, sizeComd,
                                        nbPtsPerProcs, nbPts,
                                        dis1.begin());

  dis = dis1;
  E_Float evalp = evalp1;
  /*
  E_Float evalp = K_DISTRIBUTOR2::eval(nb, NProc, meanPtsPerProc,
                                       solver, latence,
                                       comSpeed, com, comd, sizeComd,
                                       nbPtsPerProcs, nbPts,
                                       dis.begin());
  if (evalp1 < evalp)
  { dis = dis1; evalp1 = evalp; }
  */

  //printf("Adaptation init: %f\n", evalp);

  if (param == 1)
  {
    // Shuffle
    E_LONG idum = -1;
    E_Float ca = -1./K_FUNC::E_max(nbPtsMax-nbPtsMin, 1.);
    E_Float cb = 0.5;
    E_Float diff; E_Int temp;
    for (E_Int b = 0; b < nb; b++)
    {
      for (E_Int k = b+1; k < nb; k++)
      {
        if (dis[b] != dis[k] &&
            setBlocks[b] < 0 && setBlocks[k] < 0)
        {
          diff = K_FUNC::E_abs(nbPts[b] - nbPts[k]);
          // chance de muter: inv. proportionelle diff
          if (K_NOISE::stdRand(&idum) < ca*diff+cb)
          {
            temp = dis[b];
            dis[b] = dis[k];
            dis[k] = temp;
          }
        }
      }
    }
    evalp = K_DISTRIBUTOR2::eval(nb, NProc, meanPtsPerProc,
                                 solver, latence,
                                 comSpeed, com, comd, sizeComd,
                                 nbPtsPerProcs, nbPts,
                                 dis.begin());
    //printf("Adaptation randomized: %f\n", evalp);
  }

  const E_Int nitMax = 2;

  //----------------------------------------------------------
  // Gradient
  //----------------------------------------------------------
  E_Int nbItWithFail = 0;
  E_Int temp; E_Float evaln;
  E_Int* disp = dis.begin();
  while (nbItWithFail < nitMax)
  {
    for (E_Int b = 0; b < nb; b++)
    {
      if (setBlocks[b] < 0)
      {
        for (E_Int k = b+1; k < nb; k++)
        {
          if (setBlocks[k] < 0)
          {
            temp = disp[b];
            disp[b] = disp[k];
            disp[k] = temp;

            // Recherche du bloc swap qui baisse F
            evaln = K_DISTRIBUTOR2::eval(nb, NProc, meanPtsPerProc,
                                         solver, latence,
                                         comSpeed, com, comd, sizeComd,
                                         nbPtsPerProcs, nbPts,
                                         disp);
            if (evaln < evalp)
            {
              evalp = evaln; nbItWithFail = 0;
            }
            else // revert
            {
              temp = disp[b];
              disp[b] = disp[k];
              disp[k] = temp;
            }
          }
        }
      }
    }
    //printf("Adaptation: %f\n", evalp);
    nbItWithFail++;
  }

  // Sortie
  for (E_Int i = 0; i < nb; i++) out[i] = dis[i];

  bestAdapt = evalp;
  //printf("Adaptation: %f\n", bestAdapt);

  // external stats
  E_Int empty;
  stats(nbPts, NProc, com, comd, sizeComd, out, empty,
        varMin, varMax, varRMS, volRatio);

}
