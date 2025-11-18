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
# include <stdlib.h>

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
void K_DISTRIBUTOR2::genetic(
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
  //for (E_Int i = 0; i < nb; i++) printf("bloc %d : %g\n", i, nbPts[i]);

  // Nombre de blocs a equilibrer: nbeq
  E_Int nbeq = 0;
  for (E_Int i = 0; i < nb; i++)
  {
    if (setBlocks[i] < 0) nbeq++;
  }

  // Nb de noeuds des maillages imposes par processeurs: nbNodePerProc
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

  // Nb de noeuds moyens devant etre contenus par chaque processeur
  meanPtsPerProc = nbTot*1./NProc;

  // Bloc le plus petit
  E_Float nbPtsMin = K_CONST::E_MAX_FLOAT;
  for (E_Int i = 0; i < nb; i++) nbPtsMin = K_FUNC::E_min(nbPts[i], nbPtsMin);
  // Bloc le plus grand
  E_Float nbPtsMax = 0;
  for (E_Int i = 0; i < nb; i++) nbPtsMax = K_FUNC::E_max(nbPts[i], nbPtsMax);

  //----------------------------
  // Data for genetic algorithm
  //----------------------------
  // Size of population
  E_Int sizeOfPopulation = K_FUNC::E_max(E_Int(30), 2*nb);
  sizeOfPopulation = K_FUNC::E_min(sizeOfPopulation, 300);
  // Number max of generations
  E_Int nitMax = 50;
  if (param == 1) { nitMax = 0; sizeOfPopulation = 4; } // fast

  // Local data
  E_LONG idum = -1;

  //----------------------------------------------------------
  // Genetic algorithm
  //----------------------------------------------------------
  E_Int nbItWithFail = 0;
  E_Float bestEval = 1.E+12;
  FldArrayI pop(nb*sizeOfPopulation);
  E_Int* popp = pop.begin();
  FldArrayF evalPop(sizeOfPopulation);
  E_Float* evalp = evalPop.begin();
  FldArrayI survivors(sizeOfPopulation);
  E_Int* survivorsp = survivors.begin();

  // Initialization:
  // Pour le premier,
  // on essaie de mettre un truc pas trop con.
  FldArrayF nbPtsPerProcs(nbNodePerProc);
  nbPtsPerProcs.setAllValuesAtNull();
  E_Float* nbPtsPerProcsp = nbPtsPerProcs.begin();
  E_Int indBest = 0;
  E_Int* pop1 = popp;
  for (E_Int i = 0; i < nb; i++)
  {
    if (setBlocks[i] < 0)
    {
      E_Float bestEval = meanPtsPerProc;
      E_Float nbPtsForMesh = nbPts[i];
      for (E_Int p = 0; p < NProc; p++)
      {
        E_Float eval =
          K_FUNC::E_abs(meanPtsPerProc-nbPtsForMesh-nbPtsPerProcsp[p]);
        if (eval < bestEval)
        {
          bestEval = eval;
          indBest  = p;
        }
      }
      nbPtsPerProcsp[indBest] += nbPtsForMesh;
      pop1[i] = indBest;
    }
    else pop1[i] = setBlocks[i];
  }

  // Pour le deuxieme, on repartit les blocs les uns apres les autres
  E_Int* pop2 = popp+nb;
  for (E_Int i = 0; i < nb; i++)
  {
    if (setBlocks[i] < 0) pop2[i] = i%NProc;
    else pop2[i] = setBlocks[i];
  }

  // Pour le troisieme, on repartit les blocs par taille decroissante
  // Algo First Fit
  E_Int* largest = new E_Int [nb];
  vector<E_Float> npts(nbPts);
  E_Int* pop3 = popp+2*nb;
  for (E_Int i = 0; i < nb; i++)
  {
    E_Float size = -1; E_Int jlarge = 0;
    for (E_Int j = 0; j < nb; j++)
    {
      if (npts[j] > size) { jlarge = j; size = npts[j]; }
    }
    largest[i] = jlarge; npts[jlarge] = -1;
  }
  for (E_Int i = 0; i < nb; i++)
  {
    E_Int j = largest[i];
    if (setBlocks[j] < 0) pop3[j] = i%NProc;
    else pop3[j] = setBlocks[j];
  }

  // Pour le quatrieme, on repartit les blocs par taille decroissante
  // Puis sur les procs les moins charges (Best fit)
  nbPtsPerProcs = nbNodePerProc;
  nbPtsPerProcsp = nbPtsPerProcs.begin();
  nbPtsPerProcs.setAllValuesAtNull();
  E_Int* pop4 = popp+3*nb;
  for (E_Int i = 0; i < nb; i++)
  {
    E_Int j = largest[i];
    if (setBlocks[j] < 0)
    {
      if (i < NProc)
      { pop4[j] = i; nbPtsPerProcs[i] += nbPts[j]; }
      else
      {
        E_Int kless = 0; E_Float minProc = K_CONST::E_MAX_FLOAT;
        for (E_Int k = 0; k < NProc; k++)
        {
          if (nbPtsPerProcs[k] < minProc)
          { kless = k; minProc = nbPtsPerProcs[k]; }
        }
        pop4[j] = kless;
        nbPtsPerProcs[kless] += nbPts[j];
      }
    }
    else
    {
      pop4[j] = setBlocks[j];
      nbPtsPerProcs[setBlocks[j]] += nbPts[j];
    }
  }

  // Pour le cinquieme, best fit avec utilisation des coms
  /*
  E_Int* pop5 = popp+4*nb;
  E_Int* alreadySet = new E_Int [nb];
  for (E_Int i = 0; i < nb; i++) alreadySet[i] = 0;
  nbPtsPerProcs.setAllValuesAtNull();

  for (E_Int i = 0; i < nb; i++)
  {
    if (setBlocks[i] >= 0)
    {
      alreadySet[i] = 1; pop5[i] = setBlocks[i];
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
      pop5[j] = kless;
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
          pop5[kmax] = kless;
          nbPtsPerProcs[kless] += nbPts[kmax];
        }
      }
    }
  }
  delete [] alreadySet;
  */

  //for (E_Int i = 0; i < nb; i++)
  //  printf("%d %d %d\n", i, pop(i, 3), nbPts[i]);
  delete [] largest;

  // Pour les autres, c'est completement au hasard:
  for (E_Int j = 5; j <= sizeOfPopulation; j++)
  {
    E_Int* popj = popp+(j-1)*nb;
    for (E_Int i = 0; i < nb; i++)
    {
      if (setBlocks[i] < 0)
      {
        E_Int targetProc = E_Int(NProc*K_NOISE::stdRand(&idum));
        popj[i] = targetProc;
      }
      else popj[i] = setBlocks[i];
    }
  }

  // First evaluation for nice people:
//#pragma omp parallel for
  for (E_Int j = 1; j <= sizeOfPopulation; j++)
  {
    evalp[j-1] = K_DISTRIBUTOR2::eval(nb, NProc, meanPtsPerProc,
                                      solver, latence,
                                      comSpeed, com, comd, sizeComd,
                                      nbPtsPerProcs, nbPts,
                                      popp+(j-1)*nb);
  }

  //for (E_Int j = 1; j <= sizeOfPopulation; j++) printf("evalp=%f\n", evalp[j-1]);

  //printf("Adaptation init1: %f\n", evalp[0]);
  //printf("Adaptation init2: %f\n", evalp[1]);
  //printf("Adaptation init3: %f\n", evalp[2]);
  //printf("Adaptation init4: %f\n", evalp[3]);
//   for (E_Int j = 1; j < sizeOfPopulation; j++)
//   {
//     for (E_Int i = 0; i < nb; i++) printf("%d ", pop(i,j));
//     printf(" : %f\n", evalp[j-1]);
//   }

  bestEval = K_CONST::E_MAX_FLOAT;
  E_Int jBest;
  jBest = 1;
  for (E_Int j = 1; j <= sizeOfPopulation; j++)
    jBest = (evalp[j-1] < evalp[jBest-1] ? j : jBest);
  //printf("Adaptation: %f\n", evalp[jBest-1]);

  while (nbItWithFail < nitMax)
  {
    //-------------------------------------------------------------------------
    // Phase of selection:
    //-------------------------------------------------------------------------
    // On tire au hasard pour voir quelle partie de la population survit:
    // La probabilite est de: evalBest/evalCurr. Le personne la mieux
    // adaptee a sa survie assuree:
    E_Float oldBest = bestEval;
    bestEval = evalp[jBest-1];

    if (oldBest <= bestEval)
    {
      // Dans ce cas, on a degenere la population
      nbItWithFail++;
    }
    else
    {
      // Ah, on a une meilleur population:
      nbItWithFail = 0;
    }
    E_Int nbSurvivors = 0;
    if (K_FUNC::fEqualZero(bestEval)) break;
    for (E_Int j = 1; j <= sizeOfPopulation; j++)
    {
      if (K_NOISE::stdRand(&idum) > bestEval/evalp[j-1])
      {
        evalp[j-1] = -1;
      }
      else
      {
        survivorsp[nbSurvivors++] = j;
      }
    }
    //printf("pourcentage de selection: %f ", nbSurvivors*1./nb);

    // Les participants vires sont remplaces par des clones
    // des survivants:
    for (E_Int j = 1; j <= sizeOfPopulation; j++)
    {
      if (evalp[j-1] == -1)
      {
       E_Int indClone = survivorsp[E_Int(nbSurvivors*K_NOISE::stdRand(&idum))];
       for (E_Int i = 0; i < nb; i++)
       {
         popp[i+(j-1)*nb] = popp[i+(indClone-1)*nb];
       }
      }
    }

//     FldArrayI uglyPeople(pop);
//     //-------------------------------------------------------------------------
//     // Croisement de la population de nice people
//     //-------------------------------------------------------------------------
//     for ( E_Int j = 1; j <= sizeOfPopulation; j++ )
//     {
//       if ( j != jBest )
//       {
//         // tirage au hasard de papa et de maman
// 	E_Int indPapa  =
//           survivorsp[E_Int((1.*nbSurvivors*lrand48())/(MAXRAND+1.0))];
// 	E_Int indMaman =
//           survivorsp[E_Int((1.*nbSurvivors*lrand48())/(MAXRAND+1.0))];
// 	for ( E_Int i = 0; i < nb; i++ )
//         {
//           E_Int nbPts1 = nbPts[i];

//           // Recherche du bloc ayant le nb de pt le + proche
//           E_Int kCross = 0; E_Int diff = 1.e6;
//           E_Int nbPts2 = 0; E_Int diffl;
//           for ( E_Int k = 0; k < nb; k++ )
//           {
//             if (k != i)
//             {
//               nbPts2 = nbPts[k];
//               diffl = K_FUNC::E_abs(nbPts1 - nbPts2);
//               if (diffl < diff) {kCross = k; diff = diffl;}
//             }
//           }
//           diff = 1.e6;
// 	  if ( drand48() >= 1.- diff/K_FUNC::E_max(nbPts1, nbPts2) )
// 	  {
// 	    pop(i, j) = uglyPeople(i, indPapa);
// 	  }
// 	  else
// 	  {
// 	    pop(i, j) = uglyPeople(kCross, indMaman);
// 	  }
//         }
//       }
//       else
//       {
// 	for ( E_Int i = 0; i < nb; i++ )
// 	{
// 	  pop(i, j) = uglyPeople(i, j);
// 	}
//       }
//     }

    //-------------------------------------------------------------------------
    // Aggressive mutation
    //-------------------------------------------------------------------------
    if (nbItWithFail > 10)
    {
      E_Float a = -0.8/(nbPtsMax-nbPtsMin);
      E_Float b = 0.9*nbPtsMax - 0.1*nbPtsMin;
      for (E_Int j = 1; j <= sizeOfPopulation; j++)
      {
        if (j != jBest)
        {
          for (E_Int i = 0; i < nb; i++)
          {
            // chance de muter: inv. proportionelle a la taille
            if (K_NOISE::stdRand(&idum) < a*nbPts[i]+b)
            {
              E_Int targetProc = E_Int(1.*NProc*K_NOISE::stdRand(&idum));
              //printf("block %d move to %d\n", nbPts[i], targetProc);
              popp[i+(j-1)*nb] = targetProc;
            }
          }
        }
      }
    }

    //-------------------------------------------------------------------------
    // Swap
    //-------------------------------------------------------------------------
    for (E_Int j = 1; j <= sizeOfPopulation; j++)
    {
      if (j != jBest)
      {
        E_Int b = E_Int(1.*nb*K_NOISE::stdRand(&idum));
        // Recherche du bloc ayant le nb de pt le + proche
        E_Float nbPts1 = nbPts[b];
        E_Int kCross = 0; E_Float diff = K_CONST::E_MAX_FLOAT;
        E_Float nbPts2 = 0; E_Float diffl;
        for (E_Int k = 0; k < nb; k++)
        {
          if (k != b && popp[k+(j-1)*nb] != popp[b+(j-1)*nb])
          {
            nbPts2 = nbPts[k];
            diffl = K_FUNC::E_abs(nbPts1 - nbPts2);
            if (diffl < diff) {kCross = k; diff = diffl;}
          }
        }
        E_Int temp = popp[b+(j-1)*nb];
        popp[b+(j-1)*nb] = popp[kCross+(j-1)*nb];
        popp[kCross+(j-1)*nb] = temp;
      }
    }

    //-------------------------------------------------------------------------
    // Aggressive Swap
    //-------------------------------------------------------------------------
    if (nbItWithFail > 5)
    {
      for (E_Int j = 1; j <= sizeOfPopulation; j++)
      {
        if (j != jBest)
        {
          for (E_Int b = 0; b < nb; b++)
          {
            // Recherche du bloc ayant le nb de pt le + proche
            E_Float nbPts1 = nbPts[b];
            E_Int kCross = 0; E_Float diff = K_CONST::E_MAX_FLOAT;
            E_Float nbPts2 = 0; E_Float diffl;
            for (E_Int k = 0; k < nb; k++)
            {
              if (k != b && popp[k+(j-1)*nb] != popp[b+(j-1)*nb])
              {
                nbPts2 = nbPts[k];
                diffl = K_FUNC::E_abs(nbPts1 - nbPts2);
                if (diffl < diff) {kCross = k; diff = diffl;}
              }
            }
            if (K_NOISE::stdRand(&idum) >= 1.- diff/K_FUNC::E_max(nbPts1, nbPts2) )
            {
              E_Int temp = popp[b+(j-1)*nb];
              popp[b+(j-1)*nb] = popp[kCross+(j-1)*nb];
              popp[kCross+(j-1)*nb] = temp;
            }
          }
        }
      }
    }

    // On impose les contraintes
    for (E_Int j = 1; j <= sizeOfPopulation; j++ )
      for (E_Int i = 0; i < nb; i++)
        if (setBlocks[i] >= 0) popp[i+(j-1)*nb] = setBlocks[i];

    // On reevalue l'adaptation des survivants:
//#pragma omp parallel for
    for (E_Int j = 1; j <= sizeOfPopulation; j++)
    {
      evalp[j-1] = K_DISTRIBUTOR2::eval(nb, NProc, meanPtsPerProc,
                                        solver, latence,
                                        comSpeed, com, comd, sizeComd,
                                        nbPtsPerProcs, nbPts,
                                        popp+(j-1)*nb);
    }

    // Find best one
    jBest = 1;
    for (E_Int j = 1; j <= sizeOfPopulation; j++)
      jBest = (evalp[j-1] < evalp[jBest-1] ? j : jBest);
    //printf("Adaptation: %f\n", evalp[jBest-1]);

  } // loop

  // Sortie
  //printf("jbest=%d\n", jBest);
  for (E_Int i = 0; i < nb; i++) out[i] = popp[i+(jBest-1)*nb];

  bestAdapt = evalp[jBest-1];
  //printf("Adaptation: %f\n", bestAdapt);

  // external stats
  E_Int empty;
  stats(nbPts, NProc, com, comd, sizeComd, out, empty,
        varMin, varMax, varRMS, volRatio);

}
