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

# include "Metis/metis.h"

using namespace std;
using namespace K_FLD;


//========================================================================
// Find empty proc
//========================================================================
E_Int findEmptyProc(E_Int NProc, FldArrayF& nbNodePerProc)
{
  for (E_Int i = 0; i < NProc; i++)
  {
    //printf("Proc %d: %g pts\n", i, nbNodePerProc[i]);
    if (K_FUNC::E_abs(nbNodePerProc[i]) < 1.e-6) return i;
  }
  return -1;
}

//=======================================================================
// Fill empty proc with another block
//=======================================================================
E_Int fillEmpty(E_Int iEmpty, E_Int NProc, E_Int nb, std::vector<E_Int>& out, 
                FldArrayF& nbNodePerProc, FldArrayI& nbBlockPerProc,
                std::vector<E_Float>& nbPts)
{
    // Find fullest proc
    E_Float maxPts = 0; E_Int iMax = -1;
    for (E_Int i = 0; i < NProc; i++)
    {
        if (nbNodePerProc[i] > maxPts && nbBlockPerProc[i] >= 2)
        {
            maxPts = nbNodePerProc[i];
            iMax = i;
        }
    }
    if (iMax == -1) return -1; // fail
    // find first block on imax
    for (E_Int i = 0; i < nb; i++)
    {
        if (out[i] == iMax)
        {
            out[i] = iEmpty;
            nbNodePerProc[iMax] -= nbPts[i];
            nbNodePerProc[iEmpty] += nbPts[i];
            nbBlockPerProc[iMax] -= 1;
            nbBlockPerProc[iEmpty] += 1;
            return 1; // success
        }
    }
    return -1; // fail
}

//=============================================================================
// IN: nbPts: pour chaque bloc, son nombre de pts
// IN: setBlocks: pour chaque bloc, -1 si ce bloc est a equilibrer,
// si c'est une autre valeur, on considere ce bloc comme deja
// place sur ce proc
// IN: NProc: nombre de processeurs
// OUT: out: pour chaque bloc, son no de proc
//=============================================================================
void K_DISTRIBUTOR2::graph(
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

  // Enforce com graph symetry (necessary for metis but not necessarily ensured by IBM)
  if (com != NULL)
  {
    for (E_Int i = 0; i < nb; i++)
      for (E_Int j = 0; j < nb; j++)
      {
        if (com[i+j*nb] != com[j+i*nb])
        {
          //printf("No Symetry - forced: %d %d (%d %d)\n", i,j, com[i+j*nb], com[j+i*nb]);
          //com[i+j*nb] = com[j+i*nb];
          if (com[j+i*nb] > 0) com[i+j*nb] = com[j+i*nb];
          else com[j+i*nb] = com[i+j*nb];
        }
      }
  }
  // Force symetrie (matrice inferieure)
  if (comd != NULL)
  {
    E_Int b1,b2,bb1,bb2;
    bool exist;
    for (E_Int i = 0; i < sizeComd/2; i++)
    {
      b2 = comd[2*i]/nb;
      b1 = comd[2*i]-b2*nb;
      if (b1 > b2) // retourne que si il n'existe pas deja dans la partie inferieure de la matrice
      {  
        exist = false;  
        for (E_Int n = 0; n < sizeComd/2; n++)
        {
          bb2 = comd[2*n]/nb;
          bb1 = comd[2*n]-bb2*nb;
          if (n != i && bb1 == b2 && bb2 == b1) { exist = true; break; }
        }
        if (exist == false) comd[2*i] = b2 + b1*nb; // inverse 
      }
    }
  }

  /*
  for (E_Int i = 0; i < nb; i++)
  {
    for (E_Int j = 0; j < nb; j++)
    {
      printf("%d ",com[i+j*nb]);
    }
    printf("\n");
  }
  */
  // Graph : les vertex du graph sont les blocs, 
  // le poids vertex est le nbre de pts
  // les edges du graph sont les coms
  // le poids de edges est le volume de com
  
  // taille des adj
  E_Int size = 0;
  if (com != NULL)
  {
    for (E_Int i = 0; i < nb; i++)
      for (E_Int j = 0; j < nb; j++)
      {
        if (com[i+j*nb] > 0 && i != j) size++;
      }
  }
  if (comd != NULL)
  {
    //size = sizeComd/2;
    E_Int b1,b2;
    for (E_Int i = 0; i < nb; i++)
    {
      for (E_Int n = 0; n < sizeComd/2; n++)
      {
        b2 = comd[2*n]/nb;
        b1 = comd[2*n]-b2*nb;
        if (b1 < b2 && b1 == i) size++;
        else if (b1 < b2 && b2 == i) size++;
      }
    }
  }
  
  //printf("size of adj %d\n", size);
  
  idx_t* adj = new idx_t [size];
  idx_t* xadj = new idx_t [nb+1];
  idx_t* vweight = new idx_t [nb];
  idx_t* adjweight = new idx_t [size];
  
  // remplissage des poids des blocs
  //printf("weight=");
  for (E_Int i = 0; i < nb; i++)
  {
    vweight[i] = nbPts[i];
    //printf("%d ", vweight[i]);
  }
  //printf("\n");

  // com relative weight
  E_Float rel = 1.;

  // remplissage adj + xadj a partir de com
  if (com != NULL)
  {
    size = 0;
    for (E_Int i = 0; i < nb; i++)
    {
      xadj[i] = size;
      for (E_Int j = 0; j < nb; j++)
      {
        if (com[i+j*nb] > 0 && i != j)
        {
          adj[size] = j;
          //printf("%d %d\n",com[i+j*nb],com[j+i*nb]);
          //if (i < j) adjweight[size] = rel*com[i+j*nb];
          //else adjweight[size] = rel*com[j+i*nb]; // force symetrie
          adjweight[size] = rel*com[i+j*nb]+1;
          size++;
        }
      }
      xadj[nb] = size;
    }
  }
  
  // remplissage adj + xadj a partir de comd (symetrisation ici)
  if (comd != NULL)
  {
    size = 0;
    E_Int b1,b2;
    for (E_Int i = 0; i < nb; i++)
    {
      xadj[i] = size;

      for (E_Int n = 0; n < sizeComd/2; n++)
      {
        b2 = comd[2*n]/nb;
        b1 = comd[2*n]-b2*nb;
        if (b1 < b2 && b1 == i)
        {
          adj[size] = b2;
          adjweight[size] = rel*comd[2*n+1]+1;
          size++;
        }
        else if (b1 < b2 && b2 == i)
        {
          adj[size] = b1;
          adjweight[size] = rel*comd[2*n+1]+1;
          size++;
        }
      }
      xadj[nb] = size;
    }
  }

  //printf("xadj = ");
  //for (E_Int i = 0; i < nb; i++) printf("%d ", xadj[i]); 
  //printf("\n");
  //printf("size of adj %d\n", size);
  //for (E_Int i = 0; i < size; i++) printf("%d adj=%d adjw=%d\n", i, adj[i], adjweight[i]);

  E_Int objval = 0; // retour
  E_Int ncon = 1;
  
  //idx_t options[METIS_NOPTIONS];
  //METIS_SetDefaultOptions(options);
  //options[METIS_OPTION_CONTIG] = 1; // force contiguite
  //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
  //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
  //options[METIS_OPTION_MINCONN] = 1; // force min connectivite externe
  //options[METIS_OPTION_UFACTOR] = 1.; // imbalance tol

  idx_t* parts = new idx_t [nb];
  METIS_PartGraphKway(&nb, &ncon, xadj, adj, vweight, NULL, adjweight, 
                      &NProc, NULL, NULL, NULL, &objval, parts);
  //METIS_PartGraphKway(&nb, &ncon, xadj, adj, vweight, NULL, NULL, 
  //                     &NProc, NULL, NULL, options, &objval, parts);
  //METIS_PartGraphRecursive(&nb, &ncon, xadj, adj, vweight, NULL, adjweight,
  //                         &NProc, NULL, NULL, options, &objval, parts);
  
  delete [] adj; delete [] xadj;
  delete [] vweight; delete [] adjweight;

  // Sortie
  for (E_Int i = 0; i < nb; i++) out[i] = parts[i];
  delete [] parts;

  // Calcul du nombre de pts par processeurs
  nbNodePerProc.setAllValuesAtNull();
  for (E_Int i = 0; i < nb; i++)
  {
    proc = out[i];
    nbNodePerProc[proc] += nbPts[i];
  }
  //for (E_Int i = 0; i < NProc; i++) printf("%d %g\n",i,nbNodePerProc[i]);

  FldArrayI nbBlockPerProc(NProc);
  nbBlockPerProc.setAllValuesAtNull();
  for (E_Int i = 0; i < nb; i++)
  {
    proc = out[i];
    nbBlockPerProc[proc] += 1;
  }
  //printf("Info: Nb de pts moyen par proc: %d\n", int(meanPtsPerProc));
 
  // Si processeur empty, prendre un bloc sur le plus chargé iterativement
  E_Int iEmpty;
  iEmpty = findEmptyProc(NProc, nbNodePerProc);
  while (iEmpty >= 0)
  {
    iEmpty = fillEmpty(iEmpty, NProc, nb, out, nbNodePerProc, nbBlockPerProc, nbPts);
    if (iEmpty == -1) break;
    iEmpty = findEmptyProc(NProc, nbNodePerProc);
  }

  bestAdapt = 1.;
  //printf("Adaptation: %f\n", bestAdapt);

  // external stats
  E_Int empty;
  stats(nbPts, NProc, com, comd, sizeComd, out, empty, 
        varMin, varMax, varRMS, volRatio);

}
