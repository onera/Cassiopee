/*    
    Copyright 2013-2019 Onera.

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

#include "Connect/connect.h"
#include "parallel.h"
#include <algorithm>
#include <string.h>
#include "String/kstring.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/*
  Change a Elts-Vertex connectivity to a Element-Element neighbours 
  connectivity. 
  IN: cEV: Elts-Vertex connectivity. For each elt, give vertices index.
  IN: nv: nombre de noeuds dans le maillage
  IN: eltType: type d'element (TRI, QUAD,...)
  OUT: cEEN: Elts-Element neighbours connectivity. For each element, 
  give the element index of its neighbours.
  cEEN doit deja etre alloue au nombre d'elements.
  Retourne 0 (echec), 1 (succes)
  ! Marche pour les maillages non structures a elements basiques
  ! Marche uniquement pour les maillage conformes
  ! Warning: par definition, un voisin a une facette commune avec un elt
  Algo: pour chaque facette d'un element, on cherche a la matcher a un
  voisin. 
*/
//=============================================================================
E_Int K_CONNECT::connectEV2EENbrs(const char* eltType, E_Int nv, 
                                  FldArrayI& cEV,
                                  vector< vector<E_Int> >& cEEN)
{
  E_Int nelts = cEV.getSize();
  E_Int nvertex = cEV.getNfld(); // nb de noeuds par element

  if (nelts == 0) return 0;
  
  // Tableau de facettes (conforme a CGNS)
  E_Int* f; E_Int nfpe; E_Int nnpf;
  if (K_STRING::cmp(eltType, "TRI") == 0 || 
      K_STRING::cmp(eltType, "TRI*") == 0)
  {
    nfpe = 3; nnpf = 2;
    f = new E_Int[nfpe * nnpf]; 
    f[0 + 0*nfpe] = 1; f[0 + 1*nfpe] = 2;
    f[1 + 0*nfpe] = 2; f[1 + 1*nfpe] = 3;
    f[2 + 0*nfpe] = 3; f[2 + 1*nfpe] = 1;
  }
  else if (K_STRING::cmp(eltType, "QUAD") == 0 || 
           K_STRING::cmp(eltType, "QUAD*") == 0)
  {
    nfpe = 4; nnpf = 2;
    f = new E_Int[nfpe * nnpf]; 
    f[0 + 0*nfpe] = 1; f[0 + 1*nfpe] = 2;
    f[1 + 0*nfpe] = 2; f[1 + 1*nfpe] = 3;
    f[2 + 0*nfpe] = 3; f[2 + 1*nfpe] = 4;
    f[3 + 0*nfpe] = 4; f[3 + 1*nfpe] = 1;
  }
  else if (K_STRING::cmp(eltType, "TETRA") == 0 || 
           K_STRING::cmp(eltType, "TETRA*") == 0)
  {
    nfpe = 4; nnpf = 3;
    f = new E_Int[nfpe * nnpf];
    f[0 + 0*nfpe] = 1; f[0 + 1*nfpe] = 2; f[0 + 2*nfpe] = 3;
    f[1 + 0*nfpe] = 1; f[1 + 1*nfpe] = 2; f[1 + 2*nfpe] = 4;
    f[2 + 0*nfpe] = 2; f[2 + 1*nfpe] = 3; f[2 + 2*nfpe] = 4;
    f[3 + 0*nfpe] = 3; f[3 + 1*nfpe] = 1; f[3 + 2*nfpe] = 4;
  }
  else if (K_STRING::cmp(eltType, "HEXA") == 0 || 
           K_STRING::cmp(eltType, "HEXA*") == 0) 
  {
    nfpe = 6; nnpf = 4;
    f = new E_Int[nfpe * nnpf];
    f[0 + 0*nfpe] = 1; f[0 + 1*nfpe] = 4; f[0 + 2*nfpe] = 3; f[0 + 3*nfpe] = 2;
    f[1 + 0*nfpe] = 1; f[1 + 1*nfpe] = 2; f[1 + 2*nfpe] = 6; f[1 + 3*nfpe] = 5;
    f[2 + 0*nfpe] = 2; f[2 + 1*nfpe] = 3; f[2 + 2*nfpe] = 7; f[2 + 3*nfpe] = 6;
    f[3 + 0*nfpe] = 3; f[3 + 1*nfpe] = 4; f[3 + 2*nfpe] = 8; f[3 + 3*nfpe] = 7;
    f[4 + 0*nfpe] = 1; f[4 + 1*nfpe] = 5; f[4 + 2*nfpe] = 8; f[4 + 3*nfpe] = 4;
    f[5 + 0*nfpe] = 5; f[5 + 1*nfpe] = 6; f[5 + 2*nfpe] = 7; f[5 + 3*nfpe] = 8;
  }
  else if (K_STRING::cmp(eltType, "BAR") == 0 || 
           K_STRING::cmp(eltType, "BAR*") == 0)
  {
    nfpe = 2; nnpf = 1;
    f = new E_Int[nfpe * nnpf]; 
    f[0 + 0*nfpe] = 1; 
    f[1 + 0*nfpe] = 2;
  }
  else if (K_STRING::cmp(eltType, "PYRA") == 0 || 
           K_STRING::cmp(eltType, "PYRA*") == 0)
  {
    nfpe = 5; nnpf = 4;
    f = new E_Int[nfpe * nnpf];
    f[0 + 0*nfpe] = 1; f[0 + 1*nfpe] = 4; f[0 + 2*nfpe] = 3; f[0 + 3*nfpe] = 2;
    f[1 + 0*nfpe] = 1; f[1 + 1*nfpe] = 2; f[1 + 2*nfpe] = 5; f[1 + 3*nfpe] = 1;
    f[2 + 0*nfpe] = 2; f[2 + 1*nfpe] = 3; f[2 + 2*nfpe] = 5; f[2 + 3*nfpe] = 2;
    f[3 + 0*nfpe] = 3; f[3 + 1*nfpe] = 4; f[3 + 2*nfpe] = 5; f[3 + 3*nfpe] = 3;
    f[4 + 0*nfpe] = 4; f[4 + 1*nfpe] = 1; f[4 + 2*nfpe] = 5; f[4 + 3*nfpe] = 4;
  }
  else if (K_STRING::cmp(eltType, "PENTA") == 0 || 
           K_STRING::cmp(eltType, "PENTA*") == 0)
  {
    nfpe = 5; nnpf = 4;
    f = new E_Int[nfpe * nnpf];
    f[0 + 0*nfpe] = 1; f[0 + 1*nfpe] = 2; f[0 + 2*nfpe] = 5; f[0 + 3*nfpe] = 4;
    f[1 + 0*nfpe] = 2; f[1 + 1*nfpe] = 3; f[1 + 2*nfpe] = 6; f[1 + 3*nfpe] = 5;
    f[2 + 0*nfpe] = 3; f[2 + 1*nfpe] = 1; f[2 + 2*nfpe] = 4; f[2 + 3*nfpe] = 6;
    f[3 + 0*nfpe] = 1; f[3 + 1*nfpe] = 3; f[3 + 2*nfpe] = 2; f[3 + 3*nfpe] = 1;
    f[4 + 0*nfpe] = 4; f[4 + 1*nfpe] = 5; f[4 + 2*nfpe] = 6; f[4 + 3*nfpe] = 4;
  }
  else return 0;

  // Elts voisin d'un noeud
  vector< vector<E_Int> > cVE(nv);
  K_CONNECT::connectEV2VE(cEV, cVE);

#pragma omp parallel for default(shared) if (nelts > __MIN_SIZE_MEAN__)  
  for (E_Int et = 0; et < nelts; et++)
  {
    E_Int size, match, ind, ind1, ind2, eltVoisin, node;
    vector<E_Int>& cEEN1 = cEEN[et];
    cEEN1.reserve(nfpe);

    // pour chaque facette de l'element et
    for (E_Int ff = 0; ff < nfpe; ff++)
    {
      // Prend le premier noeud de la facette
      ind = cEV(et, f[ff + 0*nfpe])-1;
      
      // Recupere tous les elts ayant ce noeud
      vector<E_Int>& cVE1 = cVE[ind];
      size = cVE1.size();

      // Parcourt tout ses voisins pour savoir lequel a la facette
      // en commun
      for (E_Int v = 0; v < size; v++)
      {
        eltVoisin = cVE1[v];
        if (eltVoisin != et)
        {
          match = 0;
          for (E_Int k = 0; k < nnpf; k++)
          {
            ind1 = cEV(et, f[ff + k*nfpe]);
            for (node = 1; node <= nvertex; node++)
            {
              ind2 = cEV(eltVoisin, node);
              if (ind1 == ind2) { match++; break; }
            }
          }
          if (match == nnpf) { cEEN1.push_back(eltVoisin); goto next; }
        }
      }
      next: ;
    }
  }
  delete [] f;
  return 1;
}

//=============================================================================
/* Recherche de l'element voisin (PENTA ou HEXA) contenant aussi la facette 
   QUAD ABCD dans cEEN */
//=============================================================================
E_Int K_CONNECT::getNbrForQuadFace(
  E_Int indA, E_Int indB, E_Int indC, E_Int indD, 
  FldArrayI& cn, vector<E_Int>& cEEN)
{
  E_Int nvoisins = cEEN.size();
  E_Int ind;
  E_Int nfld = cn.getNfld();
  for (E_Int noet = 0; noet < nvoisins; noet++)
  {
    E_Int et = cEEN[noet];
    E_Int c = 0;
    for (E_Int i = 1; i <= nfld; i++)
    {
      ind = cn(et,i)-1;
      if (ind == indA || ind == indB || ind == indC || ind == indD) c++; 
    }
    if (c == 4) return et;
  }
  return -1;
}
