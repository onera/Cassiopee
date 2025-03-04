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

#include "kcore.h"
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
  // Acces universel sur BE/ME
  E_Int nc = cEV.getNConnect();
  // Acces universel aux eltTypes
  vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Number of elements and node per element for each connectivity
  vector<E_Int> nelts(nc);
  vector<E_Int> nvertex(nc);
  // Number of face per element and node per face for each connectivity
  vector<E_Int> nfpe(nc);
  vector<E_Int> nnpf(nc);
  vector<vector<E_Int> > f(nc);

  E_Int ierr = 1; // error index, 1 is nominal

  // Boucle sur toutes les connectivites pour remplir face et pre-calculer
  // le nombre de faces connectees a chaque noeud
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cEV.getConnect(ic));
    char* eltTypConn = eltTypes[ic];
    nelts[ic] = cm.getSize();
    nvertex[ic] = cm.getNfld(); // nb de noeuds par element

    if (ierr == 0) continue; // error - skip the rest of the connectivities
    if (nelts[ic] == 0) ierr = 0;
  
    // Tableau de facettes (conforme a CGNS)
    if (K_STRING::cmp(eltTypConn, "BAR") == 0 || 
             K_STRING::cmp(eltTypConn, "BAR*") == 0)
    {
      nfpe[ic] = 2; nnpf[ic] = 1;
      f[ic].reserve(nfpe[ic] * nnpf[ic]);
      f[ic][0 + 0*nfpe[ic]] = 1; 
      f[ic][1 + 0*nfpe[ic]] = 2;
    }
    else if (K_STRING::cmp(eltTypConn, "TRI") == 0 || 
        K_STRING::cmp(eltTypConn, "TRI*") == 0)
    {
      nfpe[ic] = 3; nnpf[ic] = 2;
      f[ic].reserve(nfpe[ic] * nnpf[ic]);
      f[ic][0 + 0*nfpe[ic]] = 1; f[ic][0 + 1*nfpe[ic]] = 2;
      f[ic][1 + 0*nfpe[ic]] = 2; f[ic][1 + 1*nfpe[ic]] = 3;
      f[ic][2 + 0*nfpe[ic]] = 3; f[ic][2 + 1*nfpe[ic]] = 1;
    }
    else if (K_STRING::cmp(eltTypConn, "QUAD") == 0 || 
             K_STRING::cmp(eltTypConn, "QUAD*") == 0)
    {
      nfpe[ic] = 4; nnpf[ic] = 2;
      f[ic].reserve(nfpe[ic] * nnpf[ic]);
      f[ic][0 + 0*nfpe[ic]] = 1; f[ic][0 + 1*nfpe[ic]] = 2;
      f[ic][1 + 0*nfpe[ic]] = 2; f[ic][1 + 1*nfpe[ic]] = 3;
      f[ic][2 + 0*nfpe[ic]] = 3; f[ic][2 + 1*nfpe[ic]] = 4;
      f[ic][3 + 0*nfpe[ic]] = 4; f[ic][3 + 1*nfpe[ic]] = 1;
    }
    else if (K_STRING::cmp(eltTypConn, "TETRA") == 0 || 
             K_STRING::cmp(eltTypConn, "TETRA*") == 0)
    {
      nfpe[ic] = 4; nnpf[ic] = 3;
      f[ic].reserve(nfpe[ic] * nnpf[ic]);
      f[ic][0 + 0*nfpe[ic]] = 1; f[ic][0 + 1*nfpe[ic]] = 2; f[ic][0 + 2*nfpe[ic]] = 3;
      f[ic][1 + 0*nfpe[ic]] = 1; f[ic][1 + 1*nfpe[ic]] = 2; f[ic][1 + 2*nfpe[ic]] = 4;
      f[ic][2 + 0*nfpe[ic]] = 2; f[ic][2 + 1*nfpe[ic]] = 3; f[ic][2 + 2*nfpe[ic]] = 4;
      f[ic][3 + 0*nfpe[ic]] = 3; f[ic][3 + 1*nfpe[ic]] = 1; f[ic][3 + 2*nfpe[ic]] = 4;
    }
    else if (K_STRING::cmp(eltTypConn, "PYRA") == 0 || 
             K_STRING::cmp(eltTypConn, "PYRA*") == 0)
    {
      nfpe[ic] = 5; nnpf[ic] = 4;
      f[ic].reserve(nfpe[ic] * nnpf[ic]);
      f[ic][0 + 0*nfpe[ic]] = 1; f[ic][0 + 1*nfpe[ic]] = 4; f[ic][0 + 2*nfpe[ic]] = 3; f[ic][0 + 3*nfpe[ic]] = 2;
      f[ic][1 + 0*nfpe[ic]] = 1; f[ic][1 + 1*nfpe[ic]] = 2; f[ic][1 + 2*nfpe[ic]] = 5; f[ic][1 + 3*nfpe[ic]] = 1;
      f[ic][2 + 0*nfpe[ic]] = 2; f[ic][2 + 1*nfpe[ic]] = 3; f[ic][2 + 2*nfpe[ic]] = 5; f[ic][2 + 3*nfpe[ic]] = 2;
      f[ic][3 + 0*nfpe[ic]] = 3; f[ic][3 + 1*nfpe[ic]] = 4; f[ic][3 + 2*nfpe[ic]] = 5; f[ic][3 + 3*nfpe[ic]] = 3;
      f[ic][4 + 0*nfpe[ic]] = 4; f[ic][4 + 1*nfpe[ic]] = 1; f[ic][4 + 2*nfpe[ic]] = 5; f[ic][4 + 3*nfpe[ic]] = 4;
    }
    else if (K_STRING::cmp(eltTypConn, "PENTA") == 0 || 
             K_STRING::cmp(eltTypConn, "PENTA*") == 0)
    {
      nfpe[ic] = 5; nnpf[ic] = 4; // TRI degen
      f[ic].reserve(nfpe[ic] * nnpf[ic]);
      f[ic][0 + 0*nfpe[ic]] = 1; f[ic][0 + 1*nfpe[ic]] = 2; f[ic][0 + 2*nfpe[ic]] = 5; f[ic][0 + 3*nfpe[ic]] = 4;
      f[ic][1 + 0*nfpe[ic]] = 2; f[ic][1 + 1*nfpe[ic]] = 3; f[ic][1 + 2*nfpe[ic]] = 6; f[ic][1 + 3*nfpe[ic]] = 5;
      f[ic][2 + 0*nfpe[ic]] = 3; f[ic][2 + 1*nfpe[ic]] = 1; f[ic][2 + 2*nfpe[ic]] = 4; f[ic][2 + 3*nfpe[ic]] = 6;
      f[ic][3 + 0*nfpe[ic]] = 1; f[ic][3 + 1*nfpe[ic]] = 3; f[ic][3 + 2*nfpe[ic]] = 2; f[ic][3 + 3*nfpe[ic]] = 1;
      f[ic][4 + 0*nfpe[ic]] = 4; f[ic][4 + 1*nfpe[ic]] = 5; f[ic][4 + 2*nfpe[ic]] = 6; f[ic][4 + 3*nfpe[ic]] = 4;
    }
    else if (K_STRING::cmp(eltTypConn, "HEXA") == 0 || 
             K_STRING::cmp(eltTypConn, "HEXA*") == 0) 
    {
      nfpe[ic] = 6; nnpf[ic] = 4;
      f[ic].reserve(nfpe[ic] * nnpf[ic]);
      f[ic][0 + 0*nfpe[ic]] = 1; f[ic][0 + 1*nfpe[ic]] = 4; f[ic][0 + 2*nfpe[ic]] = 3; f[ic][0 + 3*nfpe[ic]] = 2;
      f[ic][1 + 0*nfpe[ic]] = 1; f[ic][1 + 1*nfpe[ic]] = 2; f[ic][1 + 2*nfpe[ic]] = 6; f[ic][1 + 3*nfpe[ic]] = 5;
      f[ic][2 + 0*nfpe[ic]] = 2; f[ic][2 + 1*nfpe[ic]] = 3; f[ic][2 + 2*nfpe[ic]] = 7; f[ic][2 + 3*nfpe[ic]] = 6;
      f[ic][3 + 0*nfpe[ic]] = 3; f[ic][3 + 1*nfpe[ic]] = 4; f[ic][3 + 2*nfpe[ic]] = 8; f[ic][3 + 3*nfpe[ic]] = 7;
      f[ic][4 + 0*nfpe[ic]] = 1; f[ic][4 + 1*nfpe[ic]] = 5; f[ic][4 + 2*nfpe[ic]] = 8; f[ic][4 + 3*nfpe[ic]] = 4;
      f[ic][5 + 0*nfpe[ic]] = 5; f[ic][5 + 1*nfpe[ic]] = 6; f[ic][5 + 2*nfpe[ic]] = 7; f[ic][5 + 3*nfpe[ic]] = 8;
    }
    else ierr = 0;
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
  if (ierr == 0) return ierr;

  // Elts voisin d'un noeud
  vector< vector<E_Int> > cVE(nv);
  K_CONNECT::connectEV2VE(cEV, cVE);

  // Boucle sur toutes les connectivites pour remplir cEEN1
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cEV.getConnect(ic));

#pragma omp parallel for default(shared) if (nelts[ic] > __MIN_SIZE_MEAN__)  
    for (E_Int et = 0; et < nelts[ic]; et++)
    {
      E_Int match, ind, ind1, ind2, eltVoisin, node;
      vector<E_Int>& cEEN1 = cEEN[et];
      cEEN1.reserve(nfpe[ic]);

      // pour chaque facette de l'element et
      for (E_Int ff = 0; ff < nfpe[ic]; ff++)
      {
        // Prend le premier noeud de la facette
        ind = cm(et, f[ic][ff + 0*nfpe[ic]])-1;
        
        // Pour les facettes degen
        //if (f[ff + 0*nfpe[ic]] == f[ff + (nnpf[ic]-1)*nfpe[ic]]) nnpfmatch = nnpf[ic]-1;
        //else nnpfmatch = nnpf[ic];

        // Recupere tous les elts ayant ce noeud
        const vector<E_Int>& cVE1 = cVE[ind];

        // Parcourt tout ses voisins pour savoir lequel a la facette
        // en commun
        for (size_t v = 0; v < cVE1.size(); v++)
        {
          eltVoisin = cVE1[v];
          if (eltVoisin != et)
          {
            match = 0;
            for (E_Int k = 0; k < nnpf[ic]; k++)
            {
              ind1 = cm(et, f[ic][ff + k*nfpe[ic]]);
              for (node = 1; node <= nvertex[ic]; node++)
              {
                ind2 = cm(eltVoisin, node);
                if (ind1 == ind2) { match++; break; }
              }
            }
            if (match == nnpf[ic]) { cEEN1.push_back(eltVoisin); goto next; }
          }
        }
        next: ;
      }
    }
  }

  return ierr;
}

// identique mais retourne aussi le no local de la face commune
E_Int K_CONNECT::connectEV2EENbrs(const char* eltType, E_Int nv, 
                                  FldArrayI& cEV,
                                  vector< vector<E_Int> >& cEEN,
                                  vector< vector<E_Int> >& commonFace)
{
  // Acces universel sur BE/ME
  E_Int nc = cEV.getNConnect();
  // Acces universel aux eltTypes
  vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Number of elements and node per element for each connectivity
  vector<E_Int> nelts(nc);
  vector<E_Int> nvertex(nc);
  // Number of face per element and node per face for each connectivity
  vector<E_Int> nfpe(nc);
  vector<E_Int> nnpf(nc);
  vector<vector<E_Int> > f(nc);

  E_Int ierr = 1; // error index, 1 is nominal

  // Boucle sur toutes les connectivites pour remplir face et pre-calculer
  // le nombre de faces connectees a chaque noeud
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cEV.getConnect(ic));
    char* eltTypConn = eltTypes[ic];
    nelts[ic] = cm.getSize();
    nvertex[ic] = cm.getNfld(); // nb de noeuds par element

    if (ierr == 0) continue; // error - skip the rest of the connectivities
    if (nelts[ic] == 0) ierr = 0;
  
    // Tableau de facettes (conforme a CGNS)
    if (K_STRING::cmp(eltTypConn, "TRI") == 0 || 
        K_STRING::cmp(eltTypConn, "TRI*") == 0)
    {
      nfpe[ic] = 3; nnpf[ic] = 2;
      f[ic].reserve(nfpe[ic] * nnpf[ic]);
      f[ic][0 + 0*nfpe[ic]] = 1; f[ic][0 + 1*nfpe[ic]] = 2;
      f[ic][1 + 0*nfpe[ic]] = 2; f[ic][1 + 1*nfpe[ic]] = 3;
      f[ic][2 + 0*nfpe[ic]] = 3; f[ic][2 + 1*nfpe[ic]] = 1;
    }
    else if (K_STRING::cmp(eltTypConn, "QUAD") == 0 || 
             K_STRING::cmp(eltTypConn, "QUAD*") == 0)
    {
      nfpe[ic] = 4; nnpf[ic] = 2;
      f[ic].reserve(nfpe[ic] * nnpf[ic]);
      f[ic][0 + 0*nfpe[ic]] = 1; f[ic][0 + 1*nfpe[ic]] = 2;
      f[ic][1 + 0*nfpe[ic]] = 2; f[ic][1 + 1*nfpe[ic]] = 3;
      f[ic][2 + 0*nfpe[ic]] = 3; f[ic][2 + 1*nfpe[ic]] = 4;
      f[ic][3 + 0*nfpe[ic]] = 4; f[ic][3 + 1*nfpe[ic]] = 1;
    }
    else if (K_STRING::cmp(eltTypConn, "TETRA") == 0 || 
             K_STRING::cmp(eltTypConn, "TETRA*") == 0)
    {
      nfpe[ic] = 4; nnpf[ic] = 3;
      f[ic].reserve(nfpe[ic] * nnpf[ic]);
      f[ic][0 + 0*nfpe[ic]] = 1; f[ic][0 + 1*nfpe[ic]] = 2; f[ic][0 + 2*nfpe[ic]] = 3;
      f[ic][1 + 0*nfpe[ic]] = 1; f[ic][1 + 1*nfpe[ic]] = 2; f[ic][1 + 2*nfpe[ic]] = 4;
      f[ic][2 + 0*nfpe[ic]] = 2; f[ic][2 + 1*nfpe[ic]] = 3; f[ic][2 + 2*nfpe[ic]] = 4;
      f[ic][3 + 0*nfpe[ic]] = 3; f[ic][3 + 1*nfpe[ic]] = 1; f[ic][3 + 2*nfpe[ic]] = 4;
    }
    else if (K_STRING::cmp(eltTypConn, "HEXA") == 0 || 
             K_STRING::cmp(eltTypConn, "HEXA*") == 0) 
    {
      nfpe[ic] = 6; nnpf[ic] = 4;
      f[ic].reserve(nfpe[ic] * nnpf[ic]);
      f[ic][0 + 0*nfpe[ic]] = 1; f[ic][0 + 1*nfpe[ic]] = 4; f[ic][0 + 2*nfpe[ic]] = 3; f[ic][0 + 3*nfpe[ic]] = 2;
      f[ic][1 + 0*nfpe[ic]] = 1; f[ic][1 + 1*nfpe[ic]] = 2; f[ic][1 + 2*nfpe[ic]] = 6; f[ic][1 + 3*nfpe[ic]] = 5;
      f[ic][2 + 0*nfpe[ic]] = 2; f[ic][2 + 1*nfpe[ic]] = 3; f[ic][2 + 2*nfpe[ic]] = 7; f[ic][2 + 3*nfpe[ic]] = 6;
      f[ic][3 + 0*nfpe[ic]] = 3; f[ic][3 + 1*nfpe[ic]] = 4; f[ic][3 + 2*nfpe[ic]] = 8; f[ic][3 + 3*nfpe[ic]] = 7;
      f[ic][4 + 0*nfpe[ic]] = 1; f[ic][4 + 1*nfpe[ic]] = 5; f[ic][4 + 2*nfpe[ic]] = 8; f[ic][4 + 3*nfpe[ic]] = 4;
      f[ic][5 + 0*nfpe[ic]] = 5; f[ic][5 + 1*nfpe[ic]] = 6; f[ic][5 + 2*nfpe[ic]] = 7; f[ic][5 + 3*nfpe[ic]] = 8;
    }
    else if (K_STRING::cmp(eltTypConn, "BAR") == 0 || 
             K_STRING::cmp(eltTypConn, "BAR*") == 0)
    {
      nfpe[ic] = 2; nnpf[ic] = 1;
      f[ic].reserve(nfpe[ic] * nnpf[ic]);
      f[ic][0 + 0*nfpe[ic]] = 1; 
      f[ic][1 + 0*nfpe[ic]] = 2;
    }
    else if (K_STRING::cmp(eltTypConn, "PYRA") == 0 || 
             K_STRING::cmp(eltTypConn, "PYRA*") == 0)
    {
      nfpe[ic] = 5; nnpf[ic] = 4;
      f[ic].reserve(nfpe[ic] * nnpf[ic]);
      f[ic][0 + 0*nfpe[ic]] = 1; f[ic][0 + 1*nfpe[ic]] = 4; f[ic][0 + 2*nfpe[ic]] = 3; f[ic][0 + 3*nfpe[ic]] = 2;
      f[ic][1 + 0*nfpe[ic]] = 1; f[ic][1 + 1*nfpe[ic]] = 2; f[ic][1 + 2*nfpe[ic]] = 5; f[ic][1 + 3*nfpe[ic]] = 1;
      f[ic][2 + 0*nfpe[ic]] = 2; f[ic][2 + 1*nfpe[ic]] = 3; f[ic][2 + 2*nfpe[ic]] = 5; f[ic][2 + 3*nfpe[ic]] = 2;
      f[ic][3 + 0*nfpe[ic]] = 3; f[ic][3 + 1*nfpe[ic]] = 4; f[ic][3 + 2*nfpe[ic]] = 5; f[ic][3 + 3*nfpe[ic]] = 3;
      f[ic][4 + 0*nfpe[ic]] = 4; f[ic][4 + 1*nfpe[ic]] = 1; f[ic][4 + 2*nfpe[ic]] = 5; f[ic][4 + 3*nfpe[ic]] = 4;
    }
    else if (K_STRING::cmp(eltTypConn, "PENTA") == 0 || 
             K_STRING::cmp(eltTypConn, "PENTA*") == 0)
    {
      nfpe[ic] = 5; nnpf[ic] = 4;
      f[ic].reserve(nfpe[ic] * nnpf[ic]);
      f[ic][0 + 0*nfpe[ic]] = 1; f[ic][0 + 1*nfpe[ic]] = 2; f[ic][0 + 2*nfpe[ic]] = 5; f[ic][0 + 3*nfpe[ic]] = 4;
      f[ic][1 + 0*nfpe[ic]] = 2; f[ic][1 + 1*nfpe[ic]] = 3; f[ic][1 + 2*nfpe[ic]] = 6; f[ic][1 + 3*nfpe[ic]] = 5;
      f[ic][2 + 0*nfpe[ic]] = 3; f[ic][2 + 1*nfpe[ic]] = 1; f[ic][2 + 2*nfpe[ic]] = 4; f[ic][2 + 3*nfpe[ic]] = 6;
      f[ic][3 + 0*nfpe[ic]] = 1; f[ic][3 + 1*nfpe[ic]] = 3; f[ic][3 + 2*nfpe[ic]] = 2; f[ic][3 + 3*nfpe[ic]] = 1;
      f[ic][4 + 0*nfpe[ic]] = 4; f[ic][4 + 1*nfpe[ic]] = 5; f[ic][4 + 2*nfpe[ic]] = 6; f[ic][4 + 3*nfpe[ic]] = 4;
    }
    else ierr = 0;
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
  if (ierr == 0) return ierr;

  // Elts voisin d'un noeud
  vector< vector<E_Int> > cVE(nv);
  K_CONNECT::connectEV2VE(cEV, cVE);

  // Boucle sur toutes les connectivites pour remplir cEEN1
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cEV.getConnect(ic));

#pragma omp parallel for default(shared) if (nelts[ic] > __MIN_SIZE_MEAN__)  
    for (E_Int et = 0; et < nelts[ic]; et++)
    {
      E_Int match, ind, ind1, ind2, eltVoisin, node;
      vector<E_Int>& cEEN1 = cEEN[et];
      cEEN1.reserve(nfpe[ic]);
      vector<E_Int>& commonFace1 = commonFace[et];
      commonFace1.reserve(nfpe[ic]);

      // pour chaque facette de l'element et
      for (E_Int ff = 0; ff < nfpe[ic]; ff++)
      {
        // Prend le premier noeud de la facette
        ind = cm(et, f[ic][ff + 0*nfpe[ic]])-1;
        
        // Pour les facettes degen
        //if (f[ff + 0*nfpe[ic]] == f[ff + (nnpf[ic]-1)*nfpe[ic]]) nnpfmatch = nnpf[ic]-1;
        //else nnpfmatch = nnpf[ic];

        // Recupere tous les elts ayant ce noeud
        const vector<E_Int>& cVE1 = cVE[ind];

        // Parcourt tout ses voisins pour savoir lequel a la facette
        // en commun
        for (size_t v = 0; v < cVE1.size(); v++)
        {
          eltVoisin = cVE1[v];
          if (eltVoisin != et)
          {
            match = 0;
            for (E_Int k = 0; k < nnpf[ic]; k++)
            {
              ind1 = cm(et, f[ic][ff + k*nfpe[ic]]);
              for (node = 1; node <= nvertex[ic]; node++)
              {
                ind2 = cm(eltVoisin, node);
                if (ind1 == ind2) { match++; break; }
              }
            }
            if (match == nnpf[ic]) 
            { 
              cEEN1.push_back(eltVoisin);
              commonFace1.push_back(ff); // common face
              goto next; 
            }
          }
        }
        next: ;
      }
    }
  }

  return ierr;
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
