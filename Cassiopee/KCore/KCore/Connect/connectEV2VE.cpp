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
#include <iostream>
#include "Connect/connect.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/*
  Change a Elts-Vertex connectivity to a Vertex-Elts connectivity
  IN: cEV: Elts-Vertex connectivity. For each elt, give vertices index.
  OUT: cVE: Vertex-Elts connectivity. For each vertex, give elts index.
  Les indices de vertex commencent a 1
  Les indices d'elements commencent a 0
  Attention: cVE doit deja etre alloue au nombre de noeuds.
*/
//=============================================================================
void K_CONNECT::connectEV2VE(FldArrayI& cEV, vector< vector<E_Int> >& cVE)
{
  // Acces universel sur BE/ME
  E_Int nc = cEV.getNConnect();

  E_Int nv = cVE.size(); // Nombre de points du maillage
  // Calcul du nombre d'elements attaches a chaque noeud
  FldArrayI nve(nv); nve.setAllValuesAtNull();
  E_Int* nvep = nve.begin();
  E_Int vertex;
  E_Int el_offset = 0; // decalage

  // Boucle sur toutes les connectivites pour pre-calculer la valence de
  // chaque noeud
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cEV.getConnect(ic));
    // Nbre elements de cette connectivite
    E_Int ne = cm.getSize();
    // Nombre de points par elements de cette connectivite
    E_Int ns = cm.getNfld();

    for (E_Int j = 1; j <= ns; j++)
    {
      for (E_Int i = 0; i < ne; i++)
      {
        vertex = cm(i, j);
        nvep[vertex-1]++;
      }
    }
  }

#pragma omp parallel for default (shared)
  for (E_Int i = 0; i < nv; i++) cVE[i].reserve(nvep[i]);

  // Boucle sur toutes les connectivites pour remplir cVE
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cEV.getConnect(ic));
    E_Int ne = cm.getSize();
    E_Int ns = cm.getNfld();
    
    for (E_Int j = 1; j <= ns; j++)
    {
      for (E_Int i = 0; i < ne; i++)
      {
        vertex = cm(i, j);
        cVE[vertex-1].push_back(el_offset + i);
      }
    }

    el_offset += ne; // increment element offset
  }
}

std::pair<std::vector<E_Int>,std::vector<E_Int> >
K_CONNECT::connectEV2VE(K_FLD::FldArrayI& cEV)
{
  E_Int ne       = cEV.getSize();// Nbre elements contenus dans le maillage
  E_Int ns       = cEV.getNfld();// Nombre de points par elements
  E_Int stride   = cEV.getStride();

  E_Int nv = 0;
  // On prepare un tableau de pointeur pour chaque composante de la connectivité
  std::vector<E_Int*> pt_nevs(ns);
  for (E_Int j_vert = 0; j_vert < ns; ++j_vert)
    pt_nevs[j_vert] = cEV.begin(j_vert+1);

  // On compte le nombre d'elements contenant chaque sommet
  // On en profite pour compter le nombre de sommets du maillage
  // a l'aide de sa numerotation.
  std::vector<E_Int> counter(ns*ne, 0);
  E_Int* pt_cter = counter.data();
  E_Int vertex;
  E_Int* pt_nev;
  for (E_Int j_vert = 0; j_vert < ns; ++j_vert)
  {
    pt_nev = pt_nevs[j_vert];
    for (E_Int i_elt = 0; i_elt < ne; ++i_elt)
    {
      vertex = pt_nev[i_elt*stride]-1;
      nv = (nv <= vertex ? vertex+1 : nv);
      pt_cter[vertex] += 1;
    }
  }
  // On reserve le nombre de sommets + 1 pour le premier tableau
  std::vector<E_Int> beg_vert2elts(nv+1);
  // 1er element commence au debut du tableau vert2elts :
  beg_vert2elts[0] = 0;
  // Les suivants, commenceront n valeurs plus loin du debut du sommet precedent, où n est
  // le nombre d'elements contenant le sommet precedent.
  for (E_Int iv = 0; iv < nv; ++iv)
    beg_vert2elts[iv+1] = beg_vert2elts[iv] + pt_cter[iv];
  E_Int* pt_beg = beg_vert2elts.data();
  // Le tableau v2e est de taille la valeur du dernier element de beg_vert2elts :
  std::vector<E_Int> vert2elts(beg_vert2elts[nv]);
  E_Int* pt_v2e = vert2elts.data();
  // On remet a zero le compteur.
  std::fill(counter.begin(), counter.end(), 0);
  // Puis on remplit le tableau vert2elts :
  // Pour chaque element
  for (E_Int i_elt = 0; i_elt < ne; ++i_elt)
  {
    // On regarde les sommets qu'il contient
    for (E_Int i_vert = 0; i_vert < ns; ++i_vert)
    {
      vertex = pt_nevs[i_vert][i_elt*stride]-1;
      // Et on rajoute cet element aux elements contenant ce sommet
      pt_v2e[pt_beg[vertex] + pt_cter[vertex]] = i_elt;
      pt_cter[vertex] += 1;
    }
  }
  return {std::move(beg_vert2elts), std::move(vert2elts)};
}
