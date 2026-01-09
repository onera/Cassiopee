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
#include "Connect/connect.h"
#include <stdlib.h>
#include <vector>
#include <algorithm>    // std::reverse
#include <stdio.h>
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Reordonne la numerotation d'un array NGON
   Sa connectivite doit etre propre (cleanConnectivity).
   Retourne 1 si OK. */
// ============================================================================
E_Int K_CONNECT::reorderNGON(FldArrayF& f, FldArrayI& cn, E_Int dir)
{
  // Acces non universel sur le ptrs
  E_Int* ngon = cn.getNGon();
  E_Int* nface = cn.getNFace();
  E_Int* indPG = cn.getIndPG();
  E_Int* indPH = cn.getIndPH();
  // Acces universel nbre d'elements
  E_Int nelts = cn.getNElts();
  // Connectivite ets/ets voisins
  vector< vector<E_Int> > cEEN(nelts);
  FldArrayI cFE; connectNG2FE(cn, cFE);
  connectFE2EENbrs(cFE, cEEN);

  E_Int nev = 0; // nbre d'elements deja visites
  char* isVisited = (char*)calloc(nelts, sizeof(char)); // elt deja visite?
  E_Int* mustBeVisited = (E_Int*)malloc(nelts * sizeof(E_Int));
  E_Int mbv = 0;
  E_Int p, size, ie, iv, ienei, nbFaces, nbFacesNei, dummy;
  vector<E_Int> indices;

  while (nev < nelts)
  {
    // Recherche du premier elt pas encore visite
    for (p = 0; (isVisited[p] != 0); p++);
    // appartient a un nouveau composant connexe
    // alors on recupere les vertex dans l'ordre rotatif et on ordonne l'elt de depart
    K_CONNECT::getVertexIndices(cn, ngon, nface, indPG, indPH, p, indices);

/*    printf("vertices initiaux\n");
    for (E_Int i = 0; i < indices.size(); i++) printf("%d ", indices[i]);
    printf("\n");
*/
    if (dir == -1) reverse(indices.begin(),indices.end());
    orderNGONElement(p, indices, ngon, nface, indPG, indPH, cn);

    indices.clear();

    mustBeVisited[mbv] = p;
    mbv++; nev++;
    isVisited[p] = 1;

    while (mbv > 0)
    {
      mbv--;
      ie = mustBeVisited[mbv];
      size = cEEN[ie].size();
      E_Int* elt = cn.getElt(ie, nbFaces, nface, indPH);
      for (iv = 0; iv < size; iv++)
      {
        ienei = cEEN[ie][iv]; // index de l'elt voisin: ienei
        if (isVisited[ienei] == 0)
        {
          // Acces universel element voisin ienei
          E_Int* eltn = cn.getElt(ienei, nbFacesNei, nface, indPH);
          indices.clear();

          // Trouve 1 face commune entre les deux elements
          for (E_Int nof1 = 0; nof1 < nbFaces; nof1++)
            for (E_Int nof2 = 0; nof2 < nbFacesNei; nof2++)
            {
              if (elt[nof1] == eltn[nof2])
              {
                // Une fois la face commune trouvee, on reconstruit le cycle des sommets
                // en accord avec la circulation de la cellule voisine mbv

                /* Notations :

                                          indv1 (resp indv2)=>first
                __________________________________X________________________________
                |                                 |                               |
                |           <<--------            |                               |
                |                                 |face commune                   |
                |                                 |                               |
                |                           ^     |                               |
                |    |                      ^     |                               |
                |    |                      |     |                               |
                |    |                      |     |                               |
                |    |         mbv          |     |            ienei              |
                |    v  dont la circulation |     |                               |
                |    v      est connue      |     |                               |
                |                                 |                               |
                |                                 |                               |
                |           -------->>            |                               |
                |              pface              |          nextface             |
                |_________________________________|_______________________________|
                X                                 X
                pn2 (resp pn1)              indv2 (resp indv1)=>next
                                           et pn1 (resp pn2)

                */
                // Acces universel aux faces
                E_Int* faceNei = cn.getFace(eltn[nof2]-1, dummy, ngon, indPG);
                // Noeuds de la face commune
                E_Int indv1 = faceNei[0];
                E_Int indv2 = faceNei[1];
                // Recherche de la face precedente (notee pface et reperee par ipface) sur la cellule d'indice mbv
                E_Int pface, ipface;
                // Remarque : si la face commune == 0 alors on boucle et on retombe sur nbFaces-1 pour ipface
                if (nof1 == 0)
                {
                  ipface = nbFaces-1;
                }
                else
                {
                  ipface = nof1-1;
                }
                pface = elt[ipface]-1;
                // Noeuds de la face precedente (pn1 = previous node 1, idem pour pn2)
                E_Int* facep = cn.getFace(pface, dummy, ngon, indPG);
                E_Int pn1 = facep[0];
                E_Int pn2 = facep[1];
                // Trouver le noeud commun entre la face precedente et la face commune
                // alors l'autre neoud de la face commune devient le first du cycle
                // le noeud commun entre la face commune et la precedente devient le second noeud du cycle
                E_Int drawn=0;
                E_Int prev, j, first, next, n1, n2;
                if (indv1==pn1 || indv1==pn2)
                {
                  indices.push_back(indv2);indices.push_back(indv1);
                  prev = indv2 ; next = indv1 ; first = indv2;
                }
                else
                {
                  indices.push_back(indv1);indices.push_back(indv2);
                  prev = indv1 ; next = indv2 ; first = indv1;
                }
                drawn++;
                // Cherche
                while (drawn < nbFacesNei)
                {
                  for (j = 0; j < nbFacesNei; j++)
                  {
                    E_Int* faceNext = cn.getFace(eltn[j]-1, dummy, ngon, indPG);
                    n1 = faceNext[0];
                    n2 = faceNext[1];
                    if (n1 == next && n2 != prev && n2 != first)
                    { indices.push_back(n2);  prev = n1; next = n2; drawn++; break; }
                    else if (n2 == next && n1 != prev && n1 != first)
                    { indices.push_back(n1); prev = n2; next = n1; drawn++; break; }
                  }
                  if (j == nbFacesNei) drawn++; // pour eviter les boucles infinies
                }
                indices.push_back(indices[0]); // pour boucler
                goto orderie;
              }
            }

          orderie:;
          // Face commune trouvee: ordonne la connectivite pour l'elt ienei
          orderNGONElement(ienei, indices, ngon, nface, indPG, indPH, cn);

          // incremente
          mustBeVisited[mbv] = ienei;
          mbv++; nev++;
          isVisited[ienei] = 1;
        }//test isVisited[ienei] = 0
      }//fin boucle sur les voisins
    }// fin boucle while mbv > 0
  }
  free(isVisited);
  free(mustBeVisited);
  return 1;
}

//=============================================================================
/* Reordonne un elt NGON surfacique de la maniere suivante:
   les sommets sont mis dans l'ordre sous forme d'un cycle (ind1,ind2,ind3)
   Ensuite, les faces sont reordonnees:
   si une face est (ind2,ind1) on la passe en (ind1,ind2)
   Enfin les faces sont triees dans la connectivite Elts/Faces
   IN: noe: numero de l'element
   IN: indices: indices des vertex de l'element a imposer
   IN: posElts, posFaces: indirections sur cNG
   IN/OUT: cNG modifie (retourne)
*/
//=============================================================================
void K_CONNECT::orderNGONElement(E_Int noe, vector<E_Int>& indices,
                                 E_Int* ngon, E_Int* nface, E_Int* indPG,
                                 E_Int* indPH, FldArrayI& cn)
{
  E_Int nvert = indices.size();
  E_Int nf, indm, indp, dummy, nov = 0;
  // Acces universel element noe
  E_Int* elt = cn.getElt(noe, nf, nface, indPH);
  vector<E_Int> sortedFaces;

  while (nov < nvert)
  {
    indm = indices[nov];
    if (nov < nvert-1) indp = indices[nov+1];
    else indp = indices[0];

    for (E_Int nof = 0; nof < nf; nof++)
    {
      // Acces universel face iface
      E_Int iface = elt[nof]-1;
      E_Int* face = cn.getFace(iface, dummy, ngon, indPG);
      E_Int ind1 = face[0]; E_Int ind2 = face[1];

      if ((ind1 == indm && ind2 == indp)||(ind1 == indp && ind2 == indm))
      {
        sortedFaces.push_back(iface+1);
        break;
      }
    }
    nov++;
  }
  E_Int nfaces = sortedFaces.size();
  if (nfaces != nf) {printf("Warning: reorderNGON: NGON connectivity is not clean.");}
  // trier selon sortedFaces les faces ds cEF
  for (E_Int nof = 0; nof < nf; nof++)
  {
    elt[nof] = sortedFaces[nof];
  }
}

void K_CONNECT::orderNGONElement(E_Int noe, vector<E_Int>& indices,
                                 FldArrayI& posElts, FldArrayI& posFaces,
                                 FldArrayI& cNG)
{
  E_Int nvert = indices.size();
  E_Int* ptre = &cNG[posElts[noe]];
  E_Int nf = ptre[0];
  vector<E_Int> sortedFaces;
  E_Int indm, indp;
  E_Int nov = 0;

  while (nov < nvert)
  {
    indm = indices[nov];
    if (nov < nvert-1) indp = indices[nov+1];
    else indp = indices[0];

    for (E_Int nof = 0; nof < nf; nof++)
    {
      E_Int face = ptre[nof+1]-1;
      E_Int* ptrface = &cNG[posFaces[face]];
      E_Int ind1 = ptrface[1]; E_Int ind2 = ptrface[2];

      if ((ind1 == indm && ind2 == indp)||(ind1 == indp && ind2 == indm))
      {
        sortedFaces.push_back(face+1);
        break;
      }
    }
    nov++;
  }
  E_Int nfaces = sortedFaces.size();
  if (nfaces != nf) {printf("Warning: reorderNGON: NGON connectivity is not clean.");}
  // trier selon sortedFaces les faces ds cEF
  for (E_Int nof = 0; nof < nf; nof++)
  {
    ptre[nof+1] = sortedFaces[nof];
  }
}