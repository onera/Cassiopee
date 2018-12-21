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
  E_Int* cnp = cn.begin();
  E_Int sizeFN = cnp[1];
  E_Int nelts = cnp[sizeFN+2];
  //connectivite ets/ets voisins
  vector< vector<E_Int> > cEEN(nelts);
  FldArrayI cFE; connectNG2FE(cn, cFE);

  connectFE2EENbrs(cFE, cEEN);
  FldArrayI posElts; K_CONNECT::getPosElts(cn, posElts);
  FldArrayI posFaces; K_CONNECT::getPosFaces(cn, posFaces);

  E_Int nev = 0; // nbre d'elements deja visites
  char* isVisited = (char*)calloc(nelts, sizeof(char)); // elt deja visite?
  E_Int* mustBeVisited = (E_Int*)malloc(nelts * sizeof(E_Int));
  E_Int mbv = 0;
  E_Int p, size, elt,iv, ie;
  vector<E_Int> indices;

  while (nev < nelts)
  {
    // Recherche du premier elt pas encore visite
    for (p = 0; (isVisited[p] != 0); p++);
    // appartient a un nouveau composant connexe
    // alors on recupere les vertex dans l'ordre rotatif et on ordonne l'elt de depart
    K_CONNECT::getVertexIndices(cn.begin(), posFaces.begin(), posElts[p], indices);

/*    printf("vertices initiaux\n");
    for (E_Int i = 0; i < indices.size(); i++) printf("%d ", indices[i]);
    printf("\n");
*/
    if (dir == -1) reverse(indices.begin(),indices.end());
    orderNGONElement(p, indices, posElts, posFaces, cn);

/*    // DBX
    K_CONNECT::getVertexIndices(cn.begin(), posFaces.begin(), posElts[p], indices);
    printf("vertices finaux\n");
    for (E_Int i = 0; i < indices.size(); i++) printf("%d ", indices[i]);
    printf("\n");
    // ENDDBX
*/
    indices.clear();

    mustBeVisited[mbv] = p;
    mbv++; nev++;
    isVisited[p] = 1;

    while (mbv > 0)
    {
      mbv--;
      elt = mustBeVisited[mbv];
      size = cEEN[elt].size();
      E_Int* ptre = &cn[posElts[elt]];
      for (iv = 0; iv < size; iv++)
      {
        ie = cEEN[elt][iv]; // no de l'elt voisin
        if (isVisited[ie] == 0)
        {
          E_Int* ptren = &cn[posElts[ie]];
          indices.clear();

          // Trouve 1 face commune entre les deux elements
          for (E_Int nof1 = 1; nof1 <= ptre[0]; nof1++)
            for (E_Int nof2 = 1; nof2 <= ptren[0]; nof2++)
            {
              if (ptre[nof1] == ptren[nof2]) // face commune
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
                |    |         mbv          |     |              ie               |
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
                E_Int face = ptren[nof2]-1;
                E_Int* ptrface = &cn[posFaces[face]];
                // Noeuds de la face commune
                E_Int indv1 = ptrface[1];//demarre a 1
                E_Int indv2 = ptrface[2];
                // Recherche de la face precedente (notee pface et reperee par ipface) sur la cellule d'indice mbv
                E_Int pface,ipface;
                // Remarque : si la face commune == 1 alors on boucle et on retombe sur nbface pour ipface
                if (nof1==1)
                {
                  ipface = ptre[0];
                }
                else
                {
                  ipface = nof1-1;
                }
                pface = ptre[ipface]-1;
                // Noeuds de la face precedente (pn1 = previous node 1, idem pour pn2)
                E_Int* ptrpface = &cn[posFaces[pface]];
                E_Int pn1       = ptrpface[1];
                E_Int pn2       = ptrpface[2];
                // Trouver le noeud commun entre la face precedente et la face commune
                // alors l'autre neoud de la face commune devient le first du cycle
                // le noeud commun entre la face commune et la precedente devient le second noeud du cycle
                E_Int drawn=0;
                E_Int prev, j, first, next, n1, n2, nextface;
                const E_Int* ptrnextface;
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
                while (drawn < ptren[0])
                {
                  for (j = 1; j <= ptren[0]; j++)
                  {
                    nextface = ptren[j]-1;
                    ptrnextface = &cn[posFaces[nextface]];
                    n1 = ptrnextface[1];
                    n2 = ptrnextface[2];
                    if (n1 == next && n2 != prev && n2 != first)
                    { indices.push_back(n2);  prev = n1; next = n2; drawn++; break; }
                    else if (n2 == next && n1 != prev && n1 != first)
                    { indices.push_back(n1); prev = n2; next = n1; drawn++; break; }
                  }
                  if (j == ptren[0]+1) drawn++; // pour eviter les boucles infinies
                }
                indices.push_back(indices[0]); // pour boucler
                goto orderie;
              }
            }

          orderie:;
          // Face commune trouvee: ordonne la connectivite pour l'elt ie
          orderNGONElement(ie, indices, posElts, posFaces, cn);

          // incremente
          mustBeVisited[mbv] = ie;
          mbv++; nev++;
          isVisited[ie] = 1;
        }//test isVisited[ie] = 0
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
