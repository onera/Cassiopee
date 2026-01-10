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
#include <algorithm>
using namespace K_FLD;
using namespace std;

#define BEGINFACESARR 0
// Nbre de faces
#define NFACES(cn) cn[0]
#define BEGINELTSARR 2+cn[1]
// Nbre d'elements
#define NELTS(cn) cn[BEGINELTSARR]
// Pointeur sur le debut des faces
#define PTRFACES(cn) &cn[2]
// Position du debut des faces dans connect
#define POSFACES(cn) 2
// Pointeur sur le debut des elements
#define PTRELTS(cn) &cn[4+cn[1]]
// Position du debut des elements
#define POSELTS(cn) (4+cn[1])

//=============================================================================
/*
  Calcul la position des faces dans la connectivite generale NGon.
  IN: cn: connectivite NGon
  OUT: posFace: pour chaque face, sa position dans cn (alloue ici).
  Pour array2, retourne une copie. Il vaut mieux utiliser getIndPG.
*/
//=============================================================================
E_Int K_CONNECT::getPosFaces(FldArrayI& cn, FldArrayI& posFaces)
{
  if (cn.getNGonType() == 2) // array2
    return getPosFacets(cn.getNGon(), 0, cn.getNFaces(), posFaces);
  else // array1
    return getPosFacets(cn.begin(), BEGINFACESARR+2, cn[BEGINFACESARR], posFaces);
}

//==============================================================================
E_Int K_CONNECT::getIndex(const E_Int* data, E_Int facetsStart, 
                          E_Int nbFacets, E_Int* posFacets)
{
  for (E_Int i = 0; i < nbFacets; i++)
  {
    posFacets[i] = facetsStart;
    facetsStart += data[facetsStart]+1;
  }
  return 1;
}
//=============================================================================
/*
  Calcul la position des elements dans la connectivite generale NGon.
  IN: cn: connectivite NGon
  OUT: posElt: pour chaque element, sa position dans cn.
  Pour array2, retourne une copie (il vaut mieux utiliser getIndPH)
*/
//=============================================================================
E_Int K_CONNECT::getPosElts(FldArrayI& cn, FldArrayI& posElts)
{
  if (cn.getNGonType() == 2) // array2
    return getPosFacets(cn.getNFace(), 0, cn.getNElts(), posElts);
  else // array1
    return getPosFacets(cn.begin(), BEGINELTSARR+2, cn[BEGINELTSARR], posElts);
}

//=============================================================================
/* Retourne un tableau donnant pour chaque element sa dimension
   IN: cNG: connectivite NGON: Faces/Noeuds et Elts/Faces
   OUT: dimElts: tableau donnant pour chaque element sa dimension (1,2 ou 3) */
//=============================================================================
void K_CONNECT::getDimElts(FldArrayI& cNG, FldArrayI& dimElts)
{
  // Acces non universel sur le ptrs
  E_Int* ngon = cNG.getNGon();
  E_Int* nface = cNG.getNFace();
  E_Int* indPG = cNG.getIndPG();
  E_Int* indPH = cNG.getIndPH();

  // Acces universel nbre d'elements
  E_Int nelts = cNG.getNElts();

  // dimensionnement du tableau dimElts
  dimElts.malloc(nelts);
  E_Int* dimEltsp = dimElts.begin();

  #pragma omp parallel
  {
    E_Int dim;                  // dimension de l element
    E_Int nbFaces;              // nombre de faces pour un element donne
    E_Int nbNodes;              // nombre de noeuds pour une face donnee
    
    #pragma omp for
    for (E_Int i = 0; i < nelts; i++)
    {
      dim = 0;
      // Acces universel element i
      E_Int* elt = cNG.getElt(i, nbFaces, nface, indPH);
      for (E_Int f = 0; f < nbFaces; f++)
      {
        // Acces universel face elt[f]-1
        cNG.getFace(elt[f]-1, nbNodes, ngon, indPG);
        dim = max(nbNodes, dim);
      }
      dim = min(dim, E_Int(3));
      dimEltsp[i] = dim;
    }
  }
}

//=============================================================================
/* Retourne un tableau donnant pour chaque element sa dimension
   IN: cNG: connectivite NGON: Faces/Noeuds et Elts/Faces
   IN: posFaces: position of faces in cNG
   OUT: dimElts: tableau donnant pour chaque element sa dimension (1,2 ou 3) */
//=============================================================================

void K_CONNECT::getDimElts(FldArrayI& cNG, FldArrayI& posFaces,
                           FldArrayI& dimElts)
{
  E_Int* cnp = cNG.begin();
  E_Int sizeFN = cnp[1];      // taille de la connectivite face/noeuds
  E_Int* cEFp = cnp+sizeFN+4; // debut connectivite EF
  E_Int nelts = cnp[sizeFN+2];// nombre d elements
  E_Int dim;                  // dimension de l element
  E_Int nbFaces;              // nombre de faces pour un element donne
  E_Int nbNodes;              // nombre de noeuds pour une face donnee
  E_Int pos;                  // position de la face dans la connectivite
  E_Int* posFacesp = posFaces.begin();

  // dimensionnement du tableau dimElts
  dimElts.malloc(nelts);
  E_Int* dimEltsp = dimElts.begin();

  for (E_Int i = 0; i < nelts; i++)
  {
    nbFaces = cEFp[0]; // nbre de faces pour l'element i
    dim = 0;
    for (E_Int f = 1; f <= nbFaces; f++)
    {
      pos = posFacesp[cEFp[f]-1];
      nbNodes = cnp[pos];
      dim = max(nbNodes, dim);
    }
    dim = min(dim, E_Int(3));
    cEFp += nbFaces+1;
    dimEltsp[i] = dim;
  }
}

// valid pour array et array2 + openMP
void K_CONNECT::getDimElts(FldArrayI& cNG, E_Int* indPG,
                           E_Int* indPH, FldArrayI& dimElts)
{
  E_Int nelts = cNG.getNElts();
  E_Int* ptrf = cNG.getNGon();
  E_Int* ptre = cNG.getNFace();

  E_Int dim;                  // dimension de l'element
  E_Int nbFaces;              // nombre de faces pour un element donne
  E_Int nbNodes;              // nombre de noeuds pour une face donnee
  E_Int pos;                  // position de la face dans la connectivite
  E_Int p;

  // dimensionnement du tableau dimElts
  dimElts.malloc(nelts);
  E_Int* dimEltsp = dimElts.begin();

  #pragma omp parallel default(shared)
  {
    #pragma omp for
    for (E_Int i = 0; i < nelts; i++)
    {
      p = indPH[i];
      nbFaces = ptre[p]; // nbre de faces pour l'element i
      dim = 0;
      for (E_Int f = 1; f <= nbFaces; f++)
      {
        pos = indPG[ptre[p+f]-1];
        nbNodes = ptrf[pos];
        dim = max(nbNodes, dim);
      }
      dim = min(dim, E_Int(3));
      dimEltsp[i] = dim;
    }
  }
}

//=============================================================================
/* Prend un element NGON, rend la liste des vertex dans l'ordre rotatif
   si l'element est 2D, les indices commencent a 1 */
//=============================================================================
E_Int K_CONNECT::getVertexIndices(FldArrayI& cNG, E_Int* ngon, E_Int* nface,
                                  E_Int* indPG, E_Int* indPH, E_Int eltPos,
                                  vector<E_Int>& ind)
{
  // Acces universel element eltPos
  E_Int nf;
  E_Int* elt = cNG.getElt(eltPos, nf, nface, indPH);
  ind.clear();

  // dimension de l'elt
  E_Int nvertex, dummy;
  E_Int dim = 0; E_Int totVertex = 0;
  for (E_Int i = 0; i < nf; i++)
  {
    // Acces universel face elt[i]-1
    cNG.getFace(elt[i]-1, nvertex, ngon, indPG);
    totVertex += nvertex;
    dim = max(nvertex, dim);
  }
  dim = min(dim, E_Int(3));

  if (dim == 1)
  {
    // recupere vertices direct
    ind.reserve(nf);
    for (E_Int i = 0; i < nf; i++)
    {
      E_Int* face = cNG.getFace(elt[i]-1, dummy, ngon, indPG);
      ind.push_back(face[0]);
    }
  }
  else if (dim == 2)
  {
    E_Int drawn = 0;
    E_Int prev, j, first, n1, n2, next;
    E_Int n1f2, n2f2;

    // premiere face, recupere les 2 premiers noeuds
    E_Int* face1 = cNG.getFace(elt[0]-1, dummy, ngon, indPG);
    n1 = face1[0]; n2 = face1[1];
    // deuxieme face, recupere les 2 premiers noeuds
    E_Int* face2 = cNG.getFace(elt[1]-1, dummy, ngon, indPG);
    n1f2 = face2[0]; n2f2 = face2[1];
    
    if ((n1==n1f2)||(n1==n2f2))
    {
      // on trouve deux faces consecutives dans la description nface de l'element avec le noeud n1 en commun. 
      // Ces deux faces donnent l'ordre du parcours des faces sans accorder d'importance a l'orientation
      // de chacune de ces deux faces
      ind.push_back(n2); ind.push_back(n1);
      prev = n2; next = n1; first = n2;
    }
    else
    {
      // Suite des noeuds de l'element en respectant l'ordre donne par la premiere face (en suivant l'ordre de ses deux noeuds)
      
      // Remarque : le traitement du cas ou face1 et face2 sont consecutives et n2 est le noeud commun revient au traitement realise
      //            si l'on suit le sens fournit par n1->n2 de la face1
      ind.push_back(n1); ind.push_back(n2);
      prev = n1; next = n2; first = n1;
    }
    drawn++;

    // Cherche
    while (drawn <= nf)
    {
      for (j = 1; j < nf; j++)
      {
        E_Int* face = cNG.getFace(elt[j]-1, dummy, ngon, indPG);
        n1 = face[0]; n2 = face[1];
        if (n1 == next && n2 != prev && n2 != first)
        { ind.push_back(n2);  prev = n1; next = n2; drawn++; break; }
        else if (n2 == next && n1 != prev && n1 != first)
        { ind.push_back(n1); prev = n2; next = n1; drawn++; break; }
      }
      if (j == nf) drawn++; // pour eviter les boucles infinies
    }
    //if (next != first) ind.push_back(first); // force close
  }
  else // dim == 3
  {
    // direct + unique
    // recupere vertices direct
    ind.reserve(totVertex);
    for (E_Int i = 0; i < nf; i++)
    {
      E_Int* face = cNG.getFace(elt[i]-1, nvertex, ngon, indPG);
      for (E_Int j=0; j < nvertex; j++)
      {
        ind.push_back(face[j]);
      }
    }
    sort(ind.begin(), ind.end());
    ind.erase(unique(ind.begin(), ind.end()), ind.end());
  }
  return 1;
}

E_Int K_CONNECT::getVertexIndices(const E_Int* connect, const E_Int* posFaces,
                                  E_Int eltPos,
                                  vector<E_Int>& ind)
{
  const E_Int* ptrelt = &connect[eltPos];
  E_Int nf = ptrelt[0];
  //std::cout << "nb faces : " << nf << std::endl;
  E_Int face; const E_Int* ptrface;
  E_Int index, nvertex;

  ind.clear();

  // dimension de l'elt
  E_Int dim = 0; E_Int totVertex = 0;
  for (E_Int i = 0; i < nf; i++)
  {
    face = ptrelt[i+1]-1;
    //std::cout << "face : " << face <<std::endl;
    ptrface = &connect[posFaces[face]];
    nvertex = ptrface[0]; totVertex += nvertex;
    dim = max(nvertex, dim);
  }
  dim = min(dim, E_Int(3));

  if (dim == 1)
  {
    // recupere vertices direct
    ind.reserve(nf);
    for (E_Int i = 0; i < nf; i++)
    {
      face = ptrelt[i+1]-1;
      ptrface = &connect[posFaces[face]];
      index = ptrface[1];
      ind.push_back(index);
    }
  }
  else if (dim == 2)
  {
    int drawn=0;
    int prev, j, first, n1, n2, next;
    int n1f2,n2f2,face2;

    // premiere face
    face = ptrelt[1]-1;
    ptrface = &connect[posFaces[face]];
    n1 = ptrface[1]; 
    n2 = ptrface[2];
    // deuxieme face
    face2 = ptrelt[2]-1;
    const E_Int* ptrface2;
    ptrface2 = &connect[posFaces[face2]];
    n1f2 = ptrface2[1]; n2f2 = ptrface2[2];
    
    if ((n1==n1f2)||(n1==n2f2))
    {
      // on trouve deux faces consecutives dans la description nface de l'element avec le noeud n1 en commun. 
      // Ces deux faces donnent l'ordre du parcours des faces sans accorder d'importance a l'orientation
      // de chacune de ces deux faces
      ind.push_back(n2); ind.push_back(n1);
      prev = n2; next = n1; first = n2;
    }
    else
    {
      // Suite des noeuds de l'element en respectant l'ordre donne par la premiere face (en suivant l'ordre de ses deux noeuds)
      
      // Remarque : le traitement du cas ou face1 et face2 sont consecutives et n2 est le noeud commun revient au traitement realise
      //            si l'on suit le sens fournit par n1->n2 de la face1
      ind.push_back(n1); ind.push_back(n2);
      prev = n1; next = n2;first = n1;
    }
    drawn++;

    // Cherche
    while (drawn < nf)
    {
      for (j = 2; j <= nf; j++)
      {
        face = ptrelt[j]-1;
        ptrface = &connect[posFaces[face]];
        n1 = ptrface[1];
        n2 = ptrface[2];
        if (n1 == next && n2 != prev && n2 != first)
        { ind.push_back(n2);  prev = n1; next = n2; drawn++; break; }
        else if (n2 == next && n1 != prev && n1 != first)
        { ind.push_back(n1); prev = n2; next = n1; drawn++; break; }
      }
      if (j == nf+1) drawn++; // pour eviter les boucles infinies
    }
    //if (next != first) ind.push_back(first); // force close
  }
  else // dim == 3
  {
    // direct + unique
    // recupere vertices direct
    ind.reserve(totVertex);
    for (E_Int i = 0; i < nf; i++)
    {
      face = ptrelt[i+1]-1;
      ptrface = &connect[posFaces[face]];
      nvertex = ptrface[0];
      for (E_Int j=0; j < nvertex; j++)
      {
        index = ptrface[j+1];
        ind.push_back(index);
      }
    }
    sort(ind.begin(), ind.end());
    ind.erase(unique(ind.begin(), ind.end()), ind.end());
  }
  return 1;
}