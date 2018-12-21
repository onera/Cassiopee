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
#ifndef _CONNECTOR_BLKINTERP_H_
#define _CONNECTOR_BLKINTERP_H_

# include "Def/DefTypes.h"
# include "Def/DefFunction.h"
# include "Fld/FldArray.h"
# include "KMesh.h"
# include "BlkIntTreeNode.h"
# include "Array/Array.h"
# include <vector>
# include <list>
# define FldArrayF K_FLD::FldArrayF
# define FldArrayI K_FLD::FldArrayI
# define FldArrayIS K_FLD::FldArrayIS

namespace K_KINTERP
{
#include "BlkInterpData.h"
#include "BlkInterpWithKMesh.h"
#include "BlkInterpAdt.h"

  void boundingBox(std::vector<KMesh*>& KMeshList,
                   E_Float& xmax, E_Float& ymax, E_Float& zmax, 
                   E_Float& xmin, E_Float& ymin, E_Float& zmin);
  
/* Calcul de la cellule d'interpolation pour une liste de blocs
   Retourne le no du bloc d interpolation > 0 si trouvé
                                          = 0 si pas de donneur trouvé
   Si structuré : noblk : no du bloc associés aux interpData structurées interpDatas
   Si donneur non structuré : noblk shifté de la taille de interpDatas 

   IN: x,y,z coordonnees du point a interpoler
   IN: interpDatas: liste des interpData structurees
   IN: fields: liste des champs pour les blocs d'interpolation: 
               contiennent x,y,z et eventuellement celln.
   IN: nis, njs, nks: tableaux des dimensions des grilles associees a fields
   IN: posxs: position de x dans chaque element de fields              
   IN: posys: position de y dans chaque element de fields              
   IN: poszs: position de z dans chaque element de fields              
   IN: poscs: position de celln dans chaque element de fields. 
              si n existe pas, posc vaut 0              
   IN: interpMeshType: maillages d'interpolation en noeuds ou centres etendus
   IN: interpType: type d interpolation (ordre et stockage des coefs)
   IN: nature : si nature=0, alors une cellule est valide si ne contient pas de pt masque
                si nature=1, une cellule donneuse est valide si ne contient que des cellN=1
   IN: penalty : penalisation sur le volume de la cellule donneuse si est sur un bord
   OUT: indi: indices des sommets de la cellule d interpolation
   OUT: cf: coefficients d interpolation
   Attention : dans le cas ou penalty = 1, la connectivite elts/elts Voisins est construite systematiquement.
   A AMELIORER
 */
  short getInterpolationCell(
    E_Float x, E_Float y, E_Float z,
    std::vector<BlkInterpData*>& interpDatas,
    std::vector<FldArrayF*>& fields,
    std::vector<E_Int>& nis, std::vector<E_Int>& njs,
    std::vector<E_Int>& nks,
    std::vector<E_Int>& posxs, std::vector<E_Int>& posys, 
    std::vector<E_Int>& poszs, std::vector<E_Int>& poscs,
    std::vector<BlkInterpData*>& interpDatau,
    std::vector<FldArrayF*>& fieldu, std::vector<FldArrayI*>& connectu,
    std::vector<E_Int>& posxu, std::vector<E_Int>& posyu, 
    std::vector<E_Int>& poszu, std::vector<E_Int>& poscu,
    FldArrayI& indi, FldArrayF& cf,
    BlkInterpData::InterpMeshType interpMeshType,
    BlkInterpData::InterpolationType interpType,
    E_Int nature=0, E_Int penalty=0);

/* MEME CHOSE MAIS RETOURNE LE VOLUME DE LA CELLULE D'INTERPOLATION 
   OUT: voli: volume de la cellule d'interpolation
 */
  short getInterpolationCell(
    E_Float x, E_Float y, E_Float z,
    std::vector<BlkInterpData*>& interpDatas,
    std::vector<FldArrayF*>& fields,
    std::vector<E_Int>& nis, std::vector<E_Int>& njs,
    std::vector<E_Int>& nks,
    std::vector<E_Int>& posxs, std::vector<E_Int>& posys, 
    std::vector<E_Int>& poszs, std::vector<E_Int>& poscs,
    std::vector<BlkInterpData*>& interpDatau,
    std::vector<FldArrayF*>& fieldu, std::vector<FldArrayI*>& connectu,
    std::vector<E_Int>& posxu, std::vector<E_Int>& posyu, 
    std::vector<E_Int>& poszu, std::vector<E_Int>& poscu,
    E_Float& voli, FldArrayI& indi, FldArrayF& cf,
    BlkInterpData::InterpMeshType interpMeshType,
    BlkInterpData::InterpolationType interpType,
    E_Int nature=0, E_Int penalty=0);

/* Calcul de la cellule d'interpolation pour une liste de blocs structures
   Retourne le no du bloc d interpolation (demarre a 1): 0 si pas trouve
   IN: x,y,z coordonnees du point a interpoler
   IN: interpDatas: liste des interpData structurees
   IN: fields: liste des champs pour les blocs d interpolation: 
                 contiennent x,y,z et eventuellement celln.
   IN: nis, njs, nks : tableaux des dimensions des grilles associees a fields
   IN: posxs: position de x dans chaque element de fields              
   IN: posys: position de y dans chaque element de fields              
   IN: poszs: position de z dans chaque element de fields              
   IN poscs: position de celln dans chaque element de fields. 
             si n existe pas, posc vaut 0              
   IN: interpMeshType: maillages d'interpolation en noeuds ou centres etendus
   IN: interpType: type d'interpolation (ordre et stockage des coefs)
   OUT: indi: indices des sommets de la cellule d interpolation
   OUT: cf: coefficients d interpolation
 */
  short getInterpolationCell(
    E_Float x, E_Float y, E_Float z,
    std::vector<BlkInterpData*>& interpDatas,
    std::vector<FldArrayF*>& fields,
    std::vector<E_Int>& nis, std::vector<E_Int>& njs,
    std::vector<E_Int>& nks,
    std::vector<E_Int>& posxs, std::vector<E_Int>& posys, 
    std::vector<E_Int>& poszs, std::vector<E_Int>& poscs,
    FldArrayI& indis, FldArrayF& cfs,
    BlkInterpData::InterpMeshType interpMeshType,
    BlkInterpData::InterpolationType interpType,
    E_Int nature=0, E_Int penalty=0);
  
/* Recherche de cellule d interpolation a partir de maillages non structures*/
  short getInterpolationCell(
    E_Float x, E_Float y, E_Float z,
    std::vector<BlkInterpData*>& interpDatau,
    std::vector<FldArrayF*>& fieldu, std::vector<FldArrayI*>& connectu,
    std::vector<E_Int>& posxu, std::vector<E_Int>& posyu, 
    std::vector<E_Int>& poszu, std::vector<E_Int>& poscu,
    FldArrayI& indiu, FldArrayF& cfu,
    E_Int nature=0, E_Int penalty=0);

/*  Calcul de la cellule d'extrapolation pour une liste de blocs structures 
    et/ou non-structures
    Retourne le no du bloc d interpolation > 0 si trouvé
                                           = 0 si pas de donneur trouvé
   Si structuré : noblk : no du bloc associés aux interpData structurées interpDatas
   Si donneur non structuré : noblk shifté de la taille de interpDatas 

   IN: x,y,z coordonnees du point a interpoler
   IN: interpDatas: liste des interpData structurees
   IN: fields: liste des champs pour les blocs d'extrapolation structures: 
               contiennent x,y,z et eventuellement celln.
   IN: nis, njs, nks: tableaux des dimensions des grilles associees a fields
   IN: posxs: position de x dans chaque element de fields              
   IN: posys: position de y dans chaque element de fields              
   IN: poszs: position de z dans chaque element de fields              
   IN: poscs: position de celln dans chaque element de fields. 
              si n existe pas, posc vaut 0              
   IN: interpDatau: liste des interpData non-structurees
   IN: eltType: type d'elements
   IN: fieldu: liste des champs pour les blocs d extrapolation non-structures: 
               contiennent x,y,z et eventuellement celln.
   IN: connectu: liste des connectivites pour les blocs d extrapolation non-structures
   IN: posxu: position de x dans chaque element de fieldu              
   IN: posyu: position de y dans chaque element de fieldu              
   IN: poszu: position de z dans chaque element de fieldu              
   IN: poscu: position de celln dans chaque element de fieldu. 
              si n existe pas, posc vaut 0                 
   IN: interpMeshType: maillages d'interpolation en noeuds ou centres etendus
   IN: interpType: type d interpolation (ordre et stockage des coefs)
   OUT: voli: volume de la cellule d interpolation
   OUT: indi: indices des sommets de la cellule d interpolation
   OUT: cf: coefficients d interpolation
 */
  short getExtrapolationCell(
    E_Float x, E_Float y, E_Float z,
    std::vector<BlkInterpData*>& interpDatas,
    std::vector<FldArrayF*>& fields,
    std::vector<E_Int>& nis, std::vector<E_Int>& njs,
    std::vector<E_Int>& nks,
    std::vector<E_Int>& posxs, std::vector<E_Int>& posys, 
    std::vector<E_Int>& poszs, std::vector<E_Int>& poscs,
    std::vector<BlkInterpData*>& interpDatau,
    std::vector<char*> eltType,
    std::vector<FldArrayF*>& fieldu, std::vector<FldArrayI*>& connectu,
    std::vector<E_Int>& posxu, std::vector<E_Int>& posyu, 
    std::vector<E_Int>& poszu, std::vector<E_Int>& poscu,
    E_Float& voli, FldArrayI& indi, FldArrayF& cf, E_Float cfMax,
    BlkInterpData::InterpMeshType interpMeshType,
    BlkInterpData::InterpolationType interpType);


/* Creation de la liste des interpData pour des maillages non structures
   IN: fieldsIn: liste des champs contenant les infos sur les grilles
   IN: varStringIn: liste des varstring correspondantes
   IN: cnIn: connectivite 
   IN: eltsIn: liste des types d elements 
   OUT: fieldsOut: liste des champs contenant les infos sur les grilles
   OUT: varStringOut: liste des varstring correspondantes
   IN: cnOut: connectivite 
   IN: eltTypesOut: liste des types d elements 
   OUT: posxt, posyt, poszt, posct : positions de x,y,z,celln (eventuellt)
   OUT: listOfUnstrMeshes: c est pour pouvoir les detruire 
   OUT: interpDatas: liste des interpDatas crees basees sur listOfMeshes
*/
  void buildListOfUnstrInterpData(
    std::vector<FldArrayI*>& connectIn,
    std::vector<FldArrayF*>& fieldsIn,
    std::vector<char*>& varStringIn,
    std::vector<char*>& eltsIn, 
    std::vector<FldArrayI*>& connectOut,
    std::vector<FldArrayF*>& fieldsOut,
    std::vector<char*>& varStringOut,
    std::vector<char*>& eltsOut, 
    std::vector<E_Int>& posxt, 
    std::vector<E_Int>& posyt,
    std::vector<E_Int>& poszt,
    std::vector<E_Int>& posct,
    std::vector<KMesh*>& listOfUnstrMeshes,
    std::vector<BlkInterpData*>& interpDatau);

 /* Creation de la liste des interpData pour des maillages structures
   IN : fieldsIn : liste des champs contenant les infos sur les grilles
   IN : varStringIn : liste des varstring correspondantes
   IN : nitin, njtin, nktin : dimensions des grilles
   OUT : fieldsOut : liste des champs contenant les infos sur les grilles
   OUT : varStringOut : liste des varstring correspondantes
   OUT : nitout, njtout, nktout : dimensions des grilles
   OUT : posxt, posyt, poszt, posct : positions de x,y,z,celln (eventuellt)
   OUT : listOfStructMeshes : c est pour pouvoir les detruire 
   OUT : interpDatas : liste des interpDatas crees basees sur listOfMeshes
*/
  void buildListOfStructInterpData(
    std::vector<FldArrayF*>& fieldsIn,
    std::vector<char*>& varStringIn,
    std::vector<E_Int>& nitin, 
    std::vector<E_Int>& njtin,
    std::vector<E_Int>& nktin,
    std::vector<FldArrayF*>& fieldsOut,
    std::vector<char*>& varStringOut,
    std::vector<E_Int>& nitout, 
    std::vector<E_Int>& njtout,
    std::vector<E_Int>& nktout,
    std::vector<E_Int>& posxt, 
    std::vector<E_Int>& posyt,
    std::vector<E_Int>& poszt,
    std::vector<E_Int>& posct,
    std::vector<KMesh*>& listOfStructMeshes,
    std::vector<BlkInterpData*>& interpDatas);
  
                              
/* Calcul de la cellule d interpolation la plus petite ds les grilles structurees
   nature = 0 : traitement adapte a extractMesh : un donneur est valide si pas de pt masque
   nature = 1 : traitement type setInterpolations : un donneur est valide si pas de pt masque ni interpole
   penalty = 1: penalité sur le volume de la cellule donneuse si elle appartient a une frontiere*/
  E_Float selectBestStructuredInterpolationCell( 
    E_Float x, E_Float y, E_Float z,
    std::vector<BlkInterpData*>& interpDatas,
    std::vector<FldArrayF*>& fields,
    std::vector<E_Int>& nis, std::vector<E_Int>& njs, std::vector<E_Int>& nks,
    std::vector<E_Int>& posx, std::vector<E_Int>& posy, 
    std::vector<E_Int>& posz, std::vector<E_Int>& posc,
    FldArrayI& indi, 
    FldArrayF& cf, E_Int& noblks,
    BlkInterpData::InterpMeshType interpMeshType,
    BlkInterpData::InterpolationType interpType,
    E_Int nature = 0, E_Int penalty=0);

  /* Calcul de la cellule d interpolation la plus petite ds les grilles non structurees
     nature = 0 : traitement adapte a extractMesh : un donneur est valide si pas de pt masque
     nature = 1 : traitement type setInterpolations : un donneur est valide si pas de pt masque ni interpole
     penalty = 1: penalité sur le volume de la cellule donneuse si elle appartient a une frontiere 
     Attention : dans le cas ou penalty = 1, la connectivite elts/elts Voisins est construite systematiquement.
     A ameliorer
  */
  E_Float selectBestUnstructuredInterpolationCell( 
    E_Float x, E_Float y, E_Float z,
    std::vector<BlkInterpData*>& interpDatau,
    std::vector<FldArrayF*>& fieldu, std::vector<FldArrayI*>& connectu,
    std::vector<E_Int>& posx, std::vector<E_Int>& posy, 
    std::vector<E_Int>& posz, std::vector<E_Int>& posc,
    FldArrayI& indi, FldArrayF& cf, E_Int& noblku,
    E_Int nature=0, E_Int penalty=0);

/* Calcul de la cellule d extrapolation la plus petite ds les grilles structurees*/
  E_Float selectBestStructuredExtrapolationCell( 
    E_Float x, E_Float y, E_Float z,
    std::vector<BlkInterpData*>& interpDatas,
    std::vector<FldArrayF*>& fields,
    std::vector<E_Int>& nis, std::vector<E_Int>& njs, std::vector<E_Int>& nks,
    std::vector<E_Int>& posx, std::vector<E_Int>& posy, 
    std::vector<E_Int>& posz, std::vector<E_Int>& posc,
    FldArrayI& indi, 
    FldArrayF& cf, E_Int& noblks, E_Float cfMax,
    BlkInterpData::InterpMeshType interpMeshType,
    BlkInterpData::InterpolationType interpType);

  /* Calcul de la cellule d extrapolation la plus petite ds les grilles non structurees*/
  E_Float selectBestUnstructuredExtrapolationCell( 
    E_Float x, E_Float y, E_Float z,
    std::vector<BlkInterpData*>& interpDatau, std::vector<char*> eltType,
    std::vector<FldArrayF*>& fieldu, std::vector<FldArrayI*>& connectu,
    std::vector<E_Int>& posx, std::vector<E_Int>& posy, 
    std::vector<E_Int>& posz, std::vector<E_Int>& posc,
    FldArrayI& indi, FldArrayF& cf, E_Int& noblku);

  /* Interpole une liste de champs f0 dans f pour un point a interpoler d indice ind
     Le type d interpolation effectuee est dans interpType
     IN : ni, nj, nk : dimensions de f0 dans le cas structure
     -1, -1, -1 si le donneur est non structure
     IN : cn0 : connectivite elts/noeuds associee a f0
     IN : f0 : champ du bloc donneur pour l interpolation
     IN : indi : indices des pts de la molecule d interpolation
                 peuvent être definis par des sommets de la molécule donneuse
     IN : indiSize : taille de indi
     IN : cf : coefs d'interpolation 
     IN : ind  : indice du pt a interpoler 
     IN : interpType : permet de determiner la formule appliquee
     OUT : f : champs interpoles au point d indice ind 
     Retourne -1 si erreur */
  short compInterpolatedValues(E_Int ni, E_Int nj, E_Int nk,  
                               E_Int indiSize, E_Int* indi, FldArrayF& cf,
                               FldArrayF& f0,  FldArrayI& cn0, 
                               E_Int ind, FldArrayF& f,
                               E_Int interpType);


 /* Calcul d'une valeur interpolee en un pt a partir du champ nofld de f0
   IN : ni, nj, nk : dimensions de la grille d interpolation structuree
   IN : indiSize : taille de indi (indices des pts +1 pour le type)
   IN : indi : type + indices des pts d interpolation
   IN : cf : coefs d interpolation
   IN : cellNp : cellnaturefield  du bloc donneur
   OUT : valeur interpolee du champ cellN en ind */
  E_Float compInterpolatedNature(E_Int ni, E_Int nj, E_Int nk,
                                 E_Int indiSize, E_Int* indi, 
                                 FldArrayF& cf, E_Int* cellNp);  
 /* Calcul d'une valeur interpolee en un pt a partir du champ nofld de f0
   IN : ni, nj, nk : dimensions de la grille d interpolation structuree
   IN : indiSize : taille de indi (indices des pts +1 pour le type)
   IN : indi : type + indices des pts d interpolation
   IN : cf : coefs d interpolation
   IN : cellNp : cellnaturefield  du bloc donneur
   OUT : valeur interpolee du champ cellN en ind */
  E_Float compInterpolatedNatureEX(E_Int ni, E_Int nj, E_Int nk,
                                   E_Int indiSize, E_Int* indi, 
                                   FldArrayF& cf, E_Int* cellNp);  
}
#undef FldArrayF
#undef FldArrayI
#undef FldArrayIS 
#endif
//========================KCore/Interp/BlkInterp.h===========================
