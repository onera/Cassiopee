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

#ifndef _KCORE_CONNECT_H
#define _KCORE_CONNECT_H
# include <utility>
# include "kPython.h"

# include "Def/DefTypes.h"
# include "Fld/FldArray.h"
# include "Def/DefFunction.h"

# include "topologyMapping.h"
# include "hashFunctions.h"

namespace K_CONNECT
{
  /*-------------*/
  /* - General - */
  /*-------------*/

  /* Suppress identical points in coord */
  void supIdPoints(K_FLD::FldArrayF& coord, 
                   E_Int posx, E_Int posy, E_Int posz,
                   E_Float tol=1.e-14);

  /* Nettoyage de la connectivite de maillage non-structures */
  void cleanConnectivity(E_Int posx, E_Int posy, E_Int posz, 
                         E_Float eps, const char* eltType, 
                         K_FLD::FldArrayF& f, K_FLD::FldArrayI& cEV,
                         bool remove_degen=false,
                         bool ordered_merge=true);

  /* Nettoyage de la connectivite de maillage non-structures 
     (openmp coarse grain) */
  void cleanConnectivity_opt(E_Int posx, E_Int posy, E_Int posz, 
                             E_Float eps, const char* eltType, 
                             K_FLD::FldArrayF& f, K_FLD::FldArrayI& cEV,
                             bool remove_degen=false,
                             bool ordered_merge=true);

  /*-------------*/
  /*- Structure -*/
  /*-------------*/
  /* Reorder a structured field defined by f. Result in f.
   IN/OUT: f: field to be reordered
   IN/OUT: im, jm, km: dimensions of the field in each dir
   IN: oi, oj, ok: new numerotation of the field  
  */
  void reorderStructField(E_Int& im, E_Int& jm, E_Int& km, 
                          K_FLD::FldArrayF& f,
                          E_Int oi, E_Int oj, E_Int ok);

  /* Reorder a structured field defined by f. Result in fout.
   IN: f: field to be reordered
   OUT: fout: reordered field
   IN/OUT: im, jm, km: dimensions of the field in each dir
   IN: oi, oj, ok: new numerotation of the field  
  */
  void reorderStructField(E_Int& im, E_Int& jm, E_Int& km, 
                          K_FLD::FldArrayF& f,
                          K_FLD::FldArrayF& fout,
                          E_Int oi, E_Int oj, E_Int ok); 
  
    /* Detecte le couple de faces coincidentes d'arrays 1 et 2
     IN: im1,jm1,km1: dimensions de l'array 1
     IN: f1: champs de1 : coordonnees incluses
     IN: posx1, posy1, posz1: positions des coordonnees ds f1
     IN: im2,jm2,km2: dimensions de l'array 2
     IN: f2: champs de 2: coordonnees incluses
     IN: posx2, posy2, posz2: positions des coordonnees ds f2
     OUT: nof1: numero de la face de 1 coincidente avec une face de 2 
     OUT: nof2: numero de la face de 2 coincidente avec une face de 1 
     nof1 et nof2 valent -1 initialement
     RETOURNE 1 si il y a coincidence
     0 si pas de face coincidente */
  E_Int detectMatchInterface(E_Int im1, E_Int jm1, E_Int km1, 
                             E_Int posx1, E_Int posy1, E_Int posz1,
                             E_Int im2, E_Int jm2, E_Int km2,
                             E_Int posx2, E_Int posy2, E_Int posz2,
                             K_FLD::FldArrayF& f1, K_FLD::FldArrayF& f2,
                             E_Int& nof1, E_Int& nof2,
                             E_Float eps=1.e-12);
  // version lente
  E_Int detectMatchInterface2(E_Int im1, E_Int jm1, E_Int km1, 
                              E_Int posx1, E_Int posy1, E_Int posz1,
                              E_Int im2, E_Int jm2, E_Int km2,
                              E_Int posx2, E_Int posy2, E_Int posz2,
                              K_FLD::FldArrayF& f1, K_FLD::FldArrayF& f2,
                              E_Int& nof1, E_Int& nof2,
                              E_Float eps=1.e-12);

  /* Determine le coin coincident de f2 avec le coin (i1,j1,k1) de f1 
     attention: les indices demarrent a 1 en i j et k */
  E_Int getCoincidentVertexIndices(E_Int i1, E_Int j1, E_Int k1, 
                                   E_Int im1, E_Int jm1, E_Int km1,  
                                   E_Int im2, E_Int jm2, E_Int km2, 
                                   E_Int posx1, E_Int posy1, E_Int posz1, 
                                   E_Int posx2, E_Int posy2, E_Int posz2,
                                   K_FLD::FldArrayF& f1, K_FLD::FldArrayF& f2, 
                                   E_Int& i2, E_Int& j2, E_Int& k2,
                                   E_Float eps=1.e-12);

  /* Cree un maillage en centres etendus */
  void createExtendedCentersMesh(E_Int imo, E_Int jmo, E_Int kmo, 
                                 E_Float* xt, E_Float* yt, E_Float* zt,
                                 E_Int& ime, E_Int& jme, E_Int& kme, 
                                 K_FLD::FldArrayF& extendedCentersMesh);

  /*---------------------*/
  /* - specifique BARS - */
  /*---------------------*/
  /*
    Supprime les elements definis 2 fois dans une BAR
    IN: xt, yt,zt: coords
    IN: cEV: connectivite element/vertex BAR
    IN: cEVout: connectivite element/vertex sans elements doublons.
  */
  void removeDoubleBARElts(E_Float* xt, E_Float* yt, E_Float* zt, 
                           K_FLD::FldArrayI& cEV, K_FLD::FldArrayI& cEVout);

  /* Ordonnancement des elements de type bar
   IN: cEV: connectivite elt->vertex non ordonnee
   IN: field: champ defini aux noeuds  non ordonnee
   OUT: fieldout: field ordonne
   OUT: cEVout: connectivite elt->vertex ordonnee
  */
  void orderBAR(K_FLD::FldArrayF& field, K_FLD::FldArrayI& cEV,
                K_FLD::FldArrayF& fieldout, K_FLD::FldArrayI& cEVout);

  /* orderBAR2Struct: ordonne une BAR selon un i-array structure. La BAR
   ne doit pas etre ramifiee et ne doit pas contenir de doublons (i.e.
   cleanConnectivity doit etre faite avant). 
   Attention: fout doit etre allouee en externe
   IN: posx, posy, posz: position de x,y,z dans f
   IN: f: tableau a reordonner
   IN: cEV: connectivite Elts/Vertex en BAR
   OUT: fout: tableau reordonne -> i-array */
  void orderBAR2Struct(E_Int posx, E_Int posy, E_Int posz, 
                       K_FLD::FldArrayF& f, K_FLD::FldArrayI& cEV, 
                       K_FLD::FldArrayF& fout);

  /*----------------------------------*/
  /* - Connectivite element basique - */
  /*----------------------------------*/
  // Get dimensionality of a BE/ME connectivity based on its element types
  // Shall ultimately be moved to FldArrayI once _ngon=0 means BE/ME
  E_Int getDimME(const char* eltType);
  E_Int getDimME(std::vector<char*> eltTypes);

  // Get the number of facets per element type of a Multiple Element
  // connectivity. If expandToLowerDim is set to True, 'faces' of 1D and 2D
  // elements are vertices and edges, respectively.
  E_Int getNFPE(std::vector<E_Int>& nfpe, const char* eltType,
                E_Bool expandToLowerDim=true);

  /* Get all facets of a basic element*/
  E_Int getEVFacets(std::vector<std::vector<E_Int> >& facets,
                    const char* eltType, E_Bool allow_degenerated=true);
  
  /* Change a Elts-Vertex connectivity to a Vertex-Elts connectivity.
     cVE doit deja etre alloue au nombre de noeuds. */
  void connectEV2VE(K_FLD::FldArrayI& cEV,
                    std::vector< std::vector<E_Int> >& cVE);

  /* Change a Elts-Vertex connectivity to a Vertex-Vertex neighbours 
     connectivity.
     cVN doit deja etre alloue au nombre de noeuds.
     IN: corners: 0: prend les vertex voisins de V par un edge, 
         corners: 1: prend les vertex de tous les elements voisins de V.
  */
  void connectEV2VNbrs(K_FLD::FldArrayI& cEV,
                       std::vector< std::vector<E_Int> >& cVN, 
                       E_Int corners=0);

  /* Change a Elts-Vertex connectivity to an Element-Element neighbours 
     connectivity
     cEEN doit deja etre alloue au nombre d'elements.
     Attention: les voisins sont ceux ayant une facette commune avec 
     l'element initial. Maillages conformes.
  */
  E_Int connectEV2EENbrs(const char* eltType, E_Int nv, K_FLD::FldArrayI& cEV,
                         std::vector<std::vector<E_Int> >& cEEN);
  E_Int connectEV2EENbrs(const char* eltType, E_Int nv, K_FLD::FldArrayI& cEV,
                         std::vector<std::vector<E_Int> >& cEEN,
                         std::vector<std::vector<E_Int> >& commonFace);
  // Same as connectEV2EENbrs but only returns the number of neighbours
  E_Int connectEV2NNbrs(const char* eltType, E_Int nv, FldArrayI& cEV,
                        std::vector<E_Int>& cENN); 

  /* Change un connectivite Elts-Vertex (basic elements) en une connectivite
   Faces->Vertex. L'indice des faces est global, soit : nof + nelt*nfaces
   ou nfaces est le nbre de faces de l'elements, nelt le no de l'element
   (commencant a 0) et nof la numero local de la face 0,1,2,...
  */
  void connectEV2FV(K_FLD::FldArrayI& cEV, const char* eltType, 
                    K_FLD::FldArrayI& cFV);
  /* Change un connectivite Elts-Vertex (basic elements) en une connectivite
   Vertex->Faces. L'indice des faces est global, soit : nof + nelt*nfaces
   ou nfaces est le nbre de faces de l'elements, nelt le no de l'element
   (commencant a 0) et nof la numero local de la face 0,1,2,...
  */
  void connectEV2VF(K_FLD::FldArrayI& cEV, const char* eltType,
                    std::vector< std::vector<E_Int> >& cVF);
  /* Change HO EV connectivity to LO EV connectivity.
     mode=0: sub-select, mode=1: tesselate */
  E_Int connectHO2LO(const char* eltTypeHO, K_FLD::FldArrayI& cEVHO,
                     K_FLD::FldArrayI& cEVLO, E_Int mode);
  /* identifyFace */
  E_Int identifyFace(E_Int* inds, E_Int n, 
                     std::vector< std::vector<E_Int> >& cVF);

  /* 
     Determine element neighbour for a given element, with a QUAD face of
     indices indA, indB, indC, indD. For HEXA and PENTA elements. Not working
     for non conformal elements.
     IN: indA, indB, indC, indD: indices of the vertices
     IN: cEV: connectivity elemts/vertex
     IN: cEEN: connectivity elements/elements neighbours for the current 
     HEXA/PENTA element
     OUT: number of the neigbour element, -1 if not found
  */
  E_Int getNbrForQuadFace(
    E_Int indA, E_Int indB, E_Int indC, E_Int indD, 
    K_FLD::FldArrayI& cEV, std::vector<E_Int>& cEEN);

  /* Reorder un TRI array suivant dir (+1 ou -1)
     Les triangles sont tous orientes de la meme facon.
     IN: f: field
     IN/OUT: cEV: connectivite elt/vertex TRI
     IN: dir: +1 ou -1
  */
  E_Int reorderQuadTriField(K_FLD::FldArrayF& f, K_FLD::FldArrayI& cEV, 
                            E_Int dir);

  /* 
     Creation de la connectivite Elements -> Noeuds
     IN: cEV: connectivite elt->vertex non nettoyee
     IN: newId: tableau d'indirection entre les noeuds "merges" et les noeuds avant "merge"
     IN: indirp: tableau d'indirection
     OUT: cEVout: connectivite elt->vertex nettoyee
  */
  void createConnectEV(K_FLD::FldArrayI& cEV, std::vector<E_Int> newId, 
                       E_Int* indirp, K_FLD::FldArrayI& cEVout);

  /*
     Version openmp corse grain de createConnectEV
  */
  void createConnectEV_opt(K_FLD::FldArrayI& cEV, std::vector<E_Int> newId, 
                           E_Int* indirp, K_FLD::FldArrayI& cEVout);
  /* 
     Suppression des points doubles et des elements 
     degeneres (surface nulle) dans un array non-structure (f, cEV).
     IN: posx, posy, posz: position des coord. dans f
     IN: eps: tolerance pour eliminer les doublons
     IN: eltType: type d'element du array.
     IN: f, cEV: array non structure a element basique.
  */
  void cleanConnectivityBasic(E_Int posx, E_Int posy, E_Int posz, 
                              E_Float eps, const char* eltType, 
                              K_FLD::FldArrayF& f, K_FLD::FldArrayI& cEV,
                              bool ordered_merge = true);

  /*
     Version openmp corse grain de cleanConnectivityBasic
  */
  void cleanConnectivityBasic_opt(E_Int posx, E_Int posy, E_Int posz, 
                                  E_Float eps, const char* eltType, 
                                  K_FLD::FldArrayF& f, K_FLD::FldArrayI& cEV,
                                  bool ordered_merge = true);

  /* Elimine les vertex non references dans une connectivite basique */
  void cleanUnreferencedVertexBasic(K_FLD::FldArrayF& f, K_FLD::FldArrayI& cn,
                                    K_FLD::FldArrayF& fout, K_FLD::FldArrayI& cnout);

  /*-----------------------*/
  /* - Connectivite NGon - */
  /*-----------------------*/

  /*
    Calcule la connectivite Elts/noeuds cEV a partir de la 
    connectivite NGON cNG. 
    cEV doit etre alloue au nombre d'elements
  */
  void connectNG2EV(K_FLD::FldArrayI& cNG, std::vector< std::vector<E_Int> >& cEV);

  /* Calcule la connectivite Faces/Elts a partir de la connectivite 
     NGON conforme cNG */
  void connectNG2FE(K_FLD::FldArrayI& cNG, K_FLD::FldArrayI& cFE);

  /* Calcule la connectivite Elts/Faces a partir de la connectivite Face/Elts
   pour un NGON conforme */
  void connectFE2EF(K_FLD::FldArrayI& cFE, E_Int nelts, K_FLD::FldArrayI& cEF);

  /* Calcule la connectivite Elts/Elts voisins a partir de la connectivite 
     Face/Elts pour un NGON conforme.
     cEEN doit etre deja dimensionne au nombre d'elements. */
  void connectFE2EENbrs(K_FLD::FldArrayI& cFE, 
                        std::vector< std::vector<E_Int> >& cEEN);

  /* Calcule la connectivite Vertex/Faces a partir de la connectivite NGON 
     pour un NGON conforme. cVF doit etre deja alloue au nb de noeuds. */
  void connectNG2VF(K_FLD::FldArrayI& cNG, 
                    std::vector< std::vector<E_Int> >& cVF);

  /* Calcule la connectivite Vertex/Vertex voisins a partir de la 
     connectivite NGON pour un NGON conforme.
     cVF doit etre alloue au nb de noeuds */
  void connectNG2VNbrs(K_FLD::FldArrayI& cNG, 
                    std::vector< std::vector<E_Int> >& cVF);

  /* Calcule la connectivite Vertex/Elements a partir de la connectivite
      NGON pour un NGON conforme. */
  //void connectNG2ENbrs(E_Int indV, std::vector< std::vector<E_Int> >& cVF, 
  //                     K_FLD::FldArrayI& cFE, std::vector<E_Int>& ENbrs);

  /* Calcule la connectivite NFace (elts->Faces) a partir d'une
     connectivite FE */
  void connectFE2NFace(K_FLD::FldArrayI& cFE, K_FLD::FldArrayI& cNFace, 
                       E_Int& nelts);

  void connectFE2NFace3(K_FLD::FldArrayI& cFE, K_FLD::FldArrayI& cNFace, 
                       K_FLD::FldArrayI& off, E_Int& nelts);

  void connectFE2NFace4(K_FLD::FldArrayI& cFE, K_FLD::FldArrayI& cNFace, 
                       K_FLD::FldArrayI& off, E_Int& nelts);

  /* Calcul des connectivites elements a partir d'une
    connectivite mix */
  void connectMix2EV(K_FLD::FldArrayI& cMIX,
                     K_FLD::FldArrayI& cBAR, K_FLD::FldArrayI& cTRI,
                     K_FLD::FldArrayI& cQUAD, K_FLD::FldArrayI& cTETRA,
                     K_FLD::FldArrayI& cPYRA, K_FLD::FldArrayI& cPENTA,
                     K_FLD::FldArrayI& cHEXA);

  /*
    Calcul la position des faces dans la connectivite generale NGon.
    posFace[0] est la position de la face 0 dans cNG.
    IN: cNG: connectivite NGon.
    OUT: posFaces: position des faces (alloue ici).
  */
  E_Int getPosFaces(K_FLD::FldArrayI& cNG, K_FLD::FldArrayI& posFaces);

  /*
    Calcul la position des elements dans la connectivite generale NGon.
    IN: cNG: connectivite NGon
    OUT: posElts: pour chaque element, sa position dans cNG (alloue ici).
  */
  E_Int getPosElts(K_FLD::FldArrayI& cNG, K_FLD::FldArrayI& posElts);
  
   /*
    Calcul la position des facettes (face pour les elements, noeuds pour les faces) dans la connectivite generale NGon.
    IN: cNG: connectivite NGon
    OUT: posSubs: pour chaque sous entite, sa position dans cNG
  */
  template <typename Container>
  E_Int getPosFacets(const E_Int* data, E_Int facets_start, E_Int nb_facets/*nb of or PG or nodes*/, Container& posFacets)
  {
    posFacets.clear();
    posFacets.resize(nb_facets, 1);
    for (E_Int i = 0; i < nb_facets; i++)
    {
      posFacets[i] = facets_start;
      facets_start += data[facets_start]+1;
    }
    return 1;
  }

  E_Int getIndex(const E_Int* data, E_Int facets_start, 
                 E_Int nb_facets, E_Int* posFacets);

  /* Ordonne un NGON en cyclant les vertices par element */
  E_Int reorderNGON(K_FLD::FldArrayF& f, K_FLD::FldArrayI& cn, E_Int dir);

  /* Reordonne un elt NGON surfacique de la maniere suivante: 
     les sommets sont mis dans l'ordre sous forme d'un cycle (ind1,ind2, ind3)
     Ensuite, les faces sont reordonnees: 
     si une face est (ind2,ind1) on la passe en (ind1,ind2)
     Enfin les faces sont triees dans la connectivite Elts/Faces 
  IN: noe: no de l'elt
  IN: indices: indices des sommets de l elt tries selon une boucle
  IN: posElts, posFaces: pour acceder a l elt ou a la face directement
  IN/OUT: connectivite NGON */
  void orderNGONElement(E_Int noe, std::vector<E_Int>& indices, 
                        K_FLD::FldArrayI& posElts, K_FLD::FldArrayI& posFaces,
                        K_FLD::FldArrayI& cNG);

  void orderNGONElement(E_Int noe, std::vector<E_Int>& indices, 
                        E_Int* ngon, E_Int* nface, E_Int* indPG,
                        E_Int* indPH, K_FLD::FldArrayI& cn);

  // --- Nettoyage de connectivites non-structurees ---
  /* 
     Suppression des points doubles, des faces doubles et des elements 
     degeneres, des faces degenerees  dans un array NGON (f, cNG).
     IN: posx, posy, posz: position des coord. dans f
     IN: eps: tolerance pour eliminer les doublons
     IN: f, cNG: array NGON.
     IN: remove_degen: suppression des elts degeneres compte tenu de la dimension du NGON
  */
  void cleanConnectivityNGon(E_Int posx, E_Int posy, E_Int posz, 
                             E_Float eps, K_FLD::FldArrayF& f, 
                             K_FLD::FldArrayI& cNG,
                             bool remove_degen = false,
                             bool ordered_merge = true);

  /* Nettoyage des faces non referencees dans la connectivity NGon */
  void cleanUnreferencedFacesNGon(K_FLD::FldArrayI& cn);

  PyObject* V_cleanConnectivity(
    const char* varString, K_FLD::FldArrayF& f,
    K_FLD::FldArrayI& cn, const char* eltType,
    E_Float tol=0., E_Bool rmOverlappingPts=true, E_Bool rmOrphanPts=true,
    E_Bool rmDuplicatedFaces=true, E_Bool rmDuplicatedElts=true,
    E_Bool rmDegeneratedFaces=true, E_Bool rmDegeneratedElts=true,
    E_Bool exportIndirPts=false);

  // Clean connectivity - NGON
  PyObject* V_cleanConnectivityNGon(
    E_Int posx, E_Int posy, E_Int posz, const char* varString,
    K_FLD::FldArrayF& f, K_FLD::FldArrayI& cn,
    E_Float tol=0., E_Bool rmOverlappingPts=true, E_Bool rmOrphanPts=true,
    E_Bool rmDuplicatedFaces=true, E_Bool rmDuplicatedElts=true,
    E_Bool rmDegeneratedFaces=true, E_Bool rmDegeneratedElts=true,
    E_Bool exportIndirPts=false);

  E_Int V_identifyDirtyPoints(
    E_Int posx, E_Int posy, E_Int posz, 
    K_FLD::FldArrayF& f, E_Float tol,
    std::vector<E_Int>& indir, E_Bool rmOverlappingPts=true);

  E_Int V_identifyDirtyFacesNGon(
    E_Int dim, K_FLD::FldArrayF& f, K_FLD::FldArrayI& cn,
    std::vector<E_Int>& indirPG, E_Bool rmDegeneratedFaces=true);

  E_Int V_identifyDirtyFacesNGon(
    E_Int dim, K_FLD::FldArrayI& cn, E_Int* ngon, E_Int* indPG,
    std::vector<E_Int>& indirPG, E_Bool rmDegeneratedFaces=true);

  E_Int V_identifyDirtyElementsNGon(
    E_Int dim, K_FLD::FldArrayF& f, K_FLD::FldArrayI& cn,
    std::vector<E_Int>& indirPH, E_Bool rmDegeneratedFaces=true);

  E_Int V_identifyDirtyElementsNGon(
    E_Int dim, K_FLD::FldArrayI& cn, E_Int* nface, E_Int* indPH,
    std::vector<E_Int>& indirPH, E_Bool rmDegeneratedElts=true);

  // Clean connectivity - ME
  PyObject* V_cleanConnectivityME(
    E_Int posx, E_Int posy, E_Int posz, const char* varString,
    K_FLD::FldArrayF& f, K_FLD::FldArrayI& cn, const char* eltType,
    E_Float tol=0., E_Bool rmOverlappingPts=true, E_Bool rmOrphanPts=true,
    E_Bool rmDuplicatedElts=true, E_Bool rmDegeneratedElts=true,
    E_Bool exportIndirPts=false);  

  E_Int V_identifyDirtyElementsME(
    E_Int dim, K_FLD::FldArrayI& cn, std::vector<E_Int>& indir,
    std::vector<E_Int>& nuniqueElts, E_Int neltsTot=0,
    E_Bool rmDegeneratedElts=true);
  
  /* Retourne l'image de vert0 issu de la face0 dans l'element et0 
     IN: cNG: connectivite NGON: Faces/Noeuds et Elts/Faces
     IN: et0: numero de l'element considere dans cNG: demarre a 0 
     IN: vert0: pt de et0, de face0 dont on cherche l'image ds et0: demarre a 1
     IN: face0: numero de la face a laquelle appartient le pt vert0: demarre a 1
     IN: vertices: liste des vertices appartenant a face0
     OUT: vert1: pt image: demarre a 1 */
  E_Int image(E_Int vert0, E_Int face0, E_Int et0,
              std::vector<E_Int>& vertices,
              K_FLD::FldArrayI& cNG);
  E_Int image(E_Int vert0, E_Int face0, E_Int et0,
              std::vector<E_Int>& vertices, K_FLD::FldArrayI& cNG,
              E_Int* ngon, E_Int* nface, E_Int* indPG, E_Int* indPH);

  /* Retourne un tableau donnant pour chaque element sa dimension
     IN: cNG: connectivite NGON: Faces/Noeuds et Elts/Faces
     IN: posFaces: position of each face in cNG
     OUT: dimElts: tableau donnant pour chaque element sa dimension (1,2 ou 3)
  */
  void getDimElts(K_FLD::FldArrayI& cNG, K_FLD::FldArrayI& dimElts);
  void getDimElts(K_FLD::FldArrayI& cNG, K_FLD::FldArrayI& posFaces, 
                  K_FLD::FldArrayI& dimElts);
  void getDimElts(K_FLD::FldArrayI& cNG, E_Int* indPG, E_Int* indPH, 
                  K_FLD::FldArrayI& dimElts);

  /* Pour un element 2D repere par eltPos dans cNG, retourne les indices 
     des vertex de l'element dans l'ordre tournant. 
     les indices commencent a 1 */
  E_Int getVertexIndices(const E_Int* connect, const E_Int* posFaces,
                         E_Int eltPos, 
                         std::vector<E_Int>& ind);
  E_Int getVertexIndices(K_FLD::FldArrayI& cNG, E_Int* ngon, E_Int* nface,
                         E_Int* indPG, E_Int* indPH, E_Int eltPos,
                         std::vector<E_Int>& ind);

  /* ngon tools */
  E_Int check_open_cells(K_FLD::FldArrayI &cn, E_Int *is_cell_open);
  E_Int check_overlapping_cells(K_FLD::FldArrayI &cn);
  E_Int orient_boundary_ngon(E_Float *x, E_Float *y, E_Float *z,
    K_FLD::FldArrayI &cn);
  E_Int build_parent_elements_ngon(K_FLD::FldArrayI &cn, E_Int *owner,
    E_Int *neigh);
  void reversi_connex(E_Int *, E_Int *, E_Int, E_Int *,
    E_Int, std::vector<E_Int> &);
  void build_face_neighbourhood(std::vector<E_Int> &, std::vector<E_Int> &,
    std::vector<E_Int> &);
  E_Int colorConnexParts(E_Int *, E_Int *, E_Int, E_Int *);

  /* Miscellenous */
  // Perform an exclusive prefix sum on an array that is a mask comprised solely
  // of zeros and ones. Return the total number of ones, that is the total number
  // of tagged elements.
  E_Int prefixSum(std::vector<E_Int>& a);
}
#endif
