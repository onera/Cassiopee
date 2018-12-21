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
// Interface to pytree data

#ifndef _KCORE_PYTREE_H_
#define _KCORE_PYTREE_H_
#include "kPython.h"
#include "Def/DefTypes.h"
#include "Def/DefCplusPlusConst.h"

#include "Numpy/importNumpy.h"
#include <numpy/arrayobject.h>
#include "Fld/FldArray.h"
#include "Fld/DynArray.h"
#include <vector>

#define FldArrayF K_FLD::FldArrayF
#define FldArrayI K_FLD::FldArrayI

// Release shared zone
#define RELEASESHAREDZ(h, varString, eltType) {delete [] varString; delete [] eltType; for (unsigned int i = 0; i < h.size(); i++) Py_DECREF(h[i]); }
// Release hook
#define RELEASEHOOK(h) {for (unsigned int i = 0; i < h.size(); i++) Py_DECREF(h[i]); h.clear(); }

namespace K_PYTREE
{
  /* Recupere les donnees d'une zone pyTree (partagees avec python)
     IN: o: la zone python
     IN: xyz: 1: on prend les coords
     IN: loc: 0: en noeuds, 1: en centres, 2: les deux
     IN: GridCoordinates: nom du container de coordonnees
     IN: FlowSolutionNodes: nom du container de champs en noeuds
     IN: FlowSolutionCenters: nom du container de champs en centres
     OUT: ni,nj,nk: les ni, nj,nk de la zone
     OUT: cn: connectivite si zone non structuree
     OUT: csize: taille du tableau cn
     OUT: cnfld: nbre de champs dans le tableau cn
     OUT: fields: un vecteur des pointeurs des champs (partage avec numpy)
     OUT: locs: localisation des champs locs[n]=0 (en noeuds), =1 (en centres)
     Retourne 1: zone structuree valide
     Retourne 2: zone non structuree valide
     Retourne 0: zone invalide
  */
#define GETFROMZONEDATA \
  char* varString;                        \
  std::vector<E_Float*> fields;           \
  std::vector<E_Int> locs;                \
  E_Int ni, nj, nk;                       \
  std::vector<E_Int*> cn;                 \
  E_Int csize, cfld;                      \
  char* eltType;                          
  E_Int getFromZone(PyObject* o, E_Int xyz, E_Int loc,
                    char*& varString,
                    std::vector<E_Float*>& fields,
                    std::vector<E_Int>& locs,
                    E_Int& ni, E_Int& nj, E_Int& nk,
                    std::vector<E_Int*>& cn, 
                    E_Int& csize, E_Int& cnfld,
                    char*& eltType,
                    std::vector<PyArrayObject*>& hook,
                    char* GridCoordinates=NULL,
                    char* FlowSolutionNodes=NULL,
                    char* FlowSolutionCenters=NULL);
  /* Recherche par nom d'un seul niveau, retourne un seul noeud.
     IN: o: objet representant un noeud de pyTree
     IN: name: le nom du noeud
     OUT: retourne le noeud trouve, sinon retourne NULL. */
  PyObject* getNodeFromName1(PyObject* o, const char* name);
  /* Recherche par type d'un seul niveau, retourne un vecteur de noeuds.
     IN: o: objet representant un noeud de pyTree
     IN: type: le type du noeud
     OUT: out: la liste des noeuds trouves. */
  void getNodesFromType1(PyObject* o, const char* type, 
                         std::vector<PyObject*>& out);

  /* Recherche par path, retourne un seul noeud.
     IN: o: objet representant un noeud de pyTree (depart)
     IN: path: le chemin
     OUT: retourne le nouveau noeud si trouve, sinon retourne NULL */
  PyObject* getNodeFromPath(PyObject* o, const char* path);

  /* Retourne le nom d'un noeud (partage avec python) */
  char* getNodeName(PyObject* o, std::vector<PyArrayObject*>& hook);

  /* Retourne le type d'un noeud (partage avec python) */
  char* getNodeType(PyObject* o, std::vector<PyArrayObject*>& hook);

  /* Recupere la chaine de char si le noeud o contient une string 
     ou un numpy char. Le tableau de char retourne est partage avec python.
     Attention: dans le cas de numpy, il n'y a pas de \0 a la fin.
     IN: o: noeud du pyTree
     OUT: hook: contient les numpy char
     OUT: retourne le ptr sur la chaine partagee. */
  char* getValueS(PyObject* o, std::vector<PyArrayObject*>& hook);
  /* Recupere la chaine de char si le noeud o contient une string ou 
     un numpy char. Le tableau de char retourne est partage avec python.
     Attention: dans le cas de numpy, il n'y a pas de \0 a la fin.
     IN: o: noeud du pyTree
     OUT: s: taille de la chaine retournee
     OUT: hook: contient les numpy char
     OUT: retourne le ptr sur la chaine partagee. */
  char* getValueS(PyObject* o, E_Int& s, std::vector<PyArrayObject*>& hook);
  /* Retourne un entier contenu dans un noeud
     IN: o: noeud de pyTree
     OUT: entier */
  E_Int getValueI(PyObject* o);
  /* Retourne un float contenu dans un noeud
     IN: o: noeud de pyTree
     OUT: float */
  E_Float getValueF(PyObject* o);
  /* Retourne un pointeur sur un tableau d'entier si le noeud o contient
     un numpy d'entier
     IN: o: noeud de pyTree
     OUT: hook: contient les numpy d'entiers
     OUT: retourne le ptr sur le tableau d'entiers. */
  E_Int* getValueAI(PyObject* o, std::vector<PyArrayObject*>& hook);
  /* Retourne un pointeur sur un tableau d'entier si le noeud o contient
     un numpy d'entiers
     IN: o: noeud de pyTree
     OUT: s0, s1: taille du tableau d'entiers (s0,s1)
     OUT: hook: contient les numpy d'entiers
     OUT: retourne le ptr sur le tableau d'entiers. */  
  E_Int* getValueAI(PyObject* o, E_Int& s0, E_Int& s1, 
                    std::vector<PyArrayObject*>& hook);
  /* Retourne un pointeur sur un tableau de doubles si le noeud o contient
     un numpy de doubles
     IN: o: noeud de pyTree
     OUT: hook: contient les numpy de doubles
     OUT: retoure le ptr sur le tableau de doubles. */
  E_Float* getValueAF(PyObject* o, std::vector<PyArrayObject*>& hook);

  /* Ajoute un noeud a o
     IN: o: objet representant un noeud de pyTree
     IN: name: le nom du noeud a creer
     IN: type: le type du noeud a creer
     IN: value: la valeur a lui donner
     IN: size, nfld: dimension du tableau value, si necessaire */
  PyObject* createChild(PyObject* o, const char* name, 
                        const char* type, 
                        E_Float* value, E_Int size, E_Int nfld,
                        E_Int pos=-1);
  PyObject* createChild(PyObject* o, const char* name, 
                        const char* type, 
                        E_Int* value, E_Int size, E_Int nfld,
                        E_Int pos=-1);
  PyObject* createChild(PyObject* o, const char* name, 
                        const char* type, 
                        PyObject* value, 
                        E_Int pos=-1);

  /* Ajoute variable a la varString
     Si varString n'existe pas, il faut envoyer varString=NULL. */
  void add2VarString(char*& varString, E_Int& size, const char* variable);
}

#undef FldArrayF
#undef FldArrayI
#endif

