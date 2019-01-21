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
// Interface to Cassiopee arrays

#ifndef _KCORE_ARRAY_H_
#define _KCORE_ARRAY_H_
#include "Fld/DynArray.h"
#include "kPython.h"
#include "Def/DefTypes.h"
#include "Def/DefCplusPlusConst.h"

#include "Numpy/importNumpy.h"
#include <numpy/arrayobject.h>
#include "Fld/FldArray.h"
#include <vector>

#define FldArrayF K_FLD::FldArrayF
#define FldArrayI K_FLD::FldArrayI

// Release shared structured array
#define RELEASESHAREDS(array, f) {delete f; Py_DECREF(PyList_GetItem(array,1));}
// Release shared unstructured array
#define RELEASESHAREDU(array, f, cn) {delete f; delete cn;              \
    Py_DECREF(PyList_GetItem(array,1)); Py_DECREF(PyList_GetItem(array,2));}
// Release shared structured/unstructured array depending on res
#define RELEASESHAREDB(res, array, f, cn) {delete f;                  \
  Py_DECREF(PyList_GetItem(array,1));                                 \
  if (res == 2) { delete cn; Py_DECREF(PyList_GetItem(array,2)); }}
// Release shared generic array
#define RELEASESHAREDA(res, array, f, a2, a3, a4) {delete f;            \
    if (res == 1) {delete (E_Int*)a2; delete (E_Int*)a3; delete (E_Int*)a4;} \
    Py_DECREF(PyList_GetItem(array,1));                                 \
    if (res == 2) {delete (FldArrayI*)a2; Py_DECREF(PyList_GetItem(array,2));}}

#if defined E_DOUBLEREAL && defined E_DOUBLEINT
#define PYPARSETUPLE(args, format1, format2, format3, format4, ...) PyArg_ParseTuple(args, format1, __VA_ARGS__)
#elif defined E_DOUBLEREAL && !defined E_DOUBLEINT
#define PYPARSETUPLE(args, format1, format2, format3, format4, ...) PyArg_ParseTuple(args, format2, __VA_ARGS__)
#elif !defined E_DOUBLEREAL && defined E_DOUBLEINT
#define PYPARSETUPLE(args, format1, format2, format3, format4, ...) PyArg_ParseTuple(args, format3, __VA_ARGS__)
#else
#define PYPARSETUPLE(args, format1, format2, format3, format4, ...) PyArg_ParseTuple(args, format4, __VA_ARGS__)   
#endif

#if defined E_DOUBLEINT
#define PYPARSETUPLEI(args, format1, format2, ...) PyArg_ParseTuple(args, format1, __VA_ARGS__)
#else
#define PYPARSETUPLEI(args, format1, format2, ...) PyArg_ParseTuple(args, format2, __VA_ARGS__)
#endif

#if defined E_DOUBLEREAL
#define PYPARSETUPLEF(args, format1, format2, ...) PyArg_ParseTuple(args, format1, __VA_ARGS__)
#else
#define PYPARSETUPLEF(args, format1, format2, ...) PyArg_ParseTuple(args, format2, __VA_ARGS__)
#endif

namespace K_ARRAY
{
  /* Taille max de la varString pour les alloc. statiques */
  const E_Int VARSTRINGLENGTH=4096;
  /* Taille max d'un nom de variable */
  const E_Int VARNAMELENGTH=56;

  /* Retourne le nombre de variables dans varString (style "a,b,c").
     IN: varString: la chaine de variables. */
  E_Int getNumberOfVariables(char* varString);

  /* Is name present in string?
     IN: name: nom a rechercher dans string
     IN: string: la chaine recherchee.
     Retourne -1 si name n'existe pas dans string.
     Return i: si name existe et est en position i. La position etant 
     determinee suivant les virgules. */
  E_Int isNamePresent(const char* name, char* string);
  E_Int isCoordinateXPresent(char* string);
  E_Int isCoordinateYPresent(char* string);
  E_Int isCoordinateZPresent(char* string);
  E_Int isCellNatureField1Present(char* string);
  E_Int isCellNatureField2Present(char* string);
  E_Int isDensityPresent(char* string);
  E_Int isMomentumXPresent(char* string);
  E_Int isMomentumYPresent(char* string);
  E_Int isMomentumZPresent(char* string);
  E_Int isEnergyStagnationDensityPresent(char* string);
  E_Int isVelocityXPresent(char* varString);
  E_Int isVelocityYPresent(char* varString);
  E_Int isVelocityZPresent(char* varString);
  E_Int isPressurePresent(char* varString);
  E_Int isTemperaturePresent(char* varString);
  E_Int isTimePresent(char* varString);

  /* A partir de 2 varStrings, recherche si des variables communes existent et 
     retourne les vecteurs des positions des variables correspondantes dans les
     arrays (1er element : i=1).
     Les noms des variables communes sont dans varString.
     IN: varString1, varString2: les deux varStrings
     OUT: pos1: position des variables communes dans varString1
     OUT: pos2: position des variables communes dans varString2
     OUT: varString: nom des variables communes style ('a,b,c')
     Retourne:
     0 si les chaines varString1 et varString2 sont differentes
     1 si les deux chaines sont identiques.
  */
  E_Int getPosition(char* varString1, char* varString2,
                    std::vector<E_Int>& pos1, std::vector<E_Int>& pos2,
                    char* varString);
 
  /* Construit une liste des variables a partir de la chaine des variables.
     IN: varString: chaine des variables
     OUT: vars: vecteur contenant chaque chaine de variable.
     C'est la responsabilite de l'appelant de liberer la memoire des vars. */
  void extractVars(char* varString, std::vector<char*>& vars);

  /* Construit la chaine de variables varString a partir d'une liste python
     de noms de variables.
     IN: varName: liste python des noms de variables
     OUT: varString: chaine des variables (allouee par l'appelant). 
     retourne le nombre de variables. */
  E_Int getVarName(PyObject* varNames, char* varString);

  /* Analyse element string. For instance: "TRI_6" will return "TRI" and 6 */
  E_Int eltString2TypeId(char* eltString, char* eltType, E_Int& nvpe, E_Int& loc, E_Int& typeId);
  /* Analyse typeId. Retourne eltString, nvpe */  
  E_Int typeId2eltString(E_Int typeId, E_Int loc, char* eltString, E_Int& nvpe);

  /* Get data pointers dans un PyObject array.
     Pas de verification ici. */
  E_Float* getFieldPtr(PyObject* array);
  E_Int* getConnectPtr(PyObject* array);

  /* Compare 2 varStrings. 
     Retourne 0 si elles sont identiques, -1 sinon. 
     Les variables doivent etre positionnees au meme endroit */
  E_Int compareVarStrings(char* varString1, char* varString2);

  /* Ajoute un prefixe a tous les noms de variables contenu dans varString. 
     IN: varString: la chaine des variables
     IN: prefix: la chaine de prefix a ajouter aux noms des variables
     OUT: varStringOut: le resultat. */
  void addPrefixInVarString(char* varString, const char* prefix,
                            char* varStringOut);

  /* Extrait les donnees utiles d'un objet python struct array
     defini par: [ 'vars', a, ni, nj, nk ]
     ou d'un objet python unstruct array
     defini par: [ 'vars', a, c, "ELTTYPE"].
     ou ELTTYPE vaut: 
     NODE, BAR, TRI, QUAD, TETRA, PYRA, PENTA, HEXA, NGON.
     return 1: valid struct array
     f: field en stockage Fld (champs par variable)
     ni, nj, nk: number of points
     varString
     return 2: valid unstruct array
     f: field en stockage Fld (champs par variable)
     c: connectivity (champs des indices de noeuds commencants a 1)
     eltType: type of elements
     varString
     return -1: given object is not a list.
     return -2: not a valid number of elts in list.
     return -3: first element is not a var string.
     return -4: a is not a valid numpy array.
     return -5: array is structured but ni, nj, nk unvalid.
     return -6: array is unstructured but connectivity is unvalid.
     return -7: array is unstructured but elt type is unknown.
     Si shared=true, le tableau f et c partagent la memoire avec les 
     tableaux numpy. La reference sur o est incrementee.
     C'est la responsabilite de l'appelant de liberer la memoire de f et 
     eventuellement de c. Dans le cas shared=true, il faut utiliser les
     macros RELEASESHARED. */
  E_Int getFromArray(PyObject* o,
                     char*& varString,
                     FldArrayF*& f,
                     E_Int& ni, E_Int& nj, E_Int& nk,
                     FldArrayI*& c,
                     char*& eltType,
                     E_Boolean shared=false);
  E_Int getFromArray2(PyObject* o,
                      char*& varString,
                      FldArrayF*& f,
                      E_Int& ni, E_Int& nj, E_Int& nk,
                      FldArrayI*& c,
                      char*& eltType);
  E_Int getFromArray2(PyObject* o,
                      char*& varString,
                      FldArrayF*& f,
                      FldArrayI*& c,
                      char*& eltType);
  E_Int getFromArray2(PyObject* o,
                      FldArrayF*& f,
                      FldArrayI*& c);
  /* Retourne uniquement un FldArrayF (shared) sur les champs et la varstring 
     Il faut utiliser la macro RELEASESHAREDS */
  E_Int getFromArray2(PyObject* o,
                      char*& varString,
                      FldArrayF*& f);
  E_Int getFromArray2(PyObject* o, FldArrayF*& f);

  /* Extrait les donnees utiles d'un objet python struct array 
     defini par: [ 'vars', a, ni, nj, nk ]
     ou d'un objet python unstruct array
     defini par: [ 'vars', a, c, "ELTTYPE"].
     ou ELTTYPE vaut: NODE, BAR, TRI, QUAD, TETRA, PYRA, PENTA, HEXA, NGON.
     return 1: valid struct array
     f: field en stockage Dyn (champ de v1,v2,v3,...)
     ni, nj, nk: number of points
     varString
     return 2 : valid unstruct array
     f: field en stockage en stockage Dyn (champ de v1,v2,v3,...)
     c: connectivity en stockage Dyn (champ c1,c2,c3 avec indices
         commencants a zero).
     eltType: type of elements
     varString
     return -1: given object is not a list.
     return -2: not a valid number of elts in list.
     return -3: first element is not a var string.
     return -4: a is not a valid numpy array.
     return -5: array is structured but ni, nj, nk unvalid.
     return -6: array is unstructured but connectivity is unvalid.
     return -7: array is unstructured but elt type is unknown.
     C'est la responsabilite de l'appelant de liberer la memoire de f et 
     eventuellement de c. */
  E_Int getFromArray(PyObject* o,
                     char*& varString,
                     K_FLD::DynArray<E_Float>*& f,
                     E_Int& ni, E_Int& nj, E_Int& nk,
                     K_FLD::DynArray<E_Int>*& c,
                     char*& eltType);

  /* Extrait les donnees utiles d'une liste d'objets pythons.
     IN: o: liste d'objets pythons correspondant a des arrays
     OUT: res: resultat pour chaque array de arrays (1:struct, 2:unstruct, 
                                                     0: invalid, -1: skipped)
     OUT: structVarString: chaine des variables pour les arrays structures
     OUT: unstructVarString: chaine des variables pour les arrays 
                             non-structures
     OUT: structF: vecteurs des champs structures
     OUT: unstructF: vecteurs des champs non-structures
     OUT: ni,nj,nk: vecteur des dimensions pour les champ structures
     OUT: c, eltType: vecteur des connectivites et des types d'elements
                      pour les champs non-structures
     IN: skipDiffVars: selectionne les arrays ayant les memes 
                       variables que celles du premier array
     IN: skipNoCoord: selectionne seulement les array avec x,y,z
     IN: skipStructured: selectionne les arrays structures
     IN: skipUnstructured: selectionne les arrays non structures
     
     Retourne -1 si pb
     Si shared=true, tous les tableaux retournes partagent la memoire
     avec les tableaux numpy. Tous les array de o ont une reference
     incrementee.
     C'est la responsabilite de l'appelant de liberer la memoire de structF,
     unstructF et c. Dans le cas shared=true, il faut utiliser la macro
     RELEASESHARED pour chaque array de arrays ayant un res=1 ou 2.
   */
  E_Int getFromArrays(PyObject* o,
                      std::vector<E_Int>& res,
                      std::vector<char*>& structVarString,
                      std::vector<char*>& unstructVarString,
                      std::vector<FldArrayF*>& structF,
                      std::vector<FldArrayF*>& unstructF,
                      std::vector<E_Int>& ni, 
                      std::vector<E_Int>& nj, 
                      std::vector<E_Int>& nk,
                      std::vector<FldArrayI*>& c,
                      std::vector<char*>& eltType,
                      std::vector<PyObject*>& objs,
                      std::vector<PyObject*>& obju,
                      E_Boolean skipDiffVars=false,
                      E_Boolean skipNoCoord=false,
                      E_Boolean skipStructured=false,
                      E_Boolean skipUnstructured=false,
                      E_Boolean shared=false);
  /* Extrait les donnees utiles d'une liste d'objets pythons.
   IN: o: liste d'objets pythons correspondant a des arrays
   OUT: res: type de l'array dans arrays (1:struct, 2:unstruct, 0:invalid, -1: skipped). Si l'array est invalide ou skipped, il n'y pas de F, a2,...
   correspondant.
   OUT: varString: chaine des variables pour chaque array
   OUT: F: vecteurs des champs
   OUT: a2,a3,a4,: vecteur d'information pour chaque array
   Pour un array structure: a2=&ni (int*), a3=&nj (int*), a4=&nk (int*)
   Pour un array non-structure: a2=cn (FldArrayI*), a3=eltType (char*), a4=NULL
   OUT: obj: liste des objets python correspondant aux arrays
   IN: skipNoCoord: rejette les arrays sans x,y,z
   IN: skipStructured: rejette les arrays structures
   IN: skipUnstructured: rejette les arrays non structures
   IN: skipDiffVars: selectionne les arrays ayant les memes variables 
                     que celles du premier array
                     les variables sont positionnees de la meme maniere
                     pour tous les arrays
   Retourne -1 si pb
   C'est la responsabilite de l'appelant de liberer la memoire de F,
   a2,a3,a4.
  */
  E_Int getFromArrays(PyObject* o,
                      std::vector<E_Int>& res,
                      std::vector<char*>& varString,
                      std::vector<FldArrayF*>& F,
                      std::vector<void*>& a2,
                      std::vector<void*>& a3, 
                      std::vector<void*>& a4,
                      std::vector<PyObject*>& obj,
                      E_Boolean skipDiffVars,
                      E_Boolean skipNoCoord,
                      E_Boolean skipStructured,
                      E_Boolean skipUnstructured,
                      E_Boolean shared);

  /* Retourne uniquement la varString d'un array.
     varString est partage avec python.
     Retourne -1, -2 ou -3 si echec.
     IN: o: python object de l'array
     OUT: varString: chaine des variables.
  */
  E_Int getVarStringFromArray(PyObject* o, char*& varString);

  /*
    Retourne le nbre de pts (size) d'un array.
    IN: o: python object de l'array
    OUT: size: nbre de pts de l'array
    OUT: type: 1=structure, 2=non structure.
    Retourne -1 si echec.
  */
  E_Int getSizeFromArray(PyObject* o, E_Int& size, E_Int& type);

  /* 
   Retourne toutes les infos contenues dans un array,
   a l'exception des tableaux eux-memes.
   Retourne le type (1: structure, 2: non structure).
   Si type=1, ni, nj, nk sont renseignes.
   Si type=2, nvertex, nelt, sizeConnect et eltType sont renseignes.
   Cette routine ne fait pas de verification.
  */
  E_Int getInfoFromArray(PyObject* o,  char*& varString,
                         E_Int& ni, E_Int& nj, E_Int& nk,
                         E_Int& nvertex, E_Int& nelt, 
                         E_Int& sizeConect, char*& eltType);

  /* Construit un array structure a partir d'un FldArray
     IN: field: Fld champ structure
     IN: varString: variable string
     IN: ni,nj,nk: nombre de points dans field
     OUT: PyObject cree. */
  PyObject* buildArray(FldArrayF& field, const char* varString,
                       E_Int ni, E_Int nj, E_Int nk);
  /* Construit un array structure vide 
     IN: nfld: nombre de champ dans varString
     IN: varString: variable string
     IN: ni,nj,nk: nombre de points dans le champ
     OUT: PyObject cree. */
  PyObject* buildArray(E_Int nfld, const char* varString,
                       E_Int ni, E_Int nj, E_Int nk);
  /* Construit un array structure vide suivant les differentes api
     IN: nfld: nombre de champ dans varString
     IN: varString: variable string
     IN: ni,nj,nk: nombre de points dans le champ
     IN: api: 1 (array), 2 (array2)
     OUT: PyObject cree. */
  PyObject* buildArray2(E_Int nfld, const char* varString,
                        E_Int ni, E_Int nj, E_Int nk, E_Int api=1);
  PyObject* buildArray2(FldArrayF& f, const char* varString,
                        E_Int ni, E_Int nj, E_Int nk, E_Int api=1);
  /* Construit un array non structure a partir d'un FldArray
     IN: field: Fld champ non structure
     IN: varString: variable string
     IN: c: connectivite elements->noeuds (commence a 1)
     IN: et: type d'elements: 0 (NODE), 1 (BAR), 2 (TRI), 3 (QUAD)
     4 (TETRA), 5 (PYRA), 6 (PENTA), 7 (HEXA), 8 (NGON) ou star.
     IN: etString: si et=-1, utilise pour le type d'element.
     IN: center: mis a true si field est localise sur les centres,
         sinon false.
     OUT: PyObject cree. */
  PyObject* buildArray(FldArrayF& field, const char* varString,
                       FldArrayI& c, E_Int et, 
                       const char* etString=NULL, 
                       E_Boolean center=false);
  /* Construit un array non structure vide 
     IN: nfld: nombre de champs dans varString
     IN: varString: variable string
     IN: nvertex: nombre de noeuds dans le maillage
     IN: nelt: nombre d'elements dans le maillage
     IN: c: connectivite elements->noeuds (commence a 1)
     IN: et: type d'elements: 0 (NODE), 1 (BAR), 2 (TRI), 3 (QUAD),
     4 (TETRA), 5 (PYRA), 6 (PENTA), 7 (HEXA), 8 (NGON) ou star.
     IN: etString: si et=-1, utilise pour le type d'element.
     IN: center: mis a true si field est localise sur les centres,
         sinon false.
     IN: sizeConnect: taille de la connectivite, a specifier seult pour 
     les NGONS.
     OUT: PyObject cree. */
  PyObject* buildArray(E_Int nfld, const char* varString,
                       E_Int nvertex, E_Int nelt, 
                       E_Int et, const char* etString, 
                       E_Boolean center=false, E_Int sizeConnect=1); 
 PyObject* buildArray2(E_Int nfld, const char* varString,
                       E_Int nvertex, E_Int nelt, 
                       E_Int et, const char* etString, 
                       E_Boolean center=false, 
                       E_Int sizeNGon=1, E_Int sizeNFace=1,
                       E_Int nface=1, E_Int api=1);
 PyObject* buildArray2(FldArrayF& f, const char* varString, 
                       FldArrayI& cn, const char* eltType, E_Int api=1);

  /* Construit un array structure a partir d'un DynArray
     IN: field: Dyn champ structure
     IN: varString: variable string
     IN: ni,nj,nk: nombre de points dans field
     OUT: PyObject cree. */
  PyObject* buildArray(K_FLD::DynArray<E_Float>& field, const char* varString,
                       E_Int ni, E_Int nj, E_Int nk);
  
  /* Construit un array non structure a partir d'un DynArray
     IN: field: Dyn champ non structure
     IN: varString: variable string
     IN: c: Dyn connectivite elements->noeuds (commence a 0)
     IN: et: type d'elements: 0 (NODE), 1 (BAR), 2 (TRI), 3 (QUAD)
     4 (TETRA), 5 (PYRA), 6 (PENTA), 7 (HEXA) ou star.
     IN: eltString: if et=-1, utiliser pour le type d'element.
     IN: center: mis a true si field est localise sur les centres,
     sinon false.
     OUT: PyObject cree. */
  PyObject* buildArray(K_FLD::DynArray<E_Float>& field, const char* varString,
                       K_FLD::DynArray<E_Int>& c, E_Int et, 
                       const char* etString=NULL, 
                       E_Boolean center=false);

  /* Add field in array1 or array2 (in place) 
    IN: a: array
    IN: varName */
  void addFieldInArray(PyObject* a, char* varName);

  /* Extrait les donnees utiles d'un objet python o "entiers"
     Cet objet python peut etre:
     - une liste python d'entiers
     - un numpy array d'entiers
     OUT: out: FldArrayI alloue et rempli.
     Retourne 0 (FAIL), 1 (SUCCESS) */
  E_Int getFromList(PyObject* o, FldArrayI& out);

  /* Extrait les donnees utiles d'un objet python o "float"
     Cet objet python peut etre:
     - une liste python de double
     - un numpy array de double
     OUT: out: FldArrayF alloue et rempli.
     Retourne 0 (FAIL), 1 (SUCCESS) */
  E_Int getFromList(PyObject* o, FldArrayF& out);

  /* Interne */
  void cleanStructFields(std::vector<FldArrayF*>& structF);
  
  void cleanUnstrFields(std::vector<FldArrayF*>& unstrF, 
                        std::vector<FldArrayI*> cnt,
                        std::vector<char*>& eltType);
}

#undef FldArrayF
#undef FldArrayI
#endif
