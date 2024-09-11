/*    
    Copyright 2013-2024 Onera.

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

// compute stream surfaces
# include "stream.h"
# include <math.h>
# include <stdio.h>
# include <string.h>

using namespace K_FLD;
using namespace std;

#define LENGTH(dx,dy,dz) ((dx)*(dx) + (dy)*(dy) + (dz)*(dz))
//=============================================================================
/* Cree une nappe de courant a  partir d'une ligne (BAR) et d'une liste
   de grilles definies par des arrays. Les grilles d'interpolation sont celles
   qui sont structurees, contenant les infos sur la vitesse, et ayant toutes
   les memes variables dans le meme ordre.
*/
//=============================================================================
PyObject* K_POST::compStreamSurf(PyObject* self, PyObject* args)
{
  PyObject* arrays; PyObject* arrayBAR;
  PyObject* vectorNames;
  E_Int npts; E_Float signe;

  if (!PYPARSETUPLE_(args, OOO_ R_ I_,
                    &arrays, &arrayBAR, &vectorNames, &signe, &npts))
  {
      return NULL;
  }

  // Check every array in arrays
  if (PyList_Check(arrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "streamSurf: first argument must be a list.");
    return NULL;
  }
  //extraction of the 3 components of the vector used in the streamline computation 
  vector<char*> vnames;
  if (PyList_Check(vectorNames) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "streamSurf: 3 variables defined as component of the streamline vector must be defined.");
    return NULL; 
  }
  E_Int sizeOfVector = PyList_Size(vectorNames);
  if (sizeOfVector != 3 ) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "streamSurf: vector must be defined by 3 components.");
    return NULL;
  }
  for (int i = 0; i < PyList_Size(vectorNames); i++)
  {
    PyObject* tpl0 = PyList_GetItem(vectorNames, i);
    if (PyString_Check(tpl0))
    {
      char* str = PyString_AsString(tpl0);
      vnames.push_back(str);
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(tpl0)) 
    {
      char* str = (char*)PyUnicode_AsUTF8(tpl0);
      vnames.push_back(str);
    }
#endif
    else  
    {
      PyErr_SetString(PyExc_TypeError,
                      "streamSurf: vector component name must be a string.");
      return NULL;
    }
  }
  // Extract BAR
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cnBAR;
  char* eltTypeBAR; char* varStringBAR;
  E_Int res = K_ARRAY::getFromArray(arrayBAR, varStringBAR, f, im, jm, km, cnBAR, eltTypeBAR, true); 
  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(arrayBAR, f);
    PyErr_SetString(PyExc_TypeError, 
                    "streamSurf: 2nd arg is not a valid array.");
    return NULL;
  }
  if (strcmp(eltTypeBAR,"BAR")!= 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "streamSurf: 2nd arg must be a BAR array");
    RELEASESHAREDU(arrayBAR,f, cnBAR);
    return NULL;
  }
  // Extract coordinates from BAR-array
  E_Int posx = K_ARRAY::isCoordinateXPresent(varStringBAR);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varStringBAR);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varStringBAR);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "streamSurf: coordinates not found in BAR-array.");
    RELEASESHAREDU(arrayBAR,f, cnBAR);
    return NULL;
  }
  posx++; posy++; posz++;
  E_Float* xBAR = f->begin(posx);
  E_Float* yBAR = f->begin(posy);
  E_Float* zBAR = f->begin(posz);
  E_Int sizeBAR = f->getSize();

  // Extract infos from arrays
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit0; vector<E_Int> njt0; vector<E_Int> nkt0;
  vector<FldArrayI*> cnt;
  vector<char*> eltType;
  vector<PyObject*> objs0, obju0;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false; 
  E_Boolean skipDiffVars = true;
  E_Int nfld = -1;
  
  E_Int isOk = K_ARRAY::getFromArrays(arrays, resl, structVarString, unstrVarString,
                                      structF, unstrF, nit0, njt0, nkt0, 
                                      cnt, eltType, objs0, obju0, skipDiffVars,
                                      skipNoCoord, skipStructured, skipUnstructured, true);

  char* varStringOut;
  if (structVarString.size() > 0) 
  {
    varStringOut = new char [strlen(structVarString[0])+1];
    strcpy(varStringOut, structVarString[0]);
  }
  else if (unstrVarString.size() > 0) 
  {
    varStringOut = new char [strlen(unstrVarString[0])+1];
    strcpy(varStringOut, unstrVarString[0]);
  }
  else
  {
    varStringOut = new char [2];
    varStringOut[0] = '\0';
  }
  nfld = K_ARRAY::getNumberOfVariables(varStringOut);
  
  if (isOk == -1 || nfld == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "streamSurf: invalid list of arrays.");
    RELEASESHAREDU(arrayBAR,f, cnBAR);
    for (size_t nos = 0; nos < objs0.size(); nos++)
      RELEASESHAREDS(objs0[nos], structF[nos]);
    for (size_t nos = 0; nos < obju0.size(); nos++)
      RELEASESHAREDU(obju0[nos], unstrF[nos], cnt[nos]);
    return NULL;
  }
   // Build interpData 
   E_Int nzonesS = structF.size();
   E_Int nzonesU = unstrF.size();
   // InterpData structuree
   vector<E_Int> posxs1; vector<E_Int> posys1; vector<E_Int> poszs1; vector<E_Int> poscs1;
   vector<K_INTERP::InterpData*> structInterpDatas1;
   vector<FldArrayF*> structF1;
   vector<E_Int> nis1; vector<E_Int> njs1; vector<E_Int> nks1;
   vector<char*> structVarStrings1;
   E_Int isBuilt;
   for (E_Int no = 0; no < nzonesS; no++)
   {
     E_Int posx = K_ARRAY::isCoordinateXPresent(structVarString[no]); posx++;
     E_Int posy = K_ARRAY::isCoordinateYPresent(structVarString[no]); posy++;
     E_Int posz = K_ARRAY::isCoordinateZPresent(structVarString[no]); posz++;
     E_Int posc = K_ARRAY::isCellNatureField2Present(structVarString[no]); posc++;
     posxs1.push_back(posx); posys1.push_back(posy); poszs1.push_back(posz); poscs1.push_back(posc); 
     K_INTERP::InterpAdt* adt = new K_INTERP::InterpAdt(
       structF[no]->getSize(), 
       structF[no]->begin(posx),
       structF[no]->begin(posy),
       structF[no]->begin(posz),
       &nit0[no], &njt0[no], &nkt0[no], isBuilt);
     nis1.push_back(nit0[no]);
     njs1.push_back(njt0[no]);
     nks1.push_back(nkt0[no]);
     structF1.push_back(structF[no]);
     structInterpDatas1.push_back(adt);
     structVarStrings1.push_back(structVarString[no]);
   }
   // InterpData non structuree
  vector<E_Int> posxu2; vector<E_Int> posyu2; vector<E_Int> poszu2; 
  vector<E_Int> poscu2;
  vector<K_INTERP::InterpData*> unstrInterpDatas2;
  vector<FldArrayI*> cnt2;
  vector<FldArrayF*> unstrF2;
  vector<char*> unstrVarString2;
  vector<char*> eltType2;
  for (E_Int no = 0; no < nzonesU; no++)
  {
    if (strcmp(eltType[no],"TETRA")==0)
    {
      E_Int posx = K_ARRAY::isCoordinateXPresent(unstrVarString[no]); posx++;
      E_Int posy = K_ARRAY::isCoordinateYPresent(unstrVarString[no]); posy++;
      E_Int posz = K_ARRAY::isCoordinateZPresent(unstrVarString[no]); posz++;
      E_Int posc = K_ARRAY::isCellNatureField2Present(unstrVarString[no]); posc++;
      posxu2.push_back(posx); posyu2.push_back(posy); poszu2.push_back(posz); poscu2.push_back(posc); 
      K_INTERP::InterpAdt* adt = new K_INTERP::InterpAdt(unstrF[no]->getSize(), 
                                                         unstrF[no]->begin(posx),
                                                         unstrF[no]->begin(posy),
                                                         unstrF[no]->begin(posz),
                                                         cnt[no], NULL, NULL, isBuilt);
      unstrF2.push_back(unstrF[no]); cnt2.push_back(cnt[no]);
      unstrInterpDatas2.push_back(adt);
      unstrVarString2.push_back(unstrVarString[no]);
      eltType2.push_back(eltType[no]);
    }     
    else 
    {
      printf("Warning: streamSurf: unstructured element type is %s (must be TETRA).\n", eltType[no]);
      printf("Zone " SF_D_ " not taken into account to build the stream surface\n", no);
    }
  } 
  E_Int structSize =  structInterpDatas1.size();
  E_Int unstrSize = unstrInterpDatas2.size();
  E_Int interpDatasSize = structSize + unstrSize;
  
  if (interpDatasSize == 0)
  {
    RELEASESHAREDU(arrayBAR, f, cnBAR);
    PyErr_SetString(PyExc_ValueError,
                    "streamSurf: no interpData built.");
    for (size_t nos = 0; nos < objs0.size(); nos++)
      RELEASESHAREDS(objs0[nos], structF[nos]);
    for (size_t nos = 0; nos < obju0.size(); nos++)
      RELEASESHAREDU(obju0[nos], unstrF[nos], cnt[nos]);
    for (size_t nos = 0; nos < structInterpDatas1.size(); nos++)
      delete structInterpDatas1[nos];
    for (size_t nos = 0; nos < unstrInterpDatas2.size(); nos++)
      delete unstrInterpDatas2[nos];
    return NULL;
  } 
  nit0.clear(); njt0.clear(); nkt0.clear();
 
  // declaration de donnees
  // structure
  vector<FldArrayF*> structVector;
  vector<E_Int> nis; vector<E_Int> njs; vector<E_Int> nks;
  vector<E_Int> posxs; vector<E_Int> posys; vector<E_Int> poszs;
  vector<E_Int> poscs;
  vector<char*> structVarStrings;
  vector<FldArrayF*> structFields;
  vector<K_INTERP::InterpData*> structInterpDatas;
  // non structure
  vector<FldArrayF*> unstrVector;
  vector<E_Int> posxu; vector<E_Int> posyu; vector<E_Int> poszu; 
  vector<E_Int> poscu;
  vector<char*> unstrVarStrings;
  vector<FldArrayF*> unstrFields;
  vector<FldArrayI*> connectu;
  vector<char*> eltTypes;
  vector<K_INTERP::InterpData*> unstrInterpDatas;


  // seuls sont pris en compte les champs ayant les variables du vecteur
  // ts les arrays traites doivent avoir le meme nb de champs
  if (structSize > 0) 
  {
    E_Int found = extractVectorFromStructArrays(signe, nis1, njs1, nks1,
                                                posxs1, posys1, poszs1, poscs1,
                                                structVarStrings1, structF1, 
                                                structInterpDatas1,
                                                nis, njs, nks, 
                                                posxs, posys, poszs, poscs,
                                                structVarStrings, structFields,
                                                structInterpDatas,
                                                structVector, vnames);
    if (found != 1)
    {
      RELEASESHAREDU(arrayBAR,f, cnBAR);
      for (size_t nos = 0; nos < objs0.size(); nos++)
        RELEASESHAREDS(objs0[nos], structF[nos]);
      for (size_t nos = 0; nos < obju0.size(); nos++)
        RELEASESHAREDU(obju0[nos], unstrF[nos], cnt[nos]);
      for (size_t nos = 0; nos < structInterpDatas1.size(); nos++)
        delete structInterpDatas1[nos];
      for (size_t nos = 0; nos < unstrInterpDatas2.size(); nos++)
        delete unstrInterpDatas2[nos];
      
      if (found == -1) 
        PyErr_SetString(PyExc_ValueError,
                        "streamSurf: no field corresponding to vector component.");
      else // found = -2 uniquement pour l instant
        PyErr_SetString(PyExc_ValueError,
                        "streamSurf: only TETRA type is valid for unstructured mesh.");
      return NULL;
    }
  }
  // extract des variables du vecteur sur les maillages non structures 
  if (unstrSize > 0) 
  {
    E_Int found = extractVectorFromUnstrArrays(signe, posxu2, posyu2, poszu2, poscu2,
                                               unstrVarString2, unstrF2, cnt2,
                                               eltType2, unstrInterpDatas2,
                                               posxu, posyu, poszu, poscu,
                                               unstrVarStrings, unstrFields, connectu,
                                               eltTypes, unstrInterpDatas,
                                               unstrVector, vnames); 
    if (found != 1)
    {
      RELEASESHAREDU(arrayBAR, f, cnBAR);
      for (size_t nos = 0; nos < objs0.size(); nos++)
        RELEASESHAREDS(objs0[nos], structF[nos]);
      for (size_t nos = 0; nos < obju0.size(); nos++)
        RELEASESHAREDU(obju0[nos], unstrF[nos], cnt[nos]);
      for (size_t nos = 0; nos < structInterpDatas1.size(); nos++)
        delete structInterpDatas1[nos];
      for (size_t nos = 0; nos < unstrInterpDatas2.size(); nos++)
        delete unstrInterpDatas2[nos];
      if (found == -1) 
        PyErr_SetString(PyExc_ValueError,
                        "streamSurf: no field corresponding to vector component.");
      else // found = -2 uniquement pour l instant
        PyErr_SetString(PyExc_ValueError,
                        "streamSurf: only TETRA type is valid for unstructured mesh.");
      return NULL;
    }
  }

  // Petit nettoyage intermediaire
  nis1.clear(); njs1.clear(); nks1.clear();
  posxs1.clear(); posys1.clear(); poszs1.clear(); poscs1.clear();
  posxu2.clear(); posyu2.clear(); poszu2.clear(); poscu2.clear();

  // Construction du front
  // et determination des extremites (tracer tleft et tright)
  vector<tracer*> front;
  tracer *tleft = NULL;
  tracer *tright = NULL;
  createFront(xBAR, yBAR, zBAR, *cnBAR, front, sizeBAR);

  // Affectation des traceurs tleft et tright
  for (size_t v = 0; v < front.size(); v++)
  {
    if (front[v]->left == NULL) tleft = front[v];
    if (front[v]->right == NULL) tright = front[v];
  }
  if (tleft == NULL || tright == NULL) 
  { 
    size_t last = front.size()-1;
    tleft = front[0]; tright = front[last]; 
    front[0]->left = NULL; front[last]->right = NULL; 
  }
  
  // Creation des tableaux (champ et connectivite) pour la nappe de courant
  FldArrayF* field = new FldArrayF(6*sizeBAR*20*npts,nfld);
  FldArrayI* cn = new FldArrayI(4*sizeBAR*20*npts,3);

  // Compteur sur le nombre de triangles constituant la nappe de courant
  E_Int nt = 0;
  // Avance du front
  for (E_Int n = 0; n < npts; n++)
    advanceFront(front, tleft, tright, npts, nt, field, cn,
                 structInterpDatas, structFields, structVector,
                 nis, njs, nks, posxs, posys, poszs, poscs, 
                 unstrInterpDatas, unstrFields, unstrVector,
                 connectu, posxu, posyu, poszu, poscu);

  // Reallocation des tableaux (champ et connectivite) pour la nappe de courant
  field->reAllocMat(3*nt, nfld);
  cn->reAllocMat(nt,3);

  //little cleaning
  E_Int structVectorSize = structVector.size();
  for (E_Int v = 0; v < structVectorSize; v++)
    delete structVector[v];
  E_Int unstrVectorSize = unstrVector.size();
  for (E_Int v = 0; v < unstrVectorSize; v++)
    delete unstrVector[v];

  // Build array
  PyObject* tpl = K_ARRAY::buildArray(*field, varStringOut, *cn , -1 , "TRI");
  delete [] varStringOut; delete field; delete cn;
  
  for (size_t i = 0; i < front.size(); i++) { delete front[i]; }
  front.clear();

  RELEASESHAREDU(arrayBAR,f, cnBAR);
  for (size_t nos = 0; nos < objs0.size(); nos++)
    RELEASESHAREDS(objs0[nos], structF[nos]);
  for (size_t nos = 0; nos < obju0.size(); nos++)
    RELEASESHAREDU(obju0[nos], unstrF[nos], cnt[nos]);
  for (size_t nos = 0; nos < structInterpDatas1.size(); nos++)
    delete structInterpDatas1[nos];
  for (size_t nos = 0; nos < unstrInterpDatas2.size(); nos++)
    delete unstrInterpDatas2[nos];
      
  return tpl;
}
//=============================================================================
// Creation du front a partir du BAR-array
// Un front est defini par les traceurs gauche et droit des "ribbons"
//=============================================================================
void K_POST::createFront(E_Float* xBAR, E_Float* yBAR, E_Float* zBAR, 
                         FldArrayI& cnBAR, vector<tracer*>& front, E_Int npts)
{

  // Initialisation d'un front de taille npts (taille du BAR-array)
  for (E_Int i = 0; i < npts; i++)
  {
    tracer* t = new tracer;
    t->x = xBAR[i]; 
    t->y = yBAR[i];
    t->z = zBAR[i]; 
    t->left = NULL;
    t->right = NULL;
    front.push_back(t);
  }

  // Remplissage du front avec les elements du BAR-array
  E_Int ne = cnBAR.getSize(); // nbre d'elements
  for (E_Int e = 0; e < ne; e++)
  {
    E_Int i1 = cnBAR[e]-1;
    E_Int i2 = cnBAR[e];
    front[i1]->right = front[i2];
    front[i2]->left  = front[i1];    
  }
}

//=============================================================================
// Progression du front
//=============================================================================
void K_POST::advanceFront(
  vector<tracer*> front, tracer* tleft, tracer* tright,
  E_Int npts, E_Int& nt, FldArrayF* field, FldArrayI* cn,
  vector<K_INTERP::InterpData*>& listOfStructInterpData, 
  vector<FldArrayF*>& listOfStructFields,
  vector<FldArrayF*>& listOfStructVelocities,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks, 
  vector<E_Int>& posxs, vector<E_Int>& posys, vector<E_Int>& poszs, 
  vector<E_Int>& poscs,
  vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
  vector<FldArrayF*>& listOfUnstrFields,
  vector<FldArrayF*>& listOfUnstrVelocities,
  vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, vector<E_Int>& poszu, 
  vector<E_Int>& poscu)
{
  // Progression du front de la gauche vers la droite 
 advanceRibbonLeft(tleft, npts, nt, field, cn,
                    listOfStructInterpData, listOfStructFields, 
                    listOfStructVelocities, nis, njs, nks, 
                    posxs, posys, poszs, poscs,
                    listOfUnstrInterpData, listOfUnstrFields,
                    listOfUnstrVelocities, connectu,
                    posxu, posyu, poszu, poscu);

//   // Progression du front de la droite vers la gauche 
//   advanceRibbonRight(tright, npts, nt, field, cn,
//                      listOfStructInterpData, listOfStructFields, 
//                      listOfStructVelocities, nis, njs, nks, 
//                      posxs, posys, poszs, poscs,
//                      listOfUnstrInterpData, listOfUnstrFields,
//                      listOfUnstrVelocities, connectu,
//                      posxu, posyu, poszu, poscu);
}

//=============================================================================
// Avance des "ribbons" de gauche a droite
//=============================================================================
void K_POST::advanceRibbonLeft(
  tracer* t, E_Int npts, E_Int& nt, FldArrayF* field, FldArrayI* cn,
  vector<K_INTERP::InterpData*>& listOfStructInterpData, 
  vector<FldArrayF*>& listOfStructFields,
  vector<FldArrayF*>& listOfStructVelocities,
  vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks, 
  vector<E_Int>& posxs, vector<E_Int>& posys, vector<E_Int>& poszs, 
  vector<E_Int>& poscs,
  vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
  vector<FldArrayF*>& listOfUnstrFields,
  vector<FldArrayF*>& listOfUnstrVelocities,
  vector<FldArrayI*>& connectu,
  vector<E_Int>& posxu, vector<E_Int>& posyu, vector<E_Int>& poszu, 
  vector<E_Int>& poscu)
{

  // Nombre de champs
  E_Int nfld = field->getNfld();
  FldArrayF fieldL(2, nfld);
  FldArrayF fieldR(2, nfld);

  // Booleen indiquant si le rattrapage est effectue ou non 
  E_Boolean caught_up = false;
  
  // initialisation prev_diag
  E_Float prev_diag = 1e6;
    
  // Donnees pour l'appel de la methode computeStreamLineElts
  FldArrayI* cnSurf = NULL;
  E_Float* xSurf = NULL;
  E_Float* ySurf = NULL;
  E_Float* zSurf = NULL;
  E_Int sizeSurf = 0;

  // Boucle recursive sur les "ribbon"
  while(1)
  {
    if (t->right == NULL) return;

    // Coordonnees des traceurs  (points L0 et R0)
    E_Float xL0 = t->x;
    E_Float yL0 = t->y;
    E_Float zL0 = t->z;

    E_Float xR0=  t->right->x;
    E_Float yR0=  t->right->y;
    E_Float zR0=  t->right->z;

    // Coordonnees a l instant (n+1) (points L1 et R1)
    E_Float xL1, yL1, zL1;
    E_Float xR1, yR1, zR1;

    // isinterp = 1 si le point est interpolable, 0 sinon
    E_Int isinterp = 1;
    isinterp = computeStreamLineElts(
      xL0, yL0, zL0, 
      listOfStructInterpData, listOfStructFields, 
      listOfStructVelocities, nis, njs, nks, 
      posxs, posys, poszs, poscs,
      listOfUnstrInterpData, listOfUnstrFields,
      listOfUnstrVelocities, connectu,
      posxu, posyu, poszu, poscu,
      *cnSurf, xSurf, ySurf, zSurf, sizeSurf,
      fieldL);
    if (!isinterp) return;
  
    isinterp = computeStreamLineElts(
      xR0, yR0, zR0, 
      listOfStructInterpData, listOfStructFields, 
      listOfStructVelocities, nis, njs, nks, 
      posxs, posys, poszs, poscs,
      listOfUnstrInterpData, listOfUnstrFields,
      listOfUnstrVelocities, connectu,
      posxu, posyu, poszu, poscu,
      *cnSurf, xSurf, ySurf, zSurf, sizeSurf,
      fieldR);
    if (!isinterp) return;
 
    // Affectation des points L1 et R1
    xL1 = fieldL(1,1);
    yL1 = fieldL(1,2);
    zL1 = fieldL(1,3);
    xR1 = fieldR(1,1);
    yR1 = fieldR(1,2);
    zR1 = fieldR(1,3);

    // Longueur des diagonales L1R0 et L0R1 du quadrilatere 
    E_Float left_diag = LENGTH(xL1-xR0,yL1-yR0,zL1-zR0);
    E_Float right_diag = LENGTH(xL0-xR1,yL0-yR1,zL0-zR1);

    // min_diag : plus petite diagonale 
    E_Float min_diag = K_FUNC::E_min(left_diag,right_diag);


//     // RAFFINEMENT
//     // -----------------------------------------------------
//     // Test pour ajouter ou non un point (largeur / hauteur > 2)
//     E_Float ratioLeft  = LENGTH(xL0-xR0,yL0-yR0,zL0-zR0)/LENGTH(xL0-xL1,yL0-yL1,zL0-zL1);
//     E_Float ratioRight = LENGTH(xL0-xR0,yL0-yR0,zL0-zR0)/LENGTH(xR0-xR1,yR0-yR1,zR0-zR1);

//     if ((ratioLeft > 2)&& (ratioRight > 2))
//     {
//       // ajout d'un nouveau traceur
//       tracer* tnew = new tracer;
//       tnew->x    = (xL0+xR0)/2.; 
//       tnew->y    = (yL0+yR0)/2.; 
//       tnew->z    = (zL0+zR0)/2.;
//       tnew->left = t;
//       tnew->right = t->right;
//       xR0 = tnew->x;
//       yR0 = tnew->y;
//       zR0 = tnew->z;
//       isinterp = computeStreamLineElts(xR0, yR0, zR0, 
//                                        listOfStructInterpData, listOfStructFields, 
//                                        listOfStructVelocities, nis, njs, nks, 
//                                        posxs, posys, poszs, poscs,
//                                        listOfUnstrInterpData, listOfUnstrFields,
//                                        listOfUnstrVelocities, connectu,
//                                        posxu, posyu, poszu, poscu,
//                                        *cnSurf, xSurf, ySurf, zSurf, sizeSurf,
//                                        fieldR);
//       if (!isinterp) return;
      
//       xR1 = fieldR(1,1);
//       yR1 = fieldR(1,2);
//       zR1 = fieldR(1,3);
      
//       // modifications des voisins des traceurs t et t->right
//       tracer* temp = t->right;
//       t->right = tnew;
//       temp->left = tnew;

//       // Recalcul des longueur des diagonales L1R0 et L0R1 du quadrilatere 
//       left_diag = LENGTH(xL1-xR0,yL1-yR0,zL1-zR0);
//       right_diag = LENGTH(xL0-xR1,yL0-yR1,zL0-zR1);
      
//       // Recalcul de min_diag : plus petite diagonale 
//       min_diag = K_FUNC::E_min(left_diag,right_diag);

//       // Recalcul de prev_diag
//       //prev_diag = prev_diag/2.;
//       prev_diag = 1e6;
//     }

//     // DERAFFINEMENT
//     // -----------------------------------------------------
//     // Test pour supprimer un point (largeur / hauteur < 0.5)
//     else if ((ratioLeft < 0.5)&& (ratioRight < 0.5))
//     {
//       if (t->right->right != NULL) // points du ribbon de droite necessaires 
//       {
//         // Definition de RR0 : point a droite de R0
//         E_Float xRR0=  t->right->right->x;
//         E_Float yRR0=  t->right->right->y;
//         E_Float zRR0=  t->right->right->z;
//         // Calcul du points RR1 
//         FldArrayF fieldRR(2, nfld);
//         isinterp = computeStreamLineElts(xRR0, yRR0, zRR0, 
//                                          listOfStructInterpData, listOfStructFields, 
//                                          listOfStructVelocities, nis, njs, nks, 
//                                          posxs, posys, poszs, poscs,
//                                          listOfUnstrInterpData, listOfUnstrFields,
//                                          listOfUnstrVelocities, connectu,
//                                          posxu, posyu, poszu, poscu,
//                                          *cnSurf, xSurf, ySurf, zSurf, sizeSurf,
//                                          fieldRR);
//         E_Float xRR1 = fieldRR(1,1);
//         E_Float yRR1 = fieldRR(1,2);
//         E_Float zRR1 = fieldRR(1,3);
//         // Ecriture de 3 triangles :
//         // (L0, R0 , L1) 
//         writeTriangle(*field,*cn,fieldL,fieldR,fieldL,nt);
//         // (R0, RR0, R1)
//         writeTriangle(*field,*cn,fieldR,fieldRR,fieldRR,nt);
//         // (L1, R0 , RR1)
//         fieldL(0,1) = xL1; 
//         fieldL(0,2) = yL1; 
//         fieldL(0,3) = zL1; 
//         writeTriangle(*field,*cn,fieldL,fieldR,fieldRR,nt);

//         // L1 devient L0
//         xL0 = xL1;
//         yL0 = yL1;
//         zL0 = zL1;
//         // RR1 devient R0
//         xR0 = xRR1;
//         yR0 = yRR1;
//         zR0 = zRR1;
//         fieldR(0,1) = xRR1;
//         fieldR(0,2) = yRR1;
//         fieldR(0,3) = zRR1;
//         // Calcul de L1 et R1
//         isinterp = computeStreamLineElts(xL0, yL0, zL0, 
//                                          listOfStructInterpData, listOfStructFields, 
//                                          listOfStructVelocities, nis, njs, nks, 
//                                          posxs, posys, poszs, poscs,
//                                          listOfUnstrInterpData, listOfUnstrFields,
//                                          listOfUnstrVelocities, connectu,
//                                          posxu, posyu, poszu, poscu,
//                                          *cnSurf, xSurf, ySurf, zSurf, sizeSurf,
//                                          fieldL);
//         if (!isinterp) return;
        
//         xL1 = fieldL(1,1);
//         yL1 = fieldL(1,2);
//         zL1 = fieldL(1,3);
//         isinterp = computeStreamLineElts(xR0, yR0, zR0, 
//                                          listOfStructInterpData, listOfStructFields, 
//                                          listOfStructVelocities, nis, njs, nks, 
//                                          posxs, posys, poszs, poscs,
//                                          listOfUnstrInterpData, listOfUnstrFields,
//                                          listOfUnstrVelocities, connectu,
//                                          posxu, posyu, poszu, poscu,
//                                          *cnSurf, xSurf, ySurf, zSurf, sizeSurf,
//                                          fieldR);
//         if (!isinterp) return;
        
//         xR1 = fieldR(1,1);
//         yR1 = fieldR(1,2);
//         zR1 = fieldR(1,3);

//         // Gestion des traceurs
//         // --------------------
//         //  gestion de t->right
//         tracer* tr = t->right;
//         tr->x =  xRR1;
//         tr->y =  yRR1;
//         tr->z =  zRR1;
//         tr->left = t;
//         tr->right = t->right->right;
//         //delete t->right;
//         // gestion de t 
//         t->right = tr; 
//         // gestion de t ->right->right
//         tracer* trr = t->right->right;
//         trr->left = t->right;

//         // Mise a jour des coordonnees des nouveaux traceurs
//         t->x = xL0; t->y = yL0; t->z = zL0;

//         // Recalcul des longueur des diagonales L1R0 et L0R1 du quadrilatere 
//         left_diag = LENGTH(xL1-xR0,yL1-yR0,zL1-zR0);
//         right_diag = LENGTH(xL0-xR1,yL0-yR1,zL0-zR1);
        
//         // Recalcul de min_diag : plus petite diagonale 
//         min_diag = K_FUNC::E_min(left_diag,right_diag);
        
//         // Recalcul de prev_diag
//         //prev_diag = prev_diag/2.;
//         prev_diag = 1e6;
//       }
//     }

    // Entier indiquant si on avance a gauche (valeur 1) ou non (valeur 0)
    E_Int advance_on_left = (K_FUNC::E_abs(min_diag - left_diag) < K_CONST::E_GEOM_CUTOFF) ? 1 : 0;

    // Conditions de sortie de la boucle :
    // rattrapage effectue
    // avance a gauche et diagonale droite > diagonale precedente
    if (caught_up && ((advance_on_left)||(right_diag > prev_diag)))
      return;
    
    if (advance_on_left)
    {
      // Ecriture du triangle (coordonnees et connectivites)
      writeTriangle(*field,*cn,fieldL,fieldR,fieldL,nt);

      // Mise a jour du traceur
      // On remplace L0 par L1
      t->x = xL1; t->y = yL1; t->z = zL1;

      caught_up = true;
    }
    else
    {
      // Ecriture du triangle (coordonnees et connectivites)
      writeTriangle(*field,*cn,fieldL,fieldR,fieldR,nt);
      // avance du "ribbon" voisin de droite du "ribbon" actuel
      if (t->right != NULL)
      {
        advanceRibbonLeft(t->right, npts, nt, field, cn,
                          listOfStructInterpData, listOfStructFields, 
                          listOfStructVelocities, nis, njs, nks, 
                          posxs, posys, poszs, poscs,
                          listOfUnstrInterpData, listOfUnstrFields,
                          listOfUnstrVelocities, connectu,
                          posxu, posyu, poszu, poscu);
        // Mise a jour du traceur
        // On remplace R0 par R1
        t->right->x = xR1; t->right->y = yR1; t->right->z = zR1;
      }
      else
        return;
    }

    // Mise a jour de prev_diag
    prev_diag = min_diag;
  }
}

//=============================================================================
// Avance des "ribbons" de droite a gauche
//=============================================================================
void K_POST::advanceRibbonRight(tracer* t, E_Int npts, E_Int& nt, FldArrayF* field, FldArrayI* cn,
                                vector<K_INTERP::InterpData*>& listOfStructInterpData, 
                                vector<FldArrayF*>& listOfStructFields,
                                vector<FldArrayF*>& listOfStructVelocities,
                                vector<E_Int>& nis, vector<E_Int>& njs, vector<E_Int>& nks, 
                                vector<E_Int>& posxs, vector<E_Int>& posys, vector<E_Int>& poszs, 
                                vector<E_Int>& poscs,
                                vector<K_INTERP::InterpData*>& listOfUnstrInterpData, 
                                vector<FldArrayF*>& listOfUnstrFields,
                                vector<FldArrayF*>& listOfUnstrVelocities,
                                vector<FldArrayI*>& connectu,
                                vector<E_Int>& posxu, vector<E_Int>& posyu, vector<E_Int>& poszu, 
                                vector<E_Int>& poscu)
{
  // Nombre de champs
  E_Int nfld = field->getNfld();
  FldArrayF fieldL(2, nfld);
  FldArrayF fieldR(2, nfld);

  // Booleen indiquant si le rattrapage est effectue ou non 
  E_Boolean caught_up = false;
    
  // Donnees pour l'appel de la methode computeStreamLineElts
  FldArrayI* cnSurf = NULL;
  E_Float* xSurf = NULL;
  E_Float* ySurf = NULL;
  E_Float* zSurf = NULL;
  E_Int sizeSurf = 0;
  
  // DBG : quelle initialisation ??? 
  E_Float prev_diag = 1e6;

  // Boucle recursive sur les "ribbons"
  while(1)
  {

    if (t->left == NULL) return;

    // Coordonnees des traceurs  (points L0 et R0)
    E_Float xL0 = t->left->x;
    E_Float yL0 = t->left->y;
    E_Float zL0 = t->left->z;
    E_Float xR0=  t->x;
    E_Float yR0=  t->y;
    E_Float zR0=  t->z;
    // Coordonnees a l instant (n+1) (points L1 et R1)
    E_Float xL1, yL1, zL1;
    E_Float xR1, yR1, zR1;

    // isinterp = 1 si le point est interpolable, 0 sinon
    E_Int isinterp = 1;
    isinterp = computeStreamLineElts(xL0, yL0, zL0, 
                                     listOfStructInterpData, listOfStructFields, 
                                     listOfStructVelocities, nis, njs, nks, 
                                     posxs, posys, poszs, poscs,
                                     listOfUnstrInterpData, listOfUnstrFields,
                                     listOfUnstrVelocities, connectu,
                                     posxu, posyu, poszu, poscu,
                                     *cnSurf, xSurf, ySurf, zSurf, sizeSurf,
                                     fieldL);
    if (!isinterp) return;

    xL1 = fieldL(1,1);
    yL1 = fieldL(1,2);
    zL1 = fieldL(1,3);
    isinterp = computeStreamLineElts(xR0, yR0, zR0, 
                                     listOfStructInterpData, listOfStructFields, 
                                     listOfStructVelocities, nis, njs, nks, 
                                     posxs, posys, poszs, poscs,
                                     listOfUnstrInterpData, listOfUnstrFields,
                                     listOfUnstrVelocities, connectu,
                                     posxu, posyu, poszu, poscu,
                                     *cnSurf, xSurf, ySurf, zSurf, sizeSurf,
                                     fieldR);
    if (!isinterp) return;
    xR1 = fieldR(1,1);
    yR1 = fieldR(1,2);
    zR1 = fieldR(1,3);

    // Longueur des diagonales L1R0 et L0R1 du quadrilatere 
    E_Float left_diag = LENGTH(xL1-xR0,yL1-yR0,zL1-zR0);
    E_Float right_diag = LENGTH(xL0-xR1,yL0-yR1,zL0-zR1);

    // min_diag : plus petite diagonale 
    E_Float min_diag = K_FUNC::E_min(left_diag,right_diag);

    // Test pour ajouter ou non un point (largeur / hauteur > 2)
    // Entier indiquant si on avance a droite (valeur 1) ou non (valeur 0)
    E_Float ratioLeft  = LENGTH(xL0-xR0,yL0-yR0,zL0-zR0)/LENGTH(xL0-xL1,yL0-yL1,zL0-zL1);
    E_Float ratioRight = LENGTH(xL0-xR0,yL0-yR0,zL0-zR0)/LENGTH(xR0-xR1,yR0-yR1,zR0-zR1);
    if ((ratioLeft > 2)&& (ratioRight > 2))
    {
    // ajout d'un nouveau traceur
      tracer* tnew = new tracer;
      tnew->x    = (xL0+xR0)/2.; 
      tnew->y    = (yL0+yR0)/2.; 
      tnew->z    = (zL0+zR0)/2.;
      tnew->right = t;
      tnew->left = t->left;
      xL0 = tnew->x;
      yL0 = tnew->y;
      zL0 = tnew->z;
      isinterp = computeStreamLineElts(xL0, yL0, zL0, 
                                       listOfStructInterpData, listOfStructFields, 
                                       listOfStructVelocities, nis, njs, nks, 
                                       posxs, posys, poszs, poscs,
                                       listOfUnstrInterpData, listOfUnstrFields,
                                       listOfUnstrVelocities, connectu,
                                       posxu, posyu, poszu, poscu,
                                       *cnSurf, xSurf, ySurf, zSurf, sizeSurf,
                                       fieldL);
      if (!isinterp) return;
      
      xL1 = fieldL(1,1);
      yL1 = fieldL(1,2);
      zL1 = fieldL(1,3);
      
      // modifications des voisins des traceurs t et t->left
      tracer* temp = t->left;
      t->left = tnew;
      temp->right = tnew;

      // Recalcul des longueur des diagonales L1R0 et L0R1 du quadrilatere 
      left_diag = LENGTH(xL1-xR0,yL1-yR0,zL1-zR0);
      right_diag = LENGTH(xL0-xR1,yL0-yR1,zL0-zR1);
      
      // Recalcul de min_diag : plus petite diagonale 
      min_diag = K_FUNC::E_min(left_diag,right_diag);

      // Recalcul de prev_diag
      //prev_diag = prev_diag/2.;
      prev_diag = 1e6;
    }

    E_Int advance_on_right = (K_FUNC::E_abs(min_diag - right_diag) < K_CONST::E_GEOM_CUTOFF) ? 1 : 0;

    // Conditions de sortie de la boucle :
    // rattrapage effectue
    // avance a droite et diagonale gauche > diagonale precedente
    if (caught_up && ((advance_on_right)||(left_diag > prev_diag)))
      return;

    if (advance_on_right)
    {
      // Ecriture du triangle (coordonnees et connectivites)
      writeTriangle(*field,*cn,fieldL,fieldR,fieldR,nt);
      // Mise a jour du traceur
      // On remplace R0 par R1
      t->x = xR1; t->y = yR1; t->z = zR1;
      caught_up = true;
    }
    else
    {
      // Ecriture du triangle (coordonnees et connectivites)
      writeTriangle(*field,*cn,fieldL,fieldR,fieldL,nt);
      // avance du "ribbon" voisin de gauche du "ribbon" actuel
      if (t->left != NULL)
      {
        advanceRibbonRight(t->left, npts, nt, field, cn,
                           listOfStructInterpData, listOfStructFields, 
                           listOfStructVelocities, nis, njs, nks, 
                           posxs, posys, poszs, poscs,
                           listOfUnstrInterpData, listOfUnstrFields,
                           listOfUnstrVelocities, connectu,
                           posxu, posyu, poszu, poscu);
        // Mise a jour du traceur
        // On remplace L0 par L1
        t->left->x = xL1; t->left->y = yL1; t->left->z = zL1;
      }
      else
        return;
    }

    // Mise a jour de prev_diag
    prev_diag = min_diag;
  }
}

//=============================================================================
// Ecriture d'un triangle ABC (champs et connectivites) 
// constituant la nappe de courant
//=============================================================================
void K_POST::writeTriangle(FldArrayF& field, FldArrayI& cn, 
                           FldArrayF& fL0, FldArrayF& fR0, FldArrayF& fLR1,
                           E_Int& nt)
{
  // Reallocation des champs si necessaire
  if (field.getSize() <= 3*nt+2) field.reAllocMat(3*nt+2 + 1000,field.getNfld());
  if (cn.getSize() <= 3*nt+2) cn.reAllocMat(3*nt+2 + 1000,3);
  // Affectation des champs
  for (E_Int n = 1; n <= field.getNfld(); n++)
  {
    field(3*nt  ,n) = fL0(0,n);
    field(3*nt+1,n) = fR0(0,n);
    field(3*nt+2,n) = fLR1(1,n);
  }
  // Affectation des connectivites
  cn(nt,1)= 3*nt+1;
  cn(nt,2)= 3*nt+2;
  cn(nt,3)= 3*nt+3;

  // incrementation du compteur sur les triangles
  nt++;
}
