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

# include "transform.h"
# include <list>
# include <stdlib.h>
# include "Nuga/include/KdTree.h"
# include "Nuga/include/ArrayAccessor.h"

using namespace std;
using namespace K_FLD;
namespace K_TRANSFORM
{
  struct Block;
  struct BlockEdge
  {
      E_Int _grade;// grade de la facette selon la methode Weakest Descent
      E_Int _externEdge;// la face a-t-elle un pt exterieur
      Block* _blk1;
      Block* _blk2;
      E_Int _dirBlk1;
      E_Int _dirBlk2;
      E_Float _ptA[3];
      E_Float _ptB[3];
      E_Float _ptC[3];
      E_Float _ptD[3];
  };

  struct Block
  {
      E_Int _nmerging;
      list<BlockEdge*> _listOfEdges;
      list<Block*> _voisins;
      E_Int _ni;
      E_Int _nj;
      E_Int _nk;
      FldArrayF* _field;
      FldArrayF* _fieldc;
  };
  E_Int matchingEdges2D(E_Float xA1, E_Float yA1, E_Float zA1,
                        E_Float xB1, E_Float yB1, E_Float zB1,
                        E_Float xA2, E_Float yA2, E_Float zA2,
                        E_Float xB2, E_Float yB2, E_Float zB2,
                        E_Float tol=1.e-10);

  E_Int matchingEdges3D(E_Float xA1, E_Float yA1, E_Float zA1,
                        E_Float xB1, E_Float yB1, E_Float zB1,
                        E_Float xC1, E_Float yC1, E_Float zC1,
                        E_Float xD1, E_Float yD1, E_Float zD1,
                        E_Float xA2, E_Float yA2, E_Float zA2,
                        E_Float xB2, E_Float yB2, E_Float zB2,
                        E_Float xC2, E_Float yC2, E_Float zC2,
                        E_Float xD2, E_Float yD2, E_Float zD2,
                        E_Float tol=1.e-10);

  void buildListOfBlks(
    vector<E_Int>& nit, vector<E_Int>& njt, vector<E_Int>& nkt,
    vector<FldArrayF*>& coords, vector<FldArrayF*>& fieldsc,
    list<Block*>& listOfBlks);

  void buildListOfMergeablesEdges(
    E_Int dim, E_Float alphaRef, E_Int posx, E_Int posy, E_Int posz,
    Block* blk1, Block* blk2,
    list<BlockEdge*>& listOfMergeableEdges,
    E_Float tol=1.e-10);
  void buildListOfMergeablesEdges2D(
    E_Float alphaRef,  E_Int posx, E_Int posy, E_Int posz,
    Block* blk1, Block* blk2,
    list<BlockEdge*>& listOfMergeableEdges,
    E_Float tol=1.e-10);
  void buildListOfMergeablesEdges3D(
    E_Int posx, E_Int posy, E_Int posz,
    Block* blk1, Block* blk2,
    list<BlockEdge*>& listOfMergeableEdges,
    E_Float tol=1.e-10);

  void gradeEdge(E_Int dim, BlockEdge* edge, list<Block*>& listOfBlks,
                 list<BlockEdge*>& listOfMergeableEdges,
                 E_Float tol=1.e-10);

  void gradeEdge2D(BlockEdge* edge, list<Block*>& listOfBlks,
                   list<BlockEdge*>& listOfMergeableEdges,
                   E_Float tol=1.e-10);

  void gradeEdge3D(BlockEdge* edge, list<Block*>& listOfBlks,
                   list<BlockEdge*>& listOfMergeableEdges,
                   E_Float tol=1.e-10);

  void deleteListOfEdgesForBlk(Block* block1, BlockEdge* faceMin,
                               list<BlockEdge*>& listOfMergeableEdges);

  void deleteBlkInfo(Block* blk, list<Block*>& listOfBlks);
}

//=============================================================================
/* Fusion des grilles surfaciques 2D (k=1) ou volumiques selon l'algorithme
   de Rigby
   Ne suppose pour l'instant que le maillage est coincident (pas de raccords
   1 pt sur 2 ) a cause de la connectivite
   Attention : pas de partage de memoire entre python et le c car les tableaux
   sont modifies/supprimes */
//=============================================================================
PyObject* K_TRANSFORM::mergeStructGrids(PyObject* self, PyObject* args)
{
  E_Float tol=1.e-10;// tolerance match
  E_Float alphaRef=180.;
  E_Int dircons=0;
  PyObject *arrays, *arraysc;
  E_Int sizeMax;

  if (!PYPARSETUPLE_(args, OO_ II_ RR_,
                    &arrays, &arraysc, &sizeMax, &dircons, &tol, &alphaRef))
  {
      return NULL;
  }
  // Check every arrays
  if (PyList_Check(arrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "merge: arrays argument must be a list.");
    return NULL;
  }
  // Check every arrays
  if (PyList_Check(arraysc) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "merge: arraysc argument must be a list.");
    return NULL;
  }

  // Extract infos from arrays
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF;
  vector<FldArrayF*> unstrF;
  vector<E_Int> nit;  vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltType;
  vector<PyObject*> objst, objut;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = false;
  E_Bool skipUnstructured = true;
  E_Bool skipDiffVars = true;
  vector<E_Int> rest;
  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, rest, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltType, objst, objut,
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "merge: arrays is not valid.");
    return NULL;
  }
  E_Int api = -1;
  E_Int nzones = structF.size();
  E_Int dim = 2;
  for (E_Int v = 0; v < nzones; v++)
  {
    if (nkt[v] != 1) { dim = 3; break; }
  }

   // Extract infos from arrays defining fields at centers
  vector<char*> structVarStringc; vector<char*> unstrVarStringc;
  vector<FldArrayF*> structFc; vector<FldArrayF*> unstrFc;
  vector<E_Int> nitc; vector<E_Int> njtc; vector<E_Int> nktc;
  vector<FldArrayI*> cntc;
  vector<char*> eltTypec;
  vector<PyObject*> objstc, objutc;
  vector<E_Int> resc;
  E_Int nzonesc = 0;
  skipNoCoord = false;
  skipStructured = false;
  skipUnstructured = true;
  skipDiffVars = true;
  isOk = K_ARRAY::getFromArrays(
    arraysc, resc, structVarStringc, unstrVarStringc,
    structFc, unstrFc, nitc, njtc, nktc, cntc, eltTypec,
    objstc, objutc,
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  if (isOk == -1)
  {
    for (E_Int v1 = 0; v1 < nzones; v1++) RELEASESHAREDS(objst[v1], structF[v1]);
    PyErr_SetString(PyExc_TypeError,
                    "merge: arrays defining fields at centers are not valid.");
    return NULL;
  }

  nzonesc = structFc.size();
  E_Int both = 0;// both nodes and centers are merged -> both=1
  if (nzonesc > 0) both = 1;
  if (both == 1)
  {
    // check number of zones
    if (nzonesc != nzones)
    {
      for (E_Int v1 = 0; v1 < nzones; v1++) RELEASESHAREDS(objst[v1], structF[v1]);
      for (E_Int v1 = 0; v1 < nzonesc; v1++) RELEASESHAREDS(objstc[v1], structFc[v1]);

      PyErr_SetString(PyExc_TypeError,
                      "merge: number of arrays at nodes and centers must be equal.");
      return NULL;
    }
    // check dimensions
    for (E_Int v = 0; v < nzonesc; v++)
    {
      if (nitc[v] != K_FUNC::E_max(1,nit[v]-1) ||
          njtc[v] != K_FUNC::E_max(1,njt[v]-1) ||
          nktc[v] != K_FUNC::E_max(1,nkt[v]-1))
      {
        for (E_Int v1 = 0; v1< nzones; v1++) RELEASESHAREDS(objst[v1], structF[v1]);
        for (E_Int v1 = 0; v1< nzones; v1++) RELEASESHAREDS(objstc[v1], structFc[v1]);

        PyErr_SetString(PyExc_TypeError,
                        "merge: inconsistent dimensions between arrays at nodes and centers.");
        return NULL;
      }
    }
  }//fin des checks sur les centres

  //========================================================
  // Coordonnees deja verifiees dans getFromArrays
  // imposer ici que les champs sont tous dans le meme ordre
  //========================================================
  E_Int posx = K_ARRAY::isCoordinateXPresent(structVarString[0]); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(structVarString[0]); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(structVarString[0]); posz++;
  for (E_Int v = 0; v < nzones; v++)
  {
    if (api == -1) api = structF[v]->getApi();
    char* varString = structVarString[v];
    E_Int posxi = K_ARRAY::isCoordinateXPresent(varString); posxi++;
    E_Int posyi = K_ARRAY::isCoordinateYPresent(varString); posyi++;
    E_Int poszi = K_ARRAY::isCoordinateZPresent(varString); poszi++;
    if (posxi != posx || posyi != posy || poszi != posz)
    {
      PyErr_SetString(PyExc_TypeError,
                      "merge: coordinates are not located at same position.");
      for (E_Int v1 = 0; v1 < nzones; v1++) RELEASESHAREDS(objst[v1], structF[v1]);
      return NULL;
    }
  }
  if (api == -1) api = 1;

  // calcul de la normale aux centres sur la premiere couche en k
  E_Int nptsAll = 0;
  for (E_Int v = 0; v < nzones; v++) nptsAll += (nit[v]-1)*(njt[v]-1);
  FldArrayF* mergedField = NULL;
  FldArrayF* mergedFieldc = NULL;

  //---------
  E_Int oi = -1; E_Int oj = 2; E_Int ok = 3;
  vector<E_Int> pos1; vector<E_Int> pos2;
  vector<E_Int> posc1; vector<E_Int> posc2;

  char* varString = new char [strlen(structVarString[0])*2+4];
  K_ARRAY::getPosition(structVarString[0], structVarString[0], pos1, pos2, varString);
  delete [] varString;
  if (nzonesc > 0)
  {
    char* varStringc = new char [strlen(structVarStringc[0])*2+4];
    K_ARRAY::getPosition(structVarStringc[0], structVarStringc[0], posc1, posc2, varStringc);
    delete [] varStringc;
  }

  //E_Int rank = 0; E_Int idum = -1; E_Int ncandidates = 0;

  // -----------------
  /*0 - pour chaque bloc : construction des facettes et de la liste des voisins par facette */
  E_Int nio, njo, nko, nioc, njoc, nkoc;
  Block* blk1; Block* blk2;
  list<Block*> listOfBlks;
  list<BlockEdge*> listOfMergeableEdges;
  buildListOfBlks(nit, njt, nkt, structF, structFc, listOfBlks);
  E_Int compt = 0;
  E_Int res, res2;
  /* 1 - determination des facettes mergeables */
  nzones = listOfBlks.size();
  list<Block*>::iterator itr;
  list<Block*>::iterator itr2;
  list<BlockEdge*>::iterator itrf;
  list<BlockEdge*>::iterator itrf3;
//   list<BlockEdge*> candidates;
  E_Int nmergeMin;
  E_Int found = 0;
  E_Int nic1, njc1, nkc1, nic2, njc2, nkc2;
  for (itr = listOfBlks.begin(); itr != listOfBlks.end(); itr++)
  {
    for (itr2 = itr; itr2 != listOfBlks.end(); itr2++)
    {
      if (*itr2 != *itr)
        buildListOfMergeablesEdges(dim, alphaRef, posx, posy, posz,
                                   *itr, *itr2, listOfMergeableEdges, tol);
    }
  }

  /* 2 - determination du grade */
  for (itrf = listOfMergeableEdges.begin();
       itrf != listOfMergeableEdges.end(); itrf++)
    gradeEdge(dim, *itrf, listOfBlks, listOfMergeableEdges, tol);

  /* 3 - si pas de facettes mergeables : fin */
  E_Int gradeMin = 1000000; BlockEdge* faceMin = NULL;
  if (listOfMergeableEdges.size() == 0 ) goto end;
  restart:;
  /* 4 - recuperation de la facette de plus petit grade */
  gradeMin = 1000000; faceMin = NULL; nmergeMin = 1000000;
  for (itrf = listOfMergeableEdges.begin();
       itrf != listOfMergeableEdges.end(); itrf++)
  {
    if ( (*itrf)->_grade <= gradeMin && (*itrf)->_grade >= 0)
    {
      gradeMin = (*itrf)->_grade;
    }
  }
  // the first edge of grade min is kept
  // priority is given to a face with a low nb of merging operations already done
  if (dircons == 0)
  {
    for (itrf = listOfMergeableEdges.begin(); itrf != listOfMergeableEdges.end(); itrf++)
    {
      if (( (*itrf)->_blk1->_nmerging < nmergeMin || (*itrf)->_blk2->_nmerging < nmergeMin) )
      {
        nmergeMin = K_FUNC::E_min((*itrf)->_blk1->_nmerging,(*itrf)->_blk2->_nmerging);
        faceMin = *itrf;
      }
    }
  }
  else
  {
    for (itrf = listOfMergeableEdges.begin();
         itrf != listOfMergeableEdges.end(); itrf++)
    {
      if (K_FUNC::E_abs((*itrf)->_dirBlk1) == dircons && K_FUNC::E_abs((*itrf)->_dirBlk2) == dircons)
      {(*itrf)->_grade = -10000;  gradeMin = -10000; }
    }

    found = 0; nmergeMin = 1000000;
    for (itrf = listOfMergeableEdges.begin(); itrf != listOfMergeableEdges.end(); itrf++)
    {
      if ((*itrf)->_grade == -10000)
      {
        faceMin = *itrf; found = 1;
        break;
      }
    }
    if (found == 0)  goto end;
  }
//   for (itrf = listOfMergeableEdges.begin();
//        itrf != listOfMergeableEdges.end(); itrf++)
//   {if ( (*itrf)->_grade == gradeMin ) candidates.push_back(*itrf);}

//   ncandidates = candidates.size();
//   rank  = E_Int(K_NOISE::longRand(&idum) * ncandidates);
//   ncandidates = 0;
//   for (itrf = candidates.begin();
//        itrf != candidates.end(); itrf++)
//   {
//     if ( ncandidates == rank )
//     {faceMin = *itrf; break;}
//     ncandidates++;
//   }
//   candidates.clear();

  /* 7 - join blk1 et blk2 -> mergedBlk */
  blk1 = faceMin->_blk1; blk2 = faceMin->_blk2;
  mergedField = new FldArrayF();
  mergedFieldc = new FldArrayF();

  if (both == 0)
  {
    res = joinStructured(*(blk1->_field), blk1->_ni, blk1->_nj, blk1->_nk, posx, posy, posz,
                         *(blk2->_field), blk2->_ni, blk2->_nj, blk2->_nk, posx, posy, posz,
                         pos1, pos2, *mergedField, nio, njo, nko, tol);
    if (res != 0)
    {
      if (dircons == 1 && gradeMin == -10000) K_CONNECT::reorderStructField(nio, njo, nko, *mergedField, -1, -2, 3);
      else if (dircons == 2 && gradeMin == -10000) K_CONNECT::reorderStructField(nio, njo, nko, *mergedField, 2,-1, 3);
      else if (dircons == 3 && gradeMin == -10000) K_CONNECT::reorderStructField(nio, njo, nko, *mergedField, 3, 1, 2);

      //check du volume
      res2 = checkNegativeVolumeCells(dim, nio,njo,nko,*mergedField);
      if (res2 == 0) K_CONNECT::reorderStructField(nio, njo, nko, *mergedField, oi,oj,ok);
    }
  }
  else
  {
    nic1 = K_FUNC::E_max(1,blk1->_ni-1); nic2 = K_FUNC::E_max(1,blk2->_ni-1);
    njc1 = K_FUNC::E_max(1,blk1->_nj-1); njc2 = K_FUNC::E_max(1,blk2->_nj-1);
    nkc1 = K_FUNC::E_max(1,blk1->_nk-1); nkc2 = K_FUNC::E_max(1,blk2->_nk-1);
    res = joinBothStructured(*(blk1->_field), blk1->_ni, blk1->_nj, blk1->_nk, posx, posy, posz,
                             *(blk2->_field), blk2->_ni, blk2->_nj, blk2->_nk, posx, posy, posz,
                             *(blk1->_fieldc), nic1, njc1, nkc1,
                             *(blk2->_fieldc), nic2, njc2, nkc2,
                             pos1, pos2, posc1, posc2,
                             *mergedField, nio, njo, nko,
                             *mergedFieldc, nioc, njoc, nkoc, tol);
    if (res != 0)
    {
      if (dircons == 1 && gradeMin == -10000)
      {
        K_CONNECT::reorderStructField(nio, njo, nko, *mergedField, -1, -2, 3);
        K_CONNECT::reorderStructField(nioc, njoc, nkoc, *mergedFieldc, -1, -2, 3);
      }
      else if (dircons == 2 && gradeMin == -10000)
      {
        K_CONNECT::reorderStructField(nio, njo, nko, *mergedField, 2,-1, 3);
        K_CONNECT::reorderStructField(nioc, njoc, nkoc, *mergedFieldc, 2,-1, 3);
      }
      else if (dircons == 3 && gradeMin == -10000)
      {
        K_CONNECT::reorderStructField(nio, njo, nko, *mergedField, 3, 1, 2);
        K_CONNECT::reorderStructField(nioc, njoc, nkoc, *mergedFieldc, 3, 1, 2);
      }

      //check du volume
      res2 = checkNegativeVolumeCells(dim,nio,njo,nko,*mergedField);
      if (res2 == 0)
      {
        K_CONNECT::reorderStructField(nio, njo, nko, *mergedField, oi,oj,ok);
        K_CONNECT::reorderStructField(nioc, njoc, nkoc, *mergedFieldc, oi,oj,ok);
      }
    }
  } // fin des joins

  /* + s'assurer que la grille ne soit pas trop grande */
  if (nio*njo*nko > sizeMax || res == 0)
  {
    /* On enleve la facette choisie de la liste des facettes mergeables */
    itrf = blk1->_listOfEdges.begin();
    while (itrf != blk1->_listOfEdges.end())
    {
      itrf3 = itrf; itrf3++;
      if ((*itrf) == faceMin) blk1->_listOfEdges.erase(itrf);
      itrf = itrf3;
    }
    itrf = blk2->_listOfEdges.begin();
    while (itrf != blk2->_listOfEdges.end())
    {
      itrf3 = itrf; itrf3++;
      if ((*itrf) == faceMin) blk2->_listOfEdges.erase(itrf);
      itrf = itrf3;
    }
    itrf = listOfMergeableEdges.begin();
    while ((itrf != listOfMergeableEdges.end())&&((*itrf) != faceMin)) itrf++;
    listOfMergeableEdges.erase(itrf);
    delete faceMin;

    if (mergedField != NULL) delete mergedField;
    if (mergedFieldc != NULL) delete mergedFieldc;

    if (listOfMergeableEdges.size() > 0) goto restart;
    else goto end;
  }
  /* 5 - tag des facettes condamnees */
  deleteListOfEdgesForBlk(blk1, faceMin, listOfMergeableEdges);
  deleteListOfEdgesForBlk(blk2, faceMin, listOfMergeableEdges);

  /* 6 - enlever la facette faceMin de la liste des mergeables */
  itrf = listOfMergeableEdges.begin();
  while ((itrf != listOfMergeableEdges.end())&&((*itrf) != faceMin)) itrf++;
  listOfMergeableEdges.erase(itrf);
  delete faceMin;
  if (res != 0)
  {
    delete blk1->_field; delete blk1->_fieldc;
    blk1->_field = mergedField;
    if (both != 0) blk1->_fieldc = mergedFieldc;
    blk1->_ni = nio; blk1->_nj = njo; blk1->_nk = nko;
    blk1->_listOfEdges.clear();
    blk1->_voisins.clear();
    blk1->_nmerging++;
    /*8- delete blk2 */
    delete blk2->_field; delete blk2->_fieldc;
    deleteBlkInfo(blk2, listOfBlks); delete blk2;
  }
  /* 9 - nouvelles facettes mergeables pour blk1 : comme pt 1 */
  for (itr = listOfBlks.begin(); itr != listOfBlks.end(); itr++)
  {
    if ( *itr!= blk1)
      buildListOfMergeablesEdges(dim, alphaRef, posx, posy, posz, blk1, *itr, listOfMergeableEdges);
  }
  /* 10- grade de mergedBlk : comme pt2 */
  for (itrf = blk1->_listOfEdges.begin();
       itrf != blk1->_listOfEdges.end(); itrf++)
    gradeEdge(dim, *itrf, listOfBlks, listOfMergeableEdges, tol);

  /* 11- si liste des mergeables non vide : restart*/
  compt++;
  if (both == 0) { delete mergedFieldc; mergedFieldc = NULL; }
  if (listOfMergeableEdges.size() != 0) goto restart;
  end:;

  if (both == 0) // pas de centres
  {
    PyObject* l = PyList_New(0);
    for (itr = listOfBlks.begin(); itr != listOfBlks.end(); itr++)
    {
      nio = (*itr)->_ni; njo = (*itr)->_nj; nko = (*itr)->_nk;
      FldArrayF* field = (*itr)->_field;
      PyObject* tpl = K_ARRAY::buildArray3(*field, structVarString[0], nio, njo, nko, api);
      delete field;
      PyList_Append(l, tpl); Py_DECREF(tpl);
    }
    // nettoyages...
    if (dircons == 0)
    {
      for (itr = listOfBlks.begin(); itr != listOfBlks.end(); itr++)
      {
        for (itrf = (*itr)->_listOfEdges.begin(); itrf != (*itr)->_listOfEdges.end(); itrf++)
          delete *itrf;
        (*itr)->_listOfEdges.clear();
        (*itr)->_voisins.clear();
        delete *itr;
      }
    }
    else
    {
      for (itr = listOfBlks.begin(); itr != listOfBlks.end(); itr++)
      {
        for (itrf = (*itr)->_listOfEdges.begin(); itrf != (*itr)->_listOfEdges.end(); itrf++)
        {
          found = 0;
          for (itrf3 = listOfMergeableEdges.begin(); itrf3 != listOfMergeableEdges.end(); itrf3++)
          {
            if (*itrf == *itrf3) {found = 1; break;}
          }

          if (found == 0) delete *itrf;
        }
        (*itr)->_listOfEdges.clear();
        (*itr)->_voisins.clear();
        delete *itr;
      }
    }
    for (itrf = listOfMergeableEdges.begin(); itrf != listOfMergeableEdges.end(); itrf++)
      delete *itrf;
    listOfMergeableEdges.clear();
    listOfBlks.clear();
    return l;
  }
  else // traitement des centres et noeuds : retourne [An,Ac]
  {
    PyObject* lnodes = PyList_New(0); //noeuds
    PyObject* lcenters = PyList_New(0); //centres
    for (itr = listOfBlks.begin(); itr != listOfBlks.end(); itr++)
    {
      nio = (*itr)->_ni; njo = (*itr)->_nj; nko = (*itr)->_nk;
      FldArrayF* field = (*itr)->_field;
      PyObject* tpl = K_ARRAY::buildArray3(*field, structVarString[0], nio, njo, nko, api);
      PyList_Append(lnodes, tpl); Py_DECREF(tpl);
      delete field;
      nioc = K_FUNC::E_max(1,nio-1); njoc = K_FUNC::E_max(1,njo-1); nkoc = K_FUNC::E_max(1,nko-1);
      FldArrayF* fieldc = (*itr)->_fieldc;
      PyObject* tplc = K_ARRAY::buildArray3(*fieldc, structVarStringc[0], nioc, njoc, nkoc, api);
      PyList_Append(lcenters, tplc); Py_DECREF(tplc);
      delete fieldc;
    }
    // nettoyages...
    if (dircons == 0)
    {
      for (itr = listOfBlks.begin(); itr != listOfBlks.end(); itr++)
      {
        for (itrf = (*itr)->_listOfEdges.begin(); itrf != (*itr)->_listOfEdges.end(); itrf++)
          delete *itrf;
        (*itr)->_listOfEdges.clear();
        (*itr)->_voisins.clear();
        delete *itr;
      }
    }
    else
    {
      for (itr = listOfBlks.begin(); itr != listOfBlks.end(); itr++)
      {
        for (itrf = (*itr)->_listOfEdges.begin(); itrf != (*itr)->_listOfEdges.end(); itrf++)
        {
          found = 0;
          for (itrf3 = listOfMergeableEdges.begin(); itrf3 != listOfMergeableEdges.end(); itrf3++)
          {
            if (*itrf == *itrf3) {found = 1; break;}
          }

          if (found == 0) delete *itrf;
        }
        (*itr)->_listOfEdges.clear();
        (*itr)->_voisins.clear();
        delete *itr;
      }
    }
    for (itrf = listOfMergeableEdges.begin(); itrf != listOfMergeableEdges.end(); itrf++)
      delete *itrf;
    listOfMergeableEdges.clear();

    listOfBlks.clear();
    return Py_BuildValue("(OO)", lnodes, lcenters);
  }
}

//=============================================================================
void K_TRANSFORM::deleteBlkInfo(Block* blk, list<Block*>& listOfBlks)
{
  list<Block*>::iterator itr;
  blk->_listOfEdges.clear();
  blk->_voisins.clear();

  itr = listOfBlks.begin();
  while ((*itr) != blk) itr++;
  listOfBlks.erase(itr);
}
//=============================================================================
void K_TRANSFORM::deleteListOfEdgesForBlk(
  Block* block1, BlockEdge* faceMin,
  list<BlockEdge*>& listOfMergeableEdges)
{
  list<BlockEdge*>::iterator itrf;
  list<BlockEdge*>::iterator itrf3;
  list<BlockEdge*>::iterator itrf2;

  itrf = block1->_listOfEdges.begin();

  while (itrf != block1->_listOfEdges.end())
  {
    itrf3 = itrf; itrf3++;
    if ((*itrf) != faceMin)
    {
      itrf2 = listOfMergeableEdges.begin();
      while ((*itrf2) != (*itrf)) itrf2++;
      listOfMergeableEdges.erase(itrf2);
      if ((*itrf)->_blk1 != block1)
      {
        itrf2 = (*itrf)->_blk1->_listOfEdges.begin();
        while ((*itrf2) != (*itrf)) itrf2++;
        (*itrf)->_blk1->_listOfEdges.erase(itrf2);
      }
      else
      {
        itrf2 = (*itrf)->_blk2->_listOfEdges.begin();
        while ((*itrf2) != (*itrf)) itrf2++;
        (*itrf)->_blk2->_listOfEdges.erase(itrf2);
      }
      delete (*itrf);
      block1->_listOfEdges.erase(itrf);
    }
    itrf = itrf3;
  }
}

//=============================================================================
/* Construction des blocs ainsi que la detection des voisins.
   fieldsc sont les champs en centres
   Ici les facettes ne sont pas construites */
//=============================================================================
void K_TRANSFORM::buildListOfBlks(
  vector<E_Int>& nit, vector<E_Int>& njt, vector<E_Int>& nkt,
  vector<FldArrayF*>& coords, vector<FldArrayF*>& fieldsc,
  list<Block*>& listOfBlks)
{
  E_Int nzones = coords.size();
  FldArrayF bboxes(nzones,6);
  if (fieldsc.size() == 0)
  {
    for (E_Int v = 0; v < nzones; v++)
    {
      Block* blk = new Block;
      blk->_ni = nit[v]; blk->_nj = njt[v]; blk->_nk = nkt[v];
      blk->_field = coords[v]; blk->_nmerging = 0; blk->_fieldc = NULL;
      listOfBlks.push_back(blk);
    }
  }
  else
  {
    for (E_Int v = 0; v < nzones; v++)
    {
      Block* blk = new Block;
      blk->_ni = nit[v]; blk->_nj = njt[v]; blk->_nk = nkt[v];
      blk->_field = coords[v]; blk->_nmerging = 0; blk->_fieldc = fieldsc[v];
      listOfBlks.push_back(blk);
    }
  }
}
//=============================================================================
/* Determination des facettes mergeables: on regarde si deux aretes coincident
   entre le blk1 et blk2. Pas de raccord partiel pris en compte. */
//=============================================================================
void K_TRANSFORM::buildListOfMergeablesEdges(
  E_Int dim, E_Float alphaRef, E_Int posx, E_Int posy, E_Int posz,
  Block* blk1, Block* blk2,
  list<BlockEdge*>& listOfMergeableEdges, E_Float tol)
{
  if (dim == 2)
    return buildListOfMergeablesEdges2D(alphaRef, posx, posy, posz, blk1, blk2,
                                        listOfMergeableEdges, tol);
  else return buildListOfMergeablesEdges3D(posx, posy, posz, blk1, blk2,
                                           listOfMergeableEdges, tol);
}
//=============================================================================
void K_TRANSFORM::buildListOfMergeablesEdges3D(
  E_Int posx, E_Int posy, E_Int posz,
  Block* blk1, Block* blk2,
  list<BlockEdge*>& listOfMergeableEdges, E_Float tol)
{
  BlockEdge* newEdge = NULL;
  FldArrayF* f1 = blk1->_field;
  E_Float* xt1 = f1->begin(posx);
  E_Float* yt1 = f1->begin(posy);
  E_Float* zt1 = f1->begin(posz);
  FldArrayF* f2 = blk2->_field;
  E_Float* xt2 = f2->begin(posx);
  E_Float* yt2 = f2->begin(posy);
  E_Float* zt2 = f2->begin(posz);
  E_Int dirt1[6]; E_Int dirt2[6];
  dirt1[0] =-1; dirt1[1] = 1; dirt1[2] =-2; dirt1[3] = 2; dirt1[4] =-3; dirt1[5] = 3;
  dirt2[0] = 1; dirt2[1] =-1; dirt2[2] = 2; dirt2[3] =-2; dirt2[4] =-3; dirt2[5] = 3;
  E_Int indA1=-1, indB1=-1, indC1=-1, indD1=-1;
  E_Int indA2=-1, indB2=-1, indC2=-1, indD2=-1;
  E_Float xA1, yA1, zA1, xB1, yB1, zB1, xC1, yC1, zC1, xD1, yD1, zD1;
  E_Float xA2, yA2, zA2, xB2, yB2, zB2, xC2, yC2, zC2, xD2, yD2, zD2;
  E_Int dir1, dir2, ndir1=-1, ndir2=-1;
  E_Int ni1 = blk1->_ni; E_Int nj1 = blk1->_nj; E_Int nk1 = blk1->_nk;
  E_Int ni2 = blk2->_ni; E_Int nj2 = blk2->_nj; E_Int nk2 = blk2->_nk;
  if (blk2 == blk1) return;

  E_Int shiftk1 = (nk1-1)*ni1*nj1;
  E_Int shiftk2 = (nk2-1)*ni2*nj2;
  E_Int shiftj1 = (nj1-1)*ni1;
  E_Int shiftj2 = (nj2-1)*ni2;
  E_Int shifti1 = ni1-1;
  E_Int shifti2 = ni2-1;

  //detection si aretes coincidentes
  for (E_Int i1 = 0; i1 < 6; i1++)
  {
    dir1 = dirt1[i1];
    if (dir1 ==-1)//imin
    {
      indA1 = 0; indB1 = shiftj1; indC1 = indB1+shiftk1; indD1 = indA1+shiftk1;
    }
    else if (dir1 == 1)//imax
    {
      indA1 = shifti1; indB1 = indA1+shiftj1; indC1 = indB1+shiftk1; indD1 = indA1+shiftk1;
    }
    else if (dir1 ==-2) //jmin
    {
      indA1 = 0; indB1 = shifti1; indC1 = indB1+shiftk1; indD1 = indA1+shiftk1;
    }
    else if (dir1 == 2)
    {
      indA1 = shiftj1; indB1 = indA1+shifti1; indC1 = indB1+shiftk1; indD1 = indA1+shiftk1;
    }
    else if (dir1 ==-3)
    {
      indA1 = 0; indB1 = shifti1; indC1 = indB1 + shiftj1; indD1 = indA1 + shiftj1;
    }
    else if (dir1 == 3)
    {
      indA1 = shiftk1; indB1 = indA1+shifti1; indC1 = indB1 + shiftj1; indD1 = indA1 + shiftj1;
    }

    for (E_Int i2 = 0; i2 < 6; i2++)
    {
      dir2 = dirt2[i2];
      if (dir2 ==-1)//imin
      {
        indA2 = 0; indB2 = shiftj2; indC2 = indB2+shiftk2; indD2 = indA2+shiftk2;
      }
      else if (dir2 == 1)//imax
      {
        indA2 = shifti2; indB2 = indA2+shiftj2; indC2 = indB2+shiftk2; indD2 = indA2+shiftk2;
      }
      else if (dir2 ==-2) //jmin
      {
        indA2 = 0; indB2 = shifti2; indC2 = indB2+shiftk2; indD2 = indA2+shiftk2;
      }
      else if (dir2 == 2)
      {
        indA2 = shiftj2; indB2 = indA2+shifti2; indC2 = indB2+shiftk2; indD2 = indA2+shiftk2;
      }
      else if (dir2 ==-3)
      {
        indA2 = 0; indB2 = shifti2; indC2 = indB2 + shiftj2; indD2 = indA2 + shiftj2;
      }
      else if (dir2 == 3)
      {
        indA2 = shiftk2; indB2 = indA2+shifti2; indC2 = indB2 + shiftj2; indD2 = indA2 + shiftj2;
      }
      xA1 = xt1[indA1]; yA1 = yt1[indA1]; zA1 = zt1[indA1];
      xB1 = xt1[indB1]; yB1 = yt1[indB1]; zB1 = zt1[indB1];
      xC1 = xt1[indC1]; yC1 = yt1[indC1]; zC1 = zt1[indC1];
      xD1 = xt1[indD1]; yD1 = yt1[indD1]; zD1 = zt1[indD1];
      xA2 = xt2[indA2]; yA2 = yt2[indA2]; zA2 = zt2[indA2];
      xB2 = xt2[indB2]; yB2 = yt2[indB2]; zB2 = zt2[indB2];
      xC2 = xt2[indC2]; yC2 = yt2[indC2]; zC2 = zt2[indC2];
      xD2 = xt2[indD2]; yD2 = yt2[indD2]; zD2 = zt2[indD2];

      if (dir1 == 1 || dir1 == -1) ndir1 = nj1*nk1;
      else if (dir1 == 2 || dir1 == -2) ndir1 = ni1*nk1;
      else if (dir1 == 3 || dir1 == -3) ndir1 = ni1*nj1;

      if (dir2 == 1 || dir2 == -1) ndir2 = nj2*nk2;
      else if (dir2 == 2 || dir2 == -2) ndir2 = ni2*nk2;
      else if (dir2 == 3 || dir2 == -3) ndir2 = ni2*nj2;

      if (ndir1 == ndir2)
      {
        E_Int match = matchingEdges3D(xA1, yA1, zA1, xB1, yB1, zB1, xC1, yC1, zC1, xD1, yD1, zD1,
                                      xA2, yA2, zA2, xB2, yB2, zB2, xC2, yC2, zC2, xD2, yD2, zD2,
                                      tol);
        if (match == 4)
        {
          blk1->_voisins.push_back(blk2); blk2->_voisins.push_back(blk1);
          newEdge = new BlockEdge;
          newEdge->_grade = 0; newEdge->_externEdge = 0;
          newEdge->_blk1 = blk1; newEdge->_blk2 = blk2;
          newEdge->_dirBlk1 = dirt1[i1]; newEdge->_dirBlk2 = dirt2[i2];
          newEdge->_ptA[0] = xA1; newEdge->_ptA[1] = yA1; newEdge->_ptA[2] = zA1;
          newEdge->_ptB[0] = xB1; newEdge->_ptB[1] = yB1; newEdge->_ptB[2] = zB1;
          newEdge->_ptC[0] = xC1; newEdge->_ptC[1] = yC1; newEdge->_ptC[2] = zC1;
          newEdge->_ptD[0] = xD1; newEdge->_ptD[1] = yD1; newEdge->_ptD[2] = zD1;
          listOfMergeableEdges.push_back(newEdge);
          blk1->_listOfEdges.push_back(newEdge);
          blk2->_listOfEdges.push_back(newEdge);
          goto next;
        }
      }
    }
    next:;
  }
}
//=============================================================================
void K_TRANSFORM::buildListOfMergeablesEdges2D(
  E_Float alphaRef, E_Int posx, E_Int posy, E_Int posz,
  Block* blk1, Block* blk2,
  list<BlockEdge*>& listOfMergeableEdges, E_Float tol)
{
  BlockEdge* newEdge = NULL;
  FldArrayF* f1 = blk1->_field;
  E_Float* xt1 = f1->begin(posx);
  E_Float* yt1 = f1->begin(posy);
  E_Float* zt1 = f1->begin(posz);
  FldArrayF* f2 = blk2->_field;
  E_Float* xt2 = f2->begin(posx);
  E_Float* yt2 = f2->begin(posy);
  E_Float* zt2 = f2->begin(posz);
  E_Int dirt1[4]; E_Int dirt2[4];
  dirt1[0] =-1; dirt1[1] = 1; dirt1[2] =-2; dirt1[3] = 2;
  dirt2[0] = 1; dirt2[1] =-1; dirt2[2] = 2; dirt2[3] =-2;
  E_Int indA10, indB10, indC10, indD10, indA20, indB20, indC20, indD20;
  E_Int indA1=-1, indB1=-1, indC1=-1, indD1=-1, indA2=-1, indB2=-1, indC2=-1, indD2=-1;
  E_Float xA1, yA1, zA1, xB1, yB1, zB1, xA2, yA2, zA2, xB2, yB2, zB2;
  E_Int dir1, dir2, ndir1, ndir2;
  E_Int ni1 = blk1->_ni; E_Int nj1 = blk1->_nj;
  E_Int ni2 = blk2->_ni; E_Int nj2 = blk2->_nj;
  E_Float dx1, dy1, dz1, dx2, dy2, dz2;
  E_Int ie1=0, ie2=0, inc1, inc2, ishift1, ishift2, nmatch;
  E_Int grade = 0;
  if (blk2 == blk1) return;
  E_Int freezeSE = 1;
  if (K_FUNC::E_abs(alphaRef-180.) < tol) freezeSE = 0;
  indA10 = 0; indB10 = ni1-1; indC10 = indB10+(nj1-1)*ni1; indD10 = indA10 + (nj1-1)*ni1;
  indA20 = 0; indB20 = ni2-1; indC20 = indB20+(nj2-1)*ni2; indD20 = indA20 + (nj2-1)*ni2;
  E_Float dalpha;
  //detection si aretes coincidentes
  for (E_Int i1 = 0; i1 < 4; i1++)
  {
    dir1 = dirt1[i1];
    for (E_Int i2 = 0; i2 < 4; i2++)
    {
      dir2 = dirt2[i2];
      if (dir1 ==-1) {indA1 = indD10; indB1 = indA10; indC1 = indB10; indD1 = indC10;}//imin
      else if (dir1 == 1){indA1 = indB10; indB1 = indC10; indC1 = indD10; indD1 = indA10;}//imax
      else if (dir1 ==-2){indA1 = indA10; indB1 = indB10; indC1 = indC10; indD1 = indD10;}//jmin
      else if (dir1 == 2){indA1 = indC10; indB1 = indD10; indC1 = indA10; indD1 = indB10;}//jmax

      if (dir2 ==-1) {indA2 = indD20; indB2 = indA20; indC2 = indB20; indD2 = indC20;}//imin
      else if (dir2 == 1){indA2 = indB20; indB2 = indC20; indC2 = indD20; indD2 = indA20;}//imax
      else if (dir2 ==-2){indA2 = indA20; indB2 = indB20; indC2 = indC20; indD2 = indD20;}//jmin
      else if (dir2 == 2){indA2 = indC20; indB2 = indD20; indC2 = indA20; indD2 = indB20;}//jmax

      xA1 = xt1[indA1]; yA1 = yt1[indA1]; zA1 = zt1[indA1];
      xB1 = xt1[indB1]; yB1 = yt1[indB1]; zB1 = zt1[indB1];
      xA2 = xt2[indA2]; yA2 = yt2[indA2]; zA2 = zt2[indA2];
      xB2 = xt2[indB2]; yB2 = yt2[indB2]; zB2 = zt2[indB2];

      if (dir1 == 1 || dir1 == -1) ndir1 = nj1;
      else  ndir1 = ni1;
      if (dir2 == 1 || dir2 == -1) ndir2 = nj2;
      else ndir2 = ni2;
      grade = 0;

      if (ndir1 == ndir2)
      {
        E_Int match = matchingEdges2D(xA1, yA1, zA1, xB1, yB1, zB1, xA2, yA2, zA2, xB2, yB2, zB2, tol);
        if (match == 0)
        {
          // detection maillage en C
          dx1 = xB1-xA1; dx2 = xB2-xA2;
          dy1 = yB1-yA1; dy2 = yB2-yA2;
          dz1 = zB1-zA1; dz2 = zB2-zA2;
          if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol &&
              K_FUNC::E_abs(dx2)<tol && K_FUNC::E_abs(dy2)<tol && K_FUNC::E_abs(dz2)<tol )
          {
            nmatch = 0; inc1 = 0; inc2 = 0; ishift1 = 1; ishift2 = 1;

            switch (dir1)
            {
              case -1:
                ie1 = nj1; ishift1 = ni1-1;
                break;
              case  1:
                inc1 = (ni1-1); ie1 = nj1; ishift1 = ni1-1;
              break;
              case -2:
                ie1 = ni1;
                break;
              case 2:
                inc1 = (nj1-1)*ni1; ie1 = ni1;
                break;
            }
            switch (dir2)
            {
              case -1:
                ie2 = nj2; ishift2 = ni2-1;
                break;
              case 1:
                inc2 = (ni1-1); ie2 = nj2; ishift2 = ni2-1;
              break;
              case -2:
                ie2 = ni2;
                break;
              case 2:
                inc2 = (nj2-1)*ni2; ie2 = ni2;
                break;
            }
            for (E_Int ip1 = 0; ip1 < ie1; ip1++)
            {
              E_Int ind1 = ip1*ishift1 + inc1;
              for (E_Int ip2 = 0; ip2 < ie2; ip2++)
              {
                E_Int ind2 = ip2*ishift2 + inc2;
                dx1 = xt1[ind1]-xt2[ind2]; dy1 = yt1[ind1]-yt2[ind2]; dz1 = zt1[ind1]-zt2[ind2];
                if ( K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol )
                {nmatch++; goto suivant;}
              }
              suivant:;
            }

            //verifie que les pts sont coincidents 2 a 2
            if (nmatch == ndir1) match = 2;
          }
        }

        if (match == 2)
        {
          dalpha = 0.;
          if (freezeSE == 1) dalpha = compAngleBetweenZones(dir1, ni1, nj1, xt1, yt1, zt1, dir2, ni2, nj2, xt2, yt2, zt2);

          if (K_FUNC::E_abs(dalpha) < alphaRef)
          {
            increaseGradeForRotation(indA1, indB1, indC1, indD1, xt1, yt1, zt1,
                                     indA2, indB2, indC2, indD2, xt2, yt2, zt2, grade);
            blk1->_voisins.push_back(blk2); blk2->_voisins.push_back(blk1);
            newEdge = new BlockEdge;
            newEdge->_grade = grade; newEdge->_externEdge = 0;
            newEdge->_blk1 = blk1; newEdge->_blk2 = blk2;
            newEdge->_dirBlk1 = dirt1[i1]; newEdge->_dirBlk2 = dirt2[i2];
            newEdge->_ptA[0] = xA1; newEdge->_ptA[1] = yA1; newEdge->_ptA[2] = zA1;
            newEdge->_ptB[0] = xB1; newEdge->_ptB[1] = yB1; newEdge->_ptB[2] = zB1;
            listOfMergeableEdges.push_back(newEdge);
            blk1->_listOfEdges.push_back(newEdge);
            blk2->_listOfEdges.push_back(newEdge);
          }
        }
      }
    }
  }
}
//==============================================================================
/**/
//==============================================================================
E_Float K_TRANSFORM::compAngleBetweenZones(
  E_Int dir1, E_Int ni1, E_Int nj1,
  E_Float* xt1, E_Float* yt1, E_Float* zt1,
  E_Int dir2, E_Int ni2, E_Int nj2,
  E_Float* xt2, E_Float* yt2, E_Float* zt2)
{
  E_Float tol = 1.e-10;
  E_Float dalpha = 0., alpha0;
  E_Float ptA1[3], ptA2[3], ptB1[3], ptB2[3], ptC1[3], ptC2[3], ptD1[3], ptD2[3];
  E_Int shift1=-1, shift2=-1, shiftmax1=-1, inci1=-1, incj1=-1, inci2=-1, incj2=-1;
  E_Int indA10=-1, indB10=-1, indC10=-1, indD10=-1, indA20=-1, indB20=-1, indC20=-1, indD20=-1;
  E_Int ind1=-1, ind2=-1;
  E_Int reverse = 1; // reverse = 1 : no reverse, -1 : reverse

  switch (dir1)
  {
    case -1:
      inci1 = 1; incj1 = ni1; indA10 = 0; shift1 = ni1; shiftmax1 = nj1-1; ind1 = 0; ind2 = shift1;
      break;
    case 1:
      inci1 = 1; incj1 = ni1; indA10 = ni1-2; shift1 = ni1; shiftmax1 = nj1-1; ind1 = ni1-1; ind2 = ind1+shift1;
      break;
    case -2:
      inci1 = 1; incj1 = ni1; indA10 = 0; shift1 = 1; shiftmax1 = ni1-1; ind1 = 0; ind2 = ind1+ shift1;
      break;
    case 2:
      inci1 = 1; incj1 = ni1; indA10 = (nj1-2)*ni1; shift1 = 1; shiftmax1 = ni1-1; ind1 = (nj1-1)*ni1; ind2 = ind1 + shift1;
      break;
  }
  switch ( dir2 )
  {
    case -1:
      inci2 = 1; incj2 = ni2; indA20 = 0; shift2 = ni2; ind2 = 0;
      if (K_FUNC::E_abs(xt1[ind1]-xt2[ind2]) > tol ||
          K_FUNC::E_abs(yt1[ind1]-yt2[ind2]) > tol ||
          K_FUNC::E_abs(zt1[ind1]-zt2[ind2]) > tol)//reverse
      {
        indA20 = indA20+(nj2-2)*ni2; shift2 = -ni2;
      }
      break;
    case 1:
      inci2 = 1; incj2 = ni2; indA20 = ni2-2; shift2 = ni2; ind2 = ni2-1;
       if (K_FUNC::E_abs(xt1[ind1]-xt2[ind2]) > tol ||
           K_FUNC::E_abs(yt1[ind1]-yt2[ind2]) > tol ||
           K_FUNC::E_abs(zt1[ind1]-zt2[ind2]) > tol)//reverse
       {

         indA20 = indA20+(nj2-2)*ni2; shift2 = -ni2;
       }
      break;
    case -2:
      inci2 = 1; incj2 = ni2; indA20 = 0; shift2 = 1; ind2 = 0;
      if (K_FUNC::E_abs(xt1[ind1]-xt2[ind2]) > tol ||
          K_FUNC::E_abs(yt1[ind1]-yt2[ind2]) > tol ||
          K_FUNC::E_abs(zt1[ind1]-zt2[ind2]) > tol)//reverse
      {

        indA20 = indA20+ni2-2; shift2 = -1;
      }
      break;
    case 2:
      inci2 = 1; incj2 = ni2; indA20 = (nj2-2)*ni2; shift2 = 1; ind2 = (nj2-1)*ni2;
      if (K_FUNC::E_abs(xt1[ind1]-xt2[ind2]) > tol ||
          K_FUNC::E_abs(yt1[ind1]-yt2[ind2]) > tol ||
          K_FUNC::E_abs(zt1[ind1]-zt2[ind2]) > tol)//reverse
      {
        indA20 = indA20 +(ni2-2); shift2 = -1;
      }
      break;
  }

  for (E_Int i1 = 0; i1 < shiftmax1; i1++)
  {
    indB10 = indA10 + inci1; indC10 = indA10 + inci1 + incj1; indD10 = indA10 + incj1;
    indB20 = indA20 + inci2; indC20 = indA20 + inci2 + incj2; indD20 = indA20 + incj2;
    ptA1[0] = xt1[indA10]; ptB1[0] = xt1[indB10]; ptC1[0] = xt1[indC10]; ptD1[0] = xt1[indD10];
    ptA1[1] = yt1[indA10]; ptB1[1] = yt1[indB10]; ptC1[1] = yt1[indC10]; ptD1[1] = yt1[indD10];
    ptA1[2] = zt1[indA10]; ptB1[2] = zt1[indB10]; ptC1[2] = zt1[indC10]; ptD1[2] = zt1[indD10];
    ptA2[0] = xt2[indA20]; ptB2[0] = xt2[indB20]; ptC2[0] = xt2[indC20]; ptD2[0] = xt2[indD20];
    ptA2[1] = yt2[indA20]; ptB2[1] = yt2[indB20]; ptC2[1] = yt2[indC20]; ptD2[1] = yt2[indD20];
    ptA2[2] = zt2[indA20]; ptB2[2] = zt2[indB20]; ptC2[2] = zt2[indC20]; ptD2[2] = zt2[indD20];

    ptA1[0] = xt1[indA10]; ptB1[0] = xt1[indB10]; ptC1[0] = xt1[indC10]; ptD1[0] = xt1[indD10];
    ptA1[1] = yt1[indA10]; ptB1[1] = yt1[indB10]; ptC1[1] = yt1[indC10]; ptD1[1] = yt1[indD10];
    ptA1[2] = zt1[indA10]; ptB1[2] = zt1[indB10]; ptC1[2] = zt1[indC10]; ptD1[2] = zt1[indD10];

    if ( reverse == 1) alpha0 = K_COMPGEOM::getAlphaAngleBetweenQuads(ptA1, ptB1, ptC1, ptD1,
                                                                      ptA2, ptB2, ptC2, ptD2);
    else alpha0 = K_COMPGEOM::getAlphaAngleBetweenQuads(ptA1, ptB1, ptC1, ptD1,
                                                        ptA2, ptD2, ptC2, ptB2);
    if ( alpha0 != -1000.)
    {
      dalpha = K_FUNC::E_max(dalpha,K_FUNC::E_abs(180.-alpha0));
    }
    indA10 = indA10 + shift1; indA20 = indA20 + shift2;
  }
  return dalpha;
}
//=================================================================================
/* Increase grade for edge A1B1=A2B2 if their merging will lead to rotating cells*/
//=================================================================================
void K_TRANSFORM::increaseGradeForRotation(E_Int indA1, E_Int indB1, E_Int indC1, E_Int indD1,
                                           E_Float* xt1, E_Float* yt1, E_Float* zt1,
                                           E_Int indA2, E_Int indB2, E_Int indC2, E_Int indD2,
                                           E_Float* xt2, E_Float* yt2, E_Float* zt2,
                                           E_Int& grade)
{
  // test triangle A1B1C1 puis A1B1D1
  E_Float xAB1 = xt1[indB1]-xt1[indA1];
  E_Float yAB1 = yt1[indB1]-yt1[indA1];
  E_Float zAB1 = zt1[indB1]-zt1[indA1];

  E_Float xAC1 = xt1[indC1]-xt1[indA1];
  E_Float yAC1 = yt1[indC1]-yt1[indA1];
  E_Float zAC1 = zt1[indC1]-zt1[indA1];

  E_Float xBC1 = xt1[indC1]-xt1[indB1];
  E_Float yBC1 = yt1[indC1]-yt1[indB1];
  E_Float zBC1 = zt1[indC1]-zt1[indB1];

  E_Float xAD1 = xt1[indD1]-xt1[indA1];
  E_Float yAD1 = yt1[indD1]-yt1[indA1];
  E_Float zAD1 = zt1[indD1]-zt1[indA1];

  E_Float xBD1 = xt1[indD1]-xt1[indB1];
  E_Float yBD1 = yt1[indD1]-yt1[indB1];
  E_Float zBD1 = zt1[indD1]-zt1[indB1];

  E_Float d2AB1 = xAB1*xAB1+yAB1*yAB1+zAB1*zAB1; E_Float dAB1 = sqrt(d2AB1);
  E_Float d2AC1 = xAC1*xAC1+yAC1*yAC1+zAC1*zAC1;
  E_Float d2BC1 = xBC1*xBC1+yBC1*yBC1+zBC1*zBC1; E_Float dBC1 = sqrt(d2BC1);
  E_Float d2AD1 = xAD1*xAD1+yAD1*yAD1+zAD1*zAD1; E_Float dAD1 = sqrt(d2AD1);
  E_Float d2BD1 = xBD1*xBD1+yBD1*yBD1+zBD1*zBD1;

  E_Float cosDAB1 = 0.5*(-d2BD1+d2AD1+d2AB1)/(dAB1*dAD1);
  E_Float cosABC1 = 0.5*(-d2AC1+d2AB1+d2BC1)/(dAB1*dBC1);

  // test triangle A2B2C2 puis A2B2D2
  E_Float xAB2 = xt2[indB2]-xt2[indA2];
  E_Float yAB2 = yt2[indB2]-yt2[indA2];
  E_Float zAB2 = zt2[indB2]-zt2[indA2];

  E_Float xAC2 = xt2[indC2]-xt2[indA2];
  E_Float yAC2 = yt2[indC2]-yt2[indA2];
  E_Float zAC2 = zt2[indC2]-zt2[indA2];

  E_Float xBC2 = xt2[indC2]-xt2[indB2];
  E_Float yBC2 = yt2[indC2]-yt2[indB2];
  E_Float zBC2 = zt2[indC2]-zt2[indB2];

  E_Float xAD2 = xt2[indD2]-xt2[indA2];
  E_Float yAD2 = yt2[indD2]-yt2[indA2];
  E_Float zAD2 = zt2[indD2]-zt2[indA2];

  E_Float xBD2 = xt2[indD2]-xt2[indB2];
  E_Float yBD2 = yt2[indD2]-yt2[indB2];
  E_Float zBD2 = zt2[indD2]-zt2[indB2];

  E_Float d2AB2 = xAB2*xAB2+yAB2*yAB2+zAB2*zAB2; E_Float dAB2 = sqrt(d2AB2);
  E_Float d2AC2 = xAC2*xAC2+yAC2*yAC2+zAC2*zAC2;
  E_Float d2BC2 = xBC2*xBC2+yBC2*yBC2+zBC2*zBC2; E_Float dBC2 = sqrt(d2BC2);
  E_Float d2AD2 = xAD2*xAD2+yAD2*yAD2+zAD2*zAD2; E_Float dAD2 = sqrt(d2AD2);
  E_Float d2BD2 = xBD2*xBD2+yBD2*yBD2+zBD2*zBD2;

  E_Float cosDAB2 = 0.5*(-d2BD2+d2AD2+d2AB2)/(dAB2*dAD2);
  E_Float cosABC2 = 0.5*(-d2AC2+d2AB2+d2BC2)/(dAB2*dBC2);
  E_Float cosmax1 = K_FUNC::E_max(K_FUNC::E_abs(cosABC1),K_FUNC::E_abs(cosDAB1));
  E_Float cosmax2 = K_FUNC::E_max(K_FUNC::E_abs(cosABC2),K_FUNC::E_abs(cosDAB2));
  E_Float cos60 = cos(K_CONST::E_PI/3.); E_Float cos30 = cos(K_CONST::E_PI/6.); E_Float cos45 = cos(K_CONST::E_PI/4.);

  if (cosmax1 > cos60 || cosmax2 > cos60 ) grade = grade+2;
  if (cosmax1 > cos45 || cosmax2 > cos45 ) grade = grade+2;
  if (cosmax1 > cos30 || cosmax2 > cos30 ) grade = grade+10005;
  return;
}

//=============================================================================
/* Retourne le nb de pts coincident entre 2 faces */
//=============================================================================
E_Int K_TRANSFORM::matchingEdges3D(E_Float xA1, E_Float yA1, E_Float zA1,
                                   E_Float xB1, E_Float yB1, E_Float zB1,
                                   E_Float xC1, E_Float yC1, E_Float zC1,
                                   E_Float xD1, E_Float yD1, E_Float zD1,
                                   E_Float xA2, E_Float yA2, E_Float zA2,
                                   E_Float xB2, E_Float yB2, E_Float zB2,
                                   E_Float xC2, E_Float yC2, E_Float zC2,
                                   E_Float xD2, E_Float yD2, E_Float zD2,
                                   E_Float tol)
{
  E_Float dx1, dy1, dz1;
  E_Int nmatch = 0;
  dx1 = xA2-xA1; dy1 = yA2-yA1; dz1 = zA2-zA1;
  if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol) {nmatch++; goto pt2;}
  dx1 = xB2-xA1; dy1 = yB2-yA1; dz1 = zB2-zA1;
  if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol) {nmatch++; goto pt2;}
  dx1 = xC2-xA1; dy1 = yC2-yA1; dz1 = zC2-zA1;
  if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol) {nmatch++; goto pt2;}
  dx1 = xD2-xA1; dy1 = yD2-yA1; dz1 = zD2-zA1;
  if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol) {nmatch++; goto pt2;}
  return 0;
  pt2:;
  dx1 = xA2-xB1; dy1 = yA2-yB1; dz1 = zA2-zB1;
  if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol) {nmatch++; goto pt3;}
  dx1 = xB2-xB1; dy1 = yB2-yB1; dz1 = zB2-zB1;
  if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol) {nmatch++; goto pt3;}
  dx1 = xC2-xB1; dy1 = yC2-yB1; dz1 = zC2-zB1;
  if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol) {nmatch++; goto pt3;}
  dx1 = xD2-xB1; dy1 = yD2-yB1; dz1 = zD2-zB1;
  if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol) {nmatch++; goto pt3;}
  return 0;
  pt3:;
  dx1 = xA2-xC1; dy1 = yA2-yC1; dz1 = zA2-zC1;
  if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol) {nmatch++; goto pt4;}
  dx1 = xB2-xC1; dy1 = yB2-yC1; dz1 = zB2-zC1;
  if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol) {nmatch++; goto pt4;}
  dx1 = xC2-xC1; dy1 = yC2-yC1; dz1 = zC2-zC1;
  if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol) {nmatch++; goto pt4;}
  dx1 = xD2-xC1; dy1 = yD2-yC1; dz1 = zD2-zC1;
  if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol) {nmatch++; goto pt4;}
  return 0;
  pt4:;
  dx1 = xA2-xD1; dy1 = yA2-yD1; dz1 = zA2-zD1;
  if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol) {nmatch++; goto end;}
  dx1 = xB2-xD1; dy1 = yB2-yD1; dz1 = zB2-zD1;
  if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol) {nmatch++; goto end;}
  dx1 = xC2-xD1; dy1 = yC2-yD1; dz1 = zC2-zD1;
  if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol) {nmatch++; goto end;}
  dx1 = xD2-xD1; dy1 = yD2-yD1; dz1 = zD2-zD1;
  if (K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol) {nmatch++; goto end;}
  end:;
  return nmatch;
}
//=============================================================================
/* Retourne 1 si les 2 edges coincident pour leur deux sommets
   0 sinon */
//=============================================================================
E_Int K_TRANSFORM::matchingEdges2D(E_Float xA1, E_Float yA1, E_Float zA1,
                                   E_Float xB1, E_Float yB1, E_Float zB1,
                                   E_Float xA2, E_Float yA2, E_Float zA2,
                                   E_Float xB2, E_Float yB2, E_Float zB2,
                                   E_Float tol)
{
  E_Float dx1, dy1, dz1, dx2, dy2, dz2;

  //test arete A1B1 avec A2B2
  dx1 = xA2-xA1; dx2 = xB2-xB1;
  dy1 = yA2-yA1; dy2 = yB2-yB1;
  dz1 = zA2-zA1; dz2 = zB2-zB1;
  if( K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol)
  {
    // verif B1 et B2 coincidents
    if( K_FUNC::E_abs(dx2)<tol && K_FUNC::E_abs(dy2)<tol && K_FUNC::E_abs(dz2)<tol)
      return 2;
    else return 1;
  }

  //test arete A1B1 avec B2A2
  dx1 = xB2-xA1; dx2 = xA2-xB1;
  dy1 = yB2-yA1; dy2 = yA2-yB1;
  dz1 = zB2-zA1; dz2 = zA2-zB1;
  if( K_FUNC::E_abs(dx1)<tol && K_FUNC::E_abs(dy1)<tol && K_FUNC::E_abs(dz1)<tol)
  {
    // verif B1 et B2 coincidents
    if( K_FUNC::E_abs(dx2)<tol && K_FUNC::E_abs(dy2)<tol && K_FUNC::E_abs(dz2)<tol)
      return 2;
    else return 1;
  }
  return 0;
}
//=============================================================================
void K_TRANSFORM::gradeEdge(E_Int dim,
                            BlockEdge* edge, list<Block*>& listOfBlks,
                            list<BlockEdge*>& listOfMergeableEdges,
                            E_Float tol)
{
  if (dim == 2) return gradeEdge2D(edge, listOfBlks, listOfMergeableEdges, tol);
  else return gradeEdge3D(edge, listOfBlks, listOfMergeableEdges, tol);
}
//=============================================================================
void K_TRANSFORM::gradeEdge3D(BlockEdge* edge, list<Block*>& listOfBlks,
                              list<BlockEdge*>& listOfMergeableEdges,
                              E_Float tol)
{
  Block* blk1 = edge->_blk1;
  Block* blk2 = edge->_blk2;
  edge->_grade = K_FUNC::E_max(blk1->_nmerging, blk2->_nmerging);
  // test si facette contient un pt exterieur
  if (edge->_externEdge == 1) edge->_grade +=2;
  E_Float pt11[3]; E_Float pt12[3]; E_Float pt21[3]; E_Float pt22[3];
  E_Int diff = 1;
  list<Block*>::iterator itr2;
  list<BlockEdge*>::iterator itrf1;
  list<BlockEdge*>::iterator itrf2;
  list<BlockEdge*>::iterator itrf3;
  E_Int dirt2[6];
  E_Int match0, match;
  E_Int indA2 =-1, indB2 =-1, indC2 =-1, indD2 = -1;
  E_Float ptA2[3]; E_Float ptB2[3]; E_Float ptC2[3]; E_Float ptD2[3];
  E_Float ptA3[3]; E_Float ptB3[3]; E_Float ptC3[3]; E_Float ptD3[3];

  dirt2[0] = 1; dirt2[1] =-1; dirt2[2] = 2; dirt2[3] =-2; dirt2[4] = 3; dirt2[5] =-3;
  //blk1
  for (itrf2 = blk1->_listOfEdges.begin(); itrf2 != blk1->_listOfEdges.end(); itrf2++)
  {
    diff = 0;
    if ((*itrf2)->_blk1 == blk1 ) diff = (*itrf2)->_dirBlk1-edge->_dirBlk1;
    else if  ((*itrf2)->_blk2 == blk1 ) diff = (*itrf2)->_dirBlk2-edge->_dirBlk1;
    if (diff != 0 && diff != 2 && diff != 4 && diff != 6)
    {
      //edge->_grade++;
      itrf1 = listOfMergeableEdges.begin();
      while ((itrf1 != listOfMergeableEdges.end())&&((*itrf1) != (*itrf2))) itrf1++;
      if ((*itrf1) == (*itrf2)) edge->_grade++;
      /* Recherche si creation d'une nouvelle facette mergeable */
      for (itrf3 = blk2->_listOfEdges.begin();
           itrf3 != blk2->_listOfEdges.end();
           itrf3++)
      {
        // recuperation des extremites de itrf2 et itrf3
        for ( E_Int i = 0; i < 3; i++)
        {
          ptA2[i] = (*itrf2)->_ptA[i]; ptB2[i] = (*itrf2)->_ptB[i];
          ptC2[i] = (*itrf2)->_ptC[i]; ptD2[i] = (*itrf2)->_ptD[i];
          ptA3[i] = (*itrf3)->_ptA[i]; ptB3[i] = (*itrf3)->_ptB[i];
          ptC3[i] = (*itrf3)->_ptC[i]; ptD3[i] = (*itrf3)->_ptD[i];
        }
        // Trouver l arete commune (found = 1)
        E_Int found = 0;
        //--A2B2 et A3B3
        match0 = matchingEdges2D(ptA2[0], ptA2[1], ptA2[2], ptB2[0], ptB2[1], ptB2[2],
                                 ptA3[0], ptA3[1], ptA3[2], ptB3[0], ptB3[1], ptB3[2], tol);
        if (match0 != 0 )
        {
          for ( E_Int i = 0; i < 3; i++)
          {
            pt11[i] = ptC2[i]; pt12[i] = ptD2[i]; pt21[i] = ptC3[i]; pt22[i] = ptD3[i];
          }
          found = 1; goto next;
        }
        //A2B2 et B3C3
        match0 = matchingEdges2D(ptA2[0], ptA2[1], ptA2[2], ptB2[0], ptB2[1], ptB2[2],
                                 ptB3[0], ptB3[1], ptB3[2], ptC3[0], ptC3[1], ptC3[2], tol);
        if (match0 != 0 )
        {
          for ( E_Int i = 0; i < 3; i++)
          {
            pt11[i] = ptC2[i]; pt12[i] = ptD2[i]; pt21[i] = ptD3[i]; pt22[i] = ptA3[i];
          }
          found = 1; goto next;
        }
        //A2B2 et C3D3
        match0 = matchingEdges2D(ptA2[0], ptA2[1], ptA2[2], ptB2[0], ptB2[1], ptB2[2],
                                 ptC3[0], ptC3[1], ptC3[2], ptD3[0], ptD3[1], ptD3[2], tol);
        if (match0 != 0 )
        {
          for ( E_Int i = 0; i < 3; i++)
          {
            pt11[i] = ptC2[i]; pt12[i] = ptD2[i]; pt21[i] = ptA3[i]; pt22[i] = ptB3[i];
          }
          found = 1; goto next;
        }
        //A2B2 et D3A3
        match0 = matchingEdges2D(ptA2[0], ptA2[1], ptA2[2], ptB2[0], ptB2[1], ptB2[2],
                                 ptD3[0], ptD3[1], ptD3[2], ptA3[0], ptA3[1], ptA3[2], tol);
        if (match0 != 0 )
        {
          for ( E_Int i = 0; i < 3; i++)
          {
            pt11[i] = ptC2[i]; pt12[i] = ptD2[i]; pt21[i] = ptB3[i]; pt22[i] = ptC3[i];
          }
          found = 1; goto next;
        }

        //--B2C2 et A3B3
        match0 = matchingEdges2D(ptB2[0], ptB2[1], ptB2[2], ptC2[0], ptC2[1], ptC2[2],
                                 ptA3[0], ptA3[1], ptA3[2], ptB3[0], ptB3[1], ptB3[2], tol);
        if (match0 != 0 )
        {
          for (E_Int i = 0; i < 3; i++)
          {
            pt11[i] = ptA2[i]; pt12[i] = ptD2[i]; pt21[i] = ptC3[i]; pt22[i] = ptD3[i];
          }
          found = 1; goto next;
        }
        //B2C2 et B3C3
        match0 = matchingEdges2D(ptB2[0], ptB2[1], ptB2[2], ptC2[0], ptC2[1], ptC2[2],
                                 ptB3[0], ptB3[1], ptB3[2], ptC3[0], ptC3[1], ptC3[2], tol);
        if (match0 != 0 )
        {
          for (E_Int i = 0; i < 3; i++)
          {
            pt11[i] = ptA2[i]; pt12[i] = ptD2[i]; pt21[i] = ptD3[i]; pt22[i] = ptA3[i];
          }
          found = 1; goto next;
        }
        //B2C2 et C3D3
        match0 = matchingEdges2D(ptB2[0], ptB2[1], ptB2[2], ptC2[0], ptC2[1], ptC2[2],
                                 ptC3[0], ptC3[1], ptC3[2], ptD3[0], ptD3[1], ptD3[2], tol);
        if (match0 != 0 )
        {
          for (E_Int i = 0; i < 3; i++)
          {
            pt11[i] = ptA2[i]; pt12[i] = ptD2[i]; pt21[i] = ptA3[i]; pt22[i] = ptB3[i];
          }
          found = 1; goto next;
        }
        //B2C2 et D3A3
        match0 = matchingEdges2D(ptB2[0], ptB2[1], ptB2[2], ptC2[0], ptC2[1], ptC2[2],
                                 ptD3[0], ptD3[1], ptD3[2], ptA3[0], ptA3[1], ptA3[2], tol);
        if (match0 != 0 )
        {
          for (E_Int i = 0; i < 3; i++)
          {
            pt11[i] = ptA2[i]; pt12[i] = ptD2[i]; pt21[i] = ptC3[i]; pt22[i] = ptB3[i];
          }
          found = 1; goto next;
        }

        //--C2D2 et A3B3
        match0 = matchingEdges2D(ptC2[0], ptC2[1], ptC2[2], ptD2[0], ptD2[1], ptD2[2],
                                 ptA3[0], ptA3[1], ptA3[2], ptB3[0], ptB3[1], ptB3[2], tol);
        if (match0 != 0 )
        {
          for (E_Int i = 0; i < 3; i++)
          {
            pt11[i] = ptA2[i]; pt12[i] = ptB2[i]; pt21[i] = ptC3[i]; pt22[i] = ptD3[i];
          }
          found = 1; goto next;
        }
        //C2D2 et B3C3
        match0 = matchingEdges2D(ptC2[0], ptC2[1], ptC2[2], ptD2[0], ptD2[1], ptD2[2],
                                 ptB3[0], ptB3[1], ptB3[2], ptC3[0], ptC3[1], ptC3[2], tol);
        if (match0 != 0 )
        {
          for (E_Int i = 0; i < 3; i++)
          {
            pt11[i] = ptA2[i]; pt12[i] = ptB2[i]; pt21[i] = ptA3[i]; pt22[i] = ptD3[i];}
          found = 1; goto next;
        }
        //C2D2 et C3D3
        match0 = matchingEdges2D(ptC2[0], ptC2[1], ptC2[2], ptD2[0], ptD2[1], ptD2[2],
                                 ptC3[0], ptC3[1], ptC3[2], ptD3[0], ptD3[1], ptD3[2], tol);
        if (match0 != 0 )
        {
          for (E_Int i = 0; i < 3; i++)
          {
            pt11[i] = ptA2[i]; pt12[i] = ptB2[i]; pt21[i] = ptA3[i]; pt22[i] = ptB3[i];
          }
          found = 1; goto next;
        }
        //C2D2 et D3A3
        match0 = matchingEdges2D(ptC2[0], ptC2[1], ptC2[2], ptD2[0], ptD2[1], ptD2[2],
                                 ptD3[0], ptD3[1], ptD3[2], ptA3[0], ptA3[1], ptA3[2], tol);
        if (match0 != 0 )
        {
          for (E_Int i = 0; i < 3; i++)
          {
            pt11[i] = ptA2[i]; pt12[i] = ptB2[i]; pt21[i] = ptC3[i]; pt22[i] = ptB3[i];
          }
          found = 1; goto next;
        }

        //--D2A2 et A3B3
        match0 = matchingEdges2D(ptD2[0], ptD2[1], ptD2[2], ptA2[0], ptA2[1], ptA2[2],
                                 ptA3[0], ptA3[1], ptA3[2], ptB3[0], ptB3[1], ptB3[2], tol);
        if (match0 != 0 )
        {
          for (E_Int i = 0; i < 3; i++)
          {
            pt11[i] = ptC2[i]; pt12[i] = ptB2[i]; pt21[i] = ptC3[i]; pt22[i] = ptD3[i];
          }
          found = 1; goto next;
        }
        //D2A2 et B3C3
        match0 = matchingEdges2D(ptD2[0], ptD2[1], ptD2[2], ptA2[0], ptA2[1], ptA2[2],
                                 ptB3[0], ptB3[1], ptB3[2], ptC3[0], ptC3[1], ptC3[2], tol);
        if (match0 != 0 )
        {
          for (E_Int i = 0; i < 3; i++)
          {
            pt11[i] = ptC2[i]; pt12[i] = ptB2[i]; pt21[i] = ptA3[i]; pt22[i] = ptD3[i];
          }
          found = 1; goto next;
        }
        //D2A2 et C3D3
        match0 = matchingEdges2D(ptD2[0], ptD2[1], ptD2[2], ptA2[0], ptA2[1], ptA2[2],
                                 ptC3[0], ptC3[1], ptC3[2], ptD3[0], ptD3[1], ptD3[2], tol);
        if (match0 != 0 )
        {
          for (E_Int i = 0; i < 3; i++)
          {
            pt11[i] = ptC2[i]; pt12[i] = ptB2[i]; pt21[i] = ptA3[i]; pt22[i] = ptB3[i];
          }
          found = 1; goto next;
        }
        //D2A2 et D3A3
        match0 = matchingEdges2D(ptD2[0], ptD2[1], ptD2[2], ptA2[0], ptA2[1], ptA2[2],
                                 ptD3[0], ptD3[1], ptD3[2], ptA3[0], ptA3[1], ptA3[2], tol);
        if (match0 != 0 )
        {
          for (E_Int i = 0; i < 3; i++)
          {
            pt11[i] = ptC2[i]; pt12[i] = ptB2[i]; pt21[i] = ptB3[i]; pt22[i] = ptC3[i];
          }
          found = 1; goto next;
        }

        next:;
        // si oui : recherche d un bloc dont une face est coincidente avec la face constituee de pt11,pt12,pt21,pt22
        if ( found != 0 )
        {
          for (itr2 = listOfBlks.begin(); itr2 != listOfBlks.end(); itr2++)
          {
            E_Int ni2 = (*itr2)->_ni; E_Int nj2 = (*itr2)->_nj;
            FldArrayF* f2 = (*itr2)->_field;
            E_Float* xt2 = f2->begin(1);
            E_Float* yt2 = f2->begin(2);
            E_Float* zt2 = f2->begin(3);
            E_Int shiftk2 = ni2*nj2;
            E_Int shiftj2 = (nj2-1)*ni2;
            E_Int shifti2 = ni2-1;
            for (E_Int i2 = 0; i2 < 6; i2++)
            {
              E_Int dir2 = dirt2[i2];

              switch (dir2)
              {
                case -1:
                  indA2 = 0; indB2 = shiftj2; indC2 = indB2+shiftk2; indD2 = indA2+shiftk2;
                  break;
                case  1:
                  indA2 = shifti2; indB2 = indA2+shiftj2; indC2 = indB2+shiftk2; indD2 = indA2+shiftk2;
                  break;
                case -2:
                  indA2 = 0; indB2 = shifti2; indC2 = indB2+shiftk2; indD2 = indA2+shiftk2;
                  break;
                case  2:
                  indA2 = shiftj2; indB2 = indA2+shifti2; indC2 = indB2+shiftk2; indD2 = indA2+shiftk2;
                  break;
                case -3:
                  indA2 = 0; indB2 = shifti2; indC2 = indB2 + shiftj2; indD2 = indA2 + shiftj2;
                  break;
                case  3:
                  indA2 = shiftk2; indB2 = indA2+shifti2; indC2 = indB2 + shiftj2; indD2 = indA2 + shiftj2;
                  break;
              }
              match = matchingEdges3D(pt11[0], pt11[1], pt11[2], pt12[0], pt12[1], pt12[2],
                                      pt21[0], pt21[1], pt21[2], pt22[0], pt22[1], pt22[2],
                                      xt2[indA2], yt2[indA2], zt2[indA2], xt2[indB2], yt2[indB2], zt2[indB2],
                                      xt2[indC2], yt2[indC2], zt2[indC2], xt2[indD2], yt2[indD2], zt2[indD2],
                                      tol);

              if (match == 4) edge->_grade--;
            }
          }
        }
      }
    }
  }
  //blk2
  for (itrf2 = blk2->_listOfEdges.begin();
       itrf2 != blk2->_listOfEdges.end();
       itrf2++)
  {
    diff = 0;
    if ((*itrf2)->_blk1 == blk2 ) diff = (*itrf2)->_dirBlk1-edge->_dirBlk2;
    else if  ((*itrf2)->_blk2 == blk2 ) diff = (*itrf2)->_dirBlk2-edge->_dirBlk2;
    if ( diff != 0 && diff != 2 && diff != 4 && diff != 6)
    {
      //edge->_grade++;
      itrf1 = listOfMergeableEdges.begin();
      while ((itrf1 != listOfMergeableEdges.end())&&((*itrf1) != (*itrf2))) itrf1++;
      if ((*itrf1) == (*itrf2)) edge->_grade++;
    }
  }

  return;
}
//=============================================================================
void K_TRANSFORM::gradeEdge2D(BlockEdge* edge, list<Block*>& listOfBlks,
                              list<BlockEdge*>& listOfMergeableEdges,
                              E_Float tol)
{
  Block* blk1 = edge->_blk1;
  Block* blk2 = edge->_blk2;
  //edge->_grade = K_FUNC::E_max(blk1->_nmerging, blk2->_nmerging);

  // test si facette contient un pt exterieur
  if (edge->_externEdge == 1) edge->_grade +=2;
  E_Float dx1, dy1, dz1;
  E_Float xA, yA, zA, xB, yB, zB;
  E_Int diff = 1;
  list<Block*>::iterator itr2;
  list<BlockEdge*>::iterator itrf1;
  list<BlockEdge*>::iterator itrf2;
  list<BlockEdge*>::iterator itrf3;
  E_Int dirt2[4];
  E_Int indA2 =-1, indB2 =-1;
  dirt2[0] = 1; dirt2[1] =-1; dirt2[2] = 2; dirt2[3] =-2;
  E_Float xA2, yA2, zA2, xB2, yB2, zB2;

  //blk1
  for (itrf2 = blk1->_listOfEdges.begin(); itrf2 != blk1->_listOfEdges.end(); itrf2++)
  {
    diff = 0;
    if ((*itrf2)->_blk1 == blk1 ) diff = (*itrf2)->_dirBlk1-edge->_dirBlk1;
    else if  ((*itrf2)->_blk2 == blk1 ) diff = (*itrf2)->_dirBlk2-edge->_dirBlk1;
    if ( diff != 0 && diff != 2 && diff != 4)
    {
      //edge->_grade++;
      itrf1 = listOfMergeableEdges.begin();
      while ((itrf1 != listOfMergeableEdges.end())&&((*itrf1) != (*itrf2))) itrf1++;
      if ((*itrf1) == (*itrf2)) edge->_grade++;

      /* Recherche si creation d'une nouvelle facette mergeable */
      for (itrf3 = blk2->_listOfEdges.begin();
           itrf3 != blk2->_listOfEdges.end();
           itrf3++)
      {
        // recuperation des extremites de itrf2 et itrf3
        E_Float* ptA2 = (*itrf2)->_ptA; E_Float* ptB2 = (*itrf2)->_ptB;
        E_Float* ptA3 = (*itrf3)->_ptA; E_Float* ptB3 = (*itrf3)->_ptB;

        // sont ils coincidents en un point et lequel
        E_Int found = 0;
        //test A2 A3
        dx1 = ptA2[0]-ptA3[0]; dy1 = ptA2[1]-ptA3[1]; dz1 = ptA2[2]-ptA3[2];
        if (K_FUNC::E_abs(dx1) < tol && K_FUNC::E_abs(dy1) < tol && K_FUNC::E_abs(dz1) < tol)
        {
          xA = ptB2[0]; yA = ptB2[1]; zA = ptB2[2];
          xB = ptB3[0]; yB = ptB3[1]; zB = ptB3[2];
          found = 1; goto next;
        }
        // test A2 B3
        dx1 = ptA2[0]-ptB3[0]; dy1 = ptA2[1]-ptB3[1]; dz1 = ptA2[2]-ptB3[2];
        if (K_FUNC::E_abs(dx1) < tol && K_FUNC::E_abs(dy1) < tol && K_FUNC::E_abs(dz1) < tol)
        {
          xA = ptB2[0]; yA = ptB2[1]; zA = ptB2[2];
          xB = ptA3[0]; yB = ptA3[1]; zB = ptA3[2];
          found = 1; goto next;
        }
        // test B2A3
        dx1 = ptB2[0]-ptA3[0]; dy1 = ptB2[1]-ptA3[1]; dz1 = ptB2[2]-ptA3[2];
        if (K_FUNC::E_abs(dx1) < tol && K_FUNC::E_abs(dy1) < tol && K_FUNC::E_abs(dz1) < tol)
        {
          xA = ptA2[0]; yA = ptA2[1]; zA = ptA2[2];
          xB = ptB3[0]; yB = ptB3[1]; zB = ptB3[2];
          found = 1; goto next;
        }
        // test B2B3
        dx1 = ptB2[0]-ptB3[0]; dy1 = ptB2[1]-ptB3[1]; dz1 = ptB2[2]-ptB3[2];
        if (K_FUNC::E_abs(dx1) < tol && K_FUNC::E_abs(dy1) < tol && K_FUNC::E_abs(dz1) < tol)
        {
          xA = ptA2[0]; yA = ptA2[1]; zA = ptA2[2];
          xB = ptA3[0]; yB = ptA3[1]; zB = ptA3[2];
          found = 1; goto next;
        }
        // si oui : recherche d un bloc dont une arete est coincidente avec les 2 autres extremites
        next:;
        if ( found == 1 )
        {
          for (itr2 = listOfBlks.begin(); itr2 != listOfBlks.end(); itr2++)
          {
            E_Int ni2 = (*itr2)->_ni; E_Int nj2 = (*itr2)->_nj;
            FldArrayF* f2 = (*itr2)->_field;
            E_Float* xt2 = f2->begin(1);
            E_Float* yt2 = f2->begin(2);
            E_Float* zt2 = f2->begin(3);

            for (E_Int i2 = 0; i2 < 4; i2++)
            {
              E_Int dir2 = dirt2[i2];
              if ( dir2 ==-1){indA2 = 0; indB2 = (nj2-1)*ni2; }//imin
              else if ( dir2 == 1){indA2 = ni2-1; indB2 = ni2-1+ (nj2-1)*ni2; }//imax
              else if ( dir2 ==-2) {indA2 = 0; indB2 = ni2-1;}//jmin
              else if ( dir2 == 2) {indA2 = ni2-1+(nj2-1)*ni2; indB2 = (nj2-1)*ni2;}//jmax
              xA2 = xt2[indA2]; yA2 = yt2[indA2]; zA2 = zt2[indA2];
              xB2 = xt2[indB2]; yB2 = yt2[indB2]; zB2 = zt2[indB2];
              E_Int match = matchingEdges2D(xA, yA, zA, xB, yB, zB,
                                            xA2, yA2, zA2, xB2, yB2, zB2, tol);
              if ( match == 2 )
              {
                edge->_grade--;
              }
            }
          }
        }
      }
    }
  }

  //blk2
  for (itrf2 = blk2->_listOfEdges.begin();
       itrf2 != blk2->_listOfEdges.end();
       itrf2++)
  {
    diff = 0;
    if ((*itrf2)->_blk1 == blk2 ) diff = (*itrf2)->_dirBlk1-edge->_dirBlk2;
    else if  ((*itrf2)->_blk2 == blk2 ) diff = (*itrf2)->_dirBlk2-edge->_dirBlk2;
    if ( diff != 0 && diff != 2 && diff != 4)
    {
      //edge->_grade++;
      itrf1 = listOfMergeableEdges.begin();
      while ((itrf1 != listOfMergeableEdges.end())&&((*itrf1) != (*itrf2))) itrf1++;
      if ((*itrf1) == (*itrf2)) edge->_grade++;
    }
  }

}
//=============================================================================
E_Int K_TRANSFORM::checkNegativeVolumeCells(
  E_Int dim, E_Int im, E_Int jm, E_Int km,
  FldArrayF& coords)
{
  E_Int im1 = K_FUNC::E_max(im-1,1);
  E_Int jm1 = K_FUNC::E_max(jm-1,1);
  E_Int km1 = K_FUNC::E_max(km-1,1);
  E_Int ncells = im1*jm1*km1;
  E_Int ninti = im*jm1*km1;
  E_Int nintj = im1*jm*km1;
  E_Int nintk = im1*jm1*km;
  E_Int nint =  ninti + nintj + nintk;
  FldArrayF vol(ncells);
  FldArrayF surf(nint, 3);
  FldArrayF snorm(nint);
  FldArrayF centerInt(nint, 3);

  if (dim == 2)
    K_METRIC::compSurfStruct2D(
      im, jm, km,
      coords.begin(1), coords.begin(2), coords.begin(3), vol.begin());
  else
    K_METRIC::compMetricStruct(
      im, jm, km, ninti, nintj, nintk,
      coords.begin(1), coords.begin(2), coords.begin(3),
      vol.begin(), surf.begin(1), surf.begin(2), surf.begin(3),
      snorm.begin(),
      centerInt.begin(1), centerInt.begin(2), centerInt.begin(3));
  for (E_Int ind = 0; ind < ncells; ind++)
  {
    if (vol[ind] < 0.) return 0;
  }
  return 1;
}
