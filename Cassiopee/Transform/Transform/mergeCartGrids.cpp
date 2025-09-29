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
# include <list>
# include <stdlib.h>
#include "transform.h"
using namespace std;
using namespace K_FLD;
namespace K_TRANSFORM
{
  struct CartBlock;
  struct CartBlockFace
  {
      E_Int _i1, _i2, _j1, _j2, _k1, _k2;// definit la facette sur la grille
      E_Int _grade;// grade de la facette selon la methode Weakest Descent
      //E_Int _externFace;// la face a-t-elle une arete externe ?
      E_Float _xmin, _ymin, _zmin, _xmax,_ymax,_zmax;// bounding box de la facette
      CartBlock* _blk1;//bloc courant
      CartBlock* _blk2;//bloc oppose
  };
  struct CartBlock
  {
      E_Int _ni;
      E_Int _nj;
      E_Int _nk;//indices du bloc dans les trois directions
      E_Float _xmin, _ymin, _zmin, _xmax, _ymax, _zmax;
      E_Float _dh;// pas d'espace : pour definir le niveau immediatement
      list<CartBlockFace*> _listOfFaces;
  };
}

using namespace K_TRANSFORM;

//=============================================================================
/* Fusion des blocs cartesiens selon l algorithme de Rigby */
//=============================================================================
PyObject* K_TRANSFORM::mergeCartGrids(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  E_Int sizeMax;
  E_Float eps = 1.e-10;// tolerance match

  if (!PYPARSETUPLE_(args, O_ I_ R_,
                    &arrays, &sizeMax, &eps))
  {
      return NULL;
  }
  // Check every arrays
  if (PyList_Check(arrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "mergeCart: arrays argument must be a list.");
    return NULL;
  }
  // Extract infos from arrays
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt; vector<char*> eltType;
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
                    "mergeCart: arrays is not valid.");
    return NULL;
  }

  // coordonnees deja verifiees dans getFromArrays
  E_Int api = -1;
  E_Int nzones = structF.size();
  vector<E_Int> posxt(nzones); vector<E_Int> posyt(nzones); vector<E_Int> poszt(nzones);
  for (E_Int v = 0; v < nzones; v++)
  {
    char* varString = structVarString[v];
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;
    posxt[v] = posx; posyt[v] = posy; poszt[v] = posz;
    if (api == -1) api = structF[v]->getApi();
  }
  //printf("Running mergecart\n"); fflush(stdout);

  /* Determination de la bounding box des blocs */
  E_Float xmino, ymino, zmino, xmaxo, ymaxo, zmaxo;
  K_COMPGEOM::globalBoundingBox(posxt, posyt, poszt, structF,
                                xmino, ymino, zmino, xmaxo, ymaxo, zmaxo);
  list<CartBlock*> listOfBlocks;
  /* Construction de listOfBlocks */
  for (E_Int noz = 0; noz < nzones; noz++)
  {
    CartBlock* blk = new CartBlock;
    blk->_ni = nit[noz]; blk->_nj = njt[noz]; blk->_nk = nkt[noz];
    K_COMPGEOM::boundingBoxStruct(nit[noz], njt[noz], nkt[noz],
                                  structF[noz]->begin(posxt[noz]),
                                  structF[noz]->begin(posyt[noz]),
                                  structF[noz]->begin(poszt[noz]),
                                  blk->_xmin, blk->_ymin, blk->_zmin,
                                  blk->_xmax, blk->_ymax, blk->_zmax);

    blk->_dh = (blk->_xmax-blk->_xmin)/(nit[noz]-1);
    listOfBlocks.push_back(blk);
  }
  for (E_Int i = 0; i < nzones; i++) RELEASESHAREDS(objst[i], structF[i]);

  //printf("step1\n"); fflush(stdout);

  /* compte les facettes mergeables */
  /* construit la liste des facettes mergeables */
  list<CartBlock*>::iterator itr;
  list<CartBlock*>::iterator itr2;
  list<CartBlock*>::iterator itr3;

  vector<CartBlockFace*>::iterator vitrf;
  vector<CartBlockFace*>::iterator vitrf1;
  vector<CartBlockFace*>::iterator vitrf2;

  list<CartBlockFace*>::iterator itrf;
  list<CartBlockFace*>::iterator itrf1;
  list<CartBlockFace*>::iterator itrf2;
  list<CartBlockFace*>::iterator itrf3;
  list<CartBlockFace*> candidates;
  E_Int rank = 0; E_Int ncandidates = 0; E_LONG idum = -1;

  CartBlock *block1, *block2;

  vector<CartBlockFace*> listOfMergeableFace;  // liste des facettes mergeables (candidates au merge)
  CartBlockFace* face;
  E_Int gradeMin;
  CartBlockFace* faceMin;
  CartBlock* blockMerged;
  E_Int itrfconst, itrf2const, itrf3const;
  //E_Int fusi1, fusi2, fusj1, fusj2, fusk1, fusk2;
  E_Int size1, size2;
  E_Float dhfus, xfusmin, yfusmin, zfusmin, xfusmax, yfusmax, zfusmax;
  E_Int ninew, njnew, nknew;
  E_Float xmin1, ymin1, zmin1, xmax1, ymax1, zmax1, dh1;
  E_Float xmin2, ymin2, zmin2, xmax2, ymax2, zmax2, dh2;
  //E_Int ni1, nj1, nk1, ni2, nj2, nk2, dirp;
  //E_Int dir;
  E_Int nio, njo, nko;
  E_Float findex1, findex2;

  /* Premiere construction de la liste des faces fusionnables */
  for (itr = listOfBlocks.begin(); itr != listOfBlocks.end(); itr++)
  {
    //ni1 = (*itr)->_ni; nj1 = (*itr)->_nj; nk1 = (*itr)->_nk;
    dh1 = (*itr)->_dh;
    xmin1 = (*itr)->_xmin; xmax1 = (*itr)->_xmax;
    ymin1 = (*itr)->_ymin; ymax1 = (*itr)->_ymax;
    zmin1 = (*itr)->_zmin; zmax1 = (*itr)->_zmax;

    for (itr2 = itr; itr2 != listOfBlocks.end(); itr2++)
    {
      //ni2 = (*itr2)->_ni; nj2 = (*itr2)->_nj; nk2 = (*itr2)->_nk;
      dh2 = (*itr2)->_dh;
      xmin2 = (*itr2)->_xmin; xmax2 = (*itr2)->_xmax;
      ymin2 = (*itr2)->_ymin; ymax2 = (*itr2)->_ymax;
      zmin2 = (*itr2)->_zmin; zmax2 = (*itr2)->_zmax;

      if ((*itr) != (*itr2))
      {
        /* compare cas 1 */
        if(K_FUNC::fEqualZero(xmin1-xmin2, eps) &&
           K_FUNC::fEqualZero(xmax1-xmax2, eps) &&
           K_FUNC::fEqualZero(ymin1-ymin2, eps) &&
           K_FUNC::fEqualZero(ymax1-ymax2, eps) &&
           K_FUNC::fEqualZero(dh1-dh2, eps))
        {
          if (K_FUNC::fEqualZero(zmin1-zmax2, eps))//possible en 3D seulmt
          {
            /* cette facette est mergeable */
            face = new CartBlockFace;
            face->_i1 = 1; face->_i2 = (*itr)->_ni;
            face->_j1 = 1; face->_j2 = (*itr)->_nj;
            face->_k1 = 1; face->_k2 = 1;
            face->_xmin = xmin1; face->_ymin = ymin1; face->_zmin = zmin1;
            face->_xmax = xmax1; face->_ymax = ymax1; face->_zmax = zmin1;
            face->_grade = 0;
            face->_blk1 = (*itr); face->_blk2 = (*itr2);
            listOfMergeableFace.push_back(face);
            (*itr)->_listOfFaces.push_back(face);
            (*itr2)->_listOfFaces.push_back(face);
          }
          else if (K_FUNC::fEqualZero(zmax1-zmin2, eps))//possible en 3D seulmt
          {
            /* cette facette est mergeable */
            face = new CartBlockFace;
            face->_i1 = 1; face->_i2 = (*itr)->_ni;
            face->_j1 = 1; face->_j2 = (*itr)->_nj;
            face->_k1 = (*itr)->_nk; face->_k2 = (*itr)->_nk;
            face->_xmin = xmin1; face->_ymin = ymin1; face->_zmin = zmax1;
            face->_xmax = xmax1; face->_ymax = ymax1; face->_zmax = zmax1;
            face->_grade = 0;
            face->_blk1 = (*itr); face->_blk2 = (*itr2);
            listOfMergeableFace.push_back(face);
            (*itr)->_listOfFaces.push_back(face);
            (*itr2)->_listOfFaces.push_back(face);
          }
        }// cas 1

        /* compare cas 2 */
        else if (K_FUNC::fEqualZero(ymin1-ymin2, eps) &&
                 K_FUNC::fEqualZero(ymax1-ymax2, eps) &&
                 K_FUNC::fEqualZero(zmin1-zmin2, eps) && // Toujours vrai en 2D
                 K_FUNC::fEqualZero(zmax1-zmax2, eps) &&
                 K_FUNC::fEqualZero(dh1-dh2, eps))
        {
          if (K_FUNC::fEqualZero(xmin1-xmax2, eps))
          {
            /* cette facette est mergeable */
            face = new CartBlockFace;
            face->_i1 = 1; face->_i2 = 1;
            face->_j1 = 1; face->_j2 = (*itr)->_nj;
            face->_k1 = 1; face->_k2 = (*itr)->_nk;
            face->_xmin = xmin1; face->_ymin = ymin1; face->_zmin = zmin1;
            face->_xmax = xmin1; face->_ymax = ymax1; face->_zmax = zmax1;
            face->_grade = 0;
            face->_blk1 = (*itr);
            face->_blk2 = (*itr2);
            listOfMergeableFace.push_back(face);
            (*itr)->_listOfFaces.push_back(face);
            (*itr2)->_listOfFaces.push_back(face);
          }
          else if (K_FUNC::fEqualZero(xmax1-xmin2, eps))
          {
            /* cette facette est mergeable */
            face = new CartBlockFace;
            face->_i1 = (*itr)->_ni; face->_i2 = (*itr)->_ni;
            face->_j1 = 1; face->_j2 = (*itr)->_nj;
            face->_k1 = 1; face->_k2 = (*itr)->_nk;
            face->_xmin = xmax1; face->_ymin = ymin1; face->_zmin = zmin1;
            face->_xmax = xmax1; face->_ymax = ymax1; face->_zmax = zmax1;
            face->_grade = 0;
            face->_blk1 = (*itr); face->_blk2 = (*itr2);
            listOfMergeableFace.push_back(face);
            (*itr)->_listOfFaces.push_back(face);
            (*itr2)->_listOfFaces.push_back(face);
          }
        }// cas 2

        /* compare cas 3 */
        else if (K_FUNC::fEqualZero(xmin1-xmin2, eps) &&
                 K_FUNC::fEqualZero(xmax1-xmax2, eps) &&
                 K_FUNC::fEqualZero(zmin1-zmin2, eps) && // Toujours vrai en 2D
                 K_FUNC::fEqualZero(zmax1-zmax2, eps) &&
                 K_FUNC::fEqualZero(dh1-dh2,eps))
        {
          if (K_FUNC::fEqualZero(ymin1-ymax2, eps))
          {
            /* cette facette est mergeable */
            face = new CartBlockFace;
            face->_i1 = 1; face->_i2 = (*itr)->_ni;
            face->_j1 = 1; face->_j2 = 1;
            face->_k1 = 1; face->_k2 = (*itr)->_nk;
            face->_xmin = xmin1; face->_ymin = ymin1; face->_zmin = zmin1;
            face->_xmax = xmax1; face->_ymax = ymin1; face->_zmax = zmax1;
            face->_grade = 0;
            face->_blk1 = (*itr);
            face->_blk2 = (*itr2);
            listOfMergeableFace.push_back(face);
            (*itr)->_listOfFaces.push_back(face);
            (*itr2)->_listOfFaces.push_back(face);
          }
          else if (K_FUNC::fEqualZero(ymax1-ymin2, eps))
          {
            /* cette facette est mergeable */
            face = new CartBlockFace;
            face->_i1 = 1; face->_i2 = (*itr)->_ni;
            face->_j1 = (*itr)->_nj; face->_j2 = (*itr)->_nj;
            face->_k1 = 1; face->_k2 = (*itr)->_nk;
            face->_xmin = xmin1; face->_ymin = ymax1; face->_zmin = zmin1;
            face->_xmax = xmax1; face->_ymax = ymax1; face->_zmax = zmax1;
            face->_grade = 0;
            face->_blk1 = (*itr); face->_blk2 = (*itr2);
            listOfMergeableFace.push_back(face);
            (*itr)->_listOfFaces.push_back(face);
            (*itr2)->_listOfFaces.push_back(face);
          }
        }// cas 3

      } // test itr != itr2
    }// parcours itr2
  }// parcours itr

  //printf("step2\n"); fflush(stdout);
  CartBlockFace* mface;

  /* Recherche du grade de chaque facette candidate */
  for (vitrf = listOfMergeableFace.begin();
       vitrf != listOfMergeableFace.end();
       vitrf++)
  {
    mface = *vitrf;
    mface->_grade = 1;
    block1 = mface->_blk1; block2 = mface->_blk2;

    if (mface->_i1 == mface->_i2) itrfconst = 1;
    else if (mface->_j1 == mface->_j2) itrfconst = 2;
    else itrfconst = 3;

    // cette facette contient-elle une arete externe ?
    //if ( (*itrf)->_externFace == 1 ) (*itrf)->_grade = (*itrf)->_grade+2;

    /* block 1 */
    for (itrf2 = block1->_listOfFaces.begin();
         itrf2 != block1->_listOfFaces.end();
         itrf2++)
    {
      if ((*itrf2)->_i1 == (*itrf2)->_i2) { itrf2const = 1; findex1 = (*itrf2)->_xmin; }
      else if ((*itrf2)->_j1 == (*itrf2)->_j2) { itrf2const = 2; findex1 = (*itrf2)->_ymin; }
      else { itrf2const = 3; findex1 = (*itrf2)->_zmin; }

      if (itrfconst != itrf2const)
      {
        // Trouve *itrf2
        vitrf1 = listOfMergeableFace.begin();
        while ((vitrf1 != listOfMergeableFace.end())&&((*vitrf1) != (*itrf2))) vitrf1++; // !!
        if ((*vitrf1) == (*itrf2)) mface->_grade++;

        /* Recherche si creation d'une nouvelle facette mergeable */
        for (itrf3 = block2->_listOfFaces.begin();
             itrf3 != block2->_listOfFaces.end();
             itrf3++)
        {
          if ((*itrf3)->_i1 == (*itrf3)->_i2) { itrf3const = 1; findex2 = (*itrf3)->_xmin; }
          else if ((*itrf3)->_j1 == (*itrf3)->_j2) { itrf3const = 2; findex2 = (*itrf3)->_ymin; }
          else { itrf3const = 3; findex2 = (*itrf3)->_zmin; }

          if (itrf3const == itrf2const && K_FUNC::fEqualZero(findex1-findex2,eps))
          {
	    /*
            if (itrfconst == 1)
            {
              fusi1 = 1; fusi2 = (*itrf2)->_i2+(*itrf3)->_i2-1;
              fusj1 = 1; fusj2 = (*itrf2)->_j2;
              fusk1 = 1; fusk2 = (*itrf2)->_k2;
            }
            else if (itrfconst == 2)
            {
              fusi1 = 1; fusi2 = (*itrf2)->_i2;
              fusj1 = 1; fusj2 = (*itrf2)->_j2+(*itrf3)->_j2-1;
              fusk1 = 1; fusk2 = (*itrf2)->_k2;
            }
            else
            {
              fusi1 = 1; fusi2 = (*itrf2)->_i2;
              fusj1 = 1; fusj2 = (*itrf2)->_j2;
              fusk1 = 1; fusk2 = (*itrf2)->_k2+(*itrf3)->_k2-1;
            }
	    */
            dhfus = block1->_dh;
            xfusmin = K_FUNC::E_min(block1->_xmin,block2->_xmin);
            yfusmin = K_FUNC::E_min(block1->_ymin,block2->_ymin);
            zfusmin = K_FUNC::E_min(block1->_zmin,block2->_zmin);
            xfusmax = K_FUNC::E_max(block1->_xmax,block2->_xmax);
            yfusmax = K_FUNC::E_max(block1->_ymax,block2->_ymax);
            zfusmax = K_FUNC::E_max(block1->_zmax,block2->_zmax);

            for (itr2 = listOfBlocks.begin(); itr2 != listOfBlocks.end(); itr2++)
            {
              /* compare cas 1 */
              if (itrfconst != 3 &&
                  K_FUNC::fEqualZero(xfusmin-(*itr2)->_xmin, eps) &&
                  K_FUNC::fEqualZero(xfusmax-(*itr2)->_xmax, eps) &&
                  K_FUNC::fEqualZero(yfusmin-(*itr2)->_ymin, eps) &&
                  K_FUNC::fEqualZero(yfusmax-(*itr2)->_ymax, eps) &&
                  K_FUNC::fEqualZero(dhfus-(*itr2)->_dh,eps))
              {
                if (K_FUNC::fEqualZero(zfusmin-(*itr2)->_zmax, eps) ||
                    K_FUNC::fEqualZero(zfusmax-(*itr2)->_zmin, eps) )
                {
                  mface->_grade--;
                  mface->_grade--;
                }
              }
              /* compare cas 2 */
              if (itrfconst != 1 &&
                  K_FUNC::fEqualZero(yfusmin-(*itr2)->_ymin, eps) &&
                  K_FUNC::fEqualZero(yfusmax-(*itr2)->_ymax, eps) &&
                  K_FUNC::fEqualZero(zfusmin-(*itr2)->_zmin, eps) &&
                  K_FUNC::fEqualZero(zfusmax-(*itr2)->_zmax, eps) &&
                  K_FUNC::fEqualZero(dhfus-(*itr2)->_dh,eps))
              {
                if (K_FUNC::fEqualZero(xfusmin-(*itr2)->_xmax, eps) ||
                    K_FUNC::fEqualZero(xfusmax-(*itr2)->_xmin, eps) )
                {
                  mface->_grade--;
                  mface->_grade--;
                }
              }
              /* compare cas 3 */
              if (itrfconst != 2 &&
                  K_FUNC::fEqualZero(zfusmin-(*itr2)->_zmin, eps) &&
                  K_FUNC::fEqualZero(zfusmax-(*itr2)->_zmax, eps) &&
                  K_FUNC::fEqualZero(xfusmin-(*itr2)->_xmin, eps) &&
                  K_FUNC::fEqualZero(xfusmax-(*itr2)->_xmax, eps) &&
                  K_FUNC::fEqualZero(dhfus-(*itr2)->_dh,eps))
              {
                if (K_FUNC::fEqualZero(yfusmin-(*itr2)->_ymax, eps) ||
                    K_FUNC::fEqualZero(yfusmax-(*itr2)->_ymin, eps) )
                {
                  mface->_grade--;
                  mface->_grade--;
                }
              }
            }// fin boucle itr2 de listOfBlocks
          }// fin test itrf2const = itrf3const
        }// fin boucle itrf3 de blk2.listOfFaces
        /* fin de la recherche */
      }// fin itrfconst != itrf2const
    }

    for (itrf2 = block2->_listOfFaces.begin();
         itrf2 != block2->_listOfFaces.end();
         itrf2++)
    {
      if ((*itrf2)->_i1 == (*itrf2)->_i2) itrf2const = 1;
      else if ((*itrf2)->_j1 == (*itrf2)->_j2) itrf2const = 2;
      else itrf2const = 3;

      // Trouve itrf2
      vitrf1 = listOfMergeableFace.begin();
      while ((vitrf1 != listOfMergeableFace.end())&&((*vitrf1) != (*itrf2))) vitrf1++;
      if (itrfconst != itrf2const && (*vitrf1) == (*itrf2)) mface->_grade++;
    }
  } // boucle sur itrf

  if (listOfMergeableFace.size() <= 0) goto end;
  //printf("stepO\n"); fflush(stdout);

  truc:;

  /* Determine la facette de plus petit grade */
  gradeMin = 1000000; faceMin = NULL;
  for (vitrf = listOfMergeableFace.begin();
       vitrf != listOfMergeableFace.end(); vitrf++)
  {
    mface = *vitrf;
    if (mface->_grade <= gradeMin) gradeMin = mface->_grade;
  }
  for (vitrf = listOfMergeableFace.begin();
       vitrf != listOfMergeableFace.end(); vitrf++)
  { if ((*vitrf)->_grade == gradeMin) candidates.push_back(*vitrf); }

  ncandidates = candidates.size();
  rank = E_Int(K_NOISE::stdRand(&idum) * ncandidates);
  //printf("rank=%d %d %f\n", rank, ncandidates);
  ncandidates = 0;
  for (itrf = candidates.begin();
       itrf != candidates.end(); itrf++)
  {
    if (ncandidates == rank) { faceMin = *itrf; break; }
    ncandidates++;
  }
  candidates.clear();

  block1 = faceMin->_blk1; block2 = faceMin->_blk2;

  //printf("stepH\n"); fflush(stdout);

  if (K_FUNC::fEqualZero(faceMin->_xmin-faceMin->_xmax,eps))
  {
    ninew = block1->_ni+block2->_ni-1; njnew = block1->_nj; nknew = block1->_nk; //dir = 1;
  }
  else if (K_FUNC::fEqualZero(faceMin->_ymin-faceMin->_ymax,eps))
  {
    ninew = block1->_ni; njnew = block1->_nj+block2->_nj-1; nknew = block1->_nk; //dir = 2;
  }
  else
  {
    ninew = block1->_ni; njnew = block1->_nj; nknew = block1->_nk+block2->_nk-1; //dir = 3;
  }

  /* Assure que la grille creee ne soit pas trop grande */
  /* Utile pour le parallele */
  size1 = block1->_ni*block1->_nj*block1->_nk;
  size2 = block2->_ni*block2->_nj*block2->_nk;
  if (size1 + size2 > sizeMax)
  {
    /* On enleve la facette choisie de la liste des facettes mergeables */
    itrf = block1->_listOfFaces.begin();
    while (itrf != block1->_listOfFaces.end())
    {
      itrf3 = itrf;
      itrf3++;
      if ((*itrf) == faceMin) block1->_listOfFaces.erase(itrf);
      itrf = itrf3;
    }
    itrf = block2->_listOfFaces.begin();
    while (itrf != block2->_listOfFaces.end())
    {
      itrf3 = itrf;
      itrf3++;
      if ((*itrf) == faceMin) block2->_listOfFaces.erase(itrf);
      itrf = itrf3;
    }
    vitrf = listOfMergeableFace.begin();
    // Efface faceMin
    while ((vitrf != listOfMergeableFace.end())&&((*vitrf) != faceMin)) vitrf++; // !!
    listOfMergeableFace.erase(vitrf); // check
    delete faceMin;
    if (listOfMergeableFace.size() > 0) goto truc;
    else goto end;
  }

  //printf("la\n"); fflush(stdout);

  /* Tagger les facettes condamnees */
  /* les supprimer et les fusionner */
  /* cette fusion cree-t-elle de nouvelle facettes? */
  itrf = block1->_listOfFaces.begin();
  while (itrf != block1->_listOfFaces.end())
  {
    itrf3 = itrf; itrf3++;
    if ((*itrf) != faceMin)
    {
      // Efface *itrf
      vitrf2 = listOfMergeableFace.begin();
      while ((*vitrf2) != (*itrf)) vitrf2++; // !!
      listOfMergeableFace.erase(vitrf2); // check
      if ((*itrf)->_blk1 != block1)
      {
        itrf2 = (*itrf)->_blk1->_listOfFaces.begin();
        while ((*itrf2) != (*itrf)) itrf2++;
        (*itrf)->_blk1->_listOfFaces.erase(itrf2);
      }
      else
      {
        itrf2 = (*itrf)->_blk2->_listOfFaces.begin();
        while ((*itrf2) != (*itrf)) itrf2++;
        (*itrf)->_blk2->_listOfFaces.erase(itrf2);
      }
      delete (*itrf);
      block1->_listOfFaces.erase(itrf);
    }
    itrf = itrf3;
  }

  itrf = block2->_listOfFaces.begin();
  while (itrf != block2->_listOfFaces.end())
  {
    itrf3 = itrf; itrf3++;
    if ((*itrf) != faceMin)
    {
      // Efface *itrf
      vitrf2 = listOfMergeableFace.begin();
      while ((*vitrf2) != (*itrf)) vitrf2++; // !!
      listOfMergeableFace.erase(vitrf2); // check
      if ((*itrf)->_blk1 != block2)
      {
        itrf2 = (*itrf)->_blk1->_listOfFaces.begin();
        while ((*itrf2) != (*itrf)) itrf2++;
        (*itrf)->_blk1->_listOfFaces.erase(itrf2);
      }
      else
      {
        itrf2 = (*itrf)->_blk2->_listOfFaces.begin();
        while ((*itrf2) != (*itrf)) itrf2++;
        (*itrf)->_blk2->_listOfFaces.erase(itrf2);
      }
      delete (*itrf);
      block2->_listOfFaces.erase(itrf);
    }
    itrf = itrf3;
  }

  //printf("la la\n"); fflush(stdout);

  /* On enleve la facette choisie de la liste des facettes mergeables */
  // Efface faceMin
  vitrf = listOfMergeableFace.begin();
  while ((vitrf != listOfMergeableFace.end())&&((*vitrf) != faceMin)) vitrf++; // !!
  listOfMergeableFace.erase(vitrf); // check
  delete faceMin;

  /* Merger les blocks */
  /*
  if (dir == 1)
  {
    if (block1->_xmin < block2->_xmin) dirp = 1;
    else dirp =-1;
  }
  else if (dir == 2)
  {
    if (block1->_ymin < block2->_ymin) dirp = 2;
    else dirp = -2;
  }
  else
  {
    if (block1->_zmin < block2->_zmin) dirp = 3;
    else dirp = -3;
  }
  */
  block1->_ni = ninew; block1->_nj = njnew;  block1->_nk = nknew;
  block1->_listOfFaces.clear();
  block1->_xmin = K_FUNC::E_min(block1->_xmin, block2->_xmin);
  block1->_xmax = K_FUNC::E_max(block1->_xmax, block2->_xmax);
  block1->_ymin = K_FUNC::E_min(block1->_ymin, block2->_ymin);
  block1->_ymax = K_FUNC::E_max(block1->_ymax, block2->_ymax);
  block1->_zmin = K_FUNC::E_min(block1->_zmin, block2->_zmin);
  block1->_zmax = K_FUNC::E_max(block1->_zmax, block2->_zmax);

  // delete block2
  block2->_listOfFaces.clear();
  itr = listOfBlocks.begin();
  while ((*itr) != block2) itr++;
  listOfBlocks.erase(itr);
  delete block2;

  xmin1 = block1->_xmin; ymin1 = block1->_ymin; zmin1 = block1->_zmin;
  xmax1 = block1->_xmax; ymax1 = block1->_ymax; zmax1 = block1->_zmax;
  dh1 = block1->_dh;

  /* Trouver les nouvelles facettes mergeables pour le nouveau bloc */
  for (itr2 = listOfBlocks.begin(); itr2 != listOfBlocks.end(); itr2++)
  {
    if (block1 != (*itr2))
    {
      xmin2 = (*itr2)->_xmin; xmax2 = (*itr2)->_xmax;
      ymin2 = (*itr2)->_ymin; ymax2 = (*itr2)->_ymax;
      zmin2 = (*itr2)->_zmin; zmax2 = (*itr2)->_zmax;
      dh2 = (*itr2)->_dh;

      /* compare cas 1 */
      if( K_FUNC::fEqualZero(xmin1-xmin2, eps) &&
          K_FUNC::fEqualZero(xmax1-xmax2, eps) &&
          K_FUNC::fEqualZero(ymin1-ymin2, eps) &&
          K_FUNC::fEqualZero(ymax1-ymax2, eps) &&
          K_FUNC::fEqualZero(dh1-dh2,eps) )
      {
        if (K_FUNC::fEqualZero(zmin1-zmax2, eps) )//possible en 3D seulmt
        {
          /* cette facette est mergeable */
          face = new CartBlockFace;
          face->_i1 = 1; face->_i2 = block1->_ni;
          face->_j1 = 1; face->_j2 = block1->_nj;
          face->_k1 = 1; face->_k2 = 1;
          face->_xmin = xmin1; face->_ymin = ymin1; face->_zmin = zmin1;
          face->_xmax = xmax1; face->_ymax = ymax1; face->_zmax = zmin1;
          face->_grade = 0;
          face->_blk1 = block1; face->_blk2 = (*itr2);
          listOfMergeableFace.push_back(face);
          block1->_listOfFaces.push_back(face);
          (*itr2)->_listOfFaces.push_back(face);
        }
        else if (K_FUNC::fEqualZero(zmax1-zmin2, eps))//possible en 3D seulmt
        {
          /* cette facette est mergeable */
          face = new CartBlockFace;
          face->_i1 = 1; face->_i2 = block1->_ni;
          face->_j1 = 1; face->_j2 = block1->_nj;
          face->_k1 = block1->_nk; face->_k2 = block1->_nk;
          face->_xmin = xmin1; face->_ymin = ymin1; face->_zmin = zmax1;
          face->_xmax = xmax1; face->_ymax = ymax1; face->_zmax = zmax1;
          face->_grade = 0;
          face->_blk1 = block1; face->_blk2 = (*itr2);
          listOfMergeableFace.push_back(face);
          block1->_listOfFaces.push_back(face);
          (*itr2)->_listOfFaces.push_back(face);
        }
      }// fin cas 1

      /* compare cas 2 */
      else if (K_FUNC::fEqualZero(ymin1-ymin2, eps) &&
               K_FUNC::fEqualZero(ymax1-ymax2, eps) &&
               K_FUNC::fEqualZero(zmin1-zmin2, eps) && // Toujours vrai en 2D
               K_FUNC::fEqualZero(zmax1-zmax2, eps) &&
               K_FUNC::fEqualZero(dh1-dh2,eps))
      {
        if (K_FUNC::fEqualZero(xmin1-xmax2, eps))
        {
          /* cette facette est mergeable */
          face = new CartBlockFace;
          face->_i1 = 1; face->_i2 = 1;
          face->_j1 = 1; face->_j2 = block1->_nj;
          face->_k1 = 1; face->_k2 = block1->_nk;
          face->_xmin = xmin1; face->_ymin = ymin1; face->_zmin = zmin1;
          face->_xmax = xmin1; face->_ymax = ymax1; face->_zmax = zmax1;
          face->_grade = 0;
          face->_blk1 = block1; face->_blk2 = (*itr2);
          listOfMergeableFace.push_back(face);
          block1->_listOfFaces.push_back(face);
          (*itr2)->_listOfFaces.push_back(face);
        }
        else if (K_FUNC::fEqualZero(xmax1-xmin2, eps))
        {
          /* cette facette est mergeable */
          face = new CartBlockFace;
          face->_i1 = block1->_ni; face->_i2 = block1->_ni;
          face->_j1 = 1; face->_j2 = block1->_nj;
          face->_k1 = block1->_nk; face->_k2 = block1->_nk;
          face->_xmin = xmax1; face->_ymin = ymin1; face->_zmin = zmin1;
          face->_xmax = xmax1; face->_ymax = ymax1; face->_zmax = zmax1;
          face->_grade = 0;
          face->_blk1 = block1; face->_blk2 = (*itr2);
          listOfMergeableFace.push_back(face);
          block1->_listOfFaces.push_back(face);
          (*itr2)->_listOfFaces.push_back(face);
        }
      }// fin cas 2
      /* compare cas 3 */
      else if (K_FUNC::fEqualZero(xmin1-xmin2, eps) &&
               K_FUNC::fEqualZero(xmax1-xmax2, eps) &&
               K_FUNC::fEqualZero(zmin1-zmin2, eps) && // Toujours vrai en 2D
               K_FUNC::fEqualZero(zmax1-zmax2, eps) &&
               K_FUNC::fEqualZero(dh1-dh2,eps))
      {
        if (K_FUNC::fEqualZero(ymin1-ymax2, eps))
        {
          /* cette facette est mergeable */
          face = new CartBlockFace;
          face->_i1 = 1; face->_i2 = block1->_ni;
          face->_j1 = 1; face->_j2 = 1;
          face->_k1 = 1; face->_k2 = block1->_nk;
          face->_xmin = xmin1; face->_ymin = ymin1; face->_zmin = zmin1;
          face->_xmax = xmax1; face->_ymax = ymin1; face->_zmax = zmax1;
          face->_grade = 0;
          face->_blk1 = block1;
          face->_blk2 = (*itr2);
          listOfMergeableFace.push_back(face);
          block1->_listOfFaces.push_back(face);
          (*itr2)->_listOfFaces.push_back(face);
        }
        else if (K_FUNC::fEqualZero(ymax1-ymin2, eps))
        {
          /* cette facette est mergeable */
          face = new CartBlockFace;
          face->_i1 = 1; face->_i2 = block1->_ni;
          face->_j1 = block1->_nj; face->_j2 = block1->_nj;
          face->_k1 = 1; face->_k2 = block1->_nk;
          face->_xmin = xmin1; face->_ymin = ymax1; face->_zmin = zmin1;
          face->_xmax = xmax1; face->_ymax = ymax1; face->_zmax = zmax1;
          face->_grade = 0;
          face->_blk1 = block1; face->_blk2 = (*itr2);
          listOfMergeableFace.push_back(face);
          block1->_listOfFaces.push_back(face);
          (*itr2)->_listOfFaces.push_back(face);
        }
      }// cas 3
    }
  }

  //printf("step3\n"); fflush(stdout);
  if (listOfMergeableFace.size() <= 0) goto end;

  blockMerged = block1;
  /* Grader le nouveau bloc et ses voisins */
  for (itr2 = listOfBlocks.begin(); itr2 != listOfBlocks.end(); itr2++)
  {
    if (blockMerged != (*itr2))
    {
      xmin2 = (*itr2)->_xmin; xmax2 = (*itr2)->_xmax;
      ymin2 = (*itr2)->_ymin; ymax2 = (*itr2)->_ymax;
      zmin2 = (*itr2)->_zmin; zmax2 = (*itr2)->_zmax;
      dh2 = (*itr2)->_dh;
      if (K_FUNC::fEqualZero(dh1-dh2,eps))
      {
        if (((K_FUNC::fEqualZero(zmin1-zmax2, eps) ||
              K_FUNC::fEqualZero(zmax1-zmin2, eps)) &&
             xmin2 < xmax1 && xmax2 > xmin1 && ymin2 < ymax1 && ymax2 > ymin1) ||
            ((K_FUNC::fEqualZero(xmin1-xmax2, eps) ||
              K_FUNC::fEqualZero(xmax1-xmin2, eps)) &&
             ymin2 < ymax1 && ymax2 > ymin1 &&
             ((zmin2 < zmax1 && zmax2 > zmin1) ||                // 3D
              K_FUNC::fEqualZero(zmin1-zmax1, eps))) ||  // 2D
            ((K_FUNC::fEqualZero(ymin1-ymax2, eps) ||
              K_FUNC::fEqualZero(ymax1-ymin2, eps)) &&
             xmin2 < xmax1 && xmax2 > xmin1 &&
             ((zmin2 < zmax1 && zmax2 > zmin1) ||             // 3D
              K_FUNC::fEqualZero(zmin1-zmax1, eps)))) // 2D
        {
          for (itrf = (*itr2)->_listOfFaces.begin();
               itrf != (*itr2)->_listOfFaces.end();
               itrf++)
          {
            // Trouve *itrf
            vitrf1 = listOfMergeableFace.begin();
            while ((vitrf1 != listOfMergeableFace.end())&&((*vitrf1) != (*itrf))) vitrf1++; // !!
            if ((*vitrf1) == (*itrf))
            {
              (*itrf)->_grade = 1;
              block1 = (*itrf)->_blk1; block2 = (*itrf)->_blk2;
              if ((*itrf)->_i1 == (*itrf)->_i2) itrfconst = 1;
              else if ((*itrf)->_j1 == (*itrf)->_j2) itrfconst = 2;
              else itrfconst = 3;

              // cette facette contient-elle une arete externe ?
              //if ( (*itrf)->_externFace == 1 ) (*itrf)->_grade = (*itrf)->_grade+2;

              /* block 1 */
              for (itrf2 = block1->_listOfFaces.begin();
                   itrf2 != block1->_listOfFaces.end();
                   itrf2++)
              {
                if ((*itrf2)->_i1 == (*itrf2)->_i2) { itrf2const = 1; findex1 = (*itrf2)->_xmin; }
                else if ((*itrf2)->_j1 == (*itrf2)->_j2) { itrf2const = 2; findex1 = (*itrf2)->_ymin; }
                else {itrf2const = 3; findex1 = (*itrf2)->_zmin;}

                if (itrfconst != itrf2const)
                {
                  // Trouve *itrf2
                  vitrf1 = listOfMergeableFace.begin();
                  while ((vitrf1 != listOfMergeableFace.end())&&((*vitrf1) != (*itrf2))) vitrf1++; // !!
                  if ((*vitrf1) == (*itrf2)) (*itrf)->_grade++;

                  /* Recherche si creation d'une nouvelle facette mergeable */
                  for (itrf3 = block2->_listOfFaces.begin(); itrf3 != block2->_listOfFaces.end();itrf3++)
                  {
                    if ((*itrf3)->_i1 == (*itrf3)->_i2) { itrf3const = 1; findex2 = (*itrf3)->_xmin; }
                    else if ((*itrf3)->_j1 == (*itrf3)->_j2) { itrf3const = 2; findex2 = (*itrf3)->_ymin; }
                    else { itrf3const = 3; findex2 = (*itrf3)->_zmin; }

                    if (itrf3const == itrf2const && K_FUNC::fEqualZero(findex1-findex2,eps)== true)
                    {
		      /*
                      if (itrfconst == 1)
                      {
                        fusi1 = 1; fusi2 = (*itrf2)->_i2+(*itrf3)->_i2-1;
                        fusj1 = 1; fusj2 = (*itrf2)->_j2;
                        fusk1 = 1; fusk2 = (*itrf2)->_k2;
                      }
                      else if (itrfconst == 2)
                      {
                        fusi1 = 1; fusi2 = (*itrf2)->_i2;
                        fusj1 = 1; fusj2 = (*itrf2)->_j2+(*itrf3)->_j2-1;
                        fusk1 = 1; fusk2 = (*itrf2)->_k2;
                      }
                      else
                      {
                        fusi1 = 1; fusi2 = (*itrf2)->_i2;
                        fusj1 = 1; fusj2 = (*itrf2)->_j2;
                        fusk1 = 1; fusk2 = (*itrf2)->_k2+(*itrf3)->_k2-1;
                      }
		      */
                      dhfus = block1->_dh;
                      xfusmin = K_FUNC::E_min(block1->_xmin,block2->_xmin);
                      yfusmin = K_FUNC::E_min(block1->_ymin,block2->_ymin);
                      zfusmin = K_FUNC::E_min(block1->_zmin,block2->_zmin);
                      xfusmax = K_FUNC::E_max(block1->_xmax,block2->_xmax);
                      yfusmax = K_FUNC::E_max(block1->_ymax,block2->_ymax);
                      zfusmax = K_FUNC::E_max(block1->_zmax,block2->_zmax);

                      for (itr3 = listOfBlocks.begin();
                           itr3 != listOfBlocks.end();
                           itr3++)
                      {
                        /* compare cas 1 */
                        if (itrfconst != 3 &&
                            K_FUNC::fEqualZero(xfusmin-(*itr3)->_xmin, eps) &&
                            K_FUNC::fEqualZero(xfusmax-(*itr3)->_xmax, eps) &&
                            K_FUNC::fEqualZero(yfusmin-(*itr3)->_ymin, eps) &&
                            K_FUNC::fEqualZero(yfusmax-(*itr3)->_ymax, eps) &&
                            K_FUNC::fEqualZero(dhfus-(*itr3)->_dh,eps))
                        {
                          if (K_FUNC::fEqualZero(zfusmin-(*itr3)->_zmax, eps) ||
                              K_FUNC::fEqualZero(zfusmax-(*itr3)->_zmin, eps) )
                          {
                            (*itrf)->_grade--;
                            (*itrf)->_grade--;
                          }
                        }
                        /* compare cas 2 */
                        if (itrfconst != 1 &&
                            K_FUNC::fEqualZero(yfusmin-(*itr3)->_ymin, eps) &&
                            K_FUNC::fEqualZero(yfusmax-(*itr3)->_ymax, eps) &&
                            K_FUNC::fEqualZero(zfusmin-(*itr3)->_zmin, eps) &&
                            K_FUNC::fEqualZero(zfusmax-(*itr3)->_zmax, eps) &&
                            K_FUNC::fEqualZero(dhfus-(*itr3)->_dh,eps))
                        {
                          if (K_FUNC::fEqualZero(xfusmin-(*itr3)->_xmax, eps) ||
                              K_FUNC::fEqualZero(xfusmax-(*itr3)->_xmin, eps) )
                          {
                            (*itrf)->_grade--;
                            (*itrf)->_grade--;
                          }
                        }
                        /* compare cas 3 */
                        if (itrfconst != 2 &&
                            K_FUNC::fEqualZero(zfusmin-(*itr3)->_zmin, eps) &&
                            K_FUNC::fEqualZero(zfusmax-(*itr3)->_zmax, eps) &&
                            K_FUNC::fEqualZero(xfusmin-(*itr3)->_xmin, eps) &&
                            K_FUNC::fEqualZero(xfusmax-(*itr3)->_xmax, eps) &&
                            K_FUNC::fEqualZero(dhfus-(*itr3)->_dh,eps))
                        {
                          if (K_FUNC::fEqualZero(yfusmin-(*itr3)->_ymax, eps) ||
                              K_FUNC::fEqualZero(yfusmax-(*itr3)->_ymin, eps) )
                          {
                            (*itrf)->_grade--;
                            (*itrf)->_grade--;
                          }
                        }
                      }// for itr3
                    }
                  }
                  /* fin de la recherche */
                }// fin itrfconst != itrf2const
              }
              for (itrf2 = block2->_listOfFaces.begin();
                   itrf2 != block2->_listOfFaces.end();
                   itrf2++)
              {
                if ((*itrf2)->_i1 == (*itrf2)->_i2) itrf2const = 1;
                else if ((*itrf2)->_j1 == (*itrf2)->_j2) itrf2const = 2;
                else itrf2const = 3;

                // Trouve *itrf2
                vitrf1 = listOfMergeableFace.begin();
                while ((vitrf1 != listOfMergeableFace.end())&&((*vitrf1) != (*itrf2))) vitrf1++; // !!
                if (itrfconst != itrf2const && (*vitrf1) == (*itrf2)) (*itrf)->_grade++;
              }
            }//fin if itrf face mergeable
          }// boucle sur itrf de itr2
        }//fin if itr2 voisin de blockMerged
      }//dh1 == dh2
    }//fin if (blockMerged != (*itr2))
  }// boucle sur itr2 des blocks

  //printf("size of mergeable %d\n", listOfMergeableFace.size()); fflush(stdout);
  if (listOfMergeableFace.size() > 0) goto truc;
  end:;

  PyObject* l = PyList_New(0);
  PyObject* tpl;
  FldArrayF* coordso;
  E_Float *xto, *yto, *zto;
  E_Int nionjo;
  E_Float xmin, ymin, zmin, dh;
  for (itr = listOfBlocks.begin(); itr != listOfBlocks.end(); itr++)
  {
    nio = (*itr)->_ni; njo = (*itr)->_nj; nko = (*itr)->_nk;
    tpl = K_ARRAY::buildArray3(3, structVarString[0], nio, njo, nko, api);
    K_ARRAY::getFromArray3(tpl, coordso);
    xto = coordso->begin(1);
    yto = coordso->begin(2);
    zto = coordso->begin(3);
    xmin = (*itr)->_xmin;
    ymin = (*itr)->_ymin;
    zmin = (*itr)->_zmin;
    dh = (*itr)->_dh;
    nionjo = nio*njo;
    #pragma omp parallel
    {
      E_Int ind;
      #pragma omp for collapse(3)
      for (E_Int k = 0; k < nko; k++)
      for (E_Int j = 0; j < njo; j++)
      for (E_Int i = 0; i < nio; i++)
      {
        ind = i + j * nio + k * nionjo;
        xto[ind] = xmin + i * dh;
        yto[ind] = ymin + j * dh;
        zto[ind] = zmin + k * dh;
      }
    }
    RELEASESHAREDS(tpl, coordso);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }

  // nettoyages...
  for (itr = listOfBlocks.begin(); itr != listOfBlocks.end(); itr++)
  {
    list<CartBlockFace*>& listOfFaces = (*itr)->_listOfFaces;
    for (itrf = listOfFaces.begin(); itrf != listOfFaces.end(); itrf++)
      delete *itrf;
    delete *itr;
  }
  listOfBlocks.clear();

  return l;
}
