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
# include "stdio.h"
# include "transform.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Redresseur de normales : etant donnee une liste de blocs paroi, oriente les
   blocs de telle sorte que les normales soient toutes dans le meme sens
   L element de reference est le premier element de la liste
*/
//=============================================================================
PyObject* K_TRANSFORM::reorderAll(PyObject* self, PyObject* args)
{
  // Load block arrays
  PyObject* listBlks;
  E_Int dir=1; // direction of the normals
  if (!PYPARSETUPLE_(args, O_ I_, &listBlks, &dir))
  {
    return NULL;
  }
  // Check dir
  if (dir != 1 && dir != -1)
  {
    printf("Warning: reorderAll: direction is invalid. Set dir = 1.\n");
    dir = 1;
  }

  // Check every array in listBlks
  if (PyList_Check(listBlks) == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "reorderAll: argument must be a list.");
    return NULL;
  }

  E_Int nzone = PyList_Size(listBlks);
  vector<FldArrayF*> vectOfFields;
  vector<E_Int> nis;// initial ni, nj, nk
  vector<E_Int> njs;
  vector<E_Int> nks;
  vector<E_Int> nit;// modified ni, nj, nk to set nk = 1
  vector<E_Int> njt;
  vector<E_Int> nkt;
  vector<E_Int> posxt;
  vector<E_Int> posyt;
  vector<E_Int> poszt;
  E_Int posx, posy, posz;
  char* varString;
  E_Int api = -1;
  FldArrayF bbox(nzone, 6);

  // Extraction des infos pour chaque bloc
  for (int i = 0; i < nzone; i++)
  {
    E_Int nil, njl, nkl;
    PyObject* tpl = PyList_GetItem(listBlks, i);

    char* eltType;
    FldArrayF* f; FldArrayI* cn;
    E_Int res =
      K_ARRAY::getFromArray3(tpl, varString, f, nil, njl, nkl, cn, eltType);

    nis.push_back(nil); njs.push_back(njl); nks.push_back(nkl);
    if (api == -1) api = f->getApi();

    if (res != 1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "reorderAll: array is not structured.");
      return NULL;
    }
    if (nil == 1)
    {
      nil = njl;
      njl = nkl;
    }
    else if (njl == 1)
    {
      njl = nkl;
      nkl = 1;
    }
    else if (nkl == 1)
    {
      ;
    }
    else
    {
      PyErr_SetString(PyExc_TypeError,
                      "reorderAll: must be a surface." );
      return NULL;
    }
    //check if coordinates exist
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);

    if (posx == -1 || posy == -1 || posz == -1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "reorderAll: coordinates not found.");
      return NULL;
    }
    posx++; posy++; posz++;
    nit.push_back(nil);
    njt.push_back(njl);
    posxt.push_back(posx);
    posyt.push_back(posy);
    poszt.push_back(posz);

    //calcul des bounding boxes
    K_COMPGEOM::boundingBoxStruct(nis[i], njs[i], nks[i],
                                  f->begin(posx), f->begin(posy), f->begin(posz),
                                  bbox(i,1), bbox(i,2), bbox(i,3),
                                  bbox(i,4), bbox(i,5), bbox(i,6));

    vectOfFields.push_back(f);
  }//parcours de toutes les zones

  if (api == -1) api = 1;

  E_Int vectOfFieldsSize = vectOfFields.size();
  if (vectOfFieldsSize != nzone )
  {
    PyErr_SetString(PyExc_TypeError,
                    "reorderAll: wrong number of fields.");
    return NULL;
  }

  /*---------------------------------------------*/
  /* 1-Calcul de la distance dij entre les blocs */
  /*---------------------------------------------*/
  FldArrayF distMat(nzone, nzone);
  distMat.setAllValuesAt(K_CONST::E_MAX_FLOAT);
  FldArrayI indMat(nzone, nzone);
  indMat.setAllValuesAt(-1);

  E_Float dmin;
  E_Int ni1, nj1, ni2, nj2;
  E_Int ind1, ind2;
  E_Int posx1, posy1, posz1, posx2, posy2, posz2;

  for (E_Int v1 = 0; v1 < nzone; v1++)
  {
    ni1 = nit[v1];
    nj1 = njt[v1];
    posx1 = posxt[v1];
    posy1 = posyt[v1];
    posz1 = poszt[v1];
    FldArrayF* field1 = vectOfFields[v1];

    for (E_Int v2 = 0; v2 < nzone; v2++)
    {
      if ( v1 != v2 && indMat(v1, v2+1) == -1 )
      {
        E_Int intersect = testBBIntersection(v1, v2, bbox);
        if ( intersect == 1 )
        {
          ni2 = nit[v2];
          nj2 = njt[v2];
          posx2 = posxt[v2];
          posy2 = posyt[v2];
          posz2 = poszt[v2];
          FldArrayF* field2 = vectOfFields[v2];

          K_COMPGEOM::compMeanDist(ni1, nj1, 
          field1->begin(posx1), field1->begin(posy1), field1->begin(posz1),
          ni2, nj2, 
          field2->begin(posx2), field2->begin(posy2), field2->begin(posz2),
          ind1, ind2, dmin
          );

          distMat(v1,v2+1) = dmin;
          distMat(v2,v1+1) = dmin;
          indMat(v1,v2+1) = ind1;
          indMat(v2,v1+1) = ind2;
        }
      }
    } // parcours des blocs v2
  }// parcours des blocs v1

  // tri des blocs relies selon la distance minimale

  vector<FldArrayI*> listOfSortedBlks;
  sortBlocks(distMat, listOfSortedBlks);

  /*--------------------------------------------*/
  /* Algorithme de parcours du graphe des blocs */
  /*--------------------------------------------*/
  E_Int isOk = 0;
  FldArrayI rel(nzone, nzone);
  FldArrayI tagOpp(nzone);  //determine les blocs de sens oppose
  FldArrayI dejaVu(nzone);  //tableau de travail dejaVu
  E_Float eps = 1.e-12;

  E_Int max = 10;
  E_Int cnt = 0;
  while ( isOk == 0 && cnt < max )
  {
    // Determination des relations entre blocs
    rel.setAllValuesAt(-1);

    for (E_Int v1 = 0; v1 < nzone; v1++)
      for (E_Int v2 = 0; v2 < nzone; v2++)
      {
        if ( v1 != v2 && rel(v1, v2+1) == -1 )
        {
          if ( distMat(v1, v2+1) < eps  )
          {
            rel(v1,v2+1) = indMat(v1,v2+1);
            rel(v2,v1+1) = indMat(v2,v1+1);
          }
        }
      }

    tagOpp.setAllValuesAt(0);
    tagOpp[0] = 1;
    dejaVu.setAllValuesAtNull();

    // Parcours du graphe pour determiner le sens des normales
    graphPathSearch(-1, 0, vectOfFields, nit, njt,
                    posxt, posyt, poszt, listOfSortedBlks,
                    rel, distMat, dejaVu, tagOpp);

    isOk = 1;
    for (E_Int i = 0; i < nzone; i++ )
    {
      if (dejaVu[i] == 0 )
      {
        isOk = 0;
        eps = eps * 10;
        cnt++;
        break;
      }
    }
  }
  if ( cnt == max )
  {
    printf("Warning: reorderAll: can't find a single path to link all the blocks...\n");
    printf("...Check if some blocks are not in contact.\n");
  }
  else
  {
    /*--------------------------------------------*/
    /* redressage des blocs ayant un tagOpp = -1 */
    /*--------------------------------------------*/
    for (E_Int v = 0; v < nzone; v++)
    {
      if (tagOpp[v] == -1 && dir == 1)
        reorder( nit[v], njt[v], *vectOfFields[v]);
      if (tagOpp[v] == 1 && dir == -1)
      reorder( nit[v], njt[v], *vectOfFields[v]);
    }
  }
  /*--------------*/
  /* build arrays */
  /*--------------*/
  PyObject* tpl;
  PyObject* l = PyList_New(0);

  for (E_Int i = 0; i < nzone; i++)
  {
    tpl = K_ARRAY::buildArray3(*vectOfFields[i], varString,
                               nis[i], njs[i], nks[i], api);
    RELEASESHAREDS(PyList_GetItem(listBlks, i), vectOfFields[i]);
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }


  E_Int sortedBlksSize = listOfSortedBlks.size();
  for (E_Int v = 0 ; v < sortedBlksSize; v++)
    delete listOfSortedBlks[v];

  return l;
}
//=========================================================================
/* Algorithme recursif de parcours de graphe : Depth First Search
   in : nb_father :  no du bloc pere
   in : nb_cur : no du bloc courant
   in : rel : relation entre les blocs : distance < eps
   rel(i,j+1) : fournit le numero du noeud pere i verifiant cette relation
   rel(j,i+1) : fournit le numero du noeud fils j verifiant cette relation
   in : distMat : distance entre les points les plus proches entre 2 blocs
   in/out : dejaVu
   in/out : tagOpp : sens de la normale par rapport a la reference
*/
//=========================================================================
void K_TRANSFORM::graphPathSearch(E_Int nb_father, E_Int nb_cur,
                                  vector<FldArrayF*>& listOfFields,
                                  vector<E_Int>& nit, vector<E_Int>& njt,
                                  vector<E_Int>& posxt, vector<E_Int>& posyt,
                                  vector<E_Int>& poszt,
                                  vector<FldArrayI*>& listOfSortedBlks,
                                  FldArrayI& rel, FldArrayF& distMat,
                                  FldArrayI& dejaVu,
                                  FldArrayI& tagOpp)
{
  E_Int isOpp;
  //E_Int notvalid = 0;

  if ( nb_father != -1 )
  {
    E_Int ni1 = nit[nb_father];
    E_Int nj1 = njt[nb_father];
    E_Int ni2 = nit[nb_cur];
    E_Int nj2 = njt[nb_cur];
    FldArrayF* field1 = listOfFields[nb_father];
    FldArrayF* field2 = listOfFields[nb_cur];
    E_Int posx1 = posxt[nb_father];
    E_Int posy1 = posyt[nb_father];
    E_Int posz1 = poszt[nb_father];
    E_Int posx2 = posxt[nb_cur];
    E_Int posy2 = posyt[nb_cur];
    E_Int posz2 = poszt[nb_cur];
    E_Int ind1 = rel(nb_father, nb_cur+1);
    E_Int ind2 = rel(nb_cur, nb_father+1);
    E_Float dist = distMat(nb_father, nb_cur+1);

    K_COMPGEOM::rectifyNormals(ni1, nj1, ind1, 
      field1->begin(posx1), field1->begin(posy1), field1->begin(posz1),
      ni2, nj2, ind2, 
      field2->begin(posx2), field2->begin(posy2), field2->begin(posz2),
      dist, isOpp
    );

    tagOpp[nb_cur] = isOpp * tagOpp[nb_father];
  }
  // marquage du noeud
  dejaVu[nb_cur] = 1;
//  if ( notvalid == 1 )
//    dejaVu[nb_cur] = 0;

  FldArrayI& listOfBlks = *listOfSortedBlks[nb_cur];

  for (E_Int i0 = 0; i0 < listOfBlks.getSize(); i0++)
  {
    E_Int i = listOfBlks[i0];
    //parcours des voisins non deja traites du noeud courant
    if ( dejaVu[i] == 0 && rel(nb_cur, i+1) != -1 )
    {
      graphPathSearch(nb_cur, i,  listOfFields, nit, njt,
                      posxt, posyt, poszt, listOfSortedBlks,
                      rel, distMat, dejaVu, tagOpp);
    }
  }
}
//=========================================================================
/* Reorder the numerotation of mesh */
//=========================================================================
void K_TRANSFORM::reorder( const E_Int im, const E_Int jm,
                           FldArrayF& field)
{
  E_Int nfld = field.getNfld();
  FldArrayF fnew(field.getSize(), nfld);

  // reordering
  E_Int beta;
  E_Int ind, ind2;
  for (E_Int i = 0; i < im; i++)
    for (E_Int j = 0; j < jm; j++)
      {
        beta = j*im;
        ind = i + beta;
        ind2 = im-i-1 + beta;
        for (E_Int n = 1; n <= nfld; n++)
          fnew(ind2,n) = field(ind,n);

      }
  field = fnew;
}

//==========================================================================
/* Return 1 if bbox of blocks intersect, 0 elsewhere */
//==========================================================================
E_Int K_TRANSFORM::testBBIntersection(E_Int noblk1, E_Int noblk2,
                                      FldArrayF& bbox)
{
  E_Float tol = 1.e-5;
  E_Float xmax1, ymax1, zmax1, xmin1, ymin1, zmin1;
  E_Float xmin2, xmax2, ymin2, ymax2, zmin2, zmax2;

  xmin1 = bbox(noblk1,1);
  ymin1 = bbox(noblk1,2);
  zmin1 = bbox(noblk1,3);
  xmax1 = bbox(noblk1,4);
  ymax1 = bbox(noblk1,5);
  zmax1 = bbox(noblk1,6);

  xmin2 = bbox(noblk2,1);
  ymin2 = bbox(noblk2,2);
  zmin2 = bbox(noblk2,3);
  xmax2 = bbox(noblk2,4);
  ymax2 = bbox(noblk2,5);
  zmax2 = bbox(noblk2,6);

  if ( xmin1  <=  xmax2 + tol && xmax1  >=  xmin2 -tol &&
       ymin1  <=  ymax2 + tol && ymax1  >=  ymin2 -tol &&
       zmin1  <=  zmax2 + tol && zmax1  >=  zmin2 -tol)
    return 1;

  else
    return 0;
}


//=========================================================================
/* Tri les blocs etant en relation selon leur distance minimale
   Rmq : la distance est recopiee car modifiee dans l algorithme de tri */
//=========================================================================
void K_TRANSFORM::sortBlocks(FldArrayF& distMat,
                             vector<FldArrayI*>& listOfSortedBlks)
{
  E_Int nzone = distMat.getSize();

  FldArrayI tmpBlks(nzone);
  FldArrayF tmpDist(nzone);

  for (E_Int i = 0; i < nzone; i++)
  {
    E_Int cnt = 0;
    for (E_Int j = 0; j < nzone; j++)
    {
      if (distMat(i,j+1) < K_CONST::E_MAX_FLOAT && i != j)
      {
        tmpBlks[cnt] = j;
        tmpDist[cnt] = distMat(i,j+1);
        cnt++;
      }
    }

    /* tri des blocs en relation selon leur distance */
    E_Int bg = 0;
    E_Int bd = cnt-1;
    quickSortBlks(bg, bd, tmpBlks, tmpDist);

    FldArrayI* sortedBlks = new FldArrayI(cnt);
    for (E_Int j = 0; j < cnt; j++)
      (*sortedBlks)[j] = tmpBlks[j];

    listOfSortedBlks.push_back(sortedBlks);
  }
}


//=============================================================================
/* Sort the array of coordinates (x,y,z) with respect to the norm array,
   from bound bg to bound bd  */
//=============================================================================
void K_TRANSFORM::quickSortBlks(E_Int bg, E_Int bd, FldArrayI& sortedBlks,
                                FldArrayF& dist)
{
  if ( bg < bd )
  {
    E_Int ipiv = pivotingBlks(bg, bd, sortedBlks, dist);
    quickSortBlks(bg, ipiv-1, sortedBlks, dist);
    quickSortBlks(ipiv+1, bd, sortedBlks, dist);
  }
}
//=============================================================================
/* Make a splitting of the list of blks between bg and bd and return the pivot*/
//=============================================================================
E_Int K_TRANSFORM::pivotingBlks(E_Int bg, E_Int bd, FldArrayI& sortedBlks,
                                FldArrayF& dist)
{
  E_Float piv = dist[bg];
  E_Int l = bg+1;
  E_Int r = bd;

  E_Float tmp;
  E_Int tmpi;

  while ( l <= r)
  {
    while ( l <= bd && dist[l] <= piv) l++;
    while ( r >=  0 && dist[r] > piv) r--;
    if ( l < r )
    {
      // swap norm(indr) and norm(indl)
      tmp = dist[r];
      dist[r] = dist[l];
      dist[l] = tmp;

      // swap indl and indr blks
      tmpi = sortedBlks[r];
      sortedBlks[r] = sortedBlks[l];
      sortedBlks[l] = tmpi;

      r--;
      l++;
    }
  }

  // swap norm(indr) and norm(bg)
  tmp = dist[r];
  dist[r] = dist[bg];
  dist[bg] = tmp;

  // swap indl and bg blks
  tmpi = sortedBlks[r];
  sortedBlks[r] =  sortedBlks[bg];
  sortedBlks[bg] = tmpi;

  return r;
}
