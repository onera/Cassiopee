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

// Gather degenerated windows
# include "connector.h"

using namespace K_FUNC;
using namespace std;
using namespace K_FLD;

//============================================================================
/* Gather degenerate windows */
//=============================================================================
PyObject* K_CONNECTOR::gatherDegenerated(PyObject* self, PyObject* args)
{
  PyObject *listOfAllTags, *listOfBlks, *listOfWinTypes;
  PyObject *listOfNI, *listOfNJ, *listOfNK;
  E_Int dimPb;
  if (!PYPARSETUPLEI(args,
                    "OOOOOOl", "OOOOOOi",
                    &listOfAllTags, &listOfWinTypes, &listOfBlks, &listOfNI, &listOfNJ, &listOfNK, &dimPb))
  {
      return NULL;
  }
  if (PyList_Check(listOfAllTags) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "gatherDegenerated: 1st argument must be a list.");
    return NULL;
  }
  if (PyList_Check(listOfBlks) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "gatherDegenerated: 2nd argument must be a list.");
    return NULL;
  }
  if (PyList_Check(listOfWinTypes) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "gatherDegenerated: 3rd argument must be a list.");
    return NULL;
  }
  if (PyList_Check(listOfNI) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "gatherDegenerated: 4th argument must be a list.");
    return NULL;
  }
  if (PyList_Check(listOfNJ) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "gatherDegenerated: 5th argument must be a list.");
    return NULL;
  }
  if (PyList_Check(listOfNK) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "gatherDegenerated: 6th argument must be a list.");
    return NULL;
  }
  /* Extract info on tag1 and tag2 */
  vector<E_Int> rest;
  vector<char*> structVarStringt; vector<char*> unstrVarStringt;
  vector<FldArrayF*> structFt; vector<FldArrayF*> unstrFt;
  vector<E_Int> nitt; vector<E_Int> njtt; vector<E_Int> nktt;
  vector<FldArrayI*> cntt;
  vector<char*> eltTypett;
  vector<PyObject*> objstt, objutt;
  E_Boolean skipNoCoord = false;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = true;
  E_Boolean skipDiffVars = true;

  E_Int isOk = K_ARRAY::getFromArrays(
    listOfAllTags, rest, structVarStringt, unstrVarStringt,
    structFt, unstrFt, nitt, njtt, nktt, cntt, eltTypett, objstt, objutt, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nwinstot = structFt.size();
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError, "gatherDegenerated: 2nd list of arrays is not valid.");
    for (E_Int is = 0; is < nwinstot; is++)
    {
      RELEASESHAREDS(objstt[is], structFt[is]);
    }
    return NULL;
  } 

  vector<E_Int> post1; 
  E_Int postag1;
  for (E_Int i = 0; i < nwinstot; i++)
  {   
    postag1 = K_ARRAY::isNamePresent("tag1", structVarStringt[i]); postag1++;
    if (postag1 < 1) 
    {
      PyErr_SetString(PyExc_TypeError, "gatherDegenerated: tag1 must be defined in listOfAllWins.");
      for (E_Int is = 0; is < nwinstot; is++)
      {
        RELEASESHAREDS(objstt[is], structFt[is]);
      }
      return NULL;
    } 
    post1.push_back(postag1); 
  }
  /* Extract window type and corresponding original block */
  E_Int nwins2 = PyList_Size(listOfWinTypes);
  if (nwins2 != nwinstot)
  {
    PyErr_SetString(PyExc_TypeError, "gatherDegenerated: listOfAllWins and listOfWinTypes must be of same size.");
    for (E_Int is = 0; is < nwinstot; is++)
    {
      RELEASESHAREDS(objstt[is], structFt[is]);
    }
    return NULL;
  }
  E_Int nwins3 = PyList_Size(listOfBlks);
  if (nwins3 != nwinstot)
  {
    PyErr_SetString(PyExc_TypeError, "gatherDegenerated: listOfAllWins and listOfBlks must be of same size.");
    for (E_Int is = 0; is < nwinstot; is++)
    {
      RELEASESHAREDS(objstt[is], structFt[is]);
    }
    return NULL;
  }
  if (nwins2 != nwins3) 
  {
    PyErr_SetString(PyExc_TypeError, "gatherDegenerated: listOfWinTypes and listOfBlks must be of same size.");
    for (E_Int is = 0; is < nwinstot; is++)
    {
      RELEASESHAREDS(objstt[is], structFt[is]);
    }
    return NULL;
  }
  vector<E_Int> winTypes(nwinstot); vector<E_Int> noBlks(nwinstot);
  for (int now1 = 0; now1 < nwinstot; now1++)
  {
    PyObject* tplw = PyList_GetItem(listOfWinTypes,now1);
    if (PyLong_Check(tplw) == 0 && PyInt_Check(tplw) == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "gatherDegenerated: listOfWinTypes elements must be integers.");
      for (E_Int is = 0; is < nwinstot; is++)
      {
        RELEASESHAREDS(objstt[is], structFt[is]);
      }
      return NULL;
    }
    winTypes[now1] = PyLong_AsLong(tplw);
    PyObject* tplb = PyList_GetItem(listOfBlks,now1);
    if (PyLong_Check(tplb) == 0 && PyInt_Check(tplb) == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "gatherDegenerated: listOfBlks elements must be integers.");
      for (E_Int is = 0; is < nwinstot; is++)
      {
        RELEASESHAREDS(objstt[is], structFt[is]);
      }
      return NULL;
    }
    noBlks[now1] = PyLong_AsLong(tplb);
  }

  /* Extract size Ni,Nj,Nk of blks */
  E_Int sizeNI = PyList_Size(listOfNI);
  E_Int sizeNJ = PyList_Size(listOfNJ);
  E_Int sizeNK = PyList_Size(listOfNK);
  E_Int nblks = sizeNI;
  if (sizeNI != sizeNJ || sizeNI != sizeNK || sizeNJ != sizeNK)
  {
    PyErr_SetString(PyExc_TypeError,
                    "gatherDegenerated: listOfNI, listOfNJ, listOfNK must be of same size.");
    for (E_Int is = 0; is < nwinstot; is++)
    {
      RELEASESHAREDS(objstt[is], structFt[is]);
    }
    return NULL;
  }
  vector<E_Int> NIB(nblks);  vector<E_Int> NJB(nblks); vector<E_Int> NKB(nblks);
  for (int nob1 = 0; nob1 < nblks; nob1++)
  {
    PyObject* tpli = PyList_GetItem(listOfNI,nob1);
    PyObject* tplj = PyList_GetItem(listOfNJ,nob1);
    PyObject* tplk = PyList_GetItem(listOfNK,nob1);
    if ( (PyLong_Check(tpli) == 0 && PyInt_Check(tpli) == 0) || 
         (PyLong_Check(tplj) == 0 && PyInt_Check(tplj) == 0) ||
         (PyLong_Check(tplk) == 0 && PyInt_Check(tplk) == 0) )
    {
      PyErr_SetString(PyExc_TypeError,
                      "gatherDegenerated: listOfNI, listOfNJ and listOfNK elements must all be integers.");
      for (E_Int is = 0; is < nwinstot; is++)
      {
        RELEASESHAREDS(objstt[is], structFt[is]);
      }
      return NULL;
    }    
    NIB[nob1] = PyLong_AsLong(tpli);   
    NJB[nob1] = PyLong_AsLong(tplj);   
    NKB[nob1] = PyLong_AsLong(tplk);   

  }
  /*-------------------- End of checks --------------------------------------*/
  vector<E_Int> imin1; vector<E_Int> jmin1; vector<E_Int> kmin1; 
  vector<E_Int> imax1; vector<E_Int> jmax1; vector<E_Int> kmax1;  
  vector<E_Int> rcvBlk; 
  E_Int indinit1, inds1, jsw1, isw1, iew1, jew1, iwmax1, indw1;
  E_Int ieloc1, foundwin, deltaw1, istart1, jstart1, iend1, jend1;
  E_Int iminloc, imaxloc, jminloc, jmaxloc, kminloc, kmaxloc;
  E_Int im1, jm1, km1, imcw1, jmcw1;
  for (E_Int now1 = 0; now1 < nwinstot; now1++)
  {
    E_Int ncellsw1 =  structFt[now1]->getSize();
    E_Float* tagp1 = structFt[now1]->begin(post1[now1]);
    imcw1 = nitt[now1]; jmcw1 = njtt[now1]; // size of window in centers
    
    //E_Int inde1=-1;
    E_Int noblk1 = noBlks[now1]; 
    E_Int typewin1 = winTypes[now1];
    im1 = NIB[noblk1]; jm1 = NJB[noblk1]; km1 = NKB[noblk1]; 

    indinit1 = 0;    
    restart:;
    for (E_Int ind1 = indinit1; ind1 < ncellsw1; ind1++)
    {
      if (tagp1[ind1] == -3.)// degenerated 
      {
        inds1 = ind1;
        jsw1 = inds1/imcw1; isw1 = inds1-jsw1*imcw1; iew1 = imcw1; jew1 = jsw1; 
        iwmax1 = imcw1;
        for (E_Int jw1 = jsw1; jw1 < jmcw1; jw1++)
        {
          ieloc1 = isw1; deltaw1 = 0;
          for (E_Int iw1 = isw1; iw1 < iwmax1; iw1++)
          {
            indw1 = iw1 + jw1*imcw1;
            if (tagp1[indw1] != -3.) break;
            ieloc1 = iw1; deltaw1++;
          }
          if (deltaw1 == 0) break; // a j constant, pas de pt valide trouve en i
          else 
          {
            foundwin = 1;
            iew1 = K_FUNC::E_min(imcw1-1,ieloc1);
            jew1 = jw1;
            iwmax1 = iew1+1;
          }
        }
        if (foundwin == 1) 
        {
          //inde1 = iew1+jew1*imcw1;
          istart1 = K_FUNC::E_min(isw1+1,iew1+1); iend1 = K_FUNC::E_max(isw1+1,iew1+1);
          jstart1 = K_FUNC::E_min(jsw1+1,jew1+1); jend1 = K_FUNC::E_max(jsw1+1,jew1+1);
          getIndicesInBlk(istart1, iend1, jstart1, jend1, imcw1, jmcw1, typewin1, im1, jm1, km1, 
                          iminloc, imaxloc, jminloc, jmaxloc, kminloc, kmaxloc);
          imaxloc = K_FUNC::E_min(imaxloc,im1);
          jmaxloc = K_FUNC::E_min(jmaxloc,jm1);
          kmaxloc = K_FUNC::E_min(kmaxloc,km1);
          imin1.push_back(iminloc); jmin1.push_back(jminloc); kmin1.push_back(kminloc); 
          imax1.push_back(imaxloc); jmax1.push_back(jmaxloc); kmax1.push_back(kmaxloc);
          rcvBlk.push_back(noblk1);
          // update fields oppositeWins and oppositePts
          for (E_Int j1 = jstart1-1; j1 < jend1; j1++)
            for (E_Int i1 = istart1-1; i1 < iend1; i1++)
            {
              E_Int ind1 = i1 + j1*imcw1;
              tagp1[ind1] = -2.;
            }
        }
        indinit1 = iend1+(jstart1-1)*imcw1;
        goto restart;
      }// 1st pt of tag -3      
    }//parcours des pts
  }
  for (E_Int is = 0; is < nwinstot; is++)
  {
    RELEASESHAREDS(objstt[is], structFt[is]);
  }
  
  PyObject* listwins = PyList_New(0);  
  E_Int sizeL = imin1.size();
  for (E_Int v = 0; v < sizeL; v++)
  {
#define K_LONG long int
    K_LONG noRcvBlk = rcvBlk[v];
    K_LONG im1 = imin1[v]; 
    K_LONG jm1 = jmin1[v];
    K_LONG km1 = kmin1[v];
    K_LONG ip1 = imax1[v]; 
    K_LONG jp1 = jmax1[v]; 
    K_LONG kp1 = kmax1[v]; 
    PyObject* wins = Py_BuildValue("[l,[llllll]]",
                                   noRcvBlk,im1, ip1, jm1, jp1, km1, kp1);
    PyList_Append(listwins, wins);
    Py_DECREF(wins);
  }
  return listwins;
}
