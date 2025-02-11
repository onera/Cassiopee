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

// Identify matching points in windows

# include "connector.h"

using namespace K_FUNC;
using namespace std;
using namespace K_FLD;

//============================================================================
// Gather matching points in windows in structured patches
// IN: listOfAllWins: exterior faces located at centers defined as zones
//     listOfWinTypes: number of the corresponding wins (1 <=> i=1; 2 <=> i=ni 
//     listOfWinBlks: number of the block defining current window
//     listOfNI: size in direction i of all the blks
//     listOfNJ:                   j
//     listOfNK:                   k
// OUT: returns the list of indices [imin,imax,jmin,jmax,kmin,kmax], 
//                          corresponding opposite indices
//                          trirac
//=============================================================================
PyObject* K_CONNECTOR::gatherMatching(PyObject* self, PyObject* args)
{
  PyObject *listOfAllWins, *listOfAllTags, *listOfBlks, *listOfWinTypes, 
    *listOfNI, *listOfNJ, *listOfNK;
  E_Int dimPb;
  E_Float tol;
  if (!PYPARSETUPLE_(args, OOOO_ OOO_ I_ R_,
                    &listOfAllWins, &listOfAllTags, 
                    &listOfWinTypes, &listOfBlks, &listOfNI, &listOfNJ, 
                    &listOfNK, &dimPb, &tol))
  {
      return NULL;
  }

  if (PyList_Check(listOfAllWins) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "gatherMatching: 1st argument must be a list.");
    return NULL;
  }
  if (PyList_Check(listOfAllTags) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "gatherMatching: 2nd argument must be a list.");
    return NULL;
  }
  if (PyList_Check(listOfBlks) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "gatherMatching: 3rd argument must be a list.");
    return NULL;
  }
  if (PyList_Check(listOfWinTypes) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "gatherMatching: 4th argument must be a list.");
    return NULL;
  }
  if (PyList_Check(listOfNI) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "gatherMatching: 5th argument must be a list.");
    return NULL;
  }
  if (PyList_Check(listOfNJ) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "gatherMatching: 6th argument must be a list.");
    return NULL;
  }
  if (PyList_Check(listOfNK) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "gatherMatching: 7th argument must be a list.");
    return NULL;
  }
  /* Extract info on windows located at nodes */
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltTypet;
  vector<PyObject*> objst, objut;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = true;
  E_Boolean skipDiffVars = false;//true;

  E_Int isOk = K_ARRAY::getFromArrays(
    listOfAllWins, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltTypet, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nwinstot = structF.size();
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError, "gatherMatching: 1st list of arrays is not valid.");
    for (E_Int is = 0; is < nwinstot; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    return NULL;
  } 
  if (nwinstot == 0) 
  {
    PyErr_SetString(PyExc_TypeError, "gatherMatching: 1st list of arrays does not contain valid structured zones.");
    for (E_Int is = 0; is < nwinstot; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    return NULL;
  }
  vector<E_Int> posxt; vector<E_Int> posyt; vector<E_Int> poszt; 
  E_Int posxi, posyi, poszi;
  for (E_Int i = 0; i < nwinstot; i++)
  {
    posxi = K_ARRAY::isCoordinateXPresent(structVarString[i]); posxi++; 
    posyi = K_ARRAY::isCoordinateYPresent(structVarString[i]); posyi++;
    poszi = K_ARRAY::isCoordinateZPresent(structVarString[i]); poszi++;    
    posxt.push_back(posxi); posyt.push_back(posyi); poszt.push_back(poszi);
  }

  /* Extract info on tag1 and tag2 */
  vector<E_Int> rest;
  vector<char*> structVarStringt; vector<char*> unstrVarStringt;
  vector<FldArrayF*> structFt; vector<FldArrayF*> unstrFt;
  vector<E_Int> nitt; vector<E_Int> njtt; vector<E_Int> nktt;
  vector<FldArrayI*> cntt;
  vector<char*> eltTypett;
  vector<PyObject*> objstt, objutt;
  skipNoCoord = false;
  skipStructured = false;
  skipUnstructured = true;
  skipDiffVars = true;

  isOk = K_ARRAY::getFromArrays(
    listOfAllTags, rest, structVarStringt, unstrVarStringt,
    structFt, unstrFt, nitt, njtt, nktt, cntt, eltTypett, objstt, objutt, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int ntagstot = structFt.size();
  if (ntagstot != nwinstot) 
  {
    PyErr_SetString(PyExc_TypeError, "gatherMatching: 1st and 2nd lists of arrays must be of same size.");
    for (E_Int is = 0; is < nwinstot; is++)
      RELEASESHAREDS(objst[is], structF[is]);
    for (E_Int is = 0; is < ntagstot; is++)
      RELEASESHAREDS(objstt[is], structFt[is]);
    return NULL;
  }
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError, "gatherMatching: 2nd list of arrays is not valid.");
    for (E_Int is = 0; is < nwinstot; is++)
    {
      RELEASESHAREDS(objst[is], structF[is]);
      RELEASESHAREDS(objstt[is], structFt[is]);
    }
    return NULL;
  } 
  if (nwinstot == 0) 
  {
    PyErr_SetString(PyExc_TypeError, "gatherMatching: 2nd list of arrays does not contain valid structured zones.");
    for (E_Int is = 0; is < nwinstot; is++)
    {
      RELEASESHAREDS(objst[is], structF[is]);
      RELEASESHAREDS(objstt[is], structFt[is]);
    }
    return NULL;
  }
  vector<E_Int> post1; vector<E_Int> post2;
  E_Int postag1, postag2;
  for (E_Int i = 0; i < nwinstot; i++)
  {   
    postag1 = K_ARRAY::isNamePresent("tag1", structVarStringt[i]); postag1++;
    postag2 = K_ARRAY::isNamePresent("tag2", structVarStringt[i]); postag2++;
    if (postag1 < 1 || postag2 < 1) 
    {
      PyErr_SetString(PyExc_TypeError, "gatherMatching: tag1 and tag2 must be defined in listOfAllWins.");
      for (E_Int is = 0; is < nwinstot; is++)
      {
        RELEASESHAREDS(objstt[is], structFt[is]);
        RELEASESHAREDS(objst[is], structF[is]);
      }
      return NULL;
    } 
    post1.push_back(postag1); post2.push_back(postag2); 
  }
  
  /* Extract window type and corresponding original block */
  E_Int nwins2 = PyList_Size(listOfWinTypes);
  if (nwins2 != nwinstot)
  {
    PyErr_SetString(PyExc_TypeError, "gatherMatching: listOfAllWins and listOfWinTypes must be of same size.");
    for (E_Int is = 0; is < nwinstot; is++)
    {
      RELEASESHAREDS(objst[is], structF[is]);
      RELEASESHAREDS(objstt[is], structFt[is]);
    }
    return NULL;
  }
  E_Int nwins3 = PyList_Size(listOfBlks);
  if (nwins3 != nwinstot)
  {
    PyErr_SetString(PyExc_TypeError, "gatherMatching: listOfAllWins and listOfBlks must be of same size.");
    for (E_Int is = 0; is < nwinstot; is++)
    {
      RELEASESHAREDS(objstt[is], structFt[is]);
      RELEASESHAREDS(objst[is], structF[is]);
    }
    return NULL;
  }
  if (nwins2 != nwins3) 
  {
    PyErr_SetString(PyExc_TypeError, "gatherMatching: listOfWinTypes and listOfBlks must be of same size.");
    for (E_Int is = 0; is < nwinstot; is++)
    {
      RELEASESHAREDS(objst[is], structF[is]);
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
                      "gatherMatching: listOfWinTypes elements must be integers.");
      for (E_Int is = 0; is < nwinstot; is++)
      {
        RELEASESHAREDS(objst[is], structF[is]);
        RELEASESHAREDS(objstt[is], structFt[is]);
      }
      return NULL;
    }
    winTypes[now1] = PyLong_AsLong(tplw);
    PyObject* tplb = PyList_GetItem(listOfBlks,now1);
    if (PyLong_Check(tplb) == 0 && PyInt_Check(tplb) == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "gatherMatching: listOfBlks elements must be integers.");
      for (E_Int is = 0; is < nwinstot; is++)
      {
        RELEASESHAREDS(objst[is], structF[is]);
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
                    "gatherMatching: listOfNI, listOfNJ, listOfNK must be of same size.");
    for (E_Int is = 0; is < nwinstot; is++)
    {
      RELEASESHAREDS(objst[is], structF[is]);
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
    if ((PyLong_Check(tpli) == 0 && PyInt_Check(tpli) == 0) || 
        (PyLong_Check(tplj) == 0 && PyInt_Check(tplj) == 0) ||
        (PyLong_Check(tplk) == 0 && PyInt_Check(tplk) == 0))
    {
      PyErr_SetString(PyExc_TypeError,
                      "gatherMatching: listOfNI, listOfNJ and listOfNK elements must all be integers.");
      for (E_Int is = 0; is < nwinstot; is++)
      {
        RELEASESHAREDS(objst[is], structF[is]);
        RELEASESHAREDS(objstt[is], structFt[is]);
      }
      return NULL;
    }    
    NIB[nob1] = PyLong_AsLong(tpli);   
    NJB[nob1] = PyLong_AsLong(tplj);   
    NKB[nob1] = PyLong_AsLong(tplk);   

  }
  /*-------------------- End of checks --------------------------------------*/
  E_Int nmatch = 4;
  if (dimPb == 2) nmatch = 2;
  vector<E_Int> imin1; vector<E_Int> imin2;
  vector<E_Int> jmin1; vector<E_Int> jmin2;
  vector<E_Int> kmin1; vector<E_Int> kmin2;
  vector<E_Int> imax1; vector<E_Int> imax2;
  vector<E_Int> jmax1; vector<E_Int> jmax2;
  vector<E_Int> kmax1; vector<E_Int> kmax2;  
  vector<E_Int> rac1; vector<E_Int> rac2; vector<E_Int> rac3;
  vector<E_Int> rcvBlk; vector<E_Int> dnrBlk;
  vector<E_Int> rcvWinNbr; vector<E_Int> dnrWinNbr;

  E_Int im1, jm1, km1, im2, jm2, km2;
  E_Int i1, j1, i2, j2;
  E_Int imw1, jmw1, imw2, imcw1, jmcw1, imcw2, jmcw2;
  E_Int iminloc, imaxloc, jminloc, jmaxloc, kminloc, kmaxloc;
  E_Int inci, incj, inciopp, incjopp;
  E_Int incci, inccj, incciopp, inccjopp;
  E_Int isw1, jsw1, iew1, jew1, isw2, jsw2, iew2, jew2, indw1;
  E_Int ieloc1, deltaw1, deltaw2, foundwin, indinit1, iwmax1;
  E_Int istart1, iend1, istart2, iend2, jstart1, jend1, jstart2, jend2;
  E_Int isMatch;
  vector<E_Int> indTab1(nmatch);
  vector<E_Int> indTab2(nmatch);
  vector<E_Int> noOpp(nmatch);// opposite matching node in indtab2 of node indtab1[i] 
  for (E_Int now1 = 0; now1 < nwinstot; now1++)
  {
    //E_Int nptsw1 =  structF[now1]->getSize();
    E_Int ncellsw1 =  structFt[now1]->getSize();
    E_Float* oppositeWins = structFt[now1]->begin(post1[now1]);
    E_Float* oppositePts = structFt[now1]->begin(post2[now1]);
    E_Float* xn1 = structF[now1]->begin(posxt[now1]);
    E_Float* yn1 = structF[now1]->begin(posyt[now1]);
    E_Float* zn1 = structF[now1]->begin(poszt[now1]);

    imw1 = nit[now1]; jmw1 = njt[now1]; // size of window in nodes
    imcw1 = nitt[now1]; jmcw1 = njtt[now1]; // size of window in centers
    //1. Search for the 1st pt of tag1 > -1
    E_Int inds1 = -1; E_Int noblk2 = -1; E_Int inds2 = -1;
    E_Int typewin2 = -1; E_Int now2 = -1;
    E_Int inde2=-1;
    E_Int noblk1 = noBlks[now1]; 
    E_Int typewin1 = winTypes[now1];
    im1 = NIB[noblk1]; jm1 = NJB[noblk1]; km1 = NKB[noblk1]; 
    indinit1 = 0;
    
    restart:;
    inds1 = -1; 
    for (E_Int ind1 = indinit1; ind1 < ncellsw1; ind1++)
    {
      if (oppositeWins[ind1] >-1.) 
      { 
        // Get the indices (i,j) of first node 
        j1 = ind1/imcw1; i1 = ind1-j1*imcw1;
        if (i1 == imw1-1 || (j1 == jmw1-1 && dimPb == 3)) goto next;  

        indTab1[0] =  i1+j1*imw1;
        indTab1[1] =  i1+1+j1*imw1;
        if (nmatch == 4)
        {
          indTab1[2] =  i1+(j1+1)*imw1;
          indTab1[3] =  i1+1+(j1+1)*imw1;
        }
        isMatch = 0;
        // Check the four vertices of current cell
        // Do they match the four vertices of the opposite cell ?
        now2 = E_Int(oppositeWins[ind1]); inds2 = E_Int(oppositePts[ind1]);
        E_Float* xn2 = structF[now2]->begin(posxt[now2]);
        E_Float* yn2 = structF[now2]->begin(posyt[now2]);
        E_Float* zn2 = structF[now2]->begin(poszt[now2]);
        noblk2 = noBlks[now2]; typewin2 = winTypes[now2]; 
        im2 = NIB[noblk2]; jm2 = NJB[noblk2]; km2 = NKB[noblk2]; 
        imw2 = nit[now2]; // size of window in nodes
        imcw2 = nitt[now2]; jmcw2 = njtt[now2]; // size of window in centers
        j2 = inds2/imcw2; i2 = inds2-j2*imcw2;
        indTab2[0] = i2+j2*imw2; indTab2[1] = i2+1+j2*imw2;

        if (nmatch == 4)
        { 
          indTab2[2] = i2  +(j2+1)*imw2;
          indTab2[3] = i2+1+(j2+1)*imw2;
        }
        
        for (E_Int noi1 = 0; noi1 < nmatch; noi1++)
          for (E_Int noi2 = 0; noi2 < nmatch; noi2++)
          {
            if (K_FUNC::E_abs(xn1[indTab1[noi1]]-xn2[indTab2[noi2]]) < tol && 
                K_FUNC::E_abs(yn1[indTab1[noi1]]-yn2[indTab2[noi2]]) < tol && 
                K_FUNC::E_abs(zn1[indTab1[noi1]]-zn2[indTab2[noi2]]) < tol ) 
            {
              noOpp[noi1] = noi2; isMatch++;             
            }
          }
        if (isMatch >= nmatch)
        {
          inds1 = ind1;
          break;
        }
      }
      next:;
    }

    if (inds1 == -1 || isMatch < nmatch) goto nextwin;

    // We know inds1/inds2, now1/now2, noblk1/noblk2, typewin1/typewin2 and indTab1/indTab2
    // Determine the increments
    inci = 1; incj = imw1;
    inciopp = indTab2[noOpp[1]]-indTab2[noOpp[0]]; incjopp = 0;
    if (nmatch == 4)
    {incjopp = indTab2[noOpp[2]]-indTab2[noOpp[0]];}
    
    // Degenerated cases: we assume trirac 1,2,3
    if (isMatch > nmatch) 
    {
      inciopp = 1; incjopp = imw2;
    }
    // Build subwindows
    compIncrement(inds1, imcw1, oppositeWins, oppositePts, dimPb,
                  incci, inccj, incciopp, inccjopp);
    
    jsw1 = inds1/imcw1; isw1 = inds1-jsw1*imcw1; iew1 = imcw1; jew1 = jsw1; 
    deltaw2 = 0; iwmax1 = imcw1;
    for (E_Int jw1 = jsw1; jw1 < jmcw1; jw1++)
    {
      ieloc1 = isw1; deltaw1 = 0;
      for (E_Int iw1 = isw1; iw1 < iwmax1; iw1++)
      {
        indw1 = iw1 + jw1 * imcw1;
        if (oppositeWins[indw1] < 0.)
        {
          break;
        }
        else 
        {
          if (E_Int(oppositeWins[indw1]) != now2) 
          {
            break;
          }
          if (E_Int(oppositePts[indw1]) != inds2+incciopp*deltaw1+inccjopp*deltaw2)
          {
            break;
          }
          deltaw1++; ieloc1 = iw1;
        }
      }
      if (deltaw1 == 0) break; // a j constant, pas de pt valide trouve en i
      else 
      {
        foundwin = 1;

        iew1 = K_FUNC::E_min(imcw1-1,ieloc1);
        jew1 = jw1; deltaw2++;
        iwmax1 = iew1+1;
      }
    }

    // Matching patch found
    if (foundwin == 1) 
    {
      //inde1 = iew1+jew1*imcw1;
      // rcvWin [isw1,iew1,jsw1,jew1] found
      istart1 = K_FUNC::E_min(isw1+1,iew1+1); iend1 = K_FUNC::E_max(isw1+1,iew1+1);
      jstart1 = K_FUNC::E_min(jsw1+1,jew1+1); jend1 = K_FUNC::E_max(jsw1+1,jew1+1);
      getIndicesInBlk(istart1, iend1, jstart1, jend1, imcw1, jmcw1, typewin1, im1, jm1, km1, 
                      iminloc, imaxloc, jminloc, jmaxloc, kminloc, kmaxloc);
      imaxloc = K_FUNC::E_min(imaxloc,im1);
      jmaxloc = K_FUNC::E_min(jmaxloc,jm1);
      kmaxloc = K_FUNC::E_min(kmaxloc,km1);
      imin1.push_back(iminloc); jmin1.push_back(jminloc); kmin1.push_back(kminloc); 
      imax1.push_back(imaxloc); jmax1.push_back(jmaxloc); kmax1.push_back(kmaxloc); 
      // calcul trirac            
      if (dimPb == 3) compTrirac(im1, jm1, im2, jm2, typewin1, inci, incj, typewin2, inciopp, incjopp, rac1, rac2, rac3);
      else compTrirac2D(im1, jm1, im2, jm2, typewin1, typewin2, inci, inciopp, rac1, rac2, rac3);
      rcvBlk.push_back(noblk1); dnrBlk.push_back(noblk2);
      // dnrWin to be computed
      inde2 = E_Int(oppositePts[iew1+jew1*imcw1]); // max pt
      jsw2 = inds2/imcw2; isw2 = inds2-jsw2*imcw2;      
      jew2 = inde2/imcw2; iew2 = inde2-jew2*imcw2;
      istart2 = K_FUNC::E_min(isw2,iew2)+1; iend2 = K_FUNC::E_max(isw2,iew2)+1;
      jstart2 = K_FUNC::E_min(jsw2,jew2)+1; jend2 = K_FUNC::E_max(jsw2,jew2)+1;

      getIndicesInBlk(istart2, iend2, jstart2, jend2, imcw2, jmcw2, typewin2, im2, jm2, km2, 
                      iminloc, imaxloc, jminloc, jmaxloc, kminloc, kmaxloc);
      imaxloc = K_FUNC::E_min(imaxloc,im2);
      jmaxloc = K_FUNC::E_min(jmaxloc,jm2);
      kmaxloc = K_FUNC::E_min(kmaxloc,km2);
      imin2.push_back(iminloc); jmin2.push_back(jminloc); kmin2.push_back(kminloc); 
      imax2.push_back(imaxloc); jmax2.push_back(jmaxloc); kmax2.push_back(kmaxloc); 

      rcvWinNbr.push_back(now1); dnrWinNbr.push_back(now2);
      // update fields oppositeWins and oppositePts
      for (E_Int j1 = jstart1-1; j1 < jend1; j1++)
        for (E_Int i1 = istart1-1; i1 < iend1; i1++)
        {
          E_Int ind1 = i1 + j1*imcw1;
          oppositePts[ind1] = -2.;
          oppositeWins[ind1] = -2.;
        }
      E_Float* oppositeWins2 = structFt[now2]->begin(post1[now2]);
      E_Float* oppositePts2 = structFt[now2]->begin(post2[now2]); 
      for (E_Int j2 = jstart2-1; j2 < jend2; j2++)
        for (E_Int i2 = istart2-1; i2 < iend2; i2++)
        {
          E_Int ind2 = i2 + j2*imcw2;
          oppositePts2[ind2] = -2.;
          oppositeWins2[ind2] = -2.;
        }
    }  
    indinit1 = iend1+(jstart1-1)*imcw1;
    goto restart;
    nextwin:;
  }
  for (E_Int is = 0; is < nwinstot; is++)
  {
    RELEASESHAREDS(objst[is], structF[is]);
    RELEASESHAREDS(objstt[is], structFt[is]);
  }

  PyObject* listwins = PyList_New(0);  
  E_Int sizeL = imin1.size();
  for (E_Int v = 0; v < sizeL; v++)
  {
    // E_LONG ne marche pas
#define K_LONG long int
    K_LONG noRcvBlk = rcvBlk[v];
    K_LONG noDnrBlk = dnrBlk[v];
    K_LONG noRcvWin = rcvWinNbr[v];
    K_LONG noDnrWin = dnrWinNbr[v];
    K_LONG im1 = imin1[v]; K_LONG im2 = imin2[v];
    K_LONG jm1 = jmin1[v]; K_LONG jm2 = jmin2[v];
    K_LONG km1 = kmin1[v]; K_LONG km2 = kmin2[v];
    K_LONG ip1 = imax1[v]; K_LONG ip2 = imax2[v];
    K_LONG jp1 = jmax1[v]; K_LONG jp2 = jmax2[v];
    K_LONG kp1 = kmax1[v]; K_LONG kp2 = kmax2[v];
    K_LONG trirac11 = rac1[v]; 
    K_LONG trirac12 = rac2[v]; 
    K_LONG trirac13 = rac3[v]; 
    PyObject* wins = Py_BuildValue("[[ll],[llllll],[llllll],[lll],[ll]]",
                                   noRcvBlk, noDnrBlk,
                                   im1, ip1, jm1, jp1, km1, kp1,
                                   im2, ip2, jm2, jp2, km2, kp2,
                                   trirac11, trirac12, trirac13,
                                   noRcvWin, noDnrWin);
    PyList_Append(listwins, wins);
    Py_DECREF(wins);
  }
  return listwins;
}
//==============================================================================
void K_CONNECTOR::getIndicesInBlk(
  E_Int isw1, E_Int iew1, E_Int jsw1, E_Int jew1, 
  E_Int imw1, E_Int jmw1, E_Int typewin1, 
  E_Int im1, E_Int jm1, E_Int km1,
  E_Int& imin, E_Int& imax, E_Int& jmin, E_Int& jmax, E_Int& kmin, E_Int& kmax)
{
  switch (typewin1)
  {
    case 1:
      imin = 1; imax = 1;

      if (isw1 == 0) jmin = 1;
      else jmin = isw1;
      if (jsw1 == 0) kmin = 1;
      else kmin = jsw1;

      if (isw1 == jm1-1) jmax = jm1;
      else jmax = iew1+1; 
      if (jsw1 == km1-1) kmax = km1;
      else kmax = jew1+1; 

      break;

    case 2: 
      imin = im1; imax = im1;
      if (isw1 == 0) jmin = 1;
      else jmin = isw1;
      if (jsw1 == 0) kmin = 1;
      else kmin = jsw1;

      if (isw1 == jm1-1) jmax = jm1;
      else jmax = iew1+1; 
      if (jsw1 == km1-1) kmax = km1;
      else kmax = jew1+1; 
      break;

    case 3:
      jmin = 1; jmax = 1;

      if (isw1 == 0) imin = 1;
      else imin = isw1;
      if (jsw1 == 0) kmin = 1;
      else kmin = jsw1;

      if (isw1 == im1-1) imax = im1;
      else imax = iew1+1; 
      if (jsw1 == km1-1) kmax = km1;
      else kmax = jew1+1; 
      break;

    case 4:
      jmin = jm1; jmax = jm1;
      if (isw1 == 0) imin = 1;
      else imin = isw1;
      if (jsw1 == 0) kmin = 1;
      else kmin = jsw1;

      if (isw1 == im1-1) imax = im1;
      else imax = iew1+1; 
      if (jsw1 == km1-1) kmax = km1;
      else kmax = jew1+1;
      break;

    case 5:
      if (isw1 == 0) imin = 1;
      else imin = isw1;
      if (jsw1 == 0) jmin = 1;
      else jmin = jsw1;

      if (isw1 == im1-1) imax = im1;
      else imax = iew1+1; 
      if (jsw1 == jm1-1) jmax = jm1;
      else jmax = jew1+1;
      kmin = 1; kmax = 1;
      break;

    case 6:
      if (isw1 == 0) imin = 1;
      else imin = isw1;
      if (jsw1 == 0) jmin = 1;
      else jmin = jsw1;

      if (isw1 == im1-1) imax = im1;
      else imax = iew1+1; 
      if (jsw1 == jm1-1) jmax = jm1;
      else jmax = jew1+1;
      kmin = km1; kmax = km1;
      break;
  }
}

//=============================================================================
/* Calcul des trirac */
//=============================================================================
void K_CONNECTOR::compTrirac(E_Int im1, E_Int jm1, E_Int im2, E_Int jm2, 
                             E_Int typewin1, E_Int inc1, E_Int inc2, 
                             E_Int typewin2, E_Int incm1, E_Int incm2,
                             vector<E_Int>& rac1, vector<E_Int>& rac2, vector<E_Int>& rac3)
{
  switch (typewin1)
  {
    case 1:
    case 2:
      if (inc1 > 0) inc1 = im1; 
      else inc1 = -im1;
      if (inc2 > 0) inc2 = im1*jm1; 
      else inc2 = -im1*jm1;
      break;
    case 3:
    case 4:
      if (inc1 > 0) inc1 = 1; 
      else inc1 = -1;
      if (inc2 > 0) inc2 = im1*jm1; 
      else inc2 = -im1*jm1;
      break;

    case 5:
    case 6:
      if (inc1 > 0) inc1 = 1; 
      else inc1 = -1;
      if (inc2 > 0) inc2 = im1; 
      else inc2 = -im1;
      break;
  }

  switch (typewin2)
  {
    case 1:
    case 2: // (i,j) devient (j,k)

      if (incm1 == 1 ) incm1=im2; 
      else if (incm1 == -1) incm1=-im2;
      else if (incm1 > 1) incm1=im2*jm2;
      else incm1 = -im2*jm2;
      
      if (incm2 == 1) incm2=im2; 
      else if (incm2 == -1) incm2=-im2;
      else if (incm2 > 1) incm2=im2*jm2;
      else incm2 = -im2*jm2;
      
      break;

    case 3:
    case 4:// (i,j) devient (i,k)

      if (incm1 == 1 ) incm1=1; 
      else if (incm1 == -1) incm1=-1;
      else if (incm1 > 1) incm1=im2*jm2;
      else incm1=-im2*jm2;
      
      if (incm2 == 1) incm2=1; 
      else if (incm2 == -1) incm2=-1;
      else if (incm2 > 1) incm2=im2*jm2;
      else incm2=-im2*jm2;

      break;

    case 5:
    case 6:// (i,j) devient (i,j)

      if (incm1 == 1) incm1=1; 
      else if (incm1 == -1) incm1=-1;
      else if (incm1 > 1) incm1=im2;
      else incm1=-im2;
      
      if (incm2 == 1) incm2=1; 
      else if (incm2 == -1) incm2=-1;
      else if (incm2 > 1) incm2=im2;
      else incm2=-im2;

      break;
  
  }
 
  E_Int r1 = 0;
  E_Int r2 = 0;
  E_Int r3 = 0;

  if (inc1 == 1)
  {
    if (incm1 == 1)    r1 =  1;
    else if (incm1 == -1)   r1 = -1;
    else if (incm1 == im2)  r1 =  2;
    else if (incm1 == -im2) r1 = -2;
    else if (incm1 == im2*jm2)  r1 =  3;
    else if (incm1 == -im2*jm2) r1 = -3;
  }
  else if (inc1 == im1)
  {
    if (incm1 == 1)    r2 =  1;
    else if (incm1 == -1)   r2 = -1;
    else if (incm1 == im2)  r2 =  2;
    else if (incm1 == -im2) r2 = -2;
    else if (incm1 == im2*jm2)  r2 =  3;
    else if (incm1 == -im2*jm2) r2 = -3;
  }
  else if (inc1 == im1*jm1)
  {
    if (incm1 == 1)    r3 =  1;
    else if (incm1 == -1)   r3 = -1;
    else if (incm1 == im2)  r3 =  2;
    else if (incm1 == -im2) r3 = -2;
    else if (incm1 == im2*jm2)  r3 =  3;
    else if (incm1 == -im2*jm2) r3 = -3;
  }

  if (inc2 == 1)
  {
    if (incm2 == 1)    r1 =  1;
    else if (incm2 == -1)   r1 = -1;
    else if (incm2 == im2)  r1 =  2;
    else if (incm2 == -im2) r1 = -2;
    else if (incm2 == im2*jm2)  r1 =  3;
    else if (incm2 == -im2*jm2) r1 = -3;
  }
  else if (inc2 == im1)
  {
    if (incm2 == 1)    r2 =  1;
    else if (incm2 == -1)   r2 = -1;
    else if (incm2 == im2)  r2 =  2;
    else if (incm2 == -im2) r2 = -2;
    else if (incm2 == im2*jm2)  r2 =  3;
    else if (incm2 == -im2*jm2) r2 = -3;
  }
  else if (inc2 == im1*jm1)
  {
    if (incm2 == 1)    r3 =  1;
    else if (incm2 == -1)   r3 = -1;
    else if (incm2 == im2)  r3 =  2;
    else if (incm2 == -im2) r3 = -2;
    else if (incm2 == im2*jm2)  r3 =  3;
    else if (incm2 == -im2*jm2) r3 = -3;
  }
  E_Int sum = K_FUNC::E_abs(r1) + K_FUNC::E_abs(r2) + K_FUNC::E_abs(r3);
  if (r1 == 0)
  {
    if (sum == 3) r1 = signature(3, r2, r3, 1);
    if (sum == 4) r1 = signature(2, r2, r3, 1);
    if (sum == 5) r1 = signature(1, r2, r3, 1);
  }
  else if (r2 == 0)
  {
    if (sum == 3) r2 = signature(r1, 3, r3, 2);
    if (sum == 4) r2 = signature(r1, 2, r3, 2);
    if (sum == 5) r2 = signature(r1, 1, r3, 2);
  }
  else if (r3 == 0)
  {
    if (sum == 3) r3 = signature(r1, r2, 3, 3);
    if (sum == 4) r3 = signature(r1, r2, 2, 3);
    if (sum == 5) r3 = signature(r1, r2, 1, 3);
  }
  rac1.push_back(r1); rac2.push_back(r2); rac3.push_back(r3);
}

//=============================================================================
// Soit un trirac donnee (r1,r2,r3), modifie le trirac no var
// pour que le triedre soit direct. Retourne le trirac modifie
//=============================================================================
E_Int K_CONNECTOR::signature(E_Int r1, E_Int r2, E_Int r3, E_Int var)
{
  E_Float s;
  E_Int n = 0;
  if (r1 < 0) n++;
  if (r2 < 0) n++;
  if (r3 < 0) n++;
  E_Float ar1 = K_FUNC::E_abs(r1);
  E_Float ar2 = K_FUNC::E_abs(r2);
  E_Float ar3 = K_FUNC::E_abs(r3);
  s = (ar2 - ar1)*(ar3 - ar1)*(ar3 - ar2)*0.5;
  //s = s * pow(-1.,n);
  if (2*(n/2) - n != 0) s = -s;

  if (var == 1)
  {
    if (s < 0) return -r1;
    else return r1;
  }
  else if (var == 2)
  {
    if (s < 0) return -r2;
    else return r2;
  }
  else if (var == 3)
  {
    if (s < 0) return -r3;
    else return r3;
  }
  return 0;
}
//=============================================================================
/* IN: indA1: index of first matching point of win now1
   IN: imw1: size of win now1
   IN: oppositeWins, oppositePts: infos for opposite matching pts: win number and index
   IN: dimPb: problem dimension 
   OUT: inci, incj: valid increments for matching points in now1
   OUT: inciopp, incjopp: valid increments for matching opposite points in now2
   If no valid increment in one direction, set to 0 */
//==============================================================================
void K_CONNECTOR::compIncrement(
  E_Int indA1, E_Int imw1, E_Float* oppositeWins, E_Float* oppositePts, E_Int dimPb,
  E_Int& inci, E_Int& incj, E_Int& inciopp, E_Int& incjopp)
{
  E_Int indB1 = indA1+1; 
  E_Int indD1 = indA1+imw1;
  inci = 0; incj = 0;
  inciopp = 0; incjopp = 0;  
  if (oppositeWins[indB1] == oppositeWins[indA1])
  {
    inci = 1;
    inciopp = E_Int(oppositePts[indB1])-E_Int(oppositePts[indA1]);    
  }
  
  if (dimPb == 3 && oppositeWins[indD1] == oppositeWins[indA1])
  {
    incj = imw1;
    incjopp = E_Int(oppositePts[indD1])-E_Int(oppositePts[indA1]);    
  }  
}
  
//=============================================================================
/* Calcul des trirac */
//=============================================================================
void K_CONNECTOR::compTrirac2D(
  E_Int im1, E_Int jm1, E_Int im2, E_Int jm2, 
  E_Int typewin1, E_Int typewin2, E_Int inc1, E_Int incm1,
  vector<E_Int>& rac1, vector<E_Int>& rac2, vector<E_Int>& rac3)
{  
  switch (typewin1)
  {
    case 1:
    case 2:
      if (inc1 > 0) inc1 = im1; 
      else inc1 = -im1;
      break;
    case 3:
    case 4:
      if (inc1 > 0) inc1 = 1; 
      else inc1 = -1;
      break;
  }
  switch (typewin2)
  {
    case 1:
    case 2:
      if (incm1 > 0) incm1 = im2; 
      else incm1 = -im2;
      break;
    case 3:
    case 4:
      if (incm1 > 0) incm1 = 1; 
      else incm1 = -1;
      break;

  }
  E_Int r1 = 0;
  E_Int r2 = 0;

  if (inc1 == 1)
  {
    if (incm1 == 1)    r1 = 1;
    else if (incm1 == -1)   r1 = -1;
    else if (incm1 == im2)  r1 =  2;
    else if (incm1 == -im2) r1 = -2;
  }
  else if (inc1 == im1)
  {
    if (incm1 == 1)    r2 =  1;
    else if (incm1 == -1)   r2 = -1;
    else if (incm1 == im2)  r2 =  2;
    else if (incm1 == -im2) r2 = -2;
  }

  if (r1 == 0)
  {
    if (r2 == 1) r1 =-2;
    else if (r2 ==-1) r1 = 2;
    else if (r2 == 2) r1 = 1;
    else if (r2 ==-2) r1 =-1;
  }
  else if (r2 == 0)
  {
    if (r1 == 1) r2 = 2;
    else if (r1 ==-1) r2 =-2;
    else if (r1 == 2) r2 =-1;
    else if (r1 ==-2) r2 = 1; 
  }

  rac1.push_back(r1); rac2.push_back(r2); rac3.push_back(3);
}
