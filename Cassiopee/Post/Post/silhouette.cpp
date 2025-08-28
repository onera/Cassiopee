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

#include "post.h"

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

//=============================================================================
/* Detect shape edges on a surface defined by a TRI */
//=============================================================================
PyObject* K_POST::silhouette(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float vx, vy, vz;

  if (!PYPARSETUPLE_(args, O_ TRRR_, &array, &vx, &vy, &vz))
  {
    return NULL;
  }

  /*-----------------------------------------------*/
  /* Extract datas from surface mesh               */ 
  /*-----------------------------------------------*/
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* eltType; char* varString;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, nil, njl, nkl, 
                                     cn, eltType);
  // check array contains TRI
  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "silhouette: array must be unstructured.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }

  // element type
  E_Int type = 0;
  if (strcmp(eltType, "BAR") == 0) type = 2;
  else if (strcmp(eltType, "TRI") == 0) type = 3;
  else if (strcmp(eltType, "QUAD") == 0) type = 4;
  if (type == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "silhouette: array must be TRI or QUAD.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }  

  // get coordinates
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "silhouette: array must contain coordinates.");
    RELEASESHAREDU(array, f, cn); return NULL;  
  }  

  // get connectivity
  E_Int nelts = cn->getSize();

  // get Element-Element neighbours connectivity
  vector< vector<E_Int> > cEEN(nelts);
  E_Int npts = f->getSize();
  K_CONNECT::connectEV2EENbrs(eltType, npts, *cn, cEEN); 

  // prodscal: scalar product
  E_Float prodscal;
  // prodsign: sign of scalar product.
  E_Int prodsign1, prodsign2;
  // pointers on field f
  E_Int nfld = f->getNfld();
  vector<E_Float*> fp(nfld);
  for (E_Int p=0;p<nfld;p++) fp[p] = f->begin(p+1);
  // cnnew: new BAR connectivity (storing shape edges)
  FldArrayI* cnnew = new FldArrayI(nelts, 2);
  E_Int* cnnew1 = cnnew->begin(1); E_Int* cnnew2 = cnnew->begin(2);
  // fnew: new field associated with cnnew
  FldArrayF* fnew = new FldArrayF(npts,nfld);
  vector<E_Float*> fnewp(nfld);
  for (E_Int p = 0; p < nfld; p++) fnewp[p] = fnew->begin(p+1);
  // index for new connectivity
  E_Int noe = 0; E_Int nop = 0;
  // found indices of matching edge for new BAR connectivity
  E_Int found1 = 0; E_Int found2 = 1;
  // index of neighbour elements
  E_Int elt2;
  // surface arrays and data for surface computation
  FldArrayF nsurf(nelts,3);
  E_Float* nsurfx = nsurf.begin(1);E_Float* nsurfy = nsurf.begin(2); E_Float* nsurfz = nsurf.begin(3);
  FldArrayF surf(nelts, 1);

  // Compute sign of scalar product nsurf*v of elements and compare with neighbours
  // When sign is different between two neighbours, store matching edge
  if (type == 2) // BAR
  {
    E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2);
    // indA,indB: indices of BAR element
    E_Int indA, indB;
    // return all the elements
    for (E_Int elt1 = 0; elt1 < nelts; elt1++)
    {
      indA = cn1[elt1]-1; indB = cn2[elt1]-1;
      for (E_Int p=0;p<nfld;p++)
      {
        fnewp[p][indA] = fp[p][indA]; fnewp[p][indB] = fp[p][indB];
        cnnew1[elt1] = cn1[elt1]; cnnew2[elt1] = cn2[elt1];
      }
    }
  }
  else if (type == 3) // TRI
  {
    E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2); E_Int* cn3 = cn->begin(3);
    // indAi,indBi,indCi: indices of a TRI element i
    E_Int indA1, indB1, indC1, indA2, indB2, indC2;
    // ptA1,ptB1,ptC1: pointers on f for points A1,B1,C1
    vector<E_Float> ptA1(nfld); vector<E_Float> ptB1(nfld); vector<E_Float> ptC1(nfld);

    // Compute surface vector of triangles
    K_METRIC::compUnstructSurf(
      *cn, "TRI", fp[posx], fp[posy], fp[posz], 
      nsurfx, nsurfy, nsurfz, surf.begin());

    // Loop on all elements
    for (E_Int elt1 = 0; elt1 < nelts; elt1++)
    {
      indA1 = cn1[elt1]-1; indB1 = cn2[elt1]-1; indC1 = cn3[elt1]-1;
      for (E_Int p=0;p<nfld;p++){ptA1[p] = fp[p][indA1];ptB1[p] = fp[p][indB1];ptC1[p] = fp[p][indC1];}
      // sign of scalar product for element elt1
      prodsign1=0;
      prodscal = nsurfx[elt1]*vx + nsurfy[elt1]*vy + nsurfz[elt1]*vz;
      if (prodscal < 0.) prodsign1 = 1;
      vector<E_Int>& eltsVoisins = cEEN[elt1]; E_Int nvoisins = eltsVoisins.size();
      // Loop on all neighbours
      for (E_Int et2 = 0; et2 < nvoisins; et2++)
      {
        elt2 = eltsVoisins[et2];
        indA2 = cn1[elt2]-1; indB2 = cn2[elt2]-1; indC2 = cn3[elt2]-1;
        // sign of scalar product for neighbours
        prodsign2 = 0;
        prodscal = nsurfx[elt2]*vx + nsurfy[elt2]*vy + nsurfz[elt2]*vz;
        if (prodscal < 0.) prodsign2 = 1;
        // store matching edge when prodsign of elt1 and elt2 is different
        if (prodsign1 != prodsign2)
        {
          found1 = 0; found2 = 0;
          if (indA1 == indA2 || indA1 == indB2 || indA1 == indC2)
            found1 = 1;
          if (indB1 == indA2 || indB1 == indB2 || indB1 == indC2)
          {
            if (found1 == 0) found1 = 2;
            else found2 = 2;
          }
          if (indC1 == indA2 || indC1 == indB2 || indC1 == indC2)
            found2 = 3;

          if (found1 == 1) {for (E_Int p=0;p<nfld;p++) fnewp[p][nop] = ptA1[p];}
          else {for (E_Int p=0;p<nfld;p++) fnewp[p][nop] = ptB1[p];}
          if (found2 == 2) {for (E_Int p=0;p<nfld;p++) fnewp[p][nop+1] = ptB1[p];}
          else {for (E_Int p=0;p<nfld;p++) fnewp[p][nop+1] = ptC1[p];}

          cnnew1[noe] = nop+1; cnnew2[noe] = nop+2;
          noe++; nop = nop+2;
          // reallocation if neccessary
          if (nop+2 >fnew->getSize()) 
          {fnew->reAllocMat(2*fnew->getSize(),nfld); for (E_Int p=0;p<nfld;p++) fnewp[p] = fnew->begin(p+1);}
          if (noe+2 > cnnew->getSize()) 
          {cnnew->reAllocMat(2*cnnew->getSize(),2); cnnew1 = cnnew->begin(1); cnnew2 = cnnew->begin(2);}
        }
      }
    }
  }
  else if (type == 4) // QUAD 
  {
    E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2);
    E_Int* cn3 = cn->begin(3); E_Int* cn4 = cn->begin(4);
    // indAi,indBi,indCi,indDi: indices of a QUAD element i
    E_Int indA1, indB1, indC1, indD1, indA2, indB2, indC2, indD2;
    // ptA1,ptB1,ptC1,ptD1: pointers on f for points A1,B1,C1,D1
    vector<E_Float> ptA1(nfld); vector<E_Float> ptB1(nfld);
    vector<E_Float> ptC1(nfld); vector<E_Float> ptD1(nfld);

    // Compute surface vector of triangles
    K_METRIC::compUnstructSurf(*cn, "TRI", fp[posx], fp[posy], fp[posz], 
                               nsurfx, nsurfy, nsurfz, surf.begin());

    // Loop on all elements
    for (E_Int elt1 = 0; elt1 < nelts; elt1++)
    {
      indA1 = cn1[elt1]-1; indB1 = cn2[elt1]-1;
      indC1 = cn3[elt1]-1; indD1 = cn4[elt1]-1;
      for (E_Int p=0;p<nfld;p++)
      {
        ptA1[p] = fp[p][indA1]; ptB1[p] = fp[p][indB1];
        ptC1[p] = fp[p][indC1]; ptD1[p] = fp[p][indD1];
      }
      // sign of scalar product for element elt1
      prodsign1=0;
      prodscal = nsurfx[elt1]*vx + nsurfy[elt1]*vy + nsurfz[elt1]*vz;
      if (prodscal < 0.) prodsign1 = 1;
      vector<E_Int>& eltsVoisins = cEEN[elt1]; E_Int nvoisins = eltsVoisins.size();
      // Loop on all neighbours
      for (E_Int et2 = 0; et2 < nvoisins; et2++)
      {
        elt2 = eltsVoisins[et2];
        indA2 = cn1[elt2]-1; indB2 = cn2[elt2]-1; indC2 = cn3[elt2]-1; indD2 = cn4[elt2]-1;
        // sign of scalar product for neighbours
        prodsign2=0;
        prodscal = nsurfx[elt2]*vx + nsurfy[elt2]*vy + nsurfz[elt2]*vz;
        if (prodscal < 0.) prodsign2 = 1;
        // store matching edge when prodsign of elt1 and elt2 is different
        if (prodsign1 != prodsign2)
        {
          found1 = 0; found2 = 0;
          if (indA1 == indA2 || indA1 == indB2 || indA1 == indC2 || indA1 == indD2) 
            found1 = 1;
          if (indB1 == indA2 || indB1 == indB2 || indB1 == indC2 || indB1 == indD2)
          {
            if ( found1 == 0) found1 = 2;
            else found2 = 2;
          }
          if (indC1 == indA2 || indC1 == indB2 || indC1 == indC2 || indC1 == indD2)
          {
            if ( found1 == 0 ) found1 = 3; 
            else found2 = 3;
          }
          if (indD1 == indA2 || indD1 == indB2 || indD1 == indC2 || indD1 == indD2) 
            found2 = 4;

          if (found1 == 1) {for (E_Int p=0;p<nfld;p++) fnewp[p][nop] = ptA1[p];}
          else if (found1 == 2) {for (E_Int p=0;p<nfld;p++) fnewp[p][nop] = ptB1[p];}
          else {for (E_Int p=0;p<nfld;p++) fnewp[p][nop] = ptC1[p];}
          if (found2 == 2) {for (E_Int p=0;p<nfld;p++) fnewp[p][nop+1] = ptB1[p];}
          else if (found2 == 3)  {for (E_Int p=0;p<nfld;p++) fnewp[p][nop+1] = ptC1[p];}
          else {for (E_Int p=0;p<nfld;p++) fnewp[p][nop+1] = ptD1[p];}

          cnnew1[noe] = nop+1; cnnew2[noe] = nop+2;
          noe++; nop=nop+2;
          // reallocation if neccessary
          if (nop+2 > fnew->getSize()) 
          {fnew->reAllocMat(2*fnew->getSize(),nfld); for (E_Int p=0;p<nfld;p++) fnewp[p] = fnew->begin(p+1);}
          if (noe+2 > cnnew->getSize()) 
          {cnnew->reAllocMat(2*cnnew->getSize(),2); cnnew1 = cnnew->begin(1); cnnew2 = cnnew->begin(2);}
        }
      }
    }
  }

  RELEASESHAREDU(array, f, cn);

  // reallocate fnew andcnnew
  if (type != 2)
  {
    fnew->reAllocMat(nop,nfld);
    cnnew->reAllocMat(noe, 2);
  
    // clean connectivity: suppress multiple nodes, unused nodes and degenerate elements
    posx++; posy++; posz++;
    K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-10, eltType, 
                                 *fnew, *cnnew);
  }

  if (cnnew->getSize() == 0)
  {
    if (cnnew != NULL) delete cnnew;
    if (fnew != NULL) delete fnew;
    Py_INCREF(Py_None);
    return Py_None;
  }

  PyObject* tpl;
  tpl = K_ARRAY::buildArray(*fnew, varString, *cnnew, -1, "BAR");
  delete fnew; delete cnnew;
  return tpl;
}
