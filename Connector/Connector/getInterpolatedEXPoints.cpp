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

#include "connector.h"
using namespace K_FLD;
using namespace std;
using namespace K_FUNC;

//=============================================================================
/* Determine les points EX interpoles a partir du cellN */
//=============================================================================
PyObject* K_CONNECTOR::getEXPoints(PyObject* self, PyObject* args)
{
  PyObject *coordArray, *cellNArray; // domaine a interpoler
  if (!PyArg_ParseTuple(args, "OO",
                        &coordArray, &cellNArray))
  {
      return NULL;
  }
  
  /*--------------------------------------------------*/
  /* Extraction des infos sur le domaine a interpoler */
  /*--------------------------------------------------*/
  E_Int im, jm, km;
  FldArrayF* field; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(coordArray, varString, field, 
                                    im, jm, km, cn, eltType, true); 
  if (res != 1)
  {
    if (res == 2) RELEASESHAREDU(coordArray, field,cn);
    PyErr_SetString(PyExc_TypeError,
                    "getEXPoints: first argument must be structured.");
    return NULL;
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getEXPoints: first array must contain (x,y,z).");
    RELEASESHAREDS(coordArray, field);
    return NULL;
  }
  posx++; posy++; posz++;
  E_Int imc, jmc, kmc;
  FldArrayF* field1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  res = K_ARRAY::getFromArray(cellNArray, varString1, field1, 
                              imc, jmc, kmc, cn1, eltType1, true); 
  if (res != 1)
  {
    RELEASESHAREDS(coordArray, field);
    RELEASESHAREDB(res, cellNArray, field1, cn1);
    PyErr_SetString(PyExc_TypeError,
                    "getEXPoints: 2nd argument must be structured.");
    return NULL;
  }
  E_Int posc = K_ARRAY::isCellNatureField2Present(varString1);
  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getEXPoints: 2nd array must contain cellN.");
    RELEASESHAREDS(coordArray, field);
    RELEASESHAREDS(cellNArray, field1);
    return NULL;
  }
  posc++;
  if (E_max(1,im-1) != imc || E_max(1,jm-1) != jmc || E_max(1,km-1) != kmc)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getEXPoints: 1st and 2nd arrays must be located at nodes and centers respectively.");
    RELEASESHAREDS(coordArray, field);
    RELEASESHAREDS(cellNArray, field1);
    return NULL;   
  }
  /* fin des verifs*/
  E_Float* xn = field->begin(posx);
  E_Float* yn = field->begin(posy);
  E_Float* zn = field->begin(posz);

  /* Determination des pts interpoles */
  E_Int ncells = field1->getSize();
  E_Float* cellNp = field1->begin(posc);
  FldArrayI interpolatedCellArray(ncells);
  E_Int c = 0;
  for (E_Int ind = 0; ind < ncells; ind++)
  {
    if (cellNp[ind] == 2.) {interpolatedCellArray[c] = ind; c++;}
  }
  interpolatedCellArray.resize(c);

  E_Int nbEX = 0;
  E_Int inci, incj, inck;
  E_Int ind1, ind2, ind3, ind4;
  E_Int indc1, indc2, indc3, indc4, indc5, indc6;
  E_Int imjm = im*jm;
  E_Int imcjmc = imc*jmc;
  E_Int im1, ip1, jm1, jp1, km1, kp1;
  // coordEX
  FldArrayF* coordEX = new FldArrayF(6*ncells,7);
  E_Float* xt = coordEX->begin(1);
  E_Float* yt = coordEX->begin(2);
  E_Float* zt = coordEX->begin(3);
  E_Float* indirEX = coordEX->begin(4);
  E_Float* indirEX2 = coordEX->begin(5);
  E_Float* indirNode = coordEX->begin(6);
  E_Float* dirEX = coordEX->begin(7);
  
  // determination des centres d interfaces EX pres des pts masques
  for (E_Int nocell = 0; nocell < interpolatedCellArray.getSize(); nocell++)
  {
    E_Int indcell = interpolatedCellArray[nocell];
    E_Int kcell = indcell/imcjmc;
    E_Int jcell = (indcell - kcell*imcjmc)/imc;
    E_Int icell = indcell - jcell * imc - kcell * imcjmc; 
    //recherche d un voisin sur les 6 qui est masque
    im1 = K_FUNC::E_max(0,icell-1); ip1 = K_FUNC::E_min(icell+1,imc-1);
    jm1 = K_FUNC::E_max(0,jcell-1); jp1 = K_FUNC::E_min(jcell+1,jmc-1);
    km1 = K_FUNC::E_max(0,kcell-1); kp1 = K_FUNC::E_min(kcell+1,kmc-1);
    indc1 = im1 + jcell * imc + kcell * imcjmc;
    indc2 = ip1 + jcell * imc + kcell * imcjmc;
    indc3 = icell + jm1 * imc + kcell * imcjmc;
    indc4 = icell + jp1 * imc + kcell * imcjmc;
    indc5 = icell + jcell * imc + km1 * imcjmc;
    indc6 = icell + jcell * imc + kp1 * imcjmc;    

    if ((cellNp[indc1] == 0. && cellNp[indc2] == 1.) || 
        (icell == 0 && cellNp[indc2] == 1.))// voisin i-1 masque ou pt frontiere
    {
      inci = icell;
      jm1 = jcell; jp1 = jm1+1; km1 = kcell; kp1 = km1+1;
      ind1 = inci + jm1*im + km1*imjm;
      ind2 = inci + jp1*im + km1*imjm;
      ind3 = inci + jp1*im + kp1*imjm;
      ind4 = inci + jm1*im + kp1*imjm;
      xt[nbEX] = 0.25*(xn[ind1]+xn[ind2]+xn[ind3]+xn[ind4]);
      yt[nbEX] = 0.25*(yn[ind1]+yn[ind2]+yn[ind3]+yn[ind4]);
      zt[nbEX] = 0.25*(zn[ind1]+zn[ind2]+zn[ind3]+zn[ind4]);
      indirEX[nbEX] = E_Float(indcell); indirEX2[nbEX] = E_Float(indc1);
      indirNode[nbEX] = ind1; dirEX[nbEX] = 1;
      nbEX++;
    }
    if ((cellNp[indc2] == 0. && cellNp[indc1] == 1.) ||
        (icell == imc-1 && cellNp[indc1] == 1.))// voisin i+1 masque ou pt frontiere
    {
      inci = icell+1;
      jm1 = jcell; jp1 = jm1+1; km1 = kcell; kp1 = km1+1;
      ind1 = inci + jm1*im + km1*imjm;
      ind2 = inci + jp1*im + km1*imjm;
      ind3 = inci + jp1*im + kp1*imjm;
      ind4 = inci + jm1*im + kp1*imjm;
      xt[nbEX] = 0.25*(xn[ind1]+xn[ind2]+xn[ind3]+xn[ind4]);
      yt[nbEX] = 0.25*(yn[ind1]+yn[ind2]+yn[ind3]+yn[ind4]);
      zt[nbEX] = 0.25*(zn[ind1]+zn[ind2]+zn[ind3]+zn[ind4]);
      indirEX[nbEX] = E_Float(indcell); indirEX2[nbEX] = E_Float(indc2);
      indirNode[nbEX] = ind1; dirEX[nbEX] = 0;
      nbEX++;
    }
    if ((cellNp[indc3] == 0. && cellNp[indc4] == 1.)  || 
        (jcell == 0 && cellNp[indc4] == 1.))// voisin j-1 masque ou frontiere
    {
      incj = jcell*im;
      im1 = icell; ip1 = im1+1; km1 = kcell; kp1 = km1+1;
      ind1 = im1 + km1*imjm + incj;
      ind2 = ip1 + km1*imjm + incj;
      ind3 = im1 + kp1*imjm + incj;
      ind4 = ip1 + kp1*imjm + incj;
      xt[nbEX] = 0.25*(xn[ind1]+xn[ind2]+xn[ind3]+xn[ind4]);
      yt[nbEX] = 0.25*(yn[ind1]+yn[ind2]+yn[ind3]+yn[ind4]);
      zt[nbEX] = 0.25*(zn[ind1]+zn[ind2]+zn[ind3]+zn[ind4]);
      indirEX[nbEX] = E_Float(indcell);indirEX2[nbEX] = E_Float(indc3);
      indirNode[nbEX] = ind1; dirEX[nbEX] = 3;
      nbEX++;
    }

    if ((cellNp[indc4] == 0.  && cellNp[indc3] == 1.) || 
        (jcell == jmc-1 && cellNp[indc3] == 1.))// voisin j+1 masque ou frontiere
    {
      incj = (jcell+1)*im;
      im1 = icell; ip1 = im1+1; km1 = kcell; kp1 = km1+1;
      ind1 = im1 + km1*imjm + incj;
      ind2 = ip1 + km1*imjm + incj;
      ind3 = im1 + kp1*imjm + incj;
      ind4 = ip1 + kp1*imjm + incj;
      xt[nbEX] = 0.25*(xn[ind1]+xn[ind2]+xn[ind3]+xn[ind4]);
      yt[nbEX] = 0.25*(yn[ind1]+yn[ind2]+yn[ind3]+yn[ind4]);
      zt[nbEX] = 0.25*(zn[ind1]+zn[ind2]+zn[ind3]+zn[ind4]);
      indirEX[nbEX] = E_Float(indcell);indirEX2[nbEX] = E_Float(indc4);
      indirNode[nbEX] = ind1; dirEX[nbEX] = 2;

      nbEX++;
    }
    if ((cellNp[indc5] == 0.  && cellNp[indc6] == 1.) || 
        (kcell == 0 && cellNp[indc6] == 1.))// voisin k-1 masque ou frontiere
    {
      inck = kcell*imjm;
      im1 = icell; ip1 = im1+1; jm1 = jcell; jp1 = jm1+1;
      ind1 = im1 + jm1*im + inck;
      ind2 = ip1 + jm1*im + inck;
      ind3 = im1 + jp1*im + inck;
      ind4 = ip1 + jp1*im + inck; 
      xt[nbEX] = 0.25*(xn[ind1]+xn[ind2]+xn[ind3]+xn[ind4]);
      yt[nbEX] = 0.25*(yn[ind1]+yn[ind2]+yn[ind3]+yn[ind4]);
      zt[nbEX] = 0.25*(zn[ind1]+zn[ind2]+zn[ind3]+zn[ind4]);
      indirEX[nbEX] = E_Float(indcell);indirEX2[nbEX] = E_Float(indc5);
      indirNode[nbEX] = ind1; dirEX[nbEX] = 5;
      nbEX++;
    }
    if ((cellNp[indc6] == 0.  && cellNp[indc5] == 1.) || 
        (kcell == kmc-1 && cellNp[indc5] == 1.) )// voisin k+1 masque ou frontiere
    {
      inck = (kcell+1)*imjm;
      im1 = icell; ip1 = im1+1; jm1 = jcell; jp1 = jm1+1;
      ind1 = im1 + jm1*im + inck;
      ind2 = ip1 + jm1*im + inck;
      ind3 = im1 + jp1*im + inck;
      ind4 = ip1 + jp1*im + inck;        
      xt[nbEX] = 0.25*(xn[ind1]+xn[ind2]+xn[ind3]+xn[ind4]);
      yt[nbEX] = 0.25*(yn[ind1]+yn[ind2]+yn[ind3]+yn[ind4]);
      zt[nbEX] = 0.25*(zn[ind1]+zn[ind2]+zn[ind3]+zn[ind4]);
      indirEX[nbEX] = E_Float(indcell);indirEX2[nbEX] = E_Float(indc6);
      indirNode[nbEX] = ind1; dirEX[nbEX] = 4;
      nbEX++;
    }
  }
  coordEX->reAllocMat(nbEX,7);

  RELEASESHAREDS(coordArray, field);
  RELEASESHAREDS(cellNArray, field1); 
  FldArrayI* cnl = new FldArrayI(0);
  PyObject* tpl = K_ARRAY::buildArray(*coordEX, 
                                      "x,y,z,indcell1,indcell2,nodemin,EXdir", 
                                      *cnl, -1, "NODE", false);
  delete coordEX; delete cnl;
  return tpl;
}
