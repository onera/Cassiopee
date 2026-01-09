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
# include "converter.h"

using namespace K_FLD;

//=============================================================================
/* Convert a PointList to a list of PointLists (PLs) one for each face of 
   a structured grid.
   Global index of faces for structured grids start 0 and contains i interfaces
   then j interfaces, then k interfaces */
//=============================================================================
PyObject* K_CONVERTER::pointList2SPL(PyObject* self, PyObject* args)
{
  PyObject* PLarray; PyObject* PLDarray;
  E_Int ni, nj, nk;
  if (!PYPARSETUPLE_(args, OO_ III_, &PLarray, &PLDarray, &ni, &nj, &nk)) return NULL;

  // Check numpy (pointlist)
  FldArrayI* PL;
  E_Int res1 = K_NUMPY::getFromPointList(PLarray, PL);

  if (res1 == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "pointList2SPL: PL numpy is invalid.");
    return NULL;
  }
  // Chekc numpy (pointList donor)
  FldArrayI* PLD;
  E_Int res2 = K_NUMPY::getFromPointList(PLDarray, PLD);

  if (res2 == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "pointList2SPL: PLD numpy is invalid.");
    return NULL;
  }

  E_Int nf = PL->getSize();
  E_Int* p = PL->begin();
  E_Int* pd = PLD->begin();
  E_Int ni1 = K_FUNC::E_max(ni-1,1);
  E_Int nj1 = K_FUNC::E_max(nj-1,1);
  E_Int nk1 = K_FUNC::E_max(nk-1,1);
  
  E_Int ninti = ni*nj1*nk1;
  E_Int nintj = ni1*nj*nk1;
  E_Int i,j,k,ind;

  E_Int nimin = 0; // nbre d'interfaces imin
  E_Int nimax = 0;
  E_Int njmin = 0;
  E_Int njmax = 0;
  E_Int nkmin = 0;
  E_Int nkmax = 0;

  // recherche par type d'interface
  for (E_Int n = 0; n < nf; n++)
  {
    ind = p[n];
    if (ind >= ninti+nintj) // interface k
    {
      ind = ind-ninti-nintj;
      k = ind/(ni1*nj1);
      if (k == 0) nkmin++;
      else if (k == nk1) nkmax++;
    }
    else if (ind >= ninti) // interface j
    {
      ind = ind-ninti;
      k = ind/(ni1*nj);
      j = (ind-k*ni1*nj)/ni1;
      if (j == 0) njmin++;
      else if (j == nj1) njmax++;
    }
    else
    {
      k = ind/(ni*nj1);
      j = (ind-k*ni*nj1)/ni;
      i = ind-j*ni-k*ni*nj1;
      if (i == 0) nimin++;
      else if (i == ni1) nimax++;
    }
  }

  // allocate numpys
  PyObject* pimin = NULL;
  PyObject* pimax = NULL;
  PyObject* pjmin = NULL;
  PyObject* pjmax = NULL;
  PyObject* pkmin = NULL;
  PyObject* pkmax = NULL;
  PyObject* pdimin = NULL;
  PyObject* pdimax = NULL;
  PyObject* pdjmin = NULL;
  PyObject* pdjmax = NULL;
  PyObject* pdkmin = NULL;
  PyObject* pdkmax = NULL;

  E_Int* imin=NULL; E_Int* imax=NULL;
  E_Int* jmin=NULL; E_Int* jmax=NULL;
  E_Int* kmin=NULL; E_Int* kmax=NULL;
  E_Int* dimin=NULL; E_Int* dimax=NULL;
  E_Int* djmin=NULL; E_Int* djmax=NULL;
  E_Int* dkmin=NULL; E_Int* dkmax=NULL;
  
  if (nimin > 0) 
  {
    pimin = K_NUMPY::buildNumpyArray(nimin, 1, 1, 0);
    imin = K_NUMPY::getNumpyPtrI(pimin);
    pdimin = K_NUMPY::buildNumpyArray(nimin, 1, 1, 0);
    dimin = K_NUMPY::getNumpyPtrI(pdimin);
  }
  if (nimax > 0) 
  {
    pimax = K_NUMPY::buildNumpyArray(nimax, 1, 1, 0);
    imax = K_NUMPY::getNumpyPtrI(pimax);
    pdimax = K_NUMPY::buildNumpyArray(nimax, 1, 1, 0);
    dimax = K_NUMPY::getNumpyPtrI(pdimax);
  }
  if (njmin > 0) 
  {
    pjmin = K_NUMPY::buildNumpyArray(njmin, 1, 1, 0);
    jmin = K_NUMPY::getNumpyPtrI(pjmin);
    pdjmin = K_NUMPY::buildNumpyArray(njmin, 1, 1, 0);
    djmin = K_NUMPY::getNumpyPtrI(pdjmin);
  }
  if (njmax > 0) 
  {
    pjmax = K_NUMPY::buildNumpyArray(njmax, 1, 1, 0);
    jmax = K_NUMPY::getNumpyPtrI(pjmax);
    pdjmax = K_NUMPY::buildNumpyArray(njmax, 1, 1, 0);
    djmax = K_NUMPY::getNumpyPtrI(pdjmax);
  }
  if (nkmin > 0) 
  {
    pkmin = K_NUMPY::buildNumpyArray(nkmin, 1, 1, 0);
    kmin = K_NUMPY::getNumpyPtrI(pkmin);
    pdkmin = K_NUMPY::buildNumpyArray(nkmin, 1, 1, 0);
    dkmin = K_NUMPY::getNumpyPtrI(pdkmin);
  }
  if (nkmax > 0)
  {
    pkmax = K_NUMPY::buildNumpyArray(nkmax, 1, 1, 0);
    kmax = K_NUMPY::getNumpyPtrI(pkmax);
    pdkmax = K_NUMPY::buildNumpyArray(nkmax, 1, 1, 0);
    dkmax = K_NUMPY::getNumpyPtrI(pdkmax);
  }

  //printf("%d %d %d %d %d %d\n", nimin,nimax,njmin,njmax,nkmin,nkmax);

  nimin = 0; nimax = 0; njmin = 0; njmax = 0; nkmin = 0; nkmax = 0;
  for (E_Int n = 0; n < nf; n++)
  {
    ind = p[n];
    if (ind >= ninti+nintj) // interface k
    {
      // ind = i+j*ni1+k*ni1*nj1;
      ind = ind-ninti-nintj;
      k = ind/(ni1*nj1);
      if (k == 0) { kmin[nkmin] = p[n]; dkmin[nkmin] = pd[n]; nkmin++; }
      else if (k == nk1) { kmax[nkmax] = p[n]; dkmax[nkmax] = pd[n]; nkmax++; }
    }
    else if (ind >= ninti) // interface j
    {
      // ind = i+j*ni1+k*ni1*nj;
      ind = ind-ninti;
      k = ind/(ni1*nj);
      j = (ind-k*ni1*nj)/ni1;
      if (j == 0) { jmin[njmin] = p[n]; djmin[njmin] = pd[n]; njmin++; }
      else if (j == nj1) { jmax[njmax] = p[n]; djmax[njmax] = pd[n]; njmax++; } 
    }
    else // interface i
    {
      // ind = i+j*ni+k*ni*nj1;
      k = ind/(ni*nj1);
      j = (ind-k*ni*nj1)/ni;
      i = ind-j*ni-k*ni*nj1;
      if (i == 0) { imin[nimin] = p[n]; dimin[nimin] = pd[n]; nimin++; }
      else if (i == ni1) { imax[nimax] = p[n]; dimax[nimax] = pd[n]; nimax++; }
    }
  }

  if (pimin == NULL) { Py_INCREF(Py_None); pimin = Py_None; Py_INCREF(Py_None); pdimin = Py_None;}
  if (pimax == NULL) { Py_INCREF(Py_None); pimax = Py_None; Py_INCREF(Py_None); pdimax = Py_None;}
  if (pjmin == NULL) { Py_INCREF(Py_None); pjmin = Py_None; Py_INCREF(Py_None); pdjmin = Py_None;}
  if (pjmax == NULL) { Py_INCREF(Py_None); pjmax = Py_None; Py_INCREF(Py_None); pdjmax = Py_None;}
  if (pkmin == NULL) { Py_INCREF(Py_None); pkmin = Py_None; Py_INCREF(Py_None); pdkmin = Py_None;}
  if (pkmax == NULL) { Py_INCREF(Py_None); pkmax = Py_None; Py_INCREF(Py_None); pdkmax = Py_None;}

  PyObject* tpl = Py_BuildValue("[OOOOOOOOOOOO]", pimin, pimax, pjmin, pjmax, pkmin, pkmax,
                                pdimin, pdimax, pdjmin, pdjmax, pdkmin, pdkmax);

  RELEASESHAREDN(PLarray, PL);
  RELEASESHAREDN(PLDarray, PLD);

  return tpl;
}
