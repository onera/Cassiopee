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
# include "converter.h"

using namespace K_FLD;

//=============================================================================
/* Convert a vertex range (PR of structured block) into a vertex indices 
   numpy (VL). 
   The VL starts 0. */
//=============================================================================
PyObject* K_CONVERTER::PR2VL(PyObject* self, PyObject* args)
{
  PyObject* PRo; E_Int ni, nj, nk;
  PyObject* PRdo; PyObject* trfo; E_Int ni2, nj2, nk2;

  E_Int nargs = PyTuple_Size(args);
  if (nargs == 4)
  {
    if (PYPARSETUPLE_(args, O_ III_, &PRo, &ni, &nj, &nk))
    {
      PRdo = NULL;
    }
    else return NULL;
  }
  else if (PYPARSETUPLE_(args, O_ III_ OO_ III_, 
            &PRo, &ni, &nj, &nk,
            &PRdo, &trfo, &ni2, &nj2, &nk2))
  {
    if (PRdo == Py_None) PRdo = NULL;
  }
  else return NULL;

  // Check pointRange numpy
  FldArrayI* PR;
  E_Int res = K_NUMPY::getFromNumpyArray(PRo, PR);

  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "PR2VL: PointRange numpy is invalid.");
    return NULL;
  }
  // Extract imin,jmin,...
  E_Int t;
  E_Int nf = PR->getSize() * PR->getNfld();
  E_Int* p = PR->begin();
  E_Int imin=1, imax=1, jmin=1, jmax=1, kmin=1, kmax=1;
  if (nf == 2)
  { imin = p[0]; imax = p[1]; }
  else if (nf == 4)
  { imin = p[0]; imax = p[2]; jmin = p[1]; jmax = p[3]; }
  else if (nf == 6)
  { imin = p[0]; imax = p[3]; jmin = p[1]; jmax = p[4]; kmin = p[2]; kmax = p[5]; }
  else 
  {
    RELEASESHAREDN(PRo, PR);
    PyErr_SetString(PyExc_TypeError, 
                    "PR2VL: pointRange is invalid.");
    return NULL;
  }
  if (imax < imin) { t = imax; imax = imin; imin = t; }
  if (jmax < jmin) { t = jmax; jmax = jmin; jmin = t; }
  if (kmax < kmin) { t = kmax; kmax = kmin; kmin = t; }
  printf(SF_D6_ "\n",imin,imax,jmin,jmax,kmin,kmax); fflush(stdout);

  E_Int size, ind;
  //E_Int dim = 3;
  //if (nk == 1) dim = 2; // 2D

  // Compute vertex indices corresponding to pointRange
  PyObject* o;
  E_Int ii = 0;
  
  // fenetre en i
  if (imin == imax)
  {
    size = (kmax-kmin+1)*(jmax-jmin+1);
    o = K_NUMPY::buildNumpyArray(size, 1, 1);
    E_Int* p = K_NUMPY::getNumpyPtrI(o);
    for (E_Int k = kmin-1; k < kmax; k++)
    {
      for (E_Int j = jmin-1; j < jmax; j++)
      {
        ind = imin-1 + j*ni + k*ni*nj;
        p[ii] = ind; ii++;
      }
    }
    
  }
  // fenetre en j
  else if (jmin == jmax)
  {
    size = (kmax-kmin+1)*(imax-imin+1);
    o = K_NUMPY::buildNumpyArray(size, 1, 1);
    E_Int* p = K_NUMPY::getNumpyPtrI(o);
    for (E_Int k = kmin-1; k < kmax; k++)
    {
      for (E_Int i = imin-1; i < imax; i++)
      {
        ind = i + (jmin-1)*ni + k*ni*nj;
        p[ii] = ind; ii++;
      }
    }
  }
  // fenetre en k
  else if (kmin == kmax)
  {
    size = (jmax-jmin+1)*(imax-imin+1);
    o = K_NUMPY::buildNumpyArray(size, 1, 1);
    E_Int* p = K_NUMPY::getNumpyPtrI(o);
    for (E_Int j = jmin-1; j < jmax; j++)
    {
      for (E_Int i = imin-1; i < imax; i++)
      {
        ind = i + j*ni + (kmin-1)*ni*nj;
        p[ii] = ind; ii++;
      }
    }
  }
  else
  {
    printf("Warning: PR2VL: requires a 2D range.\n");
    RELEASESHAREDN(PRo, PR);
    return NULL;
  }

  RELEASESHAREDN(PRo, PR);
  PyObject* q = NULL;
    
  if (PRdo != NULL)
  {
    // Check numpy (pointRange donor)
    FldArrayI* PRdonor;
    E_Int res = K_NUMPY::getFromNumpyArray(PRdo, PRdonor);

    if (res == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "PR2VL: PointRangeDonor numpy is invalid.");
      return NULL;
    } 
    
    // Check numpy (trf)
    FldArrayI* trf;
    res = K_NUMPY::getFromNumpyArray(trfo, trf);

    if (res == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "PR2VL: Transform numpy is invalid.");
      return NULL;
    }

    // Extract imin2,imax2,...
    nf = PRdonor->getSize() * PRdonor->getNfld();
    p = PRdonor->begin();
    E_Int imin2=1, /*imax2=1,*/ jmin2=1, /*jmax2=1,*/ kmin2=1; /*kmax2=1;*/
    if (nf == 2)
    { imin2 = p[0]; /*imax2 = p[1];*/ }
    else if (nf == 4)
    { imin2 = p[0]; /*imax2 = p[2];*/ jmin2 = p[1]; /*jmax2 = p[3];*/ }
    else if (nf == 6)
    { imin2 = p[0]; /*imax2 = p[3];*/ jmin2 = p[1]; /*jmax2 = p[4];*/ kmin2 = p[2]; /*kmax2 = p[5];*/ }
    else
    {
      RELEASESHAREDN(PRdo, PRdonor);
      PyErr_SetString(PyExc_TypeError, 
                      "PR2VL: pointRangeDonor is invalid.");
      return NULL;
    }
    // revert point range if not right order
    //if (imax2 < imin2) { t = imax2; imax2 = imin2; imin2 = t; }
    //if (jmax2 < jmin2) { t = jmax2; jmax2 = jmin2; jmin2 = t; }
    //if (kmax2 < kmin2) { t = kmax2; kmax2 = kmin2; kmin2 = t; }    

    // shift matrix
    E_Int* trfp = trf->begin();
    E_Int s = trf->getSize()*trf->getNfld();
    E_Int m[9];

    for (E_Int c = 0; c < s; c++)
    {
      if (trfp[c] == 1)
      { m[0+3*c] = 1; m[1+3*c] = 0; m[2+3*c] = 0; }
      else if (trfp[c] == -1)
      { m[0+3*c] = -1; m[1+3*c] = 0; m[2+3*c] = 0; }
      else if (trfp[c] == 2)
      { m[0+3*c] = 0; m[1+3*c] = 1; m[2+3*c] = 0; }
      else if (trfp[0] == -2)
      { m[0+3*c] = 0; m[1+3*c] = -1; m[2+3*c] = 0; }
      else if (trfp[c] == 3)
      { m[0+3*c] = 0; m[1+3*c] = 0; m[2+3*c] = 1; }
      else if (trfp[c] == -3)
      { m[0+3*c] = 0; m[1+3*c] = 0; m[2+3*c] = -1; }
    }
    if (s == 2)
    {
      { m[0+3*2] = 0; m[1+3*2] = 0; m[2+3*2] = 1; }
    }

    // Compute donor indices
    E_Int ip,jp,kp;
    ii = 0;
    
    // fenetre en i
    if (imin == imax)
    {
      size = (kmax-kmin+1)*(jmax-jmin+1);
      q = K_NUMPY::buildNumpyArray(size, 1, 1);
      E_Int* p = K_NUMPY::getNumpyPtrI(q);
      E_Int i = imin-1;
      for (E_Int k = kmin-1; k < kmax; k++)
      {
        for (E_Int j = jmin-1; j < jmax; j++)
        {
          ip = imin2 + m[0+3*0]*(i-imin)+m[0+3*1]*(j-jmin)+m[0+3*2]*(k-kmin);
          jp = jmin2 + m[1+3*0]*(i-imin)+m[1+3*1]*(j-jmin)+m[1+3*2]*(k-kmin);
          kp = kmin2 + m[2+3*0]*(i-imin)+m[2+3*1]*(j-jmin)+m[2+3*2]*(k-kmin);
          ind = ip + jp*ni2 + kp*ni2*nj2;
          p[ii] = ind; ii++;
        }
      }
    }
    // fenetre en j
    else if (jmin == jmax)
    {
      size = (kmax-kmin+1)*(imax-imin+1);
      q = K_NUMPY::buildNumpyArray(size, 1, 1);
      E_Int* p = K_NUMPY::getNumpyPtrI(q);
      //E_Int j = jmin-1;
      for (E_Int k = kmin-1; k < kmax; k++)
      {
        for (E_Int i = imin-1; i < imax; i++)
        {
          ind = i + (jmin-1)*ni + k*ni*nj;
          p[ii] = ind; ii++;
        }
      }
    }
    // fenetre en k
    else if (kmin == kmax)
    {
      size = (jmax-jmin+1)*(imax-imin+1);
      q = K_NUMPY::buildNumpyArray(size, 1, 1);
      E_Int* p = K_NUMPY::getNumpyPtrI(q);
      for (E_Int j = jmin-1; j < jmax; j++)
      {
        for (E_Int i = imin-1; i < imax; i++)
        {
          ind = i + j*ni + (kmin-1)*ni*nj;
          p[ii] = ind; ii++;
        }
      }
    }

    RELEASESHAREDN(PRdo, PRdonor);
    RELEASESHAREDN(trfo, trf);
  }

  if (PRdo == NULL) return o;
  else return Py_BuildValue("OO",o,q); 
}
