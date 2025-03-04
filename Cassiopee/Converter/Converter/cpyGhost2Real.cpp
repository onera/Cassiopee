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

// cpyGhost2Real values of an array
// cpyReal2Ghost values of an array
# include "converter.h"

//=============================================================================
/* Copy array with ghost values in an array without ghost values */
//=============================================================================
PyObject* K_CONVERTER::cpyGhost2Real(PyObject* self, PyObject* args)
{
  IMPORTNUMPY;
  PyObject* arrayR; PyObject* arrayG;
  E_Int Im=0, Jm=0, Km=0, D;
  if (!PYPARSETUPLE_(args, OO_ IIII_, &arrayR, &arrayG, &D, &Im, &Jm, &Km)) return NULL;
  
  // Ecriture de la valeur
  E_Int im = Im; // indices reels
  E_Int jm = Jm;
  E_Int km = Km;
  E_Int imjm = im*jm;
  E_Int d = D;
  E_Int dim = 3;
  E_Int im2 = im+2*d; // indices ghost
  E_Int im2jm2 = im2*(jm+2*d);
  E_Int indreal, indghost;
  PyArrayObject* arR = (PyArrayObject*)arrayR;
  PyArrayObject* arG = (PyArrayObject*)arrayG;
  E_Float* dataG = (E_Float*)PyArray_DATA(arG);
  E_Float* dataR = (E_Float*)PyArray_DATA(arR);
  if (km == 0) dim = 2; 
  if (jm == 0) dim = 1; 

  switch (dim)
  {
    case 3:
#pragma omp parallel for default(shared) private(indghost,indreal)
      for (E_Int k = d; k < km+d; k++)
        for (E_Int j = d; j < jm+d; j++)
          for (E_Int i = d; i < im+d; i++)
          {
            indghost = i+j*im2+k*im2jm2;
            indreal = (i-d)+(j-d)*im+(k-d)*imjm;
            dataR[indreal] = dataG[indghost];
          }
      break;

    case 2:
#pragma omp parallel for default(shared) private(indghost,indreal)
      for (E_Int j = d; j < jm+d; j++)
        for (E_Int i = d; i < im+d; i++)
        {
          indghost = i + j*im2;
          indreal = (i-d)+(j-d)*im;
          dataR[indreal] = dataG[indghost];
        }
      break;
  
    case 1:
      for (E_Int i = d; i < im+d; i++) dataR[i-d] = dataG[i];
      break;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
/* Copy array with ghost values in an array without ghost values */
//=============================================================================
PyObject* K_CONVERTER::cpyReal2Ghost(PyObject* self, 
                                     PyObject* args)
{
  IMPORTNUMPY;
  PyObject* arrayR;
  PyObject* arrayG;
  E_Int Im=0, Jm=0, Km=0, D;
  if (!PYPARSETUPLE_(args, OO_ IIII_, &arrayG, &arrayR, &D , &Im, &Jm, &Km)) return NULL;

  // Ecriture de la valeur
  E_Int dim = 3;
  E_Int im = Im;
  E_Int jm = Jm;
  E_Int km = Km;
  E_Int imjm = im*jm;
  E_Int d  = D;
  E_Int im2 = im+2*d;
  E_Int im2jm2 = im2*(jm+2*d);
  E_Int indreal, indghost, indadj;
  PyArrayObject* arR = (PyArrayObject*)arrayR;
  PyArrayObject* arG = (PyArrayObject*)arrayG;
  E_Float* dataG = (E_Float*)PyArray_DATA(arG);
  E_Float* dataR = (E_Float*)PyArray_DATA(arR);
  if (km == 0) dim = 2; 
  if (jm == 0) dim = 1;
  
  switch (dim)
  {
    case 3:
      // copy real values
#pragma omp parallel default(shared) private(indghost,indreal,indadj)
    {
#pragma omp for
      for (E_Int k = d; k < km+d; k++)
        for (E_Int j = d; j < jm+d; j++)
          for (E_Int i = d; i < im+d; i++)
          {
            indghost = i+j*im2+k*im2jm2;
            indreal = (i-d)+(j-d)*im+(k-d)*imjm;
            dataG[indghost] = dataR[indreal];
          }
      // fill ghost values
      // direction k
#pragma omp for
      for (E_Int j = 0; j < jm; j++)
        for (E_Int i = 0; i < im; i++)
        {
          for (E_Int k = -1; k > -d-1; k--)
          {
            indghost = i+d+(j+d)*im2+(k+d)*im2jm2;
            indadj = i+d+(j+d)*im2+(k+1+d)*im2jm2;
            dataG[indghost] = dataG[indadj];
          }
          for (E_Int k = km; k < km+d; k++)
          {
            indghost = i+d+(j+d)*im2+(k+d)*im2jm2;
            indadj = i+d+(j+d)*im2+(k-1+d)*im2jm2;
            dataG[indghost] = dataG[indadj];
          }
        }
      // direction j
#pragma omp for
      for (E_Int k = -d; k < km+d; k++)
        for (E_Int i = 0; i < im; i++)
        {
          for (E_Int j = -1; j > -d-1; j--)
          {
            indghost = i+d+(j+d)*im2+(k+d)*im2jm2;
            indadj = i+d+(j+d+1)*im2+(k+d)*im2jm2;
            dataG[indghost] = dataG[indadj];
          }
          for (E_Int j= jm; j < jm+d; j++)
          {
            indghost = i+d+(j+d)*im2+(k+d)*im2jm2;
            indadj = i+d+(j-1+d)*im2+(k+d)*im2jm2;
            dataG[indghost] = dataG[indadj];
          }
        }
      // direction i
#pragma omp for
      for (E_Int k = -d; k < km+d; k++)
        for (E_Int j = -d; j < jm+d; j++)
        {
          for (E_Int i = -1; i > -d-1; i--)
          {
            indghost = i+d+(j+d)*im2+(k+d)*im2jm2;
            indadj = i+1+d+(j+d)*im2+(k+d)*im2jm2;
            dataG[indghost] = dataG[indadj];
          }
          for (E_Int i= im; i < im+d; i++)
          {
            indghost = i+d+(j+d)*im2+(k+d)*im2jm2;
            indadj = i-1+d+(j+d)*im2+(k+d)*im2jm2;
            dataG[indghost] = dataG[indadj];
          }
        }
    }
    break;
      
    case 2:
    {
      // copy real values
#pragma omp parallel for default(shared) private(indghost,indreal)
      for (E_Int j = d; j < jm+d; j++)
        for (E_Int i = d; i < im+d; i++)
        {
          indghost = i + j*im2;
          indreal = (i-d)+(j-d)*im;
          dataG[indghost] = dataR[indreal];
        }
      // fill ghost values
      // direction j
      for(E_Int i = 0; i < im; i++)
      {
        for (E_Int j = -1; j > -d-1; j--)
        {
          indghost = i+d+(j+d)*im2;
          indadj = i+d+(j+1+d)*im2;
          dataG[indghost] = dataG[indadj];
        }
        for (E_Int j = jm; j < jm+d; j++)
        {
          indghost = i+d+(j+d)*im2;
          indadj = i+d+(j-1+d)*im2;
          dataG[indghost] = dataG[indadj];
        }
      }
      // direction i (including corners)
      for(E_Int j = -d; j < jm+d; j++)
      {
        for (E_Int i = -1; i > -d-1; i--)
        {
          indghost = i+d+(j+d)*im2;
          indadj = i+1+d+(j+d)*im2;
          dataG[indghost] = dataG[indadj];
        }
        for (E_Int i = im; i < im+d; i++)
        {
          indghost = i+d+(j+d)*im2;
          indadj = i-1+d+(j+d)*im2;
          dataG[indghost] = dataG[indadj];
        }
      }
    }
    break;
 
    case 1:
      // copy real values
      for (E_Int i = d; i < im+d; i++) dataG[i] = dataR[i-d];
      // fill ghost values
      // direction i min
      for(E_Int i = -d; i < 0; i++) dataG[i+d] = dataG[d];
      // direction i max
      for(E_Int i = im; i < im+d; i++) dataG[i+d] = dataG[d+im-1];
      break;
  }
  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
// Copy connect array to connect pyTree
//=============================================================================
PyObject* K_CONVERTER::cpyConnectP2ConnectA(PyObject* self, PyObject* args)
{
  IMPORTNUMPY;
  PyObject* cA; // connectivite array
  PyObject* cP1, *cP2; // connectivites PyTree (2 connectivite si NGON, 1 sinon)
  E_Int stype=0, ne=0; // stype: nb noeuds/elmt, ne: nb elements
  E_Int nfaces=0, nelts=0; // nombre de faces et nombre d elements (valent -1 quand element != NGON)
  if (!PYPARSETUPLE_(args, OOO_ IIII_, &cA, &cP1, &cP2, &stype, &ne, &nfaces, &nelts)) return NULL;

  PyArrayObject* crA = (PyArrayObject*)cA;
  PyArrayObject* crP1 = (PyArrayObject*)cP1;
  E_Int* dataA = (E_Int*)PyArray_DATA(crA);
  E_Int* dataP1 = (E_Int*)PyArray_DATA(crP1);
  
  if (nfaces == -1)  // Elements basiques
  {
#pragma omp parallel for default(shared)
    for (E_Int p = 0; p < ne; p++)
      for (E_Int k = 0; k < stype; k++)
        dataA[p+k*ne] = dataP1[k+stype*p];
  }
  else  // Elements NGon
  {
    PyArrayObject* crP2 = (PyArrayObject*)cP2;
    E_Int* dataP2 = (E_Int*)PyArray_DATA(crP2);
    E_Int st1 = PyArray_DIMS(crP1)[0]; // dimension de la  connectivite cP1
    E_Int st2 = PyArray_DIMS(crP2)[0]; // dimension de la  connectivite cP2

    dataA[0] = nfaces;
    dataA[1] = st1;
    dataA += 2;
#pragma omp parallel for default(shared)
    for (E_Int p = 0; p < st1; p++) dataA[p] = dataP1[p];
    dataA[st1] = nelts;
    dataA[st1+1] = st2;
    dataA += st1+2;
#pragma omp parallel for default(shared)
    for (E_Int p = 0; p < st2; p++) dataA[p] = dataP2[p];
  }
  
  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
// Copy connect array to connect pyTree for NGON2 (CGNS v4)
//=============================================================================
PyObject* K_CONVERTER::cpyConnectP2ConnectA2(PyObject* self, PyObject* args)
{
  IMPORTNUMPY;
  PyObject* cA; // connectivite array
  PyObject* cP1, *cP2; // connectivites PyTree (2 connectivite si NGON, 1 sinon)
  PyObject* off1, *off2; // offsets
  E_Int stype=0, ne=0; // stype: nb noeuds/elmt, ne: nb elements
  E_Int nfaces=0, nelts=0; // nombre de faces et nombre d elements (valent -1 quand element != NGON)
  if (!PYPARSETUPLE_(args, OOO_ IIII_ OO_, 
    &cA, &cP1, &cP2, &stype, &ne, &nfaces, &nelts, &off1, &off2)) return NULL;

  PyArrayObject* crA = (PyArrayObject*)cA;
  PyArrayObject* crP1 = (PyArrayObject*)cP1;
  E_Int* dataA = (E_Int*)PyArray_DATA(crA);
  E_Int* dataP1 = (E_Int*)PyArray_DATA(crP1);
  
  PyArrayObject* crP2 = (PyArrayObject*)cP2;
  E_Int* dataP2 = (E_Int*)PyArray_DATA(crP2);
  E_Int st1 = PyArray_DIMS(crP1)[0]; // dimension de la  connectivite cP1
  E_Int st2 = PyArray_DIMS(crP2)[0]; // dimension de la  connectivite cP2

  PyArrayObject* crOff1 = (PyArrayObject*)off1;
  PyArrayObject* crOff2 = (PyArrayObject*)off2;
  E_Int* poff1 = (E_Int*)PyArray_DATA(crOff1);
  E_Int* poff2 = (E_Int*)PyArray_DATA(crOff2);

  /*
  dataA[0] = nfaces;
  dataA[1] = st1;
  dataA += 2;
  for (E_Int p = 0; p < nfaces; p++)
  {
    size = poff1[p+1]-poff1[p];
    dataA[0] = size;
    for (E_Int j = 0; j < size; j++) dataA[j+1] = dataP1[j];
    dataA += size+1;
    dataP1 += size;
  }

  dataA[0] = nelts;
  dataA[1] = st2;
  dataA += 2;
  for (E_Int p = 0; p < nelts; p++) 
  {
    size = poff2[p+1]-poff2[p];
    dataA[0] = size;
    for (E_Int j = 0; j < size; j++) dataA[j+1] = dataP2[j];
    dataA += size+1;
    dataP2 += size;
  }
  */

  dataA[0] = nfaces;
  dataA[1] = st1+nfaces;
  dataA += 2;

#pragma omp parallel
  {
    E_Int size, p1;
    #pragma omp for
    for (E_Int p = 0; p < nfaces; p++)
    {
      p1 = poff1[p];
      size = poff1[p+1]-p1;
      dataA[p1+p] = size;
      for (E_Int j = 0; j < size; j++) dataA[p1+p+j+1] = dataP1[p1+j];
    }
  }

  dataA[st1+nfaces] = nelts;
  dataA[st1+nfaces+1] = st2+nelts;
  dataA += st1+nfaces+2;

#pragma omp parallel
  {
    E_Int size, p1;
    #pragma omp for
    for (E_Int p = 0; p < nelts; p++)
    {
      p1 = poff2[p];
      size = poff2[p+1]-p1;
      dataA[p1+p] = size;
      for (E_Int j = 0; j < size; j++) dataA[p1+p+j+1] = std::abs(dataP2[p1+j]); // rip signed face
    }
  }

  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
// Copy connect array to connect pyTree
//=============================================================================
PyObject* K_CONVERTER::cpyConnectA2ConnectP(PyObject* self, PyObject* args)
{
  IMPORTNUMPY;
  PyObject* cA; PyObject* cP;
  E_Int stype=0, ne=0;
  if (!PYPARSETUPLE_(args, OO_ II_, &cA, &cP, &stype, &ne)) return NULL;

  PyArrayObject* crA = (PyArrayObject*)cA;
  PyArrayObject* crP = (PyArrayObject*)cP;
  E_Int* dataA = (E_Int*)PyArray_DATA(crA);
  E_Int* dataP = (E_Int*)PyArray_DATA(crP);

#pragma omp parallel for default(shared)
  for (E_Int p = 0; p < ne; p++)
    for (E_Int k = 0; k < stype; k++)
      dataP[k+stype*p] = dataA[p+ne*k];

  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
// Copy le champ cB dans cA en position nf
//=============================================================================
PyObject* K_CONVERTER::cpyValueByField(PyObject* self, PyObject* args)
{
  IMPORTNUMPY;
  PyObject* cA; PyObject* cB;
  E_Int np=0, nf=0;
  if (!PYPARSETUPLE_(args, OO_ II_, &cA, &cB, &np, &nf)) return NULL;

  PyArrayObject* crA = (PyArrayObject*)cA;
  PyArrayObject* crB = (PyArrayObject*)cB;
  E_Float* dataA = (E_Float*)PyArray_DATA(crA);
  E_Float* dataB = (E_Float*)PyArray_DATA(crB);
  E_Float* pA = dataA+nf*np;

#pragma omp parallel for default(shared)
  for (E_Int p = 0; p < np; p++) pA[p] = dataB[p];

  Py_INCREF(Py_None);
  return Py_None;
}
