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
// detectEmptyBC: detect empty BC of a 2D or 3D-array

# include "converter.h"
using namespace std;
using namespace K_FLD;

//=============================================================================
/* Detect empty BC of a 2 or 3D-array */
//=============================================================================
PyObject* K_CONVERTER::detectEmptyBC(PyObject* self, PyObject* args)
{
  E_Int Ni=0, Nj=0, Nk=0, D;
  PyObject* winList;
  PyObject* NwinList;
  if (!PYPARSETUPLE_(args, O_ IIII_ O_, &winList, &Ni, &Nj, &Nk, &D, &NwinList)) return NULL;

  E_Int ni  = Ni; E_Int nj  = Nj; E_Int nk  = Nk;
  E_Int dir = D;

  FldArrayI* tag = NULL;
  vector<E_Int*> win;
  vector<E_Int*> nwins;
  PyObject* tpl;
  PyObject* range;
  E_Int rangeSize = 0;

  // check if winList is a list
  if (PyList_Check(winList) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "detectEmptyBC: win argument must be a list.");
    return NULL;
  }
  // check if NwinList is a list
  if (PyList_Check(NwinList) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "detectEmptyBC: nwins argument must be a list.");
    return NULL;
  }
  E_Int winListSize = PyList_Size(winList);
  E_Int NwinListSize = PyList_Size(NwinList);

  // Convert Python winList in vector<E_Int*> win
  for (E_Int w = 0; w < winListSize; w++)
  {
    range = PyList_GetItem(winList, w);
    if (PyList_Check (range) == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "detectEmptyBC: range argument must be a list.");
      return NULL;
    }
    rangeSize = PyList_Size(range);
    E_Int* r = new E_Int[rangeSize]; win.push_back(r);
    
    for (E_Int t = 0; t < rangeSize; t++)
    {
      tpl = PyList_GetItem(range, t);  
      win[w][t] = PyInt_AsLong(tpl);
      if (win[w][t] == -1)
      {
        PyErr_SetString(PyExc_TypeError, 
                        "detectEmptyBC: win index must be an integer.");
        for (E_Int w2 = 0; w2 < winListSize; w2++) delete [] win[w2]; 
        return NULL;
      }
    }
  }

  // Convert Python NwinList in vector<E_Int*> nwins
  for (E_Int w = 0; w < NwinListSize; w++)
  {
    range = PyList_GetItem(NwinList, w);
    if (PyList_Check (range) == 0)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "detectEmptyBC: range argument must be a list.");
      for (E_Int w2 = 0; w2 < winListSize; w2++) delete [] win[w2];
      return NULL;
    }
    rangeSize = PyList_Size(range);
    E_Int* r = new E_Int[rangeSize]; nwins.push_back(r);
    for (E_Int t = 0; t < rangeSize; t++)
    {
      tpl = PyList_GetItem(range, t);
      nwins[w][t] = PyInt_AsLong(tpl);
      if (nwins[w][t] == -1)
      {
        PyErr_SetString(PyExc_TypeError, 
                        "detectEmptyBC: nwins index must be an integer.");
        for (E_Int w2 = 0; w2 < winListSize; w2++) delete [] win[w2];
        for (E_Int w2 = 0; w2 < NwinListSize; w2++) delete [] nwins[w2];    
        return NULL;
      }
    }
  }

  rangeSize = 6;
  detectEmptyBCrec(win, ni, nj, nk, dir, nwins, tag, rangeSize);

  if (tag != NULL) delete tag;
  E_Int size = nwins.size();

  // Convert vector<E_Int*> nwins in PyObject tplnwins
  PyObject * tplnwins;
  PyObject* l = PyList_New(0);
  for (E_Int i = 0; i < size; i++)
  {
#ifdef E_DOUBLEINT
    tplnwins = Py_BuildValue("[l,l,l,l,l,l]",
                             long(nwins[i][0]),long(nwins[i][1]),long(nwins[i][2]),
                             long(nwins[i][3]),long(nwins[i][4]),long(nwins[i][5]));
#else
    tplnwins = Py_BuildValue("[i,i,i,i,i,i]",
                             nwins[i][0],nwins[i][1],nwins[i][2],
                             nwins[i][3],nwins[i][4],nwins[i][5]);
#endif
    PyList_Append(l, tplnwins);
    Py_DECREF(tplnwins);
  }

  // Deleting win and nwins
  for (E_Int i = 0; i < winListSize; i++) delete [] win[i];
  for (E_Int i = 0; i < size; i++) delete [] nwins[i];
  win.clear(); nwins.clear();

  return l;
}

//=============================================================================
/* Detect empty BC of a 3D-array: recursive C method */
//=============================================================================
void K_CONVERTER::detectEmptyBCrec(
  vector<E_Int*>& win, 
  E_Int ni, E_Int nj, E_Int nk, E_Int dir, 
  vector<E_Int*>& nwins, FldArrayI*& tag, E_Int size)
{
  // variables internes a la methode
  E_Int ni1,nj1,nk1;
  E_Int i1,i2,j1,j2,k1,k2;
  E_Int found;
  E_Int winsize;
  E_Int istart = 0; E_Int jstart = 0; E_Int kstart = 0;
  E_Int iend = 0; E_Int jend = 0; E_Int kend = 0;

  ni1 = ni-1; nj1 = nj-1; nk1 = nk-1;

  winsize = win.size();
  //nwinssize = nwins.size();

  // Direction = 1 (face i=1)
  // ------------------------
  if (dir == 1)
  {
    // init du tag
    if (tag == NULL)
    {
      tag = new FldArrayI(nj*nk);
      tag->setAllValuesAtNull();
      E_Int* tagp = tag->begin();
      for (E_Int w = 0; w < winsize; w++)
      {
        E_Int* winp = win[w];
        i1 = winp[0]-1; j1 = winp[2]-1; k1 = winp[4]-1;
        i2 = winp[1]-1; j2 = winp[3]-1; k2 = winp[5]-1;
        if ((i1 == 0) &&(i2 == 0))
        {
          for (E_Int k = k1; k < k2; k++)
            for (E_Int j = j1; j < j2; j++)
              tagp[j+k*nj] = 1;
          if (j2 == nj-1)
            for (E_Int k = k1; k < k2; k++)            
              tagp[j2+k*nj] = 1;
          if (k2 == nk-1)
            for (E_Int j = j1; j < j2; j++)            
              tagp[j+k2*nj] = 1;
        }
      }
    }

    // recherche du premier pt pas ds 1 BC
    E_Int* tagp = tag->begin();
    found = 0;
    for (E_Int k = 0; k < nk; k++)
    {
      for (E_Int j = 0; j < nj; j++)
      {
        if (tagp[j+k*nj] == 0)
        { 
          jstart = j; kstart = k; found = 1; break;
        }
      }
      if (found == 1) break;
    }
    // detection de ts les pts suivants : arret qd un pt ds une bc
    jend = jstart; kend = kstart;
    for (E_Int j = jstart; j < nj; j++)
    {
      if (tagp[j+kstart*nj] == 1)
      {
        jend = j; break;
      }
      else jend = j;
    } 
    found = 0;
    for (E_Int k = kstart; k < nk; k++)
    {
      for (E_Int j = jstart; j < jend; j++)
      {
        if (tagp[j+k*nj] == 1)
        {
          kend = k; found = 1; break;
        }
        else kend = k;
      }
      if (found == 1) break;
    }

    // mise a jour de tag
    for (E_Int k = kstart; k < kend; k++)
      for (E_Int j = jstart; j < jend; j++)
        tagp[j+k*nj] = 1;
    if (jend == nj-1)
      for (E_Int k = kstart; k < kend; k++)
        tagp[jend+k*nj] = 1;
    if (kend == nk-1)
      for (E_Int j = jstart; j < jend; j++)
        tagp[j+kend*nj] = 1;
    // restart
    if ((jstart > -1 && jend > jstart && kstart > -1 && kend > kstart) || 
        (jstart == 0 && jend == 0 && nj == 1  && kstart > -1 && kend > kstart) || 
        (jstart > -1 && jend > jstart && kstart == 0 && kend == 0 && nk == 1) )
    {
      E_Int* newrange = new E_Int[size];
      newrange[0] = 1;
      newrange[1] = 1;
      newrange[2] = jstart+1;
      newrange[3] = jend+1;
      newrange[4] = kstart+1;
      newrange[5] = kend+1;
      nwins.push_back(newrange);
      detectEmptyBCrec(win, ni, nj, nk, dir, nwins, tag, size);
    }
  
    else return;
  }
  // Direction = 2 (face i=ni)
  // ---------------
  else if (dir == 2)
  {
    // init du tag
    if (tag == NULL)
    {
      tag = new FldArrayI(nj*nk);
      tag->setAllValuesAtNull();
      E_Int* tagp = tag->begin();
      for (E_Int w = 0; w < winsize; w++)
      {
        E_Int* winp = win[w];
        i1 = winp[0]-1; j1 = winp[2]-1; k1 = winp[4]-1;
        i2 = winp[1]-1; j2 = winp[3]-1; k2 = winp[5]-1;
        if ((i1 == ni1) &&(i2 == ni1))
        {
          for (E_Int k = k1; k < k2; k++)
            for (E_Int j = j1; j < j2; j++)
              tagp[j+k*nj] = 1;
          if (j2 == nj-1)
            for (E_Int k = k1; k < k2; k++)            
              tagp[j2+k*nj] = 1;
          if (k2 == nk-1)
            for (E_Int j = j1; j < j2; j++)            
              tagp[j+k2*nj] = 1;
        }
      }
    }
    // recherche du premier pt pas ds 1 BC
    E_Int* tagp = tag->begin();
    found = 0;
    for (E_Int k = 0; k < nk; k++)
    {
      for (E_Int j = 0; j < nj; j++)
      {
        if ( tagp[j+k*nj] == 0 )
        { 
          jstart = j; kstart = k; found = 1; break;
        }
      }
      if (found == 1) break;
    }
    // detection de ts les pts suivants : arret qd un pt ds une bc 
    jend = jstart; kend = kstart;
    for (E_Int j = jstart; j < nj; j++)
    {
      if (tagp[j+kstart*nj] == 1)
      {
        jend = j; break;
      }
      else jend = j;
    }
    
    found = 0;
    for (E_Int k = kstart; k < nk; k++)
    {
      for (E_Int j = jstart; j < jend; j++)
      {
        if (tagp[j+k*nj] == 1)
        {
          kend = k; found = 1; break;
        }
        else kend = k;
      }
      if (found == 1) break;
    }
    // mise a jour de tag
    for (E_Int k = kstart; k < kend; k++)
      for (E_Int j = jstart; j < jend; j++)
        tagp[j+k*nj] = 1;
    if (jend == nj-1)
      for (E_Int k = kstart; k < kend; k++)
        tagp[jend+k*nj] = 1;
    if (kend == nk-1)
      for (E_Int j = jstart; j < jend; j++)
        tagp[j+kend*nj] = 1;
    // restart
    if ((jstart > -1 && jend > jstart && kstart > -1 && kend > kstart) || 
        (jstart == 0 && jend == 0 && nj == 1  && kstart > -1 && kend > kstart) || 
        (jstart > -1 && jend > jstart && kstart == 0 && kend == 0 && nk == 1) )
    {
      E_Int* newrange = new E_Int[size];
      newrange[0] = ni;
      newrange[1] = ni;
      newrange[2] = jstart+1;
      newrange[3] = jend+1;
      newrange[4] = kstart+1;
      newrange[5] = kend+1;
      nwins.push_back(newrange);
      detectEmptyBCrec(win, ni, nj, nk, dir, nwins, tag, size);
    }
    else return;
  }
  // Direction = 3 (face j=1)
  // ---------------
  else if (dir == 3)
  {
    // init du tag
    if (tag == NULL)
    {
      tag = new FldArrayI(ni*nk);
      tag->setAllValuesAtNull();
      E_Int* tagp = tag->begin();
      for (E_Int w = 0; w < winsize; w++)
      {
        E_Int* winp = win[w];
        i1 = winp[0]-1; j1 = winp[2]-1; k1 = winp[4]-1;
        i2 = winp[1]-1; j2 = winp[3]-1; k2 = winp[5]-1;
        if ((j1 == 0) && (j2 == 0))
        {
          for (E_Int k = k1; k < k2; k++)
            for (E_Int i = i1; i < i2; i++)
              tagp[i+k*ni] = 1;
          if (i2 == ni-1)
            for (E_Int k = k1; k < k2; k++)            
              tagp[i2+k*ni] = 1;
          if (k2 == nk-1)
            for (E_Int i = i1; i < i2; i++)            
              tagp[i+k2*ni] = 1;
        }
      }
    }
    // recherche du premier pt pas ds 1 BC
    E_Int* tagp = tag->begin();
    found = 0;
    for (E_Int k = 0; k < nk; k++)
    {
      for (E_Int i = 0; i < ni; i++)
      {
        if (tagp[i+k*ni] == 0)
        { 
          istart = i; kstart = k; found = 1; break;
        }
      }
      if (found == 1) break;
    }
    // detection de ts les pts suivants : arret qd un pt ds une bc 
    iend = istart; kend = kstart;
    for (E_Int i = istart; i < ni; i++)
    {
      if (tagp[i+kstart*ni] == 1)
      {
        iend = i; break;
      }
      else iend = i;
    }
    found = 0;
    for (E_Int k = kstart; k < nk; k++)
    {
      for (E_Int i = istart; i < iend; i++)
      {
        if (tagp[i+k*ni] == 1)
        {
          kend = k; found = 1; break;
        }
        else kend = k;
      }
      if (found == 1) break;
    }

    // mise a jour de tag
    for (E_Int k = kstart; k < kend; k++)
      for (E_Int i = istart; i < iend; i++)
        tagp[i+k*ni] = 1;
    if ( iend == ni-1 )
      for (E_Int k = kstart; k < kend; k++)
        tagp[iend+k*ni] = 1;
    if ( kend == nk-1 )
      for (E_Int i = istart; i < iend; i++)
        tagp[i+kend*ni] = 1;
    // restart
    if ((istart > -1 && iend > istart && kstart > -1 && kend > kstart) || 
        (istart == 0 && iend == 0 && ni == 1  && kstart > -1 && kend > kstart) || 
        (istart > -1 && iend > istart && kstart == 0 && kend == 0 && nk == 1) )
    {
      E_Int* newrange = new E_Int[size];
      newrange[0] = istart+1;
      newrange[1] = iend+1;
      newrange[2] = 1;
      newrange[3] = 1;
      newrange[4] = kstart+1;
      newrange[5] = kend+1;
      nwins.push_back(newrange);
      detectEmptyBCrec(win, ni, nj, nk, dir, nwins, tag, size);
    }
    else return;
  }
// Direction = 4 (face j=nj)
// ---------------
  else if (dir == 4)
  {
    // init du tag
    if (tag == NULL)
    {
      tag = new FldArrayI(ni*nk);
      tag->setAllValuesAtNull();
      E_Int* tagp = tag->begin();
      for (E_Int w = 0; w < winsize; w++)
      {
        E_Int* winp = win[w];
        i1 = winp[0]-1; j1 = winp[2]-1; k1 = winp[4]-1;
        i2 = winp[1]-1; j2 = winp[3]-1; k2 = winp[5]-1;
        if ((j1 == nj1) &&(j2 == nj1))
        {
          for (E_Int k = k1; k < k2; k++)
            for (E_Int i = i1; i < i2; i++)
              tagp[i+k*ni] = 1;
          if (i2 == ni-1)
            for (E_Int k = k1; k < k2; k++)            
              tagp[i2+k*ni] = 1;
          if (k2 == nk-1)
            for (E_Int i = i1; i < i2; i++)            
              tagp[i+k2*ni] = 1;
        }
      }
    }
    // recherche du premier pt pas ds 1 BC
    found = 0;
    E_Int* tagp = tag->begin();
    for (E_Int k = 0; k < nk; k++)
    {
      for (E_Int i = 0; i < ni; i++)
      {
        if (tagp[i+k*ni] == 0)
        { 
          istart = i; kstart = k;
          found = 1;
          break;
        }
      }
      if (found == 1) break;
    }
    // detection de ts les pts suivants : arret qd un pt ds une bc 
    iend = istart; kend = kstart;
    for (E_Int i = istart; i < ni; i++)
    {
      if (tagp[i+kstart*ni] == 1)
      {
        iend = i; break;
      }
      else
        iend = i;
    }
    found = 0;
    for (E_Int k = kstart; k < nk; k++)
    {
      for (E_Int i = istart; i < iend; i++)
      {
        if (tagp[i+k*ni] == 1)
        {
          kend = k; found = 1; break;
        }
        else
          kend = k;
      }
      if (found == 1) break;
    }
    // mise a jour de tag
    for (E_Int k = kstart; k < kend; k++)
      for (E_Int i = istart; i < iend; i++)
        tagp[i+k*ni] = 1;
    if (iend == ni-1)
      for (E_Int k = kstart; k < kend; k++)
        tagp[iend+k*ni] = 1;
    if (kend == nk-1)
      for (E_Int i = istart; i < iend; i++)
        tagp[i+kend*ni] = 1;
    // restart
    if ((istart > -1 && iend > istart && kstart > -1 && kend > kstart) || 
        (istart == 0 && iend == 0 && ni == 1  && kstart > -1 && kend > kstart) || 
        (istart > -1 && iend > istart && kstart == 0 && kend == 0 && nk == 1) )    
    {
      E_Int* newrange = new E_Int[size];
      newrange[0] = istart+1;
      newrange[1] = iend+1;
      newrange[2] = nj;
      newrange[3] = nj;
      newrange[4] = kstart+1;
      newrange[5] = kend+1;
      nwins.push_back(newrange);
      detectEmptyBCrec(win, ni, nj, nk, dir, nwins, tag, size);
    }
    else return;
  }
// Direction = 5 (face k=1)
// ---------------
  else if (dir == 5)
  {
    // init du tag
    if (tag == NULL)
    {
      tag = new FldArrayI(ni*nj);
      tag->setAllValuesAtNull();
      E_Int* tagp = tag->begin();
      for (E_Int w = 0; w < winsize; w++)
      {
        E_Int* winp = win[w];
        i1 = winp[0]-1; j1 = winp[2]-1; k1 = winp[4]-1;
        i2 = winp[1]-1; j2 = winp[3]-1; k2 = winp[5]-1;
        if ((k1 == 0) &&(k2 == 0))
        {
          for (E_Int j = j1; j < j2; j++)
            for (E_Int i = i1; i < i2; i++)
              tagp[i+j*ni] = 1;
          if (i2 == ni-1)
            for (E_Int j = j1; j < j2; j++)            
              tagp[i2+j*ni] = 1;
          if (j2 == nj-1)
            for (E_Int i = i1; i < i2; i++)            
              tagp[i+j2*ni] = 1;
        }
      }
    }
    // recherche du premier pt pas ds 1 BC
    E_Int* tagp = tag->begin();
    found = 0;
    for (E_Int j = 0; j < nj; j++)
    {
      for (E_Int i = 0; i < ni; i++)
      {
        if (tagp[i+j*ni] == 0)
        { 
          istart = i; jstart = j;
          found = 1;
          break;
        }
      }
      if (found == 1) break;
    }
    // detection de ts les pts suivants : arret qd un pt ds une bc 
    iend = istart; jend = jstart;
    for (E_Int i = istart; i < ni; i++)
    {
      if (tagp[i+jstart*ni] == 1)
      {
        iend = i; break;
      }
      else
        iend = i;
    }
    
    found = 0;
    for (E_Int j = jstart; j < nj; j++)
    {
      for (E_Int i = istart; i < iend; i++)
      {
        if (tagp[i+j*ni] == 1)
        {
          jend = j; found = 1; break;
        }
        else jend = j;
      }
      if (found == 1) break;
    }
    // mise a jour de tag
    for (E_Int j = jstart; j < jend; j++)
      for (E_Int i = istart; i < iend; i++)
        tagp[i+j*ni] = 1;
    if ( iend == ni-1 )
      for (E_Int j = jstart; j < jend; j++)
        tagp[iend+j*ni] = 1;
    if ( jend == nj-1 )
      for (E_Int i = istart; i < iend; i++)
        tagp[i+jend*ni] = 1;
    // restart
    if ((istart > -1 && iend > istart && jstart > -1 && jend > jstart ) || 
        (istart == 0 && iend == 0 && ni == 1  && jstart > -1 && jend > jstart) || 
        (istart > -1 && iend > istart && jstart == 0 && jend == 0 && nj == 1) )    
    {
      E_Int* newrange = new E_Int[size];
      newrange[0] = istart+1;
      newrange[1] = iend+1;
      newrange[2] = jstart+1;
      newrange[3] = jend+1;
      newrange[4] = 1;
      newrange[5] = 1;
      nwins.push_back(newrange);
      detectEmptyBCrec(win, ni, nj, nk, dir, nwins, tag, size);
    }
    else return;
  }
// Direction = 6 (face k=nk)
// ---------------
  else
  {
    // init du tag
    if (tag == NULL)
    {
      tag = new FldArrayI(ni*nj);
      tag->setAllValuesAtNull();
      E_Int* tagp = tag->begin();
      for (E_Int w = 0; w < winsize; w++)
      {
        E_Int* winp = win[w];
        i1 = winp[0]-1; j1 = winp[2]-1; k1 = winp[4]-1;
        i2 = winp[1]-1; j2 = winp[3]-1; k2 = winp[5]-1;

        if ((k1 == nk1) &&(k2 == nk1))
        {
          for (E_Int j = j1; j < j2; j++)
            for (E_Int i = i1; i < i2; i++)
              tagp[i+j*ni] = 1;
          if (i2 == ni-1)
            for (E_Int j = j1; j < j2; j++)            
              tagp[i2+j*ni] = 1;
          if (j2 == nj-1)
            for (E_Int i = i1; i < i2; i++)            
              tagp[i+j2*ni] = 1;
        }
      }
    }
    // recherche du premier pt pas ds 1 BC
    E_Int* tagp = tag->begin();
    found = 0;
    for (E_Int j = 0; j < nj; j++)
    {
      for (E_Int i = 0; i < ni; i++)
      {
        if ( tagp[i+j*ni] == 0 )
        { 
          istart = i; jstart = j;
          found = 1;
          break;
        }
      }
      if (found == 1) break;
    }
    // detection de ts les pts suivants : arret qd un pt ds une bc 
    iend = istart; jend = jstart;
    for (E_Int i = istart; i < ni; i++)
    {
      if (tagp[i+jstart*ni] == 1)
      {
        iend = i; break;        
      }
      else iend = i;
    }
    found = 0;
    for (E_Int j = jstart; j < nj; j++)
    {
      for (E_Int i = istart; i < iend; i++)
      {
        if (tagp[i+j*ni] == 1)
        {
          jend = j; found = 1; break;
        }
        else jend = j;
      }
      if (found == 1) break;
    }
    // mise a jour de tag
    for (E_Int j = jstart; j < jend; j++)
      for (E_Int i = istart; i < iend; i++)
        tagp[i+j*ni] = 1;
    if (iend == ni-1)
      for (E_Int j = jstart; j < jend; j++)
        tagp[iend+j*ni] = 1;
    if (jend == nj-1)
      for (E_Int i = istart; i < iend; i++)
        tagp[i+jend*ni] = 1;
    // restart
    if ((istart > -1 && iend > istart && jstart > -1 && jend > jstart) || 
        (istart == 0 && iend == 0 && ni == 1  && jstart > -1 && jend > jstart) || 
        (istart > -1 && iend > istart && jstart == 0 && jend == 0 && nj == 1) )        
    {
      E_Int* newrange = new E_Int[6];
      newrange[0] = istart+1;
      newrange[1] = iend+1;
      newrange[2] = jstart+1;
      newrange[3] = jend+1;
      newrange[4] = nk;
      newrange[5] = nk;
      nwins.push_back(newrange);
      detectEmptyBCrec(win, ni, nj, nk, dir, nwins, tag, size);
    }
    else return;
  }
}
