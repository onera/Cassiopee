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
# include "converter.h"

using namespace K_FLD;

void makeRange(E_Int NJ, std::vector< std::vector<E_Int> >& line, std::vector< std::vector<E_Int> >& range)
{
  //printf("make range %d\n", NJ);

  for (E_Int j = 0; j < NJ; j++)
  {
    std::vector<E_Int>& v = line[j];
    size_t s = v.size();
    //printf("make range in %d\n", s);
    E_Int imin = 1.e6; E_Int imax = -1.e6;
    for (size_t n = 0; n < s; n++)
    {
      imin = K_FUNC::E_min(v[n], imin);
      imax = K_FUNC::E_max(v[n], imax);
    }
    if (s > 0)
    {
      range[j].push_back(imin);
      range[j].push_back(imax);
    }
  }
}

void accumulate(E_Int ni, E_Int nj, E_Int nk,
                E_Int NJ, E_Int KP, std::vector< std::vector<E_Int> >& range, E_Int type, std::vector<E_Int>& ranges)
{
  E_Int iminL, imaxL, jminL, jmaxL;
  E_Int imin, imax;
  iminL = -1; imaxL = -1;
  jminL = -1; jmaxL = -1;
  for (E_Int j = 0; j < NJ; j++)
  {
    //printf("size %d\n", range[j].size());
    if (range[j].size() > 1)
    {
      imin = range[j][0]; imax = range[j][1];
      iminL = imin; imaxL = imax;
      jminL = j; jmaxL = j;
      //printf("type=%d, %d %d %d %d %d\n", type, iminL, jminL, imaxL, jmaxL, KP);
   
      if (type == 1)
      {
        ranges.push_back(KP+1);
        ranges.push_back(KP+1);
        ranges.push_back(iminL+1);
        ranges.push_back(K_FUNC::E_min(imaxL+2,nj));
        ranges.push_back(jmaxL+1);
        ranges.push_back(K_FUNC::E_min(jmaxL+2,nk));                
      }
      else if (type == 2)
      {
        ranges.push_back(iminL+1);
        ranges.push_back(K_FUNC::E_min(imaxL+2,ni));
        ranges.push_back(KP+1);
        ranges.push_back(KP+1);
        ranges.push_back(jminL+1);
        ranges.push_back(K_FUNC::E_min(jmaxL+2,nk));                
      }
      else
      {
        ranges.push_back(iminL+1);
        ranges.push_back(K_FUNC::E_min(imaxL+2,ni));
        ranges.push_back(jminL+1);
        ranges.push_back(K_FUNC::E_min(jmaxL+2,nj));
        ranges.push_back(KP+1);
        ranges.push_back(KP+1);
      }
    }
  }
}

void addToList(PyObject* tpl, std::vector<E_Int>& ranges)
{
  E_Int nr = ranges.size()/6;
  for (E_Int i = 0; i < 6*nr; i += 6)
  {
    PyObject* s = Py_BuildValue("[llllll]", ranges[i],
                                ranges[i+1],ranges[i+2],
                                ranges[i+3],
                                ranges[i+4],ranges[i+5]);
    PyList_Append(tpl, s); Py_DECREF(s);
  }
}

//=============================================================================
/* Convert a PointList to a list of ranges for structured grids */
//=============================================================================
PyObject* K_CONVERTER::pointList2Ranges(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int ni, nj, nk;
  if (!PYPARSETUPLEI(args, "Olll", "Oiii", &array, &ni, &nj, &nk)) return NULL;

  // Check numpy (pointlist)
  FldArrayI* PL;
  E_Int res = K_NUMPY::getFromNumpyArray(array, PL, true);

  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "pointList2Ranges: numpy is invalid.");
    return NULL;
  }
  E_Int nf = PL->getSize();
  E_Int* p = PL->begin();
  E_Int ni1 = K_FUNC::E_max(ni-1,1);
  E_Int nj1 = K_FUNC::E_max(nj-1,1);
  E_Int nk1 = K_FUNC::E_max(nk-1,1);
  
  E_Int ninti = ni*nj1*nk1;
  E_Int nintj = ni1*nj*nk1;
  //E_Int nintk = ni1*nj1*nk;
  E_Int i,j,k,ind;

  // recherche du type d'interface. Compte.
  E_Int nimin = 0; // interfaces en imin
  E_Int nimax = 0;
  E_Int njmin = 0;
  E_Int njmax = 0;
  E_Int nkmin = 0;
  E_Int nkmax = 0;

  E_Int KP;
  
  for (E_Int n = 0; n < nf; n++)
  {
    ind = p[n];
    if (ind >= ninti+nintj) // interface k
    { 
      KP = (ind-ninti-nintj)/(ni1*nj1);
      if (KP == 0) nkmin++; else nkmax++;
    }
    else if (ind >= ninti) // interface j 
    { 
      k = (ind-ninti)/(ni1*nj);
      KP = (ind-ninti-k*ni1*nj)/ni1;
      if (KP == 0) njmin++; else njmax++; 
    }
    else // interface i
    { 
      k = ind/(ni*nj1);
      j = (ind-k*ni*nj1)/ni;
      KP = ind-j*ni-k*ni*nj1;
      if (KP == 0) nimin++; else nimax++;
    }
  }

  //printf("ni,nj=%d %d - %d\n",ni,nj,ind);
  
  std::vector< std::vector<E_Int> > line1; // imin 
  if (nimin > 0) line1.resize(nk);
  std::vector< std::vector<E_Int> > line2; // imax
  if (nimax > 0) line2.resize(nk);
  std::vector< std::vector<E_Int> > line3; // jmin
  if (njmin > 0) line3.resize(nk);
  std::vector< std::vector<E_Int> > line4; // jmax
  if (njmax > 0) line4.resize(nk);
  std::vector< std::vector<E_Int> > line5; // kmin 
  if (nkmin > 0) line5.resize(nj);
  std::vector< std::vector<E_Int> > line6; // kmax 
  if (nkmax > 0) line6.resize(nj);
  
  //printf("nimin=%d nimax=%d\n", nimin, nimax);
  //printf("njmin=%d njmax=%d\n", njmin, njmax);
  //printf("nkmin=%d nkmax=%d\n", nkmin, nkmax);

  // Stockage i,j de l'interface
  for (E_Int n = 0; n < nf; n++)
  {
    ind = p[n];
    if (ind >= ninti+nintj) // interface k
    {
      // ind = i+j*ni1+k*ni1*nj1;
      ind = ind-ninti-nintj;
      k = ind/(ni1*nj1);
      j = (ind-k*ni1*nj1)/ni1;
      i = ind-j*ni1-k*ni1*nj1;
      if (k == 0) line5[j].push_back(i);
      else line6[j].push_back(i);
    }
    else if (ind >= ninti) // interface j
    {
      // ind = i+j*ni1+k*ni1*nj;
      ind = ind-ninti;
      k = ind/(ni1*nj);
      j = (ind-k*ni1*nj)/ni1;
      i = ind-j*ni1-k*ni1*nj;
      if (j == 0) line3[k].push_back(i);
      else line4[k].push_back(i);
    }
    else // interface i
    {
      // ind = i+j*ni+k*ni*nj1;
      k = ind/(ni*nj1);
      j = (ind-k*ni*nj1)/ni;
      i = ind-j*ni-k*ni*nj1;
      //printf("ijk %d %d %d\n",i,j,k);
      if (i == 0) line1[k].push_back(j);
      else line2[k].push_back(j);
    }
  }
  //printf("formation des ranges\n");

  // Formation des ranges
  std::vector< std::vector<E_Int> > range1; // imin 
  if (nimin > 0) range1.resize(nk);
  std::vector< std::vector<E_Int> > range2; // imax
  if (nimax > 0) range2.resize(nk);
  std::vector< std::vector<E_Int> > range3; // jmin
  if (njmin > 0) range3.resize(nk);
  std::vector< std::vector<E_Int> > range4; // jmax
  if (njmax > 0) range4.resize(nk);
  std::vector< std::vector<E_Int> > range5; // kmin 
  if (nkmin > 0) range5.resize(nj);
  std::vector< std::vector<E_Int> > range6; // kmax 
  if (nkmax > 0) range6.resize(nj);

  if (nimin > 0) makeRange(nk, line1, range1);
  if (nimax > 0) makeRange(nk, line2, range2); 
  if (njmin > 0) makeRange(nk, line3, range3);
  if (njmax > 0) makeRange(nk, line4, range4);
  if (nkmin > 0) makeRange(nj, line5, range5); 
  if (nkmax > 0) makeRange(nj, line6, range6); 
  
  //printf("sustained\n");
  //for (E_Int j = 0; j < NJ; j++) printf("size=%d\n", range[j].size());

  // Fusion des ranges
  std::vector<E_Int> ranges1;
  std::vector<E_Int> ranges2;
  std::vector<E_Int> ranges3;
  std::vector<E_Int> ranges4;
  std::vector<E_Int> ranges5;
  std::vector<E_Int> ranges6;
  
  if (nimin > 0) accumulate(ni,nj,nk, nk, 0, range1, 1, ranges1);
  if (nimax > 0) accumulate(ni,nj,nk, nk, ni-1, range2, 1, ranges2); 
  if (njmin > 0) accumulate(ni,nj,nk, nk, 0, range3, 2, ranges3); 
  if (njmax > 0) accumulate(ni,nj,nk, nk, nj-1, range4, 2, ranges4);
  if (nkmin > 0) accumulate(ni,nj,nk, nj, 0, range5, 3, ranges5); 
  if (nkmax > 0) accumulate(ni,nj,nk, nj, nk-1, range6, 3, ranges6); 

  //printf("prod\n");
  // Construction des ranges
  
  //printf("nr=%d\n", nr);
  PyObject* tpl = PyList_New(0);

  addToList(tpl, ranges1);
  addToList(tpl, ranges2);
  addToList(tpl, ranges3);
  addToList(tpl, ranges4);
  addToList(tpl, ranges5);
  addToList(tpl, ranges6);

  RELEASESHAREDN(array, PL);

  return tpl;
}
