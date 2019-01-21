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

// Diff 2 arrays defining the same solution

#include "converter.h"
using namespace std;

using namespace K_FUNC;
using namespace K_FLD;

#ifdef __MACH__
// For Mac OSX
#include <math.h>
//#define isnan __isnand
//#define isinf __isinf

#elif defined(_WIN32)
// rien a faire maintenant

#elif defined(_WIN64)
// For Windows
#include <float.h>
#ifndef isnan
#define isnan _isnan
#endif
#ifndef isinf
#define isinf(x) (!_finite(x))
#endif

#else
#include <cmath>
#define isnan std::isnan
#define isinf std::isinf
#endif

//=============================================================================
/* Diffing arrays: 
   if 2 input arrays: they define the same mesh + diff their solution 
   if 3 input arrays: 
     - 1st arrays: fine reference solution, 
     - 2nd arrays: ref solution in centers,
     - 3rd arrays: solution to be tested.
*/
//=============================================================================
PyObject* K_CONVERTER::diffArrays(PyObject* self, PyObject* args)
{
  PyObject* arrays1; PyObject* arrays2; PyObject* arrays3;

  E_Int narrays = 3;
  if (!PyArg_ParseTuple(args, "OOO", &arrays1, &arrays2, &arrays3))
    return NULL;

  // Check every arrays
  if (PyList_Check(arrays1) == 0 || PyList_Check(arrays2) == 0)
  {
    PyErr_SetString(
      PyExc_TypeError, 
      "diffArrays: first/second argument must be lists.");
    return NULL;
  }

  if (narrays == 3 && PyList_Check(arrays3) == 0)
  {
    PyErr_SetString(
      PyExc_TypeError, 
      "diffArrays: third argument must be a list.");
    return NULL;
  }

  if (PyList_Size(arrays3) == 0) narrays = 2; 

  if (narrays == 2) // compares solution on the same grid
    return diff2(arrays1, arrays2);

  else if (narrays == 3) // compare solution to the reference
    return diff3(arrays1, arrays2, arrays3);

  else 
  {
    PyErr_SetString(
      PyExc_TypeError, 
      "diffArrays: invalid number of arrays in the list.");
    return NULL;
  }
}

//=============================================================================
/* Diff2: for 2 arrays defining the same mesh */
//=============================================================================
PyObject* K_CONVERTER::diff2(PyObject* arrays1, PyObject* arrays2)
{
  PyObject* tpl;
  E_Int res;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk;
  char* varString1; char* eltType1;
  char* varString2; char* eltType2;

  // Extract info from arrays1
  vector<E_Int> ni1; vector<E_Int> nj1; vector<E_Int> nk1;
  vector<FldArrayF*> field1;
  vector<FldArrayI*> cn1;
  vector<char*> elt1;
  vector<PyObject*> object1;

  E_Int n1 = PyList_Size(arrays1);
  E_Int nfld1 = -1; E_Int nfld2 = -1;
  
  for (int i = 0; i < n1; i++)
  {
    tpl = PyList_GetItem(arrays1, i);
    res = K_ARRAY::getFromArray2(tpl, varString1, f, ni, nj, nk, cn, eltType1);
    object1.push_back(tpl);
    if (res == 1)
    {
      if (ni*nj*nk > 0)
      {
        ni1.push_back(ni); nj1.push_back(nj); nk1.push_back(nk);
        field1.push_back(f);
        cn1.push_back(NULL);
        elt1.push_back(NULL);
        if (nfld1 == -1) nfld1 = f->getNfld();
        else
        {
          if (f->getNfld() != nfld1)
          {
            PyErr_SetString(PyExc_ValueError,
                            "diffArrays: number of fields must be equal in 1st list.");
            return NULL;
          }
        }
      }
      else printf("Warning: diffArrays: one array is empty.\n");
    }
    else if (res == 2) 
    {
      ni1.push_back(-1); nj1.push_back(-1); nk1.push_back(-1);
      field1.push_back(f);
      cn1.push_back(cn);
      elt1.push_back(eltType1);
    }
    else
    {
      printf("Warning: diffArrays: not an array. Array skipped...\n");
    }
  }

  // Extract info from arrays2
  vector<E_Int> ni2; vector<E_Int> nj2; vector<E_Int> nk2;
  vector<FldArrayF*> field2;
  vector<PyObject*> object2;
  vector<FldArrayI*> cn2;
  E_Int n2 = PyList_Size(arrays2);

  for (int i = 0; i < n2; i++)
  {
    tpl = PyList_GetItem(arrays2, i);
    res = K_ARRAY::getFromArray2(tpl, varString2, f, ni, nj, nk, cn, eltType2);
    object2.push_back(tpl);

    if (res == 1)
    {
      if (ni*nj*nk > 0)
      {
        ni2.push_back(ni); nj2.push_back(nj); nk2.push_back(nk);
        field2.push_back(f);
        cn2.push_back(NULL);
        if (nfld2 == -1) nfld2 = f->getNfld();
        else
        {
          if (f->getNfld() != nfld2)
          {
            PyErr_SetString(PyExc_ValueError,
                            "diffArrays: number of fields must be equal in 2nd list.");
            return NULL;
          }
        }
      }
      else printf("Warning: diffArrays: one array is empty.\n");
    }
    else if (res == 2)
    {
      ni2.push_back(-1); nj2.push_back(-1); nk2.push_back(-1);
      field2.push_back(f);
      cn2.push_back(cn);
    }
    else
    {
      printf("Warning: diffArrays: not an array. Array skipped...\n");
    }
  }

  if (nfld1 != nfld2) 
    printf("Warning: diffArrays: number of fields are different. Only common fields are compared.\n");

  if (field1.size() != field2.size())
    printf("Warning: diffArrays: the number of arrays is different in arrays1 and arrays2.\n");

  if (field1.size() == 0 || field2.size() == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "diffArrays: no valid field to compare.");
    return NULL;
  }

  /*-----------------------*/
  /* Computing error field */
  /*-----------------------*/

  // Extract position of common variables: 
  // nb: x,y,z are elements 0,1,2 in pos1 and pos2
  vector<E_Int> pos1; vector<E_Int> pos2;
  char varString[K_ARRAY::VARSTRINGLENGTH];
  char varStringl[K_ARRAY::VARSTRINGLENGTH];
  E_Int pos = 
    K_ARRAY::getPosition( varString1, varString2, pos1, pos2, varStringl);
  K_ARRAY::addPrefixInVarString(varStringl, "D", varString);

  if (pos == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "diffArrays: no common variables found.");
    return NULL;
  }
  else if (pos == 0) // des variables sont differentes
  {
    printf("Warning: diffArrays: some variables are different in both arguments. Only common fields are compared.\n");
  }
  
  E_Boolean coordPresent = true;
  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
  if (posx1 == -1 || posy1 == -1 || posz1 == -1)
  {
    coordPresent = false;
  }
  posx1++; posy1++; posz1++;
  
  E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    coordPresent = false;
  }
  posx2++; posy2++; posz2++;
  
  // We now try to identify the blocks by sizes and position  
  int sizefield1 = field1.size();

  /*------------------------------*/
  /* Return errors as a py object */
  /*------------------------------*/
  PyObject* l = PyList_New(0);
  vector<E_Float*> errors;
  for (int i = 0; i < sizefield1; i++)
  {
    if (ni1[i] != -1)
    {
      tpl = K_ARRAY::buildArray2(pos1.size(), varString, 
                                 ni1[i], nj1[i], nk1[i]);
    }
    else
    {
      E_Int csize = cn1[i]->getSize()*cn1[i]->getNfld(); 
      tpl = K_ARRAY::buildArray2(pos1.size(), varString, field1[i]->getSize(),
                                 cn1[i]->getSize(), -1, elt1[i], false, csize);
      E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
      K_KCORE::memcpy__(cnnp, cn1[i]->begin(), csize);
    }
    errors.push_back(K_ARRAY::getFieldPtr(tpl));
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }

  for (int i1 = 0; i1 < sizefield1; i1++)
  {
    FldArrayF& f1 = *field1[i1];
    E_Int n1 = f1.getSize();
    FldArrayF error(n1, pos1.size(), errors[i1], true);
    error.setAllValuesAt(1.e6);
    E_Boolean found = false;
    
    found = searchField2(f1, error,
                         field2,
                         pos1, pos2,
                         posx1, posy1, posz1,
                         posx2, posy2, posz2,
                         coordPresent);

    if (found == false && coordPresent == true)
      found = searchField2(f1, error,
                           field2,
                           pos1, pos2,
                           posx1, posy1, posz1,
                           posx2, posy2, posz2,
                           false);

    if (found == false)
      printf("Warning: diffArrays: a field on a block can not be compared.\n");
  }
  
  /*------------------*/
  /* Little cleaning  */
  /*------------------*/
  int sizefield2 = field2.size();
  for (int i = 0; i < sizefield1; i++) 
  { 
    if (ni1[i] == -1) { RELEASESHAREDU(object1[i], field1[i], cn1[i]); }
    else { RELEASESHAREDS(object1[i], field1[i]); } 
  }
  for (int i = 0; i < sizefield2; i++) 
  {
    if (ni2[i] == -1) { RELEASESHAREDU(object2[i], field2[i], cn2[i]); }
    else { RELEASESHAREDS(object2[i], field2[i]); } 
  }
  field1.clear(); field2.clear(); 
  object1.clear(); object2.clear();
  ni2.clear(); nj2.clear(); nk2.clear();
  
  return l;
}

//=============================================================================
/* Diff3: */
//=============================================================================
PyObject* K_CONVERTER::diff3(PyObject* arrays1, PyObject* arrays2, PyObject* arrays3)
{
  PyObject* tpl;
  E_Int res;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk;
  char* varString1; char* varString2; char* varString3;
  char* eltType1; char* eltType2; char* eltType3;
  
  /* extract arrays1 : reference mesh in centers */
  E_Int n1 = PyList_Size(arrays1);
  vector<FldArrayF*> fieldmr;
  vector<E_Int> immr; vector<E_Int> jmmr; vector<E_Int> kmmr;

  E_Int nfld1 = -1;
  E_Int nfld2 = -1;
  E_Int nfld3 = -1;
  
  for (int i = 0; i < n1; i++)
  {
    tpl = PyList_GetItem(arrays1, i);
    res = K_ARRAY::getFromArray(tpl, varString1, f, ni, nj, nk, cn, eltType1);

    if ( res == 1 )
    {
      if ( ni*nj*nk > 0 )
      {
        immr.push_back(ni);
        jmmr.push_back(nj);
        kmmr.push_back(nk);
        fieldmr.push_back(f);

        if (nfld1 == -1) nfld1 = f->getNfld();
        else 
        {
          if (f->getSize() != nfld1)
          {
            PyErr_SetString(PyExc_TypeError,
                            "diff3: number of fields must be equal in 1st list.");
            return NULL;
          }
        }
      }
      else printf("Warning: diff3: one array is empty.\n");
    }
    else 
    {
      printf("Warning: diff3: only for structured arrays. Array skipped.\n");
    }
  }
  
  /* extract arrays2 : reference solution in centers */
  vector<E_Int> imsr; vector<E_Int> jmsr; vector<E_Int> kmsr;
  E_Int n2 = PyList_Size(arrays2);
  vector<FldArrayF*> fieldsr;

  for (int i = 0; i < n2; i++)
  {
    tpl = PyList_GetItem(arrays2, i);
    res = K_ARRAY::getFromArray(tpl, varString2, f, ni, nj, nk, cn, eltType2);
    if (res == 1)
    {
      if (ni*nj*nk > 0)
      {
        imsr.push_back(ni);
        jmsr.push_back(nj);
        kmsr.push_back(nk);
        fieldsr.push_back(f);
        if (nfld2 == -1) nfld2 = f->getNfld();
        else 
        {
          if (f->getSize() != nfld2)
          {
            PyErr_SetString(PyExc_TypeError,
                            "diff3: number of fields must be equal in 2nd list.");
            return NULL;
          }
        }
      }
      else printf("Warning: diff3: one array is empty.\n");
    }
    else 
    {
      printf("Warning: diff3: only for structured arrays. Array skipped...\n");
    }
  }  

  /* extract arrays3 : reading tested solution in centers */
  vector<E_Int> imsm; vector<E_Int> jmsm; vector<E_Int> kmsm;
  E_Int n3 = PyList_Size(arrays3);
  vector<FldArrayF*> fieldsm;

  for (int i = 0; i < n3; i++)
  {
    tpl = PyList_GetItem(arrays3, i);
    res = K_ARRAY::getFromArray(tpl, varString3, f, ni, nj, nk, cn, eltType3);
    if (res == 1)
    {
      if (ni*nj*nk > 0)
      {
        imsm.push_back(ni);
        jmsm.push_back(nj);
        kmsm.push_back(nk);
        fieldsm.push_back(f);
        if (nfld3 == -1) nfld3 = f->getNfld();
        else 
        {
          if (f->getSize() != nfld3)
          {
            PyErr_SetString(PyExc_TypeError,
                            "diff3: number of fields must be equal in 3rd list.");
            return NULL;
          }
        }
      }
      else printf("Warning: diff3: one array is empty.\n");
    }
    else 
    {
      printf("Warning: diff3: only for structured arrays. Array skipped...\n");
    }
  }  
  
  if (nfld1 != nfld2 || nfld1 != nfld3 || nfld2 != nfld3) 
    printf(" Warning: diff3: number of fields are different. Only common fields are compared.\n");
  
  /* Building the precond for all grids */
  printf("INFO: Preconditionning...");
  vector<K_INTERP::InterpAdt*> adts;
  vector<FldArrayF*> errors;
  E_Int nzone = 0;
  E_Int isBuilt;
  E_Int sizefieldmr = fieldmr.size();
  vector<E_Int> niet;  vector<E_Int> njet; vector<E_Int> nket; 
  vector<FldArrayF*> vectOfExtCenters;
  for (E_Int v = 0 ; v < sizefieldmr; v++)
  {
    FldArrayF& f = *fieldmr[v];
    E_Int ni = immr[nzone];
    E_Int nj = jmmr[nzone];
    E_Int nk = kmmr[nzone];
    FldArrayF fn(ni*nj*nk,3);

    E_Int posx = K_ARRAY::isCoordinateXPresent(varString1);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString1);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString1);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      PyErr_SetString(
        PyExc_TypeError,
        "diff3: coordinates don't exist in first array.");
    }
    posx++; posy++; posz++;
    
    for (E_Int ind = 0; ind < ni*nj*nk; ind++)
    {
      fn(ind, 1) = f(ind, posx);
      fn(ind, 2) = f(ind, posy);
      fn(ind, 3) = f(ind, posz);
    }
      
    E_Int nie = ni+2; E_Int nje = nj+2; E_Int nke = nk+2;
    niet.push_back(nie); njet.push_back(nje); nket.push_back(nke);
    
    FldArrayF* fec = new FldArrayF(nie*nje*nke,fn.getNfld());//maillage centres etendus
    K_LOC::node2ExtCenterStruct(ni, nj, nk, fn, nie, nje, nke, *fec);
    vectOfExtCenters.push_back(fec);
    //L'adt est toujours construit car le maillage est en centres etendus. Pas de verif isBuilt=1
    K_INTERP::InterpAdt* myAdt = new K_INTERP::InterpAdt(nie*nje*nke, fec->begin(1), fec->begin(2), fec->begin(3),
                                                         &nie, &nje, &nke, isBuilt);
    adts.push_back(myAdt);
    nzone++;
  }
  printf("done.\n");

  /* Interpolating from one grid to another */
  /* Beware: performed interpolation are order 2 so add an error
     to your solution */
  printf("Info: diff3: Interpolating...");
    
  // Creation des donnees pour l interpolation 
  // attention : centres etendus ici
  K_INTERP::InterpAdt::InterpolationType interpType = K_INTERP::InterpAdt::O2CF;
  
  E_Int nindi = 1;
  E_Int ncf = 8;
  FldArrayI indi(nindi);
  FldArrayF cf(ncf);
  
  // Compare variables in strings : pos2 and pos3 yield the position of 
  // common variables of strings. varString is the common variable string
  char varString[K_ARRAY::VARSTRINGLENGTH];
  vector<E_Int> pos2;
  vector<E_Int> pos3;
  E_Int npos2, npos3;

  E_Int pos = K_ARRAY::getPosition(varString2, varString3, pos2, pos3, varString);
  if (pos == -1)
    PyErr_SetString(PyExc_TypeError,
                    "diffArrays: no common variables in arrays");
  else if (pos == 0) 
    printf("Warning: diff3: some variables are different in 2nd and 3rd arguments. Only common variables will be compared.\n");
  
  E_Int posx3 = K_ARRAY::isCoordinateXPresent(varString3);
  E_Int posy3 = K_ARRAY::isCoordinateYPresent(varString3);
  E_Int posz3 = K_ARRAY::isCoordinateZPresent(varString3);
  E_Int poscelln = K_ARRAY::isCellNatureField1Present(varString3);
  poscelln++;

  E_Int sizefieldsm = fieldsm.size();
  E_Int sizepos2 = pos2.size(); 
  FldArrayI cn2(0);
  for (E_Int v = 0 ; v < sizefieldsm ; v++)
  {
    FldArrayF& f = *fieldsm[v];
    FldArrayF* error = new FldArrayF(f.getSize(), f.getNfld());
    error->setAllValuesAt(1.e6);
    E_Float* error1 = error->begin(1);
    E_Float* error2 = error->begin(2);
    E_Float* error3 = error->begin(3);

    if (posx3 == -1 || posy3 == -1 || posz3 == -1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "diff3: coordinates don't exist in third array.");
      return NULL;
    }
    posx3++; posy3++; posz3++;
    
    for (E_Int ind = 0; ind < f.getSize(); ind++)
    {
      E_Float x = f(ind,posx3);
      E_Float y = f(ind,posy3);
      E_Float z = f(ind,posz3);
      
      error1[ind] = x;
      error2[ind] = y;
      error3[ind] = z;

      E_Int no = 0;
      E_Int sizeadts = adts.size();
      for (E_Int iadt = 0; iadt < sizeadts; iadt++)
      {
        K_INTERP::InterpAdt* adt = adts[iadt];        
        E_Float voli = 0.; E_Int type = 0; E_Int noblk = 0;
        E_Int nie = niet[iadt];
        E_Int nje = njet[iadt];
        E_Int nke = nket[iadt];
        short found = K_INTERP::getInterpolationCell(x, y, z, adt, vectOfExtCenters[iadt], &nie, &nje, &nke, NULL,
                                                      1, 2, 3, 0, voli, indi, cf, type, noblk, interpType);
        
        FldArrayF& f2 = *fieldsr[no];
  
        if (found != 0)
        {          
          E_Int nerr = 4;          
          for (E_Int v = 0; v < sizepos2; v++)
          {
            npos3 = pos3[v];
            npos2 = pos2[v];
            
            //verifier que variables != coordonnee
            if (npos3 != posx3 && npos3 != posy3 && npos3 != posz3)
            {
              E_Float interp = 0.;
              K_INTERP::compOneInterpolatedValue(indi.begin(), cf, f2.begin(npos2), &cn2, NULL, NULL,
                                                  type, interp);
  
              E_Float delta = E_abs( f(ind, npos3) - interp );
              
              if (delta  < E_abs((*error)(ind,nerr)))
                (*error)(ind,nerr) = delta;
              if (poscelln != 0 && npos3 == poscelln)
                (*error)(ind, nerr) = f(ind, poscelln);
              nerr++;
            }
          }
        } // if found != 0         
        no++;
      }
    }
    errors.push_back(error);
  }
  printf("done.\n");

  /* Checking non interpolated points */
  E_Int sizeerrors = errors.size();
  for (E_Int v = 0; v < sizeerrors; v++)
  {
    FldArrayF& e = *errors[v];
    
    for (E_Int ind = 0; ind < e.getSize(); ind++)
    {
      if (e(ind,4) > 9.e5) // non interpolated point
      {
        for (E_Int v1 = 4; v1 <= e.getNfld(); v1++)
        {
          e(ind,v1) = 0.; // set to 0. and cellN to 0 also
        }
      }
    }
  }
  /* Checking solid points */
  if ( poscelln != 0) // chimera
  {
    for (E_Int v = 0; v < sizeerrors; v++)
    {
      FldArrayF& e = *errors[v]; 
      
      for (E_Int ind = 0; ind < e.getSize(); ind++)
      {
        if (e(ind,poscelln) == 0) // solid point
        {
          for (E_Int v = 4; v < e.getNfld(); v++)
          {
            e(ind,v) = 0.; // set to 0.
          }
        }
      }
    }
  }
  for (unsigned int no = 0; no < vectOfExtCenters.size(); no++)
  {
    delete vectOfExtCenters[no];
    delete adts[no];
  }
  /* Sauvegarde de errors sous forme de liste python */
  PyObject* l = PyList_New(0);

  for (int i = 0; i < sizeerrors; i++)
  {
    tpl = K_ARRAY::buildArray(*errors[i], varString, 
                              imsm[i], jmsm[i], kmsm[i]);
    delete errors[i];
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  return l;
}

//=============================================================================
E_Boolean K_CONVERTER::searchField2(FldArrayF& f1,
                                    FldArrayF& error,
                                    vector<FldArrayF*>& field2,
                                    vector<E_Int>& pos1,
                                    vector<E_Int>& pos2,
                                    E_Int posx1, E_Int posy1, E_Int posz1,
                                    E_Int posx2, E_Int posy2, E_Int posz2,
                                    E_Boolean coordPresent)
{
  const E_Float EPS = 1.e-12;
  int sizefield2 = field2.size();
  E_Int sizepos1 = pos1.size();

  // Check for invalid values in f1
  E_Int n1 = f1.getSize();
  E_Int nthreads = __NUMTHREADS__;
  E_Int *hasNan = new E_Int [nthreads]; for (E_Int i = 0; i < nthreads; i++) hasNan[i] = 0;
  E_Int* hasInf = new E_Int [nthreads]; for (E_Int i = 0; i < nthreads; i++) hasInf[i] = 0;

#pragma omp parallel default(shared)
  {
    E_Int ithread = __CURRENT_THREAD__;
    for (E_Int v = 0; v < sizepos1; v++)
    {
      E_Int npos1 = pos1[v];
      E_Float* f1p = f1.begin(npos1);
#pragma omp parallel for
      for (E_Int ind = 0; ind < n1; ind++)
      {
        if (isnan(f1p[ind])) hasNan[ithread] = 1;
        if (isinf(f1p[ind])) hasInf[ithread] = 1;
      }
    }
  }

  for (E_Int i = 0; i < nthreads; i++)
  { 
    if (hasNan[i] == 1)
    {
      delete [] hasNan; delete [] hasInf;
      PyErr_SetString(PyExc_TypeError,
                      "diffArrays: arrays1 contains Nan values.");
      return false;
    }
    if (hasInf[i] == 1)
    {
      delete [] hasNan; delete [] hasInf;
      PyErr_SetString(PyExc_TypeError,
                      "diffArrays: arrays1 contains infinite values.");
      return false;
    }
  }

  // PRe requisite
  for (int i2 = 0; i2 < sizefield2; i2++)
  {
    FldArrayF& f2 = *field2[i2];
    E_Int n2 = f2.getSize();
      
    if (n1 == 0 && n2 == 0)
    { 
      delete [] hasNan; delete [] hasInf;
      return true;
    }
  }

  E_Int* found = new E_Int [nthreads]; for (E_Int i = 0; i < nthreads; i++) found[i] = 0;

#pragma omp parallel default(shared)
  {
    E_Int ithread = __CURRENT_THREAD__;

    for (int i2 = 0; i2 < sizefield2; i2++)
    {
      FldArrayF& f2 = *field2[i2];
      E_Int n2 = f2.getSize();
      if (
        (n1 == n2 && coordPresent == false) ||
        (n1 == n2 &&
         E_abs(f2(0,posx2) - f1(0,posx1)) < EPS &&
         E_abs(f2(0,posy2) - f1(0,posy1)) < EPS &&
         E_abs(f2(0,posz2) - f1(0,posz1)) < EPS &&
         E_abs(f2(n2-1,posx2) - f1(n1-1,posx1)) < EPS &&
         E_abs(f2(n2-1,posy2) - f1(n1-1,posy1)) < EPS &&
         E_abs(f2(n2-1,posz2) - f1(n1-1,posz1)) < EPS &&
         E_abs(f2(n2/2,posx2) - f1(n1/2,posx1)) < EPS &&
         E_abs(f2(n2/2,posy2) - f1(n1/2,posy1)) < EPS &&
         E_abs(f2(n2/2,posz2) - f1(n1/2,posz1)) < EPS ))
      {
        for (E_Int v = 0; v < sizepos1; v++)
        {
          E_Int npos1 = pos1[v];
          E_Int npos2 = pos2[v];
          E_Float* f1p = f1.begin(npos1);
          E_Float* f2p = f2.begin(npos2);
          E_Float* errorp = error.begin(v+1);

#pragma omp for
          for (E_Int ind = 0; ind < n1; ind++)
          {
            if (isnan(f2p[ind])) { hasNan[ithread] = 1; }
            if (isinf(f2p[ind])) { hasInf[ithread] = 1; }
            errorp[ind] = E_abs(f1p[ind] - f2p[ind]);
          }
        }
        found[ithread] = 1;
      }
    } // loop on field2
  }

  for (E_Int i = 0; i < nthreads; i++)
  { 
    if (hasNan[i] == 1)
    {
      delete [] hasNan; delete [] hasInf; delete [] found;
      PyErr_SetString(PyExc_TypeError,
                      "diffArrays: arrays2 contains Nan values.");
      return false;
    }
    if (hasInf[i] == 1)
    {
      delete [] hasNan; delete [] hasInf; delete [] found;
      PyErr_SetString(PyExc_TypeError,
                      "diffArrays: arrays2 contains infinite values.");
      return false;
    }
  }

  for (E_Int i = 0; i < nthreads; i++)
  { 
    if (found[i] == 1) 
    { 
      delete [] hasNan; delete [] hasInf; delete [] found;
      return true;
    }
  }
  delete [] hasNan; delete [] hasInf; delete [] found;
  return false;
}
