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
# include <stdio.h>
# include <string.h>

# include "post.h"
# include <vector>

using namespace K_FLD;
using namespace std;

//=============================================================================
/* Calcul de l'integrale de la solution */
// ============================================================================
PyObject* K_POST::dinteg(PyObject* self, PyObject* args)
{
  PyObject* coordArrays;
  PyObject* FArrays;
  PyObject* ratioArrays;
  
  if (!PYPARSETUPLE_(args, OOO_,  &coordArrays, &FArrays, &ratioArrays))
    return NULL;

  // verification des listes
  if (PyList_Check(coordArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "integ: first argument must be a list.");
    return NULL;
  }
  if (PyList_Check(FArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "integ: second argument must be a list.");
    return NULL;
  }
  if (PyList_Check(ratioArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "integ: third argument must be a list.");
    return NULL;
  }
  
  E_Int nCoordArrays = PyList_Size(coordArrays);
  E_Int nFArrays = PyList_Size(FArrays);
  E_Int nRatioArrays = PyList_Size(ratioArrays);
  
  // verification du nb egal d'elements dans chq liste
  if (nRatioArrays == 0)
  { 
    if (nCoordArrays != nFArrays)
    {
      PyErr_SetString(PyExc_ValueError, 
                      "integ: number of zones in 1st and 2nd arguments must be equal.");
      return NULL;    
    }
  }
  else
  {
    if (nCoordArrays != nFArrays || 
        nCoordArrays != nRatioArrays || 
        nRatioArrays != nFArrays)
    {
      PyErr_SetString(PyExc_ValueError, 
                      "integ: number of zones in 1st, 2nd and 3rd arguments must be equal.");
      return NULL;
    }
  }

  PyObject* coordObj;
  PyObject* FObj;
  PyObject* ratioObj;
  E_Int nic, njc, nkc, nif, njf, nkf, nir, njr, nkr;
  char* varStringc; char* varStringf; char* varStringr;
  FldArrayF* fc; FldArrayF* ff;
  FldArrayF* ratio;
  E_Int nFld = -1;
  char varString0[K_ARRAY::VARSTRINGLENGTH];
  E_Int sizef = 0;
  E_Int center2node = 2;  // set to 1 if coord is in nodes and F in centers
                          // set to 0 if coord and F have the same size
  E_Int case1D;           // set to 1 if linear integration
  
  char* eltTypec; FldArrayI* cnc;
  char* eltTypef; FldArrayI* cnf;
  char* eltTyper; FldArrayI* cnr;

  E_Int resc = 0;
  E_Int resf = 0;
  E_Int resr = 0;
  FldArrayF resultat;
  E_Int res = -1;
  
  for (E_Int i = 0; i < nCoordArrays; i++)
  {
    coordObj = PyList_GetItem(coordArrays, i);
    resc = K_ARRAY::getFromArray3(coordObj, varStringc, fc, 
                                  nic, njc, nkc, cnc, eltTypec); 
    if (resc != 1 && resc != 2)
    {
      PyErr_SetString(PyExc_TypeError,
                      "integ: coord is not a valid array.");
      return NULL;
    }
    E_Int posx = K_ARRAY::isCoordinateXPresent(varStringc);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varStringc);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varStringc);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      printf("Warning: integ: coordinates not found in array " SF_D_ ".", i+1);
      printf(" Array skipped...\n");
      RELEASESHAREDB(resc, coordObj, fc, cnc);
      goto next; 
    }
    posx++; posy++; posz++;

    FObj = PyList_GetItem(FArrays, i);
    resf = K_ARRAY::getFromArray3(FObj, varStringf, ff, 
                                  nif, njf, nkf, cnf, eltTypef); 

    if (resf != 1 && resf != 2)
    {
      RELEASESHAREDB(resc, coordObj, fc, cnc);
      RELEASESHAREDS(FObj, ff);
      PyErr_SetString(PyExc_TypeError, 
                      "integ: field is not a valid array.");
      return NULL;
    }
    // Verifie le nombre de variables: doit etre le meme entre tous
    if (nFld == -1) // premier passage
    {
      strcpy(varString0, varStringf);
      nFld = ff->getNfld();
      resultat.malloc(nFld);
      resultat.setAllValuesAtNull();
    }
    else 
    {
      if (ff->getNfld() != nFld)
      {
        printf("Warning: integ: array " SF_D_ " doesnt have a valid number of variables.", i+1); 
        printf("Array skipped...\n");
        RELEASESHAREDB(resc, coordObj, fc, cnc);
        RELEASESHAREDB(resf, FObj, ff, cnf);
        goto next;
      }

      // Check is variables are ordered in the same way
      E_Int ids = K_ARRAY::compareVarStrings(varString0, varStringf);
      if (ids == -1) // varstrings are different
      {
        printf("Warning: integ: variables are in a different order than first array.");
        printf(" Array skipped...\n");
        RELEASESHAREDB(resc, coordObj, fc, cnc);
        RELEASESHAREDB(resf, FObj, ff, cnf);
        goto next;
      }
    }
    
    // Cas structure
    if (resc == 1 && resf == 1)
    { 
      sizef = nif*njf*nkf;      

      if (nic > 1 && njc > 1 && nkc > 1) 
      {
        printf("Warning: integ: 3D arrays not valid. Array skipped...\n");
        RELEASESHAREDB(resc, coordObj, fc, cnc);
        RELEASESHAREDB(resf, FObj, ff, cnf);
        goto next;
      }
      // 1D ou 2D ?
      case1D = 0;
      if ((nic == 1 && njc == 1) ||
          (nic == 1 && nkc == 1) ||
          (njc == 1 && nkc == 1))
        case1D = 1;
      
      // coord et F de meme taille ?
      if (nic == nif && njc == njf && nkc == nkf)
        center2node = 0;
      else if (nic == nif+1 && njc == njf+1 && nkc == nkf+1)
        center2node = 1;
      else
      {
        if ((nic == 1 && njc == njf+1 && nkc == nkf+1) || 
            (njc == 1 && nic == nif+1 && nkc == nkf+1) || 
            (nkc == 1 && nic == nif+1 && njc == njf+1))
          center2node = 1;
        else if ((nic == 1 && njc == 1 && nkc == nkf+1) || 
                 (nic == nif+1 && njc == 1 && nkc == 1) || 
                 (nic == 1 && njc == njf+1 && nkc == 1))
          center2node = 1;
        else 
        {
          printf("Warning: integ: coord and field arrays do not represent the same zone.");
          printf(" Array skipped...\n");
          RELEASESHAREDB(resc, coordObj, fc, cnc);
          RELEASESHAREDB(resf, FObj, ff, cnf);
          goto next;
        }
      }
      
      // Valeur du ratio ?
      if (nRatioArrays == 0)
      {
        ratio = new FldArrayF(sizef);
        ratio->setAllValuesAt(1.);
      }
      else // coord + F + r
      {
        ratioObj = PyList_GetItem(ratioArrays, i);
        resr = K_ARRAY::getFromArray3(ratioObj, varStringr, ratio, 
                                      nir, njr, nkr, cnr, eltTyper); 
        if (resr != 1)
        {
          RELEASESHAREDB(resr, ratioObj, ratio, cnr);
          printf("Warning: integ: ratio " SF_D_ " is an invalid array. Set to 1.", i+1);
          ratio = new FldArrayF(sizef);
          ratio->setAllValuesAt(1.);
        }
      }

      // integ sur chaque bloc
      res = 0;

      /*
      if (case1D == 1)
          res = dintegStruct1D(nic, njc, nkc, center2node, posx, posy, posz, 
                      *fc, *ff, *ratio, resultat);
      else
          res = dintegStruct2D(nic, njc, nkc, center2node, posx, posy, posz,
                      *fc, *ff, *ratio, resultat);
      */
      if (res == 0) 
      {
        RELEASESHAREDS(coordObj, fc);
        RELEASESHAREDS(FObj, ff);
        if (nRatioArrays == 0 || resr != 1) delete ratio;
        else RELEASESHAREDS(ratioObj, ratio);
        PyErr_SetString(PyExc_ValueError,
                        "integ: integration computation fails.");
        return NULL;
      }
      RELEASESHAREDS(coordObj, fc);
      RELEASESHAREDS(FObj, ff);
      if (nRatioArrays == 0 || resr != 1) delete ratio;
      else RELEASESHAREDS(ratioObj, ratio);
    }
    else if (resc == 2 && resf == 2) // ME
    {
      // check if field is cell or node centered
      E_Int l = strlen(eltTypef);
      if (eltTypef[l-1] == '*') center2node = 1;
      else center2node = 0;

      res = 1;
      std::vector<char*> eltTypecs, eltTypefs;
      K_ARRAY::extractVars(eltTypec, eltTypecs);

      // check if elt is valid (BAR, QUAD, TRI)
      case1D = 0;
      for (size_t ic = 0; ic < eltTypecs.size(); ic++)
      {
        if (strcmp(eltTypecs[ic], "BAR") == 0) case1D = 1;
        else if ((strcmp(eltTypecs[ic], "QUAD") == 0) || (strcmp(eltTypecs[ic], "TRI") == 0)) case1D = 0;
        else res = 0;
      }

      for (size_t ic = 0; ic < eltTypecs.size(); ic++) delete [] eltTypecs[ic];

      if (res == 0)
      {
        RELEASESHAREDU(coordObj, fc, cnc);
        RELEASESHAREDU(FObj, ff, cnf);
        PyErr_SetString(PyExc_ValueError, 
                        "integ: only TRI, QUAD, or BAR unstructured arrays are possible.");
        return NULL;
      }
      
      sizef = ff->getSize();
   
      // Valeur du ratio ?
      if (nRatioArrays == 0)
      {
        ratio = new FldArrayF(sizef);
        ratio->setAllValuesAt(1.);
      }
      else // coord + F + r
      {
        ratioObj = PyList_GetItem(ratioArrays, i);
        resr = K_ARRAY::getFromArray3(
          ratioObj, varStringr, ratio, nir, njr, nkr, cnr, eltTyper); 

        if (resr != 2)
        {
          if (resr == 1) RELEASESHAREDS(ratioObj, ratio);
          printf("Warning: integ: ratio " SF_D_ " is an invalid array. Set to 1.", i+1);
          ratio = new FldArrayF(sizef);
          ratio->setAllValuesAt(1.);
        }
      }

      // integ sur chaque bloc
      res = 0;

      /*
      if (case1D == 1)
        res = dintegUnstruct1D(center2node, posx, posy, posz, 
          *cnc, eltTypec, *fc, *ff, *ratio, resultat);
      else
        res = dintegUnstruct2D(center2node, posx, posy, posz, 
          *cnc, eltTypec, *fc, *ff, *ratio, resultat);
      */

      if (res == 0)
      {
        RELEASESHAREDU(coordObj, fc, cnc);
        RELEASESHAREDU(FObj, ff, cnf);
        if (nRatioArrays == 0 || resr != 2) delete ratio;
        else RELEASESHAREDB(resr, ratioObj, ratio, cnr);
        PyErr_SetString(PyExc_ValueError,
                        "integ: integration computation fails.");
        return NULL;
      }
      RELEASESHAREDU(coordObj, fc, cnc);
      RELEASESHAREDU(FObj, ff, cnf);
      if (nRatioArrays == 0 || resr != 2) delete ratio;
      else RELEASESHAREDB(resr, ratioObj, ratio, cnr);
    }
    else 
    {
      PyErr_SetString(PyExc_TypeError, 
                      "integ: not valid arrays.");
      return NULL;
    }
    next:;
  }

  if (res == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "integ: integration failed.");
    return NULL;
  }
  PyObject* tpl;
  PyObject* l = PyList_New(0);
  
  for (E_Int i = 0; i < nFld; i++)
  {
    tpl = Py_BuildValue(R_, resultat[i]);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  return l;
}

//=============================================================================
// Integre "surfaciquement" les grandeurs de F comme des scalaires
// Retourne 1 si succes, 0 si echec
//=============================================================================
E_Int K_POST::dintegStruct2D(E_Int ni, E_Int nj, E_Int nk, 
                             E_Int center2node, E_Int posx, E_Int posy, E_Int posz,
                             FldArrayF& coord, FldArrayF& F, 
                             FldArrayF& ratio, FldArrayF& resultat)
{
  E_Int NI, NJ;
  E_Int numberOfVariables = F.getNfld();

  if      (nk == 1) {NI = ni; NJ = nj;}
  else if (nj == 1) {NI = ni; NJ = nk;}
  else if (ni == 1) {NI = nj; NJ = nk;}
  else return 0;

  E_Int ncells = (NI-1)*(NJ-1);
    
  // point
  E_Float* xt = coord.begin(posx); // coords
  E_Float* yt = coord.begin(posy);
  E_Float* zt = coord.begin(posz);
  E_Float* dxt = coord.begin(4); // tangent
  E_Float* dyt = coord.begin(5);
  E_Float* dzt = coord.begin(6);

  auto x = new E_Float [3*NI*NJ];
  for (E_Int i = 0; i < NI*NJ; i++)
  {
    x[3*i] = xt[i];
    x[3*i+1] = yt[i];
    x[3*i+2] = zt[i];
  }

  // tangent
  auto dx = new E_Float [3*NI*NJ];
  for (E_Int i = 0; i < NI; i++)
  {
    dx[3*i] = dxt[i];
    dx[3*i+1] = dyt[i];
    dx[3*i+2] = dzt[i];
  }

  // gradient
  auto df = new E_Float [3*NI*NJ];

  // actives
  auto ax = new adouble[3*NI*NJ];

  trace_on(0);

  for (E_Int i = 0; i < 3*NI*NJ; i++) 
  {
    ax[i] <<= xt[i];
  }

  // Compute surface of each "block" i cell, with coordinates coord
  adouble* surf = new adouble [ncells];
  adouble result = 0;
  adouble* results = new adouble [numberOfVariables];

  /*
  K_METRIC::d__compSurfStruct2D(
    NI, NJ, 1,
    ax, ax+3*NI*NJ, ax+6*NI*NJ,
    surf);
  */

  if (center2node == 1) 
  {
    // Compute integral, coordinates defined in node 
    // and field F in center 
    for (E_Int n = 1; n <= numberOfVariables; n++)
    { 
      /*
      K_POST::d__integStructCellCenter2D(
        NI-1, NJ-1,
        ratio.begin(), surf, F.begin(n),
        result);
      */
      results[n-1] += result;
    }
  }
  else
  {
    // Compute integral, coordinates and field have the same size
    for (E_Int n = 1; n <= numberOfVariables; n++)
    {
      /*
      K_POST::d__integStructNodeCenter2D(
        NI, NJ,
        ratio.begin(), surf, F.begin(n),
        result);
      */
      results[n-1] += result;
    }
  }

  trace_off();
  gradient(0, 3*NI*NJ, x, df);

  resultat.malloc(NI*NJ, 3*numberOfVariables);
  for (E_Int n = 0; n < numberOfVariables; n++)
  {
    E_Float* dFxt = resultat.begin(3*n+1);
    E_Float* dFyt = resultat.begin(3*n+2);
    E_Float* dFzt = resultat.begin(3*n+3);

    for (E_Int i = 0; i < NI*NJ; i++)
    {
      dFxt[i] = df[3*n+3*i];
      dFyt[i] = df[3*n+3*i+1];
      dFzt[i] = df[3*n+3*i+2];
    }
  }

  // store df in resultat
  delete [] surf; delete [] results;
  delete [] ax; delete [] df;

  return 1;
}

//=============================================================================
// Integre "surfaciquement" les grandeurs de F comme des scalaires
// Retourne 1 si succes, 0 si echec
//=============================================================================
E_Int K_POST::dintegStruct1D(E_Int ni, E_Int nj, E_Int nk, 
                             E_Int center2node, E_Int posx, E_Int posy, E_Int posz,
                             FldArrayF& coord, FldArrayF& F, 
                             FldArrayF& ratio, FldArrayF& resultat)
{
  E_Int NI;
  E_Int numberOfVariables = F.getNfld();

  if      (ni > 1) NI = ni;
  else if (nj > 1) NI = nj;
  else if (nk > 1) NI = nk;
  else return 0;

  E_Int ncells = (NI-1);

  // point
  E_Float* xt = coord.begin(posx); // coords
  E_Float* yt = coord.begin(posy);
  E_Float* zt = coord.begin(posz);
  E_Float* dxt = coord.begin(4); // tangent
  E_Float* dyt = coord.begin(5);
  E_Float* dzt = coord.begin(6);

  auto x = new E_Float [3*NI];
  for (E_Int i = 0; i < NI; i++)
  {
    x[3*i] = xt[i];
    x[3*i+1] = yt[i];
    x[3*i+2] = zt[i];
  }

  // tangent
  auto dx = new E_Float [3*NI];
  for (E_Int i = 0; i < NI; i++)
  {
    dx[3*i] = dxt[i];
    dx[3*i+1] = dyt[i];
    dx[3*i+2] = dzt[i];
  }

  // gradient
  auto df = new E_Float [3*NI];

  // actives
  auto ax = new adouble[3*NI];

  trace_on(0);

  for (E_Int i = 0; i < 3*NI; i++) 
  {
    ax[i] <<= xt[i];
  }

  // Compute surface of each "block" i cell, with coordinates coord
  auto length = new adouble [ncells];
  adouble result = 0;
  adouble* results = new adouble [numberOfVariables];

  /*
  K_METRIC::d__compSurfStruct1D(
    NI, 1, 1,
    coord.begin(posx), coord.begin(posy), coord.begin(posz),
    length.begin());
  */
  if (center2node == 1) 
  {
    // Compute integral, coordinates defined in node 
    // and field F in center 
    for (E_Int n = 1; n <= numberOfVariables; n++)
    { 
      /*
      K_POST::d__integStructCellCenter1D(
        NI-1,
        ratio.begin(), length, F.begin(n),
        result);
      */
      results[n-1] += result;
    }
  }
  else
  {
    // Compute integral, coordinates and field have the same size
    for (E_Int n = 1; n <= numberOfVariables; n++)
    {
      /*
      K_POST::d__integStructNodeCenter1D(
        NI,
        ratio.begin(), length, F.begin(n),
        result);
      */
      results[n-1] += result;
    }
  }
  
  
  trace_off();
  gradient(0, 3*NI, x, df);

  // store df in resultat
  resultat.malloc(NI, 3*numberOfVariables);
  for (E_Int n = 0; n < numberOfVariables; n++)
  {
    E_Float* dFxt = resultat.begin(3*n+1);
    E_Float* dFyt = resultat.begin(3*n+2);
    E_Float* dFzt = resultat.begin(3*n+3);

    for (E_Int i = 0; i < NI; i++)
    {
      dFxt[i] = df[3*n+3*i];
      dFyt[i] = df[3*n+3*i+1];
      dFzt[i] = df[3*n+3*i+2];
    }
  }
  
  delete [] length;
  delete [] ax; delete [] df;

  return 1;
}

//=============================================================================
// Integre les grandeurs de F comme des scalaires
// Retourne 1 si succes, 0 si echec
// 1D and 2D
//=============================================================================
E_Int K_POST::dintegUnstruct(E_Int center2node,
                             E_Int posx, E_Int posy, E_Int posz,
                             FldArrayI& cn, const char* eltType, FldArrayF& coord, 
                             FldArrayF& F, FldArrayF& ratio, 
                             FldArrayF& resultat)
{
  E_Int numberOfVariables = F.getNfld();
  E_Int ntotElts = 0;
  E_Int nc = cn.getNConnect();
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    ntotElts += nelts;
  }

  // point
  E_Float* xt = coord.begin(posx); // coords
  E_Float* yt = coord.begin(posy);
  E_Float* zt = coord.begin(posz);
  E_Float* dxt = coord.begin(4); // tangent
  E_Float* dyt = coord.begin(5);
  E_Float* dzt = coord.begin(6);

  E_Int npts = coord.getSize();

  auto x = new E_Float [3*npts];
  for (E_Int i = 0; i < npts; i++)
  {
    x[3*i] = xt[i];
    x[3*i+1] = yt[i];
    x[3*i+2] = zt[i];
  }

  // tangent
  auto dx = new E_Float [3*npts];
  for (E_Int i = 0; i < npts; i++)
  {
    dx[3*i] = dxt[i];
    dx[3*i+1] = dyt[i];
    dx[3*i+2] = dzt[i];
  }

  // gradient
  auto df = new E_Float [3*npts];

  // actives
  auto ax = new adouble[3*npts];

  trace_on(0);

  for (E_Int i = 0; i < 3*npts; i++) 
  {
    ax[i] <<= xt[i];
  }

  adouble* surf = new adouble [ntotElts];
  adouble* snx = new adouble [ntotElts];
  adouble* sny = new adouble [ntotElts];
  adouble* snz = new adouble [ntotElts];
  adouble result = 0;
  adouble* results = new adouble [numberOfVariables];

  /*
  K_METRIC::d__compSurfUnstruct(
    cn, eltType,
    coord.begin(posx), coord.begin(posy), coord.begin(posz),
    snx, sny, snz, surf);
  */
  if (center2node == 1) 
  {
    for (E_Int n = 1; n <= numberOfVariables; n++)
    {
      // Compute integral, coordinates defined in node 
      // and field F in center
      /*
      K_POST::d__integUnstructCellCenter(
        ntotElts,
        ratio.begin(), surf, F.begin(n),
        result);
      */
      results[n-1] += result;
    }
  }
  else
  {
    // Compute integral, coordinates and field have the same size
    for (E_Int n = 1; n <= numberOfVariables; n++)
    {
      /*
      K_POST::d__integUnstructNodeCenter(
        cn,
        ratio.begin(), surf, F.begin(n),
        result
      ); */
      results[n-1] += result;
    }
  }

  trace_off();
  gradient(0, 3*npts, x, df);

  // store df in resultat
  resultat.malloc(npts, 3*numberOfVariables);
  for (E_Int n = 0; n < numberOfVariables; n++)
  {
    E_Float* dFxt = resultat.begin(3*n+1);
    E_Float* dFyt = resultat.begin(3*n+2);
    E_Float* dFzt = resultat.begin(3*n+3);

    for (E_Int i = 0; i < npts; i++)
    {
      dFxt[i] = df[3*n+3*i];
      dFyt[i] = df[3*n+3*i+1];
      dFzt[i] = df[3*n+3*i+2];
    }
  }

  delete [] surf;
  delete [] snx;
  delete [] sny;
  delete [] snz;
  delete [] ax; delete [] df;

  return 1;
}

