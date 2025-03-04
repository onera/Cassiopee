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
# include <stdio.h>
# include <string.h>

# include "post.h"
# include <vector>

using namespace K_FLD;
using namespace std;

extern "C"
{
  void k6structsurft_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Int& ncells, 
    const E_Float* xt, const E_Float* yt, const E_Float* zt, 
    E_Float* length);

  void k6structsurf1dt_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Float* xt, const E_Float* yt, const E_Float* zt, 
    E_Float* length);

  void k6unstructsurf_(E_Int& npts, E_Int& nelts, E_Int& nedges, 
                       E_Int& nnodes, E_Int* cn, 
                       E_Float* coordx, E_Float* coordy, E_Float* coordz, 
                       E_Float* snx, E_Float* sny, E_Float* snz,
                       E_Float* surface);

  void k6unstructsurf1d_(E_Int& npts, E_Int& nelts, 
                         E_Int& nnodes, E_Int* cn, 
                         E_Float* coordx, E_Float* coordy, E_Float* coordz, 
                         E_Float* length);

  void k6integstruct_(const E_Int& ni, const E_Int& nj, 
                      E_Float* ratio, E_Float* surf, 
                      E_Float* F, E_Float& result);

  void k6integstruct1d_(const E_Int& ni, E_Float* ratio, E_Float* length, 
                        E_Float* F, E_Float& result);

  void k6integstructnodecenter_(const E_Int& ni, const E_Int& nj, 
                                E_Float* ratio, E_Float* surf, 
                                E_Float* F, E_Float& result);

  void k6integunstruct_(const E_Int& nbt, const E_Int& size, 
                        E_Int* cn, E_Float* ratio, E_Float* surf, 
                        E_Float* F, E_Float& result);
  
  void k6integunstruct1d_(const E_Int& nbt, const E_Int& size, 
                          E_Int* cn, E_Float* ratio, E_Float* length, 
                          E_Float* F, E_Float& result);

  void k6integstructnodecenter1d_(const E_Int& ni, E_Float* ratio, 
                                  E_Float* length, E_Float* F, 
                                  E_Float& result);

  void k6integunsnodecenter_(const E_Int& nbt, 
                             E_Float* ratio, E_Float* surf, 
                             E_Float* F, E_Float& result);
}

//=============================================================================
/* Calcul de l'integrale de la solution */
// ============================================================================
PyObject* K_POST::integ(PyObject* self, PyObject* args)
{
  PyObject* coordArrays;
  PyObject* FArrays;
  PyObject* ratioArrays;
  
  if (!PyArg_ParseTuple(args, "OOO",  &coordArrays, &FArrays, &ratioArrays))
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

  //
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
  E_Int center2node = 2; // set to 1 if coord is in nodes and F in centers
                         // set to 0 if coord and F have the same size
  E_Int case1D;      // set to 1 if linear integration
  
  char* eltTypec; FldArrayI* cnc;
  char* eltTypef; FldArrayI* cnf;
  char* eltTyper; FldArrayI* cnr;

  E_Int resc = 0;
  E_Int resf = 0;
  E_Int resr = 0;
  FldArrayF resultat;
  E_Int res = -1;
  
  for (int i = 0; i < nCoordArrays; i++)
  {
    coordObj = PyList_GetItem(coordArrays, i);
    resc = K_ARRAY::getFromArray(coordObj, varStringc, fc, 
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
      printf("Warning: integ: coordinates not found in array %d.", i+1);
      printf(" Array skipped...\n");
      delete fc; if (resc == 2) delete cnc;
      goto next; 
    }
    posx++; posy++; posz++;

    FObj = PyList_GetItem(FArrays, i);
    resf = K_ARRAY::getFromArray(FObj, varStringf, ff, 
                                 nif, njf, nkf, cnf, eltTypef); 

    if (resf != 1 && resf != 2)
    {
      delete ff; delete fc;
      if (resc == 2) delete cnc;
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
        printf("Warning: integ: array %d doesn t have a valid number of variables.", i+1); 
        printf("Array skipped...\n");
        delete ff; if (resf == 2) delete cnf;
        delete fc; if (resc == 2) delete cnc;
        goto next;
      }

      // Check is variables are ordered in the same way
      E_Int ids = K_ARRAY::compareVarStrings(varString0, varStringf);
      if (ids == -1) // varstrings are different
      {
        printf("Warning: integ: variables are in a different order than first array.");
        printf(" Array skipped...\n");
        delete ff; if (resf == 2) delete cnf;
        delete fc; if (resc == 2) delete cnc;
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
        delete ff; if (resf == 2) delete cnf;
        delete fc; if (resc == 2) delete cnc;
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
          delete ff; delete fc;
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
        resr = K_ARRAY::getFromArray(ratioObj, varStringr, ratio, 
                                     nir, njr, nkr, cnr, eltTyper); 
        if (resr != 1)
        {
          if (resr == 2) {delete ratio; delete cnr;}
          printf("Warning: integ: ratio %d is an invalid array. Set to 1.", i+1);
          ratio = new FldArrayF(sizef);
          ratio->setAllValuesAt(1.);
        }
      }

      // integ sur chaque bloc
      res = 0;
      switch (case1D)
      {
        case 1:
          res = integ11D(nic, njc, nkc, center2node, posx, posy, posz, 
                         *fc, *ff, *ratio, resultat);
          break;
        default:
          res = integ1(nic, njc, nkc, center2node, posx, posy, posz,
                       *fc, *ff, *ratio, resultat);
          break;
      }
    
      if (res == 0) 
      {
        delete ff; delete fc; delete ratio;
        PyErr_SetString(PyExc_ValueError,
                        "integ: integration computation fails.");
        return NULL;
      }
      delete ff; delete fc; delete ratio; 
    }
    else if (resc == 2 && resf == 2) // Cas non structure
    {
      if ((strcmp(eltTypec, "TRI") == 0 && strcmp(eltTypef, "TRI") == 0) || 
          (strcmp(eltTypec, "BAR") == 0 && strcmp(eltTypef, "BAR") == 0))
        center2node = 0;
      else if ((strcmp(eltTypec, "TRI") == 0 && strcmp(eltTypef, "TRI*") == 0) || 
               (strcmp(eltTypec, "BAR") == 0 && strcmp(eltTypef, "BAR*") == 0) )
        center2node = 1;
      else
      {
        delete fc; delete cnc;
        delete ff; delete cnf;
        PyErr_SetString(PyExc_ValueError, 
                        "integ: only TRI or BAR unstructured arrays are possible.");
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
        resr = K_ARRAY::getFromArray(
          ratioObj, varStringr, ratio, nir, njr, nkr, cnr, eltTyper); 

        if (resr != 2)
        {
          if (resr == 1) delete ratio;
          printf("Warning: integ: ratio %d is an invalid array. Set to 1.", i+1);
          ratio = new FldArrayF(sizef);
          ratio->setAllValuesAt(1.);
        }
      }

      // integ sur chaque bloc
      res = 0;
      if (strcmp(eltTypec, "TRI") == 0)
        res = integUnstruct1(center2node, posx, posy, posz, 
                             *cnc, *fc, *ff, *ratio, resultat);
      else
        res = integUnstruct11D(center2node, posx, posy, posz, 
                               *cnc, *fc, *ff, *ratio, resultat);
        
      if (res == 0)
      {
        delete ff; delete fc; delete cnf; delete cnc; 
        delete ratio; if (resr == 2) delete cnr;
        PyErr_SetString(PyExc_ValueError,
                        "integ: integration computation fails.");
        return NULL;
      }
      delete ff; delete fc; delete cnf; delete cnc; delete ratio; 
      if (resr == 2) delete cnr;
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
  
#ifdef E_DOUBLEREAL
  for (E_Int i = 0; i < nFld; i++)
  {
    tpl = Py_BuildValue("d", resultat[i]);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
#else
  for (E_Int i = 0; i < nFld; i++)
  {
    tpl = Py_BuildValue("f", resultat[i]);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
#endif  
  return l;
}

//=============================================================================
// Integre "surfaciquement" les grandeurs de F comme des scalaires
// Retourne 1 si succes, 0 si echec
//=============================================================================
E_Int K_POST::integ1(E_Int niBlk, E_Int njBlk, E_Int nkBlk, 
                     E_Int center2node, E_Int posx, E_Int posy, E_Int posz,
                     FldArrayF& coordBlk, FldArrayF& FBlk, 
                     FldArrayF& ratioBlk, FldArrayF& resultat)
{
  E_Int NI, NJ;
  E_Float resultBlk = 0.;
  E_Int numberOfVariables = FBlk.getNfld();
  if (nkBlk == 1) {NI = niBlk; NJ = njBlk;}
  else if (njBlk == 1) {NI = niBlk; NJ = nkBlk;}
  else if (niBlk == 1) {NI = njBlk; NJ = nkBlk;}
  else return 0;
  
  // Compute surface of each "block" i cell, with coordinates coordBlk
  E_Int ncells = (NI-1)*(NJ-1);
  FldArrayF surfBlk(ncells);
  k6structsurft_(NI, NJ, 1, ncells, coordBlk.begin(posx), coordBlk.begin(posy), coordBlk.begin(posz), surfBlk.begin());
  
  switch (center2node) 
  {
    case 1:
      // Compute integral, coordinates defined in node 
      // and field FBlk in center 
      for (E_Int n = 1; n <= numberOfVariables; n++)
      {      
        k6integstructnodecenter_(NI-1, NJ-1, ratioBlk.begin(), surfBlk.begin(),
                                 FBlk.begin(n), resultBlk);
       
        resultat[n-1] += resultBlk;
      }
      break;
      
    default:
      // Compute integral, coordinates and field have the same size
      for (E_Int n = 1; n <= numberOfVariables; n++)
      {
        k6integstruct_(NI, NJ, ratioBlk.begin(), surfBlk.begin(), 
                       FBlk.begin(n), resultBlk);
       
        resultat[n-1] += resultBlk;
      }
      break;
  }   
  return 1;
}

//=============================================================================
// Integre "lineairement" les grandeurs de F comme des scalaires
// Retourne 1 si succes, 0 si echec
//=============================================================================
E_Int K_POST::integ11D(E_Int niBlk, E_Int njBlk, E_Int nkBlk, 
                       E_Int center2node, E_Int posx, E_Int posy, E_Int posz,
                       FldArrayF& coordBlk, FldArrayF& FBlk, 
                       FldArrayF& ratioBlk, FldArrayF& resultat)
{
  E_Int NI, NJ, NK;
  E_Float resultBlk = 0.;
  E_Int numberOfVariables = FBlk.getNfld();
    
  if (nkBlk == 1 && njBlk == 1)
  {
    NI = niBlk; NJ = njBlk; NK = nkBlk;
  }
  else if (njBlk == 1 && niBlk == 1)
  {
    NI = nkBlk; NJ = njBlk; NK = niBlk;
  }
  else if (niBlk == 1 && nkBlk == 1)
  {
    NI = njBlk; NJ = niBlk; NK = nkBlk;
  }
  else return 0;
  // Compute surface of each "block" i cell, with coordinates coordBlk
  FldArrayF lengthBlk(NI-1);
  k6structsurf1dt_(NI, NJ, NK, coordBlk.begin(posx), coordBlk.begin(posy), coordBlk.begin(posz), lengthBlk.begin());
 
  switch (center2node)
  {
    case 1 :
      for (E_Int n = 1 ; n <= numberOfVariables ; n++)
      {
        // Compute integral, coordinates defined in node 
        // and field FBlk in center 
        k6integstructnodecenter1d_(NI-1, ratioBlk.begin(), lengthBlk.begin(), 
                                   FBlk.begin(n), resultBlk);
        resultat[n-1] = resultat[n-1]+resultBlk;
      }
      break;
    default:
      // Compute integral, coordinates and field have the same size
      for (E_Int n = 1 ; n <= numberOfVariables ; n++)
      {
        k6integstruct1d_(NI, ratioBlk.begin(), lengthBlk.begin(), 
                         FBlk.begin(n), resultBlk);
        resultat[n-1] = resultat[n-1]+resultBlk;
      }
      break;
  }   
  return 1;
}

//=============================================================================
// Integre les grandeurs de F comme des scalaires
// Retourne 1 si succes, 0 si echec
// Attention : cette routine n'integre que sur des elements triangulaires
//=============================================================================
E_Int K_POST::integUnstruct1(E_Int center2node,
                             E_Int posx, E_Int posy, E_Int posz,
                             FldArrayI& cnBlk, FldArrayF& coordBlk, 
                             FldArrayF& FBlk, FldArrayF& ratioBlk, 
                             FldArrayF& resultat)
{
  E_Float resultBlk = 0.;
  E_Int numberOfVariables = FBlk.getNfld();
  E_Int size = coordBlk.getSize();
  E_Int nbT = cnBlk.getSize();
  FldArrayF surfBlk(nbT);
  E_Int nnodes = 3;
  E_Int nedges = 1;
  FldArrayF snx(nbT); // normale a la surface  
  FldArrayF sny(nbT);
  FldArrayF snz(nbT);
  k6unstructsurf_(size, nbT, nedges, nnodes, cnBlk.begin(),
                  coordBlk.begin(posx), coordBlk.begin(posy), 
                  coordBlk.begin(posz),
                  snx.begin(), sny.begin(), snz.begin(), 
                  surfBlk.begin());

  switch (center2node) 
  {
    case 1:
      for (E_Int n = 1; n <= numberOfVariables; n++)
      {
        // Compute integral, coordinates defined in node 
        // and field FBlk in center 
        k6integunsnodecenter_(nbT, ratioBlk.begin(),
                              surfBlk.begin(), 
                              FBlk.begin(n), resultBlk);
        resultat[n-1] = resultat[n-1] + resultBlk;
      }
      break;
    default:
      // Compute integral, coordinates and field have the same size
      for (E_Int n = 1; n <= numberOfVariables; n++)
      {
        k6integunstruct_(nbT, size, cnBlk.begin(), ratioBlk.begin(), 
                         surfBlk.begin(), 
                         FBlk.begin(n), resultBlk);
        
        resultat[n-1] = resultat[n-1] + resultBlk;
      }   
      break;
  }   
  return 1;
}
//=============================================================================
// Integre les grandeurs de F comme des scalaires
// Retourne 1 si succes, 0 si echec
// Attention : cette routine n'integre que sur des elements "bar"
//=============================================================================
E_Int K_POST::integUnstruct11D(E_Int center2node,
                               E_Int posx, E_Int posy, E_Int posz,
                               FldArrayI& cnBlk, FldArrayF& coordBlk, 
                               FldArrayF& FBlk, FldArrayF& ratioBlk, 
                               FldArrayF& resultat)
{
  E_Float resultBlk = 0.;
  E_Int numberOfVariables = FBlk.getNfld();
  E_Int nbT = cnBlk.getSize();
  E_Int size = coordBlk.getSize();
  FldArrayF lengthBlk(nbT);
  
  E_Int nnodes = 2;
  nnodes = 2;
  k6unstructsurf1d_(size, nbT, nnodes, cnBlk.begin(),
                    coordBlk.begin(posx), coordBlk.begin(posy), 
                    coordBlk.begin(posz), lengthBlk.begin());

  switch (center2node) 
  {
    case 1:
      for (E_Int n = 1 ; n <= numberOfVariables ; n++)
      {
        // Compute integral, coordinates defined in node 
        // and field FBlk in center 
        k6integunsnodecenter_(nbT, ratioBlk.begin(),
                              lengthBlk.begin(), 
                              FBlk.begin(n), resultBlk);
        resultat[n-1] = resultat[n-1] + resultBlk;
      }
      break;
    default :
      for (E_Int n = 1 ; n <= numberOfVariables ; n++)    
      {
        // Compute integral, coordinates and field have the same size
        k6integunstruct1d_(nbT, size, cnBlk.begin(), ratioBlk.begin(), 
                           lengthBlk.begin(), 
                           FBlk.begin(n), resultBlk);
        resultat[n-1] = resultat[n-1] + resultBlk;
      }   
      break;
  }
  return 1;
}
  
//=============================================================================
/* Calcul de l'integrale de la solution (pour le pyTree) */
// ============================================================================
PyObject* K_POST::integ2(PyObject* self, PyObject* args)
{
  PyObject* zone; char* varName;
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  if (!PyArg_ParseTuple(args, "Ossss", &zone, &varName, &GridCoordinates, 
                        &FlowSolutionNodes, &FlowSolutionCenters)) return NULL;

  E_Int ni, nj, nk, cnSize, cnNfld;
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  vector<PyArrayObject*> hook;
  E_Int res = K_PYTREE::getFromZone(zone, 0, 1, varString, 
                                    fields, locs, ni, nj, nk,
                                    cn, cnSize, cnNfld,
                                    eltType, hook,
                                    GridCoordinates, 
                                    FlowSolutionNodes, FlowSolutionCenters);

  // Get the pointers on fields
  E_Int posVol = K_ARRAY::isNamePresent("vol", varString);
  E_Int posRatio = K_ARRAY::isNamePresent("ratio", varString);
  E_Int posVar = K_ARRAY::isNamePresent(varName, varString); 

  if (posVar == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "integ2: variable doesn't exist in array.");
    RELEASESHAREDZ(hook, varString, eltType);
    return NULL;
  }

  E_Float* fp = fields[posVar];
  E_Float* v = fields[posVol];
  
  E_Int n;
  if (res == 1) n = max(ni-1,E_Int(1))*max(nj-1,E_Int(1))*max(nk-1,E_Int(1));
  else n = nj;
  //printf("%d\n", n);

  E_Float ret = 0.;
  if (posRatio == -1)
  {
#pragma omp parallel for reduction(+:ret)
    for (E_Int i = 0; i < n; i++) ret += fp[i]*v[i];
  }
  else
  {
    E_Float* r = fields[posRatio];
#pragma omp parallel for reduction(+:ret)
    for (E_Int i = 0; i < n; i++) ret += fp[i]*v[i]*r[i];
  }
  RELEASESHAREDZ(hook, varString, eltType);

  return Py_BuildValue("d", ret);
}
