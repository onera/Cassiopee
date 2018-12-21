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

// Integration

# include <string.h>
# include <stdio.h>
# include "post.h"
# include <vector>

using namespace std;
using namespace K_FLD;

extern "C"
{
  void k6normstructsurft_(const E_Int& ni, const E_Int& nj, const E_Int& npts, 
                          const E_Float* xt, const E_Float* yt, const E_Float* zt,
                          E_Float* nxt, E_Float* nyt, E_Float* nzt);

  void k6normunstructsurf_(const E_Int& nt, const E_Int& nv,
                           E_Int* cn, 
                           E_Float* coordx, E_Float* coordy, E_Float* coordz,
                           E_Float* surf);

  void k6integnormprodstruct_(const E_Int& ni, const E_Int& nj, E_Float* ratio,
                              E_Float* sx, E_Float* sy, E_Float* sz, 
                              E_Float* vx, E_Float* vy, E_Float* vz, 
                              E_Float& result);

  void k6integnormprodstructnodecenter_(
    const E_Int& ni, const E_Int& nj, E_Float* ratio,
    E_Float* sx, E_Float* sy, E_Float* sz, 
    E_Float* vx, E_Float* vy, E_Float* vz, 
    E_Float& result);

  void k6integnormprodunstruct_(const E_Int& nbt, const E_Int& size, 
                                E_Int* cn, E_Float* ratio, 
                                E_Float* sx, E_Float* sy, E_Float* sz, 
                                E_Float* vx, E_Float* vy, E_Float* vz, 
                                E_Float& result);

  void k6integnormprodunsnodecenter_(const E_Int& nbt, E_Float* ratio, 
                                     E_Float* sx, E_Float* sy, E_Float* sz, 
                                     E_Float* vx, E_Float* vy, E_Float* vz, 
                                     E_Float& result);
}
//============================================================================
/* Calcul une integrale du produit scalaire 
   de la solution*normale (vect(F).vect(n)) */
// ===========================================================================
PyObject* K_POST::integNormProduct(PyObject* self, PyObject* args)
{
  PyObject* coordArrays;
  PyObject* FArrays;
  PyObject* ratioArrays;
  if (!PyArg_ParseTuple(args, "OOO", &coordArrays, &FArrays, &ratioArrays))
  {
    return NULL;
  }

  // Check every array 
  if (PyList_Check(coordArrays) == 0)
  {
    PyErr_SetString(
      PyExc_TypeError, 
      "integNormProduct: first argument must be a list.");
    return NULL;
  }
  if (PyList_Check(FArrays) == 0)
  {
    PyErr_SetString(
      PyExc_TypeError, 
      "integNormProduct: second argument must be a list.");
    return NULL;
  }
  if (PyList_Check(ratioArrays) == 0)
  {
    PyErr_SetString(
      PyExc_TypeError, 
      "integNormProduct : third argument must be a list.");
    return NULL;
  }
  E_Int nCoordArrays = PyList_Size(coordArrays);
  E_Int nFArrays = PyList_Size(FArrays);
  E_Int nRatioArrays = PyList_Size(ratioArrays);

  if (nRatioArrays == 0) 
  { 
    if (nCoordArrays != nFArrays)
    {
      PyErr_SetString(
        PyExc_ValueError, 
        "integNormProduct: number of zones in 1st and 2nd arguments must be equal.");
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
                      "integNormProduct: number of zones in 1st, 2nd and 3rd arguments must be equal.");
      return NULL;
    }
  }

  PyObject* coordObj; PyObject* FObj; PyObject* ratioObj;
  E_Int nic, njc, nkc, nif, njf, nkf, nir, njr, nkr;
  char* varStringc; char* varStringf; char* varStringr;
  char varString0[K_ARRAY::VARSTRINGLENGTH];
  FldArrayF* fc; FldArrayF* ff; FldArrayF* ratio;
  E_Int nFld = -1;
  E_Int sizef = 0;
  E_Int center2node = 2; // set to 1 if coord is in nodes and F in centers
                         // set to 0 if coord and F have the same size

  char* eltTypec; FldArrayI* cnc;
  char* eltTypef; FldArrayI* cnf;
  char* eltTyper; FldArrayI* cnr;
  E_Int resc = 0;
  E_Int resf = 0;
  E_Int resr = 0;
  E_Float resultat = 0.;
  E_Int res = -1;

  for (int i = 0; i < nCoordArrays; i++)
  {
    coordObj = PyList_GetItem(coordArrays,i);
    FObj = PyList_GetItem(FArrays,i);
    resc = K_ARRAY::getFromArray(coordObj, varStringc, fc, 
                                 nic, njc, nkc, cnc, eltTypec); 
    
    if (resc != 1 && resc != 2)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "integNormProduct: coord is not a valid array.");
      return NULL;
    }
    
    E_Int posx = K_ARRAY::isCoordinateXPresent(varStringc);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varStringc);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varStringc);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete fc; if (resc == 2) delete cnc;
      printf("Warning: integNormProduct: coordinates not found in array %d.", i+1);
      printf(" Array skipped...\n");
      goto next; 
    }    
    posx++; posy++; posz++;

    resf = K_ARRAY::getFromArray(FObj, varStringf, ff, nif, njf, nkf, 
                                 cnf, eltTypef); 
    if (resf != 1 && resf != 2)
    {
      delete ff; delete fc;
      if (resc == 2) delete cnc;
      PyErr_SetString(PyExc_TypeError, 
                      "integNormProduct: field is not a valid array.");
      return NULL;
    }
    // check number of variables
    if ( nFld == -1 ) // premier passage
    {
      if (ff->getNfld() != 3)
      {
        printf("Warning: integNormProduct: vector array %d must be 3D.\n", i+1);
        delete ff; if (resf == 2) delete cnf;
        delete fc; if (resc == 2) delete cnc;
        goto next;
      }
      strcpy(varString0, varStringf);
      nFld = ff->getNfld();
    }
    else 
    {
      if (ff->getNfld() != nFld)
      {
        printf("Warning: integNormProduct: invalid number of variables for field array %d.", i+1);
        printf(" Array skipped...\n");
        delete ff; if (resf == 2) delete cnf;
        delete fc; if (resc == 2) delete cnc;
        goto next;
      }

      // check is variables are ordered in the same way
      E_Int ids = K_ARRAY::compareVarStrings(varString0, varStringf);
      if (ids == -1) // varstrings are different
      {
        printf("Warning: integNormProduct: variables are in a different order to the first array. Array skipped...\n");
        delete ff;if (resf == 2) delete cnf;
        delete fc; if (resc == 2) delete cnc;
        goto next;
      }
    }

    // cas structure
    if (resc == 1 && resf == 1)
    {
      sizef = nif*njf*nkf;
      if ( (nic == 1 && njc == 1) ||
           (nic == 1 && nkc == 1) ||
           (njc == 1 && nkc == 1) ||
           (nic > 1 && njc > 1 && nkc > 1))
      {
        printf("Warning: integNormProduct: coordinates array must be 2D. Array skipped...\n");
        delete ff; if (resf == 2) delete cnf;
        delete fc; if (resc == 2) delete cnc;
        goto next;
      }
      
      // coord et F de meme taille ?
      if (nic == nif && njc == njf && nkc == nkf)
        center2node = 0;
      
      else if (nic == nif+1 && njc == njf+1 && nkc == nkf+1)
        center2node = 1;
      
      else
      {
        if ( (nic == 1 && njc == njf+1 && nkc == nkf+1) || 
             (njc == 1 && nic == nif+1 && nkc == nkf+1) || 
             (nkc == 1 && nic == nif+1 && njc == njf+1) )
            center2node = 1;
        else if ( (nic == 1 && njc == 1 && nkc == nkf+1) || 
                  (nic == nif+1 && njc == 1 && nkc == 1) || 
                  (nic == 1 && njc == njf+1 && nkc == 1) )
          center2node = 1;
        else 
        {
          printf("Warning: integNormProduct: coord and f arrays do not represent the same zone. Array skipped...\n");
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
        resr = K_ARRAY::getFromArray(
          ratioObj, varStringr, ratio, nir, njr, nkr, cnr, eltTyper); 
        if (resr != 1)
        {
          if ( resr == 2 ) {delete ratio; delete cnr;}
          printf("Warning: integNormProduct: ratio %d is an invalid array. Set to 1.", i+1);
          ratio = new FldArrayF(sizef);
          ratio->setAllValuesAt(1.);
        }
      }
      // integ sur chaque bloc
      res = 0;
      res = integ3(nic, njc, nkc, center2node, posx, posy, posz,
                   *fc, *ff, *ratio, resultat);
      if ( res == 0 ) 
      {
        delete ff; delete fc;delete ratio;
        PyErr_SetString(PyExc_ValueError,
                        "integNormProduct: integration computation fails.");
        return NULL;
      }
      delete ff; delete fc; delete ratio;
    }
    else if (resc == 2 && resf == 2)//cas non structure
    {
      if ( strcmp(eltTypec,"TRI")==0 && strcmp(eltTypef,"TRI")==0) 
        center2node = 0;
      else if (strcmp(eltTypec,"TRI")==0 && strcmp(eltTypef,"TRI*")==0 )
        center2node = 1;
      else
      {
        delete fc; delete cnc;
        delete ff; delete cnf;
        PyErr_SetString(PyExc_ValueError, 
                        "integNormProduct: only TRI unstructured arrays are possible.");
        return NULL;
      }

      sizef = ff->getSize();
      
      // Valeur du ratio ?
      if (nRatioArrays == 0)
      {
        ratio = new FldArrayF(sizef);
        ratio->setAllValuesAt(1.);
      }
      else //coord + F + r
      {
        ratioObj = PyList_GetItem(ratioArrays, i);
        resr = K_ARRAY::getFromArray(
          ratioObj, varStringr, ratio, nir, njr, nkr, cnr, eltTyper); 
        if (resr != 2)
        {
          if ( resr == 1 ) delete ratio;
          printf("Warning: integNormProduct: ratio %d is an invalid array. Set to 1.", i+1);
          ratio = new FldArrayF(sizef);
          ratio->setAllValuesAt(1.);
        }
      }
      // integ sur chaque bloc
      res = 0;
      
      res = integUnstruct3(center2node, posx, posy, posz, 
                           *cnc, *fc, *ff, *ratio, resultat);
      if ( res == 0 ) 
      {
        delete ff; delete fc; delete cnf; delete cnc;
        delete ratio; 
        if (resr == 2) delete cnr;
        PyErr_SetString(PyExc_ValueError,
                        "integNormProduct: integration computation fails.");
        return NULL;
      }
      delete ff; delete fc; delete cnf; delete cnc; delete ratio; 
      if (resr == 2) delete cnr;
    }
    else 
    {
      PyErr_SetString(PyExc_TypeError, 
                      "integNormProduct: invalid arrays.");
      return NULL;
    }
    next:;
  }
   if ( res == -1 )
  {
    PyErr_SetString(PyExc_TypeError,
                    "integNormProduct: integration failed.");
    return NULL;
  }

  PyObject* tpl;
#ifdef E_DOUBLEREAL
  tpl = Py_BuildValue("d", resultat);
#else
  tpl = Py_BuildValue("f", resultat);
#endif
  return tpl;
}

//=============================================================================
// Integre les grandeurs de vect(F).vect(n)
// Retourne 1 si success, 0 si echec
//=============================================================================
E_Int K_POST::integ3(E_Int niBlk, E_Int njBlk, E_Int nkBlk, 
                     E_Int center2node,
                     E_Int posx, E_Int posy, E_Int posz,
                     FldArrayF& coordBlk, FldArrayF& FBlk, 
                     FldArrayF& ratioBlk, E_Float& resultat)
{
  E_Int NI, NJ;
  E_Float resultBlk = 0.;
  
  if (nkBlk == 1)
  { NI = niBlk; NJ = njBlk; }
  else if (njBlk == 1)
  { NI = niBlk; NJ = nkBlk; }
  else if (niBlk == 1)
  { NI = njBlk; NJ = nkBlk; }
  else return 0;
 
  // Compute surface of each "block" i cell, with coordinates coordBlk
  FldArrayF nsurfBlk((NI-1)*(NJ-1),3);

  E_Int npts = coordBlk.getSize();
  k6normstructsurft_(NI, NJ, npts, coordBlk.begin(posx), coordBlk.begin(posy), coordBlk.begin(posz), 
                     nsurfBlk.begin(1), nsurfBlk.begin(2), nsurfBlk.begin(3));

  if (center2node == 1)
  {
    // Compute integral, coordinates defined in node 
    // and field FBlk in center 
    k6integnormprodstructnodecenter_(NI-1, NJ-1, ratioBlk.begin(), 
                                     nsurfBlk.begin(1), nsurfBlk.begin(2), 
                                     nsurfBlk.begin(3), FBlk.begin(1), 
                                     FBlk.begin(2), FBlk.begin(3),
                                     resultBlk);
  }
  else
  {
    // Compute integral, coordinates and field have the same size
    k6integnormprodstruct_(NI, NJ, ratioBlk.begin(), nsurfBlk.begin(1), 
                           nsurfBlk.begin(2),  nsurfBlk.begin(3),  
                           FBlk.begin(1), FBlk.begin(2), FBlk.begin(3),
                           resultBlk);
  }
  resultat = resultat + resultBlk;   
  return 1;
}
  
//=============================================================================
// Integre les grandeurs de vect(F).vect(n)
// Retourne 1 si success, 0 si echec
//=============================================================================
E_Int K_POST::integUnstruct3(E_Int center2node,
                             E_Int posx, E_Int posy, E_Int posz,
                             FldArrayI& cnBlk, FldArrayF& coordBlk, 
                             FldArrayF& FBlk, FldArrayF& ratioBlk, 
                             E_Float& resultat)
{
  E_Float resultBlk = 0.;
  E_Int size = coordBlk.getSize();
  E_Int nbT = cnBlk.getSize();
  
  // Compute surface of each "block" i cell, with coordinates coordBlk
  FldArrayF nsurfBlk(nbT, 3);
  k6normunstructsurf_(nbT, size, cnBlk.begin(), coordBlk.begin(posx), 
                      coordBlk.begin(posy), coordBlk.begin(posz),
                      nsurfBlk.begin());  
  if (center2node == 1)
  {
    // Compute integral, coordinates defined in node 
    // and field FBlk in center 
    k6integnormprodunsnodecenter_(nbT, ratioBlk.begin(), 
                                  nsurfBlk.begin(1), nsurfBlk.begin(2),
                                  nsurfBlk.begin(3), FBlk.begin(1),
                                  FBlk.begin(2), FBlk.begin(3),
                                  resultBlk);
  }
  else
  {
    // Compute integral, coordinates and field have the same size
    k6integnormprodunstruct_(nbT, size, cnBlk.begin(), ratioBlk.begin(), 
                             nsurfBlk.begin(1), nsurfBlk.begin(2),
                             nsurfBlk.begin(3), FBlk.begin(1),
                             FBlk.begin(2), FBlk.begin(3),
                             resultBlk);
  }
  resultat = resultat+resultBlk;
  return 1;
}
