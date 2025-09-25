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

// Integration

# include <string.h>
# include <stdio.h>
# include "post.h"
# include <vector>

using namespace std;
using namespace K_FLD;

//============================================================================
/* Calcul une integrale du produit scalaire 
   de la solution*normale (vect(F).vect(n)) */
// ===========================================================================
PyObject* K_POST::integNormProduct(PyObject* self, PyObject* args)
{
  PyObject* coordArrays;
  PyObject* FArrays;
  PyObject* ratioArrays;
  if (!PYPARSETUPLE_(args, OOO_, &coordArrays, &FArrays, &ratioArrays))
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
      "integNormProduct: third argument must be a list.");
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

  for (E_Int i = 0; i < nCoordArrays; i++)
  {
    coordObj = PyList_GetItem(coordArrays,i);
    FObj = PyList_GetItem(FArrays,i);
    resc = K_ARRAY::getFromArray3(coordObj, varStringc, fc, 
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
      RELEASESHAREDB(resc, coordObj, fc, cnc);
      printf("Warning: integNormProduct: coordinates not found in array %d.", i+1);
      printf(" Array skipped...\n");
      goto next; 
    }    
    posx++; posy++; posz++;

    resf = K_ARRAY::getFromArray3(FObj, varStringf, ff, nif, njf, nkf, 
                                  cnf, eltTypef); 
    if (resf != 1 && resf != 2)
    {
      RELEASESHAREDS(FObj, ff);
      RELEASESHAREDB(resc, coordObj, fc, cnc);
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
        RELEASESHAREDB(resc, coordObj, fc, cnc);
        RELEASESHAREDB(resf, FObj, ff, cnf);
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
        RELEASESHAREDB(resc, coordObj, fc, cnc);
        RELEASESHAREDB(resf, FObj, ff, cnf);
        goto next;
      }

      // check is variables are ordered in the same way
      E_Int ids = K_ARRAY::compareVarStrings(varString0, varStringf);
      if (ids == -1) // varstrings are different
      {
        printf("Warning: integNormProduct: variables are in a different order to the first array. Array skipped...\n");
        RELEASESHAREDB(resc, coordObj, fc, cnc);
        RELEASESHAREDB(resf, FObj, ff, cnf);
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
        RELEASESHAREDB(resc, coordObj, fc, cnc);
        RELEASESHAREDB(resf, FObj, ff, cnf);
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
          RELEASESHAREDS(FObj, ff);
          RELEASESHAREDS(coordObj, ff);
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
        resr = K_ARRAY::getFromArray3(
          ratioObj, varStringr, ratio, nir, njr, nkr, cnr, eltTyper); 
        if (resr != 1)
        {
          RELEASESHAREDB(resr, ratioObj, ratio, cnr);
          printf("Warning: integNormProduct: ratio %d is an invalid array. Set to 1.", i+1);
          ratio = new FldArrayF(sizef);
          ratio->setAllValuesAt(1.);
        }
      }
      // integ sur chaque bloc
      res = 0;
      res = integNormProdStruct2D(nic, njc, nkc, center2node, posx, posy, posz,
                   *fc, *ff, *ratio, resultat);

      if (res == 0) 
      {
        RELEASESHAREDS(coordObj, fc);
        RELEASESHAREDS(FObj, ff);
        if (nRatioArrays == 0 || resr != 1) delete ratio;
        else RELEASESHAREDS(ratioObj, ratio);
        PyErr_SetString(PyExc_ValueError,
                        "integNormProduct: integration computation fails.");
        return NULL;
      }
      RELEASESHAREDS(coordObj, fc);
      RELEASESHAREDS(FObj, ff);
      if (nRatioArrays == 0 || resr != 1) delete ratio;
      else RELEASESHAREDS(ratioObj, ratio);
    }
    else if (resc == 2 && resf == 2)//cas non structure
    {
      // check if field is cell or node centered
      E_Int l = strlen(eltTypef);
      if (eltTypef[l-1] == '*') center2node = 1;
      else center2node = 0;

      res = 1;
      std::vector<char*> eltTypecs, eltTypefs;
      K_ARRAY::extractVars(eltTypec, eltTypecs);

      // check if elt is valid (QUAD, TRI)
      for (E_Int ic = 0; ic < eltTypecs.size(); ic++)
      {
        if ((strcmp(eltTypecs[ic], "QUAD") != 0) && (strcmp(eltTypecs[ic], "TRI") != 0)) res = 0;
      }

      for (size_t ic = 0; ic < eltTypecs.size(); ic++) delete [] eltTypecs[ic];

      if (res == 0)
      {
        RELEASESHAREDU(coordObj, fc, cnc);
        RELEASESHAREDU(FObj, ff, cnf);
        PyErr_SetString(PyExc_ValueError, 
                        "integNormProd: only QUAD or TRI unstructured arrays are possible.");
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
        resr = K_ARRAY::getFromArray3(
          ratioObj, varStringr, ratio, nir, njr, nkr, cnr, eltTyper); 
        if (resr != 2)
        {
          if (resr == 1) RELEASESHAREDS(ratioObj, ratio);
          printf("Warning: integNormProduct: ratio %d is an invalid array. Set to 1.", i+1);
          ratio = new FldArrayF(sizef);
          ratio->setAllValuesAt(1.);
        }
      }
      // integ sur chaque bloc
      res = 0;
      res = integNormProdUnstruct2D(center2node, posx, posy, posz, 
                           *cnc, eltTypec, *fc, *ff, *ratio, resultat);

      if (res == 0) 
      {
        RELEASESHAREDU(coordObj, fc, cnc);
        RELEASESHAREDU(FObj, ff, cnf);
        if (nRatioArrays == 0 || resr != 2) delete ratio;
        else RELEASESHAREDB(resr, ratioObj, ratio, cnr);
        PyErr_SetString(PyExc_ValueError,
                        "integNormProduct: integration computation fails.");
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
                      "integNormProduct: invalid arrays.");
      return NULL;
    }
    next:;
  }
   
  if (res == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "integNormProduct: integration failed.");
    return NULL;
  }

  PyObject* tpl = Py_BuildValue(R_, resultat);
  return tpl;
}
