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

// Moment integration

# include <stdio.h>
# include <string.h>
# include "post.h"
# include <vector>

using namespace std;
using namespace K_FLD;

extern "C"
{
  void k6integmomentstruct_(const E_Int& ni, const E_Int& nj,
                            const E_Float& cx, const E_Float& cy, 
                            const E_Float& cz, E_Float* ratio, 
                            E_Float* xt, E_Float* yt, E_Float* zt, 
                            E_Float* surf, E_Float* vx, E_Float* vy, 
                            E_Float* vz, E_Float* result);

  void k6integmomentstruct1d_(const E_Int& ni, const E_Float& cx, 
                              const E_Float& cy, const E_Float& cz, 
                              E_Float* ratio, E_Float* xt, E_Float* yt, 
                              E_Float* zt, E_Float* length, E_Float* vx, 
                              E_Float* vy, E_Float* vz, E_Float* result);

  void k6integmomentstructnodecenter_(const E_Int& ni, const E_Int& nj,
                                      const E_Float& cx, const E_Float& cy, 
                                      const E_Float& cz, E_Float* ratio, 
                                      E_Float* xt, E_Float* yt, E_Float* zt, 
                                      E_Float* surf, E_Float* vx, 
                                      E_Float* vy, E_Float* vz,
                                      E_Float* result);

  void k6integmomentstructnodecenter1d_(const E_Int& ni, const E_Float& cx, 
                                        const E_Float& cy, const E_Float& cz,
                                        E_Float* ratio, E_Float* xt, 
                                        E_Float* yt, E_Float* zt,
                                        E_Float* length, E_Float* vx, 
                                        E_Float* vy, E_Float* vz, 
                                        E_Float* result);
  
  void k6integmomentunstruct_(const E_Int& nbt, const E_Int& size, E_Int* cn,
                              const E_Float& cx, const E_Float& cy, 
                              const E_Float& cz, E_Float* ratio, 
                              E_Float* xt, E_Float* yt, E_Float* zt,
                              E_Float* surf, E_Float* vx, E_Float* vy, 
                              E_Float* vz, E_Float* result);

  void k6integmomentunstruct1d_(const E_Int& nbt, const E_Int& size, 
                                E_Int* cn, const E_Float& cx, 
                                const E_Float& cy, const E_Float& cz,
                                E_Float* ratio, E_Float* xt, E_Float* yt, 
                                E_Float* zt, E_Float* length,
                                E_Float* vx, E_Float* vy, 
                                E_Float* vz, E_Float* result);

  void k6integmomentunsnodecenter_(const E_Int& nbt, const E_Int& size, 
                                   E_Int* cn, const E_Float& cx, 
                                   const E_Float& cy, const E_Float& cz,
                                   E_Float* ratio, 
                                   E_Float* xt, E_Float* yt, E_Float* zt,
                                   E_Float* surf, E_Float* vx, E_Float* vy, 
                                   E_Float* vz, E_Float* result);

  void k6integmomentunsnodecenter1d_(const E_Int& nbt, const E_Int& size, 
                                     E_Int* cn, const E_Float& cx, 
                                     const E_Float& cy, const E_Float& cz,
                                     E_Float* ratio, 
                                     E_Float* xt, E_Float* yt, E_Float* zt,
                                     E_Float* surf,
                                     E_Float* vx, E_Float* vy, E_Float* vz,
                                     E_Float* result);
}
//=============================================================================
/* Calcul une integrale du moment (OM^F) */
// ============================================================================
PyObject* K_POST::integMoment(PyObject* self, PyObject* args)
{
  E_Float cx, cy, cz;
  PyObject* coordArrays;
  PyObject* FArrays;
  PyObject* ratioArrays;
  if (!PYPARSETUPLE_(args, OOO_ TRRR_,
                    &coordArrays, &FArrays, &ratioArrays, &cx, &cy, &cz))
  {
    return NULL;
  }
  
  // Check every array in listFields
  if (PyList_Check(coordArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "integMoment: first argument must be a list.");
    return NULL;
  }
  if (PyList_Check(FArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "integMoment: second argument must be a list.");
    return NULL;
  }
  if (PyList_Check(ratioArrays) == 0)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "integMoment: third argument must be a list.");
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
        "integMoment: number of zones in 1st and 2nd arguments must be equal.");
      return NULL;  
    }
  }
  else
  {
    if (nCoordArrays != nFArrays || 
        nCoordArrays != nRatioArrays ||
        nRatioArrays != nFArrays)
    {
      PyErr_SetString(
        PyExc_ValueError, 
        "integMoment: number of zones in 1st, 2nd and 3rd arguments must be equal.");
      return NULL;    
    }
  }

  PyObject* coordObj; PyObject* FObj; PyObject* ratioObj;
  E_Int nic, njc, nkc, nif, njf, nkf, nir, njr, nkr;
  char* varStringc; char* varStringf; char* varStringr;
  char varString0[K_ARRAY::VARSTRINGLENGTH];
  FldArrayF* fc;
  FldArrayF* ff;
  FldArrayF* ratio;
  E_Int sizef = 0;
  E_Int nFld = -1;
  E_Int center2node = 2; // set to 1 if coord is in nodes and F in centers
                         // set to 0 if coord and F have the same size
  E_Int case1D;          // set to 1 if linear integration


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
    coordObj = PyList_GetItem(coordArrays,i);
    FObj = PyList_GetItem(FArrays,i);
    resc = K_ARRAY::getFromArray3(coordObj, varStringc, 
                                  fc, nic, njc, nkc, cnc, eltTypec); 
    
    if (resc != 1 && resc != 2)
    {
      PyErr_SetString(PyExc_TypeError,
                      "integMoment: coord is not a valid array.");
      return NULL;
    }
    
    E_Int posx = K_ARRAY::isCoordinateXPresent(varStringc);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varStringc);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varStringc);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      printf("Warning: integMoment: coordinates not found in array %d. Array skipped...\n",i+1);
      RELEASESHAREDB(resc, coordObj, fc, cnc);
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
                      "integMoment: field is not a valid array.");
      return NULL;
    }
    
    // check number of variables
    if (nFld == -1) // premier passage
    {
      if (ff->getNfld() != 3)
      {
        printf("Warning: integMoment: field must be a vector.\n");
        RELEASESHAREDB(resc, coordObj, fc, cnc);
        RELEASESHAREDB(resf, FObj, ff, cnf);
        goto next;
      }
      strcpy(varString0, varStringf);
      nFld = ff->getNfld();
      resultat.malloc(nFld);
      resultat.setAllValuesAtNull();
    }
    else 
    {
      if (ff->getNfld() != nFld)
      {
        printf("Warning: integMoment: invalid number of variables for field array %d.", i+1);
        printf(" Array skipped...\n");
        RELEASESHAREDB(resc, coordObj, fc, cnc);
        RELEASESHAREDB(resf, FObj, ff, cnf);
        goto next;
      }
      // check is variables are ordered in the same way
      E_Int ids = K_ARRAY::compareVarStrings(varString0, varStringf);
      if (ids == -1) // varstrings are different
      {
        printf("Warning: integMoment: variables are in a different order to the first array. Array skipped...\n");
        RELEASESHAREDB(resc, coordObj, fc, cnc);
        RELEASESHAREDB(resf, FObj, ff, cnf);
        goto next;
      }
    }

    // cas structure
    if (resc == 1 && resf == 1)
    {
      sizef = nif*njf*nkf;
      if (nic > 1 && njc > 1 && nkc > 1)
      {
        printf("Warning: integMoment: 3D arrays not valid. Array %d skipped...\n",i+1);
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
          printf("Warning: integMoment: coord and field arrays do not represent the same zone.");
          printf("Array skipped...\n");
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
        resr = K_ARRAY::getFromArray3(ratioObj, varStringr, ratio, 
                                      nir, njr, nkr, cnr, eltTyper); 
        if (resr != 1)
        {
          RELEASESHAREDB(resr, ratioObj, ratio, cnr);
          printf("Warning: integMoment: ratio %d is an invalid array. Set to 1.", i+1);
          ratio = new FldArrayF(sizef);
          ratio->setAllValuesAt(1.);
        }
      }
      // integ sur chaque bloc
      res = 0;
      switch (case1D) 
      {
        case 1:
          res = integ41D(nic, njc, nkc, center2node, 
                         posx, posy, posz, cx, cy, cz, 
                         *fc, *ff, *ratio, resultat);
          break;
        default:
          res = integ4(nic, njc, nkc, center2node, 
                       posx, posy, posz, cx, cy, cz,
                       *fc, *ff, *ratio, resultat);
          break;
      }
      if (res == 0) 
      {
        RELEASESHAREDS(coordObj, fc);
        RELEASESHAREDS(FObj, ff);
        if (nRatioArrays == 0 || resr != 1) delete ratio;
        else RELEASESHAREDS(ratioObj, ratio);
        PyErr_SetString(PyExc_ValueError,
                        "integMoment: integration computation fails.");
        return NULL;
      }    
      RELEASESHAREDS(coordObj, fc);
      RELEASESHAREDS(FObj, ff);
      if (nRatioArrays == 0 || resr != 1) delete ratio;
      else RELEASESHAREDS(ratioObj, ratio);
    }
    
    else if (resc == 2 && resf == 2) // cas non structure
    {
      if ((strcmp(eltTypec, "TRI") == 0 && strcmp(eltTypef, "TRI") == 0) || 
          (strcmp(eltTypec, "BAR") == 0 && strcmp(eltTypef, "BAR") == 0))
        center2node = 0;
      else if ((strcmp(eltTypec, "TRI") == 0 && strcmp(eltTypef, "TRI*") == 0) || 
               (strcmp(eltTypec, "BAR") == 0 && strcmp(eltTypef, "BAR*") == 0))
        center2node = 1;
      else
      {
        RELEASESHAREDU(coordObj, fc, cnc);
        RELEASESHAREDU(FObj, ff, cnf);        
        PyErr_SetString(PyExc_ValueError, 
                        "integMoment: only TRI or BAR unstructured arrays are possible.");
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
          printf("Warning: integMoment: ratio %d is an invalid array. Set to 1.", i+1);
          ratio = new FldArrayF(sizef);
          ratio->setAllValuesAt(1.);
        }
      }

      // integ sur chaque bloc
      res = 0;
      if (strcmp(eltTypec, "TRI") == 0)
        res = integUnstruct4(center2node, posx, posy, posz, cx, cy, cz, 
                             *cnc, *fc, *ff, *ratio, resultat); 
      else
        res = integUnstruct41D(center2node, posx, posy, posz, cx, cy, cz, 
                               *cnc, *fc, *ff, *ratio, resultat);
        
      if (res == 0)
      {
        RELEASESHAREDU(coordObj, fc, cnc);
        RELEASESHAREDU(FObj, ff, cnf);
        if (nRatioArrays == 0 || resr != 2) delete ratio;
        else RELEASESHAREDB(resr, ratioObj, ratio, cnr);
        PyErr_SetString(PyExc_ValueError,
                        "integMoment: integration computation fails.");
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
                      "integMoment: not valid arrays.");
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
#ifdef E_DOUBLEREAL  
    tpl = Py_BuildValue("d", resultat[i]);
#else
    tpl = Py_BuildValue("f", resultat[i]);
#endif
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }
  return l;
}

//=============================================================================
// Integre "surfaciquement" les grandeurs de M = OM^F
// Retourne 1 si success, 0 si echec
//=============================================================================
E_Int K_POST::integ4(E_Int niBlk, E_Int njBlk, E_Int nkBlk, 
                     E_Int center2node, 
                     E_Int posx, E_Int posy, E_Int posz,
                     E_Float cx, E_Float cy, E_Float cz, 
                     FldArrayF& coordBlk, 
                     FldArrayF& FBlk, FldArrayF& ratioBlk, 
                     FldArrayF& resultat)
{
  E_Int NI, NJ;
  FldArrayF res(3);
  
  if (nkBlk == 1)
  { NI = niBlk; NJ = njBlk; }
  else if (njBlk == 1)
  { NI = niBlk; NJ = nkBlk; }
  else if (niBlk == 1)
  { NI = njBlk; NJ = nkBlk; }
  else return 0;
  
  E_Int ncells = (NI-1)*(NJ-1);
  FldArrayF surfBlk(ncells);
  
  // Compute surface of each "block" i cell, with coordinates coordBlk
  K_METRIC::compStructSurft(
    NI, NJ, 1,
    coordBlk.begin(posx), coordBlk.begin(posy), coordBlk.begin(posz),
    surfBlk.begin());

  if (center2node == 1)
  {
    // Compute integral, coordinates defined in node 
    // and field FBlk in center 
    k6integmomentstructnodecenter_(NI, NJ, cx, cy, cz, ratioBlk.begin(), 
                                   coordBlk.begin(posx), coordBlk.begin(posy),
                                   coordBlk.begin(posz), surfBlk.begin(), 
                                   FBlk.begin(1),FBlk.begin(2),FBlk.begin(3), 
                                   res.begin());
  }
  else
  {
    // Compute integral, coordinates and field have the same size
    k6integmomentstruct_(NI, NJ, cx, cy, cz, ratioBlk.begin(), 
                         coordBlk.begin(posx), coordBlk.begin(posy),
                         coordBlk.begin(posz), surfBlk.begin(), 
                         FBlk.begin(1), FBlk.begin(2), FBlk.begin(3), 
                         res.begin());
  }
  
  resultat[0] += res[0];
  resultat[1] += res[1];
  resultat[2] += res[2];   
  return 1;
}
//=============================================================================
// Integre "lineairement" les grandeurs de M = OM^F
// Retourne 1 si success, 0 si echec
//=============================================================================
E_Int K_POST::integ41D(E_Int niBlk, E_Int njBlk, E_Int nkBlk, 
                       E_Int center2node, 
                       E_Int posx, E_Int posy, E_Int posz,
                       E_Float cx, E_Float cy, E_Float cz, 
                       FldArrayF& coordBlk, 
                       FldArrayF& FBlk, FldArrayF& ratioBlk, 
                       FldArrayF& resultat)
{
  E_Int NI, NJ, NK;
  FldArrayF res(3);
  resultat.setAllValuesAtNull();

  if (nkBlk == 1 && njBlk == 1)
  { NI = niBlk; NJ = njBlk; NK = nkBlk; }
  else if (njBlk == 1 && niBlk == 1)
  { NI = nkBlk; NJ = njBlk; NK = niBlk; }
  else if (niBlk == 1 && nkBlk == 1)
  { NI = njBlk; NJ = niBlk; NK = nkBlk; }
  else return 0;
  
  FldArrayF lengthBlk(NI-1);
  // Compute surface of each "block" i cell, with coordinates coordBlk
  K_METRIC::compStructSurf1dt(
    NI, NJ, NK,
    coordBlk.begin(posx), coordBlk.begin(posy), coordBlk.begin(posz),
    lengthBlk.begin());
  
  if (center2node == 1)
  {
    // Compute integral, coordinates defined in node 
    // and field FBlk in center 
    k6integmomentstructnodecenter1d_(NI, cx, cy, cz, ratioBlk.begin(), 
                                     coordBlk.begin(posx), 
                                     coordBlk.begin(posy), 
                                     coordBlk.begin(posz), 
                                     lengthBlk.begin(), FBlk.begin(1), 
                                     FBlk.begin(2), FBlk.begin(3), 
                                     res.begin());
  }
  else
  {
    // Compute integral, coordinates and field have the same size
    k6integmomentstruct1d_(NI, cx, cy, cz, ratioBlk.begin(),
                           coordBlk.begin(posx), coordBlk.begin(posy),
                           coordBlk.begin(posz), lengthBlk.begin(),
                           FBlk.begin(1), FBlk.begin(2), FBlk.begin(3),
                           res.begin());
  }
  
  resultat[0] += res[0];
  resultat[1] += res[1];
  resultat[2] += res[2];
  
  return 1;
}
  
  
//=============================================================================
// Integre les grandeurs de M = OM^F
// Retourne 1 si success, 0 si echec
//=============================================================================
E_Int K_POST::integUnstruct4(E_Int center2node,
                             E_Int posx, E_Int posy, E_Int posz,
                             E_Float cx, E_Float cy, E_Float cz, 
                             FldArrayI& cnBlk, FldArrayF& coordBlk, 
                             FldArrayF& FBlk, FldArrayF& ratioBlk, 
                             FldArrayF& resultat)
{
  FldArrayF res(3);
  
  E_Int size = coordBlk.getSize();
  E_Int nbT = cnBlk.getSize();
  FldArrayF surfBlk(nbT);

  // Compute surface of each "block" i cell, with coordinates coordBlk
  FldArrayF snx(nbT); // normale a la surface  
  FldArrayF sny(nbT);
  FldArrayF snz(nbT);
  K_METRIC::compUnstructSurf(
    cnBlk, "TRI",
    coordBlk.begin(posx), coordBlk.begin(posy), coordBlk.begin(posz),
    snx.begin(), sny.begin(), snz.begin(), surfBlk.begin());

  switch (center2node)
  {
    case 1:
    // Compute integral, coordinates defined in node 
    // and field FBlk in center 
      k6integmomentunsnodecenter_(nbT, size, cnBlk.begin(), cx, cy, cz, 
                                  ratioBlk.begin(), coordBlk.begin(posx), 
                                  coordBlk.begin(posy),coordBlk.begin(posz),
                                  surfBlk.begin(), FBlk.begin(1), 
                                  FBlk.begin(2), FBlk.begin(3), res.begin());
      break;

    default:
      // Compute integral, coordinates and field have the same size
      k6integmomentunstruct_(nbT, size, cnBlk.begin(), cx, cy, cz, 
                             ratioBlk.begin(), coordBlk.begin(posx),
                             coordBlk.begin(posy), coordBlk.begin(posz),
                             surfBlk.begin(), FBlk.begin(1), FBlk.begin(2),
                             FBlk.begin(3), res.begin());
      break;
  }
  resultat[0] += res[0];
  resultat[1] += res[1];
  resultat[2] += res[2];
   
  return 1;
}
//=============================================================================
// Integre les grandeurs de M = OM^F
// Retourne 1 si success, 0 si echec
//=============================================================================
E_Int K_POST::integUnstruct41D(E_Int center2node,
                               E_Int posx, E_Int posy, E_Int posz,
                               E_Float cx, E_Float cy, E_Float cz, 
                               FldArrayI& cnBlk, FldArrayF& coordBlk, 
                               FldArrayF& FBlk, FldArrayF& ratioBlk, 
                               FldArrayF& resultat)
{
  FldArrayF res(3);
  
  E_Int size = coordBlk.getSize();
  E_Int nbT = cnBlk.getSize();
  FldArrayF lengthBlk(nbT);

  // Compute surface of each "block" i cell, with coordinates coordBlk
  K_METRIC::compUnstructSurf1d(
    cnBlk, "BAR",
    coordBlk.begin(posx), coordBlk.begin(posy), coordBlk.begin(posz),
    lengthBlk.begin());

  switch (center2node) 
  { 
    case 1:
      // Compute integral, coordinates defined in node 
      // and field FBlk in center 
      k6integmomentunsnodecenter1d_(
        nbT, size, cnBlk.begin(), cx, cy, cz, 
        ratioBlk.begin(), coordBlk.begin(posx),
        coordBlk.begin(posy), coordBlk.begin(posz), 
        lengthBlk.begin(), FBlk.begin(1),
        FBlk.begin(2), FBlk.begin(3),
        res.begin());
      break;
    default:
      // Compute integral, coordinates and field have the same size
      k6integmomentunstruct1d_(
        nbT, size, cnBlk.begin(), cx, cy, cz, 
        ratioBlk.begin(), coordBlk.begin(posx),
        coordBlk.begin(posy), coordBlk.begin(posz), 
        lengthBlk.begin(), FBlk.begin(1), FBlk.begin(2),
        FBlk.begin(3), res.begin());
      break;
  }

  resultat[0] += res[0];
  resultat[1] += res[1];
  resultat[2] += res[2];   

  return 1;
}
  
  
