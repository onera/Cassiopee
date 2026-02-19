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
#include "converter.h"
using namespace std;
using namespace K_FLD;

//=============================================================================
/* Convertit un tableau defini en noeuds en un tableau defini en centres */
// ============================================================================
PyObject* K_CONVERTER::center2Node(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int type;
  PyObject* BCFields; // optional indR+fields on BCs
  if (!PYPARSETUPLE_(args, O_ I_ O_, &array, &type, &BCFields)) return NULL;

  // Check BCFields
  if (BCFields != Py_None)
  {
    if (PyList_Check(BCFields) == false) 
    {
      PyErr_SetString(PyExc_TypeError, 
                       "center2Node: BCFields must be a list of indR+fields.");
      return NULL;
    }
  }

  E_Int ni, nj, nk;
  FldArrayI* cnc; FldArrayF* fc;
  char* eltType; char* varString;

  E_Int res = K_ARRAY::getFromArray3(array, varString, fc, ni, nj, nk, cnc, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "center2Node: unrecognised type of array.");
    return NULL;
  }

  E_Int nelts = fc->getSize();
  E_Int nfld = fc->getNfld();
  E_Int api = fc->getApi();
  E_Int npts;

  // Essaie de trouver la variables cellN. Les traitements sont un peu
  // differents pour cette variable.
  E_Int posCellN = K_ARRAY::isCellNatureField2Present(varString);
  if (posCellN != -1) posCellN += 1;

  // Retourne le mode de sortie du champ cellnaturefield (1 (0,1), 2 (0,1,2) 
  // ou 3 (0, 1, -interpolationblock))
  E_Int mod = 0;
  E_Float nature;
  if (posCellN != -1)
  {
    mod = 1;
    E_Float* cellNp = fc->begin(posCellN);
    for (E_Int i = 0; i < nelts; i++) 
    {
      nature = cellNp[i];
      if (K_FUNC::fEqualZero(nature-2.)) { mod = 2; break; }
      else if (nature < -0.2) { mod = 3; break; }
    }
  }

  PyObject* tpl;
  FldArrayF* fn2;
  E_Int ret;

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  
  if (res == 1)
  {
    E_Int nin, njn, nkn;
    
    if (ni == 1)
    {
      if ((nj != 1)&&(nk != 1))
      { nin = 1; njn = nj+1; nkn = nk+1; }
      else if (nj == 1)
      { nin = 1; njn = 1; nkn = nk+1; }
      else //if (nk == 1)
      { nin = 1; njn = 1; nkn = nk+1; }
    }
    else if (nj == 1)
    {
      if ((ni != 1)&&(nk != 1))
      { nin = ni+1; njn = 1; nkn = nk+1; }
      else if (ni == 1)
      { nin = 1; njn = 1; nkn = nk+1; }
      else // if (nk == 1)
      { nin = ni+1; njn = 1; nkn = 1; }
    }
    else if (nk == 1)
    {
      if ((ni != 1)&&(nj != 1))
      { nin = ni+1; njn = nj+1; nkn = 1; }
      else if (ni == 1)
      { nin = 1; njn = nj+1; nkn = 1; }
      else //if (nj == 1)
      { nin = ni+1; njn = 1; nkn = 1; }
    }
    else { nin = ni+1; njn = nj+1; nkn = nk+1; }

    tpl = K_ARRAY::buildArray3(nfld, varString, nin, njn, nkn);
    K_ARRAY::getFromArray3(tpl, fn2);

    ret = K_LOC::center2nodeStruct(*fc, ni, nj, nk, posCellN, mod, 
                                   posx, posy, posz, *fn2, nin, njn, nkn,
                                   type);
    // Boundary corrections
    if (BCFields != Py_None)
    {
      //PyObject* indR = PyList_GetItem(BCFields, 0);
      //PyObject* fields = PyList_GetItem(BCFields, 1);
      //E_Int res = K_ARRAY::getFromArray3(fields, varString, fc, 
      //                                   ni, nj, nk, cnc, eltType);
      //center2NodeStructBorder(fn2, nin, njn, nkn);
      //RELEASESHAREDB(res, fields);
    }

    RELEASESHAREDS(tpl, fn2);
    RELEASESHAREDS(array, fc);
    if (ret == 0)
    {
      PyErr_SetString(PyExc_ValueError, 
                      "center2Node: algo failed for the structured array.");
      return NULL;
    }
    return tpl;
  }
  else  // res == 2
  {
    if (strchr(eltType, '*') == NULL)
    {
      PyErr_SetString(PyExc_TypeError, 
                      "center2Node: unstructured array must be eltType*.");
      RELEASESHAREDU(array, fc, cnc); return NULL;
    }

    if (K_STRING::cmp(eltType, "NODE*") == 0)
    {
      npts = nelts;
      tpl = K_ARRAY::buildArray3(nfld, varString, npts, nelts, "NODE", true, api);
      K_ARRAY::getFromArray3(tpl, fn2);
      
      #pragma omp parallel
      {
        for (E_Int n = 1; n <= nfld; n++)
        {
          E_Float* fcp = fc->begin(n);
          E_Float* fn2p = fn2->begin(n);

          #pragma omp for nowait
          for (E_Int i = 0; i < npts; i++) fn2p[i] = fcp[i]; 
        }
      }

      RELEASESHAREDS(tpl, fn2);
      RELEASESHAREDU(array, fc, cnc);
      return tpl;
    }
      
    npts = 0;

    if (K_STRING::cmp(eltType, "NGON*") == 0)
    {
      vector<vector<E_Int> > cEV(nelts);
      K_CONNECT::connectNG2EV(*cnc, cEV);

      // Calcul de npts en prenant le numero max du vertex dans la conn NGON
      #pragma omp parallel reduction(max: npts)
      {
        E_Int indv, nv;
        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          const vector<E_Int>& vertices = cEV[i];
          nv = vertices.size();
          for (E_Int j = 0; j < nv; j++)
          {
            indv = vertices[j];
            npts = K_FUNC::E_max(npts, indv);
          }
        }
      }

      PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts,
                                           *cnc, "NGON", false, api, true);
      FldArrayF* fn2;
      K_ARRAY::getFromArray3(tpl, fn2);
      ret = K_LOC::center2nodeNGon(*fc, *cnc, cEV, *fn2, posCellN, mod, type);
      cEV.clear();

      RELEASESHAREDS(tpl, fn2);
      RELEASESHAREDU(array, fc, cnc);
      if (ret == 0)
      {
        PyErr_SetString(PyExc_ValueError, 
                        "center2Node: algo failed for the NGON array.");
        return NULL;
      }
      return tpl;
    }
    else // BE/ME
    {    
      E_Int nc = cnc->getNConnect();
      char* eltType2 = new char[K_ARRAY::VARSTRINGLENGTH];
      K_ARRAY::unstarVarString(eltType, eltType2);

      // Calcul de npts en prenant le numero max dans les conns BE
      #pragma omp parallel reduction(max: npts)
      {
        for (E_Int ic = 0; ic < nc; ic++)
        {
          FldArrayI& cmc = *(cnc->getConnect(ic));
          E_Int nelts = cmc.getSize();
          E_Int nvpe = cmc.getNfld();
          
          #pragma omp for collapse(2)
          for (E_Int i = 0; i < nelts; i++)
          for (E_Int n = 1; n <= nvpe; n++)
          {
            npts = K_FUNC::E_max(npts, cmc(i,n));
          }
        }
      }

      tpl = K_ARRAY::buildArray3(nfld, varString, npts,
                                 *cnc, eltType2, false, api, true);
      K_ARRAY::getFromArray3(tpl, fn2);

      ret = K_LOC::center2nodeUnstruct(*fc, *cnc, posCellN, mod,
                                        posx, posy, posz, *fn2, type);

      delete[] eltType2;
      RELEASESHAREDS(tpl, fn2);
      RELEASESHAREDU(array, fc, cnc);

      if (ret == 0)
      {
        PyErr_SetString(PyExc_ValueError, 
                        "center2Node: algo failed for the BE/ME array.");
        return NULL;
      }
      return tpl;
    }
  }
}
