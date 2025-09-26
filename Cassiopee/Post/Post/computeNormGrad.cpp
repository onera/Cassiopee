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
# include "post.h"

using namespace std;

//=============================================================================
/* Calcul de la norme du gradient d'un ensemble de champs definis en noeuds. 
   La norme du gradient est fournie aux centres des cellules */
//=============================================================================
PyObject* K_POST::computeNormGrad(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* varname;
  if (!PYPARSETUPLE_(args, OO_, &array, &varname)) return NULL;
  
  // Check array
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk; // number of points of array
  E_Int posx = -1; E_Int posy = -1; E_Int posz = -1;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, cn, 
                                     eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeNormGrad: invalid array.");
    return NULL;
  }
  // check varname
  char* var = NULL;
  if (PyString_Check(varname)) var = PyString_AsString(varname);
#if PY_VERSION_HEX >= 0x03000000
  else if (PyUnicode_Check(varname)) var = (char*)PyUnicode_AsUTF8(varname);
#endif
  else
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeNormGrad: varname must be a string.");
    RELEASESHAREDB(res,array,f,cn); return NULL;
  } 
  E_Int posv = K_ARRAY::isNamePresent(var, varString);
  if (posv == -1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeNormGrad: variable not found in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posv++;
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeNormGrad: coordinates not found in array.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  posx++; posy++; posz++;
  E_Int api = f->getApi();
  E_Int nfld = f->getNfld();
  if (nfld == 3)
  {
    PyErr_SetString(PyExc_TypeError, 
                    "computeNormGrad: no field to compute.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  
  PyObject* tpl = NULL;
  FldArrayF* f2;
  E_Int sizeVarStringOut = strlen(var) + 5; // grad var
  char* varStringOut = new char[sizeVarStringOut];
  strcpy(varStringOut, "grad"); strcat(varStringOut,var);

  if (res == 1)
  {
    E_Int ni1 = K_FUNC::E_max(1, ni-1);
    E_Int nj1 = K_FUNC::E_max(1, nj-1);
    E_Int nk1 = K_FUNC::E_max(1, nk-1);
    E_Int ncells = ni1*nj1*nk1;
    tpl = K_ARRAY::buildArray3(1, varStringOut, ni1, nj1, nk1, api);
    
    K_ARRAY::getFromArray3(tpl, f2);
    E_Float* fnp = f2->begin(1);
    E_Float* fyp = new E_Float[ncells];
    E_Float* fzp = new E_Float[ncells];

    computeGradStruct(ni, nj, nk, 
                      f->begin(posx), f->begin(posy), f->begin(posz), 
                      f->begin(posv),
                      fnp, fyp, fzp);

    // stockage de la norme
    #pragma omp parallel for if (ncells > __MIN_SIZE_MEAN__)
    for (E_Int i = 0; i < ncells; i++)
      fnp[i] = sqrt(fnp[i]*fnp[i] + fyp[i]*fyp[i] + fzp[i]*fzp[i]);

    delete [] fyp; delete [] fzp;
  }
  else // non structure 
  {
    E_Int npts = f->getSize();
    E_Bool center = true;
    E_Bool copyConnect = true;

    if (strcmp(eltType, "NGON") == 0)
    {
      // Build unstructured NGON array from existing connectivity & empty fields
      E_Int nelts = cn->getNElts();

      tpl = K_ARRAY::buildArray3(
        1, varStringOut, npts, *cn, "NGON",
        center, api, copyConnect
      );
      K_ARRAY::getFromArray3(tpl, f2);

      E_Float* fnp = f2->begin(1);
      E_Float* fyp = new E_Float[nelts];
      E_Float* fzp = new E_Float[nelts];
      
      E_Int ierr = computeGradNGon(
        f->begin(posx), f->begin(posy), f->begin(posz), 
        f->begin(posv), *cn, 
        fnp, fyp, fzp);

      if (ierr == 1)
      {
        PyErr_SetString(PyExc_TypeError, 
                        "computeGrad: gradient can only be computed for 3D NGONs.");
        delete [] varStringOut;
        delete [] fyp; delete [] fzp;
        RELEASESHAREDS(tpl, f2);
        RELEASESHAREDB(res, array, f, cn); return NULL;         
      }

      // stockage de la norme
      #pragma omp parallel for if (nelts > __MIN_SIZE_MEAN__)
      for (E_Int i = 0; i < nelts; i++)
        fnp[i] = sqrt(fnp[i]*fnp[i] + fyp[i]*fyp[i] + fzp[i]*fzp[i]);

      delete [] fyp; delete [] fzp;
    }
    else  // ME
    {
      E_Int nc = cn->getNConnect();
      E_Int ntotElts = 0;
      for (E_Int ic = 0; ic < nc; ic++)
      {
        K_FLD::FldArrayI& cm = *(cn->getConnect(ic));
        E_Int nelts = cm.getSize();
        ntotElts += nelts;
      }
      
      tpl = K_ARRAY::buildArray3(
        1, varStringOut, npts, *cn, eltType,
        center, api, copyConnect
      );
      K_ARRAY::getFromArray3(tpl, f2);
      E_Float* fnp = f2->begin(1);
      E_Float* fyp = new E_Float[ntotElts];
      E_Float* fzp = new E_Float[ntotElts];

      computeGradUnstruct(
        f->begin(posx), f->begin(posy), f->begin(posz), *cn, eltType,
        f->begin(posv), fnp, fyp, fzp);

      // stockage de la norme
      #pragma omp parallel for if (ntotElts > __MIN_SIZE_MEAN__)
      for (E_Int i = 0; i < ntotElts; i++)
        fnp[i] = sqrt(fnp[i]*fnp[i] + fyp[i]*fyp[i] + fzp[i]*fzp[i]);

      delete [] fyp; delete [] fzp;
    }
  }
  delete [] varStringOut;
  RELEASESHAREDS(tpl, f2);
  RELEASESHAREDB(res, array, f, cn);
  return tpl;
}
