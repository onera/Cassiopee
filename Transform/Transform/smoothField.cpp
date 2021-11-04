/*    
    Copyright 2013-2021 Onera.

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

# include "transform.h"
using namespace std;
using namespace K_FLD;
using namespace K_SEARCH;

// ============================================================================
/* Smooth a field */
// ============================================================================
PyObject* K_TRANSFORM::_smoothField(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float eps;  
  PyObject* epsl; // numpy de eps (optional)
  E_Int niter, type;
  PyObject* varList;
  if (!PYPARSETUPLE(args,
                    "OdOllO", "OdOiiO",
                    "OfOllO", "OfOiiO",
                    &array, &eps, &epsl, &niter, &type, &varList))
  {
      return NULL;
  }

  // Check varList
  if (PyList_Check(varList) == false)
  {
    PyErr_SetString(PyExc_TypeError,
                    "smoothField: varList must be a list of variables.");
    return NULL;
  }
  // nbre de variables a lisser
  E_Int nvars = PyList_Size(varList);
  
  // Check epsl
  E_Float* epsf = NULL;
  if (epsl != Py_None)
  {
    E_Int size, nfld;
    K_NUMPY::getFromNumpyArray(epsl, epsf, size, nfld, true);
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray2(array, varString, 
                                     f, im, jm, km, cn, eltType);

  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "smoothField: array must be unstructured.");
    return NULL;
  }  
  
  // position des variables a lisser
  vector<E_Int> posVars(nvars);
  for (E_Int i = 0; i < nvars; i++)
  {
      PyObject* varname = PyList_GetItem(varList, i);
      char* var = NULL;
      if (PyString_Check(varname)) var = PyString_AsString(varname);
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(varname)) var = (char*)PyUnicode_AsUTF8(varname);
#endif
      E_Int pos = K_ARRAY::isNamePresent(var, varString);
      if (pos == -1) { printf("Warning: smoothField: var doesn't exists.\n"); posVars[i] = 0; }
      else posVars[i] = pos;
  }

  // nbre de vertex du maillage
  E_Int npts = f->getSize();
  vector< vector<E_Int> > cVN(npts); // vertex/elt
  if (strcmp(eltType, "NGON") == 0) K_CONNECT::connectNG2VNbrs(*cn, cVN);
  else K_CONNECT::connectEV2VNbrs(*cn, cVN);

  // local temp arrays
  E_Float* ft = new E_Float [npts*nvars];
  if (epsl == Py_None) 
  { 
    epsf = new E_Float [npts];
    # pragma omp parallel for
    for (E_Int ind = 0; ind < npts; ind++) epsf[ind] = eps;
  }

  if (type == 0) // isotrope 
  {
#pragma omp parallel
  {
    E_Int nbV, nov;
    E_Float eps2, df;
    E_Float* fp; E_Float* tp;

    for (E_Int it = 0; it < niter; it++)
    {
        for (E_Int nv = 0; nv < nvars; nv++)
        {
            fp = f->begin(posVars[nv]+1);
            tp = ft+npts*nv;
                
            #pragma omp for
            for (E_Int ind = 0; ind < npts; ind++)
            {
                vector<E_Int>& v = cVN[ind]; // vertex voisins
                nbV = v.size();
                eps2 = epsf[ind]/nbV;
            
                df = 0.;
                for (E_Int vi = 0; vi < nbV; vi++)
                {
                    nov = v[vi]-1;
                    df += fp[nov]-fp[ind];
                }
                tp[ind] = fp[ind] + eps2*df;
            }

            // end loop update
            #pragma omp for
            for (E_Int ind = 0; ind < npts; ind++)
            {
                fp[ind] = tp[ind];
            }
        }
    }
  }
  }
  else if (type == 1) // scaled
  {
      /* a ajouter si type = 0 est ok */
  }
  delete [] ft;

  RELEASESHAREDB(res, array, f, cn);
  if (epsl != Py_None) { Py_DECREF(epsl); }
  else { delete [] epsf; }
  
  Py_INCREF(Py_None);
  return Py_None; 
}
