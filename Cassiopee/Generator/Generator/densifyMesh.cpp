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

# include "generator.h"

using namespace K_FLD;
using namespace K_CONST;
using namespace std;

//=============================================================================
// Densifie le maillage d'un i-array
// IN: i-array
// IN: pas de discretisation h
//=============================================================================
PyObject* K_GENERATOR::densifyMesh(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float h;
  if (!PYPARSETUPLE_(args, O_ R_, &array, &h))
    return NULL;
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  FldArrayF* newf;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType);
  
  E_Int len, posx, posy, posz, ni;
  E_Float l, fn;
  E_Float df;
  E_Float hinv = 1./h;

  if (res == 1)
  {      
    posx = K_ARRAY::isCoordinateXPresent( varString);
    posy = K_ARRAY::isCoordinateYPresent( varString);
    posz = K_ARRAY::isCoordinateZPresent( varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      PyErr_SetString(PyExc_TypeError,
                      "densify: can't find coordinates in array.");
      return NULL;    
    }
    posx++; posy++; posz++;
    E_Float* x = f->begin(posx);
    E_Float* y = f->begin(posy);
    E_Float* z = f->begin(posz);
    E_Int npts = f->getSize();
    FldArrayI n(npts);
    E_Int nfld = f->getNfld();

// Size of new field
// -------------------
    len = 1;

    for (E_Int i = 1; i < npts; i++)
    {
      E_Float dx = x[i]-x[i-1];
      E_Float dy = y[i]-y[i-1];
      E_Float dz = z[i]-z[i-1];
      l = sqrt(dx*dx+dy*dy+dz*dz);
      n[i-1]=(E_Int)(l*hinv);
      if ((l*hinv-(E_Float)(n[i-1])) != 0.) n[i-1]=n[i-1]+1;
      if (l == 0.) n[i-1]=n[i-1]+1;
      len=len+n[i-1];
    }

// Compute new field
// -----------------------
// declaration of new field array
    newf = new FldArrayF(len, nfld);
    newf->setAllValuesAtNull();

    for (E_Int fld = 1; fld <= nfld; fld++)
    {    
      ni=0;
      E_Float* fp = f->begin(fld);
      E_Float* newfp = newf->begin(fld);

      for (E_Int i = 1; i < npts; i++)
      {
        newfp[ni] = fp[i-1];
        ni++;
        l=fp[i]-fp[i-1];
        if (l != 0)
        {
          df = l/n[i-1];
          for (E_Int ind=1;ind<n[i-1];ind++)
          {
            fn=fp[i-1]+ind*df;
            newfp[ni] = fn;
            ni++;
          }
        }
        else
        {
          for (E_Int ind=1;ind<n[i-1];ind++)
          {
            newfp[ni] = fp[i-1];
            ni++;
          }        
        }
      }
      newfp[ni] = fp[npts-1];
    }
    
    delete f;
    // Build array
    PyObject* tpl = K_ARRAY::buildArray(*newf, varString, 
                                        len, 1, 1);
    delete newf;
    return tpl;  
  }
  else if (res == 2)
  {
    if (strcmp(eltType, "BAR") == 0) 
    {
      posx = K_ARRAY::isCoordinateXPresent(varString);
      posy = K_ARRAY::isCoordinateYPresent(varString);
      posz = K_ARRAY::isCoordinateZPresent(varString);
      if (posx == -1 || posy == -1 || posz == -1)
      {
        delete f; delete cn;
        PyErr_SetString(PyExc_TypeError,
                        "densify: coordinates not found in array.");
        return NULL;
      }
      posx++; posy++; posz++;
 
      E_Int ind1,ind2;
      E_Int nelts = cn->getSize();
      FldArrayI n(nelts);
      E_Int* cn1 = cn->begin(1);
      E_Int* cn2 = cn->begin(2);
      E_Float* x = f->begin(posx);
      E_Float* y = f->begin(posy);
      E_Float* z = f->begin(posz);
      E_Int nfld = f->getNfld();

// Size of coordinates and new connectivity array
// ----------------------------------------------
      len = 1;
      E_Int npts = f->getSize();
      for (E_Int i = 0; i < nelts; i++)
      {
        ind1 = cn1[i]-1;
        ind2 = cn2[i]-1;
        E_Float dx = x[ind2]-x[ind1];
        E_Float dy = y[ind2]-y[ind1];
        E_Float dz = z[ind2]-z[ind1];
        l = sqrt(dx*dx+dy*dy+dz*dz);
        n[i] = (E_Int)(l*hinv);
        if ((l*hinv-(E_Float)(n[i])) != 0.) n[i]=n[i]+1;
        if (l == 0.) n[i] = n[i]+1;
        len=len+n[i];
      }
      
// Compute new coordinates
// -----------------------
// declaration of new coordinate array
      newf = new FldArrayF(len, nfld);
      FldArrayI* newcn = new FldArrayI(len-1, 2);
      E_Int* newcn1 = newcn->begin(1);
      E_Int* newcn2 = newcn->begin(2);
      E_Int nc1 = 0;
      E_Int nc2 = 0;
      for (E_Int i = 0; i < nelts; i++)
      {
        ind1 = cn1[i]-1;
        ind2 = cn2[i]-1;
        newcn1[nc1] = cn1[i];
        for (E_Int j = 1; j < n[i]; j++)
        {
          newcn2[nc1] = npts+1+nc2; nc1++;
          newcn1[nc1] = npts+1+nc2; nc2++;
        }
        newcn2[nc1] = cn2[i]; nc1++;
      }

// copy old field values in new field to keep the old ordering
// -----------------------------------------------------------
      for (E_Int j = 1; j <= nfld; j++)
      {
        E_Float* fp = f->begin(j);
        E_Float* newfp = newf->begin(j);
        for (E_Int i = 0; i < npts; i++)
          newfp[i] = fp[i];
      }

      for (E_Int fld=1;fld<=nfld;fld++)
      {    
        ni = npts;
        E_Float* fp = f->begin(fld);
        E_Float* newfp = newf->begin(fld);
        for (E_Int i = 0; i < nelts; i++)
        {
          ind1 = cn1[i]-1;
          ind2 = cn2[i]-1;
          l = fp[ind2]-fp[ind1];
          if (l != 0)
          {
            df=l/n[ind1];
            for (E_Int ind=1;ind<n[i];ind++)
            {
              fn=fp[ind1]+ind*df;
              newfp[ni]=fn;
              ni++;
            }
          }
          else
          {
            for (E_Int ind=1;ind<n[i];ind++)
            {
              newfp[ni]=fp[ind1];
              ni++;
            }        
          }
        }
      }
      
      // Build array
      delete f; delete cn;
      PyObject* tpl = K_ARRAY::buildArray(*newf, varString, 
                                          *newcn, 1, NULL, false);
      delete newf; delete newcn;
      return tpl;
    }
    else
    {
      delete f; delete cn;
      PyErr_SetString(PyExc_TypeError,
                      "densify: not a valid type of elements.");
      return NULL;
    }
  }
  else
  {
    delete f; delete cn;
    PyErr_SetString(PyExc_TypeError,
                    "densify: invalid array.");
    return NULL;
  }
}
