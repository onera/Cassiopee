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

// Information on geometries

# include <string.h>
# include "geom.h"
using namespace K_FUNC;
using namespace K_FLD;
using namespace std;
using namespace K_CONST;

// ============================================================================
/* Get the length of a 1D array defining a mesh */
// ============================================================================
PyObject* K_GEOM::getLength(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);
  
  E_Int posx, posy, posz;
  E_Float length;
  E_Float dx, dy, dz;

  length = 0.;
  if (res == 1)
  {      
    if (jm != 1 || km != 1)
      printf("Warning: getLength: only line j=1, k=1 is taken into account.\n");
    
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "getLength: can't find coordinates in array.");
      return NULL;    
    }
    posx++; posy++; posz++;
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);
    
    for (E_Int i = 1; i < im; i++)
    {
      E_Int i1 = i-1;
      dx = xt[i] - xt[i1];
      dy = yt[i] - yt[i1];
      dz = zt[i] - zt[i1];
      length += sqrt(dx*dx+dy*dy+dz*dz);
    }
    RELEASESHAREDS(array, f);
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
        RELEASESHAREDU(array, f, cn);
        PyErr_SetString(PyExc_TypeError,
                        "getLength: coordinates not found in array.");
        return NULL;
      }
      posx++; posy++; posz++;
 
      E_Int* cn1 = cn->begin(1);
      E_Int* cn2 = cn->begin(2);
      E_Float* xt = f->begin(posx);
      E_Float* yt = f->begin(posy);
      E_Float* zt = f->begin(posz);

      for (E_Int i = 0; i < cn->getSize(); i++)
      {
        E_Int ind1 = cn1[i]-1;
        E_Int ind2 = cn2[i]-1;
        dx = xt[ind2] - xt[ind1];
        dy = yt[ind2] - yt[ind1];
        dz = zt[ind2] - zt[ind1];
        length += sqrt(dx*dx + dy*dy + dz*dz);
      }
      RELEASESHAREDU(array, f, cn);
    }
    else 
    {
      RELEASESHAREDU(array, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "getLength: not a valid element type.");
      return NULL;
    }
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "getLength: invalid array.");
    return NULL;
  }
#ifdef E_DOUBLEREAL
  return Py_BuildValue("d", length);
#else
  return Py_BuildValue("f", length);
#endif
}

// ============================================================================
/* Get the distant index */
// ============================================================================
PyObject* K_GEOM::getDistantIndex(PyObject* self, PyObject* args)
{
  E_Float eps = 1.e-10;
  PyObject* array;
  E_Int ind;
  E_Float l;

  if (!PYPARSETUPLE_(args, O_ I_ R_,
                    &array, &ind, &l))
  {
      return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  E_Int imjm, kind, jind;
  E_Float t, dx, dy, dz;
  E_Int is, bsup, binf;
  E_Float length;

  if (res == 1)
  {
    if (jm != 1 || km != 1)
      printf("Warning: getDistantIndex: only line j=1, k=1 is taken into account.\n");
    // Data check
    if (ind > im*jm*km || ind < 1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError, 
                      "getDistantIndex: index is out of mesh bounds.");
      return NULL;
    }
      
    E_Int posx = K_ARRAY::isCoordinateXPresent( varString );
    E_Int posy = K_ARRAY::isCoordinateYPresent( varString );
    E_Int posz = K_ARRAY::isCoordinateZPresent( varString );
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "getDistantIndex: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    is = ind-1;
    imjm = im*jm;
    kind = is / imjm;
    jind = (is - kind*imjm)/im;
      
    bsup = im-1 + jind*im + kind*imjm;
    binf = jind*im + kind*imjm;
    
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);

    length = 0.;
    if (l >= 0.)
    {
      while (length < l && is < bsup)
      {
        dx = xt[is+1] - xt[is];
        dy = yt[is+1] - yt[is];
        dz = zt[is+1] - zt[is];
        t = sqrt( dx*dx + dy*dy + dz*dz );
        length = length + t;
        is++;
      }
      if (is == bsup && length < l-eps)
      {
        RELEASESHAREDS(array, f);
        PyErr_SetString(PyExc_ValueError,
                        "getDistantIndex: max bound reached without matching length.");
        return NULL;
      }
    }
    else
    {
      while (length > l && is > binf)
      {
        dx = xt[is] - xt[is-1];
        dy = yt[is] - yt[is-1];
        dz = zt[is] - zt[is-1];
        t = sqrt( dx*dx + dy*dy + dz*dz);
        length = length - t;
        is--;
      }
      if (is == binf && length > l+eps)
      {
         RELEASESHAREDS(array, f);
        PyErr_SetString(PyExc_ValueError,
                        "getDistantIndex: min bound reached with matching length.");
        return NULL;
      }
    }

    RELEASESHAREDS(array, f);
    return Py_BuildValue(I_, is);
  }
  else if (res == 2)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "getDistantIndex: can not be used on an unstructured array.");
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "getDistantIndex: invalid array.");
    return NULL;
  }
}      

// ============================================================================
/* Return the curvilinear abscissa of a 1D array defining a mesh */
// ============================================================================
PyObject* K_GEOM::getCurvilinearAbscissa(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);
  
  E_Int posx, posy, posz;
  E_Float length;
  E_Float l;

  length = 0.;
  if (res == 1)
  {      
    if (jm != 1 || km != 1)
      printf("Warning: getCurvilinearAbscissa: only line j=1, k=1 is taken into account.\n");
    
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "getCurvilinearAbscissa: can't find coordinates in array.");
      return NULL;    
    }
    posx++; posy++; posz++;

    
    PyObject* tpl;
    tpl = K_ARRAY::buildArray(1, "s", im, jm, km);
    E_Float* ab = K_ARRAY::getFieldPtr(tpl);

    ab[0] = 0.;
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);
    
    E_Float dx, dy, dz;
    for (E_Int i = 1; i < im; i++)
    {
      dx = xt[i] - xt[i-1];
      dy = yt[i] - yt[i-1];
      dz = zt[i] - zt[i-1];
 
      l = sqrt(dx*dx + dy*dy + dz*dz);
      length += l;
      ab[i] = ab[i-1] + l;
    }

    E_Float inv = 1/length;
    for (E_Int i = 1; i < im; i++) ab[i] = ab[i] * inv;
    
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else if (res == 2)
  {
    if (strcmp(eltType, "BAR") != 0) 
    {
      RELEASESHAREDU(array, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "getCurvilinearAbscissa: only for BAR-array.");
      return NULL;
    }
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDU(array, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "getCurvilinearAbscissa: can't find coordinates in array.");
      return NULL;    
    }
    posx++; posy++; posz++;
    E_Int npts = f->getSize();
    PyObject* tpl = K_ARRAY::buildArray(1, "s", npts, cn->getSize(), 
                                        -1, eltType);
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn->begin(), cn->getSize()*cn->getNfld());

    E_Float* ab = K_ARRAY::getFieldPtr(tpl);
    
    ab[0] = 0.;
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);
    E_Int* cn1 = cn->begin(1);
    E_Int* cn2 = cn->begin(2);
    
    // pt de depart : 0
    E_Int nelts = cn->getSize();
    E_Int ind2; E_Float dx, dy, dz;
    FldArrayIS dejaVu(nelts); dejaVu.setAllValuesAtNull();
    short* dejaVup = dejaVu.begin();
    for (E_Int ind1 = 0; ind1 < npts-1; ind1++)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        if (dejaVup[et] == 0 && cn1[et]-1 == ind1)
        {
          ind2 = cn2[et]-1;
          dx = xt[ind2] - xt[ind1];
          dy = yt[ind2] - yt[ind1];
          dz = zt[ind2] - zt[ind1];
          
          l = sqrt(dx*dx + dy*dy + dz*dz);
          length += l;
          ab[ind2] = ab[ind1] + l;
          dejaVup[et] = 1;
          break;
        }
      }
    }
    
    E_Float inv = 1./length;
    for (E_Int i = 1; i < npts; i++) ab[i] = ab[i] * inv;
 
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "getCurvilinearAbscissa: invalid array.");
    return NULL;
  }
}
