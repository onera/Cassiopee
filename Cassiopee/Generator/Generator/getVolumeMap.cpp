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

// getVolumeMap

# include "generator.h"

using namespace K_FUNC;
using namespace K_CONST;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Return volume map */
// ============================================================================
PyObject* K_GENERATOR::getVolumeMapOfMesh(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int method;

  if ( !PYPARSETUPLE_(args, O_ I_, &array, &method) ) return NULL;
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int posx, posy, posz;
  E_Int res;
  res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  PyObject* tpl = NULL;
  
  if (res == 1 || res == 2)
  {
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDB(res, array, f, cn);
      PyErr_SetString(PyExc_ValueError,
                      "getVolumeMap: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;

    // Coordonnees
    E_Float* xt = f->begin(posx);
    E_Float* yt = f->begin(posy);
    E_Float* zt = f->begin(posz);

    E_Int api = f->getApi();
    E_Int npts = f->getSize();

    if (res == 1) // cas structure
    {
      E_Int dim = 3;
      E_Int im1 = im-1;
      E_Int jm1 = jm-1;
      E_Int km1 = km-1;
      if (im1*jm1*km1 == 0)
      {
        if ( (im1 && (jm1 || km1)) || (jm1 && km1) ) dim = 2;
        //if ((im1*jm1)||(im1*km1)||(jm1*km1)) dim = 2;
	      else dim = 1;
      }
      if (im == 1) im1 = 1;
      if (jm == 1) jm1 = 1;
      if (km == 1) km1 = 1;
     
      //E_Int ncells = im1*jm1*km1;
      E_Int ninti = im*jm1*km1;
      E_Int nintj = im1*jm*km1;
      E_Int nintk = im1*jm1*km;
      E_Int nint =  ninti + nintj + nintk;

      tpl = K_ARRAY::buildArray3(1, "vol", im1, jm1, km1, api);
      FldArrayF* f2;
      K_ARRAY::getFromArray3(tpl, f2);
      E_Float* volap = f2->begin(1);

      // calcul du volume
      if (dim == 1)
        K_METRIC::compSurfStruct1D(im, jm, km, xt, yt, zt, volap);
      else if (dim == 2)
        K_METRIC::compSurfStruct2D(im, jm, km, xt, yt, zt, volap);
      else
      {
        FldArrayF surf(nint,3);
        FldArrayF snorm(nint);
        FldArrayF centerInt(nint, 3);
        K_METRIC::compMetricStruct(
          im, jm, km, ninti, nintj, nintk,
          xt, yt, zt,
          volap, surf.begin(1), surf.begin(2), surf.begin(3), 
          snorm.begin(), 
          centerInt.begin(1), centerInt.begin(2), centerInt.begin(3));
      }
      RELEASESHAREDS(tpl, f2);
      RELEASESHAREDS(array, f);
    }
    else if (res == 2) // cas non-structure
    {
      // Build array contenant le volume des elements au centre des elements
      tpl = K_ARRAY::buildArray3(1, "vol", npts, *cn, eltType, true, api, true);
      FldArrayF* f2; K_ARRAY::getFromArray3(tpl, f2);
      E_Float* vol = f2->begin(1); // pointeur sur le tableau de volumes

      if (strcmp(eltType, "NGON") == 0) // Elements NGON
      { 
        E_Int ierr = 0;
        if (method == 0)
          ierr = K_METRIC::compVolNGon(xt, yt, zt, *cn, vol);
        else if (method == 1)
          ierr = K_METRIC::compVolNGonImad(xt, yt, zt, *cn, vol);
        else
        {
          PyErr_SetString(PyExc_ValueError,
                          "getVolumeMap: wrong method (should be 0 or 1).");
          RELEASESHAREDS(tpl, f2);
          RELEASESHAREDU(array, f, cn);
          return NULL;
        }

        if (ierr == 1)
        {
          PyErr_SetString(PyExc_ValueError,
                          "getVolumeMap: dimension of element unknown.");
          RELEASESHAREDS(tpl, f2);
          RELEASESHAREDU(array, f, cn);
          return NULL;
        }
      }
      else // ME
      {
        E_Int nc = cn->getNConnect();

        // Get dimensionality
        E_Int dim = K_CONNECT::getDimME(eltType);
        
        if (dim == 1)
        {
          K_FLD::FldArrayI& cm = *(cn->getConnect(0)); // TODO
          E_Int nelts = cm.getSize();
          E_Int* cn1 = cn->begin(1);
          E_Int* cn2 = cn->begin(2);
          E_Int ind1, ind2;
          E_Float dx, dy, dz;
          for (E_Int i = 0; i < nelts; i++)
          {
            ind1 = cn1[i] - 1; 
            ind2 = cn2[i] - 1;
            dx = xt[ind2] - xt[ind1];
            dy = yt[ind2] - yt[ind1];
            dz = zt[ind2] - zt[ind1];
            vol[i] = sqrt(dx*dx + dy*dy + dz*dz);
          }
        }
        else if (dim == 2)
        {
          // Compute total number of elements
          E_Int ntotElts = 0;
          for (E_Int ic = 0; ic < nc; ic++)
          {
            K_FLD::FldArrayI& cm = *(cn->getConnect(ic));
            E_Int nelts = cm.getSize();
            ntotElts += nelts;
          }

          // Allocate memory to store facet normals and their areas for all
          // connectivities
          FldArrayF snx(ntotElts), sny(ntotElts), snz(ntotElts);
          
          // Compute surface of elements
          K_METRIC::compSurfUnstruct(
            *cn, eltType, xt, yt, zt,
            snx.begin(), sny.begin(), snz.begin(), vol
          );
        }
        else if (dim == 3)
        {
          // Number of facets per element
          std::vector<E_Int> nfpe;
          E_Int ierr = K_CONNECT::getNFPE(nfpe, eltType, false);
          if (ierr != 0)
          {
            RELEASESHAREDS(tpl, f2);
            RELEASESHAREDU(array, f, cn);
            return NULL;
          }
          
          // Compute total number of facets
          E_Int ntotFacets = 0;
          for (E_Int ic = 0; ic < nc; ic++)
          {
            K_FLD::FldArrayI& cm = *(cn->getConnect(ic));
            E_Int nelts = cm.getSize();
            ntotFacets += nfpe[ic]*nelts;
          }
          
          FldArrayF snx(ntotFacets), sny(ntotFacets), snz(ntotFacets);
          FldArrayF surf(ntotFacets);
          K_METRIC::compMetricUnstruct(
            *cn, eltType,
            xt, yt, zt,
            snx.begin(), sny.begin(), snz.begin(), surf.begin(), vol
          );
        }
      }
      
      RELEASESHAREDS(tpl, f2);
      RELEASESHAREDU(array, f, cn);
    }
  }
  else
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "getVolumeMap: unknown type of array.");
    return NULL;
  }
  return tpl;
}
