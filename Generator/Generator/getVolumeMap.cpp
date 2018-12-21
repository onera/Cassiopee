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

// getVolumeMap

# include "generator.h"

using namespace K_FUNC;
using namespace K_CONST;
using namespace K_FLD;
using namespace std;

extern "C"
{
  void k6compstructmetric_(
    const E_Int& im, const E_Int& jm, const E_Int& km,
    const E_Int& nbcells, const E_Int& nintt,
    const E_Int& ninti, const E_Int& nintj, 
    const E_Int& nintk, 
    E_Float* x, E_Float* y, E_Float* z, 
    E_Float* vol, E_Float* surfx, E_Float* surfy, E_Float* surfz, 
    E_Float* snorm, E_Float* cix, E_Float* ciy, E_Float* ciz);

  void k6structsurft_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk, const E_Int& ncells, 
    const E_Float* xt, const E_Float* yt, const E_Float* zt, 
    E_Float* length);

  void k6structsurf1dt_(
    const E_Int& ni, const E_Int& nj, const E_Int& nk,
    const E_Float* xt, const E_Float* yt, const E_Float* zt, 
    E_Float* length);

  void k6compunstrmetric_(E_Int& npts, E_Int& nelts, E_Int& nedges, 
                          E_Int& nnodes, E_Int* cn, 
                          E_Float* coordx, E_Float* coordy, E_Float* coordz, 
                          E_Float* xint, E_Float* yint, E_Float* zint, 
                          E_Float* snx, E_Float* sny, E_Float* snz, 
                          E_Float* surf, E_Float* vol);
}
// ============================================================================
/* Return volume map */
// ============================================================================
PyObject* K_GENERATOR::getVolumeMapOfMesh( PyObject* self,
                                           PyObject* args )
{
  PyObject* array;

  if ( !PyArg_ParseTuple(args, "O", &array) ) return NULL;
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int posx, posy, posz;
  E_Int res;
  res = K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, 
                              eltType, true);

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
     
      E_Int ncells = im1*jm1*km1;
      E_Int ninti = im*jm1*km1;
      E_Int nintj = im1*jm*km1;
      E_Int nintk = im1*jm1*km;
      E_Int nint =  ninti + nintj + nintk;

      tpl = K_ARRAY::buildArray(1, "vol", im1, jm1, km1);
      E_Float* volap = K_ARRAY::getFieldPtr(tpl);

      // calcul du volume
      if (dim == 1)
	k6structsurf1dt_(
          im, jm , km , 
          xt, yt, zt, volap);
      else if (dim == 2)
        k6structsurft_(
          im, jm, km, ncells, 
          xt, yt, zt, volap);
      else
      {
        FldArrayF surf(nint,3);
        FldArrayF snorm(nint);
        FldArrayF centerInt(nint, 3);
        k6compstructmetric_(
          im, jm, km, ncells, nint, ninti, nintj, nintk,
          xt, yt, zt,
          volap, surf.begin(1), surf.begin(2), surf.begin(3), 
          snorm.begin(), 
          centerInt.begin(1), centerInt.begin(2), centerInt.begin(3));
      }
      RELEASESHAREDS(array, f);
    }
    else if (res == 2) // cas non structure
    {
      if (strcmp(eltType, "NGON") != 0) // Elements basiques
      { 
        E_Int nelts = cn->getSize(); // nb d elements ns
        E_Int nnodes = cn->getNfld(); // nb de noeuds ds 1 element
        E_Int nedges; // nb de facettes par element
        // recherche du nb de facettes par element
        if (strcmp(eltType, "BAR") == 0) nedges = 0; 
        else if (strcmp(eltType, "TRI") == 0) nedges = 1; 
        else if (strcmp(eltType, "QUAD") == 0 ) nedges = 1;
        else if (strcmp(eltType, "TETRA") == 0) nedges = 4;
        else if (strcmp(eltType, "HEXA") == 0) nedges = 6;
        else if (strcmp(eltType, "PENTA") == 0) nedges = 5;
        else if (strcmp(eltType, "PYRA") == 0) nedges = 5;
        else 
        {
          PyErr_SetString(PyExc_TypeError,
                          "getVolumeMap: unknown type of element.");
          RELEASESHAREDU(array, f, cn); return NULL;
        }

        // Build array
        tpl = K_ARRAY::buildArray(1, "vol",
                                  f->getSize(), nelts, -1, eltType, true);
        E_Float* fieldp = K_ARRAY::getFieldPtr(tpl);
        E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
        FldArrayI cnn(nelts, cn->getNfld(), cnnp, true); cnn = *cn;

        if (nedges == 0) // BAR element
        {
          E_Int* cn1 = cn->begin(1);
          E_Int* cn2 = cn->begin(2);
          for (E_Int i = 0; i < nelts; i++)
          {
            E_Int ind1 = cn1[i]-1; 
            E_Int ind2 = cn2[i]-1;
            E_Float dx = xt[ind2]-xt[ind1];
            E_Float dy = yt[ind2]-yt[ind1];
            E_Float dz = zt[ind2]-zt[ind1];
            fieldp[i] = sqrt(dx*dx + dy*dy + dz*dz);
          }
        }
        else
        {
          FldArrayF snx(nelts, nedges);
          FldArrayF sny(nelts, nedges);
          FldArrayF snz(nelts, nedges);
          FldArrayF surf(nelts, nedges);
          FldArrayF vol(nelts);
          //tableau local au fortran
          FldArrayF xint(nelts,nedges);
          FldArrayF yint(nelts,nedges);
          FldArrayF zint(nelts,nedges);
          //
          k6compunstrmetric_(npts, nelts, nedges, nnodes, cn->begin(), 
                             xt, yt, zt, 
                             xint.begin(), yint.begin(), zint.begin(),
                             snx.begin(), sny.begin(), 
                             snz.begin(), surf.begin(), vol.begin());

          if (nedges == 1) // surfacique
          {
            E_Float* surfp = surf.begin(1);
            for (E_Int i = 0; i < nelts; i++) fieldp[i] = surfp[i];
          }
          else // elements volumique
          {
            E_Float* volp = vol.begin();
            for (E_Int i = 0; i < nelts; i++) fieldp[i] = volp[i];
          }
        }
        RELEASESHAREDU(array, f, cn);
      }
      else // Elements NGON
      {
        E_Int* cnp = cn->begin(); // pointeur sur la connectivite NGon
        E_Int sizeFN = cnp[1]; //  taille de la connectivite Face/Noeuds
        E_Int nelts = cnp[sizeFN+2];  // nombre total d elements

        // Build array contenant le volume des elements
        tpl = K_ARRAY::buildArray(1, "vol",
                                  npts, nelts,
                                  -1, eltType, true, cn->getSize());
        E_Float* volp = K_ARRAY::getFieldPtr(tpl); // pointeur sur le tableau de volumes
        E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
        FldArrayI cnn(cn->getSize(), 1, cnnp, true); cnn = *cn;
        
        // compute array vol which store volume at element centers
        E_Int err = K_METRIC::CompNGonVol(xt,yt,zt,*cn,volp);

        // sortie si une erreur a ete trouvee
        if (err == 1)
        {
          PyErr_SetString(PyExc_ValueError,
                          "getVolumeMap: dimension of element unknown.");
          RELEASESHAREDU(array, f, cn);
          return NULL;
        }
        RELEASESHAREDU(array, f, cn);
      }
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
