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

// getOrthogonalityMap

# include "generator.h"
# include <math.h>

using namespace K_FLD;
using namespace K_FUNC;

// ============================================================================
/* Return orthogonality map */
/* angle is given in degree */
// Definition of the returned value: maximum of the difference between 
// the dihedral angle of the element and the dihedral angle for an 
// "ideal" element.
// ============================================================================
PyObject* K_GENERATOR::getOrthogonalityMap(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Int normalized;

  if (!PYPARSETUPLE_(args, O_ I_, &array, &normalized)) return NULL;
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int posx, posy, posz;
  E_Int res;
  res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, 
                               eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOrthogonalityMap: unknown type of array.");
    return NULL;
  }

  E_Int api = f->getApi();
  E_Int npts = f->getSize();
  PyObject* tpl = NULL;
  FldArrayF* f2;
  
  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_ValueError,
                    "getOrthogonalityMap: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  // pointeurs associes aux coordonnees
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);

  if (res == 1) // cas structure
  {
    // Dimension du tableau
    E_Int dim = 3;
    E_Int im1 = im-1;
    E_Int jm1 = jm-1;
    E_Int km1 = km-1;
    if (im == 1) im1 = 1;
    if (jm == 1) jm1 = 1;
    if (km == 1) km1 = 1;
    E_Int ni = im;
    E_Int nj = jm;
    E_Int nk = km;

    if (im == 1) // i constant
    {
      if (jm == 1 || km == 1) dim = 1;
      else
      {
        dim = 2;
        ni = jm; nj = km;
      }
    }
    else if (jm == 1) // i constant
    {
      if (im == 1 || km == 1) dim = 1;
      else
      {
        dim = 2;
        ni = im; nj = km;
      }
    } 
    else if (km == 1) // k constant
    {
      if (im == 1 || jm == 1) dim = 1;
      else
      {
        dim = 2;
        ni = im; nj = jm;
      }
    }

    E_Int ni1 = ni-1; E_Int nj1 = nj-1; E_Int nk1 = nk-1;
    if (ni == 1) ni1 = 1;
    if (nj == 1) nj1 = 1;
    if (nk == 1) nk1 = 1;
    
    // Construction du tableau numpy stockant les angles 
    // definissant l'orthogonalite
    tpl = K_ARRAY::buildArray3(1, "orthogonality", im1, jm1, km1, api);
    K_ARRAY::getFromArray3(tpl, f2);
    E_Float* skewness = f2->begin(1); // pointeur sur le tableau d'angle

    E_Int ncells = im1*jm1*km1;
    
    // calcul de l'orthogonalite
    #pragma omp parallel
    {
      E_Float skewness1, skewness2, skewness3;
      E_Int inext, jnext, knext, ind, ind1, ind2, ind3;
      if (dim == 1)
      {
        #pragma omp for
        for (E_Int i = 0; i < ncells; i++) skewness[i] = 0.;
      }
      else if (dim == 2)
      {
        #pragma omp for collapse(2)
        for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          inext = i+1;
          jnext = j+1;
          ind = j*ni1+i;
          // (i,j) angle
          ind1 = j*ni+inext;
          ind2 = j*ni+i;
          ind3 = jnext*ni+i;
          skewness[ind] = K_COMPGEOM::computeSkewness(x, y, z, ind1, ind2, ind3, 90., normalized);
        }
      }
      else // dim == 3
      {
        #pragma omp for collapse(3)  
        for (E_Int k = 0; k < nk1; k++)
        for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          inext = i+1;
          jnext = j+1;
          knext = k+1;
          ind  = k*ni1*nj1+j*ni1+i;
          // (i,j) angle
          ind1 = k*ni*nj+j*ni+inext;
          ind2 = k*ni*nj+j*ni+i;
          ind3 = k*ni*nj+jnext*ni+i;
          skewness1 = K_COMPGEOM::computeSkewness(x, y, z, ind1, ind2, ind3, 90., normalized);
          // (i,k) angle
          ind3 = knext*ni*nj+j*ni+i;
          skewness2 = K_COMPGEOM::computeSkewness(x, y, z, ind1, ind2, ind3, 90., normalized);
          // (j,k) angle
          ind1 = k*ni*nj+jnext*ni+i;
          skewness3 = K_COMPGEOM::computeSkewness(x, y, z, ind1, ind2, ind3, 90., normalized);
          // max skewness
          skewness[ind] = E_max(E_max(skewness1,skewness2),skewness3);
        }
      }
    }

    RELEASESHAREDS(tpl, f2);
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else // Cas non structure
  {
    tpl = K_ARRAY::buildArray3(
      1, "orthogonality", npts, *cn, eltType, true, api, true
    );
    K_ARRAY::getFromArray3(tpl, f2);
    E_Float* skewness = f2->begin(1); // pointeur sur le tableau d'angle

    if (strcmp(eltType, "NGON") == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "getOrthogonalityMap: not implemented for NGON arrays.");
      RELEASESHAREDS(tpl, f2);
      RELEASESHAREDU(array, f, cn);
      return NULL;
    }

    E_Int nc = cn->getNConnect();

    // Number of facets per element
    std::vector<E_Int> nfpe;
    E_Int ierr = K_CONNECT::getNFPE(nfpe, eltType, false);
    if (ierr != 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "getOrthogonalityMap: Error computing nfpe.");
      RELEASESHAREDS(tpl, f2);
      RELEASESHAREDU(array, f, cn);
      return NULL;
    }

    // get eltTypes
    std::vector<char*> eltTypes;
    K_ARRAY::extractVars(eltType, eltTypes);

    
    // calcul de l'orthogonalite
    #pragma omp parallel
    {
      E_Float skewness1, skewness2, skewness3, skewness4;
      E_Int ind, ind1, ind2, ind3, ind4;
      E_Int nelts, nvpf;
      E_Int elOffset = 0;
      std::vector<std::vector<E_Int> > facets;
      // loop over all connectivities
      for (E_Int ic = 0; ic < nc; ic++)
      {
        FldArrayI& cm = *(cn->getConnect(ic));
        nelts = cm.getSize();
        K_CONNECT::getEVFacets(facets, eltTypes[ic], false, false);
        
        // loop over all elements of connectivity cm
        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          ind = i + elOffset; // true element index
          skewness[ind] = 0.;

          // loop over all faces of element i
          for (E_Int f = 0; f < nfpe[ic]; f++)
          {
            nvpf = facets[f].size();  // number of vertices per facet
            if (nvpf == 4) // QUAD face
            {
              ind1 = cm(i,facets[f][0])-1;
              ind2 = cm(i,facets[f][1])-1;
              ind3 = cm(i,facets[f][2])-1;
              ind4 = cm(i,facets[f][3])-1;

              skewness1 = K_COMPGEOM::computeSkewness(x, y, z, ind4, ind1, ind2, 90., normalized); // angle 412
              skewness2 = K_COMPGEOM::computeSkewness(x, y, z, ind1, ind2, ind3, 90., normalized); // angle 123
              skewness3 = K_COMPGEOM::computeSkewness(x, y, z, ind2, ind3, ind4, 90., normalized); // angle 234
              skewness4 = K_COMPGEOM::computeSkewness(x, y, z, ind3, ind4, ind1, 90., normalized); // angle 341

              skewness[ind] = E_max(E_max(E_max(E_max(skewness1,skewness2),skewness3),skewness4),skewness[ind]);
            }
            else if (nvpf == 3) // TRI face
            {
              ind1 = cm(i,facets[f][0])-1;
              ind2 = cm(i,facets[f][1])-1;
              ind3 = cm(i,facets[f][2])-1;

              skewness1 = K_COMPGEOM::computeSkewness(x, y, z, ind3, ind1, ind2, 60., normalized); // angle 312
              skewness2 = K_COMPGEOM::computeSkewness(x, y, z, ind1, ind2, ind3, 60., normalized); // angle 123
              skewness3 = K_COMPGEOM::computeSkewness(x, y, z, ind2, ind3, ind1, 60., normalized); // angle 231

              skewness[ind] = E_max(E_max(E_max(skewness1,skewness2),skewness3),skewness[ind]);
            }
          }
        }
        elOffset += nelts;
      }
    }

    for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
    RELEASESHAREDS(tpl, f2);
    RELEASESHAREDU(array, f, cn); 
    return tpl;
  }
}
