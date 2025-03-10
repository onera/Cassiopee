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

// getOrthogonalityMap

# include "generator.h"

using namespace K_CONST;
using namespace K_FLD;
using namespace K_FUNC;

extern "C"
{
  void k6unstructsurf_(E_Int& npts, E_Int& nelts, E_Int& nedges, 
                       E_Int& nnodes, E_Int* cn, 
                       E_Float* coordx, E_Float* coordy, E_Float* coordz, 
                       E_Float* snx, E_Float* sny, E_Float* snz,
                       E_Float* surface);
}

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
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;
  
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int posx, posy, posz;
  E_Int res;
  res = K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, 
                              eltType, true);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOrthogonalityMap: unknown type of array.");
    return NULL;
  }

  PyObject* tpl = NULL;
  
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
  // valeur de pi pour exprimer les angles en degre
  E_Float pi = 4*atan(1.);
  E_Float degconst = 180.0 / pi;

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
    // Direction constante (en 2D)
    // 0: 3D
    // 1: 2D i constant
    // 2: 2D j constant
    // 3: 2D k constant
    E_Int dir = 0;
    if (im == 1)
    {
      dim = 2; dir=1;
      if ((jm ==1)||(km == 1)) dim = 1;
    }
    else if (jm == 1)
    {
      dim = 2; dir = 2;
      if ((im ==1)||(km == 1)) dim = 1;
    } 
    else if (km == 1)
    {
      dim = 2; dir = 3;
      if ((im ==1)||(jm == 1)) dim = 1;
    }
    
    // Construction du tableau numpy stockant les angles 
    // definissant l'orthogonalite
    tpl = K_ARRAY::buildArray(1, "orthogonality", im1, jm1, km1);
    // pointeur sur le tableau d'angle
    E_Float* alpha = K_ARRAY::getFieldPtr(tpl);
    E_Int ncells = im1*jm1*km1;
    FldArrayF alphamax(ncells, 1, alpha, true);
    
    // calcul de l'orthogonalite
    if (dim == 1)         // dimension = 1D
    {
      for (E_Int i = 0; i < ncells; i++) alphamax[i] = 0.;
    }
    else if (dim == 2) // dimension = 2D
    {
      // (ni,nj) : dimensions of 2D-grid 
      E_Int ni=1, nj=1;
      switch (dir)
      {
        case 3: // k constant
          ni = im; nj = jm;
          break;
        case 2: // j constant
          ni = im; nj = km;
          break;
        case 1: // i constant
          ni = jm; nj = km;
          break;
      }

      E_Int ni1 = ni - 1; E_Int nj1 = nj - 1;
      if (ni == 1) ni1 = 1;
      if (nj == 1) nj1 = 1;

      // Boucle sur les indices de la grille
      for (E_Int j = 0; j < nj1; j++)
        for (E_Int i = 0; i < ni1; i++)
        {
          // indices pour le calcul de l'angle definissant l'orthogonalite
          E_Int inext = i+1;
          E_Int jnext = j+1;
          // calcul de l'angle
          E_Int ind = j*ni1+i;
          E_Int ind1 = j*ni+i;
          E_Int ind2 = j*ni+inext;
          E_Int ind3 = jnext*ni+i;
          E_Float a2 = (x[ind2]-x[ind1])*(x[ind2]-x[ind1])+(y[ind2]-y[ind1])*(y[ind2]-y[ind1])+(z[ind2]-z[ind1])*(z[ind2]-z[ind1]);
          E_Float b2 = (x[ind3]-x[ind1])*(x[ind3]-x[ind1])+(y[ind3]-y[ind1])*(y[ind3]-y[ind1])+(z[ind3]-z[ind1])*(z[ind3]-z[ind1]);
          E_Float c2 = (x[ind3]-x[ind2])*(x[ind3]-x[ind2])+(y[ind3]-y[ind2])*(y[ind3]-y[ind2])+(z[ind3]-z[ind2])*(z[ind3]-z[ind2]);
          E_Float a = sqrt(a2);
          E_Float b = sqrt(b2);
          alphamax[ind] = E_abs(acos((a2+b2-c2)/(2.*a*b))*degconst - 90.);
        }
    }
    else  // dim == 3 (dimension = 3D)
    {	      
      // (ni,nj,nk) : dimensions of 3D-grid 
      E_Int ni = im;
      E_Int nj = jm;
      E_Int nk = km;
      E_Int ni1 = ni - 1; E_Int nj1 = nj - 1; E_Int nk1 = nk - 1;
      if (ni == 1) ni1 = 1;
      if (nj == 1) nj1 = 1;
      if (nk == 1) nk1 = 1;

      // angle by mesh direction
      // Boucle sur les indices de la grille
#pragma omp parallel
      {
        E_Float alpha1, alpha2, alpha3, a2, b2, c2, a, b;
        E_Int inext, jnext, knext, ind, ind1, ind2, ind3;

        for (E_Int k=0; k < nk1; k++)
        for (E_Int j=0; j < nj1; j++)
          #pragma omp for
          for (E_Int i=0; i < ni1; i++)
          {
            // indices pour le calcul des angles definissant l'orthogonalite
            inext = i+1;
            jnext = j+1;
            knext = k+1;
            // calcul des angles
            ind  = k*ni1*nj1+j*ni1+i;
            ind1 = k*ni*nj+j*ni+i;
            ind2 = k*ni*nj+j*ni+inext;
            ind3 = k*ni*nj+jnext*ni+i;
            a2 = (x[ind2]-x[ind1])*(x[ind2]-x[ind1])+(y[ind2]-y[ind1])*(y[ind2]-y[ind1])+(z[ind2]-z[ind1])*(z[ind2]-z[ind1]);
            b2 = (x[ind3]-x[ind1])*(x[ind3]-x[ind1])+(y[ind3]-y[ind1])*(y[ind3]-y[ind1])+(z[ind3]-z[ind1])*(z[ind3]-z[ind1]);
            c2 = (x[ind3]-x[ind2])*(x[ind3]-x[ind2])+(y[ind3]-y[ind2])*(y[ind3]-y[ind2])+(z[ind3]-z[ind2])*(z[ind3]-z[ind2]);
            a = sqrt(a2);
            b = sqrt(b2);
            // ... angle correspondant aux indices (ij)
            alpha1 = E_abs(acos((a2+b2-c2)/(2.*a*b))*degconst - 90.);
            // ... angle correspondant aux indices (ik)
            ind3 = knext*ni*nj+j*ni+i;
            b2 = (x[ind3]-x[ind1])*(x[ind3]-x[ind1])+(y[ind3]-y[ind1])*(y[ind3]-y[ind1])+(z[ind3]-z[ind1])*(z[ind3]-z[ind1]);
            c2 = (x[ind3]-x[ind2])*(x[ind3]-x[ind2])+(y[ind3]-y[ind2])*(y[ind3]-y[ind2])+(z[ind3]-z[ind2])*(z[ind3]-z[ind2]);
            b = sqrt(b2);
            alpha2 = E_abs(acos((a2+b2-c2)/(2.*a*b))*degconst - 90.);
            // ... angle correspondant aux indices (jk)
            ind2 = k*ni*nj+jnext*ni+i;
            a2 = (x[ind2]-x[ind1])*(x[ind2]-x[ind1])+(y[ind2]-y[ind1])*(y[ind2]-y[ind1])+(z[ind2]-z[ind1])*(z[ind2]-z[ind1]);
            c2 = (x[ind3]-x[ind2])*(x[ind3]-x[ind2])+(y[ind3]-y[ind2])*(y[ind3]-y[ind2])+(z[ind3]-z[ind2])*(z[ind3]-z[ind2]);
            a = sqrt(a2);
            alpha3 = E_abs(acos((a2+b2-c2)/(2.*a*b))*degconst - 90.);
            alphamax[ind] = E_max(E_max(alpha1,alpha2),alpha3);
          }
      }
    }
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else // if (res == 2)
  {
    // Cas non structure
    E_Int nelts = cn->getSize();
    E_Int nnodes = cn->getNfld(); // nb de noeuds ds 1 element
    E_Int npts = f->getSize();

    E_Float* x = f->begin(posx);
    E_Float* y = f->begin(posy);
    E_Float* z = f->begin(posz);

    PyObject* tpl = K_ARRAY::buildArray(1, "orthogonality", nelts, 
					nelts, -1, eltType, true);
    E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
    K_KCORE::memcpy__(cnnp, cn->begin(), nelts*nnodes);
    E_Float* alphamaxp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF alphamax(nelts,1, alphamaxp, true);
 
    if (strcmp(eltType, "TRI") == 0)
    {
      E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2); E_Int* cn3 = cn->begin(3);
      E_Float eps = 1.e-6;
      #pragma omp parallel
      {
        E_Float d1, d2, d3, alpha12, alpha23, alpha31, sqrtd1, sqrtd2, sqrtd3;
        E_Int ind1, ind2, ind3;
        #pragma omp for
        for (E_Int et = 0; et < nelts; et++)
        {
          ind1 = cn1[et]-1;
          ind2 = cn2[et]-1;
          ind3 = cn3[et]-1;

          d1 = (x[ind2]-x[ind1])*(x[ind2]-x[ind1])+(y[ind2]-y[ind1])*(y[ind2]-y[ind1])+(z[ind2]-z[ind1])*(z[ind2]-z[ind1]); sqrtd1 = sqrt(d1);
          d2 = (x[ind3]-x[ind2])*(x[ind3]-x[ind2])+(y[ind3]-y[ind2])*(y[ind3]-y[ind2])+(z[ind3]-z[ind2])*(z[ind3]-z[ind2]); sqrtd2 = sqrt(d2);
          d3 = (x[ind1]-x[ind3])*(x[ind1]-x[ind3])+(y[ind1]-y[ind3])*(y[ind1]-y[ind3])+(z[ind1]-z[ind3])*(z[ind1]-z[ind3]); sqrtd3 = sqrt(d3);

          if (( K_FUNC::fEqualZero(d1,eps) == true )||
              ( K_FUNC::fEqualZero(d2,eps) == true )||
              ( K_FUNC::fEqualZero(d3,eps) == true ))
          {
            // degenerated element
            alphamax[et]= 120.; // 180. - 60.
          }
          else
          {
            alpha12 = E_abs(acos((d1 + d2 - d3)/(E_max(TWO*sqrtd1*sqrtd2,E_GEOM_CUTOFF)))*degconst - 60.);
            alpha23 = E_abs(acos((d2 + d3 - d1)/(E_max(TWO*sqrtd2*sqrtd3,E_GEOM_CUTOFF)))*degconst - 60.);
            alpha31 = E_abs(acos((d3 + d1 - d2)/(E_max(TWO*sqrtd3*sqrtd1,E_GEOM_CUTOFF)))*degconst - 60.);
            alphamax[et] = E_max(E_max(alpha12,alpha23),alpha31);
          }
        }
      }
    }
    else if (strcmp(eltType, "QUAD") == 0)
    {
      E_Int* cn1 = cn->begin(1); E_Int* cn2 = cn->begin(2); E_Int* cn3 = cn->begin(3); E_Int* cn4 = cn->begin(4);
      E_Int ind1, ind2, ind3, ind4;
      E_Float d1, d2, d3, d4, diag13, diag24, alpha12, alpha14, alpha32, alpha34;
      E_Float sqrtd1, sqrtd2, sqrtd3, sqrtd4;
      #pragma omp for
      for (E_Int et = 0; et < nelts; et++)
      {
        ind1 = cn1[et]-1;
        ind2 = cn2[et]-1;
        ind3 = cn3[et]-1;
        ind4 = cn4[et]-1;

        d1     = (x[ind2]-x[ind1])*(x[ind2]-x[ind1])+(y[ind2]-y[ind1])*(y[ind2]-y[ind1])+(z[ind2]-z[ind1])*(z[ind2]-z[ind1]); sqrtd1 = sqrt(d1);
        d2     = (x[ind3]-x[ind2])*(x[ind3]-x[ind2])+(y[ind3]-y[ind2])*(y[ind3]-y[ind2])+(z[ind3]-z[ind2])*(z[ind3]-z[ind2]); sqrtd2 = sqrt(d2);
        d3     = (x[ind4]-x[ind3])*(x[ind4]-x[ind3])+(y[ind4]-y[ind3])*(y[ind4]-y[ind3])+(z[ind4]-z[ind3])*(z[ind4]-z[ind3]); sqrtd3 = sqrt(d3);
        d4     = (x[ind1]-x[ind4])*(x[ind1]-x[ind4])+(y[ind1]-y[ind4])*(y[ind1]-y[ind4])+(z[ind1]-z[ind4])*(z[ind1]-z[ind4]); sqrtd4 = sqrt(d4);
        diag13 = (x[ind1]-x[ind3])*(x[ind1]-x[ind3])+(y[ind1]-y[ind3])*(y[ind1]-y[ind3])+(z[ind1]-z[ind3])*(z[ind1]-z[ind3]);
        diag24 = (x[ind2]-x[ind4])*(x[ind2]-x[ind4])+(y[ind2]-y[ind4])*(y[ind2]-y[ind4])+(z[ind2]-z[ind4])*(z[ind2]-z[ind4]);

        alpha12 = E_abs(acos((d1 + d2 - diag13)/(E_max(TWO*sqrtd1*sqrtd2,E_GEOM_CUTOFF)))*degconst - 90.);
        alpha14 = E_abs(acos((d1 + d4 - diag24)/(E_max(TWO*sqrtd1*sqrtd4,E_GEOM_CUTOFF)))*degconst - 90.);
        alpha32 = E_abs(acos((d3 + d2 - diag24)/(E_max(TWO*sqrtd3*sqrtd2,E_GEOM_CUTOFF)))*degconst - 90.);
        alpha34 = E_abs(acos((d3 + d4 - diag13)/(E_max(TWO*sqrtd3*sqrtd4,E_GEOM_CUTOFF)))*degconst - 90.);
        alphamax[et] = E_max(E_max(alpha12,alpha14),E_max(alpha32,alpha34));
      }
    }
    else if (strcmp(eltType, "TETRA") == 0)
    {
      // Compute surface normals
      E_Int nedges = 4;
      FldArrayF nsurfx(nelts, nedges);
      FldArrayF nsurfy(nelts, nedges);
      FldArrayF nsurfz(nelts, nedges);
      FldArrayF surf(nelts, nedges);
      k6unstructsurf_(npts, nelts, nedges, nnodes, cn->begin(), 
                      f->begin(posx), f->begin(posy), f->begin(posz), 
                      nsurfx.begin(), nsurfy.begin(), nsurfz.begin(), 
                      surf.begin());
      // Compute dihedral angle
      E_Float* nsurf1x = nsurfx.begin(1); E_Float* nsurf2x = nsurfx.begin(2); E_Float* nsurf3x = nsurfx.begin(3);E_Float* nsurf4x = nsurfx.begin(4);
      E_Float* nsurf1y = nsurfy.begin(1); E_Float* nsurf2y = nsurfy.begin(2); E_Float* nsurf3y = nsurfy.begin(3);E_Float* nsurf4y = nsurfy.begin(4);
      E_Float* nsurf1z = nsurfz.begin(1); E_Float* nsurf2z = nsurfz.begin(2); E_Float* nsurf3z = nsurfz.begin(3);E_Float* nsurf4z = nsurfz.begin(4);
      E_Float* surf1 = surf.begin(1); E_Float* surf2 = surf.begin(2); 
      E_Float* surf3 = surf.begin(3); E_Float* surf4 = surf.begin(4);
    
      #pragma omp parallel
      {
        E_Float s1, s2, s3,s4, s12, s13, s23, s14, s24, s34;
        E_Float alpha12, alpha14, alpha24, alpha13, alpha23, alpha34;
        E_Float alphamax1, alphamax2, alphamax3;
        #pragma omp for
        for (E_Int et = 0; et < nelts; et++)
        {
          s1 = surf1[et]; s2 = surf2[et]; s3 = surf3[et]; s4 = surf4[et];
          s12 = s1*s2; s13 = s1*s3; s14 = s1*s4;  
          s23 = s2*s3; s24 = s2*s4; s34 = s3*s4;
          alpha12 = E_abs(acos((nsurf1x[et]*nsurf2x[et] + nsurf1y[et]*nsurf2y[et] + nsurf1z[et]*nsurf2z[et])/s12)*degconst - 120.);
          alpha14 = E_abs(acos((nsurf1x[et]*nsurf4x[et] + nsurf1y[et]*nsurf4y[et] + nsurf1z[et]*nsurf4z[et])/s14)*degconst - 120.);
          alpha24 = E_abs(acos((nsurf2x[et]*nsurf4x[et] + nsurf2y[et]*nsurf4y[et] + nsurf2z[et]*nsurf4z[et])/s24)*degconst - 120.);
          alpha13 = E_abs(acos((nsurf1x[et]*nsurf3x[et] + nsurf1y[et]*nsurf3y[et] + nsurf1z[et]*nsurf3z[et])/s13)*degconst - 120.);
          alpha23 = E_abs(acos((nsurf2x[et]*nsurf3x[et] + nsurf2y[et]*nsurf3y[et] + nsurf2z[et]*nsurf3z[et])/s23)*degconst - 120.);
          alpha34 = E_abs(acos((nsurf3x[et]*nsurf4x[et] + nsurf3y[et]*nsurf4y[et] + nsurf3z[et]*nsurf4z[et])/s34)*degconst - 120.);
          alphamax1 = E_max(alpha12,alpha14);
          alphamax2 = E_max(alpha24,alpha13);
          alphamax3 = E_max(alpha23,alpha34);
          alphamax[et] = E_max(E_max(alphamax1,alphamax2),alphamax3);
        }
      }
    }
    else if (strcmp(eltType, "PYRA") == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                    "getOrthogonalityMap: not yet implemented for PYRA.");
      return NULL;
    }
    else if (strcmp(eltType, "PENTA") == 0)
    {
      // Compute surface normals
      E_Int nedges = 5;
      FldArrayF nsurfx(nelts, nedges);
      FldArrayF nsurfy(nelts, nedges);
      FldArrayF nsurfz(nelts, nedges);
      FldArrayF surf(nelts, nedges);
      k6unstructsurf_(npts, nelts, nedges, nnodes, cn->begin(), 
                      f->begin(posx), f->begin(posy), f->begin(posz), 
                      nsurfx.begin(), nsurfy.begin(), nsurfz.begin(), 
                      surf.begin());
      // Compute dihedral angle
      E_Float* nsurf1x = nsurfx.begin(1); E_Float* nsurf2x = nsurfx.begin(2); E_Float* nsurf3x = nsurfx.begin(3);
      E_Float* nsurf4x = nsurfx.begin(4); E_Float* nsurf5x = nsurfx.begin(5);
      E_Float* nsurf1y = nsurfy.begin(1); E_Float* nsurf2y = nsurfy.begin(2); E_Float* nsurf3y = nsurfy.begin(3);
      E_Float* nsurf4y = nsurfy.begin(4); E_Float* nsurf5y = nsurfy.begin(5);
      E_Float* nsurf1z = nsurfz.begin(1); E_Float* nsurf2z = nsurfz.begin(2); E_Float* nsurf3z = nsurfz.begin(3);
      E_Float* nsurf4z = nsurfz.begin(4); E_Float* nsurf5z = nsurfz.begin(5);
      E_Float* surf1 = surf.begin(1); E_Float* surf2 = surf.begin(2); 
      E_Float* surf3 = surf.begin(3); E_Float* surf4 = surf.begin(4); E_Float* surf5 = surf.begin(5);
      
      #pragma omp parallel
      {
        E_Float s1, s2, s3, s4, s5, s13, s14, s15, s23, s24, s25, s34, s35, s45;
        E_Float alpha13, alpha14, alpha15, alpha23, alpha24, alpha25, alpha34, alpha35, alpha45;
        E_Float alphamax1, alphamax2, alphamax3;  
        #pragma omp for
        for (E_Int et = 0; et < nelts; et++)
        {
          s1 = surf1[et]; s2 = surf2[et]; s3 = surf3[et]; s4 = surf4[et]; s5 = surf5[et];
          s13 = s1*s3; s14 = s1*s4; s15 = s1*s5;  
          s23 = s2*s3; s24 = s2*s4; s25 = s2*s5;  
          s34 = s3*s4; s35 = s3*s5; s45 = s4*s5;  
          alpha13 = E_abs(acos((nsurf1x[et]*nsurf3x[et] + nsurf1y[et]*nsurf3y[et] + nsurf1z[et]*nsurf3z[et])/s13)*degconst - 90.);
          alpha14 = E_abs(acos((nsurf1x[et]*nsurf4x[et] + nsurf1y[et]*nsurf4y[et] + nsurf1z[et]*nsurf4z[et])/s14)*degconst - 90.);
          alpha15 = E_abs(acos((nsurf1x[et]*nsurf5x[et] + nsurf1y[et]*nsurf5y[et] + nsurf1z[et]*nsurf5z[et])/s15)*degconst - 90.);
          alpha23 = E_abs(acos((nsurf2x[et]*nsurf3x[et] + nsurf2y[et]*nsurf3y[et] + nsurf2z[et]*nsurf3z[et])/s23)*degconst - 90.);
          alpha24 = E_abs(acos((nsurf2x[et]*nsurf4x[et] + nsurf2y[et]*nsurf4y[et] + nsurf2z[et]*nsurf4z[et])/s24)*degconst - 90.);
          alpha25 = E_abs(acos((nsurf2x[et]*nsurf5x[et] + nsurf2y[et]*nsurf5y[et] + nsurf2z[et]*nsurf5z[et])/s25)*degconst - 90.);
          alpha34 = E_abs(acos((nsurf3x[et]*nsurf4x[et] + nsurf3y[et]*nsurf4y[et] + nsurf3z[et]*nsurf4z[et])/s34)*degconst - 120.);
          alpha35 = E_abs(acos((nsurf3x[et]*nsurf5x[et] + nsurf3y[et]*nsurf5y[et] + nsurf3z[et]*nsurf5z[et])/s35)*degconst - 120.);
          alpha45 = E_abs(acos((nsurf4x[et]*nsurf5x[et] + nsurf4y[et]*nsurf5y[et] + nsurf4z[et]*nsurf5z[et])/s45)*degconst - 120.);
          alphamax1 = E_max(E_max(alpha13,alpha14),alpha15);
          alphamax2 = E_max(E_max(alpha23,alpha24),alpha25);
          alphamax3 = E_max(E_max(alpha34,alpha35),alpha45);
          alphamax[et] = E_max(E_max(alphamax1,alphamax2),alphamax3);
        }
      }
    }
    else if (strcmp(eltType, "HEXA") == 0)
    {
      // Compute surface normals
      E_Int nedges = 6;
      FldArrayF nsurfx(nelts, nedges);
      FldArrayF nsurfy(nelts, nedges);
      FldArrayF nsurfz(nelts, nedges);
      FldArrayF surf(nelts, nedges);
      k6unstructsurf_(npts, nelts, nedges, nnodes, cn->begin(), 
                      f->begin(posx), f->begin(posy), f->begin(posz), 
                      nsurfx.begin(), nsurfy.begin(), nsurfz.begin(), 
                      surf.begin());
      // Compute dihedral angle
      E_Float* nsurf1x = nsurfx.begin(1); E_Float* nsurf2x = nsurfx.begin(2); E_Float* nsurf3x = nsurfx.begin(3);
      E_Float* nsurf4x = nsurfx.begin(4); E_Float* nsurf5x = nsurfx.begin(5); E_Float* nsurf6x = nsurfx.begin(6);
      E_Float* nsurf1y = nsurfy.begin(1); E_Float* nsurf2y = nsurfy.begin(2); E_Float* nsurf3y = nsurfy.begin(3);
      E_Float* nsurf4y = nsurfy.begin(4); E_Float* nsurf5y = nsurfy.begin(5); E_Float* nsurf6y = nsurfy.begin(6);
      E_Float* nsurf1z = nsurfz.begin(1); E_Float* nsurf2z = nsurfz.begin(2); E_Float* nsurf3z = nsurfz.begin(3);
      E_Float* nsurf4z = nsurfz.begin(4); E_Float* nsurf5z = nsurfz.begin(5); E_Float* nsurf6z = nsurfz.begin(6);
      E_Float* surf1 = surf.begin(1); E_Float* surf2 = surf.begin(2); E_Float* surf3 = surf.begin(3);
      E_Float* surf4 = surf.begin(4); E_Float* surf5 = surf.begin(5); E_Float* surf6 = surf.begin(6);
      
      #pragma omp parallel
      {
        E_Float s1, s2, s3, s4, s5, s6;
        E_Float s13, s14, s15, s16, s23, s24, s25, s26, s35, s36, s45, s46;
        E_Float alpha13, alpha14, alpha15, alpha16, alpha23, alpha24, alpha25, alpha26;
        E_Float alpha35, alpha36, alpha45, alpha46;
        E_Float alphamax1, alphamax2, alphamax3, alphamax4;
        
        #pragma omp for
        for (E_Int et = 0; et < nelts; et++)
        {
          s1 = surf1[et]; s2 = surf2[et]; s3 = surf3[et];
          s4 = surf4[et]; s5 = surf5[et]; s6 = surf6[et];
          s13 = s1*s3; s14 = s1*s4; s15 = s1*s5; s16 = s1*s6; 
          s23 = s2*s3; s24 = s2*s4; s25 = s2*s5; s26 = s2*s6; 
          s35 = s3*s5; s36 = s3*s6; s45 = s4*s5; s46 = s4*s6; 
          alpha13 = E_abs(acos((nsurf1x[et]*nsurf3x[et] + nsurf1y[et]*nsurf3y[et] + nsurf1z[et]*nsurf3z[et])/s13)*degconst - 90.);
          alpha14 = E_abs(acos((nsurf1x[et]*nsurf4x[et] + nsurf1y[et]*nsurf4y[et] + nsurf1z[et]*nsurf4z[et])/s14)*degconst - 90.);
          alpha15 = E_abs(acos((nsurf1x[et]*nsurf5x[et] + nsurf1y[et]*nsurf5y[et] + nsurf1z[et]*nsurf5z[et])/s15)*degconst - 90.);
          alpha16 = E_abs(acos((nsurf1x[et]*nsurf6x[et] + nsurf1y[et]*nsurf6y[et] + nsurf1z[et]*nsurf6z[et])/s16)*degconst - 90.);
          alpha23 = E_abs(acos((nsurf2x[et]*nsurf3x[et] + nsurf2y[et]*nsurf3y[et] + nsurf2z[et]*nsurf3z[et])/s23)*degconst - 90.);
          alpha24 = E_abs(acos((nsurf2x[et]*nsurf4x[et] + nsurf2y[et]*nsurf4y[et] + nsurf2z[et]*nsurf4z[et])/s24)*degconst - 90.);
          alpha25 = E_abs(acos((nsurf2x[et]*nsurf5x[et] + nsurf2y[et]*nsurf5y[et] + nsurf2z[et]*nsurf5z[et])/s25)*degconst - 90.);
          alpha26 = E_abs(acos((nsurf2x[et]*nsurf6x[et] + nsurf2y[et]*nsurf6y[et] + nsurf2z[et]*nsurf6z[et])/s26)*degconst - 90.);
          alpha35 = E_abs(acos((nsurf3x[et]*nsurf5x[et] + nsurf3y[et]*nsurf5y[et] + nsurf3z[et]*nsurf5z[et])/s35)*degconst - 90.);
          alpha36 = E_abs(acos((nsurf3x[et]*nsurf6x[et] + nsurf3y[et]*nsurf6y[et] + nsurf3z[et]*nsurf6z[et])/s36)*degconst - 90.);
          alpha45 = E_abs(acos((nsurf4x[et]*nsurf5x[et] + nsurf4y[et]*nsurf5y[et] + nsurf4z[et]*nsurf5z[et])/s45)*degconst - 90.);
          alpha46 = E_abs(acos((nsurf4x[et]*nsurf6x[et] + nsurf4y[et]*nsurf6y[et] + nsurf4z[et]*nsurf6z[et])/s46)*degconst - 90.);
          alphamax1 = E_max(E_max(alpha13,alpha14),alpha15);
          alphamax2 = E_max(E_max(alpha16,alpha23),alpha24);
          alphamax3 = E_max(E_max(alpha25,alpha26),alpha35);
          alphamax4 = E_max(E_max(alpha36,alpha45),alpha46);
          alphamax[et] = E_max(E_max(alphamax1,alphamax2),E_max(alphamax3,alphamax4));
        }
      }
    }
    else if (strcmp(eltType, "BAR") == 0)
    {
      for (E_Int et = 0; et < nelts; et++)
      {
        alphamax[et] = 0.;
      }
    }
    else
    {
      PyErr_SetString(PyExc_TypeError,
                      "getOrthogonalityMap: unknown type of element.");
      RELEASESHAREDU(array, f, cn);
      return NULL;
    }
    RELEASESHAREDU(array, f, cn); 
    return tpl;
  }
}
