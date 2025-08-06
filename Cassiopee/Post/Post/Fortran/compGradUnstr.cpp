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
# include "Array/Array.h"
# include <math.h>
# include <vector>

/* Calcul du gradient d'un champ defini aux noeuds d une grille non structuree
   retourne le gradient defini aux centres des elts.
   Cas 1D, 2D et 3D.
*/
E_Int K_POST::computeGradUnstr(
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  K_FLD::FldArrayI& cn, const char* eltType,
  E_Float* field,
  E_Float* gradx, E_Float* grady, E_Float* gradz
)
{
  // Get ME mesh dimensionality from the first element type
  E_Int dim = 3;
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);
  if (strcmp(eltTypes[0], "BAR") == 0) dim = 1;
  else if (strcmp(eltTypes[0], "TRI") == 0 or
           strcmp(eltTypes[0], "QUAD") == 0) dim = 2;

  E_Int nelts = 0; // TODO
  E_Int nedges = 0;

  if (dim == 3)
  {
    K_FLD::FldArrayF fieldf(nelts, nedges);
    K_FLD::FldArrayF snx(nelts, nedges);
    K_FLD::FldArrayF sny(nelts, nedges);
    K_FLD::FldArrayF snz(nelts, nedges);
    K_FLD::FldArrayF surf(nelts, nedges);
    K_FLD::FldArrayF vol(nelts);
    K_FLD::FldArrayF xint(nelts, nedges);
    K_FLD::FldArrayF yint(nelts, nedges);
    K_FLD::FldArrayF zint(nelts, nedges);
    compGradUnstr3D(
      xt, yt, zt, cn, eltType,
      field, fieldf.begin(),
      snx.begin(), sny.begin(), snz.begin(), 
      surf.begin(), vol.begin(), 
      xint.begin(), yint.begin(), zint.begin(),
      gradx, grady, gradz
    );
  }
  else if (dim == 2)
  {
    K_FLD::FldArrayF snx(nelts, 1);
    K_FLD::FldArrayF sny(nelts, 1);
    K_FLD::FldArrayF snz(nelts, 1);
    K_FLD::FldArrayF surf(nelts, 1);
    compGradUnstr2D(
      xt, yt, zt, cn, eltType,
      field,
      snx.begin(), sny.begin(), snz.begin(), surf.begin(),
      gradx, grady, gradz
    );
  }
  else // dim = 1
  {
    compGradUnstr1D(
      xt, yt, zt, cn, eltType,
      field, gradx, grady, gradz);    
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
  return 1;
}

/* Calcul du gradient d'un champ defini aux noeuds d une grille non structuree
   retourne le gradient defini aux centres des elts.
   CAS 1D.
   IN: xt, yt, zt: coordonnees x, y, z des pts de la grille
   IN: cn: connectivite elts-noeuds
   IN: eltType: list of BE element types forming the ME mesh
   IN: field: champ defini aux noeuds auquel on applique grad
   OUT: gradx, grady, gradz: gradient de field %x, %y, %z
*/
void K_POST::compGradUnstr1D(
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  K_FLD::FldArrayI& cn, const char* eltType, const E_Float* field,
  E_Float* gradx, E_Float* grady, E_Float* gradz
)
{
  E_Float eps = 1e-12;
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Pre-compute element offsets
  std::vector<E_Int> nepc(nc+1); nepc[0] = 0;
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    if (strcmp(eltTypes[ic], "BAR") != 0)
    {
      fprintf(stderr, "Error: in K_POST::compGradUnstr1D.\n");
      fprintf(stderr, "Element type must be BAR.\n");
      exit(0);
    }
    nepc[ic+1] = nepc[ic] + nelts;
  }

  #pragma omp parallel
  {
    E_Int indA, indB, pos;
    E_Float xAB, yAB, zAB, lengthAB2, fAB;
      
    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
      E_Int nelts = cm.getSize();
      E_Int elOffset = nepc[ic];

      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        indA = cm(i, 1) - 1;
        indB = cm(i, 2) - 1;
                
        xAB = xt[indB] - xt[indA];
        yAB = yt[indB] - yt[indA];
        zAB = zt[indB] - zt[indA];
        lengthAB2 = xAB * xAB + yAB * yAB + zAB * zAB;

        fAB = field[indB] - field[indA];

        pos = elOffset + i;
        gradx[pos] = fAB * xAB / K_FUNC::E_max(lengthAB2, eps);
        grady[pos] = fAB * yAB / K_FUNC::E_max(lengthAB2, eps);
        gradz[pos] = fAB * zAB / K_FUNC::E_max(lengthAB2, eps);
      }
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

/* Calcul du gradient d'un champ defini aux noeuds d une grille non structuree
   retourne le gradient defini aux centres des elts.
   CAS 2D.
   IN: xt, yt, zt: coordonnees x, y, z des pts de la grille
   IN: cn: connectivite elts-noeuds
   IN: eltType: list of BE element types forming the ME mesh
   IN: field: champ defini aux noeuds auquel on applique grad
   OUT: snx, sny, snz: normales aux facettes %x, %y, %z
   OUT: surf: aires des facettes
   OUT: gradx, grady, gradz: gradient de field %x, %y, %z
*/
void K_POST::compGradUnstr2D(
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  K_FLD::FldArrayI& cn, const char* eltType, const E_Float* field,
  E_Float* snx, E_Float* sny, E_Float* snz, E_Float* surf,
  E_Float* gradx, E_Float* grady, E_Float* gradz
)
{
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Pre-compute element offsets
  std::vector<E_Int> nepc(nc+1); nepc[0] = 0;
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    if (strcmp(eltTypes[ic], "TRI") != 0 or strcmp(eltTypes[ic], "QUAD") != 0)
    {
      fprintf(stderr, "Error: in K_POST::compGradUnstr2D.\n");
      fprintf(stderr, "Element type must be TRI or QUAD.\n");
      exit(0);
    }
    nepc[ic+1] = nepc[ic] + nelts;
  }
  
  // Compute surface of elements
  K_METRIC::compUnstructSurf(cn, eltType, xt, yt, zt,
    snx, sny, snz, surf);

  #pragma omp parallel
  {
    E_Int indA, indB, indC, indD, pos;
    E_Float nn, vinv;
    E_Float xAB, yAB, zAB, xBC, yBC, zBC, xCD, yCD, zCD, xDA, yDA, zDA;
    E_Float xCA, yCA, zCA;
    E_Float nx, n1x, n2x, n3x, n4x;
    E_Float ny, n1y, n2y, n3y, n4y;
    E_Float nz, n1z, n2z, n3z, n4z;
    E_Float fAB, fBC, fCD, fDA, fCA;
      
    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
      E_Int nelts = cm.getSize();
      E_Int elOffset = nepc[ic];

      if (strcmp(eltTypes[ic], "TRI") == 0)
      {
        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          indA = cm(i, 1) - 1;
          indB = cm(i, 2) - 1;
          indC = cm(i, 3) - 1;
          pos = elOffset + i;

          nx = snx[pos];
          ny = sny[pos];
          nz = snz[pos];
          nn = sqrt(nx * nx + ny * ny + nz * nz);
          vinv = K_CONST::TWO * surf[pos] * nn;
          vinv = K_CONST::ONE / K_FUNC::E_max(vinv, K_CONST::E_MIN_VOL);

          xAB = xt[indB] - xt[indA];
          yAB = yt[indB] - yt[indA];
          zAB = zt[indB] - zt[indA];

          xBC = xt[indC] - xt[indB];
          yBC = yt[indC] - yt[indB];
          zBC = zt[indC] - zt[indB];

          xCA = xt[indA] - xt[indC];
          yCA = yt[indA] - yt[indC];
          zCA = zt[indA] - zt[indC];

          n1x = yAB * nz - zAB * ny;
          n1y = zAB * nx - xAB * nz;
          n1z = xAB * ny - yAB * nx;

          n2x = yBC * nz - zBC * ny;
          n2y = zBC * nx - xBC * nz;
          n2z = xBC * ny - yBC * nx;

          n3x = yCA * nz - zCA * ny;
          n3y = zCA * nx - xCA * nz;
          n3z = xCA * ny - yCA * nx;

          fAB = field[indA] + field[indB];
          fBC = field[indB] + field[indC];
          fCA = field[indC] + field[indA];

          gradx[pos] = vinv * (fAB * n1x + fBC * n2x + fCA * n3x);
          grady[pos] = vinv * (fAB * n1y + fBC * n2y + fCA * n3y);
          gradz[pos] = vinv * (fAB * n1z + fBC * n2z + fCA * n3z);
        }
      }
      else  // QUAD
      {
        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          indA = cm(i, 1) - 1;
          indB = cm(i, 2) - 1;
          indC = cm(i, 3) - 1;
          indD = cm(i, 4) - 1;
          pos = elOffset + i;

          nx = snx[pos];
          ny = sny[pos];
          nz = snz[pos];

          nn = sqrt(nx * nx + ny * ny + nz * nz);
          vinv = K_CONST::TWO * surf[pos] * nn;
          vinv = K_CONST::ONE / K_FUNC::E_max(vinv, K_CONST::E_MIN_VOL);

          xAB = xt[indB] - xt[indA];
          yAB = yt[indB] - yt[indA];
          zAB = zt[indB] - zt[indA];

          xBC = xt[indC] - xt[indB];
          yBC = yt[indC] - yt[indB];
          zBC = zt[indC] - zt[indB];

          xCD = xt[indD] - xt[indC];
          yCD = yt[indD] - yt[indC];
          zCD = zt[indD] - zt[indC];

          xDA = xt[indA] - xt[indD];
          yDA = yt[indA] - yt[indD];
          zDA = zt[indA] - zt[indD];

          n1x = yAB * nz - zAB * ny;
          n1y = zAB * nx - xAB * nz;
          n1z = xAB * ny - yAB * nx;

          n2x = yBC * nz - zBC * ny;
          n2y = zBC * nx - xBC * nz;
          n2z = xBC * ny - yBC * nx;

          n3x = yCD * nz - zCD * ny;
          n3y = zCD * nx - xCD * nz;
          n3z = xCD * ny - yCD * nx;

          n4x = yDA * nz - zDA * ny;
          n4y = zDA * nx - xDA * nz;
          n4z = xDA * ny - yDA * nx;

          fAB = field[indA] + field[indB];
          fBC = field[indB] + field[indC];
          fCD = field[indC] + field[indD];
          fDA = field[indD] + field[indA];
          
          gradx[pos] = vinv * (fAB * n1x + fBC * n2x + fCD * n3x + fDA * n4x);
          grady[pos] = vinv * (fAB * n1y + fBC * n2y + fCD * n3y + fDA * n4y);
          gradz[pos] = vinv * (fAB * n1z + fBC * n2z + fCD * n3z + fDA * n4z);
        }
      }
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

/* Calcul du gradient d'un champ defini aux noeuds d une grille non structuree
   retourne le gradient defini aux centres des elts.
   CAS 3D.
   IN: xt, yt, zt: coordonnees x, y, z des pts de la grille
   IN: cn: connectivite elts-noeuds
   IN: eltType: list of BE element types forming the ME mesh
   IN: field: champ defini aux noeuds auquel on applique grad
   OUT: fieldf: ieme champ defini aux facettes des elts
   OUT: snx, sny, snz: normales aux facettes %x, %y, %z
   OUT: surf: aires des facettes
   OUT: vol: Volume of the elements
   OUT: xint, yint, zint: Coordonnees du centre des facettes
   OUT: gradx, grady, gradz: gradient de field %x, %y, %z
*/
void K_POST::compGradUnstr3D(
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  K_FLD::FldArrayI& cn, const char* eltType, const E_Float* field,
  E_Float* fieldf,
  E_Float* snx, E_Float* sny, E_Float* snz, E_Float* surf, E_Float* vol,
  E_Float* xint, E_Float* yint, E_Float* zint,
  E_Float* gradx, E_Float* grady, E_Float* gradz
)
{
  // Compute surface and volume of elements
  K_METRIC::compUnstructMetric(
    cn, eltType, xt, yt, zt, 
    xint, yint, zint, snx, sny, snz, surf, vol
  );

  // Compute field on face centers
  compUnstrNodes2Faces(cn, eltType, field, fieldf);

  // Compute gradient at element centers
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);
  
  // Pre-compute element and facet offsets
  std::vector<E_Int> nepc(nc+1), nfpc(nc+1);
  nepc[0] = 0; nfpc[0] = 0;
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    E_Int nfpe;

    if (strcmp(eltTypes[ic], "TRI") == 0) nfpe = 1;
    else if (strcmp(eltTypes[ic], "QUAD") == 0) nfpe = 1;
    else if (strcmp(eltTypes[ic], "TETRA") == 0) nfpe = 4;
    else if (strcmp(eltTypes[ic], "PYRA") == 0) nfpe = 5;
    else if (strcmp(eltTypes[ic], "PENTA") == 0) nfpe = 5;
    else if (strcmp(eltTypes[ic], "HEXA") == 0) nfpe = 6;

    nepc[ic+1] = nepc[ic] + nelts;
    nfpc[ic+1] = nfpc[ic] + nfpe*nelts;
  }
  
  #pragma omp parallel
  {
    E_Int pos, nelts, nfpe, elOffset, fctOffset;
    E_Float sumx, sumy, sumz, invvol;
    
    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
      nelts = cm.getSize();
      nfpe = (nfpc[ic+1] - nfpc[ic])/nelts;  // number of facets per element
      elOffset = nepc[ic];
      fctOffset = nfpc[ic];

      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        sumx = K_CONST::E_ZERO_FLOAT;
        sumy = K_CONST::E_ZERO_FLOAT;
        sumz = K_CONST::E_ZERO_FLOAT;
        invvol = K_CONST::ONE / K_FUNC::E_max(vol[i], K_CONST::E_MIN_VOL);

        for (E_Int fi = 0; fi < nfpe; fi++)
        {
          pos = fctOffset + i * nfpe + fi;
          sumx += snx[pos] * fieldf[pos];
          sumy += sny[pos] * fieldf[pos];
          sumz += snz[pos] * fieldf[pos];
        }

        pos = elOffset + i;
        gradx[pos] = sumx * invvol;
        grady[pos] = sumy * invvol;
        gradz[pos] = sumz * invvol;
      }
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

