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
# include "post.h"
# include "Array/Array.h"
# include <math.h>
# include <vector>

/* Calcul du gradient d'un champ defini aux noeuds d une grille non structuree
   retourne le gradient defini aux centres des elts.
   Cas 1D, 2D et 3D.
*/
E_Int K_POST::computeGradUnstruct(
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  K_FLD::FldArrayI& cn, const char* eltType,
  const E_Float* field,
  E_Float* gradx, E_Float* grady, E_Float* gradz
)
{
  E_Int dim = K_CONNECT::getDimME(eltType);
  if (dim == 3)
  {
    compGradUnstruct3D(xt, yt, zt, cn, eltType, field, gradx, grady, gradz);
  }
  else if (dim == 2)
  {
    compGradUnstruct2D(xt, yt, zt, cn, eltType, field, gradx, grady, gradz);
  }
  else  // dim = 1
  {
    compGradUnstruct1D(xt, yt, zt, cn, eltType, field, gradx, grady, gradz);
  }
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
void K_POST::compGradUnstruct1D(
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
      fprintf(stderr, "Error: in K_POST::compGradUnstruct1D.\n");
      fprintf(stderr, "Element type must be BAR, not %s.\n", eltTypes[ic]);
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
   OUT: gradx, grady, gradz: gradient de field %x, %y, %z
*/
void K_POST::compGradUnstruct2D(
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  K_FLD::FldArrayI& cn, const char* eltType, const E_Float* field,
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
    if (strcmp(eltTypes[ic], "TRI") != 0 && strcmp(eltTypes[ic], "QUAD") != 0)
    {
      fprintf(stderr, "Error: in K_POST::compGradUnstruct2D.\n");
      fprintf(stderr, "Element type must be TRI or QUAD, not %s.\n", eltTypes[ic]);
      exit(0);
    }
    nepc[ic+1] = nepc[ic] + nelts;
  }

  // Allocate memory to store facet normals and their areas for all
  // connectivities
  E_Int ntotElts = nepc[nc];
  K_FLD::FldArrayF snx(ntotElts), sny(ntotElts), snz(ntotElts), surf(ntotElts);
  
  // Compute surface of elements
  K_METRIC::compSurfUnstruct(
    cn, eltType, xt, yt, zt,
    snx.begin(), sny.begin(), snz.begin(), surf.begin()
  );

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
   OUT: gradx, grady, gradz: gradient de field %x, %y, %z
*/
void K_POST::compGradUnstruct3D(
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  K_FLD::FldArrayI& cn, const char* eltType, const E_Float* field,
  E_Float* gradx, E_Float* grady, E_Float* gradz
)
{
  // Pre-compute element and facet offsets
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  E_Int ntotFacets = 0;
  E_Int ntotElts = 0;
  std::vector<E_Int> nepc(nc+1), nfpc(nc+1);
  nepc[0] = 0; nfpc[0] = 0;

  // number of facets per element
  std::vector<E_Int> nfpe(nc);  
  E_Int ierr = K_CONNECT::getNFPE(nfpe, eltType, false);
  if (ierr != 0) return;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    nepc[ic+1] = nepc[ic] + nelts;
    nfpc[ic+1] = nfpc[ic] + nfpe[ic]*nelts;  // number of facets per connectivity
    ntotFacets += nfpe[ic]*nelts;
    ntotElts += nelts;
  }

  // Allocate memory to store facet normals and their areas for all
  // connectivities, as well as the volume of the elements and fieldf, the ith
  // field defined for each facet of the elements
  FldArrayF fieldf(ntotFacets);
  FldArrayF snx(ntotFacets), sny(ntotFacets), snz(ntotFacets);
  FldArrayF surf(ntotFacets);
  FldArrayF vol(ntotElts);

  // Compute facet areas and element volumes
  K_METRIC::compMetricUnstruct(
    cn, eltType, xt, yt, zt,
    snx.begin(), sny.begin(), snz.begin(), surf.begin(), vol.begin()
  );

  // Compute field on face centers
  compUnstrNodes2Faces(cn, eltType, field, fieldf.begin());

  // Compute gradient at element centers  
  #pragma omp parallel
  {
    E_Int pose, posf, tposf, nelts, elOffset, fctOffset;
    E_Float sumx, sumy, sumz, invvol;
    
    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
      nelts = cm.getSize();
      elOffset = nepc[ic];
      fctOffset = nfpc[ic];

      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        sumx = K_CONST::E_ZERO_FLOAT;
        sumy = K_CONST::E_ZERO_FLOAT;
        sumz = K_CONST::E_ZERO_FLOAT;
        pose = elOffset + i;
        tposf = fctOffset + i * nfpe[ic];
        invvol = K_CONST::ONE / K_FUNC::E_max(vol[pose], K_CONST::E_MIN_VOL);

        for (E_Int fi = 0; fi < nfpe[ic]; fi++)
        {
          posf = tposf + fi;
          sumx += snx[posf] * fieldf[posf];
          sumy += sny[posf] * fieldf[posf];
          sumz += snz[posf] * fieldf[posf];
        }
        
        gradx[pose] = sumx * invvol;
        grady[pose] = sumy * invvol;
        gradz[pose] = sumz * invvol;
      }
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

