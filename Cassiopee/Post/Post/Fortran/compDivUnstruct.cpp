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

// ============================================================================
// Calcul de la divergence d'un champ defini aux noeuds d'une grille
// non-structuree
// Retourne la divergence definie aux centres des elements
// ============================================================================
E_Int K_POST::computeDivUnstruct(
  FldArrayI& cn, const char* eltType,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* fieldX, const E_Float* fieldY, const E_Float* fieldZ,
  E_Float* div
)
{
  E_Int dim = K_CONNECT::getDimME(eltType);
  if (dim == 2)
  {
    compDivUnstruct2D(cn, eltType, xt, yt, zt, fieldX, fieldY, fieldZ, div);
  }
  else if (dim == 3)
  {
    compDivUnstruct3D(cn, eltType, xt, yt, zt, fieldX, fieldY, fieldZ, div);
  }
  else return -1;
  return 1;
}

// ============================================================================
// CAS 3D
// ============================================================================
void K_POST::compDivUnstruct3D(
  FldArrayI& cn, const char* eltType,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* fieldX, const E_Float* fieldY, const E_Float* fieldZ,
  E_Float* div
)
{
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Number of facets per element
  std::vector<E_Int> nfpe;
  E_Int ierr = K_CONNECT::getNFPE(nfpe, eltType, false);
  if (ierr != 0) return;
  
  // Pre-compute element and facet offsets
  E_Int ntotFacets = 0;
  E_Int ntotElts = 0;
  std::vector<E_Int> nepc(nc+1), nfpc(nc+1);
  nepc[0] = 0; nfpc[0] = 0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    nepc[ic+1] = nepc[ic] + nelts;
    nfpc[ic+1] = nfpc[ic] + nfpe[ic]*nelts;  // number of facets per connectivity
    ntotFacets += nfpe[ic]*nelts;
    ntotElts += nelts;
  }
  
  // Calcul de la surface + volume des elts
  FldArrayF snx(ntotFacets), sny(ntotFacets), snz(ntotFacets), surf(ntotFacets);
  FldArrayF fieldfx(ntotFacets), fieldfy(ntotFacets), fieldfz(ntotFacets);
  FldArrayF vol(ntotElts);
  K_METRIC::compMetricUnstruct(
    cn, eltType,
    xt, yt, zt,
    snx.begin(), sny.begin(), snz.begin(), surf.begin(), vol.begin()
  );

  // Calcul des champs aux centres des facettes
  compUnstrNodes2Faces(cn, eltType, fieldX, fieldfx.begin());
  compUnstrNodes2Faces(cn, eltType, fieldY, fieldfy.begin());
  compUnstrNodes2Faces(cn, eltType, fieldZ, fieldfz.begin());

  // Calcul de la divergence au centre des elts
  #pragma omp parallel
  {
    E_Int pose, posf, tposf, nelts, elOffset, fctOffset;
    E_Float invvol, gradxx, gradyy, gradzz;

    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
      nelts = cm.getSize();
      elOffset = nepc[ic];
      fctOffset = nfpc[ic];
      
      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        gradxx = K_CONST::E_ZERO_FLOAT;
        gradyy = K_CONST::E_ZERO_FLOAT;
        gradzz = K_CONST::E_ZERO_FLOAT;
        pose = elOffset + i;
        tposf = fctOffset + i * nfpe[ic];

        for (E_Int fi = 0; fi < nfpe[ic]; fi++)
        {
          posf = tposf + fi;
          gradxx += snx[posf] * fieldfx[posf];
          gradyy += sny[posf] * fieldfy[posf];
          gradzz += snz[posf] * fieldfz[posf];
        }

        invvol = K_CONST::ONE / K_FUNC::E_max(vol[pose], K_CONST::E_MIN_VOL);
        div[pose] = (gradxx + gradyy + gradzz) * invvol;
      }
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

// ============================================================================
// CAS 2D
// ============================================================================
void K_POST::compDivUnstruct2D(
  FldArrayI& cn, const char* eltType,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* fieldX, const E_Float* fieldY, const E_Float* fieldZ,
  E_Float* div
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
      fprintf(stderr, "Error: in K_POST::compDivUnstruct2D.\n");
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
    E_Float fxAB, fxBC, fxCD, fxDA, fxCA;
    E_Float fyAB, fyBC, fyCD, fyDA, fyCA;
    E_Float fzAB, fzBC, fzCD, fzDA, fzCA;
    E_Float gradxx, gradyy, gradzz;

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
          vinv = 2.0 * surf[pos] * nn;
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

          fxAB = fieldX[indA] + fieldX[indB];
          fxBC = fieldX[indB] + fieldX[indC];
          fxCA = fieldX[indC] + fieldX[indA];

          fyAB = fieldY[indA] + fieldY[indB];
          fyBC = fieldY[indB] + fieldY[indC];
          fyCA = fieldY[indC] + fieldY[indA];

          fzAB = fieldZ[indA] + fieldZ[indB];
          fzBC = fieldZ[indB] + fieldZ[indC];
          fzCA = fieldZ[indC] + fieldZ[indA];

          gradxx = fxAB * n1x + fxBC * n2x + fxCA * n3x;
          gradyy = fyAB * n1y + fyBC * n2y + fyCA * n3y;
          gradzz = fzAB * n1z + fzBC * n2z + fzCA * n3z;

          div[pos] = vinv * (gradxx + gradyy + gradzz);
        }
      }
      else // QUAD
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
          vinv = 2.0 * surf[pos] * nn;
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

          fxAB = fieldX[indA] + fieldX[indB];
          fxBC = fieldX[indB] + fieldX[indC];
          fxCD = fieldX[indC] + fieldX[indD];
          fxDA = fieldX[indD] + fieldX[indA];

          fyAB = fieldY[indA] + fieldY[indB];
          fyBC = fieldY[indB] + fieldY[indC];
          fyCD = fieldY[indC] + fieldY[indD];
          fyDA = fieldY[indD] + fieldY[indA];

          fzAB = fieldZ[indA] + fieldZ[indB];
          fzBC = fieldZ[indB] + fieldZ[indC];
          fzCD = fieldZ[indC] + fieldZ[indD];
          fzDA = fieldZ[indD] + fieldZ[indA];

          gradxx = fxAB * n1x + fxBC * n2x + fxCD * n3x + fxDA * n4x;
          gradyy = fyAB * n1y + fyBC * n2y + fyCD * n3y + fyDA * n4y;
          gradzz = fzAB * n1z + fzBC * n2z + fzCD * n3z + fzDA * n4z;

          div[pos] = vinv * (gradxx + gradyy + gradzz);
        }
      }
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}
