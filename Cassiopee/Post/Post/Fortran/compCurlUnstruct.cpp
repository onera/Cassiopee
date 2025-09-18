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

//==============================================================================
// Calcul du rotationnel d'un vecteur defini aux noeuds d une grille non structuree
// Retourne le rotationnel defini aux centres des cellules
//==============================================================================
E_Int K_POST::computeCurlUnstruct(
  FldArrayI& cn, const char* eltType,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* ux, const E_Float* uy, const E_Float* uz,
  E_Float* rotx, E_Float* roty, E_Float* rotz
)
{
  // Get ME mesh dimensionality from the first element type
  E_Int dim = 3;
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);
  if (strcmp(eltTypes[0], "BAR") == 0) dim = 1;
  else if (strcmp(eltTypes[0], "TRI") == 0 or
           strcmp(eltTypes[0], "QUAD") == 0) dim = 2;
  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];

  if (dim == 2)
  {
    compCurlUnstruct2D(cn, eltType, xt, yt, zt, ux, uy, uz, rotx, roty, rotz);
  }
  else if (dim == 3)
  {
    compCurlUnstruct3D(cn, eltType, xt, yt, zt, ux, uy, uz, rotx, roty, rotz);
  }
  else return -1;
  return 1;
}

// ============================================================================
//  Cas 3D
// ============================================================================
void K_POST::compCurlUnstruct3D(
  FldArrayI& cn, const char* eltType,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* ux, const E_Float* uy, const E_Float* uz,
  E_Float* rotx, E_Float* roty, E_Float* rotz
)
{
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);
  
  // Pre-compute element and facet offsets
  E_Int ntotFacets = 0;
  E_Int ntotElts = 0;
  std::vector<E_Int> nfpe(nc);
  std::vector<E_Int> nepc(nc+1), nfpc(nc+1);
  nepc[0] = 0; nfpc[0] = 0;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    if (strcmp(eltTypes[ic], "TRI") == 0) nfpe[ic] = 1;
    else if (strcmp(eltTypes[ic], "QUAD") == 0) nfpe[ic] = 1;
    else if (strcmp(eltTypes[ic], "TETRA") == 0) nfpe[ic] = 4;
    else if (strcmp(eltTypes[ic], "PYRA") == 0) nfpe[ic] = 5;
    else if (strcmp(eltTypes[ic], "PENTA") == 0) nfpe[ic] = 5;
    else if (strcmp(eltTypes[ic], "HEXA") == 0) nfpe[ic] = 6;
    else
    {
      fprintf(stderr, "Error: in K_POST::compCurlUnstruct3D.\n");
      fprintf(stderr, "Unknown type of element, %s.\n", eltTypes[ic]);
    }
    nepc[ic+1] = nepc[ic] + nelts;
    nfpc[ic+1] = nfpc[ic] + nfpe[ic]*nelts;  // number of facets per connectivity
    ntotFacets += nfpe[ic]*nelts;
    ntotElts += nelts;
  }
  
  // Allocation of temp fields
  FldArrayF snx(ntotFacets), sny(ntotFacets), snz(ntotFacets), surf(ntotFacets);
  FldArrayF uintx(ntotFacets), uinty(ntotFacets), uintz(ntotFacets);
  FldArrayF vol(ntotElts);
    
  // Calcul de la surface + volume des elts
  K_METRIC::compMetricUnstruct(
    cn, eltType,
    xt, yt, zt,
    snx.begin(), sny.begin(), snz.begin(), surf.begin(), vol.begin()
  );

  // Calcul des composantes du vecteur aux centres des facettes 
  compUnstrNodes2Faces(cn, eltType, ux, uintx.begin());
  compUnstrNodes2Faces(cn, eltType, uy, uinty.begin());
  compUnstrNodes2Faces(cn, eltType, uz, uintz.begin());

  // Calcul du rotationnel au centre des elts
  #pragma omp parallel
  {
    E_Int pose, posf, tposf, nelts, elOffset, fctOffset;
    E_Float curlx, curly, curlz, vinv;
    E_Float sxi, syi, szi;

    for (E_Int ic = 0; ic < nc; ic++)
    {
      K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
      nelts = cm.getSize();
      elOffset = nepc[ic];
      fctOffset = nfpc[ic];
      
      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        pose = elOffset + i;
        
        curlx = K_CONST::E_ZERO_FLOAT;
        curly = K_CONST::E_ZERO_FLOAT;
        curlz = K_CONST::E_ZERO_FLOAT;
        pose = elOffset + i;
        tposf = fctOffset + i * nfpe[ic];

        vinv = -K_CONST::ONE / K_FUNC::E_max(vol[pose], K_CONST::E_MIN_VOL);

        for (E_Int fi = 0; fi < nfpe[ic]; fi++)
        {
          posf = tposf + fi;
          sxi = snx[posf];
          syi = sny[posf];
          szi = snz[posf];

          curlx += uinty[posf] * szi - uintz[posf] * syi;
          curly += uintz[posf] * sxi - uintx[posf] * szi;
          curlz += uintx[posf] * syi - uinty[posf] * sxi;
        }
        
        rotx[pose] = vinv * curlx;
        roty[pose] = vinv * curly;
        rotz[pose] = vinv * curlz;
      }
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
}

// ============================================================================
//  Cas 2D
// ============================================================================
void K_POST::compCurlUnstruct2D(
  FldArrayI& cn, const char* eltType,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* ux, const E_Float* uy, const E_Float* uz,
  E_Float* rotx, E_Float* roty, E_Float* rotz
)
{
  E_Int nc = cn.getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Pre-compute element and facet offsets
  E_Int ntotElts = 0;
  std::vector<E_Int> nepc(nc+1); nepc[0] = 0;
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    nepc[ic+1] = nepc[ic] + nelts;
    ntotElts += nelts;
  }
  
  // Calcul de la surface des elts
  FldArrayF snx(ntotElts), sny(ntotElts), snz(ntotElts), surf(ntotElts);
  K_METRIC::compSurfUnstruct(
    cn, eltType,
    xt, yt, zt,
    snx.begin(), sny.begin(), snz.begin(), surf.begin()
  );

  #pragma omp parallel
  {
    E_Int indA, indB, indC, indD, pos;
    E_Float curlx, curly, curlz, vinv;
    E_Float vx, vy, vz;
    E_Float nx, ny, nz, nn, n1x, n1y, n1z;
    E_Float xAB, yAB, zAB, xBC, yBC, zBC;
    E_Float xCD, yCD, zCD, xDA, yDA, zDA;
    E_Float xCA, yCA, zCA;

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
          pos = elOffset + i;
          
          nx = snx[pos];
          ny = sny[pos];
          nz = snz[pos];
          nn = sqrt(nx * nx + ny * ny + nz * nz);

          vinv = 2.0 * surf[pos] * nn;
          vinv = -K_CONST::ONE / K_FUNC::E_max(vinv, K_CONST::E_MIN_VOL);

          indA = cm(i, 1) - 1;
          indB = cm(i, 2) - 1;
          indC = cm(i, 3) - 1;

          xAB = xt[indB] - xt[indA];
          yAB = yt[indB] - yt[indA];
          zAB = zt[indB] - zt[indA];

          xBC = xt[indC] - xt[indB];
          yBC = yt[indC] - yt[indB];
          zBC = zt[indC] - zt[indB];

          xCA = xt[indA] - xt[indC];
          yCA = yt[indA] - yt[indC];
          zCA = zt[indA] - zt[indC];

          curlx = K_CONST::E_ZERO_FLOAT;
          curly = K_CONST::E_ZERO_FLOAT;
          curlz = K_CONST::E_ZERO_FLOAT;

          n1x = yAB * nz - zAB * ny;
          n1y = zAB * nx - xAB * nz;
          n1z = xAB * ny - yAB * nx;

          vx = ux[indA] + ux[indB];
          vy = uy[indA] + uy[indB];
          vz = uz[indA] + uz[indB];

          curlx += vy * n1z - vz * n1y;
          curly += vz * n1x - vx * n1z;
          curlz += vx * n1y - vy * n1x;

          n1x = yBC * nz - zBC * ny;
          n1y = zBC * nx - xBC * nz;
          n1z = xBC * ny - yBC * nx;

          vx = ux[indC] + ux[indB];
          vy = uy[indC] + uy[indB];
          vz = uz[indC] + uz[indB];

          curlx += vy * n1z - vz * n1y;
          curly += vz * n1x - vx * n1z;
          curlz += vx * n1y - vy * n1x;

          n1x = yCA * nz - zCA * ny;
          n1y = zCA * nx - xCA * nz;
          n1z = xCA * ny - yCA * nx;

          vx = ux[indC] + ux[indA];
          vy = uy[indC] + uy[indA];
          vz = uz[indC] + uz[indA];

          curlx += vy * n1z - vz * n1y;
          curly += vz * n1x - vx * n1z;
          curlz += vx * n1y - vy * n1x;

          rotx[pos] = vinv * curlx;
          roty[pos] = vinv * curly;
          rotz[pos] = vinv * curlz;
        }
      }
      else // QUAD
      {
        #pragma omp for
        for (E_Int i = 0; i < nelts; i++)
        {
          pos = elOffset + i;
          
          indA = cm(i, 1) - 1;
          indB = cm(i, 2) - 1;
          indC = cm(i, 3) - 1;
          indD = cm(i, 4) - 1;

          nx = snx[pos];
          ny = sny[pos];
          nz = snz[pos];
          nn = sqrt(nx * nx + ny * ny + nz * nz);

          vinv = 2.0 * surf[pos] * nn;
          vinv = -K_CONST::ONE / K_FUNC::E_max(vinv, K_CONST::E_MIN_VOL);

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

          curlx = K_CONST::E_ZERO_FLOAT;
          curly = K_CONST::E_ZERO_FLOAT;
          curlz = K_CONST::E_ZERO_FLOAT;

          n1x = yAB * nz - zAB * ny;
          n1y = zAB * nx - xAB * nz;
          n1z = xAB * ny - yAB * nx;

          vx = ux[indA] + ux[indB];
          vy = uy[indA] + uy[indB];
          vz = uz[indA] + uz[indB];

          curlx += vy * n1z - vz * n1y;
          curly += vz * n1x - vx * n1z;
          curlz += vx * n1y - vy * n1x;

          n1x = yBC * nz - zBC * ny;
          n1y = zBC * nx - xBC * nz;
          n1z = xBC * ny - yBC * nx;

          vx = ux[indC] + ux[indB];
          vy = uy[indC] + uy[indB];
          vz = uz[indC] + uz[indB];

          curlx += vy * n1z - vz * n1y;
          curly += vz * n1x - vx * n1z;
          curlz += vx * n1y - vy * n1x;

          n1x = yCD * nz - zCD * ny;
          n1y = zCD * nx - xCD * nz;
          n1z = xCD * ny - yCD * nx;

          vx = ux[indC] + ux[indD];
          vy = uy[indC] + uy[indD];
          vz = uz[indC] + uz[indD];

          curlx += vy * n1z - vz * n1y;
          curly += vz * n1x - vx * n1z;
          curlz += vx * n1y - vy * n1x;

          n1x = yDA * nz - zDA * ny;
          n1y = zDA * nx - xDA * nz;
          n1z = xDA * ny - yDA * nx;

          vx = ux[indD] + ux[indA];
          vy = uy[indD] + uy[indA];
          vz = uz[indD] + uz[indA];

          curlx += vy * n1z - vz * n1y;
          curly += vz * n1x - vx * n1z;
          curlz += vx * n1y - vy * n1x;

          rotx[pos] = vinv * curlx;
          roty[pos] = vinv * curly;
          rotz[pos] = vinv * curlz;
        }
      }
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic]; 
}
