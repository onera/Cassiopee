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

//=============================================================================
// Calcul du gradient d'un champ defini aux noeuds d une grille structuree
// retourne le gradient defini aux centres des cellules
//=============================================================================
E_Int K_POST::computeGradStruct(
  const E_Int ni, const E_Int nj, const E_Int nk, 
  const E_Float* xt, const E_Float* yt, const E_Float* zt, const E_Float* field,
  E_Float* gradx, E_Float* grady, E_Float* gradz
)
{
  if (ni*nj == 1 || ni*nk == 1 || nj*nk == 1)  // 1D
  {
    E_Int tni = ni;
    if (nj != 1) tni = nj;
    else if (nk != 1) tni = nk;
    compGradStruct1D(tni, xt, yt, zt, field, gradx, grady, gradz);
  }
  else if (ni == 1 || nj == 1 || nk == 1)  // 2D
  {
    compGradStruct2D(ni, nj, nk, xt, yt, zt, field, gradx, grady, gradz);
  }
  else  // 3D
  {
    compGradStruct3D(ni, nj, nk, xt, yt, zt, field, gradx, grady, gradz);
  }
  return 1;
}

// ============================================================================
// Cas 3D
// ============================================================================
void K_POST::compGradStruct3D(
  const E_Int ni, const E_Int nj, const E_Int nk,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* field,
  E_Float* gradx, E_Float* grady, E_Float* gradz
)
{
  E_Int nfld = 1;
  E_Int ni1 = K_FUNC::E_max(1, ni-1);
  E_Int nj1 = K_FUNC::E_max(1, nj-1);
  E_Int nk1 = K_FUNC::E_max(1, nk-1);
  E_Int ncells = ni1 * nj1 * nk1;

  E_Int inci = 1;
  E_Int incj = 1;
  E_Int inck = 1;

  if (ni == 2) { inci = 0; }
  else if (nj == 2) { incj = 0; }
  else if (nk == 2) { inck = 0; }
  
  E_Int ni1nj = ni1 * nj;
  E_Int ninj1 = ni * nj1;
  E_Int ni1nj1 = ni1 * nj1;

  E_Int inti = ninj1 * nk1;
  E_Int intj = ni1nj * nk1;
  E_Int intk = ni1nj1 * nk;
  E_Int intij = inti + intj;
  E_Int nint = inti + intj + intk;

  FldArrayF surfx(nint), surfy(nint), surfz(nint);
  FldArrayF snorm(nint);
  FldArrayF centerInt(nint, 3);
  FldArrayF fieldint(nint);
  FldArrayF vol(ncells);

  K_METRIC::compMetricStruct(
    ni, nj, nk, inti, intj, intk,
    xt, yt, zt,
    vol.begin(), surfx.begin(), surfy.begin(), surfz.begin(), snorm.begin(),
    centerInt.begin(1), centerInt.begin(2), centerInt.begin(3)
  );

  compIntField(ni, nj, nk, nfld, field, fieldint.begin());

  #pragma omp parallel
  {
    E_Int i, j, k;
    E_Int indint1, indint2, indint3, indint4, indint5, indint6;
    E_Float sx1, sx2, sx3, sx4, sx5, sx6;
    E_Float sy1, sy2, sy3, sy4, sy5, sy6;
    E_Float sz1, sz2, sz3, sz4, sz5, sz6;
    E_Float vinv;

    #pragma omp for
    for (E_Int indcell = 0; indcell < ncells; indcell++)
    {
      k = indcell / ni1nj1;
      j = (indcell - k*ni1nj1) / ni1;
      i = indcell - j*ni1 - k*ni1nj1;

      indint1 = i + j*ni + k*ninj1;
      indint2 = indint1 + inci;
      indint3 = i + j*ni1 + k*ni1nj + inti;
      indint4 = indint3 + incj*ni1;
      indint5 = i + j*ni1 + k*ni1nj1 + intij;
      indint6 = indint5 + inck*ni1nj1;

      sx1 = -surfx[indint1];
      sx2 =  surfx[indint2];
      sx3 = -surfx[indint3];
      sx4 =  surfx[indint4];
      sx5 = -surfx[indint5];
      sx6 =  surfx[indint6];

      sy1 = -surfy[indint1];
      sy2 =  surfy[indint2];
      sy3 = -surfy[indint3];
      sy4 =  surfy[indint4];
      sy5 = -surfy[indint5];
      sy6 =  surfy[indint6];

      sz1 = -surfz[indint1];
      sz2 =  surfz[indint2];
      sz3 = -surfz[indint3];
      sz4 =  surfz[indint4];
      sz5 = -surfz[indint5];
      sz6 =  surfz[indint6];

      vinv = K_CONST::ONE/K_FUNC::E_max(vol[indcell], K_CONST::E_MIN_VOL);

      gradx[indcell] = vinv * (sx1 * fieldint[indint1] + sx2 * fieldint[indint2]
                     + sx3 * fieldint[indint3] + sx4 * fieldint[indint4]
                     + sx5 * fieldint[indint5] + sx6 * fieldint[indint6]);

      grady[indcell] = vinv * (sy1 * fieldint[indint1] + sy2 * fieldint[indint2]
                     + sy3 * fieldint[indint3] + sy4 * fieldint[indint4]
                     + sy5 * fieldint[indint5] + sy6 * fieldint[indint6]);

      gradz[indcell] = vinv * (sz1 * fieldint[indint1] + sz2 * fieldint[indint2]
                     + sz3 * fieldint[indint3] + sz4 * fieldint[indint4]
                     + sz5 * fieldint[indint5] + sz6 * fieldint[indint6]);
    }
  }
}

// ============================================================================
// Cas 2D
// ============================================================================
void K_POST::compGradStruct2D(
  const E_Int ini, const E_Int inj, const E_Int ink,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* field,
  E_Float* gradx, E_Float* grady, E_Float* gradz
)
{
  E_Int ni = ini;
  E_Int nj = inj;
  E_Int nk = ink;
  if (ni == 1) { ni = nj; nj = nk; nk = 2; }
  else if (nj == 1) { nj = nk; nk = 2; }
  else if (nk == 1) { nk = 2; }

  E_Int ni1 = K_FUNC::E_max(1, ni-1);
  E_Int nj1 = K_FUNC::E_max(1, nj-1);
  E_Int nk1 = K_FUNC::E_max(1, nk-1);
  E_Int ncells = ni1 * nj1 * nk1;

  // Calcul de la surface totale des cellules
  FldArrayF surf(ncells), nxt(ncells), nyt(ncells), nzt(ncells);
  K_METRIC::compSurfStruct2D(ni, nj, 1, xt, yt, zt, surf.begin());
  K_METRIC::compNormStructSurf(
    ni, nj,
    xt, yt, zt,
    nxt.begin(), nyt.begin(), nzt.begin()
  );

  #pragma omp parallel
  {
    E_Int indcell, indA, indB, indC, indD;
    E_Float xAB, yAB, zAB, xBC, yBC, zBC, xCD, yCD, zCD, xDA, yDA, zDA;
    E_Float nx, n1x, n2x, n3x, n4x;
    E_Float ny, n1y, n2y, n3y, n4y;
    E_Float nz, n1z, n2z, n3z, n4z;
    E_Float fAB, fBC, fCD, fDA;
    E_Float vinv, nn;

    #pragma omp for collapse(2)
    for(E_Int j = 0; j < nj1; j++)
    for(E_Int i = 0; i < ni1; i++)
    {
      indA = i + j * ni;
      indB = indA + 1;
      indC = indB + ni;
      indD = indA + ni;

      indcell = i + j * ni1;

      nx = nxt[indcell];
      ny = nyt[indcell];
      nz = nzt[indcell];

      nn = sqrt(nx * nx + ny * ny + nz * nz);
      vinv = 2.0 * surf[indcell] * nn;
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

      gradx[indcell] = vinv * (fAB * n1x + fBC * n2x + fCD * n3x + fDA * n4x);
      grady[indcell] = vinv * (fAB * n1y + fBC * n2y + fCD * n3y + fDA * n4y);
      gradz[indcell] = vinv * (fAB * n1z + fBC * n2z + fCD * n3z + fDA * n4z);
    }
  }
}

// ============================================================================
// Cas 1D
// ============================================================================
void K_POST::compGradStruct1D(
  const E_Int ni, const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* field,
  E_Float* gradx, E_Float* grady, E_Float* gradz
)
{
  E_Int ni1 = K_FUNC::E_max(1, ni-1);
  E_Float eps = 1.e-12;

  E_Int indA, indB;
  E_Float xAB, yAB, zAB;
  E_Float lengthAB2;
  E_Float fAB;

  for (E_Int i = 0; i < ni1; i++)
  {
    indA = i;
    indB = indA + 1;

    xAB = xt[indB] - xt[indA];
    yAB = yt[indB] - yt[indA];
    zAB = zt[indB] - zt[indA];
    lengthAB2 = xAB * xAB + yAB * yAB + zAB * zAB;

    fAB = field[indB] - field[indA];
    gradx[i] = fAB * xAB / K_FUNC::E_max(lengthAB2, eps);
    grady[i] = fAB * yAB / K_FUNC::E_max(lengthAB2, eps);
    gradz[i] = fAB * zAB / K_FUNC::E_max(lengthAB2, eps);
  }
}
