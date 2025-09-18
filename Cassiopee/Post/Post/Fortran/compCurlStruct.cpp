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
// Calcul du rotationnel d'un vecteur defini aux noeuds d une grille structuree
// Retourne le rotationnel defini aux centres des cellules
// IN: ni,nj,nk: dimensions du maillage en noeuds
// IN: ncells: nbre de cellules
// IN: xt,yt,zt: coordonnees de la grille
// IN: ux, uy, uz: vecteur dont le rotationnel est a calculer
// OUT: rotx, roty, rotz: rotationnel aux centres des cellules
//==============================================================================
E_Int K_POST::computeCurlStruct(
  const E_Int ni, const E_Int nj, const E_Int nk, 
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* ux, const E_Float* uy, const E_Float* uz,
  E_Float* rotx, E_Float* roty, E_Float* rotz
)
{
  if (ni*nj == 1 || ni*nk == 1 || nj*nk == 1) return -1;
  if (ni == 1 || nj == 1 || nk == 1)
  {
    compCurlStruct2D(ni, nj, nk, xt, yt, zt, ux, uy, uz, rotx, roty, rotz);
  }
  else
  {
    compCurlStruct3D(ni, nj, nk, xt, yt, zt, ux, uy, uz, rotx, roty, rotz);
  }
  return 1;
}

// ============================================================================
//  Cas 3D
// ============================================================================
void K_POST::compCurlStruct3D(
  const E_Int ni, const E_Int nj, const E_Int nk,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* ux, const E_Float* uy, const E_Float* uz,
  E_Float* rotx, E_Float* roty, E_Float* rotz
)
{
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
  FldArrayF uintx(nint), uinty(nint), uintz(nint);
  FldArrayF vol(ncells);

  // Attention : surf n est pas oriente : tjs positif
  K_METRIC::compMetricStruct(
    ni, nj, nk, inti, intj, intk,
    xt, yt, zt,
    vol.begin(), surfx.begin(), surfy.begin(), surfz.begin(), snorm.begin(),
    centerInt.begin(1), centerInt.begin(2), centerInt.begin(3)
  );

  compIntFieldV(
    ni, nj, nk,
    ux, uy, uz,
    uintx.begin(), uinty.begin(), uintz.begin()
  );

  #pragma omp parallel
  {
    E_Int i, j, k;
    E_Int indint1, indint2, indint3, indint4, indint5, indint6;
    E_Float curlx, curly, curlz;
    E_Float sx1, sx2, sx3, sx4, sx5, sx6;
    E_Float sy1, sy2, sy3, sy4, sy5, sy6;
    E_Float sz1, sz2, sz3, sz4, sz5, sz6;
    E_Float vinv;

    #pragma omp for
    for (E_Int indcell = 0; indcell < ncells; indcell++)
    {
      k = indcell / ni1nj1;
      j = (indcell - k * ni1nj1) / ni1;
      i = indcell - j * ni1 - k * ni1nj1;

      indint1 = i + j * ni + k * ninj1;
      indint2 = indint1 + inci;
      indint3 = i + j * ni1 + k * ni1nj + inti;
      indint4 = indint3 + incj * ni1;
      indint5 = i + j * ni1 + k * ni1nj1 + intij;
      indint6 = indint5 + inck * ni1nj1;

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

      vinv = -K_CONST::ONE / K_FUNC::E_max(vol[indcell], K_CONST::E_MIN_VOL);

      curlx = uinty[indint1] * sz1 - uintz[indint1] * sy1
            + uinty[indint2] * sz2 - uintz[indint2] * sy2
            + uinty[indint3] * sz3 - uintz[indint3] * sy3
            + uinty[indint4] * sz4 - uintz[indint4] * sy4
            + uinty[indint5] * sz5 - uintz[indint5] * sy5
            + uinty[indint6] * sz6 - uintz[indint6] * sy6;

      curly = uintz[indint1] * sx1 - uintx[indint1] * sz1
            + uintz[indint2] * sx2 - uintx[indint2] * sz2
            + uintz[indint3] * sx3 - uintx[indint3] * sz3
            + uintz[indint4] * sx4 - uintx[indint4] * sz4
            + uintz[indint5] * sx5 - uintx[indint5] * sz5
            + uintz[indint6] * sx6 - uintx[indint6] * sz6;

      curlz = uintx[indint1] * sy1 - uinty[indint1] * sx1
            + uintx[indint2] * sy2 - uinty[indint2] * sx2
            + uintx[indint3] * sy3 - uinty[indint3] * sx3
            + uintx[indint4] * sy4 - uinty[indint4] * sx4
            + uintx[indint5] * sy5 - uinty[indint5] * sx5
            + uintx[indint6] * sy6 - uinty[indint6] * sx6;

      rotx[indcell] = vinv * curlx;
      roty[indcell] = vinv * curly;
      rotz[indcell] = vinv * curlz;
    }
  }
}

// ============================================================================
// Cas 2D
// ============================================================================
void K_POST::compCurlStruct2D(
  const E_Int ini, const E_Int inj, const E_Int ink,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  const E_Float* ux, const E_Float* uy, const E_Float* uz,
  E_Float* rotux, E_Float* rotuy, E_Float* rotuz
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
    E_Int indA, indB, indC, indD, indcell;
    E_Float nx, ny, nz, nn;
    E_Float vinv, curlx, curly, curlz;
    E_Float xAB, yAB, zAB, xBC, yBC, zBC, xCD, yCD, zCD, xDA, yDA, zDA;
    E_Float n1x, n1y, n1z, vx, vy, vz;

    #pragma omp for collapse(2)
    for (E_Int j = 0; j < nj1; j++)
    for (E_Int i = 0; i < ni1; i++)
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

      rotux[indcell] = vinv * curlx;
      rotuy[indcell] = vinv * curly;
      rotuz[indcell] = vinv * curlz;
    }
  }
}

// ============================================================================
// Calcul du rotationnel moyen d un champ (u,v,w) sur une cellule
// attention cette routine est uniquement 3d
// IN: ind: indice du premier sommet de la cellule
// IN: ni,nj,nk: dimensions du maillage en noeuds
// IN: velo: vecteur dont le rotationnel est a calculer. Defini sur le maillage
// OUT: rotu,rotv,rotw: rotationnel moyen de (u,v,w) sur la cellule
// ============================================================================
void K_POST::compMeanCurlOfStructCell(
  const E_Int ind, const E_Int ni, const E_Int nj, const E_Int nk,
  const E_Float* velox, const E_Float* veloy, const E_Float* veloz,
  const E_Float* xt, const E_Float* yt, const E_Float* zt,
  E_Float& rotu, E_Float& rotv, E_Float& rotw
)
{
  E_Float cellvol;
  FldArrayF uint(6, 3), surf(6, 3);
  E_Float *uintx = uint.begin(1), *uinty = uint.begin(2), *uintz = uint.begin(3);
  E_Float *surfx = surf.begin(1), *surfy = surf.begin(2), *surfz = surf.begin(3);

  // calcul des surfaces aux interfaces
  K_METRIC::compIntSurfOfCell(ind, ni, nj, nk, xt, yt, zt, surf.begin());

  // calcul du volume de la cellule
  K_METRIC::compVolOfStructCell3D(ni, nj, nk, -1, ind, xt, yt, zt, cellvol);

  E_Int l1, l2, l3, l4;
  E_Int ninj = ni*nj;

  // interface en i
  l1 = ind;
  l2 = l1 + ni;
  l3 = l2 + ninj;
  l4 = l1 + ninj;
  uintx[0] = velox[l1] + velox[l2] + velox[l3] + velox[l4];
  uinty[0] = veloy[l1] + veloy[l2] + veloy[l3] + veloy[l4];
  uintz[0] = veloz[l1] + veloz[l2] + veloz[l3] + veloz[l4];

  // interface en i+1
  l1 = ind + 1;
  l2 = l1 + ni;
  l3 = l2 + ninj;
  l4 = l1 + ninj;
  uintx[1] = velox[l1] + velox[l2] + velox[l3] + velox[l4];
  uinty[1] = veloy[l1] + veloy[l2] + veloy[l3] + veloy[l4];
  uintz[1] = veloz[l1] + veloz[l2] + veloz[l3] + veloz[l4];

  // interface en j
  l1 = ind + 1;
  l2 = l1 - 1;
  l3 = l2 + ninj;
  l4 = l1 + ninj;
  uintx[2] = velox[l1] + velox[l2] + velox[l3] + velox[l4];
  uinty[2] = veloy[l1] + veloy[l2] + veloy[l3] + veloy[l4];
  uintz[2] = veloz[l1] + veloz[l2] + veloz[l3] + veloz[l4];

  // interface en j+1
  l1 = ind + ni + 1;
  l2 = l1 - 1;
  l3 = l2 + ninj;
  l4 = l1 + ninj;
  uintx[3] = velox[l1] + velox[l2] + velox[l3] + velox[l4];
  uinty[3] = veloy[l1] + veloy[l2] + veloy[l3] + veloy[l4];
  uintz[3] = veloz[l1] + veloz[l2] + veloz[l3] + veloz[l4];

  // interface en k
  l1 = ind;
  l2 = l1 + 1;
  l3 = l2 + ni;
  l4 = l1 + ni;
  uintx[4] = velox[l1] + velox[l2] + velox[l3] + velox[l4];
  uinty[4] = veloy[l1] + veloy[l2] + veloy[l3] + veloy[l4];
  uintz[4] = veloz[l1] + veloz[l2] + veloz[l3] + veloz[l4];

  // interface en k+1
  l1 = ind + ninj;
  l2 = l1 + 1;
  l3 = l2 + ni;
  l4 = l1 + ni;
  uintx[5] = velox[l1] + velox[l2] + velox[l3] + velox[l4];
  uinty[5] = veloy[l1] + veloy[l2] + veloy[l3] + veloy[l4];
  uintz[5] = veloz[l1] + veloz[l2] + veloz[l3] + veloz[l4];

  E_Float curlx = 0., curly = 0., curlz = 0.;
  for (E_Int l = 0; l < 6; l++)
  {
    curlx += uinty[l] * surfz[l] - uintz[l] * surfy[l];
    curly += uintz[l] * surfx[l] - uintx[l] * surfz[l];
    curlz += uintx[l] * surfy[l] - uinty[l] * surfx[l];
  }

  cellvol = -K_CONST::ONE_FOURTH / K_FUNC::E_max(cellvol, K_CONST::E_MIN_VOL);

  rotu = cellvol * curlx;
  rotv = cellvol * curly;
  rotw = cellvol * curlz;
}
