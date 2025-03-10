C
C    Copyright 2013-2025 Onera.
C
C    This file is part of Cassiopee.
C
C    Cassiopee is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    Cassiopee is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.

C Calcul de la divergence d'un champ defini aux noeuds d une grille structuree
C retourne la divergence defini aux centres des cellules
C ============================================================================
      SUBROUTINE k6compstructdiv(ni, nj, nk, nbcell, nbint,
     &                           xt, yt, zt,
     &                           fieldX, fieldY, fieldZ, div,
     &                           surf, snorm, centerInt, vol,
     &                           fldintX, fldintY, fldintZ)

      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E ni, nj, nk      ! dimensions de la grille aux noeuds
      INTEGER_E nbcell          ! nb de cellules
      INTEGER_E nbint           ! nb d interfaces
      REAL_E xt(0:ni*nj*nk-1)   ! coordonnees des noeuds de la grille
      REAL_E yt(0:ni*nj*nk-1)   ! coordonnees des noeuds de la grille
      REAL_E zt(0:ni*nj*nk-1)   ! coordonnees des noeuds de la grille
      REAL_E fieldX(0:ni*nj*nk-1) !champs (x-comp) defini aux noeuds auquel on applique div
      REAL_E fieldY(0:ni*nj*nk-1) !champs (y-comp) defini aux noeuds auquel on applique div
      REAL_E fieldZ(0:ni*nj*nk-1) !champs (z-comp) defini aux noeuds auquel on applique div
C_OUT
      REAL_E  div(0:nbcell-1)   ! divergence du champ vectoriel
      REAL_E  surf(0:nbint-1,3) ! surfaces
      REAL_E  snorm(0:nbint-1)  ! norme de surf : pas utilise
      REAL_E  centerInt(0:nbint-1,3)
      REAL_E  vol(0:nbcell-1)
      REAL_E  fldintX(0:nbint-1)
      REAL_E  fldintY(0:nbint-1)
      REAL_E  fldintZ(0:nbint-1)

C_LOCAL
      INTEGER_E inci, incj, inck
      INTEGER_E indcell, i, j, k
      INTEGER_E ninj, ni1, nj1, nk1, ni1nj, ninj1, ni1nj1
      INTEGER_E indint1, indint2, indint3, indint4, indint5, indint6
      INTEGER_E inti, intj, intij, intk
      REAL_E sx1, sx2, sx3, sx4, sx5, sx6
      REAL_E sy1, sy2, sy3, sy4, sy5, sy6
      REAL_E sz1, sz2, sz3, sz4, sz5, sz6
      REAL_E vinv              ! inverse du volume d une cellule
      REAL_E gradxx, gradyy, gradzz
      INTEGER_E nfld
C ------------------------------------------------------------------
C_CONSTANTES
      nfld = 1
      ni1 = ni-1
      nj1 = nj-1
      nk1 = nk-1
C
      inci = 1
      incj = 1
      inck = 1
      if (ni .eq. 2) then
         inci = 0
      else if (nj .eq. 2) then
         incj = 0
      else if (nk .eq. 2) then
         inck = 0
      endif
C
      ninj = ni*nj
      ni1nj = ni1*nj
      ninj1 = ni*nj1
      ni1nj1 = ni1*nj1
C
      inti = ninj1*nk1
      intj = ni1nj*nk1
      intk = ni1nj1*nk
      intij = inti + intj
C
C ------------------------------------------------------------------
C attention: surf n est pas oriente: tjs positif
      call k6compstructmetric(ni, nj, nk, nbcell, nbint,inti,intj,intk,
     &     xt, yt, zt, vol,
     &     surf, surf(0,2), surf(0,3), snorm,
     &     centerInt, centerInt(0,2), centerInt(0,3))

      call k6compintfield(ni, nj, nk, nfld, nbint, fieldX, fldintX)
      call k6compintfield(ni, nj, nk, nfld, nbint, fieldY, fldintY)
      call k6compintfield(ni, nj, nk, nfld, nbint, fieldZ, fldintZ)

!$OMP PARALLEL PRIVATE(indcell,i,j,k,indint1,indint2,indint3)
!$OMP& PRIVATE(indint4,indint5,indint6,sx1,sx2,sx3,sx4,sx5,sx6)
!$OMP& PRIVATE(sy1,sy2,sy3,sy4,sy5,sy6,sz1,sz2,sz3,sz4,sz5,sz6)
!$OMP& PRIVATE(vinv,gradxx,gradyy,gradzz)

!$OMP DO
      DO indcell = 0, nbcell-1
         k = indcell / ni1nj1     ! demarrent a 0
         j = (indcell-k*ni1nj1)/ni1
         i = indcell - j *ni1 - k * ni1nj1

         indint1 = i + j*ni + k*ninj1
         indint2 = indint1 + inci
         indint3 = i + j*ni1 + k*ni1nj + inti
         indint4 = indint3 + incj*ni1
         indint5 = i + j*ni1 + k*ni1nj1 + intij
         indint6 = indint5 + inck*ni1nj1

         sx1 = -surf(indint1,1)
         sx2 =  surf(indint2,1)
         sx3 = -surf(indint3,1)
         sx4 =  surf(indint4,1)
         sx5 = -surf(indint5,1)
         sx6 =  surf(indint6,1)

         sy1 = -surf(indint1,2)
         sy2 =  surf(indint2,2)
         sy3 = -surf(indint3,2)
         sy4 =  surf(indint4,2)
         sy5 = -surf(indint5,2)
         sy6 =  surf(indint6,2)

         sz1 = -surf(indint1,3)
         sz2 =  surf(indint2,3)
         sz3 = -surf(indint3,3)
         sz4 =  surf(indint4,3)
         sz5 = -surf(indint5,3)
         sz6 =  surf(indint6,3)

         vinv = ONE/MAX(vol(indcell), E_MIN_VOL)

         gradxx = sx1 * fldintX(indint1)+sx2 * fldintX(indint2)
     &          + sx3 * fldintX(indint3)+sx4 * fldintX(indint4)
     &          + sx5 * fldintX(indint5)+sx6 * fldintX(indint6)

         gradyy = sy1 * fldintY(indint1)+sy2 * fldintY(indint2)
     &          + sy3 * fldintY(indint3)+sy4 * fldintY(indint4)
     &          + sy5 * fldintY(indint5)+sy6 * fldintY(indint6)

         gradzz = sz1 * fldintZ(indint1)+sz2 * fldintZ(indint2)
     &          + sz3 * fldintZ(indint3)+sz4 * fldintZ(indint4)
     &          + sz5 * fldintZ(indint5)+sz6 * fldintZ(indint6)

         div(indcell) = vinv *(gradxx + gradyy + gradzz)

      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      END

C=============================================================================
