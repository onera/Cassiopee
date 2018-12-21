C  
C    Copyright 2013-2019 Onera.
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

C  ============================================================================
C  Calcul du rotationnel moyen d un champ u sur une grille
C  IN: ni,nj,nk: dimensions du maillage en noeuds
C  IN: nb./Post/Post/Fortran/CompStructCurlF.for de cellules
C  IN: xt,yt,zt: coordonnees de la grille
C  IN: u: vecteur dont le rotationnel est a calculer
C  OUT: rotu: rotationnel de u aux centres des cellules
C  ============================================================================
      SUBROUTINE k6compstructcurlt( ni, nj, nk, nbcell, nbint, 
     &     xt, yt, zt, ux, uy, uz, rotx, roty, rotz, 
     &     surf, snorm, centerInt, vol,
     &     uintx, uinty, uintz) 

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E nbcell          ! nb de cellules dans la grille
      INTEGER_E ni, nj, nk      ! dimensions de la grille contenant la cellule
      INTEGER_E nbint           ! nb d interfaces
      REAL_E xt(0:ni*nj*nk-1)   ! coordonnees aux noeuds
      REAL_E yt(0:ni*nj*nk-1)   ! coordonnees aux noeuds
      REAL_E zt(0:ni*nj*nk-1)   ! coordonnees aux noeuds
      REAL_E ux(0:ni*nj*nk-1)  ! vecteur dont le rotationnel est a calculer 
      REAL_E uy(0:ni*nj*nk-1)  ! vecteur dont le rotationnel est a calculer 
      REAL_E uz(0:ni*nj*nk-1)  ! vecteur dont le rotationnel est a calculer 

C_OUT
      REAL_E  rotx(0:nbcell-1) ! rotationnel de ux
      REAL_E  roty(0:nbcell-1) ! rotationnel de uy
      REAL_E  rotz(0:nbcell-1) ! rotationnel de uz
      REAL_E  surf(0:nbint-1,3) ! surfaces
      REAL_E  snorm(0:nbint-1)  ! norme de surf : pas utilise
      REAL_E  centerInt(0:nbint-1,3)
      REAL_E  vol(0:nbcell-1)
      REAL_E  uintx(0:nbint-1) ! vecteur u aux interfaces  
      REAL_E  uinty(0:nbint-1) ! vecteur u aux interfaces  
      REAL_E  uintz(0:nbint-1) ! vecteur u aux interfaces  

C_LOCAL
      INTEGER_E inci, incj, inck
      INTEGER_E indcell, i, j, k, eq
      INTEGER_E ninj, ni1, nj1, nk1, ni1nj, ninj1, ni1nj1
      INTEGER_E indint1, indint2, indint3, indint4, indint5, indint6
      INTEGER_E inti, intj, intij, intk
      REAL_E   curlx, curly, curlz
      REAL_E sx1, sx2, sx3, sx4, sx5, sx6
      REAL_E sy1, sy2, sy3, sy4, sy5, sy6
      REAL_E sz1, sz2, sz3, sz4, sz5, sz6
      REAL_E  vinv              ! inverse du volume d une cellule
C
C-----------------------------------------------------------------------------
C_CONSTANTES
      ni1 = ni-1
      nj1 = nj-1
      nk1 = nk-1
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
C ------------------------------------------------------------------
C attention : surf n est pas oriente : tjs positif
      call k6compstructmetric(ni, nj, nk, nbcell, nbint,inti,intj,intk,
     &     xt, yt, zt, vol,  surf, surf(0,2), surf(0,3), snorm, 
     &     centerInt, centerInt(0,2), centerInt(0,3) )
       
      call k6compintfieldv(ni, nj, nk, nbint, ux, uy, uz, 
     &                     uintx, uinty, uintz)

!$OMP PARALLEL PRIVATE(indcell,i,j,k,indint1,indint2,indint3)
!$OMP& PRIVATE(indint4,indint5,indint6,sx1,sx2,sx3,sx4,sx5,sx6)
!$OMP& PRIVATE(sy1,sy2,sy3,sy4,sy5,sy6,sz1,sz2,sz3,sz4,sz5,sz6)
!$OMP& PRIVATE(vinv,curlx,curly,curlz)
!$OMP DO
      DO indcell = 0, nbcell-1
         k = indcell / ni1nj1   !  demarrent a 0
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

         vinv = -ONE/vol(indcell)
      
         curlx = 
     &        uinty(indint1)* sz1 - uintz(indint1)* sy1 +  
     &        uinty(indint2)* sz2 - uintz(indint2)* sy2 +  
     &        uinty(indint3)* sz3 - uintz(indint3)* sy3 +
     &        uinty(indint4)* sz4 - uintz(indint4)* sy4 + 
     &        uinty(indint5)* sz5 - uintz(indint5)* sy5 +
     &        uinty(indint6)* sz6 - uintz(indint6)* sy6 

         curly = 
     &        uintz(indint1)* sx1 - uintx(indint1)* sz1 +  
     &        uintz(indint2)* sx2 - uintx(indint2)* sz2 +  
     &        uintz(indint3)* sx3 - uintx(indint3)* sz3 +
     &        uintz(indint4)* sx4 - uintx(indint4)* sz4 + 
     &        uintz(indint5)* sx5 - uintx(indint5)* sz5 +
     &        uintz(indint6)* sx6 - uintx(indint6)* sz6 

         curlz = 
     &        uintx(indint1)* sy1 - uinty(indint1)* sx1 +  
     &        uintx(indint2)* sy2 - uinty(indint2)* sx2 +  
     &        uintx(indint3)* sy3 - uinty(indint3)* sx3 +
     &        uintx(indint4)* sy4 - uinty(indint4)* sx4 + 
     &        uintx(indint5)* sy5 - uinty(indint5)* sx5 +
     &        uintx(indint6)* sy6 - uinty(indint6)* sx6 

         rotx(indcell) = vinv * curlx
         roty(indcell) = vinv * curly
         rotz(indcell) = vinv * curlz
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      END

C=============== Post/Fortran/CompStructCurlF.for =============================
