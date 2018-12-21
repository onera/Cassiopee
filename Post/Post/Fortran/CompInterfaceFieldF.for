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

C Etant donnes n champs definis aux noeuds d une grille 3D , 
C calcul des champs aux centres des interfaces de la grille  
C fint = 0.25*(fa+fb+fc+fd)
C
C ============================================================================
C
      SUBROUTINE k6compintfield(ni, nj, nk, nfld, nbint, f, fint)
     
      IMPLICIT NONE

# include "Def/DefFortranConst.h"
C==============================================================================
C_IN 
      INTEGER_E ni, nj, nk      ! dimensions aux noeuds
      INTEGER_E nfld            ! nb de champs
      INTEGER_E nbint           ! nb d interfaces total
      REAL_E    f(0:ni*nj*nk-1,nfld) ! champs aux noeuds
C_OUT
      REAL_E    fint(0:nbint-1,nfld) ! champs aux interfaces
C_LOCAL
      INTEGER_E eq
      INTEGER_E i, j, k, ninj, ni1, nj1, nk1, ni1nj1, ninj1, ni1nj
      INTEGER_E l1, l2, l3, l4
      INTEGER_E li, lj, lk
      INTEGER_E incij, incik, incijk ! increments en i
      INTEGER_E incjk, incji, incjik ! increments en j
      INTEGER_E incki, inckj, inckij ! increments en k
      INTEGER_E inti, intj, intij
C==============================================================================
C_ CONSTANTES
      ninj = ni*nj
      ni1 = ni-1
      nj1 = nj-1
      nk1 = nk-1
      ni1nj1 = ni1*nj1
      ninj1 = ni*nj1
      ni1nj = ni1*nj
C
      inti = ninj1*nk1
      intj = ni1nj*nk1 
      intij = inti + intj

C     increments pour interfaces en i
      incij  = ni           
      incik  = ninj        
      incijk = ni + ninj    

C     increments pour interfaces en j
      incjk = ninj
      incji = 1
      incjik = 1 + ninj

C     increments pour interfaces en k
      incki  = 1            
      inckj  = ni          
      inckij = 1 + ni      

C     Interfaces en i

      DO k = 0, nk-2
      DO j = 0, nj-2
      DO i = 0, ni1
         li = i + j*ni + k*ninj1
         l1 = i + j*ni + k*ninj
         l2 = l1 + incij
         l3 = l1 + incijk
         l4 = l1 + incik
         DO eq = 1, nfld 
            fint(li,eq) = ONE_FOURTH * 
     &           (f(l1,eq) + f(l2,eq) + f(l3,eq) + f(l4,eq) )
         ENDDO
      ENDDO
      ENDDO
      ENDDO

C     Interfaces en j
      DO k = 0, nk-2
      DO j = 0, nj1
      DO i = 0, ni-2
         lj = i + j*ni1 + k*ni1nj + inti
         l1 = i + j*ni + k*ninj
         l2 = l1 + incjk
         l3 = l1 + incjik
         l4 = l1 + incji
         DO eq = 1, nfld 
            fint(lj,eq) = ONE_FOURTH * 
     &           (f(l1,eq) + f(l2,eq) + f(l3,eq) + f(l4,eq) )
         ENDDO
      ENDDO
      ENDDO
      ENDDO

C     Interfaces en k
      DO k = 0, nk1
      DO j = 0, nj-2
      DO i = 0, ni-2
         lk = i + j*ni1 + k*ni1nj1 + intij
         l1 = i + j*ni + k*ninj
         l2 = l1 + incki
         l3 = l1 + inckij
         l4 = l1 + inckj
         DO eq = 1, nfld 
            fint(lk,eq) = ONE_FOURTH * 
     &           (f(l1,eq) + f(l2,eq) + f(l3,eq) + f(l4,eq) )
         ENDDO
      ENDDO 
      ENDDO 
      ENDDO
C       
      END
C=============================================================================
C Compute the values of the vector (f1,f2,f3) at interfaces
C=============================================================================
      SUBROUTINE k6compintfieldv(ni, nj, nk, nbint, f1, f2, f3,
     &                           fint1, fint2, fint3 )
     
      IMPLICIT NONE

# include "Def/DefFortranConst.h"

C_IN 
      INTEGER_E ni, nj, nk      ! dimensions aux noeuds
      INTEGER_E nbint           ! nb d interfaces total
      REAL_E    f1(0:ni*nj*nk-1) ! champs aux noeuds
      REAL_E    f2(0:ni*nj*nk-1) ! champs aux noeuds
      REAL_E    f3(0:ni*nj*nk-1) ! champs aux noeuds

C_OUT
      REAL_E    fint1(0:nbint-1) ! champs aux interfaces
      REAL_E    fint2(0:nbint-1) ! champs aux interfaces
      REAL_E    fint3(0:nbint-1) ! champs aux interfaces

C_LOCAL
      INTEGER_E eq
      INTEGER_E i, j, k, ninj, ni1, nj1, nk1, ni1nj1, ninj1, ni1nj
      INTEGER_E l1, l2, l3, l4
      INTEGER_E li, lj, lk
      INTEGER_E incij, incik, incijk ! increments en i
      INTEGER_E incjk, incji, incjik ! increments en j
      INTEGER_E incki, inckj, inckij ! increments en k
      INTEGER_E inti, intj, intij
C
C -------------------------------------------------------
C_ CONSTANTES
      ninj = ni*nj
      ni1 = ni-1
      nj1 = nj-1
      nk1 = nk-1
      ni1nj1 = ni1*nj1
      ninj1 = ni*nj1
      ni1nj = ni1*nj
C
      inti = ninj1*nk1
      intj = ni1nj*nk1 
      intij = inti + intj

C     increments pour interfaces en i
      incij  = ni           
      incik  = ninj        
      incijk = ni + ninj    

C     increments pour interfaces en j
      incjk = ninj
      incji = 1
      incjik = 1 + ninj

C     increments pour interfaces en k
      incki  = 1            
      inckj  = ni          
      inckij = 1 + ni      

C     Interfaces en i

      DO k = 0, nk-2
      DO j = 0, nj-2
      DO i = 0, ni1
         li = i + j*ni + k*ninj1
         l1 = i + j*ni + k*ninj
         l2 = l1 + incij
         l3 = l1 + incijk
         l4 = l1 + incik
         fint1(li) = ONE_FOURTH * (f1(l1)+f1(l2)+f1(l3)+f1(l4))
         fint2(li) = ONE_FOURTH * (f2(l1)+f2(l2)+f2(l3)+f2(l4))
         fint3(li) = ONE_FOURTH * (f3(l1)+f3(l2)+f3(l3)+f3(l4))
      ENDDO
      ENDDO
      ENDDO

C     Interfaces en j
      DO k = 0, nk-2
      DO j = 0, nj1
      DO i = 0, ni-2
         lj = i + j*ni1 + k*ni1nj + inti
         l1 = i + j*ni + k*ninj
         l2 = l1 + incjk
         l3 = l1 + incjik
         l4 = l1 + incji
         fint1(lj) = ONE_FOURTH * (f1(l1)+f1(l2)+f1(l3)+f1(l4))
         fint2(lj) = ONE_FOURTH * (f2(l1)+f2(l2)+f2(l3)+f2(l4))
         fint3(lj) = ONE_FOURTH * (f3(l1)+f3(l2)+f3(l3)+f3(l4))
      ENDDO
      ENDDO
      ENDDO

C     Interfaces en k
      DO k = 0, nk1
      DO j = 0, nj-2
      DO i = 0, ni-2
         lk = i + j*ni1 + k*ni1nj1 + intij
         l1 = i + j*ni + k*ninj
         l2 = l1 + incki
         l3 = l1 + inckij
         l4 = l1 + inckj
         fint1(lk) = ONE_FOURTH * (f1(l1)+f1(l2)+f1(l3)+f1(l4))
         fint2(lk) = ONE_FOURTH * (f2(l1)+f2(l2)+f2(l3)+f2(l4))
         fint3(lk) = ONE_FOURTH * (f3(l1)+f3(l2)+f3(l3)+f3(l4))
      ENDDO 
      ENDDO 
      ENDDO
C       
      END
