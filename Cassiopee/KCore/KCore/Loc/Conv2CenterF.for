C  
C    Copyright 2013-2026 ONERA.
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
C  convert field of floats defined in nodes to field defined in centers,
C  3D
C  ============================================================================
      SUBROUTINE k6conv2center1(ni, nj, nk, nfld, fieldnode, 
     &                          fieldcenter)
      IMPLICIT NONE
      
#include "Def/DefFortranConst.h"      
C==============================================================================
C_IN
      INTEGER_E ni, nj, nk, nfld
      REAL_E fieldnode(0:ni*nj*nk-1,nfld)
      
C_OUT
      REAL_E fieldcenter(0:(ni-1)*(nj-1)*(nk-1)-1,nfld)
      
C_LOCAL
      INTEGER_E i, j, k, n, indc
      INTEGER_E ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7
      INTEGER_E nic, njc, nicnjc, ninj
      REAL_E a
C==============================================================================

      a = ONE_EIGHT
      nic = ni-1
      njc = nj-1
      nicnjc = nic*njc
      ninj = ni*nj

!$OMP PARALLEL PRIVATE(i, indc, ind0, ind1, ind2, ind3, ind4, ind5,
!$OMP& ind6, ind7, n, k, j)      
      DO n = 1, nfld
!$OMP DO
        DO indc = 0, nicnjc*(nk-1)-1
            k = indc / nicnjc
            j = (indc-k*nicnjc)/nic
            i = indc - j*nic - k*nicnjc
            ind0 = i+j*ni+k*ninj
            ind1 = ind0+1
            ind2 = ind0+ni
            ind3 = ind2+1
            ind4 = ind0+ninj
            ind5 = ind4+1
            ind6 = ind4+ni
            ind7 = ind6+1
              
            fieldcenter(indc,n) = a*(fieldnode(ind0,n)+
     &                 fieldnode(ind1,n)+
     &                 fieldnode(ind2,n)+
     &                 fieldnode(ind3,n)+
     &                 fieldnode(ind4,n)+
     &                 fieldnode(ind5,n)+
     &                 fieldnode(ind6,n)+
     &                 fieldnode(ind7,n))
        ENDDO
!$OMP ENDDO
      ENDDO
!$OMP END PARALLEL

      END

C  ============================================================================
C  Convert cellnature field defined in nodes to field defined in centers
C  cellN=0 (blanked or interpolated) or 1 (normal)         
C  Si au moins un 0 ==> 0, sinon 1
C  ============================================================================
      SUBROUTINE k6conv2center21(ni, nj, nk, fieldnode, fieldcenter)
      IMPLICIT NONE

C==============================================================================
C_IN
      INTEGER_E ni, nj, nk
      REAL_E fieldnode(0:ni*nj*nk-1)
      
C_OUT
      REAL_E fieldcenter(0:(ni-1)*(nj-1)*(nk-1)-1)
      
C_LOCAL
      INTEGER_E i, j, k, indc
      INTEGER_E ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7
      INTEGER_E nic, njc, nicnjc, ninj
      REAL_E f
C==============================================================================
      nic = ni-1
      njc = nj-1
      nicnjc = nic*njc
      ninj = ni*nj
!$OMP PARALLEL PRIVATE(i, j, indc, ind0, ind1, ind2, ind3, ind4, ind5,
!$OMP& ind6, ind7, f, k)
!$OMP DO
      DO indc = 0, nicnjc*(nk-1)-1
        k = indc / nicnjc
        j = (indc-k*nicnjc)/nic
        i = indc - j*nic - k*nicnjc
            
        ind0 = i+j*ni+k*ninj
        ind1 = ind0+1
        ind2 = ind0+ni
        ind3 = ind2+1
        ind4 = ind0+ninj
        ind5 = ind4+1
        ind6 = ind4+ni
        ind7 = ind6+1
        f = MIN(fieldnode(ind0), fieldnode(ind1))
        f = MIN(f, fieldnode(ind2))
        f = MIN(f, fieldnode(ind3))
        f = MIN(f, fieldnode(ind4))
        f = MIN(f, fieldnode(ind5))
        f = MIN(f, fieldnode(ind6))
        fieldcenter(indc) = MIN(f, fieldnode(ind7))
      ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
      END

C  ============================================================================
C  convert cellnature field defined in nodes to field defined in centers
C  cellN=0 (blanked), 1 (normal) or 2 (interpolated)        
C      Dans le domaine interieur: Tous 0 ==> 0
C                                 Si au moins un 2 ==> 2
C                                 Sinon 1
C      Sur la frontiere: Tous 0 ==> 0
C                        Tous 2 ==> 2
C                        Sinon 1
C  ============================================================================
      SUBROUTINE k6conv2center22(ni, nj, nk, fieldnode, fieldcenter)
      IMPLICIT NONE
C==============================================================================
C_IN
      INTEGER_E ni, nj, nk
      REAL_E fieldnode(0:ni*nj*nk-1)
      
C_OUT
      REAL_E fieldcenter(0:(ni-1)*(nj-1)*(nk-1)-1)
      
C_LOCAL
      INTEGER_E i, j, k, indc
      INTEGER_E ind0, ind1, ind2, ind3, ind4, ind5, ind6, ind7
      INTEGER_E nic, njc, nicnjc, ninj
      REAL_E somme, somme2
C==============================================================================
      nic = ni-1
      njc = nj-1
      nicnjc = nic*njc
      ninj = ni*nj
!$OMP PARALLEL PRIVATE(i, j, indc, ind0, ind1, ind2, ind3, ind4, ind5,
!$OMP& ind6, ind7, somme, somme2, k)
!$OMP DO
      DO indc = 0, nicnjc*(nk-1)-1
        k = indc / nicnjc
        j = (indc-k*nicnjc)/nic
        i = indc - j*nic - k*nicnjc
                
        ind0 = i+j*ni+k*ninj
        ind1 = ind0+1
        ind2 = ind0+ni
        ind3 = ind2+1
        ind4 = ind0+ninj
        ind5 = ind4+1
        ind6 = ind4+ni
        ind7 = ind6+1
        somme = fieldnode(ind0)+fieldnode(ind1)
     &              +fieldnode(ind2)+fieldnode(ind3)
     &              +fieldnode(ind4)+fieldnode(ind5)
     &              +fieldnode(ind6)+fieldnode(ind7)
        IF (somme .EQ. 0.) THEN
            fieldcenter(indc) = 0.D0
        ELSE
            somme2 = MAX(fieldnode(ind0),1.)
     &                 + MAX(fieldnode(ind1),1.)
     &                 + MAX(fieldnode(ind2),1.)
     &                 + MAX(fieldnode(ind3),1.)
     &                 + MAX(fieldnode(ind4),1.)
     &                 + MAX(fieldnode(ind5),1.)
     &                 + MAX(fieldnode(ind6),1.)
     &                 + MAX(fieldnode(ind7),1.)
            IF (somme2 .GT. 8.5) THEN
                fieldcenter(indc) = 2.D0
            ELSE
                fieldcenter(indc) = 1.D0
            ENDIF
        ENDIF
      ENDDO
!$OMP ENDDO
!$OMP END PARALLEL
      
C     Traitement des fontieres
C     i = 0
      DO k = 0, nk-2
         DO j = 0, nj-2
            indc = j*nic+k*nicnjc
            ind0 = j*ni+k*ninj
            ind1 = ind0+1
            ind2 = ind0+ni
            ind3 = ind2+1
            ind4 = ind0+ninj
            ind5 = ind4+1
            ind6 = ind4+ni
            ind7 = ind6+1
            somme = fieldnode(ind0)+fieldnode(ind1)
     &           +fieldnode(ind2)+fieldnode(ind3)
     &           +fieldnode(ind4)+fieldnode(ind5)
     &           +fieldnode(ind6)+fieldnode(ind7)
            IF (somme .EQ. 0.) THEN
               fieldcenter(indc) = 0.D0
            ELSE
               IF (somme .EQ. 16.) THEN
                  fieldcenter(indc) = 2.D0
               ELSE
                  fieldcenter(indc) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDDO

C     i = ni-2
      IF (ni > 2) THEN
         i = ni-2
         DO k = 0, nk-2
            DO j = 0, nj-2
               indc = i+j*nic+k*nicnjc
               ind0 = i+j*ni+k*ninj
               ind1 = ind0+1
               ind2 = ind0+ni
               ind3 = ind2+1
               ind4 = ind0+ninj
               ind5 = ind4+1
               ind6 = ind4+ni
               ind7 = ind6+1
               somme = fieldnode(ind0)+fieldnode(ind1)
     &              +fieldnode(ind2)+fieldnode(ind3)
     &              +fieldnode(ind4)+fieldnode(ind5)
     &              +fieldnode(ind6)+fieldnode(ind7)
               IF (somme .EQ. 0.) THEN
                  fieldcenter(indc) = 0.D0
               ELSE
                  IF (somme .EQ. 16.) THEN
                     fieldcenter(indc) = 2.D0
                  ELSE
                     fieldcenter(indc) = 1.D0
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      
C     j = 0
      DO k = 0, nk-2
         DO i = 0, ni-2
            indc = i+k*nicnjc
            ind0 = i+k*ninj
            ind1 = ind0+1
            ind2 = ind0+ni
            ind3 = ind2+1
            ind4 = ind0+ninj
            ind5 = ind4+1
            ind6 = ind4+ni
            ind7 = ind6+1
            somme = fieldnode(ind0)+fieldnode(ind1)
     &           +fieldnode(ind2)+fieldnode(ind3)
     &           +fieldnode(ind4)+fieldnode(ind5)
     &           +fieldnode(ind6)+fieldnode(ind7)
            IF (somme .EQ. 0.) THEN
               fieldcenter(indc) = 0.D0
            ELSE
               IF (somme .EQ. 16.) THEN
                  fieldcenter(indc) = 2.D0
               ELSE
                  fieldcenter(indc) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      
      IF (nj > 2) THEN
C     j = nj-2
         j = nj-2
         DO k = 0, nk-2
            DO i = 0, ni-2
               indc = i+j*nic+k*nicnjc
               ind0 = i+j*ni+k*ninj
               ind1 = ind0+1
               ind2 = ind0+ni
               ind3 = ind2+1
               ind4 = ind0+ninj
               ind5 = ind4+1
               ind6 = ind4+ni
               ind7 = ind6+1
               somme = fieldnode(ind0)+fieldnode(ind1)
     &              +fieldnode(ind2)+fieldnode(ind3)
     &              +fieldnode(ind4)+fieldnode(ind5)
     &              +fieldnode(ind6)+fieldnode(ind7)
               IF (somme .EQ. 0.) THEN
                  fieldcenter(indc) = 0.D0
               ELSE
                  IF (somme .EQ. 16.) THEN
                     fieldcenter(indc) = 2.D0
                  ELSE
                     fieldcenter(indc) = 1.D0
                  ENDIF
               ENDIF
             ENDDO
         ENDDO
      ENDIF

C     k = 0
      DO j = 0, nj-2
         DO i = 0, ni-2
            indc = i+j*nic
            ind0 = i+j*ni
            ind1 = ind0+1
            ind2 = ind0+ni
            ind3 = ind2+1
            ind4 = ind0+ninj
            ind5 = ind4+1
            ind6 = ind4+ni
            ind7 = ind6+1
            somme = fieldnode(ind0)+fieldnode(ind1)
     &           +fieldnode(ind2)+fieldnode(ind3)
     &           +fieldnode(ind4)+fieldnode(ind5)
     &           +fieldnode(ind6)+fieldnode(ind7)
            IF (somme .EQ. 0.) THEN
               fieldcenter(indc) = 0.D0
            ELSE
               IF (somme .EQ. 16.) THEN
                  fieldcenter(indc) = 2.D0
               ELSE
                  fieldcenter(indc) = 1.D0
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      IF (nk > 2) THEN
         
C     k = nk-1
         k = nk-2
         DO j = 0, nj-2
            DO i = 0, ni-2
               indc = i+j*nic+k*nicnjc
               ind0 = i+j*ni+k*ninj
               ind1 = ind0+1
               ind2 = ind0+ni
               ind3 = ind2+1
               ind4 = ind0+ninj
               ind5 = ind4+1
               ind6 = ind4+ni
               ind7 = ind6+1
               somme = fieldnode(ind0)+fieldnode(ind1)
     &              +fieldnode(ind2)+fieldnode(ind3)
     &              +fieldnode(ind4)+fieldnode(ind5)
     &              +fieldnode(ind6)+fieldnode(ind7)
               IF (somme .EQ. 0.) THEN
                  fieldcenter(indc) = 0.D0
               ELSE
                  IF (somme.EQ. 16.) THEN
                     fieldcenter(indc) = 2.D0
                  ELSE
                     fieldcenter(indc) = 1.D0
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      END
C ===== KCore/Fortran/Conv2CenterF.for =====
