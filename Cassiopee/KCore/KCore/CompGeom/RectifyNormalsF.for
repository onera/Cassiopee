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

C ============================================================================
C Redresse la normale du bloc 2 dans le meme sens que le bloc 1 
C Test du produit scalaire des normales aux cellules ind1 de bloc 1 
C et ind2 de bloc 2 : 
C retourne 0 si les normales des 2 blocs sont orientees dans le meme sens 
C retourne 1 si orientations des normales opposees
C  ===========================================================================
      SUBROUTINE k6rectifynormals(ni1, nj1, ind1, x1, y1, z1,
     &                            ni2, nj2, ind2, x2, y2, z2,
     &                            distmin, notvalid, isopp)
C
      IMPLICIT NONE 
C 
#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E ni1, nj1, ind1
      INTEGER_E ni2, nj2, ind2
      INTEGER_E notvalid
      REAL_E x1(0:ni1*nj1-1)
      REAL_E y1(0:ni1*nj1-1)
      REAL_E z1(0:ni1*nj1-1)
      REAL_E x2(0:ni2*nj2-1)
      REAL_E y2(0:ni2*nj2-1)
      REAL_E z2(0:ni2*nj2-1)
      REAL_E distmin            ! distance entre ind1 et ind2
C_OUT
      INTEGER_E isopp
C_LOCAL
      INTEGER_E ind11, ind12, ind13, ind14
      INTEGER_E ind21, ind22, ind23, ind24
      REAL_E eps, ps1, ps2, matchTol
      REAL_E x13, y13, z13, x24, y24, z24
      REAL_E n1x, n1y, n1z, n2x, n2y, n2z
      REAL_E norm1, norm2, norm, normab, normac, normp
      INTEGER_E i1,j1, i2, j2, inci1, incj1, inci2, incj2
      INTEGER_E ind1p, ind2p
      INTEGER_E im1, im2, jm1, jm2, ip1, ip2, jp1, jp2
      REAL_E xab, yab, zab, xac, yac, zac
      REAL_E psmin, ps3
      REAL_E v1x, v1y, v1z, v2x, v2y, v2z
      INTEGER_E iclassical 
C==============================================================================
      notvalid = 0
      isOpp = 1
C     constante  : si <n1,n2> < eps -> a l envers ->isCorr = 1 
      eps = 1.e-12
      matchTol = 1.e-12
      psmin = 0.1 ! minimum pour ne pas considerer les blocs orthogonaux

      j1 = ind1/ni1 + 1
      i1 = ind1 - (j1-1)*ni1 + 1
      j2 = ind2/ni2 + 1
      i2 = ind2 - (j2-1)*ni2 + 1      
      inci1 = 0
      incj1 = 0
      inci2 = 0
      incj2 = 0

      if (i1 .eq. ni1) then 
         inci1 = -1
      endif

      if (j1 .eq. nj1) then 
         incj1 = -1
      endif
C     
      if (i2 .eq. ni2) then
         inci2 = -1
      endif
      if (j2 .eq. nj2) then 
         incj2 = -1
      endif   

C     Traitement de A1B1C1D1
C     indices du quad A1B1C1D1
      im1 = i1 + inci1 - 1
      ip1 = i1 + inci1 
      jm1 = j1 + incj1 - 1
      jp1 = j1 + incj1 
      
      ind11 = im1 + jm1 * ni1   ! A1(i,j) 
      ind12 = ip1 + jm1 * ni1   ! B1(i+1,j)
      ind13 = ip1 + jp1 * ni1   ! C1(i+1,j+1)
      ind14 = im1 + jp1 * ni1   ! D1(i,j+1)

C     A1C1
      x13 = x1(ind13) - x1(ind11)
      y13 = y1(ind13) - y1(ind11)
      z13 = z1(ind13) - z1(ind11)
      
C     B1D1
      x24 = x1(ind14) - x1(ind12)
      y24 = y1(ind14) - y1(ind12)
      z24 = z1(ind14) - z1(ind12)

C     A1C1 x B1D1
      
      n1x = y13 * z24 - y24 * z13
      n1y = x24 * z13 - x13 * z24
      n1z = x13 * y24 - x24 * y13
            
      norm1 = n1x*n1x + n1y*n1y + n1z*n1z
      
C     Traitement de A2B2C2D2
C     indices du quad A2B2C2D2
      im2 = i2 + inci2 - 1
      ip2 = i2 + inci2 
      jm2 = j2 + incj2 - 1
      jp2 = j2 + incj2 

      ind21 = im2 + jm2 * ni2   ! A2(i,j) 
      ind22 = ip2 + jm2 * ni2   ! B2(i+1,j)
      ind23 = ip2 + jp2 * ni2   ! C2(i+1,j+1)
      ind24 = im2 + jp2 * ni2   ! D2(i,j+1)

C     A2C2
      x13 = x2(ind23) - x2(ind21)
      y13 = y2(ind23) - y2(ind21)
      z13 = z2(ind23) - z2(ind21)
      
C     B2D2
      x24 = x2(ind24) - x2(ind22)
      y24 = y2(ind24) - y2(ind22)
      z24 = z2(ind24) - z2(ind22)

C     A2C2 x B2D2
      
      n2x = y13 * z24 - y24 * z13
      n2y = x24 * z13 - x13 * z24
      n2z = x13 * y24 - x24 * y13

      norm2 = n2x*n2x + n2y*n2y + n2z*n2z

      norm = sqrt(norm1*norm2)
      norm = ONE/norm
      ps1 = ( n1x * n2x + n1y * n2y + n1z * n2z ) * norm 

      IF (distmin .gt. matchTol) THEN ! recouvrants
         if (ps1 < -eps) then 
            isopp = -1
         endif
C-----------------------------------------------
      ELSE                      !  ou coincident ou recouvrant (si interieur)
         iclassical = 0
C     recherche des points opposes 
         IF (i1 .eq. 1 .and. j1 .ne. 1 .and. j1 .ne. nj1)  THEN 
            ind1p = ind1 + 1
         ELSE IF (i1 .eq. ni1 .and. j1 .ne. 1 .and. j1 .ne. nj1)  THEN 
            ind1p = ind1 - 1
         ELSE IF (j1 .eq. 1 .and. i1 .ne. 1 .and. i1 .ne. ni1)  THEN 
            ind1p = ind1 + ni1
         ELSE IF (j1 .eq. nj1 .and. i1 .ne. 1 .and. i1 .ne. ni1)  THEN 
            ind1p = ind1 - ni1   
         ELSE 
            iclassical = 1 ! point interieur
         ENDIF
         IF (i2 .eq. 1 .and. j2 .ne. 1 .and. j2 .ne. nj2)  THEN 
            ind2p = ind2 + 1
         ELSE IF (i2 .eq. ni2 .and. j2 .ne. 1 .and. j2 .ne. nj2)  THEN 
            ind2p = ind2 - 1
         ELSE IF (j2 .eq. 1 .and. i2 .ne. 1 .and. i2 .ne. ni2)  THEN 
            ind2p = ind2 + ni2
         ELSE IF (j2 .eq. nj2 .and. i2 .ne. 1 .and. i2 .ne. ni2)  THEN 
            ind2p = ind2 - ni2   
         ELSE 
            iclassical = 1 ! point interieur
         ENDIF

c$$$         write(*,*) i1, ni1, j1, nj1
c$$$         write(*,*) i2, ni2, j2, nj2
         
         IF (iclassical .eq. 0) THEN  ! point frontiere
C     ps2 = a0a1.a0a2 
            xab = x1(ind1p) - x1(ind1)
            yab = y1(ind1p) - y1(ind1)
            zab = z1(ind1p) - z1(ind1)
            
            xac = x2(ind2p) - x2(ind2)
            yac = y2(ind2p) - y2(ind2)
            zac = z2(ind2p) - z2(ind2)
            
c$$$            write(*,*) x1(ind1),y1(ind1),z1(ind1)
c$$$            write(*,*) x1(ind1p),y1(ind1p),z1(ind1p)
c$$$            write(*,*) x2(ind2p),y2(ind2p),z2(ind2p)
            
            normab = xab*xab + yab*yab + zab*zab
            normac = xac*xac + yac*yac + zac*zac
            normp = ONE/sqrt(normab*normac)
            ps2 = (xab*xac + yab*yac + zab*zac) * normp
            
c$$$            write(*,*) ' ps1=',ps1, ',ps2=', ps2

c$$$            if ( abs(ps2) .lt. psmin )  then 
c$$$               notvalid = 1
c$$$               return
c$$$            endif
            
            IF ((ps1 .gt. eps .and. ps2 .gt. eps) .or. 
     &          (ps1 .lt. eps .and. ps2 .lt. eps)) then
               
               v1x = yab*n1z - zab * n1y
               v1y = zab*n1x - xab * n1z
               v1z = xab*n1y - yab * n1x
               v2x = yac*n2z - zac * n2y
               v2y = zac*n2x - xac * n2z
               v2z = xac*n2y - yac * n2x

               ps3 = v1x*v2x + v1y*v2y + v1z*v2z
               ps3 = ps3 * norm * normp
               IF (abs(ps3) .gt. psmin) THEN
                  IF (ps3 .gt. 0.) THEN 
                     isOpp = -1
                  ENDIF
               ELSE             ! standard 
                  IF (ps1 < -eps) THEN 
                     isOpp = -1
                  ENDIF
               ENDIF
               
               WRITE(*,*) 'ambiguous matching treatment ', ps1, ps2, ps3

c$$$               IF (ABS(ps2) < psmin) THEN
c$$$                  WRITE(*,*) 'ambiguous matching treatment ', ps1, ps2
c$$$                  notvalid = 1
c$$$               ELSE
c$$$                  WRITE(*,*) 'correct matching treatment ', ps1, ps2
c$$$               ENDIF
c$$$               isOpp  = -1
            ENDIF
            
         ELSE                   ! traitement classique
            write(*,*)  'Certainly overset ', ps1 
            write(*,*) i1, ni1, j1, nj1
            write(*,*) i2, ni2, j2, nj2
            if ( ps1 < -eps ) then 
               isopp = -1          
            ENDIF
         ENDIF                  ! classical
      ENDIF                     ! if match
      RETURN 
      END
