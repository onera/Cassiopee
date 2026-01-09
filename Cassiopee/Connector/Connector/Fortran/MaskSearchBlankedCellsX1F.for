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

C center_in + depth = 1 cas 3D

      SUBROUTINE k6searchblankedcellsx1( ni, nj, nk, 
     &     meshX, meshY, meshZ,
     &     xmin, ymin,
     &     niray, njray,
     &     hiray, hjray,
     &     indir, nz, z, 
     &     cellNatureField, isMasked )
      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E ni, nj, nk   ! Dimension of the domain
      REAL_E    meshX(0:ni*nj*nk-1) ! Abcisses of mesh points
      REAL_E    meshY(0:ni*nj*nk-1) ! Ordinates of mesh points
      REAL_E    meshZ(0:ni*nj*nk-1) ! Height of mesh points
      REAL_E    xmin, ymin      ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray, hjray    ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      REAL_E    z(0:nz-1)       ! pierce points

C_OUT
      INTEGER_E cellNatureField(0:(ni-1)*(nj-1)*(nk-1)-1) 
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l       ! For loop in the three directions
      INTEGER_E ip, jp, kp
      INTEGER_E ind, indray
      INTEGER_E indrayc, indraycp1, indraycp2, indraycp3
      REAL_E    dx1, dy1, dz1
      REAL_E    xp, yp, zp
      REAL_E    xp1, xp2, xp3, xp4, xp5, xp6, xp7, xp8
      REAL_E    yp1, yp2, yp3, yp4, yp5, yp6, yp7, yp8
      REAL_E    zp1, zp2, zp3, zp4, zp5, zp6, zp7, zp8
      REAL_E    xmincell, ymincell, zmincell
      REAL_E    xmaxcell, ymaxcell, zmaxcell
      INTEGER_E ibeg, iend
      INTEGER_E iray, irayc, iraycp1, iraycp2, iraycp3
      INTEGER_E jray, jrayc, jraycp1, jraycp2, jraycp3
      INTEGER_E iraymin, iraymax, jraymin, jraymax
      REAL_E ibmin, ibmax, ibc1, ibc2, ibc3, ibc4
      REAL_E xrayc, yrayc, alpha, beta, xmax, ymax
      INTEGER_E npmin, npmax, np, ex, pos
      INTEGER_E cellN, nij, nijray, nic, nicnjc
      REAL_E    xi(6), yi(6), zi(6)
      REAL_E    ibi(6), ib0, z1, z2, z3, z4, zint
C
C=============================================================================
C
      nij = ni * nj
      nijray = niray * njray
      nic = ni-1
      nicnjc = nic*(nj-1)
      isMasked = 0

      xmax = xmin + (niray-1) * hiray
      ymax = ymin + (njray-1) * hjray 

      pos = 0

      DO k = 0, nk-2
         DO j = 0, nj-2
            DO i = 0, ni-2
               cellN = 0
C     
#include "../../Connector/Fortran/MaskCell3DF.for"
C--------------------------
C       Barycenter of cell
C--------------------------     
               xp = xp1 + xp2 + xp3 + xp4 + xp5 + xp6 + xp7 + xp8
               yp = yp1 + yp2 + yp3 + yp4 + yp5 + yp6 + yp7 + yp8 
               zp = zp1 + zp2 + zp3 + zp4 + zp5 + zp6 + zp7 + zp8
               xp = xp * ONE_EIGHT
               yp = yp * ONE_EIGHT
               zp = zp * ONE_EIGHT
C--------------------------  
C centers of the interfaces
C--------------------------               
C     interface  2376
               xi(1) = ONE_FOURTH * (xp2 + xp3 + xp6 + xp7) 
               yi(1) = ONE_FOURTH * (yp2 + yp3 + yp6 + yp7)
               zi(1) = ONE_FOURTH * (zp2 + zp3 + zp6 + zp7) 
               
C     interface 1485
               xi(2) = ONE_FOURTH * (xp1 + xp4 + xp5 + xp8) 
               yi(2) = ONE_FOURTH * (yp1 + yp4 + yp5 + yp8) 
               zi(2) = ONE_FOURTH * (zp1 + zp4 + zp5 + zp8)           

C     interface 1265
               xi(3) = ONE_FOURTH * (xp1 + xp2 + xp6 + xp5)
               yi(3) = ONE_FOURTH * (yp1 + yp2 + yp6 + yp5)  
               zi(3) = ONE_FOURTH * (zp1 + zp2 + zp6 + zp5)           

C     interface 4378
               xi(4) = ONE_FOURTH * (xp3 + xp4 + xp7 + xp8)
               yi(4) = ONE_FOURTH * (yp3 + yp4 + yp7 + yp8)  
               zi(4) = ONE_FOURTH * (zp3 + zp4 + zp7 + zp8)

C     interface 1234
               xi(5) = ONE_FOURTH * (xp1 + xp2 + xp3 + xp4) 
               yi(5) = ONE_FOURTH * (yp1 + yp2 + yp3 + yp4)
               zi(5) = ONE_FOURTH * (zp1 + zp2 + zp3 + zp4)
              
C     interface  5678
               xi(6) = ONE_FOURTH * (xp5 + xp6 + xp7 + xp8) 
               yi(6) = ONE_FOURTH * (yp5 + yp6 + yp7 + yp8)
               zi(6) = ONE_FOURTH * (zp5 + zp6 + zp7 + zp8)   
C          
C Ray index corresponding to point cell
               IF (xmaxcell.LT.xmin) GOTO 5
               IF (xmincell.GT.xmax) GOTO 5
               IF (ymaxcell.LT.ymin) GOTO 5
               IF (ymincell.GT.ymax) GOTO 5
               
               iraymin = (xmincell - xmin)/hiray
               iraymax = (xmaxcell - xmin)/hiray+1
               jraymin = (ymincell - ymin)/hjray
               jraymax = (ymaxcell - ymin)/hjray+1
               iraymin = MAX(0, iraymin)
               iraymax = MIN(niray-1, iraymax)   
               jraymin = MAX(0, jraymin)
               jraymax = MIN(njray-1, jraymax)
               
               DO jray = jraymin, jraymax
                  DO iray = iraymin, iraymax
                     
                     indray = iray + jray * niray
                     ibeg = indir(indray)
                     IF (indray.EQ.nijray-1) THEN
                        iend = nz-1
                     ELSE
                        iend = indir(indray+1)-1
                     ENDIF
                     npmin = 0
                     DO l = ibeg, iend
                        IF (z(l).GT.zmincell) THEN
                           GOTO 1
                        ENDIF
                        npmin = npmin+1
                     ENDDO
 1                   ibmin = npmin - (npmin/2)*2
                     npmax = 0
                     DO l = ibeg, iend
                        IF (z(l).GT.zmaxcell) THEN
                           GOTO 2
                        ENDIF
                        npmax = npmax+1
                     ENDDO
 2                   ibmax = npmax - (npmax/2)*2
                     
                     IF ((ibmin.NE.0).OR.
     &                   (ibmax.NE.0).OR.
     &                   (npmin.NE.npmax)) THEN
                        
C     Cell intersects, we now check if center is in
C     point B1: coeff (1-alpha)*(1-beta)
                        irayc = (xp - xmin)/hiray
                        jrayc = (yp - ymin)/hjray
                        IF (irayc.LT.0 .OR. irayc.GT.niray-1 .OR.
     &                      jrayc.LT.0 .OR. jrayc.GT.njray-1) THEN
                           ib0 = 0.
                           GOTO 61 
                        ENDIF
                        indrayc = irayc + jrayc * niray
                        ibeg = indir(indrayc)
                        IF (indrayc.EQ.nijray-1) THEN
                           iend = nz-1
                        ELSE
                           iend = indir(indrayc+1)-1
                        ENDIF
                        np = 0
                        DO l = ibeg, iend
                           IF (z(l).GT.zp) THEN
                              GOTO 3
                           ENDIF
                           np = np+1
                        ENDDO
                        IF (np .EQ. 0) GOTO 61

 3                      ibc1 = np - (np/2)*2
                        IF (l.EQ.ibeg) THEN
                           z1 = z(l)
                           IF (ibc1.EQ.1) pos = 1
                        ELSE IF (l.EQ.iend+1) THEN
                           z1 = z(l-1)
                           IF (ibc1.EQ.1) pos = -1
                        ELSE IF (ABS(zp-z(l)).LT.ABS(zp-z(l-1))) THEN
                           z1 = z(l)
                           IF (ibc1.EQ.1) pos = 1
                        ELSE
                           z1 = z(l-1)
                           IF (ibc1.EQ.1) pos = -1
                        ENDIF

C     point B2: coef alpha*(1-beta)
                        iraycp1 = irayc+1
                        jraycp1 = jrayc
                        iraycp1 = MIN(iraycp1, niray-1)
                        indraycp1 = iraycp1 + jraycp1 * niray
                        ibeg = indir(indraycp1)
                        IF (indraycp1.EQ.nijray-1) THEN
                           iend = nz-1
                        ELSE
                           iend = indir(indraycp1+1)-1
                        ENDIF
                        np = 0
                        DO l = ibeg, iend
                           IF (z(l).GT.zp) THEN
                              GOTO 4
                           ENDIF
                           np = np+1
                        ENDDO
                        IF ( np .EQ. 0 ) GOTO 61

 4                      ibc2 = np - (np/2)*2
                        IF (l.EQ.ibeg) THEN
                           z2 = z(l)
                           IF (ibc2.EQ.1) pos = 1
                        ELSE IF (l.EQ.iend+1) THEN
                           z2 = z(l-1)
                           IF (ibc2.EQ.1) pos = -1
                        ELSE IF (ABS(zp-z(l)).LT.ABS(zp-z(l-1))) THEN
                           z2 = z(l)
                           IF (ibc2.EQ.1) pos = 1
                        ELSE
                           z2 = z(l-1)
                           IF (ibc2.EQ.1) pos = -1
                        ENDIF

C     point B3:  coef alpha*beta
                        iraycp2 = irayc+1
                        jraycp2 = jrayc+1
                        iraycp2 = MIN(iraycp2, niray-1)
                        jraycp2 = MIN(jraycp2, njray-1)
                        indraycp2 = iraycp2 + jraycp2 * niray
                        ibeg = indir(indraycp2)
                        IF (indraycp2.EQ.nijray-1) THEN
                           iend = nz-1
                        ELSE
                           iend = indir(indraycp2+1)-1
                        ENDIF
                        np = 0
                        DO l = ibeg, iend
                           IF (z(l).GT.zp) THEN
                              GOTO 34
                           ENDIF
                           np = np+1
                        ENDDO
                        IF ( np .EQ. 0 ) GOTO 61

 34                     ibc3 = np - (np/2)*2
                        IF (l.EQ.ibeg) THEN
                           z3 = z(l)
                           IF (ibc3.EQ.1) pos = 1
                        ELSE IF (l.EQ.iend+1) THEN
                           z3 = z(l-1)
                           IF (ibc3.EQ.1) pos = -1
                        ELSE IF (ABS(zp-z(l)).LT.ABS(zp-z(l-1))) THEN
                           z3 = z(l)
                           IF (ibc3.EQ.1) pos = 1
                        ELSE
                           z3 = z(l-1)
                           IF (ibc3.EQ.1) pos = -1
                        ENDIF

C     point B4: coef (1-alpha)*beta
                        iraycp3 = irayc
                        jraycp3 = jrayc+1
                        jraycp3 = MIN(jraycp3, njray-1)
                        indraycp3 = iraycp3 + jraycp3 * niray
                        ibeg = indir(indraycp3)
                        IF (indraycp3.EQ.nijray-1) THEN
                           iend = nz-1
                        ELSE
                           iend = indir(indraycp3+1)-1
                        ENDIF
                        np = 0
                        DO l = ibeg, iend
                           IF (z(l).GT.zp) THEN
                              GOTO 44
                           ENDIF
                           np = np+1
                        ENDDO
                        IF ( np .EQ. 0 ) GOTO 61

 44                     ibc4 = np - (np/2)*2
                        IF (l.EQ.ibeg) THEN
                           z4 = z(l)
                           IF (ibc4.EQ.1) pos = 1
                        ELSE IF (l.EQ.iend+1) THEN
                           z4 = z(l-1)
                           IF (ibc4.EQ.1) pos = -1
                        ELSE IF (ABS(zp-z(l)).LT.ABS(zp-z(l-1))) THEN
                           z4 = z(l)
                           IF (ibc4.EQ.1) pos = 1
                        ELSE
                           z4 = z(l-1)
                           IF (ibc4.EQ.1) pos = -1
                        ENDIF
                        
                        xrayc = xmin + hiray * irayc
                        yrayc = ymin + hjray * jrayc
                        alpha = (xp - xrayc)/hiray
                        beta  = (yp - yrayc)/hjray
                        
C                        ib0 =(1.-alpha)*(1.-beta) * ibc +
C     &                       alpha*(1.-beta) * ibcp1 + 
C     &                       alpha*beta * ibcp2 +
C     &                       (1.-alpha)*beta * ibcp3

                        zint = (1.-alpha)*(1.-beta) * z1 +
     &                       alpha*(1.-beta) * z2 + 
     &                       alpha*beta * z3 +
     &                       (1.-alpha)*beta * z4

                        IF ( ibc1.EQ.1.AND.ibc2.EQ.1.AND.
     &                       ibc3.EQ.1.AND.ibc4.EQ.1) THEN
                           cellNatureField(i+j*nic+k*nicnjc) = -1
                           isMasked = 1
                           GOTO 5
                        ENDIF

                        IF ( ibc1.EQ.0.AND.ibc2.EQ.0.AND.
     &                       ibc3.EQ.0.AND.ibc4.EQ.0) THEN
                           GOTO 61
                        ENDIF

                        IF (pos.EQ.1.AND.zp.LE.zint) THEN
                           cellNatureField(i+j*nic+k*nicnjc) = -1
                           isMasked = 1
                           GOTO 5
                        ENDIF

                        IF (pos.EQ.-1.AND.zp.GE.zint) THEN
                           cellNatureField(i+j*nic+k*nicnjc) = -1
                           isMasked = 1
                           GOTO 5
                        ENDIF
 61                     CONTINUE
   

C     Also, check if EX points are in
                        DO ex = 1, 6
                           irayc = (xi(ex) - xmin)/hiray
                           jrayc = (yi(ex) - ymin)/hjray
                           IF (irayc.LT.0 .OR. irayc.GT.niray-1.OR.
     &                         jrayc.LT.0 .OR. jrayc.GT.njray-1) THEN 
                              ibi(ex) = 0.
                              GOTO 62
                           ENDIF
                   
C C1: ex point of coef (1-alpha)*(1-beta)
                           indrayc = irayc + jrayc * niray
                           ibeg = indir(indrayc)
                           IF (indrayc.EQ.nijray-1) THEN
                              iend = nz-1
                           ELSE
                              iend = indir(indrayc+1)-1
                           ENDIF
                           np = 0
                           DO l = ibeg, iend
                              IF (z(l).GT.zi(ex)) THEN
                                 GOTO 7
                              ENDIF
                              np = np+1
                           ENDDO
                           IF ( np .EQ. 0 ) GOTO 62

 7                         ibc1 = np - (np/2)*2
                           IF (l.EQ.ibeg) THEN
                              z1 = z(l)
                              IF (ibc1.EQ.1) pos = 1
                           ELSE IF (l.EQ.iend+1) THEN
                              z1 = z(l-1)
                              IF (ibc1.EQ.1) pos = -1  
                           ELSE IF (ABS(zi(ex)-z(l)).LT.
     &                             ABS(zi(ex)-z(l-1))) THEN
                              z1 = z(l)
                              IF (ibc1.EQ.1) pos = 1
                           ELSE
                              z1 = z(l-1)
                              IF (ibc1.EQ.1) pos = -1
                           ENDIF

C C2: ex point of coef alpha*(1-beta)
                           iraycp1 = irayc+1
                           jraycp1 = jrayc
                           iraycp1 = MIN(iraycp1, niray-1)
                           indraycp1 = iraycp1 + jraycp1 * niray
                           ibeg = indir(indraycp1)
                           IF (indraycp1.EQ.nijray-1) THEN
                              iend = nz-1
                           ELSE
                              iend = indir(indraycp1+1)-1
                           ENDIF
                           np = 0
                           DO l = ibeg, iend
                              IF (z(l).GT.zi(ex)) THEN
                                 GOTO 8
                              ENDIF
                              np = np+1
                           ENDDO
                           IF ( np .EQ. 0 ) GOTO 62

 8                         ibc2 = np - (np/2)*2
                           IF (l.EQ.ibeg) THEN
                              z2 = z(l)
                              IF (ibc2.EQ.1) pos = 1
                           ELSE IF (l.EQ.iend+1) THEN
                              z2 = z(l-1)
                              IF (ibc2.EQ.1) pos = -1
                           ELSE IF (ABS(zi(ex)-z(l)).LT.
     &                             ABS(zi(ex)-z(l-1))) THEN
                              z2 = z(l)
                              IF (ibc2.EQ.1) pos = 1
                           ELSE
                              z2 = z(l-1)
                              IF (ibc2.EQ.1) pos = -1
                           ENDIF

C C3: coef alpha*beta
                           iraycp2 = irayc+1
                           jraycp2 = jrayc+1
                           iraycp2 = MIN(iraycp2, niray-1)
                           jraycp2 = MIN(jraycp2, njray-1)
                           indraycp2 = iraycp2 + jraycp2 * niray
                           ibeg = indir(indraycp2)
                           IF (indraycp2.EQ.nijray-1) THEN
                              iend = nz-1
                           ELSE
                              iend = indir(indraycp2+1)-1
                           ENDIF
                           np = 0
                           DO l = ibeg, iend
                              IF (z(l).GT.zi(ex)) THEN
                                 GOTO 72
                              ENDIF
                              np = np+1
                           ENDDO
                           IF ( np .EQ. 0 ) GOTO 62

 72                        ibc3 = np - (np/2)*2
                           IF (l.EQ.ibeg) THEN
                              z3 = z(l)
                              IF (ibc3.EQ.1) pos = 1
                           ELSE IF (l.EQ.iend+1) THEN
                              z3 = z(l-1)
                              IF (ibc3.EQ.1) pos = -1
                           ELSE IF (ABS(zi(ex)-z(l)).LT.
     &                             ABS(zi(ex)-z(l-1))) THEN
                              z3 = z(l)
                              IF (ibc3.EQ.1) pos = 1
                           ELSE
                              z3 = z(l-1)
                              IF (ibc3.EQ.1) pos = -1
                           ENDIF

C C4: coef (1-alpha)*beta
                           iraycp3 = irayc
                           jraycp3 = jrayc+1
                           jraycp3 = MIN(jraycp3, njray-1)
                           indraycp3 = iraycp3 + jraycp3 * niray
                           ibeg = indir(indraycp3)
                           IF (indraycp3.EQ.nijray-1) THEN
                              iend = nz-1
                           ELSE
                              iend = indir(indraycp3+1)-1
                           ENDIF
                           np = 0
                           DO l = ibeg, iend
                              IF (z(l).GT.zi(ex)) THEN
                                 GOTO 82
                              ENDIF
                              np = np+1
                           ENDDO
                           IF ( np .EQ. 0 ) GOTO 62

 82                        ibc4 = np - (np/2)*2
                           IF (l.EQ.ibeg) THEN
                              z4 = z(l)
                              IF (ibc4.EQ.1) pos = 1
                           ELSE IF (l.EQ.iend+1) THEN
                              z4 = z(l-1)
                              IF (ibc4.EQ.1) pos = -1
                           ELSE IF (ABS(zi(ex)-z(l)).LT.
     &                             ABS(zi(ex)-z(l-1))) THEN
                              z4 = z(l)
                              IF (ibc4.EQ.1) pos = 1
                           ELSE
                              z4 = z(l-1)
                              IF (ibc4.EQ.1) pos = -1
                           ENDIF
                           
                           xrayc = xmin + hiray * irayc
                           yrayc = ymin + hjray * jrayc
                           alpha = (xi(ex) - xrayc)/hiray
                           beta  = (yi(ex) - yrayc)/hjray

C                           ibi(ex) = 
C     &                          (1.-alpha)*(1.-beta) * ibc +
C     &                          alpha*(1.-beta) * ibcp1 + 
C     &                          alpha*beta * ibcp2 +
C     &                          (1.-alpha)*beta * ibcp3 
                           
                           zint = (1.-alpha)*(1.-beta) * z1 +
     &                          alpha*(1.-beta) * z2 + 
     &                          alpha*beta * z3 +
     &                          (1.-alpha)*beta * z4

                           IF ( ibc1.EQ.1.AND.ibc2.EQ.1.AND.
     &                          ibc3.EQ.1.AND.ibc4.EQ.1) THEN
                              cellNatureField(i+j*nic+k*nicnjc) = -1
                              isMasked = 1
                              GOTO 5
                           ENDIF

                           IF ( ibc1.EQ.0.AND.ibc2.EQ.0.AND.
     &                          ibc3.EQ.0.AND.ibc4.EQ.0) THEN
C                              cellNatureField(i+j*nic+k*nicnjc) = 0
C                              isMasked = 1
                              GOTO 62
                           ENDIF

                           IF (pos.EQ.1.AND.zi(ex).LE.zint) THEN
                              cellNatureField(i+j*nic+k*nicnjc) = -1
                              isMasked = 1
                              GOTO 5
                           ENDIF

                           IF (pos.EQ.-1.AND.zi(ex).GE.zint) THEN
                              cellNatureField(i+j*nic+k*nicnjc) = -1
                              isMasked = 1
                              GOTO 5
                           ENDIF

 62                        CONTINUE
                        ENDDO
                        
C                        IF ( ib0.GE.0.5 .OR. 
C     &                       ibi(1).GE.0.5 .OR. ibi(2).GE.0.5 .OR.
C     &                       ibi(3).GE.0.5 .OR. ibi(4).GE.0.5 .OR.
C     &                       ibi(5).GE.0.5 .OR. ibi(6).GE.0.5) THEN
C                           cellNatureField(i+j*nic+k*nicnjc) 
C     &                          = -1
C                           isMasked = 1
C                           GOTO 5
C                        ENDIF 

 6                      isMasked = 1                        
                        cellNatureField(i+j*nic+k*nicnjc) = 0 
                        GOTO 5
                     ENDIF
                  ENDDO
               ENDDO
C     
 5             CONTINUE
               
            ENDDO
         ENDDO
      ENDDO
      END
C ===== XRay/MaskSearchBlankedCellsX1F.for =====
