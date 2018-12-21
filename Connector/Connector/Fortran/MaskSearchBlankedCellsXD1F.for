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

C Search the cells which will be blanked "within" a distance of XRay mask.     
C      Put the results in a cell nature field array.
C      A cell is blanked if the bounding box of cell intersects (+/- delta) 
C      the XRay mask (+/- delta)
C      and the center of cell (+/- delta) is in the mask (+/- delta). 
C      Depth = 1 case.

      SUBROUTINE k6searchblankedcellsxd1( ni, nj, nk, 
     &     meshX, meshY, meshZ,
     &     xmin, ymin,
     &     niray, njray,
     &     hiray, hjray,
     &     indir, nz, z, 
     &     listOfInterpolatedPoints, np,
     &     delta,
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
      INTEGER_E np              ! Number of interpolated points
      INTEGER_E listOfInterpolatedPoints(0:np-1)
      REAL_E    delta           ! Wall distance

C_OUT
      INTEGER_E cellNatureField(0:(ni-1)*(nj-1)*(nk-1)-1) 
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l       ! For loop in the three directions
      INTEGER_E ip, jp, kp
      INTEGER_E ind, inds, indray
      INTEGER_E indrayc, indraycp1, indraycp2, indraycp3
      REAL_E    dx1, dy1, dz1
      REAL_E    xp, yp, zp
      REAL_E    xp1, xp2, xp3, xp4, xp5, xp6, xp7, xp8
      REAL_E    yp1, yp2, yp3, yp4, yp5, yp6, yp7, yp8
      REAL_E    zp1, zp2, zp3, zp4, zp5, zp6, zp7, zp8
      REAL_E    xmincell, ymincell, zmincell
      REAL_E    xmaxcell, ymaxcell, zmaxcell
      INTEGER_E ibeg, iend
      INTEGER_E iray, irayc, iraycex, iraycp1, iraycp2, iraycp3
      INTEGER_E jray, jrayc, jraycex, jraycp1, jraycp2, jraycp3
      INTEGER_E iraymin, iraymax, jraymin, jraymax
      INTEGER_E irayminp, iraymaxp, jrayminp, jraymaxp
      REAL_E ibmin, ibmax, ibc, ibcp1, ibcp2, ibcp3
      REAL_E xrayc, yrayc, alpha, beta, xmax, ymax
      INTEGER_E npmin, npmax, nnp, ex
      INTEGER_E cellN, nij, nijray, nic, nicnjc
      REAL_E    xi(6), yi(6), zi(6)
      REAL_E    ibi(6), ib0
      INTEGER_E n, tmp, ntemp1, ntemp2
C
      tmp = 0
      nij = ni * nj
      nijray = niray * njray
      nic = ni-1
      nicnjc = nic*(nj-1)
      isMasked = 0

      xmax = xmin + (niray-1) * hiray
      ymax = ymin + (njray-1) * hjray 

      DO n = 0, np-1
         cellN = 0
C     
         inds = listOfInterpolatedPoints(n)
         k = inds/nicnjc
         j = (inds-k*nicnjc)/nic
         i = inds-k*nicnjc-j*nic
         
C
#include "../../Connector/Fortran/MaskCell3DF.for"
C     
C     Barycenter of cell
C     
         xp = xp1 + xp2 + xp3 + xp4 + xp5 + xp6 + xp7 + xp8
         yp = yp1 + yp2 + yp3 + yp4 + yp5 + yp6 + yp7 + yp8 
         zp = zp1 + zp2 + zp3 + zp4 + zp5 + zp6 + zp7 + zp8
         xp = xp * ONE_EIGHT
         yp = yp * ONE_EIGHT
         zp = zp * ONE_EIGHT
C     
C  
C centers of the interfaces
C--------------------------               
C     interface  2376
         xi(1) = ONE_FOURTH * ( xp2 + xp3 + xp6 + xp7 ) 
         yi(1) = ONE_FOURTH * ( yp2 + yp3 + yp6 + yp7 )
         zi(1) = ONE_FOURTH * ( zp2 + zp3 + zp6 + zp7 ) 
         
C     interface 1485
         xi(2) = ONE_FOURTH * ( xp1 + xp4 + xp5 + xp8 ) 
         yi(2) = ONE_FOURTH * ( yp1 + yp4 + yp5 + yp8 ) 
         zi(2) = ONE_FOURTH * ( zp1 + zp4 + zp5 + zp8 )           
         
C     interface 1265
         xi(3) = ONE_FOURTH * ( xp1 + xp2 + xp6 + xp5 )
         yi(3) = ONE_FOURTH * ( yp1 + yp2 + yp6 + yp5 )  
         zi(3) = ONE_FOURTH * ( zp1 + zp2 + zp6 + zp5 )           
         
C     interface 4378
         xi(4) = ONE_FOURTH * ( xp3 + xp4 + xp7 + xp8 )
         yi(4) = ONE_FOURTH * ( yp3 + yp4 + yp7 + yp8 )  
         zi(4) = ONE_FOURTH * ( zp3 + zp4 + zp7 + zp8 )
         
C     interface 1234
         xi(5) = ONE_FOURTH * ( xp1 + xp2 + xp3 + xp4 ) 
         yi(5) = ONE_FOURTH * ( yp1 + yp2 + yp3 + yp4 )
         zi(5) = ONE_FOURTH * ( zp1 + zp2 + zp3 + zp4 )
         
C     interface  5678
         xi(6) = ONE_FOURTH * ( xp5 + xp6 + xp7 + xp8 ) 
         yi(6) = ONE_FOURTH * ( yp5 + yp6 + yp7 + yp8 )
         zi(6) = ONE_FOURTH * ( zp5 + zp6 + zp7 + zp8 )   
C     
C Ray index corresponding to point cell
               
         xmincell = xmincell-delta
         ymincell = ymincell-delta
         xmaxcell = xmaxcell+delta
         ymaxcell = ymaxcell+delta

         IF (xmaxcell.LT.(xmin-delta)) GOTO 5
         IF (xmincell.GT.(xmax+delta)) GOTO 5
         IF (ymaxcell.LT.(ymin-delta)) GOTO 5
         IF (ymincell.GT.(ymax+delta)) GOTO 5

         iraymin = (xmincell - xmin)/hiray
         iraymax = (xmaxcell - xmin)/hiray+1
         jraymin = (ymincell - ymin)/hjray
         jraymax = (ymaxcell - ymin)/hjray+1
         iraymin = MAX( 0, iraymin )
         iraymax = MIN( niray-1, iraymax )   
         jraymin = MAX( 0, jraymin )
         jraymax = MIN( njray-1, jraymax )
               
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
 1             ibmin = npmin - (npmin/2)*2
               npmax = 0
               DO l = ibeg, iend
                  IF (z(l).GT.zmaxcell) THEN
                     GOTO 2
                  ENDIF
                  npmax = npmax+1
               ENDDO
 2             ibmax = npmax - (npmax/2)*2
               
               IF ( (ibmin.NE.0).OR.
     &              (ibmax.NE.0).OR.
     &              (npmin.NE.npmax) ) THEN
                        
C     Cell intersects, we now check if center is in
C     point B1: coeff (1-alpha)*(1-beta)
                  irayminp = (xp - delta - xmin)/hiray
                  iraymaxp = (xp + delta - xmin)/hiray+1
                  irayminp = MAX( 0, irayminp )
                  iraymaxp = MIN( niray-1, iraymaxp ) 
                  
                  jrayminp = (yp - delta - ymin)/hjray
                  jraymaxp = (yp + delta - ymin)/hjray+1
                  jrayminp = MAX( 0, jrayminp )
                  jraymaxp = MIN( njray-1, jraymaxp )
                  DO jrayc = jrayminp, jraymaxp
                     DO irayc = irayminp, iraymaxp
C     c$$$                        IF (irayc.LT.0 .OR. irayc.GT.niray-1 .OR.
c$$$  &                      jrayc.LT.0 .OR. jrayc.GT.njray-1) THEN
c$$$  ib0 = 0.
c$$$  GOTO 61 
c$$$  ENDIF
c$$$  IF (jrayc.LT.0) GOTO 6
c$$$  IF (irayc.GT.niray-1) GOTO 6
c$$$  IF (jrayc.GT.njray-1) GOTO 6
                        indrayc = irayc + jrayc * niray
                        ibeg = indir(indrayc)
                        IF (indrayc.EQ.nijray-1) THEN
                           iend = nz-1
                        ELSE
                           iend = indir(indrayc+1)-1
                        ENDIF
                        nnp = 0
                        DO l = ibeg, iend
                           IF (z(l).GT.zp) THEN
                              GOTO 3
                           ENDIF
                           nnp = nnp+1
                        ENDDO
                        IF ( nnp .EQ. 0 ) GOTO 61

 3                      ibc = nnp - (nnp/2)*2

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
                        nnp = 0
                        DO l = ibeg, iend
                           IF (z(l).GT.zp) THEN
                              GOTO 4
                           ENDIF
                           nnp = nnp+1
                        ENDDO
                        IF ( nnp .EQ. 0 ) GOTO 61

 4                      ibcp1 = nnp - (nnp/2)*2

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
                        nnp = 0
                        DO l = ibeg, iend
                           IF (z(l).GT.zp) THEN
                              GOTO 34
                           ENDIF
                           nnp = nnp+1
                        ENDDO
                        IF ( nnp .EQ. 0 ) GOTO 61

 34                     ibcp2 = nnp - (nnp/2)*2

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
                        nnp = 0
                        DO l = ibeg, iend
                           IF (z(l).GT.zp) THEN
                              GOTO 44
                           ENDIF
                           nnp = nnp+1
                        ENDDO
                        IF ( nnp .EQ. 0 ) GOTO 61

 44                     ibcp3 = nnp - (nnp/2)*2
                        
                        xrayc = xmin + hiray * irayc
                        yrayc = ymin + hjray * jrayc
                        ntemp1 = (xp - xrayc)/hiray
                        ntemp2 = (yp - yrayc)/hjray
                        alpha = (xp-ntemp1*hiray-xrayc)/hiray
                        beta  = (yp-ntemp2*hjray-yrayc)/hjray
                       
                        ib0 =(1.-alpha)*(1.-beta) * ibc +
     &                       alpha*(1.-beta) * ibcp1 + 
     &                       alpha*beta * ibcp2 +
     &                       (1.-alpha)*beta * ibcp3

 61                     CONTINUE
C     Check if EX points are in
                        DO ex = 1, 6
                           irayminp = (xi(ex) - delta - xmin)/hiray
                           iraymaxp = (xi(ex) + delta - xmin)/hiray+1
                           irayminp = MAX( 0, irayminp )
                           iraymaxp = MIN( niray-1, iraymaxp ) 
                           
                           jrayminp = (yi(ex) - delta - ymin)/hjray
                           jraymaxp = (yi(ex) + delta - ymin)/hjray+1
                           jrayminp = MAX( 0, jrayminp )
                           jraymaxp = MIN( njray-1, jraymaxp )
                           DO jraycex = jrayminp, jraymaxp
                              DO iraycex = irayminp, iraymaxp
c$$$  
c$$$  IF (irayc.LT.0 .OR. irayc.GT.niray-1.OR.
c$$$  &                         jrayc.LT.0 .OR. jrayc.GT.njray-1) THEN 
c$$$  ibi(ex) = 0.
c$$$  GOTO 62
c$$$  ENDIF
                                 
C     C1: ex point of coef (1-alpha)*(1-beta)
                                 indrayc = iraycex + jraycex * niray
                                 ibeg = indir(indrayc)
                                 IF (indrayc.EQ.nijray-1) THEN
                                    iend = nz-1
                                 ELSE
                                    iend = indir(indrayc+1)-1
                                 ENDIF
                                 nnp = 0
                                 DO l = ibeg, iend
                                    IF (z(l).GT.zi(ex)) THEN
                                       GOTO 7
                                    ENDIF
                                    nnp = nnp+1
                                 ENDDO
                                 IF ( nnp .EQ. 0 ) GOTO 62

 7                               ibc = nnp - (nnp/2)*2
                                 
C     C2: ex point of coef alpha*(1-beta)
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
                                 nnp = 0
                                 DO l = ibeg, iend
                                    IF (z(l).GT.zi(ex)) THEN
                                       GOTO 8
                                    ENDIF
                                    nnp = nnp+1
                                 ENDDO
                                 IF ( nnp .EQ. 0 ) GOTO 62

 8                               ibcp1 = nnp - (nnp/2)*2

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
                                 nnp = 0
                                 DO l = ibeg, iend
                                    IF (z(l).GT.zi(ex)) THEN
                                       GOTO 72
                                    ENDIF
                                    nnp = nnp+1
                                 ENDDO
                                 IF ( nnp .EQ. 0 ) GOTO 62

 72                              ibcp2 = nnp - (nnp/2)*2

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
                                 nnp = 0
                                 DO l = ibeg, iend
                                    IF (z(l).GT.zi(ex)) THEN
                                       GOTO 82
                                    ENDIF
                                    nnp = nnp+1
                                 ENDDO
                                 IF ( nnp .EQ. 0 ) GOTO 62
 82                              ibcp3 = nnp - (nnp/2)*2
                                 
                                 xrayc = xmin + hiray * iraycex
                                 yrayc = ymin + hjray * jraycex
                                 ntemp1 = (xi(ex) - xrayc)/hiray
                                 ntemp2 = (yi(ex) - yrayc)/hjray
                                 alpha = (xi(ex)
     &                                -ntemp1*hiray-xrayc)/hiray
                                 beta = (yi(ex)
     &                                -ntemp2*hjray-yrayc)/hjray
                                 ibi(ex) = 
     &                                (1.-alpha)*(1.-beta) * ibc +
     &                                alpha*(1.-beta) * ibcp1 + 
     &                                alpha*beta * ibcp2 +
     &                                (1.-alpha)*beta * ibcp3
                              ENDDO
                           ENDDO 
 62                        CONTINUE
                        ENDDO
                        
                        IF ( ib0.GE.0.5 .OR. 
     &                       ibi(1).GE.0.5 .OR. ibi(2).GE.0.5 .OR.
     &                       ibi(3).GE.0.5 .OR. ibi(4).GE.0.5 .OR.
     &                       ibi(5).GE.0.5.OR.ibi(6).GE.0.5) THEN
                           cellNatureField(inds) = -1
                           isMasked = 1
                           GOTO 5
                        ENDIF 
 6                      IF (cellNatureField(inds).NE.-1) THEN
                           cellNatureField(inds) = 0                  
C     isMasked = 1
                        ENDIF
                        
C     GOTO 5
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
C     
 5       CONTINUE
               
      ENDDO
      END
