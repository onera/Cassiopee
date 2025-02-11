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

C center_in + depth = 1 cas 2D

      SUBROUTINE k6searchblankedcellsx12d( ni, nj, nk, 
     &     meshX, meshY, 
     &     xmin,
     &     niray, njray,
     &     hiray,
     &     indir, nz, z, 
     &     cellNatureField, isMasked )
      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E ni, nj, nk   ! Dimension of the domain
      REAL_E    meshX(0:ni*nj*nk-1) ! Abcisses of mesh points
      REAL_E    meshY(0:ni*nj*nk-1) ! Ordinates of mesh points
      REAL_E    xmin            ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray           ! hi, hj of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      REAL_E    z(0:nz-1)       ! pierce points

C_OUT
      INTEGER_E cellNatureField(0:MAX(ni-1,1)*MAX(nj-1,1)*MAX(nk-1,1)-1)
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l       ! For loop in the three directions
      INTEGER_E ip, jp
      INTEGER_E ind, indray, indrayc, indraycp1
      REAL_E    dx1, dy1
      REAL_E    xp, yp
      REAL_E    xp1, yp1, xp2, yp2, xp3, yp3, xp4, yp4
      REAL_E    xmincell, ymincell
      REAL_E    xmaxcell, ymaxcell
      INTEGER_E ibeg, iend
      INTEGER_E iray, irayc, iraycp1
      INTEGER_E iraymin, iraymax
      REAL_E ibmin, ibmax, ibc, ibcp1
      REAL_E xrayc, alpha, xmax
      INTEGER_E npmin, npmax, np, ex
      INTEGER_E cellN, nij, nijray, nic, nicnjc, pos
      REAL_E    xi(4), yi(4)
      REAL_E    ibi(4), ib0, z1, z2, zi

C==============================================================================
      nij = ni * nj
      nijray = niray * njray
      nic = ni-1
      nicnjc = nic*(nj-1)
      isMasked = 0

      xmax = xmin + (niray-1) * hiray
      
      pos = 0

      DO j = 0, nj-2
         DO i = 0, nic-1
            cellN = 0
C
#include "../../Connector/Fortran/MaskCell2DF.for"            
C
C       Barycenter of cell
C              
            xp = xp1 + xp2 + xp3 + xp4
            yp = yp1 + yp2 + yp3 + yp4 
            xp = xp * ONE_FOURTH
            yp = yp * ONE_FOURTH
C
C       Interface centers
C
            xi(1) = ONE_HALF*(xp1 + xp2)
            yi(1) = ONE_HALF*(yp1 + yp2)
           
            xi(2) = ONE_HALF*(xp2 + xp3)
            yi(2) = ONE_HALF*(yp2 + yp3)          

            xi(3) = ONE_HALF*(xp3 + xp4)
            yi(3) = ONE_HALF*(yp3 + yp4)

            xi(4) = ONE_HALF*(xp1 + xp4)
            yi(4) = ONE_HALF*(yp1 + yp4)  

           
C Ray index corresponding to point cell
            IF (xmaxcell.LT.xmin) GOTO 5
            IF (xmincell.GT.xmax) GOTO 5

            iraymin = (xmincell - xmin)/hiray
            iraymax = (xmaxcell - xmin)/hiray+1
            iraymin = MAX( 0, iraymin )
            iraymax = MIN( niray-1, iraymax )   
            
            DO iray = iraymin, iraymax
             
               indray = iray
               ibeg = indir(indray)
               IF (indray.EQ.nijray-1) THEN
                  iend = nz-1
               ELSE
                  iend = indir(indray+1)-1
               ENDIF
               npmin = 0
               DO l = ibeg, iend
                  IF (z(l).GT.ymincell) THEN
                     GOTO 1
                  ENDIF
                  npmin = npmin+1
               ENDDO
 1             ibmin = npmin - (npmin/2)*2
               npmax = 0
               DO l = ibeg, iend
                  IF (z(l).GT.ymaxcell) THEN
                     GOTO 2
                  ENDIF
                  npmax = npmax+1
               ENDDO
 2             ibmax = npmax - (npmax/2)*2

               IF ( (ibmin.NE.0).OR.
     &              (ibmax.NE.0).OR.
     &              (npmin.NE.npmax) ) THEN

C      Cell intersects, we now check if center is in
                  irayc = (xp - xmin)/hiray 

                  IF (irayc.LT.0 .OR. irayc.GT.niray-1) THEN
                     ib0 = 0.
                     GOTO 61
                  ENDIF
                  indrayc = irayc
                  ibeg = indir(indrayc)
                  IF (indrayc.EQ.nijray-1) THEN
                     iend = nz-1
                  ELSE
                     iend = indir(indrayc+1)-1
                  ENDIF
                  np = 0
                  DO l = ibeg, iend
                     IF (z(l).GT.yp) THEN
                        GOTO 3
                     ENDIF
                     np = np+1
                  ENDDO
                  IF ( np.EQ.0) GOTO 61

 3                ibc = np - (np/2)*2
                  IF (l.EQ.ibeg) THEN
                     z1 = z(l)
                     IF (ibc.EQ.1) pos = 1
                  ELSE IF (l.EQ.iend+1) THEN
                     z1 = z(l-1)
                     IF (ibc.EQ.1) pos = -1
                  ELSE IF (ABS(yp-z(l)).LT.ABS(yp-z(l-1))) THEN
                     z1 = z(l)
                     IF (ibc.EQ.1) pos = 1
                  ELSE
                     z1 = z(l-1)
                     IF (ibc.EQ.1) pos = -1
                  ENDIF

                  iraycp1 = irayc+1
                  iraycp1 = MIN(iraycp1, niray-1)
                  indraycp1 = iraycp1
                  ibeg = indir(indraycp1)
                  IF (indraycp1.EQ.nijray-1) THEN
                     iend = nz-1
                  ELSE
                     iend = indir(indraycp1+1)-1
                  ENDIF
                  np = 0
                  DO l = ibeg, iend
                     IF (z(l).GT.yp) THEN
                        GOTO 4
                     ENDIF
                     np = np+1
                  ENDDO
                  IF ( np.EQ.0) GOTO 61

 4                ibcp1 = np - (np/2)*2
                  IF (l.EQ.ibeg) THEN
                     z2 = z(l)
                     IF (ibcp1.EQ.1) pos = 1
                  ELSE IF (l.EQ.iend+1) THEN
                     z2 = z(l-1)
                     IF (ibcp1.EQ.1) pos = -1
                  ELSE IF (ABS(yp-z(l)).LT.ABS(yp-z(l-1))) THEN
                     z2 = z(l)
                     IF (ibcp1.EQ.1) pos = 1
                  ELSE
                     z2 = z(l-1)
                     IF (ibcp1.EQ.1) pos = -1
                  ENDIF

                  xrayc = xmin + hiray * irayc
                  alpha = (xp - xrayc)/hiray

C                  ib0 = (1.-alpha) * ibc + alpha * ibcp1
C                  IF (ib0.GE.0.5) THEN
C                     cellNatureField(i+j*nic) = -1
C                     isMasked = 1
C                     GOTO 5
C                  ENDIF

                  zi = (1.-alpha) * z1 + alpha * z2
                  IF (ibc.EQ.1 .AND. ibcp1.EQ.1) THEN
                     cellNatureField(i+j*nic) = -1
                     isMasked = 1
                     GOTO 5
                  ENDIF
                  IF (ibc.EQ.0 .AND. ibcp1.EQ.0) THEN
c                     cellNatureField(i+j*nic) = 0
c                     isMasked = 1
                     GOTO 61
                  ENDIF
 
                  IF (pos.EQ.1 .AND. yp.LE. zi) THEN
                     cellNatureField(i+j*nic) = -1
                     isMasked = 1
                     GOTO 5
                  ENDIF
 
                  IF (pos.EQ.-1 .AND. yp.GE.zi) THEN
                     cellNatureField(i+j*nic) = -1
                     isMasked = 1
                     GOTO 5
                  ENDIF

 61               CONTINUE
C         Also, check if EX points are in
                  DO ex = 1, 4
                     irayc = (xi(ex) - xmin) / hiray
                     IF (irayc.LT.0 .OR. irayc.GT.niray-1) THEN
                        ibi(ex) = 0.
                        GOTO 62
                     ENDIF
                     indrayc = irayc
                     ibeg = indir(indrayc)
                     IF (indrayc.EQ.nijray-1) THEN
                        iend = nz-1
                     ELSE
                        iend = indir(indrayc+1)-1
                     ENDIF
                     np = 0
                     DO l = ibeg, iend
                        IF (z(l).GT.yi(ex)) THEN
                           GOTO 7
                        ENDIF
                        np = np+1
                     ENDDO
                     IF ( np.EQ.0) GOTO 62
                  
 7                   ibc = np - (np/2)*2
                     IF (l.EQ.ibeg) THEN
                        z1 = z(l)
                        IF (ibc.EQ.1) pos = 1
                     ELSE IF (l.EQ.iend+1) THEN
                        z1 = z(l-1)
                        IF (ibc.EQ.1) pos = -1
                     ELSE IF (ABS(yi(ex)-z(l)).LT.
     &                       ABS(yi(ex)-z(l-1))) THEN
                        z1 = z(l)
                        IF (ibc.EQ.1) pos = 1
                     ELSE
                        z1 = z(l-1)
                        IF (ibc.EQ.1) pos = -1
                     ENDIF

                     iraycp1 = irayc+1
                     iraycp1 = MIN(iraycp1, niray-1)
                     indraycp1 = iraycp1
                     ibeg = indir(indraycp1)
                     IF (indraycp1.EQ.nijray-1) THEN
                        iend = nz-1
                     ELSE
                        iend = indir(indraycp1+1)-1
                     ENDIF
                     np = 0
                     DO l = ibeg, iend
                        IF (z(l).GT.yi(ex)) THEN
                           GOTO 8
                        ENDIF
                        np = np+1
                     ENDDO
                     IF ( np.EQ.0) GOTO 62

 8                   ibcp1 = np - (np/2)*2
                     IF (l.EQ.ibeg) THEN
                        z2 = z(l)
                        IF (ibcp1.EQ.1) pos = 1
                     ELSE IF (l.EQ.iend+1) THEN
                        z2 = z(l-1)
                        IF (ibcp1.EQ.1) pos = -1
                     ELSE IF (ABS(yi(ex)-z(l)).LT.
     &                       ABS(yi(ex)-z(l-1))) THEN
                        z2 = z(l)
                        IF (ibcp1.EQ.1) pos = 1
                     ELSE
                        z2 = z(l-1)
                        IF (ibcp1.EQ.1) pos = -1
                     ENDIF

                     xrayc = xmin + hiray * irayc
                     alpha = (xi(ex) - xrayc)/hiray

C                     ibi(ex) = (1.-alpha) * ibc + alpha * ibcp1
C                     IF (ibi(ex).GE.0.5) THEN
C                        cellNatureField(i+j*nic) = -1
C                        isMasked = 1
C                        GOTO 5
C                     ENDIF
                  
                     zi = (1.-alpha) * z1 + alpha * z2
                     IF (ibc.EQ.1 .AND. ibcp1.EQ.1) THEN
                        cellNatureField(i+j*nic) = -1
                        isMasked = 1
                        GOTO 5
                     ENDIF
                     IF (ibc.EQ.0 .AND. ibcp1.EQ.0) THEN
C                        cellNatureField(i+j*nic) = 0
C                        isMasked = 1
                        GOTO 62
                     ENDIF
                     
                     IF (pos.EQ.1 .AND. yi(ex).LE. zi) THEN
                        cellNatureField(i+j*nic) = -1
                        isMasked = 1
                        GOTO 5
                     ENDIF
                  
                     IF (pos.EQ.-1 .AND. yi(ex).GE.zi) THEN
                        cellNatureField(i+j*nic) = -1
                        isMasked = 1
                        GOTO 5
                     ENDIF

 62               CONTINUE
               ENDDO

 6             cellNatureField(i+j*nic) = 0
               isMasked = 1
               GOTO 5

            ENDIF
         ENDDO
         
 5       CONTINUE
         
      ENDDO
      ENDDO
      END




