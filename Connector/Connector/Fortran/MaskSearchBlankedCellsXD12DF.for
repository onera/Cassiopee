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
C      A cell is blanked if the bounding box of cell (+/- delta) 
C      intersects the XRay mask (+/- delta)
C      and the center of cell (+/- delta) is in the mask (+/- delta). 
C      Depth = 1 case.
C      For 2D: the mesh must be nk=2 and in plane (x,y).

      SUBROUTINE k6searchblankedcellsxd12d( ni, nj, nk, 
     &     meshX, meshY, 
     &     xmin, 
     &     niray, njray,
     &     hiray, 
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
      REAL_E    xmin            ! Ray min
      INTEGER_E niray, njray    ! Dim of ray plane
      REAL_E    hiray           ! hi of ray plane
      INTEGER_E indir(0:niray*njray-1)
      INTEGER_E nz              ! Number of Z pierce points
      INTEGER_E np              ! Number of interpolated points
      INTEGER_E listOfInterpolatedPoints(0:np-1)
      REAL_E    delta           ! Wall distance
      REAL_E    z(0:nz-1)       ! pierce points

C_OUT
      INTEGER_E cellNatureField(0:(ni-1)*(nj-1)*(nk-1)-1) 
      INTEGER_E isMasked
C_LOCAL
      INTEGER_E i, j, k, l       ! For loop in the three directions
      INTEGER_E ip, jp
      INTEGER_E ind, inds, indray, indrayc, indraycp1
      INTEGER_E n
      REAL_E    dx1, dy1
      REAL_E    xp, yp
      REAL_E    xp1, yp1, xp2, yp2, xp3, yp3, xp4, yp4
      REAL_E    xmincell, ymincell
      REAL_E    xmaxcell, ymaxcell
      INTEGER_E ibeg, iend
      INTEGER_E iray, irayc, iraycp1, iraycex
      INTEGER_E iraymin, iraymax
      INTEGER_E irayminp, iraymaxp
      REAL_E ibmin, ibmax, ibc, ibcp1
      REAL_E xrayc, alpha, xmax
      INTEGER_E npmin, npmax, nnp, ex
      INTEGER_E cellN, nij, nijray, nic, nicnjc
      REAL_E    xi(4), yi(4)
      REAL_E    ibi(4), ib0
      INTEGER_E ntemp

C
      nij = ni * nj
      nijray = niray * njray
      nic = ni-1
      nicnjc = nic*(nj-1)
      isMasked = 0

      xmax = xmin + (niray-1) * hiray
      
      DO n = 0, np-1
         cellN = 0
         inds = listOfInterpolatedPoints(n)
         j = inds/nic
         i = inds-j*nic
C     
#include "../../Connector/Fortran/MaskCell2DF.for"
C     
C     Barycenter of cell
C     
         xp = xp1 + xp2 + xp3 + xp4
         yp = yp1 + yp2 + yp3 + yp4    
         xp = xp * ONE_FOURTH
         yp = yp * ONE_FOURTH
C     
C     Interface centers
C     
         xi(1) = ONE_HALF*(xp1 + xp2)
         yi(1) = ONE_HALF*(yp1 + yp2)
         
         xi(2) = ONE_HALF*(xp2 + xp3)
         yi(2) = ONE_HALF*(yp2 + yp3)          
         
         xi(3) = ONE_HALF*(xp3 + xp4)
         yi(3) = ONE_HALF*(yp3 + yp4)
         
         xi(4) = ONE_HALF*(xp1 + xp4)
         yi(4) = ONE_HALF*(yp1 + yp4)
         
         xmincell = xmincell-delta
         xmaxcell = xmaxcell+delta
         
C     Ray index corresponding to point cell
         IF (xmaxcell.LT.(xmin-delta)) GOTO 5
         IF (xmincell.GT.(xmax+delta)) GOTO 5
         
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
 1          ibmin = npmin - (npmin/2)*2
            npmax = 0
            DO l = ibeg, iend
               IF (z(l).GT.ymaxcell) THEN
                  GOTO 2
               ENDIF
               npmax = npmax+1
            ENDDO
 2          ibmax = npmax - (npmax/2)*2
            
            IF ( (ibmin.NE.0).OR.
     &           (ibmax.NE.0).OR.
     &           (npmin.NE.npmax) ) THEN
               
C     Cell intersects, we now check if center is in
               irayminp = (xp - delta - xmin)/hiray
               iraymaxp = (xp + delta - xmin)/hiray+1
               irayminp = MAX( 0, irayminp )
               iraymaxp = MIN( niray-1, iraymaxp )   
               
               DO irayc = irayminp, iraymaxp
c$$$  IF (irayc.LT.0 .OR. irayc.GT.niray-1) THEN
c$$$  ib0 = 0. 
c$$$  GOTO 61
c$$$  ENDIF
                  indrayc = irayc
                  ibeg = indir(indrayc)
                  IF (indrayc.EQ.nijray-1) THEN
                     iend = nz-1
                  ELSE
                     iend = indir(indrayc+1)-1
                  ENDIF
                  nnp = 0
                  DO l = ibeg, iend
                     IF (z(l).GT.yp) THEN
                        GOTO 3
                     ENDIF
                     nnp = nnp+1
                  ENDDO
                  IF ( nnp.EQ.0 ) GOTO 61

 3                ibc = nnp - (nnp/2)*2
                  iraycp1 = irayc+1
                  iraycp1 = MIN(iraycp1, niray-1)
                  indraycp1 = iraycp1
                  ibeg = indir(indraycp1)
                  IF (indraycp1.EQ.nijray-1) THEN
                     iend = nz-1
                  ELSE
                     iend = indir(indraycp1+1)-1
                  ENDIF
                  nnp = 0
                  DO l = ibeg, iend
                     IF (z(l).GT.yp) THEN
                        GOTO 4
                     ENDIF
                     nnp = nnp+1
                  ENDDO
                  IF ( nnp.EQ.0 ) GOTO 61

 4                ibcp1 = nnp - (nnp/2)*2
                  xrayc = xmin + hiray * irayc
                  ntemp = (xp - xrayc)/hiray
                  alpha = (xp-ntemp*hiray-xrayc)/hiray
                  ib0 = (1.-alpha) * ibc + alpha * ibcp1
                  
 61               CONTINUE
C     Check if EX points are in
                  DO ex = 1, 4
                     irayminp = (xi(ex) - delta - xmin)/hiray
                     iraymaxp = (xi(ex) + delta - xmin)/hiray+1
                     irayminp = MAX( 0, irayminp )
                     iraymaxp = MIN( niray-1, iraymaxp )   
                     DO iraycex = irayminp, iraymaxp
c$$$                        IF (iraycex.LT.0 .OR. iraycex.GT.niray-1) THEN
c$$$                           ibi(ex) = 0.
c$$$                           GOTO 62
c$$$                        ENDIF
                        indrayc = iraycex
                        ibeg = indir(indrayc)
                        IF (indrayc.EQ.nijray-1) THEN
                           iend = nz-1
                        ELSE
                           iend = indir(indrayc+1)-1
                        ENDIF
                        nnp = 0
                        DO l = ibeg, iend
                           IF (z(l).GT.yi(ex)) THEN
                              GOTO 7
                           ENDIF
                           nnp = nnp+1
                        ENDDO
                        IF ( nnp.EQ.0 ) GOTO 62

 7                      ibc = nnp - (nnp/2)*2
                        iraycp1 = iraycex+1
                        iraycp1 = MIN(iraycp1, niray-1)
                        indraycp1 = iraycp1
                        ibeg = indir(indraycp1)
                        IF (indraycp1.EQ.nijray-1) THEN
                           iend = nz-1
                        ELSE
                           iend = indir(indraycp1+1)-1
                        ENDIF
                        nnp = 0
                        DO l = ibeg, iend
                           IF (z(l).GT.yi(ex)) THEN
                              GOTO 8
                           ENDIF
                           nnp = nnp+1
                        ENDDO
                        IF ( nnp.EQ.0 ) GOTO 62

 8                      ibcp1 = nnp - (nnp/2)*2
                        xrayc = xmin + hiray * iraycex
                        ntemp = (xi(ex) - xrayc)/hiray
                        alpha = (xi(ex)-ntemp*hiray-xrayc)/hiray
                        ibi(ex) = (1.-alpha) * ibc + alpha * ibcp1
                     ENDDO
 62                  CONTINUE
                  ENDDO
                  
                  IF ( ib0.GE.0.5 .OR. ibi(1).GE.0.5 .OR.
     &                 ibi(2).GE.0.5 .OR. ibi(3).GE.0.5 .OR.
     &                 ibi(4).GE.0.5 ) THEN
                     cellNatureField(inds) = -1
                     isMasked = 1
                     GOTO 5
                  ENDIF 
                  
 6                IF (cellNatureField(inds).NE.-1) THEN
                     cellNatureField(inds) = 0
C     isMasked = 1
                  ENDIF
                  
C                  GOTO 5
               ENDDO
            ENDIF
         ENDDO
C     
 5       CONTINUE
         
      ENDDO
      END




