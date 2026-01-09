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

C=============================================================================
C Common part for center_in algorithms
C=============================================================================
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

               IF ((ibmin.NE.0).OR.
     &             (ibmax.NE.0).OR.
     &             (npmin.NE.npmax)) THEN

C      Cell intersects, we now check if center is in
                  irayc = (xp - xmin)/hiray
 
                  IF (irayc.LT.0) GOTO 6
                  IF (irayc.GT.niray-1) GOTO 6
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
                  IF ( nnp.EQ.0) GOTO 6
 3                ibc = nnp - (nnp/2)*2
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
                  nnp = 0
                  DO l = ibeg, iend
                     IF (z(l).GT.yp) THEN
                        GOTO 4
                     ENDIF
                     nnp = nnp+1
                  ENDDO
                  IF ( nnp.EQ.0) GOTO 6

 4                ibcp1 = nnp - (nnp/2)*2
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

C                  ib = (1.-alpha) * ibc + alpha * ibcp1
C                  IF ( ib.GE.0.5 ) THEN
C                     cellNatureField(et) = -1
C                     isMasked = 1
C                     GOTO 5
C                  ENDIF 

                  zi = (1.-alpha) * z1 + alpha * z2
                  IF (ibc.EQ.1 .AND. ibcp1.EQ.1) THEN
                     cellNatureField(et) = -1
                     isMasked = 1
                     GOTO 5
                  ENDIF
                  IF (ibc.EQ.0 .AND. ibcp1.EQ.0) THEN
c                     cellNatureField(et) = 0
                     isMasked = 1
                     GOTO 5
                  ENDIF
 
                  IF (pos.EQ.1 .AND. yp.LE. zi) THEN
                     cellNatureField(et) = -1
                     isMasked = 1
                     GOTO 5
                  ENDIF
 
                  IF (pos.EQ.-1 .AND. yp.GE.zi) THEN
                     cellNatureField(et) = -1
                     isMasked = 1
                     GOTO 5
                  ENDIF

 6                isMasked = 1
c                  cellNatureField(et) = 0
                  GOTO 5

               ENDIF
C     
            ENDDO
 5          CONTINUE
