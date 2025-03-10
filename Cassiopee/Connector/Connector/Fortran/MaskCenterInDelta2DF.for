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
                  indrayc = irayc
                  ibeg = indir(indrayc)
                  IF (indrayc.EQ.nijray-1) THEN
                     iend = nz-1
                  ELSE
                     iend = indir(indray+1)-1
                  ENDIF
                  nnp = 0
                  DO l = ibeg, iend
                     IF (z(l).GT.yp) THEN
                        GOTO 3
                     ENDIF
                     nnp = nnp+1
                  ENDDO
                  IF ( nnp .EQ. 0 ) GOTO 6

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
                  IF ( nnp .EQ. 0 ) GOTO 6

 4                ibcp1 = nnp - (nnp/2)*2
                  xrayc = xmin + hiray * irayc
                  ntemp = (xp - xrayc)/hiray
                  alpha = (xp-ntemp*hiray-xrayc)/hiray
                  ib = (1.-alpha) * ibc + alpha * ibcp1
                  
                  IF ( ib.GE.0.5 ) THEN
                     cellNatureField(et) = -1
                     isMasked = 1
C                     tmp = tmp+1
                     GOTO 5
                  ENDIF 
                  
 6                IF (cellNatureField(et).NE.-1) THEN
c                     cellNatureField(et) = 0
C                     isMasked = 1
                  END IF
                  
C     WRITE(*,*) 'cell ',i,j,ymincell, ymaxcell,
C     &                 ' is blanked:',npmin,
C     &                 npmax, ', ray=',indray, ' - ',ibeg, iend,
C     &                 ' = ',l
C                  GOTO 5
               ENDDO
            ENDIF
C     
         ENDDO
 5       CONTINUE
