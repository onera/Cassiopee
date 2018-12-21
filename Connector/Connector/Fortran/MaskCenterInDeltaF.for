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

C     Ray index corresponding to point cell

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
c$$$  IF (irayc.LT.0) GOTO 6
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
                           nnp = np+1
                        ENDDO
                        IF (nnp.EQ.0) GOTO 6

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
                        IF (nnp.EQ.0) GOTO 6

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
                        IF (nnp.EQ.0) GOTO 6

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
                        IF (nnp.EQ.0) GOTO 6

 44                     ibcp3 = nnp - (nnp/2)*2
                        
                        xrayc = xmin + hiray * irayc
                        yrayc = ymin + hjray * jrayc
                        ntemp1 = (xp - xrayc)/hiray
                        ntemp2 = (yp - yrayc)/hjray
                        alpha = (xp-ntemp1*hiray-xrayc)/hiray
                        beta  = (yp-ntemp2*hjray-yrayc)/hjray
                        
                        ib = (1.-alpha)*(1.-beta) * ibc +
     &                       alpha*(1.-beta) * ibcp1 + 
     &                       alpha*beta * ibcp2 +
     &                       (1.-alpha)*beta * ibcp3

                        IF ( ib.GE.0.5 ) THEN
                           cellNatureField(et) = -1
                           isMasked = 1
C                           tmp = tmp+1
                           GOTO 5
                        ENDIF 
 6                      IF (cellNatureField(et).NE.-1) THEN
c                           cellNatureField(et) = 0
c                           isMasked = 1
                        ENDIF
                        
C                        GOTO 5
                     ENDDO
                  ENDDO
               ENDIF
C     
            ENDDO
         ENDDO
 5       CONTINUE
