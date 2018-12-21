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
               IF (xmaxcell.LT.xmin) GOTO 5
               IF (xmincell.GT.xmax) GOTO 5
               IF (ymaxcell.LT.ymin) GOTO 5
               IF (ymincell.GT.ymax) GOTO 5
               
               iraymin = (xmincell-xmin)/hiray
               iraymax = (xmaxcell-xmin)/hiray+1
               jraymin = (ymincell-ymin)/hjray
               jraymax = (ymaxcell-ymin)/hjray+1
               iraymin = MAX(0, iraymin)
               iraymax = MIN(niray-1, iraymax)
               jraymin = MAX(0, jraymin)
               jraymax = MIN(njray-1, jraymax)
 
               DO jray = jraymin, jraymax
                  DO iray = iraymin, iraymax
                     
                     indray = iray + jray*niray
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
                        irayc = (xp-xmin)/hiray
                        jrayc = (yp-ymin)/hjray
                        IF (irayc.LT.0) GOTO 6
                        IF (jrayc.LT.0) GOTO 6
                        IF (irayc.GT.niray-1) GOTO 6
                        IF (jrayc.GT.njray-1) GOTO 6
                        indrayc = irayc + jrayc*niray
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
                        IF (np.EQ.0) GOTO 6
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
                        IF (np.EQ.0) GOTO 6
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
                        IF (np.EQ.0) GOTO 6
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
                        IF (np.EQ.0) GOTO 6
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
                        alpha = (xp-xrayc)/hiray
                        beta  = (yp-yrayc)/hjray
                        
c$$$                        ib = (1.-alpha)*(1.-beta) * ibc1 +
c$$$     &                       alpha*(1.-beta) * ibc2 + 
c$$$     &                       alpha*beta * ibc3 +
c$$$     &                       (1.-alpha)*beta * ibc4                       
c$$$                        IF ( ib.GE.0.5 ) THEN
c$$$                           cellNatureField(et )
c$$$     &                          = -1
c$$$                           isMasked = 1
c$$$                           GOTO 5
c$$$                        ENDIF 

                        zi = (1.-alpha)*(1.-beta)*z1 +
     &                       alpha*(1.-beta)*z2 + 
     &                       alpha*beta*z3 +
     &                       (1.-alpha)*beta*z4

                        IF (ibc1.EQ.1.AND.ibc2.EQ.1.AND.
     &                      ibc3.EQ.1.AND.ibc4.EQ.1) THEN
                           cellNatureField(et) = -1
                           isMasked = 1
                           GOTO 5
                        ENDIF

                        IF (ibc1.EQ.0.AND.ibc2.EQ.0.AND.
     &                      ibc3.EQ.0.AND.ibc4.EQ.0) THEN
c                           cellNatureField(et) = 0
                           isMasked = 1
                           GOTO 5
                        ENDIF

                        IF (pos.EQ.1.AND.zp.LE.zi) THEN
                           cellNatureField(et) = -1
                           isMasked = 1
                           GOTO 5
                        ENDIF

                        IF (pos.EQ.-1.AND.zp.GE.zi) THEN
                           cellNatureField(et) = -1
                           isMasked = 1
                           GOTO 5
                        ENDIF

 6                      isMasked = 1
c                        cellNatureField(et) = 0
                        GOTO 5

                     ENDIF
C     
                  ENDDO
               ENDDO
 5             CONTINUE
