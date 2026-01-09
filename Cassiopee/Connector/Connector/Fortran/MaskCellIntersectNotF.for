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

C-------------------------------------------------------------------------------
C common part for unstructured 3D blanking routines for inverted blanking 
C (isNot = 1) and in case of cell intersect blanking type
C-------------------------------------------------------------------------------
C
C     Ray index corresponding to point cell
         iraymin = (xmincell - xmin)/hiray
         iraymax = (xmaxcell - xmin)/hiray+1
         IF ((iraymin.LT.0).OR.(iraymax.GT.niray-1)) THEN
            cellNatureField(et) = -1
            isMasked = 1
            GOTO 6
         ENDIF  
         iraymin = MAX( 0, iraymin )
         iraymax = MIN( niray-1, iraymax )
         
         jraymin = (ymincell - ymin)/hjray
         jraymax = (ymaxcell - ymin)/hjray+1
         IF ((jraymin.LT.0).OR.(jraymax.GT.njray-1)) THEN
            cellNatureField(et) = -1
            isMasked = 1
            GOTO 6
         ENDIF 
         jraymin = MAX( 0, jraymin )
         jraymax = MIN( njray-1, jraymax )
         
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
                     GOTO 4
                  ENDIF
                  npmin = npmin+1
               ENDDO
 4             npmax = 0
               DO l = ibeg, iend
                  IF (z(l).GT.zmaxcell) THEN
                     GOTO 5
                  ENDIF
                  npmax = npmax+1
               ENDDO
 5             IF ( ((npmin/2)*2.EQ.npmin).OR.
     &              ((npmax/2)*2.EQ.npmax).OR.
     &              (npmin.NE.npmax) ) THEN
C     OUTSIDE
                  cellNatureField(et) = -1
                  isMasked = 1
                  GOTO 6
               ENDIF
C     
            ENDDO
         ENDDO
 6       CONTINUE
