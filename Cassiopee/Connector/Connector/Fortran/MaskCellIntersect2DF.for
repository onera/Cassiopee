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

C------------------------------------------------------------------------------
C common part for unstructured 2D blanking routines for classical blanking (isNot = 0)
C and in case of cell intersect blanking type
C------------------------------------------------------------------------------
C
C Ray index corresponding to point cell
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
 1             npmax = 0

               DO l = ibeg, iend
                  IF (z(l).GT.ymaxcell) THEN
                     GOTO 2
                  ENDIF
                  npmax = npmax+1
               ENDDO
 2             IF ( ((npmin/2)*2.NE.npmin).OR.
     &              ((npmax/2)*2.NE.npmax).OR.
     &              (npmin.NE.npmax) ) THEN
C     INSIDE
                  cellNatureField(et) = -1
                  isMasked = 1
                  GOTO 3
               ENDIF
C     
            ENDDO
 3          CONTINUE
