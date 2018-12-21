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

C  ============================================================================
C Write the generic.grd file according to USURP         
C this contains information about mesh coordinates for a given block
C  ============================================================================
      SUBROUTINE k6writegrid(nzone, nptsmax, nit, njt, nkt, xb, yb, zb)


      IMPLICIT NONE
C==============================================================================
C_IN
      INTEGER_E nzone
      INTEGER_E nptsmax
      INTEGER_E nit(nzone)
      INTEGER_E njt(nzone)
      INTEGER_E nkt(nzone)
      REAL_E xb(nptsmax,nzone)
      REAL_E yb(nptsmax,nzone)
      REAL_E zb(nptsmax,nzone)
C_LOCAL
      INTEGER_E i,n, npts 
C==============================================================================
      open(unit=3,file="generic.grd",status="unknown")

      write(unit=3,fmt=*) nzone
      write(unit=3,fmt=*)(nit(n),njt(n),nkt(n),n=1,nzone)

      do n=1,nzone
         npts = nit(n)*njt(n)*nkt(n)
         write(unit=3,fmt=*) (xb(i,n),i=1,npts), 
     &        (yb(i,n),i=1,npts), 
     &        (zb(i,n),i=1,npts)
      enddo
      close(3)
      return
      END
