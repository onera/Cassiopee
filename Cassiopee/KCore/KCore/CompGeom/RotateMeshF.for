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
C ============================================================================
C  ============================================================================
C  Make a rotation of a mesh
C  ============================================================================
      SUBROUTINE k6rotatemesh(dim, center, axis, teta, x0, y0, z0)
C
      IMPLICIT NONE
C==============================================================================
C_IN
      INTEGER_E dim                     ! mesh size
      REAL_E center(1:3)                ! center of rotation
      REAL_E axis(1:3)                  ! rotation Vector
      REAL_E teta                       ! angle
C_OUT
      REAL_E x0(1:dim), y0(1:dim), z0(1:dim) ! rotated mesh point
C_LOCAL
      INTEGER_E ind
      REAL_E unx, uny, unz
      REAL_E norm
      REAL_E px, py, pz
      REAL_E rx, ry, rz
      REAL_E e0, e1, e2, e3
      REAL_E a1, a2, steta, stetas2
C==============================================================================
C nx,ny,nz must be unit Vector
      norm = axis(1)*axis(1) + axis(2)*axis(2) + axis(3)*axis(3)
      IF (norm.LE.1.e-12) THEN
         WRITE(*,*) 'Error: k6rotatemesh: nx,ny,nz has null norm.'
         RETURN
      ENDIF
         
      norm = 1./SQRT(norm)
      unx = axis(1)*norm
      uny = axis(2)*norm
      unz = axis(3)*norm
      steta = sin(teta)
      stetas2 = sin(teta*0.5D0)
         
C quaternion
      e0 = cos(teta*0.5D0)
      e1 = -unx*stetas2
      e2 = -uny*stetas2
      e3 = -unz*stetas2
      a1 = e0*e0-e1*e1-e2*e2-e3*e3

      DO ind = 1, dim         
         rx = x0(ind)-center(1)
         ry = y0(ind)-center(2)
         rz = z0(ind)-center(3)
         a2 = e1*rx+e2*ry+e3*rz
         px = a1*rx+2*e1*a2-(ry*unz-rz*uny)*steta
         py = a1*ry+2*e2*a2-(rz*unx-rx*unz)*steta
         pz = a1*rz+2*e3*a2-(rx*uny-ry*unx)*steta

         x0(ind) = center(1) + px
         y0(ind) = center(2) + py
         z0(ind) = center(3) + pz
      ENDDO
      END

C  ============================================================================
C  Make a rotation of a mesh
C same function as previous, but different interface and omp present
C  ============================================================================
      SUBROUTINE k6rotatemesh2(npts, x, y, z,
     &                    xc, yc, zc, nx, ny, nz,
     &                    teta, 
     &                    xo, yo, zo)
C
      IMPLICIT NONE
C==============================================================================
C_IN
      INTEGER_E npts            ! mesh size
      REAL_E x(1:npts)          ! mesh coordinates
      REAL_E y(1:npts)
      REAL_E z(1:npts)
      REAL_E xc,yc,zc           ! center of rotation
      REAL_E nx,ny,nz           ! rotation vector
      REAL_E teta               ! angle
C_OUT
      REAL_E xo(1:npts)         ! rotated mesh
      REAL_E yo(1:npts)
      REAL_E zo(1:npts)
C_LOCAL
      INTEGER_E ind
      REAL_E unx, uny, unz
      REAL_E norm
      REAL_E px, py, pz
      REAL_E rx, ry, rz
      REAL_E e0, e1, e2, e3
      REAL_E a1, a2, sinteta, sinteta5
C==============================================================================
C     nx,ny,nz must be unit vector
      norm = nx*nx+ny*ny+nz*nz
      IF (norm.LE.1.e-12) THEN
         WRITE(*,*) 'Error: rotate: nx,ny,nz has null norm.'
         RETURN
      ENDIF
         
      norm = 1.D0/SQRT(norm)
      unx = nx*norm
      uny = ny*norm
      unz = nz*norm
         
      sinteta = sin(teta)
      sinteta5 = sin(teta*0.5D0)
C quaternion
      e0 = cos(teta*0.5D0)
      e1 = -unx*sinteta5
      e2 = -uny*sinteta5
      e3 = -unz*sinteta5
      a1 = e0*e0-e1*e1-e2*e2-e3*e3

!$OMP PARALLEL PRIVATE(ind, rx, ry, rz, a2, px, py, pz)
!$OMP DO 
      DO ind = 1, npts
         rx = x(ind)-xc
         ry = y(ind)-yc
         rz = z(ind)-zc
         a2 = e1*rx+e2*ry+e3*rz
         px = a1*rx+2*e1*a2-(ry*unz-rz*uny)*sinteta
         py = a1*ry+2*e2*a2-(rz*unx-rx*unz)*sinteta
         pz = a1*rz+2*e3*a2-(rx*uny-ry*unx)*sinteta
         xo(ind) = xc+px
         yo(ind) = yc+py
         zo(ind) = zc+pz
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      END

C============================= CompGeom/RotateMeshF.for ====================
