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

C Calcul du rotationnel d un vecteur defini aux noeuds d une grille surfacique
C non structuree 
C retourne le rotationnel du vecteur defini aux centres des cellules
C ============================================================================
      SUBROUTINE k6compunstrcurl2d(npts, nelts, nedges, nnodes, cn,  
     &                             xt ,yt, zt, snx, sny, snz, surf, vol,
     &                             ux, uy, uz, rotx, roty, rotz)
      
      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E npts            ! nb de points ds le maillage
      INTEGER_E nelts           ! nb d elts dans le maillage
      INTEGER_E nedges          ! nb de facettes par elt
      INTEGER_E nnodes          ! nb de noeuds par elt 
      INTEGER_E cn(0:nelts-1,nnodes) ! connectivite elt->noeud
      REAL_E ux(0:npts-1)       ! composante %x du vecteur
      REAL_E uy(0:npts-1)       ! 
      REAL_E uz(0:npts-1)       ! 
      REAL_E xt(0:npts-1)       ! coordonnees des noeuds de la grille
      REAL_E yt(0:npts-1)       ! coordonnees des noeuds de la grille
      REAL_E zt(0:npts-1)       ! coordonnees des noeuds de la grille
C_OUT
      REAL_E rotx(0:nelts-1)  ! rotationnel de ux
      REAL_E roty(0:nelts-1)  ! rotationnel de uy
      REAL_E rotz(0:nelts-1)  ! rotationnel de uz

C_LOCAL
      INTEGER_E fi, eti, nfaces
      INTEGER_E indA, indB, indC, indD
   
      REAL_E vol(0:nelts-1)
      REAL_E surf(0:nelts-1,1) ! surfaces des facettes
      REAL_E snx(0:nelts-1,1) ! normales aux facettes %x
      REAL_E sny(0:nelts-1,1) ! normales aux facettes %y
      REAL_E snz(0:nelts-1,1) ! normales aux facettes %z
      REAL_E curlx, curly, curlz, vinv
      REAL_E sxi, syi, szi, vx, vy, vz
      REAL_E nx, ny, nz, nn, n1x, n1y, n1z
      REAL_E xAB, yAB, zAB, xBC, yBC, zBC, xCD, yCD, zCD, xDA, yDA, zDA
      REAL_E xCA, yCA, zCA

C ------------------------------------------------------------------
     
C     calcul de la surface des elts
      nfaces = 1

      call k6unstructsurf(npts, nelts, nfaces, nnodes, cn, 
     &                    xt, yt, zt, snx, sny, snz, surf)
 
C     cas tri
      IF (nnodes .eq. 3) THEN 
         DO eti = 0, nelts-1
            
            nx = snx(eti,1)
            ny = sny(eti,1)
            nz = snz(eti,1)
            nn = sqrt(nx*nx+ny*ny+nz*nz)

            vinv = 2D0 * surf(eti,1) * nn 
            vinv = -ONE/MAX(vinv, E_MIN_VOL)

            indA = cn(eti,1)-1
            indB = cn(eti,2)-1
            indC = cn(eti,3)-1           

            xAB = xt(indB)-xt(indA)
            yAB = yt(indB)-yt(indA)
            zAB = zt(indB)-zt(indA)
            
            xBC = xt(indC)-xt(indB)
            yBC = yt(indC)-yt(indB)
            zBC = zt(indC)-zt(indB)
            
            xCA = xt(indA)-xt(indC)
            yCA = yt(indA)-yt(indC)
            zCA = zt(indA)-zt(indC)
           
            curlx = ZERO
            curly = ZERO
            curlz = ZERO

            n1x = yAB*nz - zAB*ny
            n1y = zAB*nx - xAB*nz
            n1z = xAB*ny - yAB*nx

            vx = ux(indA) + ux(indB)
            vy = uy(indA) + uy(indB)
            vz = uz(indA) + uz(indB)

            curlx = curlx + vy*n1z - vz*n1y
            curly = curly + vz*n1x - vx*n1z
            curlz = curlz + vx*n1y - vy*n1x
     
            n1x = yBC*nz - zBC*ny
            n1y = zBC*nx - xBC*nz
            n1z = xBC*ny - yBC*nx
            
            vx = ux(indC) + ux(indB)
            vy = uy(indC) + uy(indB)
            vz = uz(indC) + uz(indB)
            
            curlx = curlx + vy*n1z - vz*n1y
            curly = curly + vz*n1x - vx*n1z
            curlz = curlz + vx*n1y - vy*n1x
     
            n1x = yCA*nz - zCA*ny
            n1y = zCA*nx - xCA*nz
            n1z = xCA*ny - yCA*nx

            vx = ux(indC) + ux(indA)
            vy = uy(indC) + uy(indA)
            vz = uz(indC) + uz(indA)

            curlx = curlx + vy*n1z - vz*n1y
            curly = curly + vz*n1x - vx*n1z
            curlz = curlz + vx*n1y - vy*n1x
     
            rotx(eti) = vinv * curlx
            roty(eti) = vinv * curly
            rotz(eti) = vinv * curlz
         ENDDO

      ELSE IF (nnodes .eq. 4) THEN
         
         DO eti = 0, nelts-1
            indA = cn(eti,1)-1
            indB = cn(eti,2)-1
            indC = cn(eti,3)-1          
            indD = cn(eti,4)-1
            
            nx = snx(eti,1)
            ny = sny(eti,1)
            nz = snz(eti,1)
            nn = sqrt(nx*nx+ny*ny+nz*nz)

            vinv = 2 * surf(eti,1) * nn 
            vinv = -ONE/MAX(vinv, E_MIN_VOL)

            xAB = xt(indB)-xt(indA)
            yAB = yt(indB)-yt(indA)
            zAB = zt(indB)-zt(indA)
            
            xBC = xt(indC)-xt(indB)
            yBC = yt(indC)-yt(indB)
            zBC = zt(indC)-zt(indB)
         
            xCD = xt(indD)-xt(indC)
            yCD = yt(indD)-yt(indC)
            zCD = zt(indD)-zt(indC)
            
            xDA = xt(indA)-xt(indD)
            yDA = yt(indA)-yt(indD)
            zDA = zt(indA)-zt(indD)
C
            curlx = ZERO
            curly = ZERO
            curlz = ZERO
C
            n1x = yAB*nz - zAB*ny
            n1y = zAB*nx - xAB*nz
            n1z = xAB*ny - yAB*nx
            
            vx = ux(indA) + ux(indB)
            vy = uy(indA) + uy(indB)
            vz = uz(indA) + uz(indB)
            
            curlx = curlx + vy*n1z - vz*n1y
            curly = curly + vz*n1x - vx*n1z
            curlz = curlz + vx*n1y - vy*n1x
C     
            n1x = yBC*nz - zBC*ny
            n1y = zBC*nx - xBC*nz
            n1z = xBC*ny - yBC*nx
            
            vx = ux(indC) + ux(indB)
            vy = uy(indC) + uy(indB)
            vz = uz(indC) + uz(indB)
            curlx = curlx + vy*n1z - vz*n1y
            curly = curly + vz*n1x - vx*n1z
            curlz = curlz + vx*n1y - vy*n1x
C     
            n1x = yCD*nz - zCD*ny
            n1y = zCD*nx - xCD*nz
            n1z = xCD*ny - yCD*nx
            
            vx = ux(indC) + ux(indD)
            vy = uy(indC) + uy(indD)
            vz = uz(indC) + uz(indD)
            curlx = curlx + vy*n1z - vz*n1y
            curly = curly + vz*n1x - vx*n1z
            curlz = curlz + vx*n1y - vy*n1x
C     
            n1x = yDA*nz - zDA*ny
            n1y = zDA*nx - xDA*nz
            n1z = xDA*ny - yDA*nx
            
            vx = ux(indD) + ux(indA)
            vy = uy(indD) + uy(indA)
            vz = uz(indD) + uz(indA)
            
            curlx = curlx + vy*n1z - vz*n1y
            curly = curly + vz*n1x - vx*n1z
            curlz = curlz + vx*n1y - vy*n1x
C     
            rotx(eti) = vinv * curlx
            roty(eti) = vinv * curly
            rotz(eti) = vinv * curlz
         ENDDO    
      ELSE 
         WRITE(*,*) ' CompUnstrCurl2D: wrong value of nnodes'
         STOP
      ENDIF
      END

C ================== Post/Fortran/CompUnstrCurl2DF.for ===================
