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

C Calcul de la divergence d'un champ defini aux noeuds d une grille non structuree
C retourne la divergence defini aux centres des elts
C CAS 2D
C ============================================================================
      SUBROUTINE k6compunstrdiv2d(npts, nelts, nnodes,
     &                            cn, xt, yt, zt,
     &                            fieldX, fieldY, fieldZ,
     &                            snx, sny, snz, surf,
     &                            div)

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E npts            ! nb de points ds le maillage
      INTEGER_E nelts           ! nb d elts dans le maillage
      INTEGER_E nnodes          ! nb de noeuds par elt
      INTEGER_E cn(0:nelts-1,nnodes) ! connectivite elt->noeud
      REAL_E fieldX(0:npts-1)   ! champ (x-comp) defini aux noeuds
      REAL_E fieldY(0:npts-1)   ! champ (y-comp) defini aux noeuds
      REAL_E fieldZ(0:npts-1)   ! champ (z-comp) defini aux noeuds
      REAL_E xt(0:npts-1)       ! coordonnees des noeuds de la grille
      REAL_E yt(0:npts-1)       ! coordonnees des noeuds de la grille
      REAL_E zt(0:npts-1)       ! coordonnees des noeuds de la grille

C_OUT
      REAL_E div(0:nelts-1) ! divergence du champ vectoriel

C_LOCAL
      INTEGER_E nfaces, et
      INTEGER_E indA, indB, indC, indD
      REAL_E  nn, vinv
      REAL_E xAB, yAB, zAB, xBC, yBC, zBC, xCD, yCD, zCD, xDA, yDA, zDA
      REAL_E xCA, yCA, zCA
      REAL_E nx, n1x, n2x, n3x, n4x
      REAL_E ny, n1y, n2y, n3y, n4y
      REAL_E nz, n1z, n2z, n3z, n4z
      REAL_E fxAB, fxBC, fxCD, fxDA, fxCA
      REAL_E fyAB, fyBC, fyCD, fyDA, fyCA
      REAL_E fzAB, fzBC, fzCD, fzDA, fzCA
      REAL_E gradxx, gradyy, gradzz
      REAL_E snx(0:nelts-1,1)
      REAL_E sny(0:nelts-1,1)
      REAL_E snz(0:nelts-1,1)
      REAL_E surf(0:nelts-1,1) ! surface de l element (=aire)

C-------------------------------------------------------------------------
C     calcul de la surface des elts
      nfaces = 1

      call k6unstructsurf(npts, nelts, nfaces, nnodes, cn,
     &                    xt, yt, zt, snx, sny, snz, surf)

C     CAS TRI
      IF (nnodes .eq. 3) THEN
         DO et = 0, nelts-1
            indA = cn(et,1)-1
            indB = cn(et,2)-1
            indC = cn(et,3)-1

            nx = snx(et,1)
            ny = sny(et,1)
            nz = snz(et,1)
            nn = sqrt(nx*nx+ny*ny+nz*nz)
            vinv = 2D0 * surf(et,1) * nn
            vinv = ONE/MAX(vinv, E_MIN_VOL)

            xAB = xt(indB)-xt(indA)
            yAB = yt(indB)-yt(indA)
            zAB = zt(indB)-zt(indA)

            xBC = xt(indC)-xt(indB)
            yBC = yt(indC)-yt(indB)
            zBC = zt(indC)-zt(indB)

            xCA = xt(indA)-xt(indC)
            yCA = yt(indA)-yt(indC)
            zCA = zt(indA)-zt(indC)

            n1x = yAB*nz - zAB*ny
            n1y = zAB*nx - xAB*nz
            n1z = xAB*ny - yAB*nx

            n2x = yBC*nz - zBC*ny
            n2y = zBC*nx - xBC*nz
            n2z = xBC*ny - yBC*nx

            n3x = yCA*nz - zCA*ny
            n3y = zCA*nx - xCA*nz
            n3z = xCA*ny - yCA*nx

            fxAB = fieldX(indA) + fieldX(indB)
            fxBC = fieldX(indB) + fieldX(indC)
            fxCA = fieldX(indC) + fieldX(indA)

            fyAB = fieldY(indA) + fieldY(indB)
            fyBC = fieldY(indB) + fieldY(indC)
            fyCA = fieldY(indC) + fieldY(indA)

            fzAB = fieldZ(indA) + fieldZ(indB)
            fzBC = fieldZ(indB) + fieldZ(indC)
            fzCA = fieldZ(indC) + fieldZ(indA)

            gradxx = fxAB*n1x + fxBC*n2x + fxCA*n3x
            gradyy = fyAB*n1y + fyBC*n2y + fyCA*n3y
            gradzz = fzAB*n1z + fzBC*n2z + fzCA*n3z

            div(et) = vinv *( gradxx + gradyy + gradzz )

         ENDDO
      ELSE                      ! Quad

         DO et = 0, nelts-1
            indA = cn(et,1)-1
            indB = cn(et,2)-1
            indC = cn(et,3)-1
            indD = cn(et,4)-1

            nx = snx(et,1)
            ny = sny(et,1)
            nz = snz(et,1)

            nn = sqrt(nx*nx+ny*ny+nz*nz)
            vinv = 2 * surf(et,1) * nn
            vinv = ONE/MAX(vinv, E_MIN_VOL)

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

            n1x = yAB*nz - zAB*ny
            n1y = zAB*nx - xAB*nz
            n1z = xAB*ny - yAB*nx

            n2x = yBC*nz - zBC*ny
            n2y = zBC*nx - xBC*nz
            n2z = xBC*ny - yBC*nx

            n3x = yCD*nz - zCD*ny
            n3y = zCD*nx - xCD*nz
            n3z = xCD*ny - yCD*nx

            n4x = yDA*nz - zDA*ny
            n4y = zDA*nx - xDA*nz
            n4z = xDA*ny - yDA*nx

            fxAB = fieldX(indA) + fieldX(indB)
            fxBC = fieldX(indB) + fieldX(indC)
            fxCD = fieldX(indC) + fieldX(indD)
            fxDA = fieldX(indD) + fieldX(indA)

            fyAB = fieldY(indA) + fieldY(indB)
            fyBC = fieldY(indB) + fieldY(indC)
            fyCD = fieldY(indC) + fieldY(indD)
            fyDA = fieldY(indD) + fieldY(indA)

            fzAB = fieldZ(indA) + fieldZ(indB)
            fzBC = fieldZ(indB) + fieldZ(indC)
            fzCD = fieldZ(indC) + fieldZ(indD)
            fzDA = fieldZ(indD) + fieldZ(indA)

            gradxx = fxAB*n1x + fxBC*n2x + fxCD*n3x + fxDA*n4x
            gradyy = fyAB*n1y + fyBC*n2y + fyCD*n3y + fyDA*n4y
            gradzz = fzAB*n1z + fzBC*n2z + fzCD*n3z + fzDA*n4z

            div(et) = vinv *( gradxx + gradyy + gradzz )

         ENDDO
      ENDIF
      END

C=====================Post/Fortran/CompUnstrDiv2DF.for======================

