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

C Calcul du rotationnel d'un vecteur defini aux noeuds d'une grille 
C non structuree
C retourne le rotationnel defini aux centres des elts
C CAS 3D
C ============================================================================
      SUBROUTINE k6compunstrcurl(dim, npts, nelts, nedges, nnodes, cn,  
     &                           xt, yt, zt, snx, sny, snz, surf, vol,
     &                           uintx, uinty, uintz, ux, uy, uz, 
     &                           xint, yint, zint,
     &                           rotx, roty, rotz )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E dim             ! dimension du pb
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
      REAL_E rotx(0:nelts-1)    ! rotationnel de ux
      REAL_E roty(0:nelts-1)    ! rotationnel de uy
      REAL_E rotz(0:nelts-1)    ! rotationnel de uz

C_LOCAL
      REAL_E uintx(0:nelts-1,nedges) ! composante %x du vecteur aux facettes
      REAL_E uinty(0:nelts-1,nedges)
      REAL_E uintz(0:nelts-1,nedges)
      REAL_E vol(0:nelts-1)
      REAL_E surf(0:nelts-1,nedges) ! surfaces des facettes
      REAL_E snx(0:nelts-1,nedges) ! normales aux facettes %x
      REAL_E sny(0:nelts-1,nedges) ! normales aux facettes %y
      REAL_E snz(0:nelts-1,nedges) ! normales aux facettes %z
      REAL_E xint(0:nelts-1, nedges) ! coordonnees du centre des facettes 
      REAL_E yint(0:nelts-1, nedges)
      REAL_E zint(0:nelts-1, nedges)
      INTEGER_E fi, eti
      REAL_E curlx, curly, curlz, vinv
      REAL_E sxi, syi, szi

C ------------------------------------------------------------------
C     calcul de la surface + volume des elts
      call k6compunstrmetric(npts, nelts, nedges, nnodes, cn, 
     &     xt, yt, zt, xint, yint, zint, snx, sny, snz, surf, vol)
C     calcul des composantes du vecteur aux centres des facettes 
      call k6unstrnodes2faces(dim, npts, nelts, nedges, nnodes,  
     &     cn, ux, uintx)
      call k6unstrnodes2faces(dim, npts, nelts, nedges, nnodes,  
     &     cn, uy, uinty)
      call k6unstrnodes2faces(dim, npts, nelts, nedges, nnodes,  
     &     cn, uz, uintz)
         
c     calcul du rotationnel au centre des elts

!$OMP PARALLEL PRIVATE(eti,curlx,curly,curlz,vinv,fi,sxi,syi,szi)
!$OMP DO
      DO eti = 0, nelts-1
         
         curlx = ZERO
         curly = ZERO
         curlz = ZERO

         vinv = -ONE/MAX(vol(eti), E_MIN_VOL)

         DO fi = 1, nedges
            sxi = snx(eti,fi) 
            syi = sny(eti,fi) 
            szi = snz(eti,fi) 

            curlx = curlx + uinty(eti,fi)* szi - uintz(eti,fi)* syi
            curly = curly + uintz(eti,fi)* sxi - uintx(eti,fi)* szi
            curlz = curlz + uintx(eti,fi)* syi - uinty(eti,fi)* sxi
         ENDDO
         
         rotx(eti) = vinv * curlx
         roty(eti) = vinv * curly
         rotz(eti) = vinv * curlz
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
      END

C=====================Post/Fortran/CompUnstrCurlF.for=====================

