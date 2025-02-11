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

C  ============================================================================
C  Calcul du rotationnel moyen d un champ (u,v,w) sur une cellule tetra
C  IN: velo: vecteur dont le rotationnel est a calculer. Defini sur le maillage
C  OUT: rotu,rotv,rotw: rotationnel moyen de (u,v,w) sur la cellule
C  ============================================================================
      SUBROUTINE k6compmeancurloftetracell( npts, ind1, ind2, ind3,ind4,
     &                                      velox, veloy, veloz,  
     &                                      xt, yt, zt, rotu, rotv,rotw)

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E npts            ! nb de noeuds ds le maillage
      INTEGER_E ind1, ind2, ind3, ind4 
      REAL_E xt(0:npts-1)
      REAL_E yt(0:npts-1)
      REAL_E zt(0:npts-1)
      REAL_E velox(0:npts-1)    !vecteur dont le rotationnel est a calculer ds la cellule
      REAL_E veloy(0:npts-1)
      REAL_E veloz(0:npts-1)
C_OUT
      REAL_E  rotu, rotv, rotw  ! curl of u,v,w 
C_LOCAL
      INTEGER_E edge, l
      INTEGER_E cnloc(1,4)
      REAL_E snx(1,4)
      REAL_E sny(1,4)
      REAL_E snz(1,4)
      REAL_E surf(1,4)
      REAL_E xloc(4)
      REAL_E yloc(4)
      REAL_E zloc(4)
      REAL_E xint(4)
      REAL_E yint(4)
      REAL_E zint(4)
      REAL_E x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
      REAL_E vol, curlx, curly, curlz
      REAL_E ONE_THIRD
      REAL_E ux(4) ! vitesse aux interfaces %x
      REAL_E uy(4) ! vitesse aux interfaces %y 
      REAL_E uz(4) ! vitesse aux interfaces %z
C     
      ONE_THIRD = ONE/THREE
C
C-----------------------------------------------------------------------------
C
      x1 = xt(ind1)
      x2 = xt(ind2)
      x3 = xt(ind3)
      x4 = xt(ind4)
     
      y1 = yt(ind1)
      y2 = yt(ind2)
      y3 = yt(ind3)
      y4 = yt(ind4)

      z1 = zt(ind1)
      z2 = zt(ind2)
      z3 = zt(ind3)
      z4 = zt(ind4)
      
C     face 123
      xint(1) = ONE_THIRD * ( x1 + x2 + x3 )
      yint(1) = ONE_THIRD * ( y1 + y2 + y3 )
      zint(1) = ONE_THIRD * ( z1 + z2 + z3 )
C     face 124
      xint(2) = ONE_THIRD * ( x1 + x2 + x4 )
      yint(2) = ONE_THIRD * ( y1 + y2 + y4 )
      zint(2) = ONE_THIRD * ( z1 + z2 + z4 )
C     face 234
      xint(3) = ONE_THIRD * ( x2 + x3 + x4 )
      yint(3) = ONE_THIRD * ( y2 + y3 + y4 )
      zint(3) = ONE_THIRD * ( z2 + z3 + z4 )
C     face 134
      xint(4) = ONE_THIRD * ( x1 + x3 + x4 )
      yint(4) = ONE_THIRD * ( y1 + y3 + y4 )
      zint(4) = ONE_THIRD * ( z1 + z3 + z4 )

c     calcul des normales/surfaces
      xloc(1) = x1
      xloc(2) = x2
      xloc(3) = x3
      xloc(4) = x4

      yloc(1) = y1
      yloc(2) = y2
      yloc(3) = y3
      yloc(4) = y4

      zloc(1) = z1
      zloc(2) = z2
      zloc(3) = z3
      zloc(4) = z4 
     
      cnloc(1,1) = 1
      cnloc(1,2) = 2
      cnloc(1,3) = 3
      cnloc(1,4) = 4
      
      call k6comptetrasurf(4, 1, cnloc, xloc, yloc, zloc,
     &     snx, sny, snz, surf)

c     calcul du volume
      vol = ZERO

      DO edge = 1, 4
         vol = vol + xint(edge) * snx(1,edge) + 
     &        yint(edge) * sny(1,edge) + zint(edge) * snz(1,edge)       
      ENDDO

      vol = ONE_THIRD * vol

C-------------------------------------------------------------------------
C     Calcul des vecteurs vitesse aux interfaces
C     interfaces  123, 124, 234, 134 resp
c
      ux(1) = velox(ind1) + velox(ind2) + velox(ind3)
      uy(1) = veloy(ind1) + veloy(ind2) + veloy(ind3) 
      uz(1) = veloz(ind1) + veloz(ind2) + veloz(ind3)
C
      ux(2) = velox(ind1) + velox(ind2) + velox(ind4)
      uy(2) = veloy(ind1) + veloy(ind2) + veloy(ind4)
      uz(2) = veloz(ind1) + veloz(ind2) + veloz(ind4)
C
      ux(3) = velox(ind2) + velox(ind3) + velox(ind4)
      uy(3) = veloy(ind2) + veloy(ind3) + veloy(ind4)
      uz(3) = veloz(ind2) + veloz(ind3) + veloz(ind4)
C
      ux(4) = velox(ind1) + velox(ind3) + velox(ind4)
      uy(4) = veloy(ind1) + veloy(ind3) + veloy(ind4)
      uz(4) = veloz(ind1) + veloz(ind3) + veloz(ind4)
C
C
C-------------------------------------------------------------------------
      curlx = 0.
      curly = 0.
      curlz = 0.
      DO l = 1,4
         curlx = curlx + uy(l)*snz(1,l) - uz(l)*sny(1,l)
         curly = curly + uz(l)*snx(1,l) - ux(l)*snz(1,l)
         curlz = curlz + ux(l)*sny(1,l) - uy(l)*snx(1,l)
      ENDDO
           
      vol = - ONE_THIRD * vol   ! prise en compte de vmoy = (sum v)/3
      rotu = vol * curlx
      rotv = vol * curly
      rotw = vol * curlz

      END

C================== Post/Fortran/CompMeanCurlOfTetraCellF.for ============
