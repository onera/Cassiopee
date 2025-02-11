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

C depend du type d'element considere 
C les normales aux surfaces sont orientees vers l'exterieur de l'element
C ============================================================================
      SUBROUTINE k6unstructsurf(npts, nelts, nedges, nnodes, cn, 
     &                          xt, yt, zt,
     &                          surfnx, surfny, surfnz, surface) 
      IMPLICIT NONE
      
# include "Def/DefFortranConst.h"
C==============================================================================
C_IN
      INTEGER_E npts            ! nb de pts dans le maillage
      INTEGER_E nelts           ! nb d elements
      INTEGER_E nedges          ! nb de facettes par elemt
      INTEGER_E nnodes          ! nb de noeuds par elemt
      INTEGER_E cn(0:nelts-1,nnodes) ! connectivite elt->noeud
      REAL_E xt(0:npts-1)       !coordonnees des pts de la grille
      REAL_E yt(0:npts-1)       !coordonnees des pts de la grille
      REAL_E zt(0:npts-1)       !coordonnees des pts de la grille

C_OUT
      REAL_E surfnx(0:nelts-1,nedges) ! normale aux surfaces %x
      REAL_E surfny(0:nelts-1,nedges) ! normale aux surfaces %y
      REAL_E surfnz(0:nelts-1,nedges) ! normale aux surfaces %z
      REAL_E surface(0:nelts-1,nedges) ! surface des elts  
      
C==============================================================================
      IF (nedges .eq. 1 .and. nnodes .eq. 3) THEN ! TRI
         call k6comptrisurf(npts, nelts, cn, xt, yt, zt, 
     &                      surfnx, surfny, surfnz, surface)

      ELSE IF (nedges .eq. 1 .and. nnodes .eq. 4) THEN ! QUAD
         call k6compquadsurf(npts, nelts, cn, xt, yt, zt,
     &                       surfnx, surfny, surfnz, surface)

      ELSE IF (nedges .eq. 4 .and. nnodes .eq. 4) THEN ! TETRA
         call k6comptetrasurf(npts, nelts, cn, xt, yt, zt,
     &                      surfnx, surfny, surfnz, surface)

      ELSE IF (nedges .eq. 6 .and. nnodes .eq. 8) THEN ! HEXA
         call k6comphexasurf(npts, nelts, cn, xt, yt, zt,
     &                       surfnx, surfny, surfnz, surface)

      ELSE IF (nedges .eq. 5 .and. nnodes .eq. 6) THEN ! PENTA
         call k6comppentasurf(npts, nelts, cn, xt, yt, zt,
     &                        surfnx, surfny, surfnz, surface)

      ELSE IF (nedges .eq. 5 .and. nnodes .eq. 5) THEN ! PYRA
         call k6comppyrasurf(npts, nelts, cn, xt, yt, zt,
     &                        surfnx, surfny, surfnz, surface)

      ELSE 
         write(*,*) 'Error: in KCore/Metric/CompUnstrSurfF.for.'
         write(*,*) 'Unknown type of elements.'
         STOP
      ENDIF
      END

C-----------------------------------------------------------------------------
      SUBROUTINE k6comptrisurf(npts, nelts, cn, xt, yt, zt,
     &                         surfnx, surfny, surfnz, surface)
      IMPLICIT NONE
#include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E nelts
      INTEGER_E npts
      INTEGER_E cn(0:nelts-1,3)
      REAL_E xt(0:npts-1)
      REAL_E yt(0:npts-1)
      REAL_E zt(0:npts-1)

C_OUT
      REAL_E surfnx(0:nelts-1,1)
      REAL_E surfny(0:nelts-1,1)
      REAL_E surfnz(0:nelts-1,1)
      REAL_E surface(0:nelts-1,1)
C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind3
      REAL_E surf, surfx, surfy, surfz
      REAL_E l1x, l1y, l1z, l2x, l2y, l2z
C==============================================================================

      DO i = 0, nelts-1
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         ind3 = cn(i,3)-1
         l1x = xt(ind1)-xt(ind2)
         l1y = yt(ind1)-yt(ind2)
         l1z = zt(ind1)-zt(ind2)
         l2x = xt(ind1)-xt(ind3)
         l2y = yt(ind1)-yt(ind3)
         l2z = zt(ind1)-zt(ind3)
         surfx = (l1y*l2z-l1z*l2y)
         surfy = (l1z*l2x-l1x*l2z)
         surfz = (l1x*l2y-l1y*l2x)
         surf = sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
         surfnx(i,1) = ONE_HALF * surfx
         surfny(i,1) = ONE_HALF * surfy
         surfnz(i,1) = ONE_HALF * surfz
         surface(i,1) = ONE_HALF* surf
      ENDDO
      END
C-----------------------------------------------------------------------------
      SUBROUTINE k6compquadsurf(npts, nelts, cn, xt, yt, zt, 
     &                         surfnx, surfny, surfnz, surface)

      IMPLICIT NONE
#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E nelts
      INTEGER_E npts
      INTEGER_E cn(0:nelts-1,4)
      REAL_E xt(0:npts-1)
      REAL_E yt(0:npts-1)
      REAL_E zt(0:npts-1)

C_OUT
      REAL_E surfnx(0:nelts-1,1)
      REAL_E surfny(0:nelts-1,1)
      REAL_E surfnz(0:nelts-1,1)
      REAL_E surface(0:nelts-1,1)
C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind3, ind4
      REAL_E surfx, surfy, surfz
      REAL_E surf1x, surf1y, surf1z, surf2x, surf2y, surf2z
      REAL_E l1x, l1y, l1z, l2x, l2y, l2z

      DO i = 0, nelts-1
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         ind3 = cn(i,3)-1
         ind4 = cn(i,4)-1
C     AB x AC
         l1x = xt(ind2) - xt(ind1)
         l1y = yt(ind2) - yt(ind1)
         l1z = zt(ind2) - zt(ind1)
         
         l2x = xt(ind3) - xt(ind1)
         l2y = yt(ind3) - yt(ind1)
         l2z = zt(ind3) - zt(ind1)  
         
         surf1x = (l1y*l2z-l1z*l2y)
         surf1y = (l1z*l2x-l1x*l2z)
         surf1z = (l1x*l2y-l1y*l2x)
C     AC x AD
         l1x = xt(ind3) - xt(ind1)
         l1y = yt(ind3) - yt(ind1)
         l1z = zt(ind3) - zt(ind1)
            
         l2x = xt(ind4) - xt(ind1)
         l2y = yt(ind4) - yt(ind1)
         l2z = zt(ind4) - zt(ind1)  
         
         surf2x = (l1y*l2z-l1z*l2y)
         surf2y = (l1z*l2x-l1x*l2z)
         surf2z = (l1x*l2y-l1y*l2x)
            
         surfx = surf1x + surf2x
         surfy = surf1y + surf2y
         surfz = surf1z + surf2z
         surfnx(i,1) = ONE_HALF * surfx
         surfny(i,1) = ONE_HALF * surfy
         surfnz(i,1) = ONE_HALF * surfz
         surface(i,1) = ONE_HALF * 
     &        sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
      ENDDO
      END
C-----------------------------------------------------------------------------
      SUBROUTINE k6comptetrasurf(npts, nelts, cn, xt, yt, zt,
     &                           surfnx, surfny, surfnz, surface)

      IMPLICIT NONE
#include "Def/DefFortranConst.h"


C_IN
      INTEGER_E nelts
      INTEGER_E npts
      INTEGER_E cn(0:nelts-1,4)
      REAL_E xt(0:npts-1)
      REAL_E yt(0:npts-1)
      REAL_E zt(0:npts-1)
C_OUT
      REAL_E surfnx(0:nelts-1,4)
      REAL_E surfny(0:nelts-1,4)
      REAL_E surfnz(0:nelts-1,4)
      REAL_E surface(0:nelts-1,4)
C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind3, ind4
      REAL_E surf
      REAL_E surfx, surfy, surfz
      REAL_E l1x, l1y, l1z, l2x, l2y, l2z

      DO i = 0, nelts-1
         ind1 = cn(i,1)-1       ! A1
         ind2 = cn(i,2)-1       ! A2
         ind3 = cn(i,3)-1       ! A3
         ind4 = cn(i,4)-1       ! A4

C     facette A1A2A3
         l1x = xt(ind1)-xt(ind2)
         l1y = yt(ind1)-yt(ind2)
         l1z = zt(ind1)-zt(ind2)
         l2x = xt(ind3)-xt(ind2)
         l2y = yt(ind3)-yt(ind2)
         l2z = zt(ind3)-zt(ind2)
         surfx = (l1y*l2z-l1z*l2y)
         surfy = (l1z*l2x-l1x*l2z)
         surfz = (l1x*l2y-l1y*l2x)
         surf = sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
         surfnx(i,1) = ONE_HALF * surfx
         surfny(i,1) = ONE_HALF * surfy
         surfnz(i,1) = ONE_HALF * surfz
         surface(i,1) = ONE_HALF* surf
C     facette A1A2A4
         l1x = xt(ind2)-xt(ind1)
         l1y = yt(ind2)-yt(ind1)
         l1z = zt(ind2)-zt(ind1)
         l2x = xt(ind4)-xt(ind1)
         l2y = yt(ind4)-yt(ind1)
         l2z = zt(ind4)-zt(ind1)
         surfx = (l1y*l2z-l1z*l2y)
         surfy = (l1z*l2x-l1x*l2z)
         surfz = (l1x*l2y-l1y*l2x)
         surf = sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
         surfnx(i,2) = ONE_HALF * surfx
         surfny(i,2) = ONE_HALF * surfy
         surfnz(i,2) = ONE_HALF * surfz
         surface(i,2) = ONE_HALF* surf     

C     facette A2A3A4
         l1x = xt(ind3)-xt(ind2)
         l1y = yt(ind3)-yt(ind2)
         l1z = zt(ind3)-zt(ind2)
         l2x = xt(ind4)-xt(ind2)
         l2y = yt(ind4)-yt(ind2)
         l2z = zt(ind4)-zt(ind2)
         surfx = (l1y*l2z-l1z*l2y)
         surfy = (l1z*l2x-l1x*l2z)
         surfz = (l1x*l2y-l1y*l2x)
         surf = sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
         surfnx(i,3) = ONE_HALF * surfx
         surfny(i,3) = ONE_HALF * surfy
         surfnz(i,3) = ONE_HALF * surfz
         surface(i,3) = ONE_HALF* surf 
    
C     facette A1A3A4
         l1x = xt(ind4)-xt(ind1)
         l1y = yt(ind4)-yt(ind1)
         l1z = zt(ind4)-zt(ind1)
         l2x = xt(ind3)-xt(ind1)
         l2y = yt(ind3)-yt(ind1)
         l2z = zt(ind3)-zt(ind1)
         surfx = (l1y*l2z-l1z*l2y)
         surfy = (l1z*l2x-l1x*l2z)
         surfz = (l1x*l2y-l1y*l2x)
         surfnx(i,4) = ONE_HALF * surfx
         surfny(i,4) = ONE_HALF * surfy
         surfnz(i,4) = ONE_HALF * surfz
         surf = sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
         surface(i,4) = ONE_HALF* surf    
      ENDDO
      END
C-----------------------------------------------------------------------------
      SUBROUTINE k6comphexasurf(npts, nelts, cn, xt, yt, zt,
     &                          surfnx, surfny, surfnz, surface)

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E nelts
      INTEGER_E npts
      INTEGER_E cn(0:nelts-1,8) ! connect elt->noeud
      REAL_E xt(0:npts-1)
      REAL_E yt(0:npts-1)
      REAL_E zt(0:npts-1)
C_OUT
      REAL_E surfnx(0:nelts-1,6)
      REAL_E surfny(0:nelts-1,6)
      REAL_E surfnz(0:nelts-1,6)
      REAL_E surface(0:nelts-1,6) ! 6 facettes par elts
C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8
      REAL_E surfx, surfy, surfz
      REAL_E surf1x, surf1y, surf1z, surf2x, surf2y, surf2z
      REAL_E l1x, l1y, l1z, l2x, l2y, l2z

      DO i = 0, nelts-1
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         ind3 = cn(i,3)-1
         ind4 = cn(i,4)-1
         ind5 = cn(i,5)-1
         ind6 = cn(i,6)-1
         ind7 = cn(i,7)-1
         ind8 = cn(i,8)-1
C
C     premiere facette A1A2A3A4
C     A2A1 x A2A3
         l1x = xt(ind1) - xt(ind2)
         l1y = yt(ind1) - yt(ind2)
         l1z = zt(ind1) - zt(ind2)
         
         l2x = xt(ind3) - xt(ind2)
         l2y = yt(ind3) - yt(ind2)
         l2z = zt(ind3) - zt(ind2)  
         
         surf1x = (l1y*l2z-l1z*l2y)
         surf1y = (l1z*l2x-l1x*l2z)
         surf1z = (l1x*l2y-l1y*l2x)
C     A4A3 x A4A1
         l1x = xt(ind3) - xt(ind4)
         l1y = yt(ind3) - yt(ind4)
         l1z = zt(ind3) - zt(ind4)
            
         l2x = xt(ind1) - xt(ind4)
         l2y = yt(ind1) - yt(ind4)
         l2z = zt(ind1) - zt(ind4)  
         
         surf2x = (l1y*l2z-l1z*l2y)
         surf2y = (l1z*l2x-l1x*l2z)
         surf2z = (l1x*l2y-l1y*l2x)
            
         surfx = surf1x + surf2x
         surfy = surf1y + surf2y
         surfz = surf1z + surf2z

         surfnx(i,1) = ONE_HALF * surfx
         surfny(i,1) = ONE_HALF * surfy
         surfnz(i,1) = ONE_HALF * surfz
         surface(i,1) = ONE_HALF * 
     &        sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
C
C deuxieme facette A5A6A7A8
C     A5A6 x A5A7
         l1x = xt(ind6) - xt(ind5)
         l1y = yt(ind6) - yt(ind5)
         l1z = zt(ind6) - zt(ind5)
         
         l2x = xt(ind7) - xt(ind5)
         l2y = yt(ind7) - yt(ind5)
         l2z = zt(ind7) - zt(ind5)  
         
         surf1x = (l1y*l2z-l1z*l2y)
         surf1y = (l1z*l2x-l1x*l2z)
         surf1z = (l1x*l2y-l1y*l2x)
C     A5A7 x A5A8
         l1x = xt(ind7) - xt(ind5)
         l1y = yt(ind7) - yt(ind5)
         l1z = zt(ind7) - zt(ind5)
            
         l2x = xt(ind8) - xt(ind5)
         l2y = yt(ind8) - yt(ind5)
         l2z = zt(ind8) - zt(ind5)  
         
         surf2x = (l1y*l2z-l1z*l2y)
         surf2y = (l1z*l2x-l1x*l2z)
         surf2z = (l1x*l2y-l1y*l2x)
            
         surfx = surf1x + surf2x
         surfy = surf1y + surf2y
         surfz = surf1z + surf2z
         surfnx(i,2) = ONE_HALF * surfx
         surfny(i,2) = ONE_HALF * surfy
         surfnz(i,2) = ONE_HALF * surfz
         surface(i,2) = ONE_HALF * 
     &        sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
c
C troisieme facette 4158
C     A4A1 x A4A5
         l1x = xt(ind1) - xt(ind4)
         l1y = yt(ind1) - yt(ind4)
         l1z = zt(ind1) - zt(ind4)
         
         l2x = xt(ind5) - xt(ind4)
         l2y = yt(ind5) - yt(ind4)
         l2z = zt(ind5) - zt(ind4)  
         
         surf1x = (l1y*l2z-l1z*l2y)
         surf1y = (l1z*l2x-l1x*l2z)
         surf1z = (l1x*l2y-l1y*l2x)
C     A4A5 x A4A8
         l1x = xt(ind5) - xt(ind4)
         l1y = yt(ind5) - yt(ind4)
         l1z = zt(ind5) - zt(ind4)
            
         l2x = xt(ind8) - xt(ind4)
         l2y = yt(ind8) - yt(ind4)
         l2z = zt(ind8) - zt(ind4)  
         
         surf2x = (l1y*l2z-l1z*l2y)
         surf2y = (l1z*l2x-l1x*l2z)
         surf2z = (l1x*l2y-l1y*l2x)
            
         surfx = surf1x + surf2x
         surfy = surf1y + surf2y
         surfz = surf1z + surf2z
         surfnx(i,3) = ONE_HALF * surfx
         surfny(i,3) = ONE_HALF * surfy
         surfnz(i,3) = ONE_HALF * surfz
         surface(i,3) = ONE_HALF * 
     &        sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
c
C quatrieme facette A2A3A7A6
C     A2A3x A2A7
         l1x = xt(ind3) - xt(ind2)
         l1y = yt(ind3) - yt(ind2)
         l1z = zt(ind3) - zt(ind2)
         
         l2x = xt(ind7) - xt(ind2)
         l2y = yt(ind7) - yt(ind2)
         l2z = zt(ind7) - zt(ind2)  
         
         surf1x = (l1y*l2z-l1z*l2y)
         surf1y = (l1z*l2x-l1x*l2z)
         surf1z = (l1x*l2y-l1y*l2x)
C     A2A7 x A2A6
         l1x = xt(ind7) - xt(ind2)
         l1y = yt(ind7) - yt(ind2)
         l1z = zt(ind7) - zt(ind2)
            
         l2x = xt(ind6) - xt(ind2)
         l2y = yt(ind6) - yt(ind2)
         l2z = zt(ind6) - zt(ind2)  
         
         surf2x = (l1y*l2z-l1z*l2y)
         surf2y = (l1z*l2x-l1x*l2z)
         surf2z = (l1x*l2y-l1y*l2x)
            
         surfx = surf1x + surf2x
         surfy = surf1y + surf2y
         surfz = surf1z + surf2z
         surfnx(i,4) = ONE_HALF * surfx
         surfny(i,4) = ONE_HALF * surfy
         surfnz(i,4) = ONE_HALF * surfz
         surface(i,4) = ONE_HALF * 
     &        sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
c
C cinquieme facette A1A2A6A5
C     A1A2 x A1A6
         l1x = xt(ind2) - xt(ind1)
         l1y = yt(ind2) - yt(ind1)
         l1z = zt(ind2) - zt(ind1)
         
         l2x = xt(ind6) - xt(ind1)
         l2y = yt(ind6) - yt(ind1)
         l2z = zt(ind6) - zt(ind1)  
         
         surf1x = (l1y*l2z-l1z*l2y)
         surf1y = (l1z*l2x-l1x*l2z)
         surf1z = (l1x*l2y-l1y*l2x)
C     A1A6 x A1A5
         l1x = xt(ind6) - xt(ind1)
         l1y = yt(ind6) - yt(ind1)
         l1z = zt(ind6) - zt(ind1)
            
         l2x = xt(ind5) - xt(ind1)
         l2y = yt(ind5) - yt(ind1)
         l2z = zt(ind5) - zt(ind1)  
         
         surf2x = (l1y*l2z-l1z*l2y)
         surf2y = (l1z*l2x-l1x*l2z)
         surf2z = (l1x*l2y-l1y*l2x)
            
         surfx = surf1x + surf2x
         surfy = surf1y + surf2y
         surfz = surf1z + surf2z
         surfnx(i,5) = ONE_HALF * surfx
         surfny(i,5) = ONE_HALF * surfy
         surfnz(i,5) = ONE_HALF * surfz
         surface(i,5) = ONE_HALF * 
     &        sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
c    
C sixieme facette A3A4A8A7
C     A3A4 x A3A8
         l1x = xt(ind4) - xt(ind3)
         l1y = yt(ind4) - yt(ind3)
         l1z = zt(ind4) - zt(ind3)
         
         l2x = xt(ind8) - xt(ind3)
         l2y = yt(ind8) - yt(ind3)
         l2z = zt(ind8) - zt(ind3)  
         
         surf1x = (l1y*l2z-l1z*l2y)
         surf1y = (l1z*l2x-l1x*l2z)
         surf1z = (l1x*l2y-l1y*l2x)
C    A3A8 x A3A7
         l1x = xt(ind8) - xt(ind3)
         l1y = yt(ind8) - yt(ind3)
         l1z = zt(ind8) - zt(ind3)
            
         l2x = xt(ind7) - xt(ind3)
         l2y = yt(ind7) - yt(ind3)
         l2z = zt(ind7) - zt(ind3)  
         
         surf2x = (l1y*l2z-l1z*l2y)
         surf2y = (l1z*l2x-l1x*l2z)
         surf2z = (l1x*l2y-l1y*l2x)
            
         surfx = surf1x + surf2x
         surfy = surf1y + surf2y
         surfz = surf1z + surf2z

         surfnx(i,6) = ONE_HALF * surfx
         surfny(i,6) = ONE_HALF * surfy
         surfnz(i,6) = ONE_HALF * surfz
         surface(i,6) = ONE_HALF * 
     &        sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
c
      ENDDO
      END  
C-----------------------------------------------------------------------------
C mailles prismatiques
C
      SUBROUTINE k6comppentasurf(npts, nelts, cn, xt, yt, zt,
     &                           surfnx, surfny, surfnz, surface)

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E nelts
      INTEGER_E npts
      INTEGER_E cn(0:nelts-1,6)
      REAL_E xt(0:npts-1)
      REAL_E yt(0:npts-1)
      REAL_E zt(0:npts-1)
C_OUT
      REAL_E surfnx(0:nelts-1,5)
      REAL_E surfny(0:nelts-1,5)
      REAL_E surfnz(0:nelts-1,5)
      REAL_E surface(0:nelts-1,5)
C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind3, ind4, ind5, ind6
      REAL_E surf, surfx, surfy, surfz
      REAL_E surf1x, surf1y, surf1z, surf2x, surf2y, surf2z
      REAL_E l1x, l1y, l1z, l2x, l2y, l2z

      DO i = 0, nelts-1

         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         ind3 = cn(i,3)-1
         ind4 = cn(i,4)-1
         ind5 = cn(i,5)-1
         ind6 = cn(i,6)-1

C     premiere facette : triangle A1A2A3
C     A2A1 x A2A3
         l1x = xt(ind1)-xt(ind2)
         l1y = yt(ind1)-yt(ind2)
         l1z = zt(ind1)-zt(ind2)
         l2x = xt(ind3)-xt(ind2)
         l2y = yt(ind3)-yt(ind2)
         l2z = zt(ind3)-zt(ind2)
         surfx = (l1y*l2z-l1z*l2y)
         surfy = (l1z*l2x-l1x*l2z)
         surfz = (l1x*l2y-l1y*l2x)
         surf = sqrt(surfx*surfx+surfy*surfy+surfz*surfz)

         surfnx(i,1) = ONE_HALF * surfx
         surfny(i,1) = ONE_HALF * surfy
         surfnz(i,1) = ONE_HALF * surfz
         surface(i,1) = ONE_HALF* surf

C     deuxieme facette triangle A4A5A6
C     A4A5 x A4A6
         l1x = xt(ind5)-xt(ind4)
         l1y = yt(ind5)-yt(ind4)
         l1z = zt(ind5)-zt(ind4)
         l2x = xt(ind6)-xt(ind4)
         l2y = yt(ind6)-yt(ind4)
         l2z = zt(ind6)-zt(ind4)
         surfx = (l1y*l2z-l1z*l2y)
         surfy = (l1z*l2x-l1x*l2z)
         surfz = (l1x*l2y-l1y*l2x)
         surf = sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
         surfnx(i,2) = ONE_HALF * surfx
         surfny(i,2) = ONE_HALF * surfy
         surfnz(i,2) = ONE_HALF * surfz
         surface(i,2) = ONE_HALF* surf
C
C     troisieme facette quad 1254 
C     A1A2 x A1A5
         l1x = xt(ind2) - xt(ind1)
         l1y = yt(ind2) - yt(ind1)
         l1z = zt(ind2) - zt(ind1)
         
         l2x = xt(ind5) - xt(ind1)
         l2y = yt(ind5) - yt(ind1)
         l2z = zt(ind5) - zt(ind1)  
         
         surf1x = (l1y*l2z-l1z*l2y)
         surf1y = (l1z*l2x-l1x*l2z)
         surf1z = (l1x*l2y-l1y*l2x)
C     A1A5 x A1A4
         l1x = xt(ind5) - xt(ind1)
         l1y = yt(ind5) - yt(ind1)
         l1z = zt(ind5) - zt(ind1)
            
         l2x = xt(ind4) - xt(ind1)
         l2y = yt(ind4) - yt(ind1)
         l2z = zt(ind4) - zt(ind1)  
         
         surf2x = (l1y*l2z-l1z*l2y)
         surf2y = (l1z*l2x-l1x*l2z)
         surf2z = (l1x*l2y-l1y*l2x)
            
         surfx = surf1x + surf2x
         surfy = surf1y + surf2y
         surfz = surf1z + surf2z
         surfnx(i,3) = ONE_HALF * surfx
         surfny(i,3) = ONE_HALF * surfy
         surfnz(i,3) = ONE_HALF * surfz
         surface(i,3) = ONE_HALF * 
     &        sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
C
C     quatrieme facette quad 2365
C     A2A3 x A2A6
         l1x = xt(ind3) - xt(ind2)
         l1y = yt(ind3) - yt(ind2)
         l1z = zt(ind3) - zt(ind2)
         
         l2x = xt(ind6) - xt(ind2)
         l2y = yt(ind6) - yt(ind2)
         l2z = zt(ind6) - zt(ind2)  
         
         surf1x = (l1y*l2z-l1z*l2y)
         surf1y = (l1z*l2x-l1x*l2z)
         surf1z = (l1x*l2y-l1y*l2x)
C     A2A6 x A2A5
         l1x = xt(ind6) - xt(ind2)
         l1y = yt(ind6) - yt(ind2)
         l1z = zt(ind6) - zt(ind2)
            
         l2x = xt(ind5) - xt(ind2)
         l2y = yt(ind5) - yt(ind2)
         l2z = zt(ind5) - zt(ind2)  
         
         surf2x = (l1y*l2z-l1z*l2y)
         surf2y = (l1z*l2x-l1x*l2z)
         surf2z = (l1x*l2y-l1y*l2x)
            
         surfx = surf1x + surf2x
         surfy = surf1y + surf2y
         surfz = surf1z + surf2z
         surfnx(i,4) = ONE_HALF * surfx
         surfny(i,4) = ONE_HALF * surfy
         surfnz(i,4) = ONE_HALF * surfz
         surface(i,4) = ONE_HALF * 
     &        sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
C
C     cinquieme facette quad 3146 
C     A3A1 x A3A4
         l1x = xt(ind1) - xt(ind3)
         l1y = yt(ind1) - yt(ind3)
         l1z = zt(ind1) - zt(ind3)
         
         l2x = xt(ind4) - xt(ind3)
         l2y = yt(ind4) - yt(ind3)
         l2z = zt(ind4) - zt(ind3)  
         
         surf1x = (l1y*l2z-l1z*l2y)
         surf1y = (l1z*l2x-l1x*l2z)
         surf1z = (l1x*l2y-l1y*l2x)
C     A3A4 x A3A6
         l1x = xt(ind4) - xt(ind3)
         l1y = yt(ind4) - yt(ind3)
         l1z = zt(ind4) - zt(ind3)
            
         l2x = xt(ind6) - xt(ind3)
         l2y = yt(ind6) - yt(ind3)
         l2z = zt(ind6) - zt(ind3)  
         
         surf2x = (l1y*l2z-l1z*l2y)
         surf2y = (l1z*l2x-l1x*l2z)
         surf2z = (l1x*l2y-l1y*l2x)
            
         surfx = surf1x + surf2x
         surfy = surf1y + surf2y
         surfz = surf1z + surf2z
         surfnx(i,5) = ONE_HALF * surfx
         surfny(i,5) = ONE_HALF * surfy
         surfnz(i,5) = ONE_HALF * surfz
         surface(i,5) = ONE_HALF * 
     &        sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
      ENDDO
      END
C-----------------------------------------------------------------------------
C mailles pyra
C
      SUBROUTINE k6comppyrasurf(npts, nelts, cn, xt, yt, zt,
     &                          surfnx, surfny, surfnz, surface)

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E nelts
      INTEGER_E npts
      INTEGER_E cn(0:nelts-1,5)
      REAL_E xt(0:npts-1)
      REAL_E yt(0:npts-1)
      REAL_E zt(0:npts-1)
C_OUT
      REAL_E surfnx(0:nelts-1,5)
      REAL_E surfny(0:nelts-1,5)
      REAL_E surfnz(0:nelts-1,5)
      REAL_E surface(0:nelts-1,5)
C_LOCAL
      INTEGER_E i, indt1, indt2, indt3
      INTEGER_E ind1, ind2, ind3, ind4, ind5
      REAL_E surf, surfx, surfy, surfz
      REAL_E surf1x, surf1y, surf1z, surf2x, surf2y, surf2z
      REAL_E l1x, l1y, l1z, l2x, l2y, l2z

      DO i = 0, nelts-1

         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         ind3 = cn(i,3)-1
         ind4 = cn(i,4)-1
         ind5 = cn(i,5)-1
C
C     premiere facette A1A2A3A4
C     A2A1 x A2A3
         l1x = xt(ind1) - xt(ind2)
         l1y = yt(ind1) - yt(ind2)
         l1z = zt(ind1) - zt(ind2)
         
         l2x = xt(ind3) - xt(ind2)
         l2y = yt(ind3) - yt(ind2)
         l2z = zt(ind3) - zt(ind2)  
         
         surf1x = (l1y*l2z-l1z*l2y)
         surf1y = (l1z*l2x-l1x*l2z)
         surf1z = (l1x*l2y-l1y*l2x)
C     A4A3 x A4A1
         l1x = xt(ind3) - xt(ind4)
         l1y = yt(ind3) - yt(ind4)
         l1z = zt(ind3) - zt(ind4)
            
         l2x = xt(ind1) - xt(ind4)
         l2y = yt(ind1) - yt(ind4)
         l2z = zt(ind1) - zt(ind4)  
         
         surf2x = (l1y*l2z-l1z*l2y)
         surf2y = (l1z*l2x-l1x*l2z)
         surf2z = (l1x*l2y-l1y*l2x)
            
         surfx = surf1x + surf2x
         surfy = surf1y + surf2y
         surfz = surf1z + surf2z

         surfnx(i,1) = ONE_HALF * surfx
         surfny(i,1) = ONE_HALF * surfy
         surfnz(i,1) = ONE_HALF * surfz
         surface(i,1) = ONE_HALF * 
     &        sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
C
C     Deuxieme facette : triangle A1A2A5
         indt1 = ind1
         indt2 = ind2
         indt3 = ind5
         l1x = xt(indt1)-xt(indt2)
         l1y = yt(indt1)-yt(indt2)
         l1z = zt(indt1)-zt(indt2)
         l2x = xt(indt1)-xt(indt3)
         l2y = yt(indt1)-yt(indt3)
         l2z = zt(indt1)-zt(indt3)
         surfx = (l1y*l2z-l1z*l2y)
         surfy = (l1z*l2x-l1x*l2z)
         surfz = (l1x*l2y-l1y*l2x)
         surf = sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
         surfnx(i,2) = ONE_HALF * surfx
         surfny(i,2) = ONE_HALF * surfy
         surfnz(i,2) = ONE_HALF * surfz
         surface(i,2) = ONE_HALF* surf
C
C     3eme facette : triangle A2A3A5
         indt1 = ind2
         indt2 = ind3
         indt3 = ind5
         l1x = xt(indt1)-xt(indt2)
         l1y = yt(indt1)-yt(indt2)
         l1z = zt(indt1)-zt(indt2)
         l2x = xt(indt1)-xt(indt3)
         l2y = yt(indt1)-yt(indt3)
         l2z = zt(indt1)-zt(indt3)
         surfx = (l1y*l2z-l1z*l2y)
         surfy = (l1z*l2x-l1x*l2z)
         surfz = (l1x*l2y-l1y*l2x)
         surf = sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
         surfnx(i,3) = ONE_HALF * surfx
         surfny(i,3) = ONE_HALF * surfy
         surfnz(i,3) = ONE_HALF * surfz
         surface(i,3) = ONE_HALF* surf
C
C     4eme facette : triangle A3A4A5
         indt1 = ind3
         indt2 = ind4
         indt3 = ind5
         l1x = xt(indt1)-xt(indt2)
         l1y = yt(indt1)-yt(indt2)
         l1z = zt(indt1)-zt(indt2)
         l2x = xt(indt1)-xt(indt3)
         l2y = yt(indt1)-yt(indt3)
         l2z = zt(indt1)-zt(indt3)
         surfx = (l1y*l2z-l1z*l2y)
         surfy = (l1z*l2x-l1x*l2z)
         surfz = (l1x*l2y-l1y*l2x)
         surf = sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
         surfnx(i,4) = ONE_HALF * surfx
         surfny(i,4) = ONE_HALF * surfy
         surfnz(i,4) = ONE_HALF * surfz
         surface(i,4) = ONE_HALF* surf
C
C     5eme facette : triangle A4A1A5
         indt1 = ind4
         indt2 = ind1
         indt3 = ind5
         l1x = xt(indt1)-xt(indt2)
         l1y = yt(indt1)-yt(indt2)
         l1z = zt(indt1)-zt(indt2)
         l2x = xt(indt1)-xt(indt3)
         l2y = yt(indt1)-yt(indt3)
         l2z = zt(indt1)-zt(indt3)
         surfx = (l1y*l2z-l1z*l2y)
         surfy = (l1z*l2x-l1x*l2z)
         surfz = (l1x*l2y-l1y*l2x)
         surf = sqrt(surfx*surfx+surfy*surfy+surfz*surfz)
         surfnx(i,5) = ONE_HALF * surfx
         surfny(i,5) = ONE_HALF * surfy
         surfnz(i,5) = ONE_HALF * surfz
         surface(i,5) = ONE_HALF* surf
      ENDDO
      END
C  ============================================================================
C  Calcul de la longueur entre chaque sommet pour une ligne non structuree
C  ============================================================================
      SUBROUTINE k6unstructsurf1d(npts, nelts, nnodes, cn, 
     &                            xt, yt, zt,
     &                            length) 
      IMPLICIT NONE
      
# include "Def/DefFortranConst.h"
C_IN
      INTEGER_E npts            ! nb de pts dans le maillage
      INTEGER_E nelts           ! nb d elements
      INTEGER_E nnodes          ! nb de noeuds par elemt
      INTEGER_E cn(0:nelts-1,nnodes) ! connectivite elt->noeud
      REAL_E xt(0:npts-1)       ! coordonnees des pts de la grille
      REAL_E yt(0:npts-1)       ! coordonnees des pts de la grille
      REAL_E zt(0:npts-1)       ! coordonnees des pts de la grille
C_OUT
      REAL_E length(0:nelts-1) ! longueur des elts  
C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2
      REAL_E lx, ly, lz
      
      DO i = 0, nelts-1
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         lx = xt(ind1)-xt(ind2)
         ly = yt(ind1)-yt(ind2)
         lz = zt(ind1)-zt(ind2)
         length(i) = sqrt(lx*lx+ly*ly+lz*lz)
      ENDDO
      END

C ======================= KCore/Metric/CompUnstrSurfF.for ====================
