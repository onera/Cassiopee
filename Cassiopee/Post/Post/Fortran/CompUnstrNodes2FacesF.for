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

C     Calcul des valeurs d un champ aux faces des elements a partir des 
C     des valeurs aux noeuds.
C     Depend du type d'element considere,  
C     stockage : cf KCore/Metric/CompUnstructSurf
C
C ============================================================================
      SUBROUTINE k6unstrnodes2faces(dim, npts, nelts, nedges, nnodes,  
     &                              cn, fieldn, fieldf)
      IMPLICIT NONE
      
# include "Def/DefFortranConst.h"

C==============================================================================
C_IN
      INTEGER_E dim             ! 2 ou 3 
      INTEGER_E npts            ! nb de pts dans le maillage
      INTEGER_E nelts           ! nb d elements
      INTEGER_E nedges          ! nb de facettes par elemt
      INTEGER_E nnodes          ! nb de noeuds par elemt
      INTEGER_E cn(0:nelts-1,nnodes) ! connectivite elt->noeud
      REAL_E fieldn(0:npts-1)   ! champ aux noeuds
      
C_OUT
      REAL_E fieldf(0:nelts-1,nedges) ! champs sur les facettes des elts  
      
C==============================================================================
      IF (dim .eq. 2) THEN 
         IF ( nedges .eq. 3 .and. nnodes .eq. 3 ) THEN ! TRI
            call k6comptrifield(npts, nelts, cn, fieldn, fieldf) 
            
         ELSE IF ( nedges .eq. 4 .and. nnodes .eq. 4 ) THEN ! QUAD
            call k6compquadfield(npts, nelts, cn, fieldn, fieldf) 
         ENDIF
      ELSE                      ! dim 3
         IF ( nedges .eq. 4 .and. nnodes .eq. 4 ) THEN ! TETRA
            call k6comptetrafield(npts, nelts, cn, fieldn, fieldf)
            
         ELSE IF ( nedges .eq. 6 .and. nnodes .eq. 8 ) THEN ! HEXA
            call k6comphexafield(npts, nelts, cn, fieldn, fieldf)
            
         ELSE IF ( nedges .eq. 5 .and. nnodes .eq. 6 ) THEN ! PENTA
            call k6comppentafield(npts, nelts, cn, fieldn, fieldf)
         ENDIF
      ENDIF
      END

C-----------------------------------------------------------------------------
      SUBROUTINE k6comptrifield(npts, nelts, cn, fieldn, fieldf)

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E nelts
      INTEGER_E npts
      INTEGER_E cn(0:nelts-1,3)
      REAL_E fieldn(0:npts-1)

C_OUT
      REAL_E fieldf(0:nelts-1,3)
C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind3
      REAL_E f1, f2, f3

C lecture : P1P2, P2P3, P3P1
      DO i = 0, nelts-1
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         ind3 = cn(i,3)-1
         f1 = fieldn(ind1)
         f2 = fieldn(ind2)
         f3 = fieldn(ind3)
         
         fieldf(i,1) = ONE_HALF * ( f1 + f2 ) 
         fieldf(i,2) = ONE_HALF * ( f2 + f3 ) 
         fieldf(i,3) = ONE_HALF * ( f3 + f1 ) 
      ENDDO
      END
C-----------------------------------------------------------------------------
      SUBROUTINE k6compquadfield(npts, nelts, cn, fieldn, fieldf) 

      IMPLICIT NONE
#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E nelts
      INTEGER_E npts
      INTEGER_E cn(0:nelts-1,4)
      REAL_E fieldn(0:npts-1)
C_OUT
      REAL_E fieldf(0:nelts-1,4)
C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind3, ind4
      REAL_E f1, f2, f3, f4

      DO i = 0, nelts-1
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         ind3 = cn(i,3)-1
         ind4 = cn(i,4)-1
         f1 = fieldn(ind1)
         f2 = fieldn(ind2)
         f3 = fieldn(ind3)
         f4 = fieldn(ind4)
         
         fieldf(i,1) = ONE_HALF * ( f1 + f2 ) 
         fieldf(i,2) = ONE_HALF * ( f2 + f3 ) 
         fieldf(i,3) = ONE_HALF * ( f3 + f4 ) 
         fieldf(i,3) = ONE_HALF * ( f4 + f1 ) 
 
      ENDDO
      END
C-----------------------------------------------------------------------------
      SUBROUTINE k6comptetrafield(npts, nelts, cn, fieldn, fieldf)

      IMPLICIT NONE
#include "Def/DefFortranConst.h"


C_IN
      INTEGER_E nelts
      INTEGER_E npts
      INTEGER_E cn(0:nelts-1,4)
      REAL_E fieldn(0:npts-1)

C_OUT
      REAL_E fieldf(0:nelts-1,4)
C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind3, ind4
      REAL_E f1, f2, f3, f4
      REAL_E inv

      inv = ONE/THREE

      DO i = 0, nelts-1
         ind1 = cn(i,1)-1       ! A1
         ind2 = cn(i,2)-1       ! A2
         ind3 = cn(i,3)-1       ! A3
         ind4 = cn(i,4)-1       ! A4

         f1 = fieldn(ind1)
         f2 = fieldn(ind2)
         f3 = fieldn(ind3)
         f4 = fieldn(ind4)

C     facette A1A2A3
         fieldf(i,1) = inv * ( f1 + f2 + f3)
C     facette A1A2A4
         fieldf(i,2) = inv * ( f1 + f2 + f4)
C     facette A2A3A4
         fieldf(i,3) = inv * ( f2 + f3 + f4)         
C     facette A1A3A4
         fieldf(i,4) = inv * ( f1 + f3 + f4)         
      
      ENDDO
      END
C-----------------------------------------------------------------------------
      SUBROUTINE k6comphexafield(npts, nelts, cn, fieldn, fieldf )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E nelts
      INTEGER_E npts
      INTEGER_E cn(0:nelts-1,8) ! connect elt->noeud
      REAL_E fieldn(0:npts-1)

C_OUT
      REAL_E fieldf(0:nelts-1,6) ! 6 facettes par elts

C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8
      REAL_E f1, f2, f3, f4, f5, f6, f7, f8
      REAL_E inv
      
      inv = ONE_FOURTH 

      DO i = 0, nelts-1
         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         ind3 = cn(i,3)-1
         ind4 = cn(i,4)-1
         ind5 = cn(i,5)-1
         ind6 = cn(i,6)-1
         ind7 = cn(i,7)-1
         ind8 = cn(i,8)-1

         f1 = fieldn(ind1)
         f2 = fieldn(ind2)
         f3 = fieldn(ind3)
         f4 = fieldn(ind4)
         f5 = fieldn(ind5)
         f6 = fieldn(ind6)
         f7 = fieldn(ind7)
         f8 = fieldn(ind8)
C
C     premiere facette A1A2A3A4
         fieldf(i,1) = inv * (f1+f2+f3+f4)
C
C     deuxieme facette A5A6A7A8
         fieldf(i,2) = inv * (f5+f6+f7+f8)
c
C     troisieme facette 4158
         fieldf(i,3) = inv * (f4+f1+f5+f8)
c     
C     quatrieme facette A2A3A7A6
         fieldf(i,4) = inv * (f2+f3+f7+f6)
c     
C     cinquieme facette A1A2A6A5
         fieldf(i,5) = inv * (f1+f2+f6+f5)
c     
C     sixieme facette A3A4A8A7
         fieldf(i,6) = inv * (f3+f4+f8+f7)
      ENDDO
      END  
C-----------------------------------------------------------------------------
C mailles prismatiques
C
      SUBROUTINE k6comppentafield(npts, nelts, cn, fieldn, fieldf )

      IMPLICIT NONE

#include "Def/DefFortranConst.h"

C_IN
      INTEGER_E nelts
      INTEGER_E npts
      INTEGER_E cn(0:nelts-1,6)
      REAL_E fieldn(0:npts-1)

C_OUT
      REAL_E fieldf(0:nelts-1,5)

C_LOCAL
      INTEGER_E i
      INTEGER_E ind1, ind2, ind3, ind4, ind5, ind6
      REAL_E f1, f2, f3, f4, f5, f6
      REAL_E inv3

      inv3 = ONE/THREE
      
      DO i = 0, nelts-1

         ind1 = cn(i,1)-1
         ind2 = cn(i,2)-1
         ind3 = cn(i,3)-1
         ind4 = cn(i,4)-1
         ind5 = cn(i,5)-1
         ind6 = cn(i,6)-1

         f1 = fieldn(ind1)
         f2 = fieldn(ind2)
         f3 = fieldn(ind3)
         f4 = fieldn(ind4)
         f5 = fieldn(ind5)
         f6 = fieldn(ind6)

C     premiere facette : triangle A1A2A3
         fieldf(i,1) = inv3 * (f1 + f2 + f3)

C     deuxieme facette triangle A4A5A6
         fieldf(i,2) = inv3 * (f4 + f5 + f6)
C
C     troisieme facette quad 1254 
         fieldf(i,3) = ONE_FOURTH * ( f1+f2+f5+f4)
C
C     quatrieme facette quad 2365
         fieldf(i,4) = ONE_FOURTH * ( f2+f3+f6+f5)
C
C     cinquieme facette quad 3146 
         fieldf(i,5) = ONE_FOURTH * ( f3+f1+f4+f6)

      ENDDO
      END

C ================ Post/Fortran/CompUnstrNodes2FacesF.for ==============
