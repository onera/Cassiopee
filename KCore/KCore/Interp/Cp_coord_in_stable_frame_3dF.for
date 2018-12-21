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

C ============================================================================
C Routine extraite de sAbrinA par G. Desquesnes
C Calcule les coordonnees des points d interpolation et du point
C a interpoler dans un repere mieux preconditionne.
C (ex: si mailles initialement allongées, effectue une contraction etc.)
C==============================================================================
      subroutine cp_coord_in_stable_frame_3d(
     &             npts_interp_1D, npts_interp_3D,
     &             x_interp, y_interp, z_interp, 
     &             xp, yp, zp )
c***********************************************************************
c
c     ACT 
c_A    Calcul les coordonnees des points d'interpolation et du point
c_A    a interpoler dans un repere stable plus stable pour la suite.
c
c     IN 
c_I     npts_interp_1D : nombre de points d'interpolation 1D 
c_I     npts_interp_3D : nombre de points d'interpolation 3D 
c
c     OUT
c
c     I/O
c_/     x_interp,
c_/     y_interp,
c_/     z_interp : tableau de taille npts_interp_1D^3=npts_interp_3D
c_/                 contenant les adresses des points d'interpolation. 
c_/     xp, yp, zp: coord du point a interpoler
c
c***********************************************************************
c
      implicit none

# include "Interp/IndInterpF.h"
# include "Def/DefFortranConst.h"
c
c-----declaration des variables de la routine
c
      INTEGER_E  npts_interp_1D
      INTEGER_E npts_interp_3D

      REAL_E x_interp(npts_interp_3D)
      REAL_E y_interp(npts_interp_3D)
      REAL_E z_interp(npts_interp_3D)
      REAL_E xp, yp, zp
c
c-----declaration des variables locales
c
      INTEGER_E i0, j0, k0
      INTEGER_E i
      INTEGER_E di, dj, dk
      INTEGER_E id, id1, id2
      INTEGER_E idebug

      REAL_E pt_O(3), u_i(3), u_j(3), u_k(3)
      REAL_E pts(3, 4)
      REAL_E xtmp, ytmp, ztmp

c-----------------------------------------------------------------------
      idebug = 0

c.... coord du point d'interpolation de reference
      i0 = (npts_interp_1D+1)/2 
      j0 = (npts_interp_1D+1)/2
      k0 = (npts_interp_1D+1)/2
c     
c.... calcul du nouveau point d'origine
c
      DO i = 1,3
         u_i(i) = 0.D0
         u_j(i) = 0.D0
         u_k(i) = 0.D0
      ENDDO
c
      id = ind_interp( npts_interp_1D, i0, j0, k0) 

      pt_O(1) = x_interp(id)
      pt_O(2) = y_interp(id)
      pt_O(3) = z_interp(id)
c
c.... calcul des nouveaux vecteurs de base
c
      id = ind_interp( npts_interp_1D, i0+1, j0, k0) 
c
c.... calcul u_i
      do dk=0, 1
         do dj=0, 1 
            id1 = ind_interp( npts_interp_1D, i0  , j0+dj, k0+dk)
            id2 = ind_interp( npts_interp_1D, i0+1, j0+dj, k0+dk)
            u_i(1) = u_i(1) + x_interp(id2)-x_interp(id1)
            u_i(2) = u_i(2) + y_interp(id2)-y_interp(id1)
            u_i(3) = u_i(3) + z_interp(id2)-z_interp(id1)
         enddo
      enddo
      DO i = 1,3
         u_i(i) = u_i(i)* ONE_FOURTH
      ENDDO

c.... calcul de u_j
      do dk=0, 1
         do di=0, 1 
            id1 = ind_interp( npts_interp_1D, i0+di, j0  , k0+dk)
            id2 = ind_interp( npts_interp_1D, i0+di, j0+1, k0+dk)
            u_j(1) = u_j(1) + x_interp(id2)-x_interp(id1)
            u_j(2) = u_j(2) + y_interp(id2)-y_interp(id1)
            u_j(3) = u_j(3) + z_interp(id2)-z_interp(id1)
        enddo
      enddo
      DO i = 1,3
         u_j(i) = u_j(i) * ONE_FOURTH
      ENDDO


c.... calcul de u_k
      do dj=0, 1
         do di=0, 1
            id1 = ind_interp( npts_interp_1D, i0+di, j0+dj, k0)
            id2 = ind_interp( npts_interp_1D, i0+di, j0+dj, k0+1)
            u_k(1) = u_k(1) + x_interp(id2)-x_interp(id1)
            u_k(2) = u_k(2) + y_interp(id2)-y_interp(id1)
            u_k(3) = u_k(3) + z_interp(id2)-z_interp(id1)
         enddo
      enddo
      DO i = 1,3
         u_k(i) = u_k(i) * ONE_FOURTH
      ENDDO
c
c.... application du changement de base
c
      DO i = 1, 3
         pts(i,1) = pt_O(i) 
         pts(i,2) = pt_O(i) + u_i(i)
         pts(i,3) = pt_O(i) + u_j(i)
         pts(i,4) = pt_O(i) + u_k(i)
      ENDDO
      
      do i=1, npts_interp_3D
c......calcul de coordonnees du point i dans la nouvelle base
         call coord_in_ref_frame_3d( pts, 
     &        x_interp(i),y_interp(i),z_interp(i),
     &        xtmp, ytmp, ztmp )
c......mise a jour des coordonnees
        x_interp(i) = xtmp
        y_interp(i) = ytmp
        z_interp(i) = ztmp        
      enddo

c.... calcul de coordonnees du point a interpoler dans la nouvelle base
      call coord_in_ref_frame_3d( pts, xp, yp, zp, xtmp, ytmp, ztmp )
      
c.... mise a jour de ses coordonnees
      xp = xtmp
      yp = ytmp
      zp = ztmp

      end
