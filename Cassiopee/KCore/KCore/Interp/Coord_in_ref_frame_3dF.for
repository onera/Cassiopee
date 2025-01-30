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
C Routine extraite de sAbrinA par G. Desquesnes
C Calcule les coordonnees d un point donn� dans le rep de reference.
C =============================================================================
      subroutine coord_in_ref_frame_3d(pts, x1, y1, z1, x2, y2, z2 )
c***********************************************************************
c
c  ACT
c_A  Determine les coord du point (x,y,z) dans le 
c_A  repere 3D pts(.,1), pts(.,2), pts(.,3), pts(.,4)
c
c  IN
c_I  pts       : tableau reel de taille 3 x 4, pts definissant le repere
c_I  x1, y1, z1: coord du point dans le repere de depart 
c 
c  OUT
c_O  x2, y2, z2: coord du point dans le nouveau repere 
c
c***********************************************************************
c
      IMPLICIT NONE

# include "Def/DefFortranConst.h"

c
c-----declaration des variables de la routine
c
      REAL_E pts(3,4), x1, y1, z1, x2, y2, z2

c
c-----declaration des variables locales
c
      INTEGER_E idebug
      INTEGER_E i
      REAL_E t(3)
      REAL_E mat(3,3)
      REAL_E invdet, u(3), v(3), w(3)

      idebug = 0

      DO i = 1, 3
         u(i) = pts(i,2) - pts(i,1)
         v(i) = pts(i,3) - pts(i,1)
         w(i) = pts(i,4) - pts(i,1)
      ENDDO

      t(1) = x1 - pts(1,1)
      t(2) = y1 - pts(2,1)
      t(3) = z1 - pts(3,1)

      if( idebug.ge.1 ) then
        write(*,*) '  u=', (u(i),i=1,3) 
        write(*,*) '  v=', (v(i),i=1,3) 
        write(*,*) '  w=', (w(i),i=1,3)
        write(*,*) '  t=', (t(i),i=1,3) 
      endif

      DO i = 1, 3
         mat(i,1) = u(i)
         mat(i,2) = v(i)
         mat(i,3) = w(i)
      ENDDO
      
      call determinant( THREE_I, mat, invdet )
      if( idebug.ge.1 ) then
        write(*,*) '  det(u,v,w)=', invdet
      endif 

      invdet = 1.D0/invdet

c.... x dans le nouveau repere
      DO i = 1, 3
         mat(i,1) = t(i) 
         mat(i,2) = v(i) 
         mat(i,3) = w(i) 
      ENDDO

      call determinant(THREE_I, mat, x2)
      if( idebug.ge.1 ) then
        write(*,*) '  det(t, v, w)=', x2
      endif
      x2 = x2 * invdet

c.... y dans le nouveau repere
      DO i = 1, 3
         mat(i,1) = u(i)
         mat(i,2) = t(i)
         mat(i,3) = w(i)
      ENDDO
      call determinant(THREE_I, mat, y2)
      if( idebug.ge.1 ) then
        write(*,*) '  det(u, t, w)=', y2
      endif
      y2 = y2 * invdet

c.... z dans le nouveau repere
      DO i = 1, 3
         mat(i,1) = u(i)
         mat(i,2) = v(i)
         mat(i,3) = t(i)
      ENDDO
      call determinant(THREE_I, mat, z2)
      if( idebug.ge.1 ) then
        write(*,*) '  det(u, v, t)=', z2
      endif
      z2 = z2 * invdet

      if( idebug.ge.1 ) then
        write(*,*) 'x2=', x2, 'y2=', y2, 'z2=', z2
      endif

      if( idebug.ge.1 ) then
        write(*,*) '--- fin coord_in_ref_frame_3d ---'
      endif

      end

c-----------------------------------------------------------------------
 
      subroutine determinant(ni, mat, det)
c***********************************************************************
c
c     G.Desquesnes
c     
c_ACT
c_A  calcul le determinant de la matrice carree mat
c     
c_IN
c_I   ni  : dimension de la matrice
c_I   mat : matrice reelle de taille ni x ni
c     
c_OUT
c_O  det : determinat de la matrice
c
c***********************************************************************
c
      implicit none
c
c---- declaration des variables de la routine
c
      INTEGER_E ni 
      REAL_E    mat(ni, ni)
      REAL_E    det

C
      INTEGER_E ierror, idebug

c
      idebug = 0

      if( idebug.ge.1 ) then
        write(*,*) '=== determinant, debug=', idebug, ' ==='
      endif

c
c.... Verification des parametres d'entree
c
      ierror = 0

      if( ni.lt.1 ) then
        write(*,9000)
        write(*,*) '  dimension de la matrice non valide'
        if( idebug.ge.1 ) then
          write(*,*) '    ni=', ni
        endif
        ierror = 1
      endif

      if( ni.gt.3 ) then
        write(*,9000)
        write(*,*) '  le cas ni>3 n''est pas encore implemete'
        if( idebug.ge.1 ) then
          write(*,*) '    ni=', ni
        endif
        ierror = 1
      endif
 
      if( ierror.ne.0 ) then
        write(*,*) 'fin du programme pour cause d''erreurs'
        stop 1
      endif

c
c.... Algo
      if (ni.eq.3) then
        det = mat(1,1)*mat(2,2)*mat(3,3)
     &       +mat(2,1)*mat(3,2)*mat(1,3)
     &       +mat(3,1)*mat(1,2)*mat(2,3)
     &       -mat(1,1)*mat(3,2)*mat(2,3)
     &       -mat(2,1)*mat(1,2)*mat(3,3)
     &       -mat(3,1)*mat(2,2)*mat(1,3) 

      elseif(ni.eq.2 ) then
         det = mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)

      elseif( ni.eq.1 ) then
         det = mat(1,1)
        
      endif

      if( idebug.ge.1 ) then
        write(*,*) '--- fin determinant ---'
      endif

c
c.... Format
c     
9000  format('erreur determinant :')

      end 
