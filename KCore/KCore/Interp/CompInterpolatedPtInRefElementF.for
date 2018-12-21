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
C Routine extraite de coef_lag_3d de sAbrinA par G. Desquesnes
C Calcule les coordonnees du point à interpoler dans l element de 
C reference

C =============================================================================
      SUBROUTINE compinterpolatedptinrefelt(
     & xt, yt, zt, npts, ic, jc, kc, ni, nj, x0,y0,z0, 
     & npts_interp_1D, npts_interp_3D,   
     & x_interp, y_interp, z_interp, base, A1, B1, indxc, indxr, ipiv,
     & xx, yy, zz, erreur)

      IMPLICIT NONE
C==============================================================================
C_IN
      INTEGER_E npts            ! nb de points du maillage
      INTEGER_E ic, jc, kc      ! indices du premier noeud de la molecule d interpolation
      INTEGER_E ni, nj          ! taille du maillage dans chaque direction
      REAL_E xt(0:npts-1)       ! maillage
      REAL_E yt(0:npts-1)
      REAL_E zt(0:npts-1)
      REAL_E x0, y0, z0         ! coordonnees du pt a interpoler ds le repere original
      INTEGER_E npts_interp_1D  ! nb de pts d'interpolation par direction = ordre
      INTEGER_E npts_interp_3D  ! nb de pts d interpolation total 
C_OUT
      REAL_E xx, yy,zz          ! coordonnees du pt a interpoler ds le rep de reference
      INTEGER_E erreur          ! vaut 1 si erreur 
      INTEGER_E di,dj,dk        ! si erreur, permettant de savoir comment
                                ! modifier les coordonnees des cellules 
                                ! interpolatrices.
C_LOCAL
      INTEGER_E l 
      INTEGER_E id_col
      INTEGER_E ind, ind0, id_pt
      INTEGER_E idebug, i, j, k, ninj
      REAL_E xp, yp, zp         ! coord de x0,y0,z0 ds rep preconditionne
      REAL_E x_interp(npts_interp_3D)
      REAL_E y_interp(npts_interp_3D)
      REAL_E z_interp(npts_interp_3D)
      REAL_E base(npts_interp_1D)
      REAL_E A1(npts_interp_3D,npts_interp_3D)
      REAL_E thres
      REAL_E mi, mj, mk, m
      REAL_E B1(npts_interp_3D,3)
c.....pour le systeme lineaire
      INTEGER_E ok
      INTEGER_E indxc(npts_interp_3D), indxr(npts_interp_3d)
      INTEGER_E ipiv(npts_interp_3D)
      REAL_E xp1, yp1, zp1
      REAL_E eps
      INTEGER_E max 
C==============================================================================
      
      idebug = 0
c
      ninj = ni * nj
c
c     Parametres determinant la limite de la cellule
      max = npts_interp_1D/2 
      eps = 0.001D0
      thres = max + eps 
c
c.... a la fin de la routine, erreur indiquera si le calcul s'est 
c     bien passe
      erreur = 0
      di = 0
      dj = 0
      dk = 0
c     
      xp = x0
      yp = y0
      zp = z0
c      
      ind0 = ic-1 + (jc-1)* ni + (kc-1) * ninj
      id_pt = 1
      do k=1, npts_interp_1D 
         do j=1, npts_interp_1D 
            do i=1, npts_interp_1D 
               ind = ind0 + i-1 + (j-1) * ni + (k-1) * ninj
               x_interp(id_pt) = xt(ind)
               y_interp(id_pt) = yt(ind)
               z_interp(id_pt) = zt(ind)
C               if ( idebug .ge. 1) then
C                  write(*,*) x_interp(id_pt),y_interp(id_pt),
C     &                 z_interp(id_pt)
C               endif
               id_pt = id_pt+1
            ENDDO
         ENDDO
      ENDDO

c.... calcul des coordonnees des points d'interpolation et du point
c     a interpoler dans un repere plus stable pour la suite 
      call cp_coord_in_stable_frame_3d( npts_interp_1D, npts_interp_3D,
     &     x_interp, y_interp, z_interp,
     &     xp, yp, zp ) 

c$$$      xp = x0
c$$$      yp = y0
c$$$      zp = z0
c     

c.....Initialisation des coefficients de la matrice permettant de 
c     calculer la transformation vers le "maillage" de reference
      do id_pt= 1, npts_interp_3D
         id_col= 1
         zp1 = 1
         do k=1, npts_interp_1D 
            yp1 = 1
            do j=1, npts_interp_1D 
               xp1 = 1
               do i=1, npts_interp_1D 
                  A1(id_pt, id_col)= xp1*yp1*zp1              
                  id_col = id_col+1
                  xp1 = xp1*x_interp(id_pt)
               enddo
               yp1 = yp1*y_interp(id_pt)
            enddo
            zp1 = zp1*z_interp(id_pt)
         enddo
      enddo
      
c     
c.....Initialisation du vecteur image
c
      if ( npts_interp_1D . eq. 2) then 
         base(1) = 0.D0
         base(2) = 1.D0
      else if ( npts_interp_1D . eq. 3) then 
         base(1) = -1.D0
         base(2) = 0.D0
         base(3) = 1.D0
      else if ( npts_interp_1D . eq. 5) then 
         base(1) = -2.D0
         base(2) = -1.D0
         base(3) = 0.D0
         base(4) = 1.D0
         base(5) = 2.D0
      else 
         WRITE(*,*) 'CompInterpolatedPtInRefElement: ',
     &        'bad number of interpolation points'
         STOP
      endif 
      
      if( idebug.ge.1 ) then
        write(*,*) '  assemblage du vecteur image'
      endif

c.... premiere colonne pour les x
      id_pt = 1
      do k=1, npts_interp_1D
         do j=1, npts_interp_1D
            do i=1, npts_interp_1D
               B1(id_pt,1) = base(i)
               B1(id_pt,2) = base(j)
               B1(id_pt,3) = base(k)
C               if ( idebug .ge. 1) then 
C                  write(*,*) 'B =',B1(id_pt,1), B1(id_pt,2),B1(id_pt,3)
C               endif
               id_pt = id_pt + 1
            enddo 
         enddo
      enddo
      
c.... Resolution du systeme AX = B ou X est l'inconnu
      if( idebug.ge.1 ) then
        write(*,*) '  resolution du systeme'
      endif
      
      call kgaussj(A1, npts_interp_3D, B1, 3,
     &            ok, indxc, indxr, ipiv)
c
c.....coordonnees du point a interpoler dans le "maillage" de reference
      if ( idebug.ge.1 ) then
        write(*,*) 'calcul des coord. du pt a interpoler dans le '
     &        ,'repere de reference'
      endif
           
      if ( idebug .ge. 1 ) then
         DO l = 1, npts_interp_3D
            xx = 0.D0
            yy = 0.D0
            zz = 0.D0
            id_pt = 1
            
            DO k = 1, npts_interp_1D
               mk = (z_interp(l))**(k-1)
               DO j = 1, npts_interp_1D
                  mj = (y_interp(l))**(j-1)
                  DO i = 1, npts_interp_1D
                     mi =  (x_interp(l))**(i-1)
                     m  =  mi*mj*mk 
                     xx =  xx + B1(id_pt,1)*m
                     yy =  yy + B1(id_pt,2)*m
                     zz =  zz + B1(id_pt,3)*m
            
                     id_pt = id_pt + 1
                  enddo
               enddo
            enddo
         ENDDO
      ENDIF
c$$$      if ( idebug .ge. 1 ) then 
c$$$         open(unit=3, file="out.txt", status="unknown")
c$$$         write(unit=3,fmt=*)'    3'
c$$$         write(unit=3,fmt=*)'x                   6e18.8'          
c$$$         write(unit=3,fmt=*)'     1   5   5   5'
c$$$         write(unit=3,fmt=10)(x_interp(l),l=1,npts_interp_3D)
c$$$         write(unit=3,fmt=*)'y                   6e18.8'          
c$$$         write(unit=3,fmt=*)'     1   5   5   5'
c$$$         write(unit=3,fmt=10)(y_interp(l),l=1,npts_interp_3D)
c$$$         write(unit=3,fmt=*)'z                   6e18.8'          
c$$$         write(unit=3,fmt=*)'     1   5   5   5'
c$$$         write(unit=3,fmt=10)(z_interp(l),l=1,npts_interp_3D)
c$$$ 10      format(6e18.8)
c$$$         close(3)
c$$$      endif

      xx = 0.D0
      yy = 0.D0
      zz = 0.D0
      id_pt = 1
     
      mk = 1
      DO k = 1, npts_interp_1D
C         mk = zp**(k-1)
         mj = 1
         DO j = 1, npts_interp_1D
C            mj = yp**(j-1)
            mi = 1
            DO i = 1, npts_interp_1D
C               mi =  xp**(i-1)
               m  =  mi*mj*mk 
               xx =  xx + B1(id_pt,1)*m
               yy =  yy + B1(id_pt,2)*m
               zz =  zz + B1(id_pt,3)*m
               id_pt = id_pt + 1
               mi = mi * xp
            enddo
            mj = mj * yp
         enddo
         mk = mk * zp
      enddo

C a l ordre 2 : (xx,yy,zz) doivent etre compris entre 0 et 1
C a l ordre 3 :                                 entre -1 et 1
C a l ordre 5 :                                 entre -2 et 2  
      if ( npts_interp_1D .eq. 2 ) then ! ordre 2
         if ( xx .lt. eps .or. yy .lt. eps .or. zz .lt. eps .or.
     &        xx .ge. thres .or. yy .ge. thres .or. zz .ge. thres) then
            if ( idebug .ge. 1 ) then
               write(*,*) ' Warning: compLagrangeCoef: computation of 
     &              the interpolated point coordinates in reference frame 
     &           is wrong.'
               write(*,*) ' thres = ' ,thres
               write(*,*) 'x0 = ', x0, y0, z0
               write(*,*) 'xp = ', xp, yp, zp
               write(*,*) 'ksi =', xx, yy, zz
            endif
            erreur = 1
         endif 
      else                      ! ordres 3 et 5
         if(abs(xx).ge.thres .or. abs(yy).ge.thres .or.abs(zz).ge.thres)
     &        then
            if ( idebug .ge. 1 ) then
               write(*,*) ' Warning: compLagrangeCoef: computation of 
     &              the interpolated point coordinates in reference frame 
     &           is wrong.'
               write(*,*) ' thres = ' ,thres
               write(*,*) 'x0 = ', x0, y0, z0
               write(*,*) 'xp = ', xp, yp, zp
               write(*,*) 'ksi =', xx, yy, zz
            endif
            erreur = 1
         endif
      endif
      end
