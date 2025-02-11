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
C Check the quality of mesh
C  ==========================================================================
C       
C       Written by : O.P. Jacquotte and G. Coussement
C       
        SUBROUTINE k6checkmesh(ni,nj,nk, xd,yd,zd)        
C_IN
        INTEGER_E ni, nj, nk
        REAL_E xd(*), yd(*), zd(*)
C_LOCAL
        INTEGER_E nij
        PARAMETER  (noddm = 1100000)
        PARAMETER  (noddomdm = 1100000)
        PARAMETER  (ndodm = 1)
        DIMENSION x(3*noddm)
        DIMENSION vol(noddm)
        DIMENSION angl(3*noddm)
        DIMENSION xlgth(3*noddm)
        DIMENSION volm(3*noddm)
        DIMENSION work(3*noddm)
c       
        DIMENSION xi(noddomdm)
        DIMENSION ixi(noddomdm)
c       
        DIMENSION voldom(ndodm)
        DIMENSION ipnt(ndodm),imax(ndodm),jmax(ndodm),kmax(ndodm)
        
c       CALCULATION OF CELL VOLUME
        WRITE(*,3000)
        nij = ni*nj
        nbdom = 1
        nbnode = nij
        ipnt(1) = 1
        nint = -8
        imax(1) = ni
        jmax(1) = nj
        kmax(1) = 1
        iangl = 1
        
        DO i = 1, nij
           x(i) = xd(i)
           x(i+nij) = yd(i)
        ENDDO
      
        call k6jacob(nbdom,nbnode,ipnt, imax, jmax, kmax ,x, vol, 
     &               voldom, nint, work, ixi, xi)
     
c       MESH QUALITY MEASURMENTS
        WRITE(*,4000)
        call k6angl_qual(iangl,nbdom,nbnode,ipnt,imax,jmax,kmax,x,
     &                   angl,work(1),work(1+nbnode))
        WRITE(*,5000)
        call k6lgth_qual(nbdom,nbnode,ipnt,imax,jmax,kmax,x,xlgth,work)
        WRITE(*,6000)
        call k6volm_qual(nbdom,nbnode,ipnt,imax,jmax,kmax,x,vol,volm,
     &  work,xlgth,ixi,xi)
c       
 1000        format(/,1x,'Input Reading',/)
 2000        format(/,1x,'Mesh Reading',/)
 3000        format(/,1x,'Volume Control',/)
 4000        format(/,1x,'Angle Analysis',/)
 5000        format(/,1x,'Length Analysis',/)
 6000        format(/,1x,'Volume Analysis',/)
c    
        END
c
c=======================================================================
      function acosdd(angle)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c authors      : Gregory COUSSEMENT       (ONERA:OAt1)
c version     : 1.1          date        : 09/01/95 ; Gregory COUSSEMENT
c
c Copyright   : (c) 1994 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
       pi = 4.0*atan(1.0)
c
       acosdd = 180.0/pi * acos(angle)
c
       return
       end
c
c=======================================================================
      subroutine k6angl_qual(iangl,nbdom,nbnode,ipnt,imax,jmax,kmax,x
     &                      ,angl,anglmax,work)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autors      : Gregory COUSSEMENT       (ONERA:OAt1)
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT
c
c Copyright   : (c) 1992 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
      dimension work(*)
      dimension x(*)
      dimension angl(*)
      dimension anglmax(*)
      dimension ipnt(*),imax(*),jmax(*),kmax(*)
c
      do 100 n=1,nbdom
       ip       = ipnt(n)
       im       = imax(n)
       jm       = jmax(n)
       km       = kmax(n)
       im1      = max(1,im-1)
       jm1      = max(1,jm-1)
       km1      = max(1,km-1)
       ind      = 3*(ip-1) + 1
       call k6qual_angl(iangl,im,jm,km,im1,jm1,km1,x(ind),angl(ind)
     &                 ,work(1),work(1+3*8*im1*jm1))
 100  continue
c
      nbinter   = 18
      fmin      =  0.0
      fmax      = 90.0
      do 200 n=1,nbdom
       write(*,1000) n
       ip       = ipnt(n)
       im       = imax(n)
       jm       = jmax(n)
       km       = kmax(n)
c
       nt1      = 20
       nt2      = 21
       ind      = 3*(ip-1)+1+0*(im*jm*km)
       write(*,2000)
c      call histo (nbinter,fmin,fmax,angl(ind),im*jm*km,nt1,nt2)
       call k6histo2(nbinter,fmin,fmax,angl(ind),im*jm*km,nt1,nt2)
       nt1      = 22
       nt2      = 23
       ind      = 3*(ip-1)+1+1*(im*jm*km)
       write(*,3000)
c      call histo (nbinter,fmin,fmax,angl(ind),im*jm*km,nt1,nt2)
       call k6histo2(nbinter,fmin,fmax,angl(ind),im*jm*km,nt1,nt2)
       nt1      = 24
       nt2      = 25
       ind      = 3*(ip-1)+1+2*(im*jm*km)
       write(*,4000)
c      call histo (nbinter,fmin,fmax,angl(ind),im*jm*km,nt1,nt2)
       call k6histo2(nbinter,fmin,fmax,angl(ind),im*jm*km,nt1,nt2)
       nt1      = 26
       nt2      = 27
       ind      = 3*(ip-1)+1+0*(im*jm*km)
       write(*,5000)
c      call histo (nbinter,fmin,fmax,angl(ind),3*im*jm*km,nt1,nt2)
       call k6histo2(nbinter,fmin,fmax,angl(ind),3*im*jm*km,nt1,nt2)
       nbvar    = 3
       call k6maxthree(im*jm*km,nbvar,angl(ind),anglmax(ip))
       nt1      = 28
       nt2      = 29
       ind      = 3*(ip-1)+1+0*(im*jm*km)
       write(*,6000)
c      call histo (nbinter,fmin,fmax,anglmax(ind),im*jm*km,nt1,nt2)
       call k6histo2(nbinter,fmin,fmax,anglmax(ind),im*jm*km,nt1,nt2)
 200  continue
c
      nvar = 3
c     call plot(nbdom,nbnode,ipnt,imax,jmax,kmax,nvar,x,angl,'ANGL.PLT')
      nvar = 1 
c     call plot(nbdom,nbnode,ipnt,imax,jmax,kmax,nvar,x
c    &         ,anglmax,'ANGLMAX.PLT')
c
 1000 format(/,2x,'Sub. Dom.:',i5)
 2000 format(/,3x,'I Surface Angles')
 3000 format(/,3x,'J Surface Angles')
 4000 format(/,3x,'K Surface Angles')
 5000 format(/,3x,'All Surface Angles')
 6000 format(/,3x,'Maximum Angles')
c
      return
      end

c
c=======================================================================
      subroutine k6angle(im,jm,km,im1,jm1,km1,x,angl)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autors      : Gregory COUSSEMENT       (ONERA:OAt1)
c
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT
c
c Copyright   : (c) 1991 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
      dimension    x(3,8,im1*jm1)
      dimension angl(3,8,im1*jm1)
c
      small   = 1.0e-6
c
c
      do 100 ij = 1,im1*jm1
c
       x12 = x(1,2,ij) - x(1,1,ij)
       y12 = x(2,2,ij) - x(2,1,ij)
       z12 = x(3,2,ij) - x(3,1,ij)
c
       x43 = x(1,3,ij) - x(1,4,ij)
       y43 = x(2,3,ij) - x(2,4,ij)
       z43 = x(3,3,ij) - x(3,4,ij)
c
       x14 = x(1,4,ij) - x(1,1,ij)
       y14 = x(2,4,ij) - x(2,1,ij)
       z14 = x(3,4,ij) - x(3,1,ij)
c
       x23 = x(1,3,ij) - x(1,2,ij)
       y23 = x(2,3,ij) - x(2,2,ij)
       z23 = x(3,3,ij) - x(3,2,ij)
c
       x56 = x(1,6,ij) - x(1,5,ij)
       y56 = x(2,6,ij) - x(2,5,ij)
       z56 = x(3,6,ij) - x(3,5,ij)
c
       x87 = x(1,7,ij) - x(1,8,ij)
       y87 = x(2,7,ij) - x(2,8,ij)
       z87 = x(3,7,ij) - x(3,8,ij)
c
       x58 = x(1,8,ij) - x(1,5,ij)
       y58 = x(2,8,ij) - x(2,5,ij)
       z58 = x(3,8,ij) - x(3,5,ij)
c
       x67 = x(1,7,ij) - x(1,6,ij)
       y67 = x(2,7,ij) - x(2,6,ij)
       z67 = x(3,7,ij) - x(3,6,ij)
c
       x15 = x(1,5,ij) - x(1,1,ij)
       y15 = x(2,5,ij) - x(2,1,ij)
       z15 = x(3,5,ij) - x(3,1,ij)
c
       x26 = x(1,6,ij) - x(1,2,ij)
       y26 = x(2,6,ij) - x(2,2,ij)
       z26 = x(3,6,ij) - x(3,2,ij)
c
       x37 = x(1,7,ij) - x(1,3,ij)
       y37 = x(2,7,ij) - x(2,3,ij)
       z37 = x(3,7,ij) - x(3,3,ij)
c
       x48 = x(1,8,ij) - x(1,4,ij)
       y48 = x(2,8,ij) - x(2,4,ij)
       z48 = x(3,8,ij) - x(3,4,ij)
c
       xi  = x12
       yi  = y12
       zi  = z12
       xj  = x14
       yj  = y14
       zj  = z14
       xk  = x15
       yk  = y15
       zk  = z15
       xin    = sqrt(xi*xi + yi*yi + zi*zi)
       xjn    = sqrt(xj*xj + yj*yj + zj*zj)
       xkn    = sqrt(xk*xk + yk*yk + zk*zk)
       txin   = 0.5*(1.0 + sign(1.0,xin-small))
       txjn   = 0.5*(1.0 + sign(1.0,xjn-small))
       txkn   = 0.5*(1.0 + sign(1.0,xkn-small))
       txin1  = 1.0 - txin
       txjn1  = 1.0 - txjn
       txkn1  = 1.0 - txkn
       xi     = txin*xi + txin1*0.0 
       yi     = txin*yi + txin1*0.0 
       zi     = txin*zi + txin1*0.0 
       xj     = txjn*xj + txjn1*0.0 
       yj     = txjn*yj + txjn1*0.0 
       zj     = txjn*zj + txjn1*0.0 
       xk     = txkn*xk + txkn1*0.0 
       yk     = txkn*yk + txkn1*0.0 
       zk     = txkn*zk + txkn1*0.0 
       xin    = max(xin,small)
       xjn    = max(xjn,small)
       xkn    = max(xkn,small)
       xi     = xi/xin
       yi     = yi/xin
       zi     = zi/xin
       xj     = xj/xjn
       yj     = yj/xjn
       zj     = zj/xjn
       xk     = xk/xkn
       yk     = yk/xkn
       zk     = zk/xkn
       xixj   = xi*xj + yi*yj + zi*zj
       xjxk   = xj*xk + yj*yk + zj*zk
       xkxi   = xk*xi + yk*yi + zk*zi
       angl(1,1,ij)   = abs(acosdd(xjxk) - 90.0)
       angl(2,1,ij)   = abs(acosdd(xkxi) - 90.0)
       angl(3,1,ij)   = abs(acosdd(xixj) - 90.0)
c
       xi  = x12
       yi  = y12
       zi  = z12
       xj  = x23
       yj  = y23
       zj  = z23
       xk  = x26
       yk  = y26
       zk  = z26
       xin    = sqrt(xi*xi + yi*yi + zi*zi)
       xjn    = sqrt(xj*xj + yj*yj + zj*zj)
       xkn    = sqrt(xk*xk + yk*yk + zk*zk)
       txin   = 0.5*(1.0 + sign(1.0,xin-small))
       txjn   = 0.5*(1.0 + sign(1.0,xjn-small))
       txkn   = 0.5*(1.0 + sign(1.0,xkn-small))
       txin1  = 1.0 - txin
       txjn1  = 1.0 - txjn
       txkn1  = 1.0 - txkn
       xi     = txin*xi + txin1*0.0 
       yi     = txin*yi + txin1*0.0 
       zi     = txin*zi + txin1*0.0 
       xj     = txjn*xj + txjn1*0.0 
       yj     = txjn*yj + txjn1*0.0 
       zj     = txjn*zj + txjn1*0.0 
       xk     = txkn*xk + txkn1*0.0 
       yk     = txkn*yk + txkn1*0.0 
       zk     = txkn*zk + txkn1*0.0 
       xin    = max(xin,small)
       xjn    = max(xjn,small)
       xkn    = max(xkn,small)
       xi     = xi/xin
       yi     = yi/xin
       zi     = zi/xin
       xj     = xj/xjn
       yj     = yj/xjn
       zj     = zj/xjn
       xk     = xk/xkn
       yk     = yk/xkn
       zk     = zk/xkn
       xixj   = xi*xj + yi*yj + zi*zj
       xjxk   = xj*xk + yj*yk + zj*zk
       xkxi   = xk*xi + yk*yi + zk*zi
       angl(1,2,ij)   = abs(acosdd(xjxk) - 90.0)
       angl(2,2,ij)   = abs(acosdd(xkxi) - 90.0)
       angl(3,2,ij)   = abs(acosdd(xixj) - 90.0)
c
       xi  = x43
       yi  = y43
       zi  = z43
       xj  = x23
       yj  = y23
       zj  = z23
       xk  = x37
       yk  = y37
       zk  = z37
       xin    = sqrt(xi*xi + yi*yi + zi*zi)
       xjn    = sqrt(xj*xj + yj*yj + zj*zj)
       xkn    = sqrt(xk*xk + yk*yk + zk*zk)
       txin   = 0.5*(1.0 + sign(1.0,xin-small))
       txjn   = 0.5*(1.0 + sign(1.0,xjn-small))
       txkn   = 0.5*(1.0 + sign(1.0,xkn-small))
       txin1  = 1.0 - txin
       txjn1  = 1.0 - txjn
       txkn1  = 1.0 - txkn
       xi     = txin*xi + txin1*0.0 
       yi     = txin*yi + txin1*0.0 
       zi     = txin*zi + txin1*0.0 
       xj     = txjn*xj + txjn1*0.0 
       yj     = txjn*yj + txjn1*0.0 
       zj     = txjn*zj + txjn1*0.0 
       xk     = txkn*xk + txkn1*0.0 
       yk     = txkn*yk + txkn1*0.0 
       zk     = txkn*zk + txkn1*0.0 
       xin    = max(xin,small)
       xjn    = max(xjn,small)
       xkn    = max(xkn,small)
       xi     = xi/xin
       yi     = yi/xin
       zi     = zi/xin
       xj     = xj/xjn
       yj     = yj/xjn
       zj     = zj/xjn
       xk     = xk/xkn
       yk     = yk/xkn
       zk     = zk/xkn
       xixj   = xi*xj + yi*yj + zi*zj
       xjxk   = xj*xk + yj*yk + zj*zk
       xkxi   = xk*xi + yk*yi + zk*zi
       angl(1,3,ij)   = abs(acosdd(xjxk) - 90.0)
       angl(2,3,ij)   = abs(acosdd(xkxi) - 90.0)
       angl(3,3,ij)   = abs(acosdd(xixj) - 90.0)
c
       xi  = x43
       yi  = y43
       zi  = z43
       xj  = x14
       yj  = y14
       zj  = z14
       xk  = x48
       yk  = y48
       zk  = z48
       xin    = sqrt(xi*xi + yi*yi + zi*zi)
       xjn    = sqrt(xj*xj + yj*yj + zj*zj)
       xkn    = sqrt(xk*xk + yk*yk + zk*zk)
       txin   = 0.5*(1.0 + sign(1.0,xin-small))
       txjn   = 0.5*(1.0 + sign(1.0,xjn-small))
       txkn   = 0.5*(1.0 + sign(1.0,xkn-small))
       txin1  = 1.0 - txin
       txjn1  = 1.0 - txjn
       txkn1  = 1.0 - txkn
       xi     = txin*xi + txin1*0.0 
       yi     = txin*yi + txin1*0.0 
       zi     = txin*zi + txin1*0.0 
       xj     = txjn*xj + txjn1*0.0 
       yj     = txjn*yj + txjn1*0.0 
       zj     = txjn*zj + txjn1*0.0 
       xk     = txkn*xk + txkn1*0.0 
       yk     = txkn*yk + txkn1*0.0 
       zk     = txkn*zk + txkn1*0.0 
       xin    = max(xin,small)
       xjn    = max(xjn,small)
       xkn    = max(xkn,small)
       xi     = xi/xin
       yi     = yi/xin
       zi     = zi/xin
       xj     = xj/xjn
       yj     = yj/xjn
       zj     = zj/xjn
       xk     = xk/xkn
       yk     = yk/xkn
       zk     = zk/xkn
       xixj   = xi*xj + yi*yj + zi*zj
       xjxk   = xj*xk + yj*yk + zj*zk
       xkxi   = xk*xi + yk*yi + zk*zi
       angl(1,4,ij)   = abs(acosdd(xjxk) - 90.0)
       angl(2,4,ij)   = abs(acosdd(xkxi) - 90.0)
       angl(3,4,ij)   = abs(acosdd(xixj) - 90.0)
c
       xi  = x56
       yi  = y56
       zi  = z56
       xj  = x58
       yj  = y58
       zj  = z58
       xk  = x15
       yk  = y15
       zk  = z15
       xin    = sqrt(xi*xi + yi*yi + zi*zi)
       xjn    = sqrt(xj*xj + yj*yj + zj*zj)
       xkn    = sqrt(xk*xk + yk*yk + zk*zk)
       txin   = 0.5*(1.0 + sign(1.0,xin-small))
       txjn   = 0.5*(1.0 + sign(1.0,xjn-small))
       txkn   = 0.5*(1.0 + sign(1.0,xkn-small))
       txin1  = 1.0 - txin
       txjn1  = 1.0 - txjn
       txkn1  = 1.0 - txkn
       xi     = txin*xi + txin1*0.0 
       yi     = txin*yi + txin1*0.0 
       zi     = txin*zi + txin1*0.0 
       xj     = txjn*xj + txjn1*0.0 
       yj     = txjn*yj + txjn1*0.0 
       zj     = txjn*zj + txjn1*0.0 
       xk     = txkn*xk + txkn1*0.0 
       yk     = txkn*yk + txkn1*0.0 
       zk     = txkn*zk + txkn1*0.0 
       xin    = max(xin,small)
       xjn    = max(xjn,small)
       xkn    = max(xkn,small)
       xi     = xi/xin
       yi     = yi/xin
       zi     = zi/xin
       xj     = xj/xjn
       yj     = yj/xjn
       zj     = zj/xjn
       xk     = xk/xkn
       yk     = yk/xkn
       zk     = zk/xkn
       xixj   = xi*xj + yi*yj + zi*zj
       xjxk   = xj*xk + yj*yk + zj*zk
       xkxi   = xk*xi + yk*yi + zk*zi
       angl(1,5,ij)   = abs(acosdd(xjxk) - 90.0)
       angl(2,5,ij)   = abs(acosdd(xkxi) - 90.0)
       angl(3,5,ij)   = abs(acosdd(xixj) - 90.0)
c
       xi  = x56
       yi  = y56
       zi  = z56
       xj  = x67
       yj  = y67
       zj  = z67
       xk  = x26
       yk  = y26
       zk  = z26
       xin    = sqrt(xi*xi + yi*yi + zi*zi)
       xjn    = sqrt(xj*xj + yj*yj + zj*zj)
       xkn    = sqrt(xk*xk + yk*yk + zk*zk)
       txin   = 0.5*(1.0 + sign(1.0,xin-small))
       txjn   = 0.5*(1.0 + sign(1.0,xjn-small))
       txkn   = 0.5*(1.0 + sign(1.0,xkn-small))
       txin1  = 1.0 - txin
       txjn1  = 1.0 - txjn
       txkn1  = 1.0 - txkn
       xi     = txin*xi + txin1*0.0 
       yi     = txin*yi + txin1*0.0 
       zi     = txin*zi + txin1*0.0 
       xj     = txjn*xj + txjn1*0.0 
       yj     = txjn*yj + txjn1*0.0 
       zj     = txjn*zj + txjn1*0.0 
       xk     = txkn*xk + txkn1*0.0 
       yk     = txkn*yk + txkn1*0.0 
       zk     = txkn*zk + txkn1*0.0 
       xin    = max(xin,small)
       xjn    = max(xjn,small)
       xkn    = max(xkn,small)
       xi     = xi/xin
       yi     = yi/xin
       zi     = zi/xin
       xj     = xj/xjn
       yj     = yj/xjn
       zj     = zj/xjn
       xk     = xk/xkn
       yk     = yk/xkn
       zk     = zk/xkn
       xixj   = xi*xj + yi*yj + zi*zj
       xjxk   = xj*xk + yj*yk + zj*zk
       xkxi   = xk*xi + yk*yi + zk*zi
       angl(1,6,ij)   = abs(acosdd(xjxk) - 90.0)
       angl(2,6,ij)   = abs(acosdd(xkxi) - 90.0)
       angl(3,6,ij)   = abs(acosdd(xixj) - 90.0)
c
       xi  = x87
       yi  = y87
       zi  = z87
       xj  = x67
       yj  = y67
       zj  = z67
       xk  = x37
       yk  = y37
       zk  = z37
       xin    = sqrt(xi*xi + yi*yi + zi*zi)
       xjn    = sqrt(xj*xj + yj*yj + zj*zj)
       xkn    = sqrt(xk*xk + yk*yk + zk*zk)
       txin   = 0.5*(1.0 + sign(1.0,xin-small))
       txjn   = 0.5*(1.0 + sign(1.0,xjn-small))
       txkn   = 0.5*(1.0 + sign(1.0,xkn-small))
       txin1  = 1.0 - txin
       txjn1  = 1.0 - txjn
       txkn1  = 1.0 - txkn
       xi     = txin*xi + txin1*0.0 
       yi     = txin*yi + txin1*0.0 
       zi     = txin*zi + txin1*0.0 
       xj     = txjn*xj + txjn1*0.0 
       yj     = txjn*yj + txjn1*0.0 
       zj     = txjn*zj + txjn1*0.0 
       xk     = txkn*xk + txkn1*0.0 
       yk     = txkn*yk + txkn1*0.0 
       zk     = txkn*zk + txkn1*0.0 
       xin    = max(xin,small)
       xjn    = max(xjn,small)
       xkn    = max(xkn,small)
       xi     = xi/xin
       yi     = yi/xin
       zi     = zi/xin
       xj     = xj/xjn
       yj     = yj/xjn
       zj     = zj/xjn
       xk     = xk/xkn
       yk     = yk/xkn
       zk     = zk/xkn
       xixj   = xi*xj + yi*yj + zi*zj
       xjxk   = xj*xk + yj*yk + zj*zk
       xkxi   = xk*xi + yk*yi + zk*zi
       angl(1,7,ij)   = abs(acosdd(xjxk) - 90.0)
       angl(2,7,ij)   = abs(acosdd(xkxi) - 90.0)
       angl(3,7,ij)   = abs(acosdd(xixj) - 90.0)
c
       xi  = x87
       yi  = y87
       zi  = z87
       xj  = x58
       yj  = y58
       zj  = z58
       xk  = x48
       yk  = y48
       zk  = z48
       xin    = sqrt(xi*xi + yi*yi + zi*zi)
       xjn    = sqrt(xj*xj + yj*yj + zj*zj)
       xkn    = sqrt(xk*xk + yk*yk + zk*zk)
       txin   = 0.5*(1.0 + sign(1.0,xin-small))
       txjn   = 0.5*(1.0 + sign(1.0,xjn-small))
       txkn   = 0.5*(1.0 + sign(1.0,xkn-small))
       txin1  = 1.0 - txin
       txjn1  = 1.0 - txjn
       txkn1  = 1.0 - txkn
       xi     = txin*xi + txin1*0.0 
       yi     = txin*yi + txin1*0.0 
       zi     = txin*zi + txin1*0.0 
       xj     = txjn*xj + txjn1*0.0 
       yj     = txjn*yj + txjn1*0.0 
       zj     = txjn*zj + txjn1*0.0 
       xk     = txkn*xk + txkn1*0.0 
       yk     = txkn*yk + txkn1*0.0 
       zk     = txkn*zk + txkn1*0.0 
       xin    = max(xin,small)
       xjn    = max(xjn,small)
       xkn    = max(xkn,small)
       xi     = xi/xin
       yi     = yi/xin
       zi     = zi/xin
       xj     = xj/xjn
       yj     = yj/xjn
       zj     = zj/xjn
       xk     = xk/xkn
       yk     = yk/xkn
       zk     = zk/xkn
       xixj   = xi*xj + yi*yj + zi*zj
       xjxk   = xj*xk + yj*yk + zj*zk
       xkxi   = xk*xi + yk*yi + zk*zi
       angl(1,8,ij)   = abs(acosdd(xjxk) - 90.0)
       angl(2,8,ij)   = abs(acosdd(xkxi) - 90.0)
       angl(3,8,ij)   = abs(acosdd(xixj) - 90.0)
c
 100  continue
c
      return
      end

c
c=======================================================================
      subroutine k6assemb2(im,jm,km,im1,jm1,km1,x,xx,k)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autor       : Olivier-Pierre JACQUOTTE (ONERA:OAt1)
c version     : 1.1          date        : 01/01/91 ; Olivier-Pierre JACQUOTTE
c          revised           date        : 01/09/91 ; Gregory COUSSEMENT
c
c Copyright   : (c) 1991 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
c
      dimension  x(im,jm,km,3)
      dimension xx(3,8,im1,jm1)
c
      do 100 j = 1,jm1
      do 100 i = 1,im1
c
       i1      = min(i+1,im)
       j1      = min(j+1,jm)
       k1      = min(k+1,km)
c
       xx(1,1,i,j) = x(i ,j ,k ,1)
       xx(1,2,i,j) = x(i1,j ,k ,1)
       xx(1,3,i,j) = x(i1,j1,k ,1)
       xx(1,4,i,j) = x(i ,j1,k ,1)
       xx(1,5,i,j) = x(i ,j ,k1,1)
       xx(1,6,i,j) = x(i1,j ,k1,1)
       xx(1,7,i,j) = x(i1,j1,k1,1)
       xx(1,8,i,j) = x(i ,j1,k1,1)
c
       xx(2,1,i,j) = x(i ,j ,k ,2)
       xx(2,2,i,j) = x(i1,j ,k ,2)
       xx(2,3,i,j) = x(i1,j1,k ,2)
       xx(2,4,i,j) = x(i ,j1,k ,2)
       xx(2,5,i,j) = x(i ,j ,k1,2)
       xx(2,6,i,j) = x(i1,j ,k1,2)
       xx(2,7,i,j) = x(i1,j1,k1,2)
       xx(2,8,i,j) = x(i ,j1,k1,2)
c
       xx(3,1,i,j) = x(i ,j ,k ,3)
       xx(3,2,i,j) = x(i1,j ,k ,3)
       xx(3,3,i,j) = x(i1,j1,k ,3)
       xx(3,4,i,j) = x(i ,j1,k ,3)
       xx(3,5,i,j) = x(i ,j ,k1,3)
       xx(3,6,i,j) = x(i1,j ,k1,3)
       xx(3,7,i,j) = x(i1,j1,k1,3)
       xx(3,8,i,j) = x(i ,j1,k1,3)
c
 100  continue
c
      return
      end

c
c=======================================================================
      subroutine k6assemb3(im,jm,km,im1,jm1,km1,x,xx,k)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autor       : Olivier-Pierre JACQUOTTE (ONERA:OAt1)
c version     : 1.1          date        : 01/01/91 ; Olivier-Pierre JACQUOTTE
c          revised           date        : 01/09/91 ; Gregory COUSSEMENT
c
c Copyright   : (c) 1991 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
c
      dimension  x(im,jm,km,3)
      dimension xx(3,8,im1,jm1)
c
      do 100 j = 1,jm1
      do 100 i = 1,im1
c
       i1      = min(i+1,im)
       j1      = min(j+1,jm)
       k1      = min(k+1,km)
c
       x(i1,j ,k ,1) = x(i1,j ,k ,1) + xx(1,2,i,j)
       x(i1,j1,k ,1) = x(i1,j1,k ,1) + xx(1,3,i,j)
       x(i1,j ,k1,1) = x(i1,j ,k1,1) + xx(1,6,i,j)
       x(i1,j1,k1,1) = x(i1,j1,k1,1) + xx(1,7,i,j)
       x(i ,j ,k ,1) = x(i ,j ,k ,1) + xx(1,1,i,j)
       x(i ,j1,k ,1) = x(i ,j1,k ,1) + xx(1,4,i,j)
       x(i ,j ,k1,1) = x(i ,j ,k1,1) + xx(1,5,i,j)
       x(i ,j1,k1,1) = x(i ,j1,k1,1) + xx(1,8,i,j)
c
       x(i1,j ,k ,2) = x(i1,j ,k ,2) + xx(2,2,i,j)
       x(i1,j1,k ,2) = x(i1,j1,k ,2) + xx(2,3,i,j)
       x(i1,j ,k1,2) = x(i1,j ,k1,2) + xx(2,6,i,j)
       x(i1,j1,k1,2) = x(i1,j1,k1,2) + xx(2,7,i,j)
       x(i ,j ,k ,2) = x(i ,j ,k ,2) + xx(2,1,i,j)
       x(i ,j1,k ,2) = x(i ,j1,k ,2) + xx(2,4,i,j)
       x(i ,j ,k1,2) = x(i ,j ,k1,2) + xx(2,5,i,j)
       x(i ,j1,k1,2) = x(i ,j1,k1,2) + xx(2,8,i,j)
c
       x(i1,j ,k ,3) = x(i1,j ,k ,3) + xx(3,2,i,j)
       x(i1,j1,k ,3) = x(i1,j1,k ,3) + xx(3,3,i,j)
       x(i1,j ,k1,3) = x(i1,j ,k1,3) + xx(3,6,i,j)
       x(i1,j1,k1,3) = x(i1,j1,k1,3) + xx(3,7,i,j)
       x(i ,j ,k ,3) = x(i ,j ,k ,3) + xx(3,1,i,j)
       x(i ,j1,k ,3) = x(i ,j1,k ,3) + xx(3,4,i,j)
       x(i ,j ,k1,3) = x(i ,j ,k1,3) + xx(3,5,i,j)
       x(i ,j1,k1,3) = x(i ,j1,k1,3) + xx(3,8,i,j)
c
 100  continue
c
      return
      end

c
c=======================================================================
      subroutine k6assemb5(im,jm,km,im1,jm1,km1,x,xx,k)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autor       : Gregory COUSSEMENT (ONERA:OAt1)
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT
c
c Copyright   : (c) 1992 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
c
      dimension  x(im,jm,km,3)
      dimension xx(3,8,im1,jm1)
c
      do 100 j = 1,jm1
      do 100 i = 1,im1
c
       i1      = min(i+1,im)
       j1      = min(j+1,jm)
       k1      = min(k+1,km)
c
       x(i1,j ,k ,1) = max(x(i1,j ,k ,1),xx(1,2,i,j))
       x(i1,j1,k ,1) = max(x(i1,j1,k ,1),xx(1,3,i,j))
       x(i1,j ,k1,1) = max(x(i1,j ,k1,1),xx(1,6,i,j))
       x(i1,j1,k1,1) = max(x(i1,j1,k1,1),xx(1,7,i,j))
       x(i ,j ,k ,1) = max(x(i ,j ,k ,1),xx(1,1,i,j))
       x(i ,j1,k ,1) = max(x(i ,j1,k ,1),xx(1,4,i,j))
       x(i ,j ,k1,1) = max(x(i ,j ,k1,1),xx(1,5,i,j))
       x(i ,j1,k1,1) = max(x(i ,j1,k1,1),xx(1,8,i,j))
c
       x(i1,j ,k ,2) = max(x(i1,j ,k ,2),xx(2,2,i,j))
       x(i1,j1,k ,2) = max(x(i1,j1,k ,2),xx(2,3,i,j))
       x(i1,j ,k1,2) = max(x(i1,j ,k1,2),xx(2,6,i,j))
       x(i1,j1,k1,2) = max(x(i1,j1,k1,2),xx(2,7,i,j))
       x(i ,j ,k ,2) = max(x(i ,j ,k ,2),xx(2,1,i,j))
       x(i ,j1,k ,2) = max(x(i ,j1,k ,2),xx(2,4,i,j))
       x(i ,j ,k1,2) = max(x(i ,j ,k1,2),xx(2,5,i,j))
       x(i ,j1,k1,2) = max(x(i ,j1,k1,2),xx(2,8,i,j))
c
       x(i1,j ,k ,3) = max(x(i1,j ,k ,3),xx(3,2,i,j))
       x(i1,j1,k ,3) = max(x(i1,j1,k ,3),xx(3,3,i,j))
       x(i1,j ,k1,3) = max(x(i1,j ,k1,3),xx(3,6,i,j))
       x(i1,j1,k1,3) = max(x(i1,j1,k1,3),xx(3,7,i,j))
       x(i ,j ,k ,3) = max(x(i ,j ,k ,3),xx(3,1,i,j))
       x(i ,j1,k ,3) = max(x(i ,j1,k ,3),xx(3,4,i,j))
       x(i ,j ,k1,3) = max(x(i ,j ,k1,3),xx(3,5,i,j))
       x(i ,j1,k1,3) = max(x(i ,j1,k1,3),xx(3,8,i,j))
c
 100  continue
c
      return
      end
c
c=======================================================================
      subroutine k6correct(im,jm,km,angl)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autors      : Gregory COUSSEMENT       (ONERA:OAt1)
c
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT
c
c Copyright   : (c) 1991 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
      dimension angl(im,jm,km,3)
c
      do 100 k = 1,km
      do 100 j = 1,jm
      do 100 i = 1,im
       angl(i,j,k,1)  = angl(i,j,k,1)*0.125
       angl(i,j,k,2)  = angl(i,j,k,2)*0.125
       angl(i,j,k,3)  = angl(i,j,k,3)*0.125
 100  continue
      do 200 k = 1,km
      do 200 j = 1,jm
       angl(1 ,j,k,1) = angl(1 ,j,k,1)*2.0
       angl(1 ,j,k,2) = angl(1 ,j,k,2)*2.0
       angl(1 ,j,k,3) = angl(1 ,j,k,3)*2.0
       angl(im,j,k,1) = angl(im,j,k,1)*2.0
       angl(im,j,k,2) = angl(im,j,k,2)*2.0
       angl(im,j,k,3) = angl(im,j,k,3)*2.0
 200  continue
      do 300 k = 1,km
      do 300 i = 1,im
       angl(i,1 ,k,1) = angl(i,1 ,k,1)*2.0
       angl(i,1 ,k,2) = angl(i,1 ,k,2)*2.0
       angl(i,1 ,k,3) = angl(i,1 ,k,3)*2.0
       angl(i,jm,k,1) = angl(i,jm,k,1)*2.0
       angl(i,jm,k,2) = angl(i,jm,k,2)*2.0
       angl(i,jm,k,3) = angl(i,jm,k,3)*2.0
 300  continue
      do 400 j = 1,jm
      do 400 i = 1,im
       angl(i,j,1 ,1) = angl(i,j,1 ,1)*2.0
       angl(i,j,1 ,2) = angl(i,j,1 ,2)*2.0
       angl(i,j,1 ,3) = angl(i,j,1 ,3)*2.0
       angl(i,j,km,1) = angl(i,j,km,1)*2.0
       angl(i,j,km,2) = angl(i,j,km,2)*2.0
       angl(i,j,km,3) = angl(i,j,km,3)*2.0
 400  continue
c
      return
      end

c
c=======================================================================
      subroutine k6histo(nbinter,fmin,fmax,f,imax,nt1,nt2)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autor       : Gregory COUSSEMENT       (ONERA:OAt1)
c
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT 
c
c Copyright   : (c) 1991 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
        dimension ifreq(1000)
        dimension freq(1000)
        dimension f(imax)
c
        small         = 1.0e-06
c
        nmax          = nbinter+1
c
        do 100 index = 1,nmax
         ifreq(index) = 0
 100    continue
c
        do 200 i = 1,imax
         index        = (f(i) - fmin)/(fmax - fmin)*(nmax-1)+1.5
         index        = max(min(nmax,index),1)
         ifreq(index) = ifreq(index) + 1
 200    continue
c
        do 300 index = 2,nmax
         freq(index)  = fmin + (fmax-fmin)*(index-1.5)/(nmax-1)
 300    continue
        freq(1     )  = fmin
        freq(nmax+1)  = fmax
c
        write(*,999) imax
        do 400 index = 1,nmax
         write(*,1000) freq(index),freq(index+1)
     &                 ,(1.0*ifreq(index))/(1.0*imax),ifreq(index)
 400    continue
c
        write(nt1,2000) 2*nmax
        write(nt2,2000) 2*nmax
        do 500 index = 1,nmax 
         write(nt1,3000) freq(index)
         write(nt1,3000) freq(index+1)
         write(nt2,3000) (1.0*ifreq(index))/(1.0*imax)
         write(nt2,3000) (1.0*ifreq(index))/(1.0*imax)
 500    continue
c
 999    format(
     &     /,3x,'          FREQUENCY HISTOGRAM          Tot.Nb.Nodes :',
     &  i6,/,3x,'Bound Inf.   Bound Sup.        Freq.       Nb.Nodes')
 1000   format(3(1x,f12.6),9x,i6)
 1001   format(2(1x,f12.6),i6)
 2000   format(1x,i6)
 3000   format(1x,e15.7)
c
        return
        end

c
c=======================================================================
      subroutine k6histo2(nbinter,fmin,fmax,f,imax,nt1,nt2)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autor       : Gregory COUSSEMENT       (ONERA:OAt1)
c
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT 
c
c Copyright   : (c) 1991 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
        dimension ifreq(1000)
        dimension freq(1000)
        dimension f(imax)
c
        small         = 1.0e-06
c
        nmax          = nbinter + 1
c
        do 100 index = 1,nmax-1
         ifreq(index) = 0
 100    continue
c
        do 200 i = 1,imax
         index = (f(i) - fmin)/(fmax - fmin)*(nmax-1) + 1
         index = max(min(nmax-1,index),1)
         ifreq(index) = ifreq(index) + 1
 200    continue
c
        do 300 index = 1,nmax
         freq(index)  = fmin + (fmax-fmin)*(index-1)/(nmax-1)
 300    continue
c
        write(*,999) imax
        do 400 index = 1,nmax-1
         write(*,1000) freq(index),freq(index+1)
     &                 ,(1.0*ifreq(index))/(1.0*imax),ifreq(index)
 400    continue
c
 999    format(
     &     /,3x,'          FREQUENCY HISTOGRAM          Tot.Nb.Nodes :',
     &  i6,/,3x,'Bound Inf.   Bound Sup.        Freq.       Nb.Nodes')
 1000   format(3(1x,f12.6),9x,i6)
 1001   format(2(1x,f12.6),i6)
 2000   format(1x,i5)
 3000   format(1x,e15.7)
c
        return
        end

c
c=======================================================================
      subroutine k6jacob(nbdom,nbnode,ipnt,imax,jmax,kmax
     &                  ,x,vol,voldom,nint,work,ixi,xi)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autors      : Gregory COUSSEMENT       (ONERA:OAt1)
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT
c
c Copyright   : (c) 1992 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
c    DECLARATIONS
      dimension x(*)
      dimension work(*)
      dimension vol(*)
      dimension voldom(*)
      dimension ipnt(*),imax(*),jmax(*),kmax(*)
      dimension xi(*)
      dimension ixi(*)
c
      voltot    = 0.0
      do 10 i=1,nbnode
       vol(i)   = 0.0
 10   continue
c
      nbcell    = 0
c
      do 100 n=1,nbdom
       ip       = ipnt(n)
       im       = imax(n)
       jm       = jmax(n)
       km       = kmax(n)
       im1      = max(1,im-1)
       jm1      = max(1,jm-1)
       km1      = max(1,km-1)
       nbcell   = nbcell + im1*jm1*km1
       ind      = 3*(ip-1)+1
c
       call k6volum(n,im,jm,km,im1,jm1,km1,x(ind),vol(ip),voldom(n),neg
     &             ,nint,work(1),work(im))
c
       write(*,1000) n,neg,voldom(n)
       voltot    = voltot+voldom(n)
 100  continue
      write(*,2000) voltot
c
      nbinter   = 40 
      fmin      =  0.0
      fmax      = voltot/nbcell
      do 200 n=1,nbdom
       write(*,3000) n
       ip       = ipnt(n)
       im       = imax(n)
       jm       = jmax(n)
       km       = kmax(n)
       im1      = max(1,im-1)
       jm1      = max(1,jm-1)
       km1      = max(1,km-1)
       nt1      = 20
       nt2      = 21
       nmax     = im1*jm1*km1
       write(*,4000)
c      call histo (nbinter,fmin,fmax,vol(ip),nmax,nt1,nt2)
       call k6histo2(nbinter,fmin,fmax,vol(ip),nmax,nt1,nt2)
 200  continue
c
      do 300 n=1,nbdom
       ip       = ipnt(n)
       im       = imax(n)
       jm       = jmax(n)
       km       = kmax(n)
       lm       = 1
       im1      = max(1,im-1)
       jm1      = max(1,jm-1)
       km1      = max(1,km-1)
       call k6c_v3d(im,jm,km,ixi,xi)
       call k6lint3d(work(ip),im ,jm ,km ,lm
     &              ,vol(ip) ,im1,jm1,km1,lm,ixi,xi)
 300  continue
c
      nvar = 1 
c     call plot(nbdom,nbnode,ipnt,imax,jmax,kmax,nvar,x,work,'VOL.PLT')
c
 1000 format(/,' SubDom.:',i5,' Neg.Cells :',i5,'  SubDom.Vol.:',e15.7)
 2000 format(/5x,'Volume of the domain  VOL  = ',e9.3,/)
 3000 format(/,2x,'Sub. Dom.:',i5)
 4000 format(/,3x,'All Volumes')
c
      return
      end

c
c=======================================================================
      subroutine k6lecdat(nbdom,nbnode,ipnt,imax,jmax,kmax,nint,iangl
     &                   ,ndodm,noddm)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autor       : Gregory COUSSEMENT (ONERA:OAt1)
c
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT
c
c Copyright   : (c) 1992 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
c    DECLARATIONS
      character*20 texte
      dimension ipnt(*),imax(*),jmax(*),kmax(*)
c
      open(10,file='DATA',status='UNKNOWN')
c
      read(10,'(a)') texte
c
c     READ NUMBER OF SUB-DOMAINS AND DIMENSIONS
c
      nbnode     = 1
      read (10,*)    nbdom
      do 100 n=1,nbdom
       read (10,1000) imax(n),jmax(n),kmax(n)
       ipnt(n)=nbnode
       nbnode=nbnode+imax(n)*jmax(n)*kmax(n)
 100  continue
      nbnode = nbnode-1                                  ! TOTAL NUMBER OF NODES
      read(10,2000) nint
      read(10,2000) iangl
c
c     TEST THE MESH SIZE AND STOP IF NOT ENOUGH MEMORY
c
      if (nbnode.gt.noddm) then
       write(*,3000)
       stop
      end if
c
 1000 format(3i5)
 2000 format(i5)
 3000 format(2x,'Too many nodes in the mesh compare'
     &         ,' to dimensioning parameters')
c
      return
      end

c
c=======================================================================
      subroutine k6lgth_qual(nbdom,nbnode,ipnt,imax,jmax,kmax,x
     &                      ,xlgth,xlgthmax)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autors      : Gregory COUSSEMENT       (ONERA:OAt1)
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT
c
c Copyright   : (c) 1992 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
      dimension x(*)
      dimension xlgth(*)
      dimension xlgthmax(*)
      dimension ipnt(*),imax(*),jmax(*),kmax(*)
c
      do 100 n=1,nbdom
       ip       = ipnt(n)
       im       = imax(n)
       jm       = jmax(n)
       km       = kmax(n)
       ind      = 3*(ip-1)+1
       call k6qual_lgth(im,jm,km,x(ind),xlgth(ind))
 100  continue
c
      nbinter   = 40 
      fmin      =  1.0
      fmax      =  5.0
      do 200 n=1,nbdom
       write(*,1000) n
       ip       = ipnt(n)
       im       = imax(n)
       jm       = jmax(n)
       km       = kmax(n)
c
       nt1      = 20
       nt2      = 21
       ind      = 3*(ip-1)+1+0*(im*jm*km)
       write(*,2000)
c      call histo (nbinter,fmin,fmax,xlgth(ind),im*jm*km,nt1,nt2)
       call k6histo2(nbinter,fmin,fmax,xlgth(ind),im*jm*km,nt1,nt2)
       nt1      = 22
       nt2      = 23
       ind      = 3*(ip-1)+1+1*(im*jm*km)
       write(*,3000)
c      call histo (nbinter,fmin,fmax,xlgth(ind),im*jm*km,nt1,nt2)
       call k6histo2(nbinter,fmin,fmax,xlgth(ind),im*jm*km,nt1,nt2)
       nt1      = 24
       nt2      = 25
       ind      = 3*(ip-1)+1+2*(im*jm*km)
       write(*,4000)
c      call histo (nbinter,fmin,fmax,xlgth(ind),im*jm*km,nt1,nt2)
       call k6histo2(nbinter,fmin,fmax,xlgth(ind),im*jm*km,nt1,nt2)
       nt1      = 26
       nt2      = 27
       ind      = 3*(ip-1)+1+0*(im*jm*km)
       write(*,5000)
c      call histo (nbinter,fmin,fmax,xlgth(ind),3*im*jm*km,nt1,nt2)
       call k6histo2(nbinter,fmin,fmax,xlgth(ind),3*im*jm*km,nt1,nt2)
       nbvar    = 3
       call k6maxthree(im*jm*km,nbvar,xlgth(ind),xlgthmax(ip))
       nt1      = 28
       nt2      = 29
       ind      = 3*(ip-1)+1+0*(im*jm*km)
       write(*,6000)
c      call histo (nbinter,fmin,fmax,xlgthmax(ind),im*jm*km,nt1,nt2)
       call k6histo2(nbinter,fmin,fmax,xlgthmax(ind),im*jm*km,nt1,nt2)
 200  continue
c
      nvar = 3
c     call plot(nbdom,nbnode,ipnt,imax,jmax,kmax,nvar,x
c    &         ,xlgth,'LGTH.PLT')
      nvar = 1 
c     call plot(nbdom,nbnode,ipnt,imax,jmax,kmax,nvar,x
c    &         ,xlgthmax,'LGTHMAX.PLT')
c
 1000 format(/,2x,'Sub. Dom.:',i5)
 2000 format(/,3x,'I Surface Lengths')
 3000 format(/,3x,'J Surface Lengths')
 4000 format(/,3x,'K Surface Lengths')
 5000 format(/,3x,'All Surface Lengths')
 6000 format(/,3x,'Maximum Lengths')
c
      return
      end

c
c=======================================================================
      subroutine k6maxthree(imax,nbvar,var,varmax)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autors      : Gregory COUSSEMENT       (ONERA:OAt1)
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT
c
c Copyright   : (c) 1992 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
      dimension var(imax,nbvar)
      dimension varmax(imax)
c
      big = 9999.
c
      do 100 i = 1,imax
       varmax(i) = -big
 100  continue
c
      do 200 nvar= 1,nbvar
      do 200 i = 1,imax
       varmax(i) = max(varmax(i),var(i,nvar))
 200  continue
c
      return
      end
c
c=======================================================================
      subroutine k6meshin(nbdom,nbnode,ipnt,imax,jmax,kmax,x)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autor       : Olivier-Pierre JACQUOTTE (ONERA:OAt1)
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT
c
c Copyright   : (c) 1991 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
c
c    DECLARATIONS
      dimension x(*)
      dimension ipnt(*),imax(*),jmax(*),kmax(*)
      character*20 fm,zzz
c
      open(unit=15,file='MESH0',iostat=iocheck
     &    ,status='old',access='sequential')
      if (iocheck.ne.0) then
       write(*,999) 'MESH0'
       stop
      end if
       read(15,'(i5)') nfield
       if (nfield.ne.3*nbdom+1) then
        write(* ,1000)
        stop
       end if
       read(15,'(a/a/a)') zzz,zzz,tit
c
       do 100 n=1,nbdom
        read(15,'(a,a/4i6)') zzz,fm,n0,im,jm,km
        if ((im.ne.imax(n)).or.(jm.ne.jmax(n)).or.(km.ne.kmax(n))) then
         write(* ,2000)
         stop
        end if
        istart = 3*(ipnt(n)-1)+1+0*(im*jm*km)
        iend   = 3*(ipnt(n)-1)+1+1*(im*jm*km)-1
        read(15,'('//fm//')') (x(ind),ind=istart,iend)
c
        read(15,'(a,a/4i6)') zzz,fm,n0,im,jm,km
        istart = 3*(ipnt(n)-1)+1+1*(im*jm*km)
        iend   = 3*(ipnt(n)-1)+1+2*(im*jm*km)-1
        read(15,'('//fm//')') (x(ind),ind=istart,iend)
c
        read(15,'(a,a/4i6)') zzz,fm,n0,im,jm,km
        istart = 3*(ipnt(n)-1)+1+2*(im*jm*km)
        iend   = 3*(ipnt(n)-1)+1+3*(im*jm*km)-1
        read(15,'('//fm//')') (x(ind),ind=istart,iend)
 100   continue
c
      close (unit=15)
c
  999 format(1x,'ERROR : FILE ',A20,' DOES NOT EXIST !')
 1000 format(2x,'ERROR : In mesh reading not enough field') 
 2000 format(2x,'ERROR : In mesh reading, number of grid points 
     &not equal to the one declared in the DATA file') 
c
      return
      end

c
c=======================================================================
      subroutine k6minthree(imax,nbvar,var,varmin)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autors      : Gregory COUSSEMENT       (ONERA:OAt1)
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT
c
c Copyright   : (c) 1992 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
      dimension var(imax,nbvar)
      dimension varmin(imax)
c
      big = 9999.
c
      do 100 i = 1,imax
       varmin(i) = big
 100  continue
c
      do 200 nvar= 1,nbvar
      do 200 i = 1,imax
       varmin(i) = min(varmin(i),var(i,nvar))
 200  continue
c
      return
      end
 
c
c=======================================================================   
      subroutine k6plot(nbdom,nbnode,ipnt,imax,jmax,kmax,nvar,x,var
     &                 ,string)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autor       : Gregory COUSSEMENT       (ONERA:OAt1)
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT 
c
c Copyright   : (c) 1992 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
        dimension     x(*)
        dimension     var(*)
        dimension     ipnt(*),imax(*),jmax(*),kmax(*)
        dimension     nb(10)
        character*20  ch(3)
        character*20  va(6)
        data          ch/'x','y','z'/
        data          va/'va','va','va','va','va','va'/
        data          nb/1,2,3,4,5,6,7,8,9,10/          
        character*20  string
c
        nf = 40
        open(unit=nf,file=string,form='formatted',status='UNKNOWN')
c       write(nf,'(i5)') nbdom*(3+nvar)
        write(nf,'(i5)') nbdom*nvar
        do 100 n = 1,nbdom
c
         ip      = ipnt(n)
         im      = imax(n)
         jm      = jmax(n)
         km      = kmax(n)
c        do 110 m = 1,3
c         istart  = 3*(ip-1) + 1 + (m-1)*im*jm*km
c         iend    = 3*(ip-1) + 1 +   m  *im*jm*km - 1
c         write(nf,'(2a20)') ch(m),'5e16.8'
c         write(nf,'(5i6)') n,im,jm,km
c         write(nf,'(5e16.8)') (x(i),i=istart,iend)
c110         continue
c
         do 120 m = 1,nvar
          istart  = nvar*(ip-1) + 1 + (m-1)*im*jm*km
          iend    = nvar*(ip-1) + 1 +   m  *im*jm*km - 1
          write(nf,'(2a20)') va(m),'5e16.8'
          write(nf,'(5i6)') nb(m),n,im,jm,km
          write(nf,'(6e12.5)') (var(i),i=istart,iend)
 120         continue
c
 100        continue        
c
        close(unit=nf)
c
        return
        end

c
c=======================================================================
      subroutine k6qual_angl(iangl,im,jm,km,im1,jm1,km1,x,angl,xx,angll)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autors      : Gregory COUSSEMENT       (ONERA:OAt1)
c
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT
c
c Copyright   : (c) 1991 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
      dimension x(3*im*jm*km)
      dimension angl(3*im*jm*km)
      dimension xx(*)
      dimension angll(*)
c
      do 200 k=1,km1
c
       call k6assemb2(im,jm,km,im1,jm1,km1,x,xx,k)
       call k6angle  (im,jm,km,im1,jm1,km1,xx,angll)
       if (iangl.eq.0) then
        call k6assemb3(im,jm,km,im1,jm1,km1,angl,angll,k)
       else
        call k6assemb5(im,jm,km,im1,jm1,km1,angl,angll,k)
       end if
c
 200  continue
c
       if (iangl.eq.0) then
        call k6correct(im,jm,km,angl)
       end if
c
      return
      end

c
c=======================================================================
      subroutine k6qual_lgth(im,jm,km,x,xlgth)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autor       : Gregory COUSSEMENT       (ONERA:OAt1)
c
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT 
c
c Copyright   : (c) 1991 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
      dimension x(im,jm,km,3)
      dimension xlgth(im,jm,km,3)
c
        small          = 1.0e-6
c
         do 100 k= 1,km
         do 100 j= 1,jm
         do 100 i= 1,im
          xlgth(i,j,k,1) = 1.0
          xlgth(i,j,k,2) = 1.0
          xlgth(i,j,k,3) = 1.0
 100         continue
c
c
        if (im.ge.3) then
         do 200 k= 1,km
         do 200 j= 1,jm
         do 200 i= 1,im
          im1d         = min(1,iabs(i-1 ))
          ip1d         = min(1,iabs(i-im))
          i0           = i - im1d + ip1d
          im1          = i0 - 1 
          ip1          = i0 + 1 
          dxl          = x(i0,j,k,1) - x(im1,j,k,1)
          dyl          = x(i0,j,k,2) - x(im1,j,k,2)
          dzl          = x(i0,j,k,3) - x(im1,j,k,3)
          dsl          = sqrt(dxl*dxl + dyl*dyl + dzl*dzl)
          dsl          = max(dsl,small)
          dxr          = x(ip1,j,k,1) - x(i0,j,k,1)
          dyr          = x(ip1,j,k,2) - x(i0,j,k,2)
          dzr          = x(ip1,j,k,3) - x(i0,j,k,3)
          dsr          = sqrt(dxr*dxr + dyr*dyr + dzr*dzr)
          dsr          = max(dsr,small)
          xlgth(i,j,k,1)= max(dsl/dsr,dsr/dsl)
 200         continue
        end if
c
        if (jm.ge.3) then
         do 300 k= 1,km
         do 300 j= 1,jm
         do 300 i= 1,im
          jm1d         = min(1,iabs(j-1 ))
          jp1d         = min(1,iabs(j-jm))
          j0           = j - jm1d + jp1d
          jm1          = j0 - 1 
          jp1          = j0 + 1 
          dxl          = x(i,j0,k,1) - x(i,jm1,k,1)
          dyl          = x(i,j0,k,2) - x(i,jm1,k,2)
          dzl          = x(i,j0,k,3) - x(i,jm1,k,3)
          dsl          = sqrt(dxl*dxl + dyl*dyl + dzl*dzl)
          dsl          = max(dsl,small)
          dxr          = x(i,jp1,k,1) - x(i,j0,k,1)
          dyr          = x(i,jp1,k,2) - x(i,j0,k,2)
          dzr          = x(i,jp1,k,3) - x(i,j0,k,3)
          dsr          = sqrt(dxr*dxr + dyr*dyr + dzr*dzr)
          dsr          = max(dsr,small)
          xlgth(i,j,k,2)= max(dsl/dsr,dsr/dsl)
 300         continue
        end if
c
        if (km.ge.3) then
         do 400 k= 1,km
         do 400 j= 1,jm
         do 400 i= 1,im
          km1d         = min(1,iabs(k-1 ))
          kp1d         = min(1,iabs(k-km))
          k0           = k - km1d + kp1d
          km1          = k0 - 1 
          kp1          = k0 + 1 
          dxl          = x(i,j,k0,1) - x(i,j,km1,1)
          dyl          = x(i,j,k0,2) - x(i,j,km1,2)
          dzl          = x(i,j,k0,3) - x(i,j,km1,3)
          dsl          = sqrt(dxl*dxl + dyl*dyl + dzl*dzl)
          dsl          = max(dsl,small)
          dxr          = x(i,j,kp1,1) - x(i,j,k0,1)
          dyr          = x(i,j,kp1,2) - x(i,j,k0,2)
          dzr          = x(i,j,kp1,3) - x(i,j,k0,3)
          dsr          = sqrt(dxr*dxr + dyr*dyr + dzr*dzr)
          dsr          = max(dsr,small)
          xlgth(i,j,k,3)= max(dsl/dsr,dsr/dsl)
 400         continue
        end if
c
        return
        end


c
c=======================================================================
      subroutine k6qual_volm(im,jm,km,im1,jm1,km1,vol,volmi,volmj,volmk)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autor       : Gregory COUSSEMENT       (ONERA:OAt1)
c
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT 
c
c Copyright   : (c) 1991 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
      dimension vol(im1,jm1,km1)
      dimension volmi(im,jm1,km1)
      dimension volmj(im1,jm,km1)
      dimension volmk(im1,jm1,km)
c
        small          = 1.0e-6
c
        if (im.ge.3) then
         do 200 k= 1,km1
         do 200 j= 1,jm1
         do 200 i= 1,im
          i0           = min(im1,max(2,i))
          i0m1         = i0-1
          volr         = max(vol(i0  ,j,k),small)
          voll         = max(vol(i0m1,j,k),small)
          volmi(i,j,k) = max(volr/voll,voll/volr)
 200         continue
        else
         do 210 k= 1,km1
         do 210 j= 1,jm1
         do 210 i= 1,im
          volmi(i,j,k) = 1.0
 210         continue
        end if
c
        if (jm.ge.3) then
         do 300 k= 1,km1
         do 300 j= 1,jm
         do 300 i= 1,im1
          j0           = min(jm1,max(2,j))
          j0m1         = j0-1
          volr         = max(vol(i,j0  ,k),small)
          voll         = max(vol(i,j0m1,k),small)
          volmj(i,j,k) = max(volr/voll,voll/volr)
 300         continue
        else
         do 310 k= 1,km1
         do 310 j= 1,jm
         do 310 i= 1,im1
          volmj(i,j,k) = 1.0
 310         continue
        end if
c
        if (km.ge.3) then
         do 400 k= 1,km
         do 400 j= 1,jm1
         do 400 i= 1,im1
          k0           = min(km1,max(2,k))
          k0m1         = k0-1
          volr         = max(vol(i,j,k0  ),small)
          voll         = max(vol(i,j,k0m1),small)
          volmk(i,j,k) = max(volr/voll,voll/volr)
 400         continue
        else
         do 410 k= 1,km
         do 410 j= 1,jm1
         do 410 i= 1,im1
          volmk(i,j,k) = 1.0
 410         continue
        end if
c
        return
        end

c
c=======================================================================
      subroutine k6volm_qual(nbdom,nbnode,ipnt,imax,jmax,kmax,x,vol
     &                      ,volm,volmmax,wrk,ixi,xi)
c=======================================================================
c
c ----------------------------------------------------------------------
c
c autors      : Gregory COUSSEMENT       (ONERA:OAt1)
c version     : 1.1          date        : 01/09/92 ; Gregory COUSSEMENT
c
c Copyright   : (c) 1992 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
      dimension vol(*)
      dimension x(*)
      dimension volm(*)
      dimension volmmax(*)
      dimension ipnt(*),imax(*),jmax(*),kmax(*)
      dimension xi(*)
      dimension wrk(*)
      dimension ixi(*)
c
      do 100 n=1,nbdom
       ip       = ipnt(n)
       im       = imax(n)
       jm       = jmax(n)
       km       = kmax(n)
       im1      = max(1,im-1)
       jm1      = max(1,jm-1)
       km1      = max(1,km-1)
       id       = 3*(ip-1) + 1
       id1      = id
       id2      = id1 + im*jm1*km1
       id3      = id2 + im1*jm*km1
       call k6qual_volm(im,jm,km,im1,jm1,km1,vol(id)
     &                 ,wrk(id1),wrk(id2),wrk(id3))
 100  continue
c
      nbinter   = 40 
      fmin      =  1.0
      fmax      =  5.0
      do 200 n=1,nbdom
       write(*,1000) n
       ip       = ipnt(n)
       im       = imax(n)
       jm       = jmax(n)
       km       = kmax(n)
       im1      = max(1,im-1)
       jm1      = max(1,jm-1)
       km1      = max(1,km-1)
       id       = 3*(ip-1) + 1
       id1      = id
       id2      = id1 + im*jm1*km1
       id3      = id2 + im1*jm*km1
       nt1      = 20
       nt2      = 21
       nmax     = im*jm1*km1
       write(*,2000)
c      call histo (nbinter,fmin,fmax,wrk(id1),nmax,nt1,nt2)
       call k6histo2(nbinter,fmin,fmax,wrk(id1),nmax,nt1,nt2)
       nt1      = 22
       nt2      = 23
       nmax     = im1*jm*km1
       write(*,3000)
c      call histo (nbinter,fmin,fmax,wrk(id2),nmax,nt1,nt2)
       call k6histo2(nbinter,fmin,fmax,wrk(id2),nmax,nt1,nt2)
       nt1      = 24
       nt2      = 25
       nmax     = im1*jm1*km
       write(*,4000)
c      call histo (nbinter,fmin,fmax,wrk(id3),nmax,nt1,nt2)
       call k6histo2(nbinter,fmin,fmax,wrk(id3),nmax,nt1,nt2)
       nt1      = 26
       nt2      = 27
       nmax     = im*jm1*km1+im1*jm*km1+im1*jm1*km
       write(*,5000)
c      call histo (nbinter,fmin,fmax,wrk(id1),nmax,nt1,nt2)
       call k6histo2(nbinter,fmin,fmax,wrk(id1),nmax,nt1,nt2)
 200  continue
c
      do 300 n=1,nbdom
       ip       = ipnt(n)
       im       = imax(n)
       jm       = jmax(n)
       km       = kmax(n)
       lm       = 1
       im1      = max(1,im-1)
       jm1      = max(1,jm-1)
       km1      = max(1,km-1)
       id       = 3*(ip-1) + 1
       id1      = id
       id2      = id1 + im*jm1*km1
       id3      = id2 + im1*jm*km1
       jd1      = id
       jd2      = jd1 + im*jm*km
       jd3      = jd2 + im*jm*km
       call k6ci_v3d(im,jm,km,ixi,xi)
       call k6lint3d(volm(jd1),im ,jm ,km ,lm
     &            , wrk(id1),im ,jm1,km1,lm,ixi,xi)
       call k6cj_v3d(im,jm,km,ixi,xi)
       call k6lint3d(volm(jd2),im ,jm ,km ,lm
     &            , wrk(id2),im1,jm ,km1,lm,ixi,xi)
       call k6ck_v3d(im,jm,km,ixi,xi)
       call k6lint3d(volm(jd3),im ,jm ,km ,lm
     &            , wrk(id3),im1,jm1,km ,lm,ixi,xi)
       nbvar     = 3
       call k6maxthree(im*jm*km,nbvar,volm(id),volmmax(ip))
 300  continue
c
      nvar = 3
c     call plot(nbdom,nbnode,ipnt,imax,jmax,kmax,nvar,x,volm,'VOLM.PLT')
      nvar = 1 
c     call plot(nbdom,nbnode,ipnt,imax,jmax,kmax,nvar,x
c    &         ,volmmax,'VOLMMAX.PLT')
c
 1000 format(/,2x,'Sub. Dom.:',i5)
 2000 format(/,3x,'I Surface Volumes')
 3000 format(/,3x,'J Surface Volumes')
 4000 format(/,3x,'K Surface Volumes')
 5000 format(/,3x,'All Surface Volumes')
 6000 format(/,3x,'Maximum Volumes')
c
      return
      end


c
c ======================================================================
        subroutine k6volum(ndom,im,jm,km,im1,jm1,km1
     &                    ,x,vol,voldom,neg,nint,dvol,voli)
c ======================================================================
c
c       SUBROUTINE FOR CALCULATION OF CELL VOLUME
c
c ----------------------------------------------------------------------
c
c AUTOR       : GREGORY COUSSEMENT    (ONERA:OAt1)
c VERSION     : 1.1
c DATE        : 01/09/92
c
c COPYRIGHT   : (c) 1992 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
c Entry: - x(3,im,jm,km), the coordinates of the grid points mesh where the 
c          the field "phi" value is given.
c        - im,jm,km, the dimension of the mesh x in each directions.
c        - phi(im,jm,km), the field "phi" at the grid points of the mesh "x"
c        - nint, is the number of integration point used to calulate
c          the averaged derivative.
c Exit : - dudx(im,jm,km,3), the average spatial derivative on the cell
c          I,J,K formed by the vertecies i,j,k  i+1,j,k  i+1,j+1,k  i,j+1,k
c          i,j,k+1  i+1,j,k+1 i+1,j+1,k+1 and i,j+1,k+1.
c
c  ----------------------------------------------------------------------
c
c Goal : Calculate for each cell the average spatial derivative of the 
c        field "phi" given at the point i,j,k of the mesh "x".
c        The method is based on finite element theory with three-linear
c        variations and Gauss point for integration.
c
c ----------------------------------------------------------------------
c
c
        dimension    x     (im,jm,km,3)
        dimension    vol   (im1,jm1,km1)
        dimension    dvol  (im1)
        dimension    voli  (im1)
c
        dimension np(8),xsi(15),eta(15),zet(15),r(2)
        data np/1,4,5,8,3,6,7,2/
c
        eps = 1.e-10
c
        neg = 0
c
        im1 = max(1,im-1)
        jm1 = max(1,jm-1)
        km1 = max(1,km-1)
c
c     INITIALIZATION OF INTEGRATION POINT POSITION
c
        inint       = nint
        nnint       = abs(nint)
c
                        int1=1
        if (nnint.eq.1) int1=9
        if (nnint.eq.4) int1=2
        if (nnint.eq.6) int1=10
                        int2=int1+nnint-1
c
        dr          = 1.0/sqrt(3.0)
        if (inint.eq.-6) dr=1.0
        if (inint.eq.-8) dr=1.0
        r0          = 0.0
        rm          = r0 - dr
        rp          = r0 + dr
        r(1)        = rm
        r(2)        = rp
        n           = 1
        do 100 k=1,2
        do 100 j=1,2
        do 100 i=1,2
         xsi(np(n)) = r(i)
         eta(np(n)) = r(j)
         zet(np(n)) = r(k)
         n          = n+1
 100    continue
c
        xsi(9)      = r0
        eta(9)      = r0
        zet(9)      = r0
c
        xsi(10)     = rm 
        eta(10)     = r0 
        zet(10)     = r0 
        xsi(11)     = rp 
        eta(11)     = r0 
        zet(11)     = r0 
        xsi(12)     = r0 
        eta(12)     = rm 
        zet(12)     = r0 
        xsi(13)     = r0 
        eta(13)     = rp 
        zet(13)     = r0 
        xsi(14)     = r0 
        eta(14)     = r0 
        zet(14)     = rm 
        xsi(15)     = r0 
        eta(15)     = r0 
        zet(15)     = rp 
c
        coefi  = 0.5*float(1 + isign(1,1-im))
        coefj  = 0.5*float(1 + isign(1,1-jm))
        coefk  = 0.5*float(1 + isign(1,1-km))
        coefim = 1.0 - coefi
        coefjm = 1.0 - coefj
        coefkm = 1.0 - coefk
c
        voldom = 0.0
        do 200 k=1,km1
        do 200 j=1,jm1
         do 210 i=1,im1
          dvol(i)   = 0.0
          voli(i)   = 0.0
 210     continue
         do 220 int=int1,int2
          u         = xsi(int)
          v         = eta(int)
          w         = zet(int)
          up        = 1.0 + u
          um        = 1.0 - u
          vp        = 1.0 + v 
          vm        = 1.0 - v 
          wp        = 1.0 + w 
          wm        = 1.0 - w 
c 
          do 221 i=1,im1
c
           ip1      = min(i+1,im) 
           jp1      = min(j+1,jm) 
           kp1      = min(k+1,km) 
c 
c        VALUE OF DERIVATIVE AT THE COURRANT POINT
c
           dxd1     = 0.125*
     &        (wm*(vm*(-x(i  ,j  ,k  ,1)+x(ip1,j  ,k  ,1))
     &            +vp*(-x(i  ,jp1,k  ,1)+x(ip1,jp1,k  ,1)))
     &        +wp*(vm*(-x(i  ,j  ,kp1,1)+x(ip1,j  ,kp1,1))
     &            +vp*(-x(i  ,jp1,kp1,1)+x(ip1,jp1,kp1,1)))) 
           dxd2     = 0.125*
     &        (wm*(um*(-x(i  ,j  ,k  ,1)+x(i  ,jp1,k  ,1))
     &            +up*(-x(ip1,j  ,k  ,1)+x(ip1,jp1,k  ,1)))
     &        +wp*(um*(-x(i  ,j  ,kp1,1)+x(i  ,jp1,kp1,1))
     &            +up*(-x(ip1,j  ,kp1,1)+x(ip1,jp1,kp1,1))))
           dxd3     = 0.125*
     &        (vm*(um*(-x(i  ,j  ,k  ,1)+x(i  ,j  ,kp1,1))
     &            +up*(-x(ip1,j  ,k  ,1)+x(ip1,j  ,kp1,1)))
     &        +vp*(um*(-x(i  ,jp1,k  ,1)+x(i  ,jp1,kp1,1))
     &            +up*(-x(ip1,jp1,k  ,1)+x(ip1,jp1,kp1,1))))
           dyd1     = 0.125*
     &        (wm*(vm*(-x(i  ,j  ,k  ,2)+x(ip1,j  ,k  ,2))
     &            +vp*(-x(i  ,jp1,k  ,2)+x(ip1,jp1,k  ,2)))
     &        +wp*(vm*(-x(i  ,j  ,kp1,2)+x(ip1,j  ,kp1,2))
     &            +vp*(-x(i  ,jp1,kp1,2)+x(ip1,jp1,kp1,2)))) 
           dyd2     = 0.125*
     &        (wm*(um*(-x(i  ,j  ,k  ,2)+x(i  ,jp1,k  ,2))
     &            +up*(-x(ip1,j  ,k  ,2)+x(ip1,jp1,k  ,2)))
     &        +wp*(um*(-x(i  ,j  ,kp1,2)+x(i  ,jp1,kp1,2))
     &            +up*(-x(ip1,j  ,kp1,2)+x(ip1,jp1,kp1,2))))
           dyd3     = 0.125*
     &        (vm*(um*(-x(i  ,j  ,k  ,2)+x(i  ,j  ,kp1,2))
     &            +up*(-x(ip1,j  ,k  ,2)+x(ip1,j  ,kp1,2)))
     &        +vp*(um*(-x(i  ,jp1,k  ,2)+x(i  ,jp1,kp1,2))
     &            +up*(-x(ip1,jp1,k  ,2)+x(ip1,jp1,kp1,2))))
           dzd1     = 0.125*
     &        (wm*(vm*(-x(i  ,j  ,k  ,3)+x(ip1,j  ,k  ,3))
     &            +vp*(-x(i  ,jp1,k  ,3)+x(ip1,jp1,k  ,3)))
     &        +wp*(vm*(-x(i  ,j  ,kp1,3)+x(ip1,j  ,kp1,3))
     &            +vp*(-x(i  ,jp1,kp1,3)+x(ip1,jp1,kp1,3)))) 
           dzd2     = 0.125*
     &        (wm*(um*(-x(i  ,j  ,k  ,3)+x(i  ,jp1,k  ,3))
     &            +up*(-x(ip1,j  ,k  ,3)+x(ip1,jp1,k  ,3)))
     &        +wp*(um*(-x(i  ,j  ,kp1,3)+x(i  ,jp1,kp1,3))
     &            +up*(-x(ip1,j  ,kp1,3)+x(ip1,jp1,kp1,3))))
           dzd3     = 0.125*
     &        (vm*(um*(-x(i  ,j  ,k  ,3)+x(i  ,j  ,kp1,3))
     &            +up*(-x(ip1,j  ,k  ,3)+x(ip1,j  ,kp1,3)))
     &        +vp*(um*(-x(i  ,jp1,k  ,3)+x(i  ,jp1,kp1,3))
     &            +up*(-x(ip1,jp1,k  ,3)+x(ip1,jp1,kp1,3))))
c
           dxd1     = coefim*dxd1 + coefi*0.125
           dyd1     = coefim*dyd1 + coefi*0.0
           dzd1     = coefim*dzd1 + coefi*0.0
           dxd2     = coefjm*dxd2 + coefj*0.0
           dyd2     = coefjm*dyd2 + coefj*0.125
           dzd2     = coefjm*dzd2 + coefj*0.0
           dxd3     = coefkm*dxd3 + coefk*0.0
           dyd3     = coefkm*dyd3 + coefk*0.0
           dzd3     = coefkm*dzd3 + coefk*0.125
c
c        VALUE OF THE DETERMINANT
c
             delta    = dxd1*(dyd2*dzd3 - dzd2*dyd3) 
     &                    - dxd2*(dyd1*dzd3 - dzd1*dyd3) 
     &                    + dxd3*(dyd1*dzd2 - dzd1*dyd2)
c
c        ELEMENTARY VOLUME AT THE CURRENT POINT
c
           dvol(i)    = delta
           voli(i)    = voli(i)    + delta
 221      continue
          do 222 i=1,im1
           if (dvol(i).lt.-eps) then
            write(*,1000) i,j,k,ndom,xsi(int),eta(int),zet(int)
           end if
 222      continue
 220     continue
         do 230 i=1,im1
           vol(i,j,k) = 8.0*voli(i)/nnint
 230     continue
         do 240 i=1,im1
           voldom = voldom + vol(i,j,k)
           if (vol(i,j,k).lt.-eps) then
            write(*,2000) i,j,k,ndom
            neg = neg+1
           end if
 240     continue
 200    continue
c
 1000   format(1x,'WARNING : Cell ',3(i4,'*'),'of domain ',i4
     &  ,' have a local negative volume near Gauss point',3f6.3)
 2000   format(1x,'WARNING : Cell ',3(i4,'*'),'of domain ',i4
     &  ,' have a global negative volume')
c
        return
        end

c
c ======================================================================
        subroutine k6lint3d(phi,im,jm,km,lm,phi0,im0,jm0,km0,lm0,ixi,xi)
c ======================================================================
c
c       SUBROUTINE FOR INTERPOLATION PROCESS BASED ON A TRI-LINEAR
c       INTERPOLATION.
c
c ----------------------------------------------------------------------
c
c AUTOR       : GREGORY COUSSEMENT    (ONERA:OAt1)
c VERSION     : 1.1
c DATE        : 01/11/91
c
c REFERENCES  : GREGORY COUSSEMENT AND OLIVER-PIERRE JACQUOTTE
c               "Structured Mesh Adaption: Space Accuracy and
c                Interpolation Methods.", Proceeding of the 2nd Workshop
c                on Reliability and Adaptive Methods in Computational Mechanics.
c                Krakow, POLAND, October 14-16, 1991.
c
c COPYRIGHT   : (c) 1991 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
c Entry: - im,jm,km, the dimensions of the mesh x where the field phi is
c          to be interpolated.
c        - lm=lm0, the number of field value phi to interpolate.
c        - x0(im0,jm0,km0,lm0), the field values given on the point
c          i0,j0,k0 on the initial mesh x0.
c        - im0,jm0,km0, the dimensions of the initial mesh x0 in
c          each directions.
c        - ixi(3,im,jm,km), the indices in the mesh x0 where the point
c          i,j,k, of the mesh x, is contained as the ones given by
c          routine NEWTSEARCH
c        - xi(3,im,jm,km), the interpolating coefficients for the point
c          i,j,k of the mesh x as the ones given by routine NEWTSEARCH.
c Exit : - phi(im,jm,km,lm) the interpolated field value.
c
c ----------------------------------------------------------------------
c
c Aim  : This routine interpolates the value phi at point x,y,z knowing 
c        the cell containing this point or if not the nearest cell 
c       and the tri-linear interpolation coefficients on a initial gird.
c
c ----------------------------------------------------------------------
c
        dimension    phi0(    im0,jm0,km0,lm0)
        dimension     phi(    im ,jm ,km ,lm )
        dimension      xi(  3,im ,jm ,km )
        dimension     ixi(  3,im ,jm ,km )
c
c        find phi at x,y,z by tri-linear interpolation
c
        do 100 l  = 1,lm
        do 100 k  = 1,km
        do 100 j  = 1,jm
        do 100 i  = 1,im
         ic    = ixi(1,i,j,k)
         jc    = ixi(2,i,j,k)
         kc    = ixi(3,i,j,k)
         icp1  = min(ic+1,im0)
         jcp1  = min(jc+1,jm0)
         kcp1  = min(kc+1,km0)
         u     = xi(1,i,j,k)
         v     = xi(2,i,j,k) 
         w     = xi(3,i,j,k)
         up    = 1.0 + u
         um    = 1.0 - u
         vp    = 1.0 + v 
         vm    = 1.0 - v 
         wp    = 1.0 + w 
         wm    = 1.0 - w 
         phi(i,j,k,l)  = 0.125*
     &    (wm*( vm*(um*phi0(ic,jc  ,kc  ,l)+up*phi0(icp1,jc  ,kc  ,l))
     &        + vp*(um*phi0(ic,jcp1,kc  ,l)+up*phi0(icp1,jcp1,kc  ,l)))
     &    +wp*( vm*(um*phi0(ic,jc  ,kcp1,l)+up*phi0(icp1,jc  ,kcp1,l))
     &        + vp*(um*phi0(ic,jcp1,kcp1,l)+up*phi0(icp1,jcp1,kcp1,l))))
c
 100        continue
c
        return
        end

c
c ======================================================================
        subroutine k6ci_v3d(im,jm,km,ixi,xi)
c ======================================================================
c
c       SUBROUTINE FOR FACE CENTER ( ALONG I PLANE) 
c       TO CELL VERTEX INTERPOLATION USED
c       IN A TRI-LINEAR INTERPOLATION PROCESS
c
c ----------------------------------------------------------------------
c
c AUTOR       : GREGORY COUSSEMENT    (ONERA:OAt1)
c VERSION     : 1.1
c DATE        : 01/02/92
c
c REFERENCES  : GREGORY COUSSEMENT AND OLIVER-PIERRE JACQUOTTE
c               "Structured Mesh Adaption: Space Accuracy and
c                Interpolation Methods.", Proceeding of the 2nd Workshop
c                on Reliability and Adaptive Methods in Computational Mechanics.
c                Krakow, POLAND, October 14-16, 1991.
c
c COPYRIGHT   : (c) 1991 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
c Entry: - im,jm,km, the dimension of the mesh x in each directions.
c Exit : - ixi(3,im,jm,km), the indices in the mesh x where the vertex point
c          i,j,k, of the mesh x, is contained in the mesh constituated
c          of the center point I,J,K.
c        - xi(3,im,jm,km), the interpolating coefficients for the vertex point
c          i,j,k of the mesh x.
c
c  ----------------------------------------------------------------------
c
c Goal : Give, in the mesh constituated by the face center (i), the cell I,J,K 
c       containing the vertex point i,j,k or if not available the nearest cell 
c       and the corresponding tri-linear interpolating coefficients.
c
c ----------------------------------------------------------------------
c
        dimension      xi(3,im ,jm ,km )
        dimension     ixi(3,im ,jm ,km )
c
        do 100 k=1,km
        do 100 j=1,jm
        do 100 i=1,im
         ixi(1,i,j,k) = i
         ixi(2,i,j,k) = j-1
         ixi(3,i,j,k) = k-1
         xi (1,i,j,k) =-1.0
         xi (2,i,j,k) = 0.0
         xi (3,i,j,k) = 0.0
 100    continue
c
        im1           = max(im-1,1)
        jm1           = max(jm-1,1)
        km1           = max(km-1,1)
        im2           = max(im-2,1)
        jm2           = max(jm-2,1)
        km2           = max(km-2,1)
c
        i2            = min(2,im1)
        j2            = min(2,jm1)
        k2            = min(2,km1)
c
        do 200 j=1,jm
        do 200 i=1,im
         ixi(3,i,j,1 ) = 1
         xi (3,i,j,1 ) =-2.0
         ixi(3,i,j,km) = km2 
         xi (3,i,j,km) = 2.0
 200    continue
c
        do 300 i=1,im
        do 300 k=1,km
         ixi(2,i,1 ,k) = 1
         xi (2,i,1 ,k) =-2.0
         ixi(2,i,jm,k) = jm2 
         xi (2,i,jm,k) = 2.0
 300    continue
c
        do 400 k=1,km
        do 400 j=1,jm
c        ixi(1,1 ,j,k) = 1
c        xi (1,1 ,j,k) =-1.0
c        ixi(1,im,j,k) = im1
c        xi (1,im,j,k) = 1.0
         ixi(1,1 ,j,k) = i2
         xi (1,1 ,j,k) =-1.0 - 2.0*(i2-1)
         ixi(1,im,j,k) = im2
         xi (1,im,j,k) = 1.0 + 2.0*(im1-im2)
 400    continue
c
        return
        end

c
c ======================================================================
        subroutine k6cj_v3d(im,jm,km,ixi,xi)
c ======================================================================
c
c       SUBROUTINE FOR FACE CENTER ( ALONG J PLANE) 
c       TO CELL VERTEX INTERPOLATION USED
c       IN A TRI-LINEAR INTERPOLATION PROCESS
c
c ----------------------------------------------------------------------
c
c AUTOR       : GREGORY COUSSEMENT    (ONERA:OAt1)
c VERSION     : 1.1
c DATE        : 01/02/92
c
c REFERENCES  : GREGORY COUSSEMENT AND OLIVER-PIERRE JACQUOTTE
c               "Structured Mesh Adaption: Space Accuracy and
c                Interpolation Methods.", Proceeding of the 2nd Workshop
c                on Reliability and Adaptive Methods in Computational Mechanics.
c                Krakow, POLAND, October 14-16, 1991.
c
c COPYRIGHT   : (c) 1991 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
c Entry: - im,jm,km, the dimension of the mesh x in each directions.
c Exit : - ixi(3,im,jm,km), the indices in the mesh x where the vertex point
c          i,j,k, of the mesh x, is contained in the mesh constituated
c          of the center point I,J,K.
c        - xi(3,im,jm,km), the interpolating coefficients for the vertex point
c          i,j,k of the mesh x.
c
c  ----------------------------------------------------------------------
c
c Goal : Give, in the mesh constituated by the face center (j), the cell I,J,K 
c       containing the vertex point i,j,k or if not available the nearest cell 
c       and the corresponding tri-linear interpolating coefficients.
c
c ----------------------------------------------------------------------
c
        dimension      xi(3,im ,jm ,km )
        dimension     ixi(3,im ,jm ,km )
c
        do 100 k=1,km
        do 100 j=1,jm
        do 100 i=1,im
         ixi(1,i,j,k) = i-1
         ixi(2,i,j,k) = j
         ixi(3,i,j,k) = k-1
         xi (1,i,j,k) = 0.0
         xi (2,i,j,k) =-1.0
         xi (3,i,j,k) = 0.0
 100    continue
c
        im1           = max(im-1,1)
        jm1           = max(jm-1,1)
        km1           = max(km-1,1)
        im2           = max(im-2,1)
        jm2           = max(jm-2,1)
        km2           = max(km-2,1)
c
        i2            = min(2,im1)
        j2            = min(2,jm1)
        k2            = min(2,km1)
c
        do 200 j=1,jm
        do 200 i=1,im
         ixi(3,i,j,1 ) = 1
         xi (3,i,j,1 ) =-2.0
         ixi(3,i,j,km) = km2 
         xi (3,i,j,km) = 2.0
 200    continue
c
        do 300 i=1,im
        do 300 k=1,km
c        ixi(2,i,1 ,k) = 1
c        xi (2,i,1 ,k) =-1.0
c        ixi(2,i,jm,k) = jm1 
c        xi (2,i,jm,k) = 1.0
         ixi(2,i,1 ,k) = j2 
         xi (2,i,1 ,k) =-1.0 - 2.0*(j2-1)
         ixi(2,i,jm,k) = jm2 
         xi (2,i,jm,k) = 1.0 + 2.0*(jm1-jm2)
 300    continue
c
        do 400 k=1,km
        do 400 j=1,jm
         ixi(1,1 ,j,k) = 1
         xi (1,1 ,j,k) =-2.0
         ixi(1,im,j,k) = im2 
         xi (1,im,j,k) = 2.0
 400    continue
c
        return
        end

c
c ======================================================================
        subroutine k6ck_v3d(im,jm,km,ixi,xi)
c ======================================================================
c
c       SUBROUTINE FOR FACE CENTER ( ALONG K PLANE) 
c       TO CELL VERTEX INTERPOLATION USED
c       IN A TRI-LINEAR INTERPOLATION PROCESS
c
c ----------------------------------------------------------------------
c
c AUTOR       : GREGORY COUSSEMENT    (ONERA:OAt1)
c VERSION     : 1.1
c DATE        : 01/02/92
c
c REFERENCES  : GREGORY COUSSEMENT AND OLIVER-PIERRE JACQUOTTE
c               "Structured Mesh Adaption: Space Accuracy and
c                Interpolation Methods.", Proceeding of the 2nd Workshop
c                on Reliability and Adaptive Methods in Computational Mechanics.
c                Krakow, POLAND, October 14-16, 1991.
c
c COPYRIGHT   : (c) 1991 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
c Entry: - im,jm,km, the dimension of the mesh x in each directions.
c Exit : - ixi(3,im,jm,km), the indices in the mesh x where the vertex point
c          i,j,k, of the mesh x, is contained in the mesh constituated
c          of the center point I,J,K.
c        - xi(3,im,jm,km), the interpolating coefficients for the vertex point
c          i,j,k of the mesh x.
c
c  ----------------------------------------------------------------------
c
c Goal : Give, in the mesh constituated by the face center (k), the cell I,J,K 
c       containing the vertex point i,j,k or if not available the nearest cell 
c       and the corresponding tri-linear interpolating coefficients.
c
c ----------------------------------------------------------------------
c
        dimension      xi(3,im ,jm ,km )
        dimension     ixi(3,im ,jm ,km )
c
        do 100 k=1,km
        do 100 j=1,jm
        do 100 i=1,im
         ixi(1,i,j,k) = i-1
         ixi(2,i,j,k) = j-1
         ixi(3,i,j,k) = k
         xi (1,i,j,k) = 0.0
         xi (2,i,j,k) = 0.0
         xi (3,i,j,k) =-1.0
 100    continue
c
        im1           = max(im-1,1)
        jm1           = max(jm-1,1)
        km1           = max(km-1,1)
        im2           = max(im-2,1)
        jm2           = max(jm-2,1)
        km2           = max(km-2,1)
c
        i2            = min(2,im1)
        j2            = min(2,jm1)
        k2            = min(2,km1)
c
        do 200 j=1,jm
        do 200 i=1,im
c        ixi(3,i,j,1 ) = 1
c        xi (3,i,j,1 ) =-1.0
c        ixi(3,i,j,km) = km1
c        xi (3,i,j,km) = 1.0
         ixi(3,i,j,1 ) = k2
         xi (3,i,j,1 ) =-1.0 - 2.0*(k2-1)
         ixi(3,i,j,km) = km2
         xi (3,i,j,km) = 1.0 + 2.0*(km1-km2)
 200    continue
c
        do 300 i=1,im
        do 300 k=1,km
         ixi(2,i,1 ,k) = 1
         xi (2,i,1 ,k) =-2.0
         ixi(2,i,jm,k) = jm2 
         xi (2,i,jm,k) = 2.0
 300    continue
c
        do 400 k=1,km
        do 400 j=1,jm
         ixi(1,1 ,j,k) = 1
         xi (1,1 ,j,k) =-2.0
         ixi(1,im,j,k) = im2 
         xi (1,im,j,k) = 2.0
 400    continue
c
        return
        end

c
c ======================================================================
        subroutine k6c_v3d(im,jm,km,ixi,xi)
c ======================================================================
c
c       SUBROUTINE FOR CELL CENTER TO CELL VERTEX INTERPOLATION USED
c       IN A TRI-LINEAR INTERPOLATION PROCESS
c
c ----------------------------------------------------------------------
c
c AUTOR       : GREGORY COUSSEMENT    (ONERA:OAt1)
c VERSION     : 1.1
c DATE        : 01/02/92
c
c REFERENCES  : GREGORY COUSSEMENT AND OLIVER-PIERRE JACQUOTTE
c               "Structured Mesh Adaption: Space Accuracy and
c                Interpolation Methods.", Proceeding of the 2nd Workshop
c                on Reliability and Adaptive Methods in Computational Mechanics.
c                Krakow, POLAND, October 14-16, 1991.
c
c COPYRIGHT   : (c) 1992 ONERA OAt1, All Rights Reserved
c
c ----------------------------------------------------------------------
c
c Entry: - im,jm,km, the dimension of the mesh x in each directions.
c Exit : - ixi(3,im,jm,km), the indices in the mesh x where the vertex point
c          i,j,k, of the mesh x, is contained in the mesh constituated
c          of the center point I,J,K.
c        - xi(3,im,jm,km), the interpolating coefficients for the vertex point
c          i,j,k of the mesh x.
c
c  ----------------------------------------------------------------------
c
c Goal : Give, in the mesh constituated by the cell center, the cell I,J,K 
c       containing the vertex point i,j,k or if not available the nearest cell 
c       and the corresponding tri-linear interpolating coefficients.
c
c ----------------------------------------------------------------------
c
        dimension      xi(3,im ,jm ,km )
        dimension     ixi(3,im ,jm ,km )
c
        do 100 k=1,km
        do 100 j=1,jm
        do 100 i=1,im
         ixi(1,i,j,k) = i-1
         ixi(2,i,j,k) = j-1
         ixi(3,i,j,k) = k-1
         xi (1,i,j,k) = 0.0
         xi (2,i,j,k) = 0.0
         xi (3,i,j,k) = 0.0
 100    continue
c
        im2           = max0(im-2,1)
        jm2           = max0(jm-2,1)
        km2           = max0(km-2,1)
c
        do 200 j=1,jm
        do 200 i=1,im
         ixi(3,i,j,1 ) = 1
         xi (3,i,j,1 ) = -2.0
         ixi(3,i,j,km) = km2 
         xi (3,i,j,km) = 2.0
 200    continue
c
        do 300 i=1,im
        do 300 k=1,km
         ixi(2,i,1 ,k) = 1
         xi (2,i,1 ,k) = -2.0
         ixi(2,i,jm,k) = jm2 
         xi (2,i,jm,k) = 2.0
 300    continue
c
        do 400 k=1,km
        do 400 j=1,jm
         ixi(1,1 ,j,k) = 1
         xi (1,1 ,j,k) = -2.0
         ixi(1,im,j,k) = im2 
         xi (1,im,j,k) = 2.0
 400    continue
c
        return
        end
