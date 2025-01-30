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
c 
c  Helicopter rotor motion
c  Use both the absolute frame of reference and the relative frame of 
c  reference of each blade
c
      SUBROUTINE evalrotfor(time, psi0, psi0_b,
     &                      blade_axis,
     &                      mux, muy, muz,
     &                      x0, y0, z0, axis0, alp0,
     &                      x1, y1, z1, axis1, omega,
     &                      x2, y2, z2, axis2, prelag,
     &                      x3, y3, z3, axis3, precon,
     &                      x4, y4, z4, axis4, del0,nhdel,delc,dels,
     &                      x5, y5, z5, axis5, bet0,nhbet,betc,bets,
     &                      x6, y6, z6, axis6, tet0,nhtet,tetc,tets,
     &                      psideg, deldeg, betdeg, tetdeg,
     &                      r0, rot, xr,
     &                      s0, omeg)
      IMPLICIT NONE
C
# include "Def/DefFortranConst.h"
c      
C_Intrinsic Functions
C
      INTRINSIC COS, SIN, ATAN
C
C =============================================================================
C_IN
      INTEGER_E nhbet           ! number of harmonics for flapping
      INTEGER_E nhtet           ! number of harmonics for pitching
      INTEGER_E nhdel           ! number of harmonics for lead-lag
C
      REAL_E x0,y0,z0,          ! origin of rotor shaft 
     &       x1,y1,z1,          ! origin of pre-lag 
     &       x2,y2,z2,          ! origin of pre-con
     &       x3,y3,z3,          ! origin of the helicopter rotor
     &       x4,y4,z4,          ! origin of lead-lag motion
     &       x5,y5,z5,          ! origin of flapping motion
     &       x6,y6,z6           ! origin of pitching motion
      REAL_E alp0               ! rotor shaft angle
      REAL_E prelag             ! pre-lag angle
      REAL_E precon             ! pre-con angle
c
      REAL_E bet0,betc(nhbet),bets(nhbet) !harmonics for flapping
      REAL_E del0,delc(nhdel),dels(nhdel) !harmonics for lead-lag
      REAL_E tet0,tetc(nhtet),tets(nhtet) !harmonics for pitching 
C
      REAL_E time               ! physical time
      REAL_E omega              ! rotor angular velocity
C
      REAL_E psi0               ! angle for initial rotor rotation (in degrees)
      REAL_E psi0_b             ! angle for blade position with respect to leading blade (in degrees)
      REAL_E mux,muy,muz        ! advance velocity components in the absolute frame
C
      INTEGER_E blade_axis      ! blade spanwise vector
C
      INTEGER_E axis0,          ! axis of rotor shaft
     &          axis1,          ! rotor axis
     &          axis2,          ! pre-lag axis
     &          axis3,          ! pre-con axis
     &          axis4,          ! axis of lead-lag motion
     &          axis5,          ! axis of flapping motion
     &          axis6           ! axis of pitching motion

C_OUT                                                                   
      REAL_E r0(3)              ! position of xr in AF 
      REAL_E rot(3,3)           ! rotation in time for instantaneous position vector in AF
      REAL_E xr(3)              ! center of rotation in RF
C
      REAL_E omeg(3)            ! angular velocity components in relative frame
      REAL_E s0(1:3)            ! translatory velocity components in the relative frame
      REAL_E psideg             ! azimuthal angle in degrees
      REAL_E betdeg,tetdeg,deldeg ! flapping, pitching abd lead-lag angles in degrees

C_LOCAL
      REAL_E pis2               ! =  pi / 2
      REAL_E rad2deg            !   
      REAL_E deg2rad            ! 
      REAL_E psi                ! rotor angle of rotation in the absolute frame
      REAL_E psim               ! blade azimuthal angle
      REAL_E omgt               ! omega * time
C
      REAL_E psi0r              !  psi0 in radians
      REAL_E psi0_br            !  psi0_b in radians
      REAL_E alp0r              !  rotor shaft angle in radians
      REAL_E prelagr            !  pre-lag angle in radians
      REAL_E preconr            !  pre-con angle in radians
c
      REAL_E bet0r,betcr(nhbet),betsr(nhbet) !harmonics for flapping in radians
      REAL_E del0r,delcr(nhdel),delsr(nhdel) !harmonics for lead-lag in radians
      REAL_E tet0r,tetcr(nhtet),tetsr(nhtet) !harmonics for pitching in radians
C
      REAL_E cpsi0br,spsi0br    ! cos(psi0_br), sin (psi0_br)
      REAL_E calp0,salp0        ! cos(alp0) and sin(alp0)
      REAL_E cprelag,sprelag    ! cos(pre-lag) and sin(pre-lag)
      REAL_E cprecon,sprecon    ! cos(pre-con) and sin(pre-con)
      REAL_E cpsim,spsim        ! cos(psim) and sin(psim)
      REAL_E psimp              ! time derivative of psim
      REAL_E bet                ! flapping angle
      REAL_E cbet,sbet          ! cos(bet) and sin (bet)
      REAL_E betp               ! time derivative of bet
      REAL_E del                ! lead-lag angle
      REAL_E cdel,sdel          ! cos(del) and sin (del)
      REAL_E delp               ! time derivative of del
      REAL_E tet                ! lead-lag angle
      REAL_E ctet,stet          ! cos(tet) and sin (tet)
      REAL_E tetp               ! time derivative of tet
C
      REAL_E c0(3,3),
     &       c1(3,3),
     &       c2(3,3),
     &       c3(3,3),
     &       c4(3,3),
     &       c5(3,3),
     &       c6(3,3),
     &       ca(3,3)  ! rotation matrices of blade motions
      REAL_E c0t(3,3),
     &       c1t(3,3),
     &       c2t(3,3),
     &       c3t(3,3),
     &       c4t(3,3),
     &       c5t(3,3),
     &       c6t(3,3),
     &       cat(3,3) ! transposed rotation matrices of blade motions
      REAL_E dc0(3,3),
     &       dc1(3,3),
     &       dc2(3,3),
     &       dc3(3,3),
     &       dc4(3,3),
     &       dc5(3,3),
     &       dc6(3,3) !time derivatives of rotation matrices of blade motions
      REAL_E rot0(3,3),
     &       rot1(3,3),
     &       rot2(3,3),
     &       rot3(3,3),
     &       rot4(3,3),
     &       rot5(3,3),
     &       rot6(3,3),
     &       rota(3,3),
     &       rotb(3,3) !local rotation matrices
C
      REAL_E tomeg(3,3)           !local angular velocity vector
      REAL_E omeg01,omeg02,omeg03 !local angular velocity vector components
      REAL_E omega1,omega2,omega3 !local angular velocity vector components
      REAL_E omeg11,omeg12,omeg13 !local angular velocity vector components
      REAL_E omeg21,omeg22,omeg23 !local angular velocity vector components
      REAL_E omeg31,omeg32,omeg33 !local angular velocity vector components
      REAL_E omeg41,omeg42,omeg43 !local angular velocity vector components
      REAL_E omeg51,omeg52,omeg53 !local angular velocity vector components
      REAL_E omeg61,omeg62,omeg63 !local angular velocity vector components
      REAL_E omegb1,omegb2,omegb3 !local angular velocity vector components
C
      REAL_E x01,y01,z01        !x1-x0,y1-y0,z1-z0 
      REAL_E x12,y12,z12        !x2-x1,y2-y1,z2-z1
      REAL_E x23,y23,z23        !x3-x2,y3-y2,z3-z2
      REAL_E x34,y34,z34        !x4-x3,y4-y3,z4-z3
      REAL_E x45,y45,z45        !x5-x4,y5-y4,z5-z4
      REAL_E x56,y56,z56        !x6-x5,y6-y5,z6-z5
      REAL_E x61,y61,z61        !x1-x6,y1-y6,z1-z6
C
      REAL_E v0x,v0y,v0z        !d/dt [x01,y01,z01] in absolute frame
      REAL_E v1x,v1y,v1z        !d/dt [x12,y12,z12] in absolute frame
      REAL_E v2x,v2y,v2z        !d/dt [x23,y23,z23] in absolute frame
      REAL_E v3x,v3y,v3z        !d/dt [x34,y34,z34] in absolute frame
      REAL_E v4x,v4y,v4z        !d/dt [x45,y45,z45] in absolute frame
      REAL_E v5x,v5y,v5z        !d/dt [x56,y56,z56] in absolute frame
      REAL_E v6x,v6y,v6z        !d/dt [x60,y60,z60] in absolute frame
C
      REAL_E v01,v02,v03
C
      INTEGER_E ih              ! loop index for harmonics
C=============================================================================
C
c$$$      write(*,*) '------------------------------------------------'
c$$$      write(*,*) '++++++ Argument list of TxbEvalRotForF.for +++++'
c$$$      write(*,*) '------------------------------------------------'
c$$$      write(*,*) 'time = ',time,' , psi0 =',psi0,' , psi0_b =',psi0
c$$$      write(*,*) ' '
c$$$      write(*,*) 'blade_axis =',blade_axis
c$$$      write(*,*) ' '
c$$$      write(*,*) 'mux = ',mux,' , my =',muy,' , muz =',muz
c$$$      write(*,*) ' '
c$$$      write(*,*) 'x0 = ',x0,' , y0 =',y0,' , z0 =',z0
c$$$      write(*,*) 'axis0 =',axis0,' , alp0 =',alp0
c$$$      write(*,*) ' '
c$$$      write(*,*) 'x1 = ',x1,' , y1 =',y1,' , z1 =',z1
c$$$      write(*,*) 'axis1 =',axis1,' , omega =',omega
c$$$      write(*,*) ' '
c$$$      write(*,*) 'x2 = ',x2,' , y2 =',y2,' , z2 =',z2
c$$$      write(*,*) 'axis2 =',axis2,' , pre_lag =',prelag
c$$$      write(*,*) ' '
c$$$      write(*,*) 'x3 = ',x3,' , y3 =',y3,' , z3 =',z3
c$$$      write(*,*) 'axis3 =',axis3,' , pre_con =',precon
c$$$      write(*,*) ' '
c$$$      write(*,*) 'x4 = ',x4,' , y4 =',y4,' , z4 =',z4
c$$$      write(*,*) 'axis4 =',axis4,' , del0 =',del0,' , nhdel =',nhdel
c$$$      write(*,*) ' '
c$$$      write(*,*) 'x5 = ',x5,' , y5 =',y5,' , z5 =',z5
c$$$      write(*,*) 'axis5 =',axis5,' , bet0 =',bet0,' , nhbet =',nhbet
c$$$      write(*,*) ' '
c$$$      write(*,*) 'x6 = ',x6,' , y6 =',y6,' , z6 =',z6
c$$$      write(*,*) 'axis6 =',axis6,' , tet0 =',tet0,' , nhtet =',nhtet
c$$$      write(*,*) '------------------------------------------------'
C
c Angle azimutal psim tenant compte de la postion relative des blocks
c
      pis2 = E_PI/2.
      deg2rad = E_PI/180.
      rad2deg = 180./E_PI
c
c Conversion des angles du mouvement en radians
c
      psi0r   = psi0*deg2rad
      psi0_br = psi0_b*deg2rad
c
      alp0r = alp0*deg2rad
      prelagr = prelag*deg2rad
      preconr = precon*deg2rad
c
      del0r = del0*deg2rad
c
      DO ih = 1, nhdel
         delcr(ih) = delc(ih)*deg2rad
         delsr(ih) = dels(ih)*deg2rad
      END DO
c
      bet0r = bet0*deg2rad
c
      DO ih = 1, nhbet
         betcr(ih) = betc(ih)*deg2rad
         betsr(ih) = bets(ih)*deg2rad
      END DO
c
      tet0r = tet0*deg2rad
c
      DO ih = 1, nhtet
         tetcr(ih) = tetc(ih)*deg2rad
         tetsr(ih) = tets(ih)*deg2rad
      END DO
c
      omgt = omega*time
      psim = omgt + psi0r

c Cas maillage HELICE avec x comme axe de rotation de rotor et y dans la 
c direction de la pale
      IF (axis1.EQ.1.AND.blade_axis.EQ.2) THEN 
         psi = psim

c Cas maillage HELICE avec x comme axe de rotation de rotor et z dans la 
c direction de la pale
      ELSEIF (axis1.EQ.1.AND.blade_axis.EQ.3) THEN 
         psi = psim + pis2

c Cas maillage CHIMERE avec z comme axe de rotation de rotor et x dans la 
c direction de la pale
      ELSEIF (axis1.EQ.3.AND.blade_axis.EQ.1) THEN
         psi = psim
      
c Cas maillage MULTIBLOCK avec z comme axe de rotation de rotor et y dans la 
c direction de la pale
      ELSEIF (axis1.EQ.3.AND.blade_axis.EQ.2) THEN
         psi = psim + pis2
      
      ELSE
          write(*,*) '------------------------------------------------'
          write(*,*) 'WARNING: A blade spanwise axis does not fit any '
          write(*,*) '          of the standard definitions admitted  '
          write(*,*) '          for helicopter blades or propellers : '
          write(*,*) ' 1) Rotation axis = 0z, Blade spanwise axis = Ox'
          write(*,*) ' 2) Rotation axis = 0z, Blade spanwise axis = Oy'
          write(*,*) ' 3) Rotation axis = 0x, Blade spanwise axis = Oy'
          write(*,*) ' 4) Rotation axis = 0x, Blade spanwise axis = Oz'
          write(*,*) ' See fortran file /Tbx/Function/TbxEvalRotF.for '
          write(*,*) '------------------------------------------------'
          psi = psim
      ENDIF

c
c  ---------------------------------------------------------------
c  ----- Cinematique particularisee a une application rotor ------
c  ---------------------------------------------------------------
c
c Tenseur 1er mvt --> inclinaison du rotor (ind. du temps)
c les angles de l'inclinaison rotor ne sont pas relatif a la pale
c
      calp0 = COS(alp0r)
      salp0 = SIN(alp0r)
c 
c tenseur rotation --> C0(i,j)
c
      IF (axis0.EQ.1) THEN
        c0(1,1) = 1.
        c0(1,2) = 0.
        c0(1,3) = 0.
        c0(2,1) = 0.
        c0(2,2) = calp0
        c0(2,3) = -salp0
        c0(3,1) = 0.
        c0(3,2) = salp0
        c0(3,3) = calp0
c
c tenseur vitesse d/dt C0(i,j)
c
        dc0(1,1) = 0.
        dc0(1,2) = 0.
        dc0(1,3) = 0.
        dc0(2,1) = 0.
        dc0(2,2) = 0.
        dc0(2,3) = 0.
        dc0(3,1) = 0.
        dc0(3,2) = 0.
        dc0(3,3) = 0.
c
      ELSEIF (axis0.EQ.2) THEN
        c0(1,1) = calp0
        c0(1,2) = 0.
        c0(1,3) = salp0
        c0(2,1) = 0.
        c0(2,2) = 1.
        c0(2,3) = 0.
        c0(3,1) = -salp0
        c0(3,2) = 0.
        c0(3,3) = calp0
c 
c tenseur vitesse d/dt C0(i,j)
c
        dc0(1,1) = 0.
        dc0(1,2) = 0.
        dc0(1,3) = 0.
        dc0(2,1) = 0.
        dc0(2,2) = 0.
        dc0(2,3) = 0.
        dc0(3,1) = 0.
        dc0(3,2) = 0.
        dc0(3,3) = 0.
cc
      ELSEIF(axis0.EQ.3) THEN
        c0(1,1) = calp0
        c0(1,2) = -salp0
        c0(1,3) = 0.
        c0(2,1) = salp0
        c0(2,2) = calp0
        c0(2,3) = 0.
        c0(3,1) = 0.
        c0(3,2) = 0.
        c0(3,3) = 1.
c 
c tenseur vitesse d/dt C0(i,j)
c
        dc0(1,1) = 0.
        dc0(1,2) = 0.
        dc0(1,3) = 0.
        dc0(2,1) = 0.
        dc0(2,2) = 0.
        dc0(2,3) = 0.
        dc0(3,1) = 0.
        dc0(3,2) = 0.
        dc0(3,3) = 0.

      ELSE

        WRITE(*,*) 'Warning: evalrotfor: incorrect axis1. Set to 1.'
        c0(1,1) = 1.
        c0(1,2) = 0.
        c0(1,3) = 0.
        c0(2,1) = 0.
        c0(2,2) = calp0
        c0(2,3) = -salp0
        c0(3,1) = 0.
        c0(3,2) = salp0
        c0(3,3) = calp0
        dc0(1,1) = 0.
        dc0(1,2) = 0.
        dc0(1,3) = 0.
        dc0(2,1) = 0.
        dc0(2,2) = 0.
        dc0(2,3) = 0.
        dc0(3,1) = 0.
        dc0(3,2) = 0.
        dc0(3,3) = 0.

      ENDIF
c
c Passage dans le repere de la premiere pale (dans laquelle est donnee la 
c cinematique) on suppose que le centre de rotation est X1
c
      cpsi0br = COS(psi0_br)
      spsi0br = -SIN(psi0_br)

      IF (axis1.EQ.1) THEN
        ca(1,1) = 1.
        ca(1,2) = 0.
        ca(1,3) = 0.
        ca(2,1) = 0.
        ca(2,2) = cpsi0br
        ca(2,3) = -spsi0br
        ca(3,1) = 0.
        ca(3,2) = spsi0br
        ca(3,3) = cpsi0br
c
      ELSEIF (axis1.EQ.2) THEN
        ca(1,1) = cpsi0br
        ca(1,2) = 0.
        ca(1,3) = spsi0br
        ca(2,1) = 0.
        ca(2,2) = 1.
        ca(2,3) = 0.
        ca(3,1) = -spsi0br
        ca(3,2) = 0.
        ca(3,3) = cpsi0br
c
      ELSEIF (axis1.EQ.3) THEN
        ca(1,1) = cpsi0br
        ca(1,2) = -spsi0br
        ca(1,3) = 0.
        ca(2,1) = spsi0br
        ca(2,2) = cpsi0br
        ca(2,3) = 0.
        ca(3,1) = 0.
        ca(3,2) = 0.
        ca(3,3) = 1.
      ELSE
        WRITE(*,*) 'Warning: evalrotfor: incorrect axis1. Set to 1.'
        ca(1,1) = 1.
        ca(1,2) = 0.
        ca(1,3) = 0.
        ca(2,1) = 0.
        ca(2,2) = cpsi0br
        ca(2,3) = -spsi0br
        ca(3,1) = 0.
        ca(3,2) = spsi0br
        ca(3,3) = cpsi0br
      ENDIF
c
c Tenseur 2eme mvt --> rotation uniforme
c
      cpsim  = COS(psim)
      spsim  = SIN(psim)
      psimp  = omega
c 
c tenseur rotation --> C1(i,j)
c
      IF (axis1.EQ.1) THEN
        c1(1,1) = 1.
        c1(1,2) = 0.
        c1(1,3) = 0.
        c1(2,1) = 0.
        c1(2,2) = cpsim
        c1(2,3) = -spsim
        c1(3,1) = 0.
        c1(3,2) = spsim
        c1(3,3) = cpsim
c
c tenseur vitesse d/dt C1(i,j)
c
        dc1(1,1) = 0.
        dc1(1,2) = 0.
        dc1(1,3) = 0.
        dc1(2,1) = 0.
        dc1(2,2) = -psimp*spsim
        dc1(2,3) = -psimp*cpsim
        dc1(3,1) = 0.
        dc1(3,2) = psimp*cpsim
        dc1(3,3) = -psimp*spsim
c
      ELSEIF (axis1.EQ.2) THEN
        c1(1,1) = cpsim
        c1(1,2) = 0.
        c1(1,3) = spsim
        c1(2,1) = 0.
        c1(2,2) = 1.
        c1(2,3) = 0.
        c1(3,1) = -spsim
        c1(3,2) = 0.
        c1(3,3) = cpsim
c 
c tenseur vitesse d/dt C1(i,j)
c
        dc1(1,1) = -psimp*spsim
        dc1(1,2) = 0.
        dc1(1,3) = psimp*cpsim
        dc1(2,1) = 0.
        dc1(2,2) = 0.
        dc1(2,3) = 0.
        dc1(3,1) = -psimp*cpsim
        dc1(3,2) = 0.
        dc1(3,3) = -psimp*spsim
c
      ELSEIF(axis1.EQ.3) THEN
        c1(1,1) = cpsim
        c1(1,2) = -spsim
        c1(1,3) = 0.
        c1(2,1) = spsim
        c1(2,2) = cpsim
        c1(2,3) = 0.
        c1(3,1) = 0.
        c1(3,2) = 0.
        c1(3,3) = 1.
c 
c tenseur vitesse d/dt C1(i,j)
c
        dc1(1,1) = -psimp*spsim
        dc1(1,2) = -psimp*cpsim
        dc1(1,3) = 0.
        dc1(2,1) = psimp*cpsim
        dc1(2,2) = -psimp*spsim
        dc1(2,3) = 0.
        dc1(3,1) = 0.
        dc1(3,2) = 0.
        dc1(3,3) = 0.

      ELSE
        WRITE(*,*) 'Warning: evalrotfor: incorrect axis1. Set to 1.'
        c1(1,1) = 1.
        c1(1,2) = 0.
        c1(1,3) = 0.
        c1(2,1) = 0.
        c1(2,2) = cpsim
        c1(2,3) = -spsim
        c1(3,1) = 0.
        c1(3,2) = spsim
        c1(3,3) = cpsim
        dc1(1,1) = 0.
        dc1(1,2) = 0.
        dc1(1,3) = 0.
        dc1(2,1) = 0.
        dc1(2,2) = -psimp*spsim
        dc1(2,3) = -psimp*cpsim
        dc1(3,1) = 0.
        dc1(3,2) = psimp*cpsim
        dc1(3,3) = -psimp*spsim

      ENDIF
c
c Tenseur 3eme mvt --> rotation de pre-traine (ind. du temps)
c
      cprelag = COS(prelagr)
      sprelag = SIN(prelagr)
c 
c tenseur rotation --> C2(i,j)
c
      IF (axis2.EQ.1) THEN
        c2(1,1) = 1.
        c2(1,2) = 0.
        c2(1,3) = 0.
        c2(2,1) = 0.
        c2(2,2) = cprelag
        c2(2,3) = -sprelag
        c2(3,1) = 0.
        c2(3,2) = sprelag
        c2(3,3) = cprelag
c 
c tenseur vitesse d/dt C2(i,j)
c
        dc2(1,1) = 0.
        dc2(1,2) = 0.
        dc2(1,3) = 0.
        dc2(2,1) = 0.
        dc2(2,2) = 0.
        dc2(2,3) = 0.
        dc2(3,1) = 0.
        dc2(3,2) = 0.
        dc2(3,3) = 0.
c
      ELSEIF (axis2.EQ.2) THEN
        c2(1,1) = cprelag
        c2(1,2) = 0.
        c2(1,3) = sprelag
        c2(2,1) = 0.
        c2(2,2) = 1.
        c2(2,3) = 0.
        c2(3,1) = -sprelag
        c2(3,2) = 0.
        c2(3,3) = cprelag
c 
c tenseur vitesse d/dt C2(i,j)
c
        dc2(1,1) = 0.
        dc2(1,2) = 0.
        dc2(1,3) = 0.
        dc2(2,1) = 0.
        dc2(2,2) = 0.
        dc2(2,3) = 0.
        dc2(3,1) = 0.
        dc2(3,2) = 0.
        dc2(3,3) = 0.
cc
      ELSEIF (axis2.EQ.3) THEN
         c2(1,1) = cprelag
         c2(1,2) = -sprelag
         c2(1,3) = 0.
         c2(2,1) = sprelag
         c2(2,2) = cprelag
         c2(2,3) = 0.
         c2(3,1) = 0.
         c2(3,2) = 0.
         c2(3,3) = 1.
c 
c tenseur vitesse d/dt C2(i,j)
c
         dc2(1,1) = 0.
         dc2(1,2) = 0.
         dc2(1,3) = 0.
         dc2(2,1) = 0.
         dc2(2,2) = 0.
         dc2(2,3) = 0.
         dc2(3,1) = 0.
         dc2(3,2) = 0.
         dc2(3,3) = 0.
      
      ELSE

        WRITE(*,*) 'Warning: evalrotfor: incorrect axis2. Set to 1.'
        c2(1,1) = 1.
        c2(1,2) = 0.
        c2(1,3) = 0.
        c2(2,1) = 0.
        c2(2,2) = cprelag
        c2(2,3) = -sprelag
        c2(3,1) = 0.
        c2(3,2) = sprelag
        c2(3,3) = cprelag
        dc2(1,1) = 0.
        dc2(1,2) = 0.
        dc2(1,3) = 0.
        dc2(2,1) = 0.
        dc2(2,2) = 0.
        dc2(2,3) = 0.
        dc2(3,1) = 0.
        dc2(3,2) = 0.
        dc2(3,3) = 0.
      ENDIF
c
c Tenseur 4eme mvt --> rotation de pre-conicite (ind. du temps)
c
      cprecon = COS(preconr)
      sprecon = SIN(preconr)
c 
c tenseur rotation --> C3(i,j)
c
      IF (axis3.EQ.1) THEN
         c3(1,1) = 1.
         c3(1,2) = 0.
         c3(1,3) = 0.
         c3(2,1) = 0.
         c3(2,2) = cprecon
         c3(2,3) = -sprecon
         c3(3,1) = 0.
         c3(3,2) = sprecon
         c3(3,3) = cprecon
c 
c tenseur vitesse d/dt C3(i,j)
c
         dc3(1,1) = 0.
         dc3(1,2) = 0.
         dc3(1,3) = 0.
         dc3(2,1) = 0.
         dc3(2,2) = 0.
         dc3(2,3) = 0.
         dc3(3,1) = 0.
         dc3(3,2) = 0.
         dc3(3,3) = 0.
c
      ELSEIF(axis3.EQ.2) THEN
         c3(1,1) = cprecon
         c3(1,2) = 0.
         c3(1,3) = sprecon
         c3(2,1) = 0.
         c3(2,2) = 1.
         c3(2,3) = 0.
         c3(3,1) = -sprecon
         c3(3,2) = 0.
         c3(3,3) = cprecon
c 
c tenseur vitesse d/dt C3(i,j)
c
         dc3(1,1) = 0.
         dc3(1,2) = 0.
         dc3(1,3) = 0.
         dc3(2,1) = 0.
         dc3(2,2) = 0.
         dc3(2,3) = 0.
         dc3(3,1) = 0.
         dc3(3,2) = 0.
         dc3(3,3) = 0.
cc
      ELSEIF (axis3.EQ.3) THEN
         c3(1,1) = cprecon
         c3(1,2) = -sprecon
         c3(1,3) = 0.
         c3(2,1) = sprecon
         c3(2,2) = cprecon
         c3(2,3) = 0.
         c3(3,1) = 0.
         c3(3,2) = 0.
         c3(3,3) = 1.
c 
c tenseur vitesse d/dt C3(i,j)
c
         dc3(1,1) = 0.
         dc3(1,2) = 0.
         dc3(1,3) = 0.
         dc3(2,1) = 0.
         dc3(2,2) = 0.
         dc3(2,3) = 0.
         dc3(3,1) = 0.
         dc3(3,2) = 0.
         dc3(3,3) = 0.
      ELSE
         WRITE(*,*) 'Warning: evalrotfor: incorrect axis3. Set to 1.'
         c3(1,1) = 1.
         c3(1,2) = 0.
         c3(1,3) = 0.
         c3(2,1) = 0.
         c3(2,2) = cprecon
         c3(2,3) = -sprecon
         c3(3,1) = 0.
         c3(3,2) = sprecon
         c3(3,3) = cprecon
c 
c tenseur vitesse d/dt C3(i,j)
c
         dc3(1,1) = 0.
         dc3(1,2) = 0.
         dc3(1,3) = 0.
         dc3(2,1) = 0.
         dc3(2,2) = 0.
         dc3(2,3) = 0.
         dc3(3,1) = 0.
         dc3(3,2) = 0.
         dc3(3,3) = 0.

      ENDIF
c
c Les autres mouvements sont fonction de l'azimut psi+psi0_b
c
      psi = psi+psi0_br
      psideg = psi*rad2deg
c
c Tenseur 5eme mvt --> angle de trainee
c
      del = del0r
      DO ih = 1,nhdel
         del = del + delcr(ih)*COS(float(ih)*psi)
     &             + delsr(ih)*SIN(float(ih)*psi)
      END DO
c
      deldeg = del*rad2deg
      cdel = COS(del)
      sdel = SIN(del)
c
      delp = 0.
      DO ih = 1,nhdel
         delp = delp + omega * float(ih) *
     &        (-delcr(ih)*SIN(float(ih)*psi)
     &         +delsr(ih)*COS(float(ih)*psi))
      END DO
c
c tenseur --> C4(i,j)
c
      IF (axis4.EQ.1) THEN
        c4(1,1) = 1.
        c4(1,2) = 0.
        c4(1,3) = 0.
        c4(2,1) = 0.
        c4(2,2) = cdel
        c4(2,3) = -sdel
        c4(3,1) = 0.
        c4(3,2) = sdel
        c4(3,3) = cdel
c 
c tenseur vitesse d/dt C4(i,j)
c
        dc4(1,1) = 0.
        dc4(1,2) = 0.
        dc4(1,3) = 0.
        dc4(2,1) = 0.
        dc4(2,2) = -delp*sdel
        dc4(2,3) = -delp*cdel
        dc4(3,1) = 0.
        dc4(3,2) = delp*cdel
        dc4(3,3) = -delp*sdel
c
      ELSEIF (axis4.EQ.2) THEN
        c4(1,1) = cdel
        c4(1,2) = 0.
        c4(1,3) = sdel
        c4(2,1) = 0.
        c4(2,2) = 1.
        c4(2,3) = 0.
        c4(3,1) = -sdel
        c4(3,2) = 0.
        c4(3,3) = cdel
c 
c tenseur vitesse d/dt C4(i,j)
c
        dc4(1,1) = -delp*sdel
        dc4(1,2) = 0.
        dc4(1,3) = delp*cdel
        dc4(2,1) = 0.
        dc4(2,2) = 0.
        dc4(2,3) = 0.
        dc4(3,1) = -delp*cdel
        dc4(3,2) = 0.
        dc4(3,3) = -delp*sdel
c
      ELSEIF (axis4.EQ.3) THEN
        c4(1,1) = cdel
        c4(1,2) = -sdel
        c4(1,3) = 0.
        c4(2,1) = sdel
        c4(2,2) = cdel
        c4(2,3) = 0.
        c4(3,1) = 0.
        c4(3,2) = 0.
        c4(3,3) = 1.
c
c tenseur vitesse d/dt C4(i,j)
c
        dc4(1,1) = -delp*sdel
        dc4(1,2) = -delp*cdel
        dc4(1,3) = 0.
        dc4(2,1) = delp*cdel
        dc4(2,2) = -delp*sdel
        dc4(2,3) = 0.
        dc4(3,1) = 0.
        dc4(3,2) = 0.
        dc4(3,3) = 0.

      ELSE

        WRITE(*,*) 'Warning: evalrotfor: incorrect axis4. Set to 1.' 
        c4(1,1) = 1.
        c4(1,2) = 0.
        c4(1,3) = 0.
        c4(2,1) = 0.
        c4(2,2) = cdel
        c4(2,3) = -sdel
        c4(3,1) = 0.
        c4(3,2) = sdel
        c4(3,3) = cdel
        dc4(1,1) = 0.
        dc4(1,2) = 0.
        dc4(1,3) = 0.
        dc4(2,1) = 0.
        dc4(2,2) = -delp*sdel
        dc4(2,3) = -delp*cdel
        dc4(3,1) = 0.
        dc4(3,2) = delp*cdel
        dc4(3,3) = -delp*sdel

      ENDIF

c
c Tenseur 6eme mvt --> angle de battement
c
      bet = bet0r
      DO ih = 1,nhbet
         bet = bet + betcr(ih)*COS(float(ih)*psi)
     &             + betsr(ih)*SIN(float(ih)*psi)
      END DO
C
      betdeg = bet*rad2deg
      cbet = COS(bet)
      sbet = SIN(bet)
c
      betp = 0.
      DO ih = 1,nhbet
         betp = betp + omega * float(ih) *
     &        (-betcr(ih)*SIN(float(ih)*psi)
     &         +betsr(ih)*COS(float(ih)*psi))
      END DO
c
c tenseur --> C5(i,j)
c
      IF (axis5.EQ.1) THEN
        c5(1,1) = 1.
        c5(1,2) = 0.
        c5(1,3) = 0.
        c5(2,1) = 0.
        c5(2,2) = cbet
        c5(2,3) = -sbet
        c5(3,1) = 0.
        c5(3,2) = sbet
        c5(3,3) = cbet
c 
c tenseur vitesse d/dt C5(i,j)
c
        dc5(1,1) = 0.
        dc5(1,2) = 0.
        dc5(1,3) = 0.
        dc5(2,1) = 0.
        dc5(2,2) = -betp*sbet
        dc5(2,3) = -betp*cbet
        dc5(3,1) = 0.
        dc5(3,2) = betp*cbet
        dc5(3,3) = -betp*sbet
c
      ELSEIF (axis5.EQ.2) THEN
        c5(1,1) = cbet
        c5(1,2) = 0.
        c5(1,3) = sbet
        c5(2,1) = 0.
        c5(2,2) = 1.
        c5(2,3) = 0.
        c5(3,1) = -sbet
        c5(3,2) = 0.
        c5(3,3) = cbet
c 
c tenseur vitesse d/dt C5(i,j)
c
        dc5(1,1) = -betp*sbet
        dc5(1,2) = 0.
        dc5(1,3) = betp*cbet
        dc5(2,1) = 0.
        dc5(2,2) = 0.
        dc5(2,3) = 0.
        dc5(3,1) = -betp*cbet
        dc5(3,2) = 0.
        dc5(3,3) = -betp*sbet
c
      ELSEIF (axis5.EQ.3) THEN
        c5(1,1) = cbet
        c5(1,2) = -sbet
        c5(1,3) = 0.
        c5(2,1) = sbet
        c5(2,2) = cbet
        c5(2,3) = 0.
        c5(3,1) = 0.
        c5(3,2) = 0.
        c5(3,3) = 1.
c 
c tenseur vitesse d/dt C5(i,j)
c
        dc5(1,1) = -betp*sbet
        dc5(1,2) = -betp*cbet
        dc5(1,3) = 0.
        dc5(2,1) = betp*cbet
        dc5(2,2) = -betp*sbet
        dc5(2,3) = 0.
        dc5(3,1) = 0.
        dc5(3,2) = 0.
        dc5(3,3) = 0.

      ELSE
        WRITE(*,*) 'Warning: evalrotfor: incorrect axis5. Set to 1.'
        c5(1,1) = 1.
        c5(1,2) = 0.
        c5(1,3) = 0.
        c5(2,1) = 0.
        c5(2,2) = cbet
        c5(2,3) = -sbet
        c5(3,1) = 0.
        c5(3,2) = sbet
        c5(3,3) = cbet
        dc5(1,1) = 0.
        dc5(1,2) = 0.
        dc5(1,3) = 0.
        dc5(2,1) = 0.
        dc5(2,2) = -betp*sbet
        dc5(2,3) = -betp*cbet
        dc5(3,1) = 0.
        dc5(3,2) = betp*cbet
        dc5(3,3) = -betp*sbet

      ENDIF
c
c Tenseur 7eme mvt --> pas cyclique 
c
      tet = tet0r
      DO ih = 1,nhtet
         tet = tet + tetcr(ih)*COS(float(ih)*psi)
     &             + tetsr(ih)*SIN(float(ih)*psi)
      END DO
c
      tetdeg = tet*rad2deg
      ctet = COS(tet)
      stet = SIN(tet)
c
      tetp = 0.
      DO ih = 1,nhtet
         tetp = tetp + omega * float(ih) *
     &        (-tetcr(ih)*SIN(float(ih)*psi)
     &         +tetsr(ih)*COS(float(ih)*psi))
      END DO
c
c tenseur --> C6(i,j)
c
      IF (axis6.EQ.1) THEN
        c6(1,1) = 1.
        c6(1,2) = 0.
        c6(1,3) = 0.
        c6(2,1) = 0.
        c6(2,2) = ctet
        c6(2,3) = -stet
        c6(3,1) = 0.
        c6(3,2) = stet
        c6(3,3) = ctet
c
c tenseur vitesse d/dt C6(i,j)
c
        dc6(1,1) = 0.
        dc6(1,2) = 0.
        dc6(1,3) = 0.
        dc6(2,1) = 0.
        dc6(2,2) = -tetp*stet
        dc6(2,3) = -tetp*ctet
        dc6(3,1) = 0.
        dc6(3,2) = tetp*ctet
        dc6(3,3) = -tetp*stet
c
      ELSEIF (axis6.EQ.2) THEN
        c6(1,1) = ctet
        c6(1,2) = 0.
        c6(1,3) = stet
        c6(2,1) = 0.
        c6(2,2) = 1.
        c6(2,3) = 0.
        c6(3,1) = -stet
        c6(3,2) = 0.
        c6(3,3) = ctet
c
c tenseur vitesse d/dt C6(i,j)
c
        dc6(1,1) = -tetp*stet
        dc6(1,2) = 0.
        dc6(1,3) = tetp*ctet
        dc6(2,1) = 0.
        dc6(2,2) = 0.
        dc6(2,3) = 0.
        dc6(3,1) = -tetp*ctet
        dc6(3,2) = 0.
        dc6(3,3) = -tetp*stet
c
      ELSEIF (axis6.EQ.3) THEN
        c6(1,1) = ctet
        c6(1,2) = -stet
        c6(1,3) = 0.
        c6(2,1) = stet
        c6(2,2) = ctet
        c6(2,3) = 0.
        c6(3,1) = 0.
        c6(3,2) = 0.
        c6(3,3) = 1.
c
c tenseur vitesse d/dt C6(i,j)
c
        dc6(1,1) = -tetp*stet
        dc6(1,2) = -tetp*ctet
        dc6(1,3) = 0.
        dc6(2,1) = tetp*ctet
        dc6(2,2) = -tetp*stet
        dc6(2,3) = 0.
        dc6(3,1) = 0.
        dc6(3,2) = 0.
        dc6(3,3) = 0.

      ELSE

        WRITE(*,*) 'Warning: evalrotfor: incorrect axis6. Set to 1.'
        c6(1,1) = 1.
        c6(1,2) = 0.
        c6(1,3) = 0.
        c6(2,1) = 0.
        c6(2,2) = ctet
        c6(2,3) = -stet
        c6(3,1) = 0.
        c6(3,2) = stet
        c6(3,3) = ctet
        dc6(1,1) = 0.
        dc6(1,2) = 0.
        dc6(1,3) = 0.
        dc6(2,1) = 0.
        dc6(2,2) = -tetp*stet
        dc6(2,3) = -tetp*ctet
        dc6(3,1) = 0.
        dc6(3,2) = tetp*ctet
        dc6(3,3) = -tetp*stet

      ENDIF
c
c  --------------------------------------------------------
c  -------------- Fin donnees rotor -----------------------
c  --------------------------------------------------------
c
c Calculs preliminaires: transposees 
c
c C0 transposee --> C0T(i,j)
c
      c0t(1,1) = c0(1,1)
      c0t(1,2) = c0(2,1)
      c0t(1,3) = c0(3,1)
      c0t(2,1) = c0(1,2)
      c0t(2,2) = c0(2,2)
      c0t(2,3) = c0(3,2)
      c0t(3,1) = c0(1,3)
      c0t(3,2) = c0(2,3)
      c0t(3,3) = c0(3,3)
c
c CA transposee --> CAT(i,j)
c
      cat(1,1) = ca(1,1)
      cat(1,2) = ca(2,1)
      cat(1,3) = ca(3,1)
      cat(2,1) = ca(1,2)
      cat(2,2) = ca(2,2)
      cat(2,3) = ca(3,2)
      cat(3,1) = ca(1,3)
      cat(3,2) = ca(2,3)
      cat(3,3) = ca(3,3)
c
c C1 transposee --> C1T(i,j)
c
      c1t(1,1) = c1(1,1)
      c1t(1,2) = c1(2,1)
      c1t(1,3) = c1(3,1)
      c1t(2,1) = c1(1,2)
      c1t(2,2) = c1(2,2)
      c1t(2,3) = c1(3,2)
      c1t(3,1) = c1(1,3)
      c1t(3,2) = c1(2,3)
      c1t(3,3) = c1(3,3)
c
c C2 transposee --> C2T(i,j)
c
      c2t(1,1) = c2(1,1)
      c2t(1,2) = c2(2,1)
      c2t(1,3) = c2(3,1)
      c2t(2,1) = c2(1,2)
      c2t(2,2) = c2(2,2)
      c2t(2,3) = c2(3,2)
      c2t(3,1) = c2(1,3)
      c2t(3,2) = c2(2,3)
      c2t(3,3) = c2(3,3)
c
c C3 transposee --> C3T(i,j)
c
      c3t(1,1) = c3(1,1)
      c3t(1,2) = c3(2,1)
      c3t(1,3) = c3(3,1)
      c3t(2,1) = c3(1,2)
      c3t(2,2) = c3(2,2)
      c3t(2,3) = c3(3,2)
      c3t(3,1) = c3(1,3)
      c3t(3,2) = c3(2,3)
      c3t(3,3) = c3(3,3)
c
c C4 transposee --> C4T(i,j)
c
      c4t(1,1) = c4(1,1)
      c4t(1,2) = c4(2,1)
      c4t(1,3) = c4(3,1)
      c4t(2,1) = c4(1,2)
      c4t(2,2) = c4(2,2)
      c4t(2,3) = c4(3,2)
      c4t(3,1) = c4(1,3)
      c4t(3,2) = c4(2,3)
      c4t(3,3) = c4(3,3)
c
c C5 transposee --> C5T(i,j)
c
      c5t(1,1) = c5(1,1)
      c5t(1,2) = c5(2,1)
      c5t(1,3) = c5(3,1)
      c5t(2,1) = c5(1,2)
      c5t(2,2) = c5(2,2)
      c5t(2,3) = c5(3,2)
      c5t(3,1) = c5(1,3)
      c5t(3,2) = c5(2,3)
      c5t(3,3) = c5(3,3)
c
c C6 transposee --> C6T(i,j)
c
      c6t(1,1) = c6(1,1)
      c6t(1,2) = c6(2,1)
      c6t(1,3) = c6(3,1)
      c6t(2,1) = c6(1,2)
      c6t(2,2) = c6(2,2)
      c6t(2,3) = c6(3,2)
      c6t(3,1) = c6(1,3)
      c6t(3,2) = c6(2,3)
      c6t(3,3) = c6(3,3)      
c
c Calculs preliminaires: matrices de passages successives
c
      rot0(1,1) = c0(1,1)
      rot0(1,2) = c0(1,2)
      rot0(1,3) = c0(1,3)
      rot0(2,1) = c0(2,1)
      rot0(2,2) = c0(2,2)
      rot0(2,3) = c0(2,3)
      rot0(3,1) = c0(3,1)
      rot0(3,2) = c0(3,2)
      rot0(3,3) = c0(3,3)
c
c rota = rot0.cat
c
      rota(1,1)=rot0(1,1)*cat(1,1)+rot0(1,2)*cat(2,1)+rot0(1,3)*cat(3,1)
      rota(1,2)=rot0(1,1)*cat(1,2)+rot0(1,2)*cat(2,2)+rot0(1,3)*cat(3,2)
      rota(1,3)=rot0(1,1)*cat(1,3)+rot0(1,2)*cat(2,3)+rot0(1,3)*cat(3,3)
c
      rota(2,1)=rot0(2,1)*cat(1,1)+rot0(2,2)*cat(2,1)+rot0(2,3)*cat(3,1)
      rota(2,2)=rot0(2,1)*cat(1,2)+rot0(2,2)*cat(2,2)+rot0(2,3)*cat(3,2)
      rota(2,3)=rot0(2,1)*cat(1,3)+rot0(2,2)*cat(2,3)+rot0(2,3)*cat(3,3)
c
      rota(3,1)=rot0(3,1)*cat(1,1)+rot0(3,2)*cat(2,1)+rot0(3,3)*cat(3,1)
      rota(3,2)=rot0(3,1)*cat(1,2)+rot0(3,2)*cat(2,2)+rot0(3,3)*cat(3,2)
      rota(3,3)=rot0(3,1)*cat(1,3)+rot0(3,2)*cat(2,3)+rot0(3,3)*cat(3,3)
c
c rot1 = rota.c1 = c0.c1
c
      rot1(1,1) = rota(1,1)*c1(1,1)+rota(1,2)*c1(2,1)+rota(1,3)*c1(3,1)
      rot1(1,2) = rota(1,1)*c1(1,2)+rota(1,2)*c1(2,2)+rota(1,3)*c1(3,2)
      rot1(1,3) = rota(1,1)*c1(1,3)+rota(1,2)*c1(2,3)+rota(1,3)*c1(3,3)
c
      rot1(2,1) = rota(2,1)*c1(1,1)+rota(2,2)*c1(2,1)+rota(2,3)*c1(3,1)
      rot1(2,2) = rota(2,1)*c1(1,2)+rota(2,2)*c1(2,2)+rota(2,3)*c1(3,2)
      rot1(2,3) = rota(2,1)*c1(1,3)+rota(2,2)*c1(2,3)+rota(2,3)*c1(3,3)
c
      rot1(3,1) = rota(3,1)*c1(1,1)+rota(3,2)*c1(2,1)+rota(3,3)*c1(3,1)
      rot1(3,2) = rota(3,1)*c1(1,2)+rota(3,2)*c1(2,2)+rota(3,3)*c1(3,2)
      rot1(3,3) = rota(3,1)*c1(1,3)+rota(3,2)*c1(2,3)+rota(3,3)*c1(3,3)
c
c rot2 = rot1.c2 = c0.c1.c2
c
      rot2(1,1) = rot1(1,1)*c2(1,1)+rot1(1,2)*c2(2,1)+rot1(1,3)*c2(3,1)
      rot2(1,2) = rot1(1,1)*c2(1,2)+rot1(1,2)*c2(2,2)+rot1(1,3)*c2(3,2)
      rot2(1,3) = rot1(1,1)*c2(1,3)+rot1(1,2)*c2(2,3)+rot1(1,3)*c2(3,3)
c
      rot2(2,1) = rot1(2,1)*c2(1,1)+rot1(2,2)*c2(2,1)+rot1(2,3)*c2(3,1)
      rot2(2,2) = rot1(2,1)*c2(1,2)+rot1(2,2)*c2(2,2)+rot1(2,3)*c2(3,2)
      rot2(2,3) = rot1(2,1)*c2(1,3)+rot1(2,2)*c2(2,3)+rot1(2,3)*c2(3,3)
c
      rot2(3,1) = rot1(3,1)*c2(1,1)+rot1(3,2)*c2(2,1)+rot1(3,3)*c2(3,1)
      rot2(3,2) = rot1(3,1)*c2(1,2)+rot1(3,2)*c2(2,2)+rot1(3,3)*c2(3,2)
      rot2(3,3) = rot1(3,1)*c2(1,3)+rot1(3,2)*c2(2,3)+rot1(3,3)*c2(3,3)
c
c rot3 = rot2.c3 = c0.c1.c2.c3
c
      rot3(1,1) = rot2(1,1)*c3(1,1)+rot2(1,2)*c3(2,1)+rot2(1,3)*c3(3,1)
      rot3(1,2) = rot2(1,1)*c3(1,2)+rot2(1,2)*c3(2,2)+rot2(1,3)*c3(3,2)
      rot3(1,3) = rot2(1,1)*c3(1,3)+rot2(1,2)*c3(2,3)+rot2(1,3)*c3(3,3)
c
      rot3(2,1) = rot2(2,1)*c3(1,1)+rot2(2,2)*c3(2,1)+rot2(2,3)*c3(3,1)
      rot3(2,2) = rot2(2,1)*c3(1,2)+rot2(2,2)*c3(2,2)+rot2(2,3)*c3(3,2)
      rot3(2,3) = rot2(2,1)*c3(1,3)+rot2(2,2)*c3(2,3)+rot2(2,3)*c3(3,3)
c
      rot3(3,1) = rot2(3,1)*c3(1,1)+rot2(3,2)*c3(2,1)+rot2(3,3)*c3(3,1)
      rot3(3,2) = rot2(3,1)*c3(1,2)+rot2(3,2)*c3(2,2)+rot2(3,3)*c3(3,2)
      rot3(3,3) = rot2(3,1)*c3(1,3)+rot2(3,2)*c3(2,3)+rot2(3,3)*c3(3,3)
c
c rot4 = rot3.c4 = c0.c1.c2.c3.c4
c
      rot4(1,1) = rot3(1,1)*c4(1,1)+rot3(1,2)*c4(2,1)+rot3(1,3)*c4(3,1)
      rot4(1,2) = rot3(1,1)*c4(1,2)+rot3(1,2)*c4(2,2)+rot3(1,3)*c4(3,2)
      rot4(1,3) = rot3(1,1)*c4(1,3)+rot3(1,2)*c4(2,3)+rot3(1,3)*c4(3,3)
c
      rot4(2,1) = rot3(2,1)*c4(1,1)+rot3(2,2)*c4(2,1)+rot3(2,3)*c4(3,1)
      rot4(2,2) = rot3(2,1)*c4(1,2)+rot3(2,2)*c4(2,2)+rot3(2,3)*c4(3,2)
      rot4(2,3) = rot3(2,1)*c4(1,3)+rot3(2,2)*c4(2,3)+rot3(2,3)*c4(3,3)
c
      rot4(3,1) = rot3(3,1)*c4(1,1)+rot3(3,2)*c4(2,1)+rot3(3,3)*c4(3,1)
      rot4(3,2) = rot3(3,1)*c4(1,2)+rot3(3,2)*c4(2,2)+rot3(3,3)*c4(3,2)
      rot4(3,3) = rot3(3,1)*c4(1,3)+rot3(3,2)*c4(2,3)+rot3(3,3)*c4(3,3)
c
c rot5 = rot4.c5 = c0.c1.c2.c3.c4.c5
c
      rot5(1,1) = rot4(1,1)*c5(1,1)+rot4(1,2)*c5(2,1)+rot4(1,3)*c5(3,1)
      rot5(1,2) = rot4(1,1)*c5(1,2)+rot4(1,2)*c5(2,2)+rot4(1,3)*c5(3,2)
      rot5(1,3) = rot4(1,1)*c5(1,3)+rot4(1,2)*c5(2,3)+rot4(1,3)*c5(3,3)
c
      rot5(2,1) = rot4(2,1)*c5(1,1)+rot4(2,2)*c5(2,1)+rot4(2,3)*c5(3,1)
      rot5(2,2) = rot4(2,1)*c5(1,2)+rot4(2,2)*c5(2,2)+rot4(2,3)*c5(3,2)
      rot5(2,3) = rot4(2,1)*c5(1,3)+rot4(2,2)*c5(2,3)+rot4(2,3)*c5(3,3)
c
      rot5(3,1) = rot4(3,1)*c5(1,1)+rot4(3,2)*c5(2,1)+rot4(3,3)*c5(3,1)
      rot5(3,2) = rot4(3,1)*c5(1,2)+rot4(3,2)*c5(2,2)+rot4(3,3)*c5(3,2)
      rot5(3,3) = rot4(3,1)*c5(1,3)+rot4(3,2)*c5(2,3)+rot4(3,3)*c5(3,3)
c
c rot6 = rot5.c6 = c0.c1.c2.c3.c4.c5.c6
c
      rot6(1,1) = rot5(1,1)*c6(1,1)+rot5(1,2)*c6(2,1)+rot5(1,3)*c6(3,1)
      rot6(1,2) = rot5(1,1)*c6(1,2)+rot5(1,2)*c6(2,2)+rot5(1,3)*c6(3,2)
      rot6(1,3) = rot5(1,1)*c6(1,3)+rot5(1,2)*c6(2,3)+rot5(1,3)*c6(3,3)
c
      rot6(2,1) = rot5(2,1)*c6(1,1)+rot5(2,2)*c6(2,1)+rot5(2,3)*c6(3,1)
      rot6(2,2) = rot5(2,1)*c6(1,2)+rot5(2,2)*c6(2,2)+rot5(2,3)*c6(3,2)
      rot6(2,3) = rot5(2,1)*c6(1,3)+rot5(2,2)*c6(2,3)+rot5(2,3)*c6(3,3)
c
      rot6(3,1) = rot5(3,1)*c6(1,1)+rot5(3,2)*c6(2,1)+rot5(3,3)*c6(3,1)
      rot6(3,2) = rot5(3,1)*c6(1,2)+rot5(3,2)*c6(2,2)+rot5(3,3)*c6(3,2)
      rot6(3,3) = rot5(3,1)*c6(1,3)+rot5(3,2)*c6(2,3)+rot5(3,3)*c6(3,3)
c
c rotb final (retour dans le repere de depart)
c
      rotb(1,1)=rot6(1,1)*ca(1,1)+rot6(1,2)*ca(2,1)+rot6(1,3)*ca(3,1)
      rotb(1,2)=rot6(1,1)*ca(1,2)+rot6(1,2)*ca(2,2)+rot6(1,3)*ca(3,2)
      rotb(1,3)=rot6(1,1)*ca(1,3)+rot6(1,2)*ca(2,3)+rot6(1,3)*ca(3,3)
c
      rotb(2,1)=rot6(2,1)*ca(1,1)+rot6(2,2)*ca(2,1)+rot6(2,3)*ca(3,1)
      rotb(2,2)=rot6(2,1)*ca(1,2)+rot6(2,2)*ca(2,2)+rot6(2,3)*ca(3,2)
      rotb(2,3)=rot6(2,1)*ca(1,3)+rot6(2,2)*ca(2,3)+rot6(2,3)*ca(3,3)
c
      rotb(3,1)=rot6(3,1)*ca(1,1)+rot6(3,2)*ca(2,1)+rot6(3,3)*ca(3,1)
      rotb(3,2)=rot6(3,1)*ca(1,2)+rot6(3,2)*ca(2,2)+rot6(3,3)*ca(3,2)
      rotb(3,3)=rot6(3,1)*ca(1,3)+rot6(3,2)*ca(2,3)+rot6(3,3)*ca(3,3)
c
c Calculs preliminaires: tenseurs antisymetriques [omega(i)' X]
c
c vecteur rotation associe au tenseur
c C0T d/dt C0 --> omeg0'
c
c tomeg = c0t.dc0
c
      tomeg(1,3) = c0t(1,1)*dc0(1,3)+c0t(1,2)*dc0(2,3)+c0t(1,3)*dc0(3,3)
c
      tomeg(2,1) = c0t(2,1)*dc0(1,1)+c0t(2,2)*dc0(2,1)+c0t(2,3)*dc0(3,1)
c
      tomeg(3,2) = c0t(3,1)*dc0(1,2)+c0t(3,2)*dc0(2,2)+c0t(3,3)*dc0(3,2)
c
      omeg01 = tomeg(3,2)
      omeg02 = tomeg(1,3)
      omeg03 = tomeg(2,1)
c
c vecteur rotation associe au tenseur 
c CA [ omeg0 X] CAT + CA d/dt CAT --> omega' 
c
c tomeg = 0.
c
      tomeg(1,3) = 0.
c
      tomeg(2,1) = 0.
c
      tomeg(3,2) = 0.
c
      omega1 = omeg01*(cat(3,3)*cat(2,2)-cat(2,3)*cat(3,2))
     &        +omeg02*(cat(1,3)*cat(3,2)-cat(3,3)*cat(1,2))
     &        +omeg03*(cat(2,3)*cat(1,2)-cat(1,3)*cat(2,2))+tomeg(3,2)
      omega2 = omeg01*(cat(3,1)*cat(2,3)-cat(2,1)*cat(3,3))
     &        +omeg02*(cat(1,1)*cat(3,3)-cat(3,1)*cat(1,3))
     &        +omeg03*(cat(2,1)*cat(1,3)-cat(1,1)*cat(2,3))+tomeg(1,3)
      omega3 = omeg01*(cat(3,2)*cat(2,1)-cat(2,2)*cat(3,1))
     &        +omeg02*(cat(1,2)*cat(3,1)-cat(3,2)*cat(1,1))
     &        +omeg03*(cat(2,2)*cat(1,1)-cat(1,2)*cat(2,1))+tomeg(2,1)   
c
c vecteur rotation associe au tenseur 
c C1T [ omega X] C1 + C1T d/dt C1 --> omeg1' 
c
c tomeg = c1t.dc1
c
      tomeg(1,3) = c1t(1,1)*dc1(1,3)+c1t(1,2)*dc1(2,3)+c1t(1,3)*dc1(3,3)
c
      tomeg(2,1) = c1t(2,1)*dc1(1,1)+c1t(2,2)*dc1(2,1)+c1t(2,3)*dc1(3,1)
c
      tomeg(3,2) = c1t(3,1)*dc1(1,2)+c1t(3,2)*dc1(2,2)+c1t(3,3)*dc1(3,2)
c
      omeg11 = omega1*(c1(3,3)*c1(2,2)-c1(2,3)*c1(3,2))
     &        +omega2*(c1(1,3)*c1(3,2)-c1(3,3)*c1(1,2))
     &        +omega3*(c1(2,3)*c1(1,2)-c1(1,3)*c1(2,2))+tomeg(3,2)
      omeg12 = omega1*(c1(3,1)*c1(2,3)-c1(2,1)*c1(3,3))
     &        +omega2*(c1(1,1)*c1(3,3)-c1(3,1)*c1(1,3))
     &        +omega3*(c1(2,1)*c1(1,3)-c1(1,1)*c1(2,3))+tomeg(1,3)
      omeg13 = omega1*(c1(3,2)*c1(2,1)-c1(2,2)*c1(3,1))
     &        +omega2*(c1(1,2)*c1(3,1)-c1(3,2)*c1(1,1))
     &        +omega3*(c1(2,2)*c1(1,1)-c1(1,2)*c1(2,1))+tomeg(2,1)
c
c vecteur rotation associe au tenseur 
c C2T [ omeg1 X] C2 + C2T d/dt C2 --> omeg2' 
c
c tomeg = c2t.dc2
c
      tomeg(1,3) = c2t(1,1)*dc2(1,3)+c2t(1,2)*dc2(2,3)+c2t(1,3)*dc2(3,3)
c
      tomeg(2,1) = c2t(2,1)*dc2(1,1)+c2t(2,2)*dc2(2,1)+c2t(2,3)*dc2(3,1)
c
      tomeg(3,2) = c2t(3,1)*dc2(1,2)+c2t(3,2)*dc2(2,2)+c2t(3,3)*dc2(3,2)
c
      omeg21 = omeg11*(c2(3,3)*c2(2,2)-c2(2,3)*c2(3,2))
     &        +omeg12*(c2(1,3)*c2(3,2)-c2(3,3)*c2(1,2))
     &        +omeg13*(c2(2,3)*c2(1,2)-c2(1,3)*c2(2,2))+tomeg(3,2)
      omeg22 = omeg11*(c2(3,1)*c2(2,3)-c2(2,1)*c2(3,3))
     &        +omeg12*(c2(1,1)*c2(3,3)-c2(3,1)*c2(1,3))
     &        +omeg13*(c2(2,1)*c2(1,3)-c2(1,1)*c2(2,3))+tomeg(1,3)
      omeg23 = omeg11*(c2(3,2)*c2(2,1)-c2(2,2)*c2(3,1))
     &        +omeg12*(c2(1,2)*c2(3,1)-c2(3,2)*c2(1,1))
     &        +omeg13*(c2(2,2)*c2(1,1)-c2(1,2)*c2(2,1))+tomeg(2,1)
c
c vecteur rotation associe au tenseur 
c C3T [ omeg2 X] C3 + C3T d/dt C3 --> omeg3' 
c
c tomeg = c3t.dc3
c
      tomeg(1,3) = c3t(1,1)*dc3(1,3)+c3t(1,2)*dc3(2,3)+c3t(1,3)*dc3(3,3)
c
      tomeg(2,1) = c3t(2,1)*dc3(1,1)+c3t(2,2)*dc3(2,1)+c3t(2,3)*dc3(3,1)
c
      tomeg(3,2) = c3t(3,1)*dc3(1,2)+c3t(3,2)*dc3(2,2)+c3t(3,3)*dc3(3,2)
c
      omeg31 = omeg21*(c3(3,3)*c3(2,2)-c3(2,3)*c3(3,2))
     &        +omeg22*(c3(1,3)*c3(3,2)-c3(3,3)*c3(1,2))
     &        +omeg23*(c3(2,3)*c3(1,2)-c3(1,3)*c3(2,2))+tomeg(3,2)
      omeg32 = omeg21*(c3(3,1)*c3(2,3)-c3(2,1)*c3(3,3))
     &        +omeg22*(c3(1,1)*c3(3,3)-c3(3,1)*c3(1,3))
     &        +omeg23*(c3(2,1)*c3(1,3)-c3(1,1)*c3(2,3))+tomeg(1,3)
      omeg33 = omeg21*(c3(3,2)*c3(2,1)-c3(2,2)*c3(3,1))
     &        +omeg22*(c3(1,2)*c3(3,1)-c3(3,2)*c3(1,1))
     &        +omeg23*(c3(2,2)*c3(1,1)-c3(1,2)*c3(2,1))+tomeg(2,1)
c
c vecteur rotation associe au tenseur 
c C4T [ omeg3 X] C4 + C4T d/dt C4 --> omeg4' 
c
c tomeg = c4t.dc4
c
      tomeg(1,3) = c4t(1,1)*dc4(1,3)+c4t(1,2)*dc4(2,3)+c4t(1,3)*dc4(3,3)
c
      tomeg(2,1) = c4t(2,1)*dc4(1,1)+c4t(2,2)*dc4(2,1)+c4t(2,3)*dc4(3,1)
c
      tomeg(3,2) = c4t(3,1)*dc4(1,2)+c4t(3,2)*dc4(2,2)+c4t(3,3)*dc4(3,2)
c
      omeg41 = omeg31*(c4(3,3)*c4(2,2)-c4(2,3)*c4(3,2))
     &        +omeg32*(c4(1,3)*c4(3,2)-c4(3,3)*c4(1,2))
     &        +omeg33*(c4(2,3)*c4(1,2)-c4(1,3)*c4(2,2))+tomeg(3,2)
      omeg42 = omeg31*(c4(3,1)*c4(2,3)-c4(2,1)*c4(3,3))
     &        +omeg32*(c4(1,1)*c4(3,3)-c4(3,1)*c4(1,3))
     &        +omeg33*(c4(2,1)*c4(1,3)-c4(1,1)*c4(2,3))+tomeg(1,3)
      omeg43 = omeg31*(c4(3,2)*c4(2,1)-c4(2,2)*c4(3,1))
     &        +omeg32*(c4(1,2)*c4(3,1)-c4(3,2)*c4(1,1))
     &        +omeg33*(c4(2,2)*c4(1,1)-c4(1,2)*c4(2,1))+tomeg(2,1)
c
c vecteur rotation associe au tenseur 
c C5T [ omeg4 X] C5 + C5T d/dt C5 --> omeg5' 
c
c tomeg = c5t.dc5
c
      tomeg(1,3) = c5t(1,1)*dc5(1,3)+c5t(1,2)*dc5(2,3)+c5t(1,3)*dc5(3,3)
c
      tomeg(2,1) = c5t(2,1)*dc5(1,1)+c5t(2,2)*dc5(2,1)+c5t(2,3)*dc5(3,1)
c
      tomeg(3,2) = c5t(3,1)*dc5(1,2)+c5t(3,2)*dc5(2,2)+c5t(3,3)*dc5(3,2)
c
      omeg51 = omeg41*(c5(3,3)*c5(2,2)-c5(2,3)*c5(3,2))
     &        +omeg42*(c5(1,3)*c5(3,2)-c5(3,3)*c5(1,2))
     &        +omeg43*(c5(2,3)*c5(1,2)-c5(1,3)*c5(2,2))+tomeg(3,2)
      omeg52 = omeg41*(c5(3,1)*c5(2,3)-c5(2,1)*c5(3,3))
     &        +omeg42*(c5(1,1)*c5(3,3)-c5(3,1)*c5(1,3))
     &        +omeg43*(c5(2,1)*c5(1,3)-c5(1,1)*c5(2,3))+tomeg(1,3)
      omeg53 = omeg41*(c5(3,2)*c5(2,1)-c5(2,2)*c5(3,1))
     &        +omeg42*(c5(1,2)*c5(3,1)-c5(3,2)*c5(1,1))
     &        +omeg43*(c5(2,2)*c5(1,1)-c5(1,2)*c5(2,1))+tomeg(2,1)
c
c vecteur rotation associe au tenseur 
c C6T [ omeg5 X] C6 + C6T d/dt C6 --> omeg6' 
c
c tomeg = c6t.dc6
c
      tomeg(1,3) = c6t(1,1)*dc6(1,3)+c6t(1,2)*dc6(2,3)+c6t(1,3)*dc6(3,3)
c
      tomeg(2,1) = c6t(2,1)*dc6(1,1)+c6t(2,2)*dc6(2,1)+c6t(2,3)*dc6(3,1)
c
      tomeg(3,2) = c6t(3,1)*dc6(1,2)+c6t(3,2)*dc6(2,2)+c6t(3,3)*dc6(3,2)
c
      omeg61 = omeg51*(c6(3,3)*c6(2,2)-c6(2,3)*c6(3,2))
     &        +omeg52*(c6(1,3)*c6(3,2)-c6(3,3)*c6(1,2))
     &        +omeg53*(c6(2,3)*c6(1,2)-c6(1,3)*c6(2,2))+tomeg(3,2)
      omeg62 = omeg51*(c6(3,1)*c6(2,3)-c6(2,1)*c6(3,3))
     &        +omeg52*(c6(1,1)*c6(3,3)-c6(3,1)*c6(1,3))
     &        +omeg53*(c6(2,1)*c6(1,3)-c6(1,1)*c6(2,3))+tomeg(1,3)
      omeg63 = omeg51*(c6(3,2)*c6(2,1)-c6(2,2)*c6(3,1))
     &        +omeg52*(c6(1,2)*c6(3,1)-c6(3,2)*c6(1,1))
     &        +omeg53*(c6(2,2)*c6(1,1)-c6(1,2)*c6(2,1))+tomeg(2,1)
c
c vecteur rotation associe au tenseur 
c CAT [ omeg6 X] CA + CA d/dt CA --> omegb' 
c
c tomeg = 0.
c
      tomeg(1,3) = 0.
c
      tomeg(2,1) = 0.
c
      tomeg(3,2) = 0.
c
      omegb1 = omeg61*(ca(3,3)*ca(2,2)-ca(2,3)*ca(3,2))
     &        +omeg62*(ca(1,3)*ca(3,2)-ca(3,3)*ca(1,2))
     &        +omeg63*(ca(2,3)*ca(1,2)-ca(1,3)*ca(2,2))+tomeg(3,2)
      omegb2 = omeg61*(ca(3,1)*ca(2,3)-ca(2,1)*ca(3,3))
     &        +omeg62*(ca(1,1)*ca(3,3)-ca(3,1)*ca(1,3))
     &        +omeg63*(ca(2,1)*ca(1,3)-ca(1,1)*ca(2,3))+tomeg(1,3)
      omegb3 = omeg61*(ca(3,2)*ca(2,1)-ca(2,2)*ca(3,1))
     &        +omeg62*(ca(1,2)*ca(3,1)-ca(3,2)*ca(1,1))
     &        +omeg63*(ca(2,2)*ca(1,1)-ca(1,2)*ca(2,1))+tomeg(2,1)
c
c Calculs preliminaires: decalages entre reperes 
c                        dans le repere absolu
c
      x01 = x1-x0
      y01 = y1-y0
      z01 = z1-z0
c
      x12 = x2-x1
      y12 = y2-y1
      z12 = z2-z1
c
      x23 = x3-x2
      y23 = y3-y2
      z23 = z3-z2
c
      x34 = x4-x3
      y34 = y4-y3
      z34 = z4-z3
c
      x45 = x5-x4
      y45 = y5-y4
      z45 = z5-z4
c
      x56 = x6-x5
      y56 = y6-y5
      z56 = z6-z5
c
      x61 = x1-x6
      y61 = y1-y6
      z61 = z1-z6
c
c Sortie: origine du repere entraine dans RA
c
      r0(1) = x0+rot0(1,1)*x01+rot0(1,2)*y01+rot0(1,3)*z01+
     &           rot1(1,1)*x12+rot1(1,2)*y12+rot1(1,3)*z12+ 
     &           rot2(1,1)*x23+rot2(1,2)*y23+rot2(1,3)*z23+ 
     &           rot3(1,1)*x34+rot3(1,2)*y34+rot3(1,3)*z34+
     &           rot4(1,1)*x45+rot4(1,2)*y45+rot4(1,3)*z45+
     &           rot5(1,1)*x56+rot5(1,2)*y56+rot5(1,3)*z56+
     &           rot6(1,1)*x61+rot6(1,2)*y61+rot6(1,3)*z61+
     &           mux*time
      r0(2) = y0+rot0(2,1)*x01+rot0(2,2)*y01+rot0(2,3)*z01+
     &           rot1(2,1)*x12+rot1(2,2)*y12+rot1(2,3)*z12+ 
     &           rot2(2,1)*x23+rot2(2,2)*y23+rot2(2,3)*z23+ 
     &           rot3(2,1)*x34+rot3(2,2)*y34+rot3(2,3)*z34+
     &           rot4(2,1)*x45+rot4(2,2)*y45+rot4(2,3)*z45+
     &           rot5(2,1)*x56+rot5(2,2)*y56+rot5(2,3)*z56+
     &           rot6(2,1)*x61+rot6(2,2)*y61+rot6(2,3)*z61+
     &           muy*time
      r0(3) = z0+rot0(3,1)*x01+rot0(3,2)*y01+rot0(3,3)*z01+
     &           rot1(3,1)*x12+rot1(3,2)*y12+rot1(3,3)*z12+ 
     &           rot2(3,1)*x23+rot2(3,2)*y23+rot2(3,3)*z23+ 
     &           rot3(3,1)*x34+rot3(3,2)*y34+rot3(3,3)*z34+
     &           rot4(3,1)*x45+rot4(3,2)*y45+rot4(3,3)*z45+
     &           rot5(3,1)*x56+rot5(3,2)*y56+rot5(3,3)*z56+
     &           rot6(3,1)*x61+rot6(3,2)*y61+rot6(3,3)*z61+
     &           muz*time
c
c
c Sortie : matrice de passage RA/RR
c
      rot(1,1) = rotb(1,1)
      rot(1,2) = rotb(1,2)
      rot(1,3) = rotb(1,3)
      rot(2,1) = rotb(2,1)
      rot(2,2) = rotb(2,2)
      rot(2,3) = rotb(2,3)
      rot(3,1) = rotb(3,1)
      rot(3,2) = rotb(3,2)
      rot(3,3) = rotb(3,3)
c
c Sorties: vitesse du maillage dans le repere absolu
c
      v0x = rot0(1,1)*(omeg02*z01-omega3*y01)
     &     +rot0(1,2)*(omeg03*x01-omega1*z01)
     &     +rot0(1,3)*(omeg01*y01-omega2*x01)
      v0y = rot0(2,1)*(omeg02*z01-omega3*y01)
     &     +rot0(2,2)*(omeg03*x01-omega1*z01)
     &     +rot0(2,3)*(omeg01*y01-omega2*x01)
      v0z = rot0(3,1)*(omeg02*z01-omega3*y01)
     &     +rot0(3,2)*(omeg03*x01-omega1*z01)
     &     +rot0(3,3)*(omeg01*y01-omega2*x01)
c
      v1x = rot1(1,1)*(omeg12*z12-omeg13*y12)
     &     +rot1(1,2)*(omeg13*x12-omeg11*z12)
     &     +rot1(1,3)*(omeg11*y12-omeg12*x12)
      v1y = rot1(2,1)*(omeg12*z12-omeg13*y12)
     &     +rot1(2,2)*(omeg13*x12-omeg11*z12)
     &     +rot1(2,3)*(omeg11*y12-omeg12*x12)
      v1z = rot1(3,1)*(omeg12*z12-omeg13*y12)
     &     +rot1(3,2)*(omeg13*x12-omeg11*z12)
     &     +rot1(3,3)*(omeg11*y12-omeg12*x12)
c
      v2x = rot2(1,1)*(omeg22*z23-omeg23*y23)
     &     +rot2(1,2)*(omeg23*x23-omeg21*z23)
     &     +rot2(1,3)*(omeg21*y23-omeg22*x23)
      v2y = rot2(2,1)*(omeg22*z23-omeg23*y23)
     &     +rot2(2,2)*(omeg23*x23-omeg21*z23)
     &     +rot2(2,3)*(omeg21*y23-omeg22*x23)
      v2z = rot2(3,1)*(omeg22*z23-omeg23*y23)
     &     +rot2(3,2)*(omeg23*x23-omeg21*z23)
     &     +rot2(3,3)*(omeg21*y23-omeg22*x23)
c
      v3x = rot3(1,1)*(omeg32*z34-omeg33*y34)
     &     +rot3(1,2)*(omeg33*x34-omeg31*z34)
     &     +rot3(1,3)*(omeg31*y34-omeg32*x34)
      v3y = rot3(2,1)*(omeg32*z34-omeg33*y34)
     &     +rot3(2,2)*(omeg33*x34-omeg31*z34)
     &     +rot3(2,3)*(omeg31*y34-omeg32*x34)
      v3z = rot3(3,1)*(omeg32*z34-omeg33*y34)
     &     +rot3(3,2)*(omeg33*x34-omeg31*z34)
     &     +rot3(3,3)*(omeg31*y34-omeg32*x34)
c
      v4x = rot4(1,1)*(omeg42*z45-omeg43*y45)
     &     +rot4(1,2)*(omeg43*x45-omeg41*z45)
     &     +rot4(1,3)*(omeg41*y45-omeg42*x45)
      v4y = rot4(2,1)*(omeg42*z45-omeg43*y45)
     &     +rot4(2,2)*(omeg43*x45-omeg41*z45)
     &     +rot4(2,3)*(omeg41*y45-omeg42*x45)
      v4z = rot4(3,1)*(omeg42*z45-omeg43*y45)
     &     +rot4(3,2)*(omeg43*x45-omeg41*z45)
     &     +rot4(3,3)*(omeg41*y45-omeg42*x45)
c
      v5x = rot5(1,1)*(omeg52*z56-omeg53*y56)
     &     +rot5(1,2)*(omeg53*x56-omeg51*z56)
     &     +rot5(1,3)*(omeg51*y56-omeg52*x56)
      v5y = rot5(2,1)*(omeg52*z56-omeg53*y56)
     &     +rot5(2,2)*(omeg53*x56-omeg51*z56)
     &     +rot5(2,3)*(omeg51*y56-omeg52*x56)
      v5z = rot5(3,1)*(omeg52*z56-omeg53*y56)
     &     +rot5(3,2)*(omeg53*x56-omeg51*z56)
     &     +rot5(3,3)*(omeg51*y56-omeg52*x56)
c
      v6x = rot6(1,1)*(omeg62*z61-omeg63*y61)
     &     +rot6(1,2)*(omeg63*x61-omeg61*z61)
     &     +rot6(1,3)*(omeg61*y61-omeg62*x61)
      v6y = rot6(2,1)*(omeg62*z61-omeg63*y61)
     &     +rot6(2,2)*(omeg63*x61-omeg61*z61)
     &     +rot6(2,3)*(omeg61*y61-omeg62*x61)
      v6z = rot6(3,1)*(omeg62*z61-omeg63*y61)
     &     +rot6(3,2)*(omeg63*x61-omeg61*z61)
     &     +rot6(3,3)*(omeg61*y61-omeg62*x61)
c
c Attention: suppose le centre de la premiere rotation ind. du temps
c     
      v01 = v0x+v1x+v2x+v3x+v4x+v5x+v6x + mux
      v02 = v0y+v1y+v2y+v3y+v4y+v5y+v6y + muy
      v03 = v0z+v1z+v2z+v3z+v4z+v5z+v6z + muz
c
c Sortie de la vitesse d'entrainement dans le repere entraine:
c ------------------------------------------------------------
c
c      s0(1) = rotb(1,1)*v01+rotb(2,1)*v02+rotb(3,1)*v03
c      s0(2) = rotb(1,2)*v01+rotb(2,2)*v02+rotb(3,2)*v03
c      s0(3) = rotb(1,3)*v01+rotb(2,3)*v02+rotb(3,3)*v03

C
c      omeg(1) = omegb1
c      omeg(2) = omegb2
c      omeg(3) = omegb3
C

C CB :sortie de s0 et omega dans le repere absolu
      s0(1) = v01
      s0(2) = v02
      s0(3) = v03

      omeg(1) = rotb(1,1)*omegb1+rotb(1,2)*omegb2+rotb(1,3)*omegb3
      omeg(2) = rotb(2,1)*omegb1+rotb(2,2)*omegb2+rotb(2,3)*omegb3
      omeg(3) = rotb(3,1)*omegb1+rotb(3,2)*omegb2+rotb(3,3)*omegb3
C FIN CB

      xr(1) = x1
      xr(2) = y1
      xr(3) = z1
C
      END
C ============ Tbx/Motion/TbxEvalRotF.for ================
