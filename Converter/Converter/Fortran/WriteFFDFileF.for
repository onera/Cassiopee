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

C  ============================================================================
C  Write the files zone_n_ffd72 according to FFD72 file format         
C  ============================================================================
      SUBROUTINE writeffdfilef(nd, kod_int,nnode, nelmt
     & , real_state1, real_state2
     & , coord_x, coord_y, coord_z
     & , Density, MomentumX, MomentumY, MomentumZ, Energy
     & , Mutsmu, TurbulentDistance
C    & , nlimt, Density_l, MomentX_l, MomentY_l, MomentZ_l
C    & , Energy_l, Mutsmu_l
     & , nlimt, Var_l
     & , nmatch, ielmtmtch1, ielmtmtch2, izoneznzn, ielmtznzn
     & , npoint, nodmtch, kpoinmtch) 
      IMPLICIT NONE
      INTEGER_E nd, nnode, nelmt, nlimt, kod_int(4)
      REAL_E real_state1(10), real_state2(4)
      REAL_E coord_x(nnode), coord_y(nnode), coord_z(nnode)
      REAL_E Density(nelmt)
      REAL_E MomentumX(nelmt), MomentumY(nelmt), MomentumZ(nelmt)
      REAL_E Energy(nelmt), MutsMu(nelmt), TurbulentDistance(nelmt)
      REAL_E Var_l(6*nlimt)
C     REAL_E MomentX_l(nlimt), MomentY_l(nlimt), MomentZ_l(nlimt)
C     REAL_E Energy_l(nlimt), MutsMu_l(nlimt)
      INTEGER_E nmatch, ielmtmtch1(nmatch), ielmtmtch2(nmatch)
      INTEGER_E npoint, nodmtch(npoint), kpoinmtch(nmatch+1)
C     Tableau de connectivité en multizone
      INTEGER_E ielmtznzn(nlimt), izoneznzn(nlimt)
C
      LOGICAL CellCenter, Rans
C
      INTEGER_E n
      CHARACTER *40 kfilezone
      WRITE(kfilezone,'(A,I5.5,A)')'FFD72_',nd,'.bin'
      print *, "ecriture du fichier ",trim(kfilezone)
      print *, "nnode=",nnode," nelmt=",nelmt
C     print *, "coord_x(1), coord_y(1), coord_z(1)"
C     print *,  coord_x(1), coord_y(1), coord_z(1) 
C     print *, "coord_x(nnode), coord_y(nnode), coord_z(nnode)"
C     print *,  coord_x(nnode), coord_y(nnode), coord_z(nnode) 
C     print *, "MomentumX(1), MomentumY(1), MomentumZ(1)"
C     print *,  MomentumX(1), MomentumY(1), MomentumZ(1) 
C     print *, "MomentumX(nelmt), MomentumY(nelmt), MomentumZ(nelmt)"
C     print *,  MomentumX(nelmt), MomentumY(nelmt), MomentumZ(nelmt) 
C     print *, "Density(1), Energy(1), MutsMu(1)"
C     print *,  Density(1), Energy(1), MutsMu(1) 
C     print *, "Density(nelmt), Energy(nelmt), MutsMu(nelmt)"
C     print *,  Density(nelmt), Energy(nelmt), MutsMu(nelmt) 
C
C     print *, "MomentX_l(1), MomentY_l(1), MomentZ_l(1)"
C     print *,  MomentX_l(1), MomentY_l(1), MomentZ_l(1) 
C     print *, "MomentX_l(nlimt), MomentY_l(nlimt), MomentZ_l(nlimt)"
C     print *,  MomentX_l(nlimt), MomentY_l(nlimt), MomentZ_l(nlimt) 
C     print *, "Density_l(1), Energy_l(1), MutsMu_l(1)"
C     print *,  Density_l(1), Energy_l(1), MutsMu_l(1) 
C     print *, "Density_l(nlimt), Energy_l(nlimt), MutsMu_l(nlimt)"
C     print *,  Density_l(nlimt), Energy_l(nlimt), MutsMu_l(nlimt) 
C
C  int_kode = (kodcc,kodsolver,kodnst,kod2d )
C==============================================================================
C
      CELLCENTER = kod_int(1).GT.0     
      Rans       = kod_int(3).GT.0     
      print *, kod_int
      print *, real_state1
      print *, real_state2
      OPEN(72,file=trim(kfileZone),form="unformatted"
     &       ,convert='big_endian')
      WRITE(72) kod_int
      WRITE(72) real_state1
      WRITE(72) real_state2
      WRITE(72) nnode, npoint, nmatch, nlimt, nelmt 
      WRITE(72) nodmtch
      WRITE(72) kpoinmtch
      WRITE(72) ielmtmtch1
      WRITE(72) ielmtmtch2
      WRITE(72) izoneznzn
      WRITE(72) ielmtznzn
      WRITE(72)(coord_x(n),coord_y(n),coord_z(n),n=1,nnode)
      IF(Rans)THEN
      WRITE(72) (Density(n), MomentumX(n), MomentumY(n), MomentumZ(n),
     &        Energy(n), MutsMu(n),n=1,nelmt)
      ELSE
      WRITE(72) (Density(n), MomentumX(n), MomentumY(n), MomentumZ(n),
     &        Energy(n),n=1,nelmt)
      ENDIF
      WRITE(72)Var_l
      if(Rans)WRITE(72)TurbulentDistance
      CLOSE(72)
      return
      END
