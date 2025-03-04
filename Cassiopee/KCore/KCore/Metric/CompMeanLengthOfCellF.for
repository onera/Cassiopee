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
C Calcul d'une longueur caracteristique d'une cellule structuree ou tetra: 
C moyenne des longueurs des cotes de la cellule
C ============================================================================
      SUBROUTINE k6compmeanlengthofstructcell(ni, nj, nk, indA, 
     &                                        xt, yt, zt, meanl)
C
      IMPLICIT NONE
C
C==============================================================================
C_IN
      INTEGER_E ni, nj, nk
      INTEGER_E indA
      REAL_E xt(0:ni*nj*nk-1)
      REAL_E yt(0:ni*nj*nk-1)
      REAL_E zt(0:ni*nj*nk-1)
C_OUT
      REAL_E meanl
C_LOCAL
      INTEGER_E ninj, i, j ,k, alpha, beta, gamma
      INTEGER_E indB, indC, indD, indE, indF, indG, indH
      REAL_E d1, d2, d3, x21, y21, z21, x43, y43, z43, x65, y65, z65
      
C==============================================================================
C     
      ninj = ni*nj
        
C     Warning : here i,j,k start from 0
      k = indA/ninj
      j = indA/ni - k*nj
      i = indA - k*ninj - j*ni
      
C     Increments
      alpha = 1
      beta = 1
      gamma = 1

      if (i .eq. ni-1) then 
         i = i-1
      endif
      if (j .eq. nj-1) then
         j = j-1
      endif
      if (k .eq. nk-1) then 
         k = k-1
      endif

C     cas 2D                                
      if (nk .eq. 1) gamma = 0
      if (ni .eq. 1) alpha = 0
      if (nj .eq. 1) beta = 0

      indA = i + j * ni      + k * ninj
      indB = (i+alpha) + j * ni      + k * ninj !(i+1,j,k)
      indC = (i+alpha) + (j+beta)*ni + k * ninj !(i+1,j+1,k)
      indD =  i        + (j+beta)*ni + k * ninj !(i,j+1,k)
      
      indE = i         + j * ni      + (k+gamma) * ninj ! (i,j,k+1) 
      indF = (i+alpha) + j * ni      + (k+gamma) * ninj ! (i+1,j,k+1)
      indG = (i+alpha) + (j+beta)*ni + (k+gamma) * ninj ! (i+1,j+1,k+1)
      indH =  i        + (j+beta)*ni + (k+gamma) * ninj ! (i,j+1,k+1) 

C     premier centre : facette ABCD
      x21 = xt(indE) + xt(indF) + xt(indG) + xt(indH) - 
     &     (xt(indA) + xt(indB) + xt(indC) + xt(indD) ) 
      y21 = yt(indE) + yt(indF) + yt(indG) + yt(indH) - 
     &     (yt(indA) + yt(indB) + yt(indC) + yt(indD) ) 
      z21 = zt(indE) + zt(indF) + zt(indG) + zt(indH) - 
     &     (zt(indA) + zt(indB) + zt(indC) + zt(indD) )

      x43 = xt(indA) + xt(indD) + xt(indE) + xt(indH) - 
     &     (xt(indB) + xt(indC) + xt(indF) + xt(indG) ) 
      y43 = yt(indA) + yt(indD) + yt(indE) + yt(indH) - 
     &     (yt(indB) + yt(indC) + yt(indF) + yt(indG) ) 
      z43 = zt(indA) + zt(indD) + zt(indE) + zt(indH) - 
     &     (zt(indB) + zt(indC) + zt(indF) + zt(indG) ) 


      x65 = xt(indG) + xt(indH) + xt(indC) + xt(indD) -
     &     (xt(indA) + xt(indB) + xt(indF) + xt(indE) )
      y65 = yt(indG) + yt(indH) + yt(indC) + yt(indD) - 
     &     (yt(indA) + yt(indB) + yt(indF) + yt(indE) ) 
      z65 = zt(indG) + zt(indH) + zt(indC) + zt(indD) - 
     &     (zt(indA) + zt(indB) + zt(indF) + zt(indE) )

      d1 = sqrt(x21*x21 + y21*y21 + z21*z21)
      d2 = sqrt(x43*x43 + y43*y43 + z43*z43)
      d3 = sqrt(x65*x65 + y65*y65 + z65*z65)

      meanl = (d1+d2+d3)/12.D0

      END

C-----------------------------------------------------------------------------

      SUBROUTINE k6compmeanlengthoftetracell(npts,indA,indB,indC,indD, 
     &                                       xt, yt, zt, meanl)
                                     
C
      IMPLICIT NONE
# include "Def/DefFortranConst.h"
C
C_IN
      INTEGER_E npts            !  nb de pts 
      INTEGER_E indA, indB, indC, indD   ! sommets du tetraedre
      REAL_E xt(0:npts-1)
      REAL_E yt(0:npts-1)
      REAL_E zt(0:npts-1)
C_OUT
      REAL_E meanl
C_LOCAL
      REAL_E xA, xB, xC, xD
      REAL_E yA, yB, yC, yD
      REAL_E zA, zB, zC, zD
      REAL_E dAB, dAC, dDA, dBC, dBD, dCD 
      
C--------------------------------------------------------------------------
      xA = xt(indA)
      yA = yt(indA)
      zA = zt(indA)
      
      xB = xt(indB)
      yB = yt(indB)
      zB = zt(indB)

      xC = xt(indC)
      yC = yt(indC)
      zC = zt(indC)

      xD = xt(indD)
      yD = yt(indD)
      zD = zt(indD)
    
      dAB = sqrt((xB-xA)*(xB-xA) + (yB-yA)*(yB-yA) + (zB-zA)*(zB-zA))
      dAC = sqrt((xC-xA)*(xC-xA) + (yC-yA)*(yC-yA) + (zC-zA)*(zC-zA))
      dDA = sqrt((xA-xD)*(xA-xD) + (yA-yD)*(yA-yD) + (zA-zD)*(zA-zD))      
      dBC = sqrt((xC-xB)*(xC-xB) + (yC-yB)*(yC-yB) + (zC-zB)*(zC-zB))
      dBD = sqrt((xD-xB)*(xD-xB) + (yD-yB)*(yD-yB) + (zD-zB)*(zD-zB))
      dCD = sqrt((xD-xC)*(xD-xC) + (yD-yC)*(yD-yC) + (zD-zC)*(zD-zC))
      meanl =(dAB+dAC+dBC+dCD+dDA)/6.D0
      END
C================ Fortran/CompMeanLengthOfCellF.for ======================
