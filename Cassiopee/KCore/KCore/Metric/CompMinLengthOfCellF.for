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
C Calcul d'une longueur caracteristique de la cellule : minimum  des
c longueurs des cotes de la cellule
C ============================================================================
      SUBROUTINE k6compminlengthofcell(ni, nj, nk, indA, 
     &                                  xt, yt, zt, minl)
C
      IMPLICIT NONE
# include "Def/DefFortranConst.h"
C
C==============================================================================
C_IN
      INTEGER_E ni, nj, nk
      INTEGER_E indA
      REAL_E xt(0:ni*nj*nk-1)
      REAL_E yt(0:ni*nj*nk-1)
      REAL_E zt(0:ni*nj*nk-1)
C_OUT
      REAL_E minl
C_LOCAL
      INTEGER_E ninj, i, j ,k, alpha, beta, gamma
      INTEGER_E indB, indC, indD, indE, indF, indG, indH
      REAL_E xA, xB, xC, xD, xE, xF, xG, xH
      REAL_E yA, yB, yC, yD, yE, yF, yG, yH
      REAL_E zA, zB, zC, zD, zE, zF, zG, zH
      REAL_E dAB, dCd, dEF, dGH, dAE, dBF, dDH, dCG, dAD, dBC, dEH, dFG

C==============================================================================
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
      if (ni .eq. 1) alpha = 0
      if (nj .eq. 1) beta = 0
      if (nk .eq. 1) gamma = 0

      indA = i + j * ni      + k * ninj
      indB = (i+alpha) + j * ni      + k * ninj !(i+1,j,k)
      indC = (i+alpha) + (j+beta)*ni + k * ninj !(i+1,j+1,k)
      indD =  i        + (j+beta)*ni + k * ninj !(i,j+1,k)
      
      indE = i         + j * ni      + (k+gamma) * ninj ! (i,j,k+1) 
      indF = (i+alpha) + j * ni      + (k+gamma) * ninj ! (i+1,j,k+1)
      indG = (i+alpha) + (j+beta)*ni + (k+gamma) * ninj ! (i+1,j+1,k+1)
      indH =  i        + (j+beta)*ni + (k+gamma) * ninj ! (i,j+1,k+1) 

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

      xE = xt(indE)
      yE = yt(indE)
      zE = zt(indE)

      xF = xt(indF)
      yF = yt(indF)
      zF = zt(indF)
      
      xG = xt(indG)
      yG = yt(indG)
      zG = zt(indG)

      xH = xt(indH)
      yH = yt(indH)
      zH = zt(indH)

      if (alpha .eq. 0) then
         dAB = MAXFLOAT
         dCD = MAXFLOAT
         dEF = MAXFLOAT
         dGH = MAXFLOAT
      else
         dAB = sqrt((xB-xA)*(xB-xA) + (yB-yA)*(yB-yA) + (zB-zA)*(zB-zA))
         dCD = sqrt((xD-xC)*(xD-xC) + (yD-yC)*(yD-yC) + (zD-zC)*(zD-zC))
         dEF = sqrt((xF-xE)*(xF-xE) + (yF-yE)*(yF-yE) + (zF-zE)*(zF-zE))
         dGH = sqrt((xH-xG)*(xH-xG) + (yH-yG)*(yH-yG) + (zH-zG)*(zH-zG))
      endif
      
      if (beta .eq. 0) then
         dAD = MAXFLOAT
         dBC = MAXFLOAT
         dEH = MAXFLOAT
         dFG = MAXFLOAT
      else
         dAD = sqrt((xD-xA)*(xD-xA) + (yD-yA)*(yD-yA) + (zd-zA)*(zD-zA))
         dBC = sqrt((xC-xB)*(xC-xB) + (yC-yB)*(yC-yB) + (zC-zB)*(zC-zB))
         dEH = sqrt((xH-xE)*(xH-xE) + (yH-yE)*(yH-yE) + (zH-zE)*(zH-zE))
         dFG = sqrt((xG-xF)*(xG-xF) + (yG-yF)*(yG-yF) + (zG-zF)*(zG-zF))
      endif

      if (gamma .eq. 0) then
         dAE = MAXFLOAT
         dBF = MAXFLOAT
         dDH = MAXFLOAT
         dCG = MAXFLOAT
      else
         dAE = sqrt((xE-xA)*(xE-xA) + (yE-yA)*(yE-yA) + (zE-zA)*(zE-zA))
         dBF = sqrt((xF-xB)*(xF-xB) + (yF-yB)*(yF-yB) + (zF-zB)*(zF-zB))
         dDH = sqrt((xH-xD)*(xH-xD) + (yH-yD)*(yH-yD) + (zH-zD)*(zH-zD))
         dCG = sqrt((xG-xC)*(xG-xC) + (yG-yC)*(yG-yC) + (zG-zC)*(zG-zC))
      endif
      
      minl = MIN(dAB, dCD)
      minl = MIN(minl, dEF)
      minl = MIN(minl, dGH)
      minl = MIN(minl, dAD)
      minl = MIN(minl, dBC)
      minl = MIN(minl, dEH)
      minl = MIN(minl, dFG)
      minl = MIN(minl, dAE)
      minl = MIN(minl, dBF)
      minl = MIN(minl, dDH)
      minl = MIN(minl, dCG)
      END

C-----------------------------------------------------------------------------

      SUBROUTINE k6compminlengthoftetracell(npts,indA,indB,indC,indD, 
     &                                      xt, yt, zt, minl)

      IMPLICIT NONE
C_IN
      INTEGER_E npts            !  nb de pts 
      INTEGER_E indA, indB, indC, indD   ! sommets du tetraedre
      REAL_E xt(0:npts-1)
      REAL_E yt(0:npts-1)
      REAL_E zt(0:npts-1)
C_OUT
      REAL_E minl
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
   
      minl = MIN(dAB, dAC)
      minl = MIN(minl, dDA)
      minl = MIN(minl, dBC)
      minl = MIN(minl, dBD)
      minl = MIN(minl, dCD)
      END
C===================== KCore/Metric/CompMinLengthOfCellF.for ===============
