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
C     -------------------------
C     Integer Fortran constants
C     -------------------------
      INTEGER_E ZERO_I
      INTEGER_E ONE_I
      INTEGER_E TWO_I
      INTEGER_E THREE_I
      INTEGER_E FOUR_I
      INTEGER_E E_MAXITERNWT

      INTEGER_E X_AXIS_E, Y_AXIS_E, Z_AXIS_E

      PARAMETER (ZERO_I = 0)
      PARAMETER (ONE_I  = 1)
      PARAMETER (TWO_I  = 2)
      PARAMETER (THREE_I = 3)
      PARAMETER (FOUR_I = 4)
      PARAMETER (E_MAXITERNWT = 30)

      PARAMETER (X_AXIS_E=1, Y_AXIS_E=2, Z_AXIS_E=3)

C     ----------------------
C     Real Fortran constants
C     ----------------------
      REAL_E ZERO 
      REAL_E ONE_FOURTH
      REAL_E ONE_EIGHT
      REAL_E ONE_HALF
      REAL_E ONE
      REAL_E TWO
      REAL_E THREE
      REAL_E FOUR
      REAL_E E_CUTOFF
      REAL_E E_MIN_SURFACE
      REAL_E E_MIN_VOL
      REAL_E E_PI
      REAL_E E_MAXEXP
      REAL_E INFTYS
      REAL_E MAXFLOAT

#ifdef E_DOUBLEREAL
      PARAMETER (ZERO         =  0.0D0  )
      PARAMETER (ONE_FOURTH   =  0.25D0 )
      PARAMETER (ONE_EIGHT    =  0.125D0 )
      PARAMETER (ONE_HALF     =  0.5D0  )
      PARAMETER (ONE          =  1.0D0  )
      PARAMETER (TWO          =  2.0D0  )
      PARAMETER (THREE        =  3.0D0  )
      PARAMETER (FOUR         =  4.0D0  )
C     Useful 
      PARAMETER (E_CUTOFF     = 1.0D-11)
      PARAMETER (E_MIN_SURFACE= 1.0D-30)
      PARAMETER (E_MIN_VOL    = 1.0D-30)
      PARAMETER (E_PI         = 3.14159265359)
C     Avoid exponential overflow
      PARAMETER (E_MAXEXP     = 50.D0)
      PARAMETER (INFTYS       = 1.0D-30 )
      PARAMETER (MAXFLOAT     = 1.0E+30 )
#else
      PARAMETER (ZERO         = 0.0E0  )
      PARAMETER (ONE_EIGHT    = 0.125E0 )
      PARAMETER (ONE_FOURTH   = 0.25E0 )
      PARAMETER (ONE_HALF     = 0.5E0  )
      PARAMETER (ONE          = 1.0E0  )
      PARAMETER (TWO          = 2.0E0  )
      PARAMETER (THREE        = 3.0E0  )
      PARAMETER (FOUR         = 4.0E0  )
      PARAMETER (E_CUTOFF     = 1.0E-06)
      PARAMETER (E_MIN_SURFACE= 1.0E-30)
      PARAMETER (E_MIN_VOL    = 1.0D-30)
      PARAMETER (E_PI         = 3.14159265359)
      PARAMETER (E_MAXEXP     = 50.E0)
      PARAMETER (INFTYS       =  1.0E-30 )
      PARAMETER (MAXFLOAT     =  1.0E+30 )
#endif

c ===== Def/DefFortranConst.h === Last line ===
