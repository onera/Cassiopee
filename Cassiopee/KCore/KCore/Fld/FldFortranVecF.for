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

C =============================================================================
C  F90 version of FldArray::setallvalues()
C This subroutine copies a array of float to another one.
C           setallvaluesf ( N, array_out, array_in )
C
C  setallvaluesf_( _sizeLoc*_nfldLoc, _tab, table );
C
C  for (E_Int l=0; l<_sizeLoc*_nfldLoc; l++)
C  {
C    _tab[l] = table[l];
C  }
C =============================================================================
      SUBROUTINE k6setallvaluesf(length, lhs, rhs)

      IMPLICIT NONE
C_IN
      INTEGER_E           length
      REAL_E              rhs(length)
C_OUT
      REAL_E              lhs(length)

      INTEGER_E i
!$OMP PARALLEL PRIVATE(i)
!$OMP DO
      DO i = 1, length
        lhs(i) = rhs(i)
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END

C =============================================================================
C E_Int version of setallvaluesf 
C =============================================================================
      SUBROUTINE k6setallvaluesi(length, lhs, rhs)

      IMPLICIT NONE
C_IN
      INTEGER_E           length
      INTEGER_E           rhs(length)
C_OUT
      INTEGER_E           lhs(length)
C_LOC
      INTEGER_E i
!$OMP PARALLEL PRIVATE(i) 
!$OMP DO
      DO i = 1, length
        lhs(i) = rhs(i)
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END

C =============================================================================
C int version of setallvaluesf
C =============================================================================
      SUBROUTINE k6setallvaluesi4(length, lhs, rhs)

      IMPLICIT NONE
C_IN
      INTEGER_E length
      INTEGER*4 rhs(length)
C_OUT
      INTEGER*4 lhs(length)
C_LOCAL
      INTEGER_E i

      DO i = 1, length
        lhs(i) = rhs(i)
      END DO

      END

C =============================================================================
C Indirect copy (gather)
C =============================================================================
      SUBROUTINE k6fldgather(neq, length, szrhs, indic, rhs, lhs)

      IMPLICIT NONE
C_IN
      INTEGER_E           neq, length, szrhs
      INTEGER_E           indic(0:length-1)
      REAL_E              rhs(0:szrhs-1, neq)
C_OUT
      REAL_E              lhs(0:length-1, neq)
C_LOCAL
      INTEGER_E n, i

      DO n = 1, neq
         DO i = 0, length-1
            lhs(i,n) = rhs(indic(i), n)
         END DO
      END DO
      END

C =============================================================================
C Indirect copy (gather) with double indirection
C =============================================================================
      SUBROUTINE k6fldgather2(neq, length, szindic2, szrhs,
     &                      indic1, indic2, rhs, lhs)

      IMPLICIT NONE
C_IN
      INTEGER_E           neq, length, szindic2, szrhs
      INTEGER_E           indic1(0:length-1)
      INTEGER_E           indic2(0:szindic2-1)
      REAL_E              rhs(0:szrhs-1, neq)
C_OUT
      REAL_E              lhs(0:length-1, neq)
C_LOCAL
      INTEGER_E n, i

      DO n = 1, neq
         DO i = 0, length-1
            lhs(i,n) = rhs(indic2(indic1(i)), n)
         END DO
      END DO
      END

C =============================================================================
C Specialized Set method.
C This subroutine copies a part of an array to an other one.
C           setallvaluesfrompartf ( N1,N2, array_out, array_in )
C
C  setallvaluesfrompartf_( _sizeLoc*(end-begin+1),_sizeLoc*(begin-1), 
C  _tab, table );
C
C  for (E_Int l=0; l<_sizeLoc*(end-begin+1); l++)
C  {
C    _tab[l] = table[_sizeLoc*(begin-1+l];
C  }
C =============================================================================
C F90 version of FldArray::setallvaluesAt()
C This subroutine initializes an array with a constant value.
C           setallvaluesatf ( N, array, val )
C
C  setallvaluesatf_( _sizeLoc*_nfldLoc, _tab, val );
C
C  for (E_Int l=0; l<_sizeLoc*_nfldLoc; l++)
C  {
C    _tab[l] = val;
C  }
C =============================================================================
      SUBROUTINE k6setallvaluesatf(length, x, val)

      IMPLICIT NONE
C_IN
      INTEGER_E           length
      REAL_E              val
C_OUT
      REAL_E              x(1:length)

      INTEGER_E i
!$OMP PARALLEL PRIVATE(i)
!$OMP DO
      DO i = 1, length
        x(i) = val
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END

C =============================================================================
C F90 version of FldArray::operator= (const FldArrayNec& rhs)
C This subroutine initializes a subarray with a constant value.
C           operatoregalf ( array, first, last, val )
C
C  operatoregalf_( array , first , last , val );
C
C      for (i=first; i<last+1; i++)
C      {
C        array[i] = val;
C      }
C =============================================================================
      SUBROUTINE k6operatoregalf (x, first, last, val)

      IMPLICIT NONE
C_IN
      INTEGER_E           first, last
      REAL_E              val
C_OUT
      REAL_E              x(0:last-1)

      INTEGER_E i

      DO i = first, last-1
         x(i) = val
      END DO

      END

C =============================================================================
C F90 version of FldArray::operator+= (const FldArrayNec& y)
C This subroutine adds an array to an other one.
C           operatorplusf ( N, x, y )
C
C  operatorplusf_( _sizeLoc*_nfldLoc, _tab, y._tab );
C
C  for(i=0;i<_sizeLoc*_nfldLoc;i++)
C    _tab[i] += y._tab[i];
C =============================================================================
      SUBROUTINE k6operatorplusf(length, x, y)

      IMPLICIT NONE
C_IN
      INTEGER_E           length
      REAL_E              y(1:length)

C_INOUT
      REAL_E              x(1:length)

      INTEGER_E i
!$OMP PARALLEL PRIVATE(i)
!$OMP DO
      DO i = 1, length
         x(i) = x(i) + y(i)
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END

C =============================================================================
C F90 version of FldArray::operator+= (const FldArrayNec&)
C  This subroutine adds an array to an other one.
C           operatorsubf ( N, x, y )
C
C  operatorsubf_( _sizeLoc*_nfldLoc, _tab, y._tab );
C
C  for(i=0;i<_sizeLoc*_nfldLoc;i++)
C    _tab[i] -= y._tab[i];
C =============================================================================
      SUBROUTINE k6operatorsubf(length, x, y)

      IMPLICIT NONE
C_IN
      INTEGER_E           length
      REAL_E              y(1:length)

C_INOUT
      REAL_E              x(1:length)

      INTEGER_E i
!$OMP PARALLEL PRIVATE(i)
!$OMP DO
      DO i = 1, length
         x(i) = x(i) - y(i)
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END

C =============================================================================
C F90 optimization
C Linear interpolation between two arrays
C =============================================================================
      SUBROUTINE k6fldlib2a( tab, tab0 ,tab1, c0, c1, first, last)
      IMPLICIT NONE
C_IN
      INTEGER_E           first,last
      REAL_E              tab0(first:last)
      REAL_E              tab1(first:last)
      REAL_E              c0, c1

C_INOUT
      REAL_E              tab(first:last)

C_LOC
      INTEGER_E           i
!$OMP PARALLEL PRIVATE(i)
!$OMP DO
      DO i = first, last
        tab(i) = c0*tab0(i) + c1*tab1(i)
      END DO  
!$OMP END DO
!$OMP END PARALLEL
      END

C =============================================================================
C convert from REAL*4 (float) to REAL*8 (double)
C =============================================================================
      SUBROUTINE k6conv2double(size, rhs, lhs)
      IMPLICIT NONE
C_IN
      INTEGER_E size
      REAL*4    rhs(0:size-1)

C_OUT
      REAL*8    lhs(0:size-1)

C_LOCAL
      INTEGER_E i

      DO i = 0, size-1
        lhs(i) = rhs(i)
      END DO  

      END

C =============================================================================
C Convert from REAL*8 (double) to REAL*4 (float)
C =============================================================================
      SUBROUTINE k6conv2simple(size, rhs, lhs)
      IMPLICIT NONE
C_IN
      INTEGER_E size
      REAL*8    rhs(0:size-1)

C_OUT
      REAL*4    lhs(0:size-1)

C_LOCAL
      INTEGER_E i

      DO i = 0, size-1
         lhs(i) = rhs(i)
      END DO

      END

C =============================================================================
C Matrix-Vector product
C =============================================================================
      SUBROUTINE k6multmatvecf(lengthi, lenghtj, x, array)

      IMPLICIT NONE

C_IN
      INTEGER_E lengthi, lenghtj
      REAL_E array(lengthi)

C_INOUT
      REAL_E x(1:lengthi*lenghtj)

C_LOCAL
      INTEGER_E ind
      INTEGER_E i, j, t
C==============================================================================

      DO j = 1, lenghtj
         ind = (j-1)*lengthi
         DO i = 1, lengthi
            t = i+ind
            x(t) = array(i)*x(t)
         END DO
      END DO

      END
C =============================================================================
#ifdef ELSA_DOCUMENTATION
// ============================================================================
// @Name setallvaluesati
// @Memo F90 version of FldArrayI::setallvaluesAt()
/* @Text
*/
Fortran setallvaluesati();
#endif 
C------------------------------------------------------------------------------
      SUBROUTINE k6setallvaluesati(length, x, val)

      IMPLICIT NONE
C_IN
      INTEGER_E length
      INTEGER_E val

C_OUT
      INTEGER_E x(1:length)
      
C_LOCAL
      INTEGER_E i
!$OMP PARALLEL PRIVATE(i)
!$OMP DO
      DO i = 1, length
         x(i) = val
      END DO
!$OMP END DO
!$OMP END PARALLEL
      END
C ===================== End of file Fld/FldFortranVecF.for ===============
