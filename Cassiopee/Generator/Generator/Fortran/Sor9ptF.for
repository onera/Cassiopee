c
      SUBROUTINE k6sor(x, y, n, m,
     $     a11,a12,a22,b11,b12,b22,c11,c12,c22,rhsx,rhsy,xerr)
c
c     iterates on a 9-point stencil in two independent variables, which
c     results from a pair of equations of the form:
c
c     a11*x11 + b11*x12 + c11*x22 + a12*y11 + b12*y12 + c12*y22 = rhsx
c     a12*x11 + b12*x12 + c12*x22 + a22*y11 + b22*y12 + c22*y22 = rhsy
c
      IMPLICIT NONE
      INTEGER_E n,m
      REAL_E x(n,m), y(n,m)
      REAL_E a11(n,m), a12(n,m), a22(n,m)
      REAL_E b11(n,m), b12(n,m), b22(n,m)
      REAL_E c11(n,m), c12(n,m), c22(n,m)
      REAL_E rhsx(n,m), rhsy(n,m)
      REAL_E omg

      REAL_E dx, dy, tdx, tdy, dxx, dxy, dyy, xerr
      INTEGER_E ncnt,i,j
      REAL_E cf11, cf12, cf22, x12, y12, x11, y11, x22, y22, rhs1, rhs2
      REAL_E xtmp, ytmp, den
      INTEGER_E ncntmax
c
      ncntmax = 20
      omg = 0.01D0
      dx = 1.D0/n
      dy = 1.D0/m
      tdx = 2.D0*dx
      tdy = 2.D0*dy
      dxx = dx*dx
      dxy = 4.D0*dx*dy
      dyy = dy*dy

      xerr = 0.D0
      ncnt = 0
 1    ncnt = ncnt+1
c
      do j = 2, m-1
         do i = 2, n-1
c
            cf11 = 2.D0*( a11(i,j)/dxx + c11(i,j)/dyy )
            cf12 = 2.D0*( a12(i,j)/dxx + c12(i,j)/dyy )
            cf22 = 2.D0*( a22(i,j)/dxx + c22(i,j)/dyy )
c
            x12 = ( x(i+1,j+1) - x(i+1,j-1) - x(i-1,j+1) + x(i-1,j-1) )/
     $           dxy
            y12 = ( y(i+1,j+1) - y(i+1,j-1) - y(i-1,j+1) + y(i-1,j-1) )/
     $           dxy
c
            x11 = ( x(i+1,j) + x(i-1,j) )/dxx
            y11 = ( y(i+1,j) + y(i-1,j) )/dxx
            x22 = ( x(i,j+1) + x(i,j-1) )/dyy
            y22 = ( y(i,j+1) + y(i,j-1) )/dyy
c
            rhs1 = a11(i,j)*x11 + a12(i,j)*y11 + b11(i,j)*x12 + 
     $           b12(i,j)*y12 + c11(i,j)*x22 + c12(i,j)*y22 - 
     $           rhsx(i,j)
            rhs2 = a12(i,j)*x11 + a22(i,j)*y11 + b12(i,j)*x12 + 
     $           b22(i,j)*y12 + c12(i,j)*x22 + c22(i,j)*y22 -
     $           rhsy(i,j)
c
            den = cf11*cf22-cf12*cf12
            IF (den.EQ.0.) THEN
               WRITE(*,*) cf11, cf22, cf12
               WRITE(*,*) 'den = 0'
               STOP
            ENDIF
            xtmp = (cf22*rhs1-cf12*rhs2)/den
            ytmp = (cf11*rhs2-cf12*rhs1)/den
            xtmp = x(i,j) + omg*( xtmp - x(i,j) )
            ytmp = y(i,j) + omg*( ytmp - y(i,j) )
            xerr = max( xerr, abs( xtmp-x(i,j) ) )
            xerr = max( xerr, abs( ytmp-y(i,j) ) )
            x(i,j) = xtmp
            y(i,j) = ytmp
         ENDDO
      ENDDO
c
      IF (ncnt.LT.ncntmax) go to 1
c
      RETURN
      END
