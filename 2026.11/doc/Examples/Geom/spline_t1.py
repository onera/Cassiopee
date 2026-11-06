# - spline (array) -
import Geom as D
import Converter as C
import Generator as G
import KCore.test as test

# Spline 1D
a = D.polyline([(1.5,2,0.),(0,1.6,0.),(2,1.6,0.)])
b = D.spline(a, 3, N=200)
c = D.spline(a, 3, density=10.)
test.testA([b,c],1)

# Spline 2D
ni = 4; nj = 4
a = G.cart((0,0,0), (1,1,1), (ni,nj,1))

C.setValue(a, (1,1,1), [1.,1.,2.])
C.setValue(a, (1,2,1), [1.,2.,5.])
C.setValue(a, (1,3,1), [1.,3.,5.])
C.setValue(a, (1,4,1), [1.,4.,2.])
C.setValue(a, (2,1,1), [2.,1.,2.])
C.setValue(a, (2,2,1), [2.,2.,5.])
C.setValue(a, (2,3,1), [2.,3.,5.])
C.setValue(a, (2,4,1), [2.,4.,2.])
C.setValue(a, (3,1,1), [3.,1.,2.])
C.setValue(a, (3,2,1), [3.,2.,5.])
C.setValue(a, (3,3,1), [3.,3.,5.])
C.setValue(a, (3,4,1), [3.,4.,2.])
C.setValue(a, (4,1,1), [4.,1.,2.])
C.setValue(a, (4,2,1), [4.,2.,5.])
C.setValue(a, (4,3,1), [4.,3.,5.])
C.setValue(a, (4,4,1), [4.,4.,2.])
b = D.spline(a, 4, N=30, M=30)
c = D.spline(a, 4, density=10.)
test.testA([b,c],2)

test.writeCoverage(100)
