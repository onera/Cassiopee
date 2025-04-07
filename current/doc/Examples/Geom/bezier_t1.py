# - bezier -
import Geom as D
import Converter as C
import Generator as G
import KCore.test as test

# 1D
pts = C.array('x,y,z', 7, 1, 1)
x = pts[1][0]; y = pts[1][1]; z = pts[1][2]
x[0] = 6.; x[1] = 5.4; x[2] = 4.8; x[3] = 2.5; x[4] = 0.3
y[0] = 0.01; y[1] = 0.036; y[2] = 0.064; y[3] = 0.21; y[4] = 0.26; y[5] = 0.047
z[0] = 1.; z[1] = 1.; z[2] = 1.; z[3] = 1.; z[4] = 1.; z[5] = 1.; z[6] = 1.

a = D.bezier(pts)
test.testA([a],1)

# 2D
ni = 2; nj = 3
a = G.cart((0,0,0), (1,1,1), (ni,nj,1))
C.setValue(a, (1,1,1), [1.,1.,2.])
C.setValue(a, (1,2,1), [1.,2.,5.])
C.setValue(a, (1,3,1), [1.,3.,2.])
C.setValue(a, (2,1,1), [2.,1.,2.])
C.setValue(a, (2,2,1), [2.,2.,5.])
C.setValue(a, (2,3,1), [2.,3.,2.])
b = D.bezier(a, 10, 10)
test.testA([b],2)
test.writeCoverage(100)
