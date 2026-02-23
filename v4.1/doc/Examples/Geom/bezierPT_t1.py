# - bezier (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Bezier 1D
pts = D.polyline([(0.,0.,0.), (0.,1.,0.), (2.,1.,0.), (2.,0.,0.),\
                  (4.,-1.,0.), (5.,6.,0.),])
a = D.bezier(pts, 100); a[0] = 'bezier'
t = C.newPyTree(['Base',1]); t[2][1][2].append(a)
test.testT(t, 1)

# 2D
ni = 2; nj = 3
a = G.cart((0,0,0), (1,1,1), (ni,nj,1))
C.setValue(a,'GridCoordinates', (1,1,1), [1.,1.,2.])
C.setValue(a,'GridCoordinates', (1,2,1), [1.,2.,5.])
C.setValue(a,'GridCoordinates', (1,3,1), [1.,3.,2.])
C.setValue(a,'GridCoordinates', (2,1,1), [2.,1.,2.])
C.setValue(a,'GridCoordinates', (2,2,1), [2.,2.,5.])
C.setValue(a,'GridCoordinates', (2,3,1), [2.,3.,2.])
b = D.bezier(a, 10, 10)
test.testT([b],2)
test.writeCoverage(100)
