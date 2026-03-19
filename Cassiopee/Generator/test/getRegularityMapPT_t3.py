# - getRegularityMap (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Geom.PyTree as D
import Transform.PyTree as T
import KCore.test as test

h = 0.2 # minimum spacing
r = 1.1 # 10% stretching

# TEST 1D BAR with constant stretching
a = [(0.+h*r**i,0.,0.) for i in range(10)]
a = D.polyline(a)
a = C.convertArray2Hexa(a)
G._getRegularityMap(a)
test.testT(a, 1)

# Test 2D QUAD with constant stretching
a = G.cartRx3((-1,-1,0), (1,1,0), (h,h,1), (-5,-5,0), (5,5,0), (r,r,1), dim=2)
a = C.convertArray2Hexa(a)
G._getRegularityMap(a)
test.testT(a, 2)

# Test 3D HEXA with constant stretching
a = G.cartRx3((-1,-1,-1), (1,1,1), (h,h,h), (-5,-5,-5), (5,5,5), (r,r,r), dim=3)
a = C.convertArray2Hexa(a)
G._getRegularityMap(a)
test.testT(a, 3)