# - getAngleRegularityMap (array) -
import Converter as C
import Generator as G
import Transform as T
import KCore.test as test

# Test 1D structure
a = G.cart((0,0,0), (1,1,1), (10,1,1))
a = T.deformPoint(a, (0,0,0), (0.1,0.1,1.), 0.5, 0.4)
a = T.deformPoint(a, (5,0,0), (0.,1.,0.), 0.5, 0.5)
ac = C.node2Center(a)
reg = G.getAngleRegularityMap(a)
reg = C.addVars([ac,  reg])
test.testA([reg], 1)

# Test 2D structure
a = G.cart((0,0,0), (1,1,1), (10,10,1))
a = T.deformPoint(a, (0,0,0), (0.1,0.1,1.), 0.5, 0.4)
a = T.deformPoint(a, (5,5,0), (1.,1.,0.), 0.5, 0.5)
ac = C.node2Center(a)
reg = G.getAngleRegularityMap(a)
reg = C.addVars([ac,  reg])
test.testA([reg], 2)

# Test 3D structure
a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = T.deformPoint(a, (0,0,0), (0.1,0.1,1.), 0.5, 0.4)
a = T.deformPoint(a, (5,5,9), (1.,1.,1.), 0.5, 0.5)
ac = C.node2Center(a)
reg = G.getAngleRegularityMap(a)
reg = C.addVars([ac,  reg])
test.testA([reg], 3)
