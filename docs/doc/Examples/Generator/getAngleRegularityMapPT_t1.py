# - getAngleRegularityMap (pyTree) -
import Generator.PyTree as G
import Transform.PyTree as T
import KCore.test as test

# Test 1D structure
a = G.cart((0,0,0), (1,1,1), (10,1,1))
a = T.deformPoint(a, (0,0,0), (0.1,0.1,1.), 0.5, 0.4)
a = T.deformPoint(a, (5,0,0), (0.,1.,0.), 0.5, 0.5)
t = G.getAngleRegularityMap(a)
test.testT(t, 1)

# Test 2D structure
a = G.cart((0,0,0), (1,1,1), (10,10,1))
a = T.deformPoint(a, (0,0,0), (0.1,0.1,1.), 0.5, 0.4)
a = T.deformPoint(a, (5,5,0), (1.,1.,0.), 0.5, 0.5)
t = G.getAngleRegularityMap(a)
test.testT(t, 2)

# Test 3D structure
a = G.cart((0,0,0), (1,1,1), (10,10,10))
a = T.deformPoint(a, (0,0,0), (0.1,0.1,1.), 0.5, 0.4)
a = T.deformPoint(a, (5,5,9), (1.,1.,1.), 0.5, 0.5)
t = G.getAngleRegularityMap(a)
test.testT(t, 3)