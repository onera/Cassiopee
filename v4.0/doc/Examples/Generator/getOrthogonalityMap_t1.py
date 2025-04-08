# - getOrthogonalityMap (array) -
import Generator as G
import Converter as C
import KCore.test as test

# Test 3D structure
a = G.cylinder((0.,0.,0.), 0.5, 1., 360., 0., 10., (50,50,10))
ac = C.node2Center(a)
ortho = G.getOrthogonalityMap(a)
ortho = C.addVars([ac,  ortho])
test.testA([ortho], 1)

# Test 2D structure
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
ac = C.node2Center(a)
ortho = G.getOrthogonalityMap(a)
ortho = C.addVars([ac,  ortho])
test.testA([ortho], 2)

# Test 1D structure
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (1,10,1))
ac = C.node2Center(a)
ortho = G.getOrthogonalityMap(a)
ortho = C.addVars([ac,  ortho])
test.testA([ortho], 3)
