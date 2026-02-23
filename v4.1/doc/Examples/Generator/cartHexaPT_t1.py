# - cartHexa (pyTree) -
import Generator.PyTree as G
import KCore.test as test

# 1D
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,1,1))
test.testT(a, 1)
# 2D
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
test.testT(a, 2)
# 3D
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
test.testT(a, 3)
test.writeCoverage(100)
