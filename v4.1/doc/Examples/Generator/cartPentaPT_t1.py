# - cartPenta (pyTree) -
import Generator.PyTree as G
import KCore.test as test

# 3D seulement
a = G.cartPenta((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
test.testT(a,1)
test.writeCoverage(100)
