# - cartTetra (pyTree) -
import Generator.PyTree as G
import KCore.test as test

# BAR
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,1,1))
test.testT(a, 1)
# TRI
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
test.testT(a, 2)
# TETRA
a = G.cartTetra((0.,0.,0.), (1.,1.,1.), (2,2,2))
test.testT(a, 3)
test.writeCoverage(100)
