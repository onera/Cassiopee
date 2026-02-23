# - cartTetra (array) -
import Generator as G
import KCore.test as test

# BAR
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,1,1))
test.testA([a], 1)
# TRI
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
test.testA([a], 2)
# TETRA
a = G.cartTetra((0.,0.,0.), (1.,1.,1.), (2,2,2))
test.testA([a], 3)
test.writeCoverage(100)
