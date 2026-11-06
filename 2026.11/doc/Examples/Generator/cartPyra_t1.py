# - cartPyra (array) -
import Generator as G
import KCore.test as test

# Ne fonctionne qu'en 3D
a = G.cartPyra((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
test.testA([a], 1)

test.writeCoverage(100)
