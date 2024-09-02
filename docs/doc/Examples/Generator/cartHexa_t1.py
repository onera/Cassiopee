# - cartHexa (array) -
import Generator as G
import KCore.test as test

# BAR
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,1,1))
test.testA([a], 1)

# QUAD
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
test.testA([a], 2)

# HEXA
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
test.testA([a], 3)

test.writeCoverage(100)
