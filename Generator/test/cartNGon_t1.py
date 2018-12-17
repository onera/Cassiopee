# - cartNGon (array) -
import Generator as G
import KCore.test as test

# 1D
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (11,1,1))
test.testA([a],1)
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (1,11,1))
test.testA([a],2)
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (1,1,11))
test.testA([a],3)
# 2D
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (11,11,1))
test.testA([a],4)
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (11,1,11))
test.testA([a],5)
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (1,11,11))
test.testA([a],6)
# 3D
a = G.cartNGon((0.,0.,0.), (0.1,0.1,0.1), (11,11,11))
test.testA([a],7)
test.writeCoverage(100)
