# - close (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import KCore.test as test

# -- ME
a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
b = G.cartHexa((9,0,0), (1,1,1), (10,10,10))
a = C.mergeConnectivity(a, b)
a = G.close(a, 1.e-6)
test.testT(a, 1)

# -- NGON v3
# NGON 1D
a = G.cartNGon((0,0,0), (1,1,1), (10,1,1), api=1)
a = G.close(a, 1.e-6)
test.testT(a, 2)

# NGON 2D
a = G.cartNGon((0,0,0), (1,1,1), (10,10,1), api=1)
a = G.close(a, 1.e-6)
test.testT(a, 3)

# NGON 3D
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=1)
a = G.close(a, 1.e-6)
test.testT(a, 4)

# -- NGON v4
# NGON 3D
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=3)
a = G.close(a, 1.e-6)
test.testT(a, 5)
