# - resignNGonFaces (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Unstructured - NGon
a = G.cartNGon((0., 0., 0.), (1., 1., 1.), (5, 5, 5), api=1)
C._unsignNGonFaces(a)
C._resignNGonFaces(a)
test.testT(a, 1)

# Unstructured - unsigned NGon v4
a = G.cartNGon((0., 0., 0.), (1., 1., 1.), (5, 5, 5), api=3)
C._unsignNGonFaces(a)
C._resignNGonFaces(a)
test.testT(a, 2)

# Unstructured - signed NGon v4
a = G.cartNGon((0., 0., 0.), (1., 1., 1.), (5, 5, 5), api=3)
C._signNGonFaces(a)
C._unsignNGonFaces(a)
C._resignNGonFaces(a)
test.testT(a, 3)
