# - unsignNGonFaces (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# Structured
a = G.cart((0., 0., 0.), (1., 1., 1.), (5, 5, 5))
C._unsignNGonFaces(a)
test.testT(a, 1)

# Unstructured - BE
a = G.cartTetra((0., 0., 0.), (1., 1., 1.), (5, 5, 5))
C._unsignNGonFaces(a)
test.testT(a, 2)

# Unstructured - NGon
a = G.cartNGon((0., 0., 0.), (1., 1., 1.), (5, 5, 5), api=1)
C._unsignNGonFaces(a)
test.testT(a, 3)

# Unstructured - unsigned NGon v4
a = G.cartNGon((0., 0., 0.), (1., 1., 1.), (5, 5, 5), api=3)
C._unsignNGonFaces(a)
test.testT(a, 4)

# Unstructured - signed NGon v4
a = G.cartNGon((0., 0., 0.), (1., 1., 1.), (5, 5, 5), api=3)
C._signNGonFaces(a)
C._unsignNGonFaces(a)
test.testT(a, 5)
