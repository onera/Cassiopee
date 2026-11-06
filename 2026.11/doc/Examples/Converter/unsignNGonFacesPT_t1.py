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

# pyTree with several zones, signed input, restore signness
a = G.cart((0,0,0), (1,1,1), (5,5,5))
b = G.cartNGon((4,0,0), (1,1,1), (5,5,5), api=3)
c = G.cartTetra((8,0,0), (1,1,1), (5,5,5))
t = C.newPyTree(['Base', 3]); t[2][1][2].extend([a, b, c])
C._signNGonFaces(t)
C._unsignNGonFaces(t)
C._signNGonFaces(t, force=False)
test.testT(t, 6)

# pyTree with several zones, unsigned input, restore signness
a = G.cart((0,0,0), (1,1,1), (5,5,5))
b = G.cartNGon((4,0,0), (1,1,1), (5,5,5), api=3)
c = G.cartTetra((8,0,0), (1,1,1), (5,5,5))
t = C.newPyTree(['Base', 3]); t[2][1][2].extend([a, b, c])
C._unsignNGonFaces(t)
C._signNGonFaces(t, force=False)
test.testT(t, 7)

# pyTree with several zones, unsigned input, force sign
a = G.cart((0,0,0), (1,1,1), (5,5,5))
b = G.cartNGon((4,0,0), (1,1,1), (5,5,5), api=3)
c = G.cartTetra((8,0,0), (1,1,1), (5,5,5))
t = C.newPyTree(['Base', 3]); t[2][1][2].extend([a, b, c])
t = C.unsignNGonFaces(t)
t = C.signNGonFaces(t, force=True)
test.testT(t, 8)