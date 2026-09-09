# - exteriorFaces (pyTree) -
import KCore.test as test
import Post.PyTree as P
import Generator.PyTree as G
import Converter.PyTree as C

# -- STRUCT
# 1D
a = G.cart((0,0,0), (1,1,1), (10,1,1))
b = P.exteriorFaces(a)
test.testT(b, 1)

# 2D
a = G.cart((0,0,0), (1,1,1), (10,6,1))
b = P.exteriorFaces(a)
test.testT(b, 2)

# 3D
a = G.cart((0,0,0), (1,1,1), (4,4,6))
b = P.exteriorFaces(a)
test.testT(b, 3)

# -- BE
# TRI
a = G.cartTetra((0,0,0), (1,1,1), (20,3,1))
b = P.exteriorFaces(a)
test.testT(b, 4)

# QUAD
a = G.cartHexa((0,0,0), (1,1,1), (20,3,1))
b = P.exteriorFaces(a)
test.testT(b, 5)

# TETRA
a = G.cartTetra((0,0,0), (1,1,1), (3,3,3))
b = P.exteriorFaces(a)
test.testT(b, 6)

# HEXA
a = G.cartHexa((0,0,0), (1,1,1), (3,3,3))
b = P.exteriorFaces(a)
test.testT(b, 7)

# BAR
a = G.cartTetra((0,0,0), (1,1,1), (5,1,1))
b = P.exteriorFaces(a)
test.testT(b, 8)

# -- NGON
# NGON3D, api1
a = G.cartNGon((0,0,0), (1,1,1), (5,5,5), api=1)
b = P.exteriorFaces(a)
test.testT(b, 9)

# NGON2D, api1
a = G.cartNGon((0,0,0), (1,1,1), (5,5,1), api=1)
b = P.exteriorFaces(a)
test.testT(b, 10)

# NGON3D, api3
a = G.cartNGon((0,0,0), (1,1,1), (5,5,5), api=3)
b = P.exteriorFaces(a)
test.testT(b, 11)

# NGON2D, api3
a = G.cartNGon((0,0,0), (1,1,1), (5,5,1), api=3)
b = P.exteriorFaces(a)
test.testT(b, 12)

# -- ME
# 2D: tri-quad-tri
#          |
#         tri
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (5,10,1))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.2), (5,10,1))
c = G.cartTetra((0.8,0.,0.), (0.1,0.1,0.2), (5,10,1))
d = G.cartTetra((0.4,-0.9,0.), (0.1,0.1,0.2), (5,10,1))
a = C.mergeConnectivity([a, b, c, d], None)
b = P.exteriorFaces(a)
test.testT(b, 13)

# 3D: pyra - penta - hexa
a = G.cartPyra((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
b = G.cartPenta((0.4,0.,0.), (0.1,0.1,0.1), (5,5,5))
c = G.cartHexa((0.8,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.mergeConnectivity([a, b, c], None)
b = P.exteriorFaces(a)
test.testT(b, 14)
