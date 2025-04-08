# - exteriorFaces (pyTree) -
import KCore.test as test
import Post.PyTree as P
import Generator.PyTree as G

# 1D
a = G.cart((0,0,0), (1,1,1), (10,1,1))
b = P.exteriorFaces(a)
test.testT(b,1)

# 2D array
a = G.cart((0,0,0), (1,1,1), (10,6,1))
b = P.exteriorFaces(a)
test.testT(b,2)

# 3D array
a = G.cart((0,0,0), (1,1,1), (4,4,6))
b = P.exteriorFaces(a)
test.testT(b,3)

# TRI
a = G.cartTetra((0,0,0), (1,1,1), (20,3,1))
b = P.exteriorFaces(a)
test.testT(b,4)

# QUAD
a = G.cartHexa((0,0,0), (1,1,1), (20,3,1))
b = P.exteriorFaces(a)
test.testT(b,5)

# TETRA
a = G.cartTetra((0,0,0), (1,1,1), (3,3,3))
b = P.exteriorFaces(a)
test.testT(b,6)

# HEXA
a = G.cartHexa((0,0,0), (1,1,1), (3,3,3))
b = P.exteriorFaces(a)
test.testT(b,7)

# BAR
a = G.cartTetra((0,0,0), (1,1,1), (5,1,1))
b = P.exteriorFaces(a)
test.testT(b, 8)

# NGON3D
a = G.cartNGon((0,0,0), (1,1,1), (5,5,5))
b = P.exteriorFaces(a)
test.testT(b, 9)

# NGON2D
a = G.cartNGon((0,0,0), (1,1,1), (5,5,1))
b = P.exteriorFaces(a)
test.testT(b, 10)
