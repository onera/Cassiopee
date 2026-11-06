# - exteriorElts (array) -
import Post as P
import Generator as G
import Geom as D
import KCore.test as test

# Struct
a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = P.exteriorElts(a)
test.testA([b], 1)

# Hexa
a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
b = P.exteriorElts(a)
test.testA([b], 2)

# Tetra
a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
b = P.exteriorElts(a)
test.testA([b], 3)

# QUAD
a = G.cart((0,0,0), (1,1,1), (10,10,1))
b = P.exteriorElts(a)
test.testA([b], 4)

# TRI
a = G.cartTetra((0,0,0), (1,1,1), (10,10,1))
b = P.exteriorElts(a)
test.testA([b], 5)

# BAR
a = G.cartTetra((0,0,0), (1,1,1), (10,1,1))
b = P.exteriorElts(a)
test.testA([b], 6)

# Octree
s = D.sphere((0,0,0), 1., N=100); snear = 0.5
a = G.octree([s], [snear], dfar=5.)
b = P.exteriorElts(a)
test.testA([b], 7)

# NGon
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
b = P.exteriorElts(a)
test.testA([b], 8)

# Penta
a = G.cartPenta((0,0,0), (1,1,1), (10,10,10))
b = P.exteriorElts(a)
test.testA([b], 9)

# Pyra
a = G.cartPyra((0,0,0), (1,1,1), (10,10,10))
b = P.exteriorElts(a)
test.testA([b], 10)
