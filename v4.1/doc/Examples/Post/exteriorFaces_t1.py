# - exteriorFaces (array) -
import KCore.test as test
import Post as P
import Generator as G
import Geom as D
import Converter as C

# 1D array
a = G.cart((0,0,0), (1,1,1), (10,1,1))
b = P.exteriorFaces(a)
test.testA([b],1)
indices = []
b = P.exteriorFaces(a,indices)
test.testO([b,indices],21)

# 2D array
a = G.cart((0,0,0), (1,1,1), (10,6,1))
b = P.exteriorFaces(a)
test.testA([b],2)
indices = []
b = P.exteriorFaces(a,indices)
test.testO([b,indices],22)
# 3D array
a = G.cart((0,0,0), (1,1,1), (4,4,6))
b = P.exteriorFaces(a)
test.testA([b],3)
indices = []
a = G.cart((0,0,0), (1,1,1), (3,2,2))
b = P.exteriorFaces(a,indices)
test.testO([b,indices],23)

# BAR
a = G.cart((0,0,0), (1,1,1), (10,1,1))
a = C.convertArray2Tetra(a); a = G.close(a)
b = P.exteriorFaces(a)
test.testA([b],9)
indices = []
b = P.exteriorFaces(a,indices)
test.testO([b,indices],29)
a = D.circle((0,0,0),1.)
a = C.convertArray2Tetra(a); a = G.close(a)
indices=[]
try:
    b = P.exteriorFaces(a,indices)
except: pass

# TRI
a = G.cartTetra((0,0,0), (1,1,1), (20,3,1))
b = P.exteriorFaces(a)
test.testA([b],4)
indices = []
b = P.exteriorFaces(a,indices)
test.testO([b,indices],24)

# QUAD
a = G.cartHexa((0,0,0), (1,1,1), (20,3,1))
b = P.exteriorFaces(a)
test.testA([b],5)
indices = []
b = P.exteriorFaces(a)
test.testO([b,indices],25)

# TETRA
a = G.cartTetra((0,0,0), (1,1,1), (3,3,3))
b = P.exteriorFaces(a)
test.testA([b],6)
indices = []
b = P.exteriorFaces(a)
test.testO([b,indices],26)

# HEXA
a = G.cartHexa((0,0,0), (1,1,1), (3,3,3))
b = P.exteriorFaces(a)
test.testA([b],7)
indices = []
b = P.exteriorFaces(a)
test.testO([b,indices],27)

# NGON
a = G.cartNGon((0,0,0), (1,1,1), (3,3,3))
b = P.exteriorFaces(a)
test.testA([b],8)
indices = []
b = P.exteriorFaces(a)
test.testO([b,indices],28)
# EMPTY
a = D.sphere( (0,0,0), 1.)
a = C.convertArray2Tetra(a)
a = G.close(a)
try:
    b = P.exteriorFaces(a)
except: pass
indices = []
try:
    b = P.exteriorFaces(a,indices)
except: pass
