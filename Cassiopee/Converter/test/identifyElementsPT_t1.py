# - identifyElements (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Post.PyTree as P
import KCore.test as test

# structured
a = G.cart((0,0,0), (1,0,0), (10,1,1))
f = P.exteriorElts(a)
hook = C.createHook(a, function='elementCenters')
elts = C.identifyElements(hook, f)
C.freeHook(hook)
test.testO(elts, 1)

a = G.cart((0,0,0), (1,1,0), (10,10,1))
f = P.exteriorElts(a)
hook = C.createHook(a, function='elementCenters')
elts = C.identifyElements(hook, f)
C.freeHook(hook)
test.testO(elts, 2)

a = G.cart((0,0,0), (1,1,1), (10,10,10))
f = P.exteriorElts(a)
hook = C.createHook(a, function='elementCenters')
elts = C.identifyElements(hook, f)
C.freeHook(hook)
test.testO(elts, 3)

# 2D BE: tri
a = G.cartTetra((0,0,0), (1,1,0), (10,10,1))
f = P.exteriorElts(a)
hook = C.createHook(a, function='elementCenters')
elts = C.identifyElements(hook, f)
C.freeHook(hook)
test.testO(elts, 4)

# 3D BE: hexa
a = G.cartHexa((0,0,0), (1,1,1), (10,10,10))
f = P.exteriorElts(a)
hook = C.createHook(a, function='elementCenters')
elts = C.identifyElements(hook, f)
C.freeHook(hook)
test.testO(elts, 5)

# 3D ME: tri-quad
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (5,10,1))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.2), (5,10,1))
a = C.mergeConnectivity([a, b], None)
f = P.exteriorElts(a)
hook = C.createHook(a, function='elementCenters')
elts = C.identifyElements(hook, f)
C.freeHook(hook)
test.testO(elts, 6)

# 3D ME: pyra - penta - hexa
a = G.cartPyra((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
b = G.cartPenta((0.4,0.,0.), (0.1,0.1,0.1), (5,5,5))
c = G.cartHexa((0.8,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.mergeConnectivity([a, b, c], None)
f = P.exteriorElts(a)
hook = C.createHook(a, function='elementCenters')
elts = C.identifyElements(hook, f)
C.freeHook(hook)
test.testO(elts, 7)

# 2D NGon v3
a = G.cartNGon((0,0,0), (1,1,0), (10,10,1), api=1)
f = P.exteriorElts(a)
hook = C.createHook(a, function='elementCenters')
#elts = C.identifyElements(hook, f)
#print(elts)
C.freeHook(hook)
#test.testO(elts, 8)

# 3D NGon v3
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=1)
f = P.exteriorElts(a)
hook = C.createHook(a, function='elementCenters')
#elts = C.identifyElements(hook, f)
#print(elts)
C.freeHook(hook)
#test.testO(elts, 9)

# 2D NGon v4
a = G.cartNGon((0,0,0), (1,1,0), (10,10,1), api=3)
f = P.exteriorElts(a)
hook = C.createHook(a, function='elementCenters')
elts = C.identifyElements(hook, f)
C.freeHook(hook)
test.testO(elts, 10)

# 3D NGon v4
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10), api=3)
f = P.exteriorElts(a)
hook = C.createHook(a, function='elementCenters')
elts = C.identifyElements(hook, f)
C.freeHook(hook)
test.testO(elts, 11)
