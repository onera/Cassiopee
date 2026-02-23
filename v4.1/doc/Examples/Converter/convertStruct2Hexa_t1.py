# - convertArray2Hexa -
import Converter as C
import Generator as G
import KCore.test as test

# 1D: BAR
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,1,1))
a = C.addVars(a, 'F')
b = C.convertArray2Hexa(a)
test.testA([b], 11)

# 2D: QUAD
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a = C.addVars(a, 'F')
b = C.convertArray2Hexa(a)
test.testA([b], 1)

# 3D: HEXA
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
a = C.addVars(a, 'F')
b = C.convertArray2Hexa(a)
test.testA([b], 2)

# On lists
a1 = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a2 = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
l = C.convertArray2Hexa([a1,a2])
test.testA(l, 3)
