# - convertArray2Tetra (array) -
import Converter as C
import Generator as G
import KCore.test as test

# 1D : BAR
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,1,1))
a = C.initVars(a, '{F}={x}+{y}')
a = C.convertArray2Tetra(a)
test.testA([a], 11)

# 2D : TRI + champ
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a = C.addVars(a, 'F')
a = C.convertArray2Tetra(a)
test.testA([a], 1)

# 3D : TETRA
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
b = C.convertArray2Tetra(a)
test.testA([b], 2)

# On lists
a1 = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a2 = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
l = C.convertArray2Tetra([a1,a2])
test.testA(l, 3)
