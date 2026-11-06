# - convertArray2Tetra -
import Converter as C
import Generator as G
import KCore.test as test

# 2D : QUAD -> TRI
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
b = C.convertArray2Tetra(a)
test.testA([b], 1)

# 3D : HEXA -> TETRA
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
b = C.convertArray2Tetra(a)
test.testA([b], 2)

# On lists
a1 = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a2 = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
l = C.convertArray2Tetra([a1,a2])
test.testA(l, 3)
