# - convertArray2NGon(pyTree): BE, api 3-
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# 2D BE: tri
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a = C.convertArray2NGon(a, api=3)
t = C.newPyTree(['Base', a])
test.testT(t, 1)

# 2D BE: quad
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a = C.convertArray2NGon(a, api=3)
t = C.newPyTree(['Base', a])
test.testT(t, 2)

# 3D BE: tetra
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
a = C.convertArray2NGon(a, api=3)
t = C.newPyTree(['Base', a])
test.testT(t, 3)

# 3D BE: pyra
a = G.cartPyra((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
a = C.convertArray2NGon(a, api=3)
t = C.newPyTree(['Base', a])
test.testT(t, 4)

# 3D BE: penta
a = G.cartPenta((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
a = C.convertArray2NGon(a, api=3)
t = C.newPyTree(['Base', a])
test.testT(t, 5)

# 3D BE: hexa
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
a = C.convertArray2NGon(a, api=3)
t = C.newPyTree(['Base', a])
test.testT(t, 6)
