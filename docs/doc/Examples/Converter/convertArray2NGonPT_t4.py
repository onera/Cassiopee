# - convertArray2NGon(pyTree): ME, api 3-
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

# 2D ME: tri-quad-tri
#             |
#            tri
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (5,10,1))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.2), (5,10,1))
c = G.cartTetra((0.8,0.,0.), (0.1,0.1,0.2), (5,10,1))
d = G.cartTetra((0.4,-0.9,0.), (0.1,0.1,0.2), (5,10,1))
a = C.mergeConnectivity([a, b, c, d], None)
a = C.convertArray2NGon(a, api=3)
t = C.newPyTree(['Base', a])
test.testT(t, 1)

# 3D ME: tetra - pyra
#          |      |
#        penta   hexa
a = G.cartTetra((0.,0.4,0.), (0.1,0.1,0.1), (5,5,5))
b = G.cartPyra((0.4,0.4,0.), (0.1,0.1,0.1), (5,5,5))
c = G.cartPenta((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
d = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = C.mergeConnectivity([a, b, c, d], None)
a = C.convertArray2NGon(a, api=3)
t = C.newPyTree(['Base', a])
test.testT(t, 2)
