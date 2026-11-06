# - exteriorFacesStructured (pyTree) -
import Converter.PyTree as C
import Post.PyTree as P
import Generator.PyTree as G
import KCore.test as test

# 1D
a = G.cart((0,0,0), (1,1,1), (10,1,1))
A = P.exteriorFacesStructured(a)
t = C.newPyTree(['Base',0]); t[2][1][2] += A
test.testT(t,1)

# 2D
a = G.cart((0,0,0), (1,1,1), (1,6,10))
A = P.exteriorFacesStructured(a)
t = C.newPyTree(['Base',1]); t[2][1][2] += A
test.testT(t,2)

# 3D
a = G.cart((0,0,0), (1,1,1), (4,4,6))
A = P.exteriorFacesStructured(a)
t = C.newPyTree(['Base',2]); t[2][1][2] += A
test.testT(t,3)
