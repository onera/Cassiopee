# - exteriorFacesStructured (array) -
import KCore.test as test
import Post as P
import Generator as G

# 1D array
a = G.cart((0,0,0), (1,1,1), (10,1,1))
A = P.exteriorFacesStructured(a)
test.testA(A,1)

# 2D array
a = G.cart((0,0,0), (1,1,1), (1,6,10))
A = P.exteriorFacesStructured(a)
test.testA(A,2)

# 3D array
a = G.cart((0,0,0), (1,1,1), (4,4,6))
A = P.exteriorFacesStructured(a)
test.testA(A,3)
