# - dual (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Post.PyTree as P
import Intersector.PyTree as XOR
import KCore.test as test

# -- NGON, api 3
# 1D
N = 11
a = G.cartNGon((0,0,0), (1,0,0), (N,1,1), api=3)
C._initVars(a, '{F}={CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateX}')
T._dual(a)
#test.testT(a, 1) TODO WRONG REF - points missing at both ends

# 2D
a = G.cartNGon((0,0,0), (1,0,1), (N,1,N), api=3)
C._initVars(a, '{F}={CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateX}')
T._dual(a)
C._initVars(a, '{centers:H}={centers:CoordinateZ}')
test.testT(a, 2)

# 3D
a = G.cartNGon((0,0,0), (1,1,1), (N,N,N), api=3)
C._initVars(a, '{F}={CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateX}')
b = T.dual(a)
test.testT(b, 3)


## -- ME
# 2D: tri-quad-tri
a = G.cartTetra((0.,0.,0.), (0.1,0.1,0.2), (5,10,1))
b = G.cartHexa((0.4,0.,0.), (0.1,0.1,0.2), (5,10,1))
c = G.cartTetra((0.8,0.,0.), (0.1,0.1,0.2), (5,10,1))
a = T.join([a, b, c])
C._initVars(a, '{F}={CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateX}')
d = T.dual(a)
C._initVars(d, '{centers:H}={centers:CoordinateY}')
test.testT(d, 4)

# 3D: pyra - penta - hexa
indices = []
a = G.cartPyra((0.,0.,0.), (0.1,0.1,0.1), (5,5,5))
b = G.cartPenta((0.4,0.,0.), (0.1,0.1,0.1), (5,5,5))
c = G.cartHexa((0.8,0.,0.), (0.1,0.1,0.1), (5,5,5))
a = T.join([a, b, c])
C._initVars(a, '{F}={CoordinateX}')
C._initVars(a, '{centers:G}={centers:CoordinateX}')
T._dual(a)
XOR._reorient(a)
C._initVars(a, '{centers:H}={centers:CoordinateY}')
test.testT(a, 5)
